use num::Complex;

use crate::{
    BesselError, HankelKind, Scaling,
    amos::core,
    types::{BesselFloat, BesselResult},
};

/// (-1)^n sign factor for integer order reflection.
#[inline]
pub fn integer_sign<T: BesselFloat>(n: i64) -> T {
    if n % 2 == 0 { T::one() } else { -T::one() }
}

/// Check if `nu` is a non-negative integer. Returns `Some(n)` if so.
#[inline]
pub fn as_integer<T: BesselFloat>(nu: T) -> Option<i64> {
    if nu == nu.floor() {
        // Safe conversion: orders beyond i64 range are not practical
        Some(nu.to_i64().unwrap())
    } else {
        None
    }
}

/// Compute sin(π·x) with exact values at half-integers.
///
/// Reduces the argument modulo 2 first, so `sinpi(n)` is exactly 0 for
/// any integer `n`, and `sinpi(n + 0.5)` is exactly ±1. This avoids the
/// catastrophic rounding errors of `(x * PI).sin()` when x is a
/// half-integer (e.g. `sin(1.5 * PI)` = −1.837e-16 instead of 0).
///
/// Algorithm follows scipy/xsf: reduce to [0, 0.5], use symmetry.
#[inline]
pub(crate) fn sinpi<T: BesselFloat>(x: T) -> T {
    // sinpi is odd: sinpi(-x) = -sinpi(x)
    let (ax, sign) = if x < T::ZERO {
        (-x, -T::one())
    } else {
        (x, T::one())
    };

    // Reduce to [0, 2): r = ax mod 2
    let r = ax % T::TWO;

    // Exact special values
    if r == T::ZERO || r == T::one() {
        return T::ZERO;
    }
    if r == T::HALF {
        return sign;
    }
    if r == T::from_f64(1.5) {
        return -sign;
    }

    // Use symmetry to reduce to [0, 0.5]
    let s = if r < T::HALF {
        (r * T::PI()).sin()
    } else if r < T::one() {
        ((T::one() - r) * T::PI()).sin()
    } else if r < T::from_f64(1.5) {
        -((r - T::one()) * T::PI()).sin()
    } else {
        -((T::TWO - r) * T::PI()).sin()
    };

    sign * s
}

/// Compute cos(π·x) with exact values at integers and half-integers.
///
/// Reduces the argument modulo 2 first, so `cospi(n + 0.5)` is exactly 0
/// for any integer `n`, and `cospi(n)` is exactly ±1. This avoids the
/// catastrophic rounding errors of `(x * PI).cos()` when x is a
/// half-integer (e.g. `cos(1.5 * PI)` = −1.837e-16 instead of 0).
///
/// Algorithm follows scipy/xsf: reduce to [0, 0.5], use symmetry.
#[inline]
pub(crate) fn cospi<T: BesselFloat>(x: T) -> T {
    // cospi is even: cospi(-x) = cospi(x)
    let ax = x.abs();

    // Reduce to [0, 2): r = ax mod 2
    let r = ax % T::TWO;

    // Exact special values
    if r == T::ZERO {
        return T::one();
    }
    if r == T::HALF || r == T::from_f64(1.5) {
        return T::ZERO;
    }
    if r == T::one() {
        return -T::one();
    }

    // Use symmetry to reduce to [0, 0.5]
    if r < T::HALF {
        (r * T::PI()).cos()
    } else if r < T::one() {
        -((T::one() - r) * T::PI()).cos()
    } else if r < T::from_f64(1.5) {
        -((r - T::one()) * T::PI()).cos()
    } else {
        ((T::TWO - r) * T::PI()).cos()
    }
}

/// J_{-ν}(z) = cos(νπ)·J_ν(z) − sin(νπ)·Y_ν(z)  (DLMF 10.2.3)
#[inline]
pub fn reflect_j_element<T: BesselFloat>(order: T, j: Complex<T>, y: Complex<T>) -> Complex<T> {
    j * cospi(order) - y * sinpi(order)
}

/// H^(m)_{-ν}(z) = exp(±νπi)·H^(m)_ν(z)  (DLMF 10.4.6/7)
#[inline]
pub fn reflect_h_element<T: BesselFloat>(order: T, kind: HankelKind, h: Complex<T>) -> Complex<T> {
    let cos_nu_pi = cospi(order);
    let sin_nu_pi = sinpi(order);
    let rotation = match kind {
        HankelKind::First => Complex::new(cos_nu_pi, sin_nu_pi),
        HankelKind::Second => Complex::new(cos_nu_pi, -sin_nu_pi),
    };
    h * rotation
}

/// Y_{-ν}(z) = sin(νπ)·J_ν(z) + cos(νπ)·Y_ν(z)  (DLMF 10.2.3)
#[inline]
pub fn reflect_y_element<T: BesselFloat>(order: T, j: Complex<T>, y: Complex<T>) -> Complex<T> {
    j * sinpi(order) + y * cospi(order)
}

/// I_{-ν}(z) = I_ν(z) + (2/π)·sin(νπ)·K_ν(z)  (DLMF 10.27.2)
#[inline]
pub fn reflect_i_element<T: BesselFloat>(order: T, i: Complex<T>, k: Complex<T>) -> Complex<T> {
    k * (T::TWO / T::PI() * sinpi(order)) + i
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) enum UnderflowLocation {
    Start,
    End,
}

impl UnderflowLocation {
    /// Check if the element at `index` in a buffer of `len` underflowed.
    #[inline]
    pub fn is_underflow(self, index: usize, len: usize, n_zeros: usize) -> bool {
        match self {
            UnderflowLocation::Start => index < n_zeros,
            UnderflowLocation::End => index >= len.saturating_sub(n_zeros),
        }
    }

    /// Calculate how many zeros remain when taking the first `n_remaining` elements.
    #[inline]
    pub fn positive_tail_zeros(self, len: usize, n_remaining: usize, n_zeros: usize) -> usize {
        match self {
            UnderflowLocation::Start => n_zeros.min(n_remaining),
            UnderflowLocation::End => n_zeros.saturating_sub(len.saturating_sub(n_remaining)),
        }
    }
}

pub(crate) trait ReflectableBessel<T: BesselFloat> {
    /// Location where zeros appear upon underflow (Start or End).
    const UNDERFLOW_LOCATION: UnderflowLocation;

    /// The secondary function type needed for non-integer reflection (e.g. BesselY for BesselJ).
    type Secondary: ReflectableBessel<T>;

    /// Returns the secondary function instance, or None if none is needed (e.g. for K and H).
    fn secondary(&self) -> Option<Self::Secondary>;

    /// Evaluates the core Amos function for positive orders.
    fn eval(&self, z: Complex<T>, order: T, scaling: Scaling, n: usize) -> BesselResult<T>;

    /// DLMF reflection formula for non-integer orders: f_{-ν}(z) from f_ν(z) and optional g_ν(z).
    fn reflect_non_int(
        &self,
        order: T,
        primary: Complex<T>,
        secondary: Option<Complex<T>>,
    ) -> Complex<T>;

    /// DLMF reflection formula for integer orders: f_{-n}(z) from f_n(z).
    fn reflect_int(&self, order: i64, primary: Complex<T>) -> Complex<T>;
}

pub(crate) struct BesselJ;
pub(crate) struct BesselY;
pub(crate) struct BesselI;
pub(crate) struct BesselK;
pub(crate) struct Hankel(pub HankelKind);
pub(crate) struct NoSecondary;

impl<T: BesselFloat> ReflectableBessel<T> for NoSecondary {
    const UNDERFLOW_LOCATION: UnderflowLocation = UnderflowLocation::Start;
    type Secondary = NoSecondary;

    #[inline]
    fn secondary(&self) -> Option<Self::Secondary> {
        None
    }

    #[inline]
    fn eval(&self, _z: Complex<T>, _order: T, _scaling: Scaling, _n: usize) -> BesselResult<T> {
        unreachable!("NoSecondary should never be evaluated directly")
    }

    #[inline]
    fn reflect_non_int(
        &self,
        _order: T,
        _primary: Complex<T>,
        _secondary: Option<Complex<T>>,
    ) -> Complex<T> {
        unreachable!("NoSecondary has no reflection formula")
    }

    #[inline]
    fn reflect_int(&self, _order: i64, _primary: Complex<T>) -> Complex<T> {
        unreachable!("NoSecondary has no reflection formula")
    }
}

impl<T: BesselFloat> ReflectableBessel<T> for BesselJ {
    const UNDERFLOW_LOCATION: UnderflowLocation = UnderflowLocation::End;
    type Secondary = BesselY;

    #[inline]
    fn secondary(&self) -> Option<Self::Secondary> {
        Some(BesselY)
    }

    #[inline]
    fn eval(&self, z: Complex<T>, order: T, scaling: Scaling, n: usize) -> BesselResult<T> {
        core::complex_bessel_j(z, order, scaling, n)
    }

    #[inline]
    fn reflect_non_int(&self, order: T, j: Complex<T>, y: Option<Complex<T>>) -> Complex<T> {
        reflect_j_element(order, j, y.unwrap())
    }

    #[inline]
    fn reflect_int(&self, order: i64, j: Complex<T>) -> Complex<T> {
        j * integer_sign::<T>(order)
    }
}

impl<T: BesselFloat> ReflectableBessel<T> for BesselY {
    const UNDERFLOW_LOCATION: UnderflowLocation = UnderflowLocation::Start;
    type Secondary = BesselJ;

    #[inline]
    fn secondary(&self) -> Option<Self::Secondary> {
        Some(BesselJ)
    }

    #[inline]
    fn eval(&self, z: Complex<T>, order: T, scaling: Scaling, n: usize) -> BesselResult<T> {
        core::complex_bessel_y(z, order, scaling, n)
    }

    #[inline]
    fn reflect_non_int(&self, order: T, y: Complex<T>, j: Option<Complex<T>>) -> Complex<T> {
        reflect_y_element(order, j.unwrap(), y)
    }

    #[inline]
    fn reflect_int(&self, order: i64, y: Complex<T>) -> Complex<T> {
        y * integer_sign::<T>(order)
    }
}

impl<T: BesselFloat> ReflectableBessel<T> for BesselI {
    const UNDERFLOW_LOCATION: UnderflowLocation = UnderflowLocation::End;
    type Secondary = BesselK;

    #[inline]
    fn secondary(&self) -> Option<Self::Secondary> {
        Some(BesselK)
    }

    #[inline]
    fn eval(&self, z: Complex<T>, order: T, scaling: Scaling, n: usize) -> BesselResult<T> {
        core::complex_bessel_i(z, order, scaling, n)
    }

    #[inline]
    fn reflect_non_int(&self, order: T, i: Complex<T>, k: Option<Complex<T>>) -> Complex<T> {
        reflect_i_element(order, i, k.unwrap())
    }

    #[inline]
    fn reflect_int(&self, _order: i64, i: Complex<T>) -> Complex<T> {
        i
    }
}

impl<T: BesselFloat> ReflectableBessel<T> for BesselK {
    const UNDERFLOW_LOCATION: UnderflowLocation = UnderflowLocation::Start;
    type Secondary = NoSecondary;

    #[inline]
    fn secondary(&self) -> Option<Self::Secondary> {
        None
    }

    #[inline]
    fn eval(&self, z: Complex<T>, order: T, scaling: Scaling, n: usize) -> BesselResult<T> {
        core::complex_bessel_k(z, order, scaling, n)
    }

    #[inline]
    fn reflect_non_int(
        &self,
        _order: T,
        k: Complex<T>,
        _secondary: Option<Complex<T>>,
    ) -> Complex<T> {
        k
    }

    #[inline]
    fn reflect_int(&self, _order: i64, k: Complex<T>) -> Complex<T> {
        k
    }
}

impl<T: BesselFloat> ReflectableBessel<T> for Hankel {
    const UNDERFLOW_LOCATION: UnderflowLocation = UnderflowLocation::Start;
    type Secondary = NoSecondary;

    #[inline]
    fn secondary(&self) -> Option<Self::Secondary> {
        None
    }

    #[inline]
    fn eval(&self, z: Complex<T>, order: T, scaling: Scaling, n: usize) -> BesselResult<T> {
        core::complex_bessel_h(z, order, scaling, self.0, n)
    }

    #[inline]
    fn reflect_non_int(
        &self,
        order: T,
        h: Complex<T>,
        _secondary: Option<Complex<T>>,
    ) -> Complex<T> {
        reflect_h_element(order, self.0, h)
    }

    #[inline]
    fn reflect_int(&self, order: i64, h: Complex<T>) -> Complex<T> {
        h * integer_sign::<T>(order)
    }
}

pub(crate) fn reflect_orders<T: BesselFloat, Op: ReflectableBessel<T>>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    n: usize,
    op: Op,
) -> BesselResult<T> {
    if order >= T::ZERO {
        return op.eval(z, order, scaling, n);
    }

    let mut partial_loss_of_significance = false;

    let mut unwrap_plos = |result: BesselResult<T>| match result {
        Ok(vals) => Ok(vals),
        Err(BesselError::PartialLossOfSignificance { y, n_zeros }) => {
            partial_loss_of_significance = true;
            Ok((y, n_zeros))
        }
        Err(e) => Err(e),
    };

    let abs_order: T = order.abs();
    let n_order = abs_order.ceil().to_usize().unwrap();
    let n_negative = n_order.min(n);

    let (pos_result, negative_data) = if let Some(int_order) = as_integer(abs_order) {
        // if we have a negative integer order, then orders are
        // (if int_order is denoted by o)
        // -o, -o+1, -o+2, ..., 0, 1, 2, ... n-(o+1)
        // e.g for order = -3, int_order = 3, n = 5 orders are -3, -2, -1, 0, 1,
        // of course, if n < int_order, then we never reach order 0
        // -o, -o+1, -o+2, ..., -o+n
        // e.g. if n = 2, int_order = 3, orders are -3, -2,

        // now we need positive forms of all the orders we need in either
        // negative or positive form.

        let n64 = n as i64;
        let min_order = (-int_order + n64).abs().min(0);
        let max_order = (n64 - (int_order + 1)).max(int_order);
        let n_positive = ((max_order - min_order) + 1) as usize;
        let result =
            unwrap_plos(op.eval(z, T::from_isize(min_order as isize), scaling, n_positive))?;
        (result, None)
    } else {
        // General case: need both J and Y at positive |ν|
        //
        // say order = -2.7, then we need
        // -2.7, -1.7, -0.7, 0.3, 1.3, 2.3, ...

        let first_negative = order.abs() - T::from_usize(n_negative - 1);

        let primary_neg_result = unwrap_plos(op.eval(z, first_negative, scaling, n_negative))?;
        let secondary_neg_result = op
            .secondary()
            .map(|s| unwrap_plos(s.eval(z, first_negative, scaling, n_negative)))
            .transpose()?;

        let n_positive = n.saturating_sub(n_negative);
        let pos_result = if n_positive > 0 {
            let first_positive = T::ONE + order.fract();
            unwrap_plos(op.eval(z, first_positive, scaling, n_positive))?
        } else {
            (Vec::new(), 0)
        };
        (pos_result, Some((primary_neg_result, secondary_neg_result)))
    };

    let (prim_positive, n_zeros_prim_positive) = pos_result;
    let n_positive = prim_positive.len();
    let mut answer = Vec::with_capacity(n);
    let mut n_zeros = 0;

    if let Some(((prim_negative, n_zeros_prim_neg), secondary_result)) = negative_data {
        let (sec_negative, n_zeros_sec_neg) = secondary_result.unzip();

        let prim_rev = prim_negative.into_iter().rev();
        let sec_rev = sec_negative.map(|sec| sec.into_iter().rev());

        for (i, (prim_val, sec_val)) in prim_rev.zip_option(sec_rev).enumerate() {
            let current_order = order + T::from_usize(i);
            answer.push(op.reflect_non_int(current_order.abs(), prim_val, sec_val));
            let index = n_negative - 1 - i;
            if Op::UNDERFLOW_LOCATION.is_underflow(index, n_negative, n_zeros_prim_neg)
                && (op.secondary().is_none()
                    || Op::Secondary::UNDERFLOW_LOCATION.is_underflow(
                        index,
                        n_negative,
                        n_zeros_sec_neg.unwrap(),
                    ))
            {
                n_zeros += 1;
            }
        }
    } else {
        // Special case for negative integer order: J(-n, z) = (-1)^n J(n, z)
        for i in 0..n_negative {
            let current_order = order + T::from_usize(i);
            let abs_order = current_order.abs().to_usize().unwrap();
            if Op::UNDERFLOW_LOCATION.is_underflow(abs_order, n_positive, n_zeros_prim_positive) {
                n_zeros += 1;
            }
            answer.push(op.reflect_int(abs_order as i64, prim_positive[abs_order]));
        }
    }

    let n_remaining = n - n_negative;
    if n_remaining > 0 {
        let mut j_positive = prim_positive;
        answer.extend(j_positive.drain(..n_remaining));
        n_zeros += Op::UNDERFLOW_LOCATION.positive_tail_zeros(
            n_positive,
            n_remaining,
            n_zeros_prim_positive,
        );
    }
    if partial_loss_of_significance {
        Err(BesselError::PartialLossOfSignificance { y: answer, n_zeros })
    } else {
        Ok((answer, n_zeros))
    }
}

trait ZipOptionExt: Iterator + Sized {
    fn zip_option<J: Iterator>(
        self,
        mut maybe_iter: Option<J>,
    ) -> impl Iterator<Item = (Self::Item, Option<J::Item>)> {
        self.map(move |val| (val, maybe_iter.as_mut().and_then(|it| it.next())))
    }
}

impl<I: Iterator> ZipOptionExt for I {}
