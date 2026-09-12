use num::Complex;

use crate::{
    BesselError, HankelKind, Scaling,
    amos::{algorithms, validate_inputs},
    prelude::*,
    types::{BesselFloat, BesselResult},
};

/// (-1)^n sign factor for integer order reflection.
#[inline]
pub(crate) fn integer_sign<T: BesselFloat>(n: i64) -> T {
    if n % 2 == 0 { T::one() } else { -T::one() }
}

/// Check if `nu` is a non-negative integer. Returns `Some(n)` if so.
#[inline]
pub(crate) fn as_integer<T: BesselFloat>(nu: T) -> Option<i64> {
    if nu.is_finite() && nu == nu.floor() {
        // Safe conversion: orders beyond i64 range are not practical
        nu.to_i64()
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
pub(crate) fn reflect_j_element<T: BesselFloat>(
    order: T,
    j: Complex<T>,
    y: Complex<T>,
) -> Complex<T> {
    j * cospi(order) - y * sinpi(order)
}

/// H^(m)_{-ν}(z) = exp(±νπi)·H^(m)_ν(z)  (DLMF 10.4.6/7)
#[inline]
pub(crate) fn reflect_h_element<T: BesselFloat>(
    order: T,
    kind: HankelKind,
    h: Complex<T>,
) -> Complex<T> {
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
pub(crate) fn reflect_y_element<T: BesselFloat>(
    order: T,
    j: Complex<T>,
    y: Complex<T>,
) -> Complex<T> {
    j * sinpi(order) + y * cospi(order)
}

/// I_{-ν}(z) = I_ν(z) + (2/π)·sin(νπ)·K_ν(z)  (DLMF 10.27.2)
///
/// When `scaling == Scaling::Scaled`, K_ν is scaled by exp(z) while I_ν is scaled by exp(-|Re(z)|).
/// The K_ν term must be converted to the I_ν scaling frame by multiplying by exp(-|Re(z)| - z).
#[inline]
pub(crate) fn reflect_i_element<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    i: Complex<T>,
    k: Complex<T>,
) -> Complex<T> {
    let k_scaled = match scaling {
        Scaling::Unscaled => k,
        Scaling::Scaled => {
            let x = z.re;
            let r = -x.abs() - x; // <= 0 for all x
            let mc = T::MACHINE_CONSTANTS;
            if -r > mc.exponent_limit {
                Complex::<T>::ZERO
            } else {
                (k * r.exp()) * Complex::<T>::cis(-z.im)
            }
        }
    };
    k_scaled * (T::TWO / T::PI() * sinpi(order)) + i
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) enum UnderflowLocation {
    Start,
    End,
}

impl UnderflowLocation {
    #[inline]
    pub(crate) fn slice_zeros(self, len: usize, start: usize, end: usize, n_zeros: usize) -> usize {
        let slice_len = end + 1 - start;
        match self {
            UnderflowLocation::Start => n_zeros.saturating_sub(start).min(slice_len),
            UnderflowLocation::End => {
                let tail = len.saturating_sub(end + 1);
                n_zeros.saturating_sub(tail).min(slice_len)
            }
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
        z: Complex<T>,
        order: T,
        scaling: Scaling,
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
        _z: Complex<T>,
        _order: T,
        _scaling: Scaling,
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
        algorithms::complex_bessel_j(z, order, scaling, n)
    }

    #[inline]
    fn reflect_non_int(
        &self,
        _z: Complex<T>,
        order: T,
        _scaling: Scaling,
        j: Complex<T>,
        y: Option<Complex<T>>,
    ) -> Complex<T> {
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
        algorithms::complex_bessel_y(z, order, scaling, n)
    }

    #[inline]
    fn reflect_non_int(
        &self,
        _z: Complex<T>,
        order: T,
        _scaling: Scaling,
        y: Complex<T>,
        j: Option<Complex<T>>,
    ) -> Complex<T> {
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
        algorithms::complex_bessel_i(z, order, scaling, n)
    }

    #[inline]
    fn reflect_non_int(
        &self,
        z: Complex<T>,
        order: T,
        scaling: Scaling,
        i: Complex<T>,
        k: Option<Complex<T>>,
    ) -> Complex<T> {
        reflect_i_element(z, order, scaling, i, k.unwrap())
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
        algorithms::complex_bessel_k(z, order, scaling, n)
    }

    #[inline]
    fn reflect_non_int(
        &self,
        _z: Complex<T>,
        _order: T,
        _scaling: Scaling,
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
        algorithms::complex_bessel_h(z, order, scaling, self.0, n)
    }

    #[inline]
    fn reflect_non_int(
        &self,
        _z: Complex<T>,
        order: T,
        _scaling: Scaling,
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
    validate_inputs(z, order, n)?;
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

    let finish = |y: Vec<Complex<T>>, n_zeros: usize, plos: bool| {
        if plos {
            Err(BesselError::PartialLossOfSignificance { y, n_zeros })
        } else {
            Ok((y, n_zeros))
        }
    };

    let abs_order: T = order.abs();
    let n_order = abs_order.ceil().to_usize().unwrap();
    let n_negative = n_order.min(n);

    // 1. Negative integer orders: J(-n, z) = (-1)^n J(n, z)
    // Evaluated with a single positive Amos call starting at order 0.
    if let Some(int_order) = as_integer(abs_order) {
        let max_order = (n as i64 - 1 - int_order).max(int_order);
        let n_positive = (max_order + 1) as usize;
        let (mut pos_values, pos_n_zeros) = unwrap_plos(op.eval(z, T::ZERO, scaling, n_positive))?;

        let order_size = int_order as usize;
        let start_ind = order_size + 1 - n_negative;
        let mut n_zeros =
            Op::UNDERFLOW_LOCATION.slice_zeros(n_positive, start_ind, order_size, pos_n_zeros);

        let mut answer = Vec::with_capacity(n);
        for i in 0..n_negative {
            let cur_order = order_size - i;
            answer.push(op.reflect_int(cur_order as i64, pos_values[cur_order]));
        }

        let n_remaining = n - n_negative;
        if n_remaining > 0 {
            answer.extend(pos_values.drain(..n_remaining));
            n_zeros +=
                Op::UNDERFLOW_LOCATION.slice_zeros(n_positive, 0, n_remaining - 1, pos_n_zeros);
        }

        return finish(answer, n_zeros, partial_loss_of_significance);
    }

    // 2. Negative non-integer orders (DLMF reflection formulas)
    let first_negative = order.abs() - T::from_usize(n_negative - 1);
    let (prim_neg, n_zeros_prim_neg) =
        unwrap_plos(op.eval(z, first_negative, scaling, n_negative))?;
    let sec_neg_result = op
        .secondary()
        .map(|s| unwrap_plos(s.eval(z, first_negative, scaling, n_negative)))
        .transpose()?;

    let (sec_neg, n_zeros_sec_neg) = sec_neg_result.unzip();
    let secondary_neg_iter = sec_neg.map(|sec| sec.into_iter().rev());

    let mut answer = Vec::with_capacity(n);
    for (i, (prim_val, sec_val)) in prim_neg
        .into_iter()
        .rev()
        .zip_option(secondary_neg_iter)
        .enumerate()
    {
        let cur_abs_order = order.abs() - T::from_usize(i);
        answer.push(op.reflect_non_int(z, cur_abs_order, scaling, prim_val, sec_val));
    }

    let mut n_zeros = match n_zeros_sec_neg {
        Some(sec_zeros) => (n_zeros_prim_neg + sec_zeros).saturating_sub(n_negative),
        None => n_zeros_prim_neg,
    };

    // Push the remaining positive orders onto the end, if required.
    let n_remaining = n - n_negative;
    if n_remaining > 0 {
        let first_positive = order + T::from_usize(n_negative);
        let (pos_values, pos_n_zeros) =
            unwrap_plos(op.eval(z, first_positive, scaling, n_remaining))?;
        answer.extend(pos_values);
        n_zeros += pos_n_zeros;
    }

    finish(answer, n_zeros, partial_loss_of_significance)
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
