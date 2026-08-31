use num::Complex;

use crate::{
    BesselError, HankelKind, Scaling,
    types::{BesselFloat, BesselResult, BesselValues},
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

#[allow(type_alias_bounds)]
pub(crate) type BesselSig<T: BesselFloat = f64> =
    fn(Complex<T>, T, Scaling, usize) -> Result<BesselValues<T>, BesselError<T>>;

pub(crate) fn no_secondary<T: BesselFloat>() -> Option<(BesselSig<T>, UnderflowLocation)> {
    None
}

pub(crate) fn reflect_orders<
    T: BesselFloat,
    BesselSigPrim,
    BesselSigSec,
    FReflectNonInt,
    FReflectInt,
>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    n: usize,
    primary_fn: BesselSigPrim,
    primary_loc: UnderflowLocation,
    negative_fn: Option<(BesselSigSec, UnderflowLocation)>,
    reflect_non_int: FReflectNonInt,
    reflect_int: FReflectInt,
) -> BesselResult<T>
where
    BesselSigPrim: Fn(Complex<T>, T, Scaling, usize) -> Result<BesselValues<T>, BesselError<T>>,
    BesselSigSec: Fn(Complex<T>, T, Scaling, usize) -> Result<BesselValues<T>, BesselError<T>>,
    FReflectNonInt: Fn(T, Complex<T>, Option<Complex<T>>) -> Complex<T>,
    FReflectInt: Fn(i64, Complex<T>) -> Complex<T>,
{
    let abs_order: T = order.abs();
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
        let result = primary_fn(z, T::from_isize(min_order as isize), scaling, n_positive)?;
        (result, None)
    } else {
        // General case: need both J and Y at positive |ν|
        //
        // say order = -2.7, then we need
        // -2.7, -1.7, -0.7, 0.3, 1.3, 2.3, ...

        let n_negative = order.abs().ceil().abs().to_usize().unwrap();
        let first_negative = order.fract().abs();
        let primary_neg_result = primary_fn(z, first_negative, scaling, n_negative)?;
        let secondary_neg_result = if let Some((ref negative_fn, _)) = negative_fn {
            Some(negative_fn(z, first_negative, scaling, n_negative)?)
        } else {
            None
        };

        let n_positive = n.saturating_sub(n_negative);

        let first_positive = T::ONE + order.fract();
        let pos_result = if n_positive > 0 {
            primary_fn(z, first_positive, scaling, n_positive)?
        } else {
            (Vec::new(), 0)
        };
        (pos_result, Some((primary_neg_result, secondary_neg_result)))
    };

    let (prim_positive, n_zeros_prim_positive) = pos_result;
    let n_positive = prim_positive.len();
    let mut answer = Vec::with_capacity(n);
    let mut n_zeros = 0;

    let n_remaining = if negative_data.is_none() {
        // Special case for negative integer order: J(-n, z) = (-1)^n J(n, z)
        let mut n_negative = 0;
        for i in 0..n {
            let current_order = order + T::from_usize(i);
            let abs_order = current_order.abs().to_usize().unwrap();
            if primary_loc.is_underflow(abs_order, n_positive, n_zeros_prim_positive) {
                n_zeros += 1;
            }
            if current_order < T::ZERO {
                answer.push(reflect_int(abs_order as i64, prim_positive[abs_order]));
                n_negative += 1;
            } else {
                // abort on first positive order
                break;
            }
        }
        n - n_negative
    } else {
        let ((prim_negative, n_zeros_prim_neg), secondary_result) = negative_data.unwrap();
        let (sec_negative, n_zeros_sec_neg) = if let Some((vals, nz)) = secondary_result {
            (Some(vals), Some(nz))
        } else {
            (None, None)
        };
        let n_negative = prim_negative.len();
        for i in 0..n {
            let current_order = order + T::from_usize(i);
            if current_order < T::ZERO {
                let index = n_negative - 1 - i;
                answer.push(reflect_non_int(
                    current_order.abs(),
                    prim_negative[index],
                    sec_negative.as_ref().map(|y| y[index]),
                ));
                if primary_loc.is_underflow(index, n_negative, n_zeros_prim_neg)
                    && negative_fn.as_ref().is_some_and(|(_, loc)| {
                        loc.is_underflow(index, n_negative, n_zeros_sec_neg.unwrap())
                    })
                {
                    n_zeros += 1;
                }
            } else {
                // abort on first positive order
                break;
            }
        }
        n_positive
    };

    if n_remaining > 0 {
        // make j_positive mutable for draining - was immutable to this point
        let mut j_positive = prim_positive;
        answer.extend(j_positive.drain(..n_remaining));
        n_zeros += primary_loc.positive_tail_zeros(n_positive, n_remaining, n_zeros_prim_positive);
    }
    Ok((answer, n_zeros))
}
