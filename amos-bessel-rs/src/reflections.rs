use num::Complex;

use crate::{HankelKind, types::BesselFloat};

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
    let (ax, sign) = if x < T::zero() {
        (-x, -T::one())
    } else {
        (x, T::one())
    };

    // Reduce to [0, 2): r = ax mod 2
    let r = ax % T::two();

    // Exact special values
    if r == T::zero() || r == T::one() {
        return T::zero();
    }
    if r == T::half() {
        return sign;
    }
    if r == T::from_f64(1.5) {
        return -sign;
    }

    // Use symmetry to reduce to [0, 0.5]
    let s = if r < T::half() {
        (r * T::PI()).sin()
    } else if r < T::one() {
        ((T::one() - r) * T::PI()).sin()
    } else if r < T::from_f64(1.5) {
        -((r - T::one()) * T::PI()).sin()
    } else {
        -((T::two() - r) * T::PI()).sin()
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
    let r = ax % T::two();

    // Exact special values
    if r == T::zero() {
        return T::one();
    }
    if r == T::half() || r == T::from_f64(1.5) {
        return T::zero();
    }
    if r == T::one() {
        return -T::one();
    }

    // Use symmetry to reduce to [0, 0.5]
    if r < T::half() {
        (r * T::PI()).cos()
    } else if r < T::one() {
        -((T::one() - r) * T::PI()).cos()
    } else if r < T::from_f64(1.5) {
        -((r - T::one()) * T::PI()).cos()
    } else {
        ((T::two() - r) * T::PI()).cos()
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
    k * (T::two() / T::PI() * sinpi(order)) + i
}
