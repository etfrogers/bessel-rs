use std::f64::consts::PI;

use num::Complex;

use crate::HankelKind;

/// (-1)^n sign factor for integer order reflection.
#[inline]
pub fn integer_sign(n: i64) -> f64 {
    if n % 2 == 0 { 1.0 } else { -1.0 }
}

/// Check if `nu` is a non-negative integer. Returns `Some(n)` if so.
#[inline]
pub fn as_integer(nu: f64) -> Option<i64> {
    if nu == nu.floor() {
        // Safe conversion: orders beyond i64 range are not practical
        Some(nu as i64)
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
pub(crate) fn sinpi(x: f64) -> f64 {
    // sinpi is odd: sinpi(-x) = -sinpi(x)
    let (ax, sign) = if x < 0.0 { (-x, -1.0) } else { (x, 1.0) };

    // Reduce to [0, 2): r = ax mod 2
    let r = ax % 2.0;

    // Exact special values
    if r == 0.0 || r == 1.0 {
        return 0.0;
    }
    if r == 0.5 {
        return sign;
    }
    if r == 1.5 {
        return -sign;
    }

    // Use symmetry to reduce to [0, 0.5]
    let s = if r < 0.5 {
        (r * PI).sin()
    } else if r < 1.0 {
        ((1.0 - r) * PI).sin()
    } else if r < 1.5 {
        -((r - 1.0) * PI).sin()
    } else {
        -((2.0 - r) * PI).sin()
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
pub(crate) fn cospi(x: f64) -> f64 {
    // cospi is even: cospi(-x) = cospi(x)
    let ax = x.abs();

    // Reduce to [0, 2): r = ax mod 2
    let r = ax % 2.0;

    // Exact special values
    if r == 0.0 {
        return 1.0;
    }
    if r == 0.5 || r == 1.5 {
        return 0.0;
    }
    if r == 1.0 {
        return -1.0;
    }

    // Use symmetry to reduce to [0, 0.5]
    if r < 0.5 {
        (r * PI).cos()
    } else if r < 1.0 {
        -((1.0 - r) * PI).cos()
    } else if r < 1.5 {
        -((r - 1.0) * PI).cos()
    } else {
        ((2.0 - r) * PI).cos()
    }
}

/// J_{-ν}(z) = cos(νπ)·J_ν(z) − sin(νπ)·Y_ν(z)  (DLMF 10.2.3)
#[inline]
pub fn reflect_j_element(order: f64, j: Complex<f64>, y: Complex<f64>) -> Complex<f64> {
    j * cospi(order) - y * sinpi(order)
}

/// H^(m)_{-ν}(z) = exp(±νπi)·H^(m)_ν(z)  (DLMF 10.4.6/7)
#[inline]
pub fn reflect_h_element(order: f64, kind: HankelKind, h: Complex<f64>) -> Complex<f64> {
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
pub fn reflect_y_element(order: f64, j: Complex<f64>, y: Complex<f64>) -> Complex<f64> {
    j * sinpi(order) + y * cospi(order)
}

/// I_{-ν}(z) = I_ν(z) + (2/π)·sin(νπ)·K_ν(z)  (DLMF 10.27.2)
#[inline]
pub fn reflect_i_element(order: f64, i: Complex<f64>, k: Complex<f64>) -> Complex<f64> {
    k * (2.0 / PI * sinpi(order)) + i
}
