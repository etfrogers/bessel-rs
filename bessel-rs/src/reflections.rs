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

fn cospi(x: f64) -> f64 {
    (x * PI).cos()
}

fn sinpi(x: f64) -> f64 {
    (x * PI).sin()
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
