mod amos;
mod bessel_zeros;

pub use amos::{BesselError, GammaError, Scaling};
use amos::{
    complex_airy, complex_airy_b, complex_bessel_h, complex_bessel_i, complex_bessel_j,
    complex_bessel_k, complex_bessel_y,
};
use num::complex::Complex64;

use crate::amos::HankelKind;
pub use bessel_zeros::{BesselFunType, bessel_zeros};

// TODO work with abritrary bit-depth floats
// TODO bessel derivatives
// TODO negative orders

pub fn bessel_j<T: Into<Complex64>>(order: f64, z: T) -> Result<Complex64, BesselError> {
    complex_bessel_j(z.into(), order, Scaling::Unscaled, 1).map(|v| v.0[0])
}

pub fn bessel_i<T: Into<Complex64>>(order: f64, z: T) -> Result<Complex64, BesselError> {
    complex_bessel_i(z.into(), order, Scaling::Unscaled, 1).map(|v| v.0[0])
}

pub fn bessel_k<T: Into<Complex64>>(order: f64, z: T) -> Result<Complex64, BesselError> {
    complex_bessel_k(z.into(), order, Scaling::Unscaled, 1).map(|v| v.0[0])
}

pub fn bessel_y<T: Into<Complex64>>(order: f64, z: T) -> Result<Complex64, BesselError> {
    complex_bessel_y(z.into(), order, Scaling::Unscaled, 1).map(|v| v.0[0])
}

pub fn hankel<T: Into<Complex64>>(
    order: f64,
    z: T,
    kind: HankelKind,
) -> Result<Complex64, BesselError> {
    complex_bessel_h(z.into(), order, Scaling::Unscaled, kind, 1).map(|v| v.0[0])
}

pub fn airy<T: Into<Complex64>>(z: T) -> Result<Complex64, BesselError> {
    complex_airy(z.into(), false, Scaling::Unscaled).map(|v| v.0)
}

pub fn airyp<T: Into<Complex64>>(z: T) -> Result<Complex64, BesselError> {
    complex_airy(z.into(), true, Scaling::Unscaled).map(|v| v.0)
}

pub fn airy_b<T: Into<Complex64>>(z: T) -> Result<Complex64, BesselError> {
    complex_airy_b(z.into(), false, Scaling::Unscaled).map(|v| v)
}

pub fn airy_bp<T: Into<Complex64>>(z: T) -> Result<Complex64, BesselError> {
    complex_airy_b(z.into(), true, Scaling::Unscaled)
}

#[cfg(test)]
mod tests;
