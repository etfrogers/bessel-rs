#![feature(test)]
mod amos;
pub use amos::{BesselError, GammaError, Scaling};
use amos::{complex_bessel_h, complex_bessel_j, zbesi};
use num::complex::Complex64;

use crate::amos::HankelKind;

// TODO work with abritrary bit-depth floats

pub fn bessel_j<T: Into<Complex64>>(order: f64, z: T) -> Result<Complex64, BesselError> {
    complex_bessel_j(z.into(), order, Scaling::Unscaled, 1).map(|v| v.0[0])
}

pub fn bessel_i<T: Into<Complex64>>(order: f64, z: T) -> Result<Complex64, BesselError> {
    zbesi(z.into(), order, Scaling::Unscaled, 1).map(|v| v.0[0])
}

pub fn hankel<T: Into<Complex64>>(
    order: f64,
    z: T,
    kind: HankelKind,
) -> Result<Complex64, BesselError> {
    complex_bessel_h(z.into(), order, Scaling::Unscaled, kind, 1).map(|v| v.0[0])
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod bench;
