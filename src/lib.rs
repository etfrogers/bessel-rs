#![feature(float_gamma)]
mod amos;
use amos::zbesj;
pub use amos::{BesselError, GammaError, Scaling};
use num::complex::Complex64;

pub fn bessel_j<T: Into<Complex64>>(order: f64, z: T) -> Result<Complex64, BesselError> {
    zbesj(z.into(), order, Scaling::Unscaled, 1).map(|v| v.0[0])
}

pub fn add(left: u64, right: u64) -> u64 {
    left + right
}

#[cfg(test)]
mod tests;
