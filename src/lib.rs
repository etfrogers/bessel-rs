#![feature(float_gamma)]
mod amos;
pub use amos::{BesselError, Scaling, GammaError};
use amos::zbesj;
use num::complex::Complex64;

pub fn bessel_j(order: f64, z: f64) -> Result<Complex64, BesselError> {
    zbesj(z.into(), order, Scaling::Unscaled, 1).map(|v| v.0[0])
}

pub fn add(left: u64, right: u64) -> u64 {
    left + right
}

#[cfg(test)]
mod tests;
