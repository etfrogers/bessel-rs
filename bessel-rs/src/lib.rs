mod amos;
mod bessel_zeros;

pub use amos::{BesselError, GammaError, Scaling};
use amos::{
    complex_airy, complex_airy_b, complex_bessel_h, complex_bessel_i, complex_bessel_j,
    complex_bessel_k, complex_bessel_y,
};
use num::{
    Complex,
    complex::{Complex64, ComplexFloat},
};

use crate::amos::{HankelKind, MACHINE_CONSTANTS};
pub use bessel_zeros::{BesselFunType, bessel_zeros};

pub trait BackTo<T> {
    fn back_to(&self) -> T;
}

impl BackTo<Complex64> for Complex64 {
    fn back_to(&self) -> Complex64 {
        *self
    }
}

impl BackTo<f64> for f64 {
    fn back_to(&self) -> f64 {
        *self
    }
}

impl BackTo<f64> for Complex64 {
    fn back_to(&self) -> f64 {
        assert!(
            self.im() < 1000.0 * MACHINE_CONSTANTS.abs_error_tolerance,
            "Imaginary part of Bessel funtion for real input is too large: {:?}",
            self
        );
        self.re()
    }
}

// TODO work with abritrary bit-depth floats
// TODO bessel derivatives
// TODO negative orders

pub fn bessel_j<ZT, OT>(order: OT, z: ZT) -> Result<ZT, BesselError>
where
    ZT: Into<Complex<f64>>,
    Complex64: BackTo<ZT>,
    OT: Into<f64>,
{
    complex_bessel_j(z.into(), order.into(), Scaling::Unscaled, 1)
        .map(|v| v.0[0])
        .map(|v| v.back_to())
}

pub fn bessel_i<ZT, OT>(order: OT, z: ZT) -> Result<ZT, BesselError>
where
    ZT: Into<Complex<f64>> + BackTo<ZT>,
    Complex64: BackTo<ZT>,
    OT: Into<f64>,
{
    complex_bessel_i(z.into(), order.into(), Scaling::Unscaled, 1)
        .map(|v| v.0[0])
        .map(|v| v.back_to())
}

pub fn bessel_k<ZT, OT>(order: OT, z: ZT) -> Result<ZT, BesselError>
where
    ZT: Into<Complex<f64>> + BackTo<ZT>,
    Complex64: BackTo<ZT>,
    OT: Into<f64>,
{
    complex_bessel_k(z.into(), order.into(), Scaling::Unscaled, 1)
        .map(|v| v.0[0])
        .map(|v| v.back_to())
}

pub fn bessel_y<ZT, OT>(order: OT, z: ZT) -> Result<ZT, BesselError>
where
    ZT: Into<Complex<f64>> + BackTo<ZT>,
    Complex64: BackTo<ZT>,
    OT: Into<f64>,
{
    complex_bessel_y(z.into(), order.into(), Scaling::Unscaled, 1)
        .map(|v| v.0[0])
        .map(|v| v.back_to())
}

pub fn hankel<ZT, OT>(order: OT, z: ZT, kind: HankelKind) -> Result<ZT, BesselError>
where
    ZT: Into<Complex<f64>> + BackTo<ZT>,
    Complex64: BackTo<ZT>,
    OT: Into<f64>,
{
    complex_bessel_h(z.into(), order.into(), Scaling::Unscaled, kind, 1)
        .map(|v| v.0[0])
        .map(|v| v.back_to())
}

pub fn airy<ZT>(z: ZT) -> Result<ZT, BesselError>
where
    ZT: Into<Complex<f64>> + BackTo<ZT>,
    Complex64: BackTo<ZT>,
{
    complex_airy(z.into(), false, Scaling::Unscaled)
        .map(|v| v.0)
        .map(|v| v.back_to())
}

pub fn airyp<ZT>(z: ZT) -> Result<ZT, BesselError>
where
    Complex64: BackTo<ZT>,
    ZT: Into<Complex<f64>> + BackTo<ZT>,
{
    complex_airy(z.into(), true, Scaling::Unscaled)
        .map(|v| v.0)
        .map(|v| v.back_to())
}

pub fn airy_b<ZT>(z: ZT) -> Result<ZT, BesselError>
where
    Complex64: BackTo<ZT>,
    ZT: Into<Complex<f64>> + BackTo<ZT>,
{
    complex_airy_b(z.into(), false, Scaling::Unscaled)
        .map(|v| v)
        .map(|v| v.back_to())
}

pub fn airy_bp<ZT>(z: ZT) -> Result<ZT, BesselError>
where
    Complex64: BackTo<ZT>,
    ZT: Into<Complex<f64>> + BackTo<ZT>,
{
    complex_airy_b(z.into(), true, Scaling::Unscaled).map(|v| v.back_to())
}

#[cfg(test)]
mod tests;
