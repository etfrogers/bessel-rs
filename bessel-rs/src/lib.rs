//! # Bessel-Rs: Bessel functions in pure Rust
//!
//! A crate implementing pure Rust translations of [Amos' complex
//! Bessel function algorthims](https://www.netlib.org/amos/)
//!
//! The development of this crate started as a a reaction to finding that:
//! 1) no pure Rust implementations of Bessel functions existed
//! 2) all the standard library functions (in all languages) seem to use a wrapper
//!    around Amos' original fortran.
//!
//! For example,
//! - [Python - Scipy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.jv.html#scipy.special.jv)
//! - [Julia - OpenSpecFunc](https://github.com/JuliaMath/openspecfun)
//! - [Rust (the best way of calculating Bessel functions I've found)](https://crates.io/crates/complex-bessel-rs/1.2.1)
//!
//! There are other implentations in some cases for integer order and real argument,
//! but the general case is all Amos.
//!
//! The aim of this crate is to translate the Amos Fortran code
//! into idomatic Rust, while retaining full compatibility with Amos' code.
//!
//! The primary test of this crate is that it gives the same values
//! (to within approx 10 significant figures, subject to the considerations
//! below) as the Fortran code for all inputs.
//!
//! ## Note on errors - edited from the original Amos documentation
//!
//! In most complex variable computation, one must evaluate
//! elementary functions. When the magnitude of z or (effective) order is
//! large, losses of significance by argument reduction occur.
//! consequently, if either one exceeds `u1=(0.5/eps).sqrt()`, then
//! losses exceeding half of machine precision are likely and an error flag
//! [BesselError::PartialLossOfSignificance] is triggered where
//! `eps = f64::EPSILON` is double precision unit roundoff.
//! If either z or order is larger than `u2=0.5/eps`, then all significance is
//! lost and [BesselError::LossOfSignificance] is returned.
//! In order to use the int function, arguments
//! must be further restricted not to exceed the largest machine
//! integer, `u3=(i32::MAX as f64) * 0.5`[^note1].
//! Thus, the magnitude of z and (effective) order is
//! restricted by `u2.min(u3)`. With the 64 bit versions given above, u1, u2, and u3
//! are approximately  1.3e+8, 1.8e+16, 2.1e+9.[^note2]
//!
//! [^note1]: This is used for the time being to match the Fortran
//! code, but may be relaxed later.
//!
//! [^note2]: It is intended to extend the current `f64` implementation
//! to allow different size floats in the future.
//!
//! The approximate relative error in the magnitude of a complex
//! bessel function can be expressed by `eps * 10.0.pow(s)` where
//! `eps = f64::EPSILON` is the nominal precision and `10.0.pow(s)`
//! represents the increase in error due to argument reduction in the
//! elementary functions. Here,
//! `s = 1.0.max(z.abs().log10().abs().max(order.log10().abs())`
//! approximately (i.e. s = max(1, abs(exponent of
//! abs(z), abs(exponent of order))). However, the phase angle may
//! have only absolute accuracy. This is most likely to occur when
//! one component (in absolute value) is larger than the other by
//! several orders of magnitude. if one component is `10.0.pow(k)` larger
//! than the other, then one can expect only `p.log10().abs() - k).max(0)`,
//! significant digits; or, stated another way, when k exceeds
//! the exponent of p, no significant digits remain in the smaller
//! component. However, the phase angle retains absolute accuracy
//! because, in complex arithmetic with precision p, the smaller
//! component will not (as a rule) decrease below p times the
//! magnitude of the larger component. In these extreme cases,
//! the principal phase angle is on the order of +p, -p, pi/2-p,
//! or -pi/2+p.

mod amos;
mod bessel_zeros;

pub use amos::{BesselError, GammaError, HankelKind, Scaling};
pub use amos::{
    complex_airy, complex_airy_b, complex_bessel_h, complex_bessel_i, complex_bessel_j,
    complex_bessel_k, complex_bessel_y,
};
use num::{
    Complex,
    complex::{Complex64, ComplexFloat},
};

use crate::amos::MACHINE_CONSTANTS;
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
// TODO Make tolerance margin variable (particularly for partial loss of siginifcance)
// TODO ignore nans in PLOS
// TODO Overflow to pos of negative infinity
// TODO ComplexOutputForRealInput error, rather than panic

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
    complex_airy_b(z.into(), false, Scaling::Unscaled).map(|v| v.back_to())
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
