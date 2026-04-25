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
//! TODO other options
//! TODO negative approach
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
mod reflections;
mod types;

use std::ops::Mul;

pub use amos::{GammaError, HankelKind, Scaling};
pub use amos::{
    complex_airy, complex_airy_b, complex_bessel_h, complex_bessel_i, complex_bessel_j,
    complex_bessel_k, complex_bessel_y, complex_hankel1, complex_hankel2,
};
use num::{Complex, complex::Complex64};
pub use types::{BackTo, BesselError};

pub use bessel_zeros::{BesselFunType, bessel_zeros};

use crate::reflections::{
    as_integer, integer_sign, reflect_h_element, reflect_i_element, reflect_j_element,
    reflect_y_element,
};

// TODO work with abritrary bit-depth floats
// TODO bessel derivatives
// TODO Make tolerance margin variable (particularly for partial loss of siginifcance) and tighten where possible
// TODO ignore nans in PLOS
// TODO Overflow to positive or negative infinity, or zero?
// TODO Handle Vectors/ndarrays for z
// TODO Unwrap PLOS?
// TODO smaller crates for real, zeros
// TODO check validity of real functions

pub trait BesselInput: Into<Complex<f64>> + BackTo<Self> + Mul<f64, Output = Self>
where
    Complex<f64>: BackTo<Self>,
{
}

impl BesselInput for f64 {}
impl BesselInput for Complex<f64> {}

pub fn bessel_j<ZT, OT>(order: OT, z: ZT) -> Result<ZT, BesselError>
where
    ZT: BesselInput,
    Complex64: BackTo<ZT>,
    OT: Into<f64>,
{
    let order: f64 = order.into();
    let z: Complex<f64> = z.into();
    if order >= 0.0 {
        return bessel_j_single(order, z)?.back_to();
    }

    // Special case for negative integer order: J(-n, z) = (-1)^n J(n, z)
    let abs_order: f64 = order.abs();
    if let Some(n) = as_integer(abs_order) {
        return Ok(bessel_j_single(abs_order, z)?.back_to()? * integer_sign(n));
    }

    // General case: need both J and Y at positive |ν|
    let j = bessel_j_single(abs_order, z)?;
    let y = bessel_y_single(abs_order, z)?;

    Ok(reflect_j_element(abs_order, j.into(), y.into()).back_to()?)
}

pub fn bessel_i<ZT, OT>(order: OT, z: ZT) -> Result<ZT, BesselError>
where
    ZT: Into<Complex<f64>> + BackTo<ZT>,
    Complex<f64>: BackTo<ZT>,
    OT: Into<f64>,
{
    let order = order.into();
    let z = z.into();
    if order >= 0.0 {
        return bessel_i_single(order, z)?.back_to();
    }

    let abs_order = order.abs();

    // Integer shortcut: I_{-n}(z) = I_n(z)
    if as_integer(abs_order).is_some() {
        return bessel_i_single(abs_order, z)?.back_to();
    }

    // General case: need both I and K at positive |ν|
    let i = bessel_i_single(abs_order, z)?;
    let k = bessel_k_single(abs_order, z)?;

    Ok(reflect_i_element(abs_order, i, k).back_to()?)
}

pub fn bessel_k<ZT, OT>(order: OT, z: ZT) -> Result<ZT, BesselError>
where
    ZT: Into<Complex<f64>> + BackTo<ZT>,
    Complex64: BackTo<ZT>,
    OT: Into<f64>,
{
    // K_{-ν}(z) = K_ν(z) (DLMF 10.27.3)
    let abs_order = order.into().abs();
    bessel_k_single(abs_order, z.into())?.back_to()
}

pub fn bessel_y<ZT, OT>(order: OT, z: ZT) -> Result<ZT, BesselError>
where
    ZT: Into<Complex<f64>> + BackTo<ZT>,
    Complex64: BackTo<ZT>,
    OT: Into<f64>,
{
    let order = order.into();
    let z = z.into();
    if order >= 0.0 {
        return bessel_y_single(order, z)?.back_to();
    }

    let abs_order = order.abs();

    // Integer shortcut: Y_{-n}(z) = (-1)^n * Y_n(z)
    if let Some(n) = as_integer(abs_order) {
        let y = bessel_y_single(abs_order, z)?;
        return Ok((y * integer_sign(n)).back_to()?);
    }

    // General case: need both J and Y at positive |ν|
    let j = bessel_j_single(abs_order, z)?;
    let y = bessel_y_single(abs_order, z)?;

    Ok(reflect_y_element(abs_order, j, y).back_to()?)
}

pub fn hankel<ZT, OT>(order: OT, z: ZT, kind: HankelKind) -> Result<ZT, BesselError>
where
    ZT: Into<Complex<f64>> + BackTo<ZT>,
    Complex64: BackTo<ZT>,
    OT: Into<f64>,
{
    let order = order.into();
    let abs_order = order.abs();
    let mut h = match kind {
        HankelKind::First => hankel1_single(abs_order, z.into())?,
        HankelKind::Second => hankel2_single(abs_order, z.into())?,
    };

    if order < 0.0 {
        // Need to reflect the Hankel function for negative orders, but this is just a rotation, so no loss of significance.
        h = reflect_h_element(abs_order, kind, h);
    }
    Ok(h.back_to()?)
}

pub fn airy<ZT>(z: ZT) -> Result<ZT, BesselError>
where
    ZT: Into<Complex<f64>> + BackTo<ZT>,
    Complex64: BackTo<ZT>,
{
    complex_airy(z.into(), false, Scaling::Unscaled)?
        .0
        .back_to()
}

pub fn airyp<ZT>(z: ZT) -> Result<ZT, BesselError>
where
    Complex64: BackTo<ZT>,
    ZT: Into<Complex<f64>> + BackTo<ZT>,
{
    complex_airy(z.into(), true, Scaling::Unscaled)?.0.back_to()
}

pub fn airy_b<ZT>(z: ZT) -> Result<ZT, BesselError>
where
    Complex64: BackTo<ZT>,
    ZT: Into<Complex<f64>> + BackTo<ZT>,
{
    complex_airy_b(z.into(), false, Scaling::Unscaled)?.back_to()
}

pub fn airy_bp<ZT>(z: ZT) -> Result<ZT, BesselError>
where
    Complex64: BackTo<ZT>,
    ZT: Into<Complex<f64>> + BackTo<ZT>,
{
    complex_airy_b(z.into(), true, Scaling::Unscaled)?.back_to()
}

use paste::paste;
simple_bessel_wrapper!(bessel_j);
simple_bessel_wrapper!(bessel_y);
simple_bessel_wrapper!(bessel_i);
simple_bessel_wrapper!(bessel_k);
simple_bessel_wrapper!(hankel1);
simple_bessel_wrapper!(hankel2);

#[cfg(test)]
mod tests;
