#![warn(missing_docs, clippy::all)]
//! # amos-bessel-rs: Bessel functions in pure Rust
//!
//! A crate implementing pure Rust translations of [Amos' complex
//! Bessel function algorithms](https://www.netlib.org/amos/)
//!
//! The aim of this crate is to translate the Amos Fortran code
//! into idiomatic Rust, while retaining full compatibility with Amos' code.
//!
//! ## Alternatives
//!
//! To calculate Bessel functions in Rust there are now several alternatives
//!
//! - [Complex Bessel rs](https://crates.io/crates/complex-bessel-rs/) - A wrapper around the Amos' Fortran functions with a Rust API.
//!   Good if you want guarantees that answers will be the same as Fortran, but requires a Fortran compiler in your toolchain to compile.
//! - [Complex Bessel](http://docs.rs/complex-bessel/latest/complex_bessel/) - A line-by-line translation of Amos code with a very good
//!   [comparison tool](https://github.com/elgar328/complex-bessel-test) to confirm both accuracy and computational speed. Carefully optimised
//!   for accuracy and speed using detailed tools (e.g. implementation of FMA) to aid the compiler.
//! - [This crate](https://docs.rs/amos-bessel-rs/latest/amos_bessel_rs/) - A more idiomatic translation of the Fortran code: using Rust
//!   tools. Relies on the compiler to optimise as best it can. A fork of the elgar328's [comparison tool](https://github.com/etfrogers/complex-bessel-test) shows similar accuracy and
//!   execution speed.
//! - **Real Bessel** - WIP (soon to be released) crate that calculates real-only Bessel function's *J*, and *Y* for integer order. *J* takes
//!   real inputs, *Y* is restricted to positive inputs (to give real answers). This implementation is faster for these simple cases.
//!
//! The primary test of this crate is that it gives the same values
//! (to within approx 10 significant figures, subject to the considerations
//! below) as the Fortran code for all inputs.
//!
//! ## Usage
//!
//! The primary interface to this crate are the functions in the crate root:
//! [bessel_j], [bessel_i], etc.
//! These functions are implemented for any real `order` (including negative orders)
//! and real or complex argument `z`. All return a single value of the function at the
//! given order and argument. If the argument is real (f64) then the function will attempt to
//! return a real (f64) output. If the calculated answer is complex, then a [BesselError::ComplexOutputForRealInput]
//! error will be returned, which contains the complex value, in case it is needed.
//! If the argument is complex, the output will be complex.
//!
//! Other errors returned by these functions are on overflow, failure to converge, or
//! loss of significance in the answer (triggered by extreme values of inputs). See [BesselError] for more details.
//!
//! These functions implement calculation for negative orders (not in Amos' original functions)
//! using the reflection formulae from [DLMF](https://dlmf.nist.gov/10)
//!
//! ### Amos interface
//!
//! Amos' original functions are available in the [amos] module, and are called complex_\[func\].
//! These functions expose additional functionality, but at some loss of simplicity.
//!
//! #### Limitations
//!
//! - They do not implement the reflection formulae, and so are limited to positive orders (returning a [BesselError::InvalidInput]
//!   error if a negative order is requested).
//!
//! #### Inputs
//!
//! - The Amos functions take an additional `scaling` parameter, of type [Scaling], which, if set to [Scaling::Scaled],
//!   returns the scaled Bessel function value. That is, the return value is the function value multiplied by a
//!   (positive or negative) exponential factor to remove the exponential growth or decay that occurs as the argument is goes to infinity.
//!   The precise formula for the scaled function value is given in documentation for each function. The simple functions
//!   always return the unscaled function value.
//!
//! - The Amos functions take an additional argument `n` which specifies the number of orders to return. The return value is a
//!   vector of `Complex<f64>` containing the values of the function at orders `[order, order + 1, ..., order + n]`.
//!   This is implemented because, if multiple orders are required, the computation of them in a single run of the algroithm
//!   is more efficient than running it multiple times.
//!
//! #### Return values
//!
//! - They return an additional error: [BesselError::PartialLossOfSignificance], in cases where the algorithm
//!   as converged, but the result is not as accurate as normal due to loss of significance. It occurs on
//!   extreme values of inputs, and is a feature of the Amos algorithm. It is hidden from the user in the simpler functions,
//!   so that the user does not need to worry about it: if the error is returned by the underlying Amos function, then
//!   it is unwrapped and returned as `Ok(value)`. by the simpler functions.
//!
//! - The general form of the Amos functions return is a `Result<(Vec<Complex<f64>>, usize>, BesselError>`, where the `Vec` contains
//!   the values of the function at orders `[order, order + 1, ..., order + n]` and `nz` contains the number of the elements
//!   in the Vec that have been set to zero due to underflow.
//!
//!
//! ## Note on accuracy - edited from the original Amos documentation
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

/// Container for the complex_\[func\] version of the Bessel and Airy functions
/// for finer control of the calculation and results
pub mod amos;

mod reflections;
mod types;

use std::ops::Mul;

use crate::{
    amos::{
        complex_airy, complex_airy_b, complex_bessel_i, complex_bessel_j, complex_bessel_k,
        complex_bessel_y, complex_hankel1, complex_hankel2,
    },
    reflections::{
        as_integer, integer_sign, reflect_h_element, reflect_i_element, reflect_j_element,
        reflect_y_element,
    },
};
pub use amos::{HankelKind, Scaling};

use num::{Complex, complex::Complex64};
pub use types::{BackFrom, BesselError};
use types::{BesselResult, simple_bessel_wrapper};

// TODO work with abritrary bit-depth floats
// TODO bessel derivatives
// TODO Overflow to positive or negative infinity, or zero?
// TODO Handle Vectors/ndarrays for z

/// A trait for types that can be used as input to Bessel functions.
///
/// This trait is implemented for `f64` and `Complex<f64>`, allowing the Bessel functions
/// to accept both real and complex arguments.
pub trait BesselInput:
    Into<Complex<f64>>
    + BackFrom<Complex<f64>>
    + Mul<f64, Output = Self>
    + BackFrom<BesselResult<Complex<f64>>>
{
}

impl BesselInput for f64 {}
impl BesselInput for Complex<f64> {}

/// Computes the Bessel function of the first kind Jv(z).
///
/// # Arguments
/// * `order` - The order of the Bessel function (can be non-integer).
/// * `z` - The complex or real argument.
pub fn bessel_j<ZT: BesselInput, OT: Into<f64>>(order: OT, z: ZT) -> Result<ZT, BesselError> {
    let order: f64 = order.into();
    let z: Complex<f64> = z.into();
    if order >= 0.0 {
        return ZT::back_from(&bessel_j_single(order, z));
    }

    // Special case for negative integer order: J(-n, z) = (-1)^n J(n, z)
    let abs_order: f64 = order.abs();
    if let Some(n) = as_integer(abs_order) {
        return ZT::back_from(&bessel_j_single(abs_order, z)).map(|x| x * integer_sign(n));
    }

    // General case: need both J and Y at positive |ν|
    let j = bessel_j_single(abs_order, z)?;
    let y = bessel_y_single(abs_order, z)?;

    ZT::back_from(&reflect_j_element(abs_order, j, y))
}

/// Computes the modified Bessel function of the first kind Iv(z).
///
/// # Arguments
/// * `order` - The order of the Bessel function.
/// * `z` - The complex or real argument.
pub fn bessel_i<ZT: BesselInput, OT: Into<f64>>(order: OT, z: ZT) -> Result<ZT, BesselError> {
    let order = order.into();
    let z = z.into();
    if order >= 0.0 {
        return ZT::back_from(&bessel_i_single(order, z));
    }

    let abs_order = order.abs();

    // Integer shortcut: I_{-n}(z) = I_n(z)
    if as_integer(abs_order).is_some() {
        return ZT::back_from(&bessel_i_single(abs_order, z));
    }

    // General case: need both I and K at positive |ν|
    let i = bessel_i_single(abs_order, z)?;
    let k = bessel_k_single(abs_order, z)?;

    ZT::back_from(&reflect_i_element(abs_order, i, k))
}

/// Computes the modified Bessel function of the second kind Kv(z).
///
/// # Arguments
/// * `order` - The order of the Bessel function.
/// * `z` - The complex or real argument.
pub fn bessel_k<ZT: BesselInput, OT: Into<f64>>(order: OT, z: ZT) -> Result<ZT, BesselError> {
    // K_{-ν}(z) = K_ν(z) (DLMF 10.27.3)
    let abs_order = order.into().abs();
    ZT::back_from(&bessel_k_single(abs_order, z.into()))
}

/// Computes the Bessel function of the second kind Yv(z).
///
/// # Arguments
/// * `order` - The order of the Bessel function.
/// * `z` - The complex or real argument.
pub fn bessel_y<ZT: BesselInput, OT: Into<f64>>(order: OT, z: ZT) -> Result<ZT, BesselError> {
    let order = order.into();
    let z = z.into();
    if order >= 0.0 {
        return ZT::back_from(&bessel_y_single(order, z));
    }

    let abs_order = order.abs();

    // Integer shortcut: Y_{-n}(z) = (-1)^n * Y_n(z)
    if let Some(n) = as_integer(abs_order) {
        let y = bessel_y_single(abs_order, z)?;
        return ZT::back_from(&(y * integer_sign(n)));
    }

    // General case: need both J and Y at positive |ν|
    let j = bessel_j_single(abs_order, z)?;
    let y = bessel_y_single(abs_order, z)?;

    ZT::back_from(&reflect_y_element(abs_order, j, y))
}

/// Computes the Hankel function Hv(z) of the first or second kind.
///
/// # Arguments
/// * `order` - The order of the Hankel function.
/// * `z` - The complex or real argument.
/// * `kind` - The kind of Hankel function (First or Second).
pub fn hankel<ZT: BesselInput, OT: Into<f64>>(
    order: OT,
    z: ZT,
    kind: HankelKind,
) -> Result<ZT, BesselError> {
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
    ZT::back_from(&h)
}

/// Computes the Airy function Ai(z).
pub fn airy<ZT: BesselInput>(z: ZT) -> Result<ZT, BesselError> {
    ZT::back_from(&complex_airy(z.into(), false, Scaling::Unscaled).map(|x| x.0))
}

/// Computes the derivative of the Airy function Ai'(z).
pub fn airyp<ZT: BesselInput>(z: ZT) -> Result<ZT, BesselError> {
    ZT::back_from(&complex_airy(z.into(), true, Scaling::Unscaled).map(|x| x.0))
}

/// Computes the Airy function of the second kind Bi(z).
pub fn airy_b<ZT: BesselInput>(z: ZT) -> Result<ZT, BesselError> {
    ZT::back_from(&complex_airy_b(z.into(), false, Scaling::Unscaled))
}

/// Computes the derivative of the Airy function of the second kind Bi'(z).
pub fn airy_bp<ZT: BesselInput>(z: ZT) -> Result<ZT, BesselError> {
    ZT::back_from(&complex_airy_b(z.into(), true, Scaling::Unscaled))
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

#[cfg(any(test, feature = "test-utils"))]
#[doc(hidden)]
pub mod test_utils;
