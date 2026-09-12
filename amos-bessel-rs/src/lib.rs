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
//! and real or complex argument `z` (for both `f64` and `f32` via [BesselFloat]). All return a single
//! value of the function at the given order and argument. If the argument is real then the function
//! will attempt to return a real output. If the calculated answer is complex, then a [BesselError::ComplexOutputForRealInput]
//! error will be returned, which contains the complex value, in case it is needed.
//! If the argument is complex, the output will be complex.
//!
//! Other errors returned by these functions are on overflow, failure to converge, or
//! loss of significance in the answer (triggered by extreme values of inputs). See [BesselError] for more details.
//!
//! These functions implement calculation for negative orders using the reflection formulae
//! from [DLMF](https://dlmf.nist.gov/10).
//!
//! ### Amos interface
//!
//! Amos' core functions are available in the [amos] module, named `complex_[func]`.
//! These functions expose additional functionality, but at some loss of simplicity:
//!
//! #### Negative orders
//!
//! - Like the root functions, the [amos] functions natively support negative orders across all function
//!   families using the reflection identities from DLMF.
//!
//! #### Inputs
//!
//! - The Amos functions take an additional `scaling` parameter, of type [Scaling], which, if set to [Scaling::Scaled],
//!   returns the scaled Bessel function value. That is, the return value is the function value multiplied by a
//!   (positive or negative) exponential factor to remove the exponential growth or decay that occurs as the argument goes to infinity.
//!   The precise formula for the scaled function value is given in documentation for each function. The simple functions
//!   always return the unscaled function value.
//!
//! - The Amos functions take an additional argument `n` which specifies the number of orders to return. The return value is a
//!   vector of `Complex<T>` containing the values of the function at orders `[order, order + 1, ..., order + n - 1]`.
//!   This is implemented because, if multiple orders are required, computing them in a single run of the algorithm
//!   is much more efficient than running it repeatedly.
//!
//! #### Return values
//!
//! - They return an additional error variant: [BesselError::PartialLossOfSignificance], in cases where the algorithm
//!   has converged, but the result is not as accurate as normal due to loss of significance. It occurs on
//!   extreme values of inputs, and is a feature of the Amos algorithm. It is hidden from the user in the simpler functions,
//!   so that the user does not need to worry about it: if the error is returned by the underlying Amos function, then
//!   it is unwrapped and returned as `Ok(value)` by the simpler functions.
//!
//! - The general form of the Amos functions return is a `Result<(Vec<Complex<T>>, usize), BesselError<T>>`, where the `Vec` contains
//!   the values of the function at orders `[order, order + 1, ..., order + n - 1]` and `n_zeros` contains the number of elements
//!   in the `Vec` that have been set to zero due to underflow.
//!
//! ### Derivatives
//!
//! Derivatives of Bessel and Hankel functions with respect to the argument $z$ are provided in the
//! [`derivatives`] module. For each function family, both the first derivative (e.g. [`derivatives::bessel_j_p`])
//! and arbitrary $k$-th order derivatives (e.g. [`derivatives::bessel_j_derivative`]) are available:
//!
//! ```rust
//! use amos_bessel_rs::derivatives::{bessel_j_p, bessel_j_derivative};
//!
//! // First derivative J_0'(1.0):
//! let dj = bessel_j_p(0.0, 1.0).unwrap();
//!
//! // Second derivative (d/dz)^2 J_0(1.0):
//! let d2j = bessel_j_derivative(0.0, 1.0, 2).unwrap();
//! ```
//!
//! ## Note on accuracy
//!
//! When the magnitude of `z` or (effective) `order` is extremely large, losses of significance
//! by argument reduction occur in the underlying computations.
//!
//! If either one exceeds `u1 = (0.5/eps).sqrt()` (approx `1.3e8` for `f64`), losses exceeding half
//! of machine precision are likely and [BesselError::PartialLossOfSignificance] is triggered.
//! If either `z` or `order` is larger than `u2 = 0.5/eps` (approx `1.8e16` for `f64`), then all
//! significance is lost and [BesselError::LossOfSignificance] is returned.
//!
//! For a full mathematical breakdown of the relative error and phase angle accuracy
//! based on the original Amos documentation, please see the
//! [Performance & Accuracy Guide](https://etfrogers.github.io/bessel-rs/).
#![cfg_attr(not(feature = "std"), no_std)]

#[macro_use]
extern crate alloc;

use num::Complex;

pub(crate) mod prelude {
    pub use alloc::borrow::ToOwned;
    pub use alloc::string::{String, ToString};
    pub use alloc::vec::Vec;
}

/// Container for the complex_\[func\] version of the Bessel and Airy functions
/// for finer control of the calculation and results
pub mod amos;

/// Functions for computing derivatives of Bessel functions with respect to the argument $z$.
pub mod derivatives;
pub(crate) mod reflections;
mod types;

pub use amos::{HankelKind, Scaling};
use amos::{
    complex_airy, complex_airy_b, complex_bessel_i, complex_bessel_j, complex_bessel_k,
    complex_bessel_y, complex_hankel1, complex_hankel2,
};
use types::{AllowPlos, simple_bessel_wrapper};
pub use types::{BesselError, BesselFloat, BesselInput};

// TODO Overflow to positive or negative infinity, or zero?

/// Computes the Bessel function of the first kind Jv(z).
///
/// # Arguments
/// * `order` - The order of the Bessel function (can be non-integer).
/// * `z` - The complex or real argument.
pub fn bessel_j<FT: BesselFloat, ZT: BesselInput<FT>, OT: Into<FT>>(
    order: OT,
    z: ZT,
) -> Result<ZT, BesselError<FT>> {
    bessel_j_single(order.into(), z.into()).and_then(ZT::back_from)
}

/// Computes the modified Bessel function of the first kind Iv(z).
///
/// # Arguments
/// * `order` - The order of the Bessel function.
/// * `z` - The complex or real argument.
pub fn bessel_i<FT: BesselFloat, ZT: BesselInput<FT>, OT: Into<FT>>(
    order: OT,
    z: ZT,
) -> Result<ZT, BesselError<FT>> {
    bessel_i_single(order.into(), z.into()).and_then(ZT::back_from)
}

/// Computes the modified Bessel function of the second kind Kv(z).
///
/// # Arguments
/// * `order` - The order of the Bessel function.
/// * `z` - The complex or real argument.
pub fn bessel_k<FT: BesselFloat, ZT: BesselInput<FT>, OT: Into<FT>>(
    order: OT,
    z: ZT,
) -> Result<ZT, BesselError<FT>> {
    bessel_k_single(order.into(), z.into()).and_then(ZT::back_from)
}

/// Computes the Bessel function of the second kind Yv(z).
///
/// # Arguments
/// * `order` - The order of the Bessel function.
/// * `z` - The complex or real argument.
pub fn bessel_y<FT: BesselFloat, ZT: BesselInput<FT>, OT: Into<FT>>(
    order: OT,
    z: ZT,
) -> Result<ZT, BesselError<FT>> {
    bessel_y_single(order.into(), z.into()).and_then(ZT::back_from)
}

/// Computes the Hankel function Hv(z) of the first or second kind.
///
/// # Arguments
/// * `order` - The order of the Hankel function.
/// * `z` - The complex or real argument.
/// * `kind` - The kind of Hankel function (First or Second).
pub fn hankel<FT: BesselFloat, ZT: BesselInput<FT>, OT: Into<FT>>(
    order: OT,
    z: ZT,
    kind: HankelKind,
) -> Result<ZT, BesselError<FT>> {
    match kind {
        HankelKind::First => hankel1_single(order.into(), z.into()),
        HankelKind::Second => hankel2_single(order.into(), z.into()),
    }
    .and_then(ZT::back_from)
}

/// Computes the Airy function Ai(z).
pub fn airy<FT: BesselFloat, ZT: BesselInput<FT>>(z: ZT) -> Result<ZT, BesselError<FT>> {
    complex_airy(z.into(), false, Scaling::Unscaled)
        .map(|x| x.0)
        .allow_plos()
        .and_then(ZT::back_from)
}

/// Computes the derivative of the Airy function Ai'(z).
pub fn airyp<FT: BesselFloat, ZT: BesselInput<FT>>(z: ZT) -> Result<ZT, BesselError<FT>> {
    complex_airy(z.into(), true, Scaling::Unscaled)
        .map(|x| x.0)
        .allow_plos()
        .and_then(ZT::back_from)
}

/// Computes the Airy function of the second kind Bi(z).
pub fn airy_b<FT: BesselFloat, ZT: BesselInput<FT>>(z: ZT) -> Result<ZT, BesselError<FT>> {
    complex_airy_b(z.into(), false, Scaling::Unscaled)
        .allow_plos()
        .and_then(ZT::back_from)
}

/// Computes the derivative of the Airy function of the second kind Bi'(z).
pub fn airy_bp<FT: BesselFloat, ZT: BesselInput<FT>>(z: ZT) -> Result<ZT, BesselError<FT>> {
    complex_airy_b(z.into(), true, Scaling::Unscaled)
        .allow_plos()
        .and_then(ZT::back_from)
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
