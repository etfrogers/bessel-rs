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
//! - [This crate](https://docs.rs/amos-bessel-rs/latest/amos_bessel_rs/) - A modern, idiomatic pure-Rust translation of Amos' algorithms.
//!   Features zero-allocation buffer-passing APIs (`_into`), register-resident recurrence loops, and SIMD-friendly Horner polynomial evaluation.
//!   Benchmarks show it consistently matches or outpaces both Fortran AMOS and other translations while offering full `#![no_std]` support.
//! - [Complex Bessel](http://docs.rs/complex-bessel/latest/complex_bessel/) - A line-by-line translation of Amos code with a very good
//!   [comparison tool](https://github.com/elgar328/complex-bessel-test) to confirm both accuracy and computational speed. Optimised
//!   for accuracy and speed using FMA tools to aid the compiler.
//! - [Real Bessel](https://crates.io/crates/real-bessel) - A dedicated pure `core` zero-allocation crate that calculates real-only Bessel
//!   functions *J* and *Y* for integer order. Faster for these simple cases.
//! - [Complex Bessel rs](https://crates.io/crates/complex-bessel-rs/) - A wrapper around the Amos' Fortran functions with a Rust API.
//!   Good if you want guarantees that answers will be identical to Fortran, but requires a Fortran compiler in your toolchain to compile.
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
//! - Sequence functions are available in both allocating forms ([`amos::complex_bessel_j`], etc.) returning `Result<(Vec<Complex<T>>, SequenceInfo), BesselError<T>>`
//!   and zero-allocation slice-filling forms ([`amos::complex_bessel_j_into`], etc.) returning `Result<SequenceInfo, BesselError<T>>`.
//!
//! - The returned [`SequenceInfo`] provides metadata about the computation:
//!   - `n_zeros`: the number of elements set to zero due to underflow. Underflow zeroes occur at the **end** of the sequence (highest orders)
//!     for $J_\nu$ and $I_\nu$, and at the **start** of the sequence (lowest orders) for $Y_\nu$, $K_\nu$, and $H_\nu^{(m)}$.
//!   - `partial_loss_of_significance`: `true` if extreme values of $|z|$ or `order` caused argument reduction to lose more than half
//!     of machine precision. The computed values are still returned as `Ok` because the algorithm converged.
//!
//! - Errors ([`BesselError`]) implement `Copy` with zero heap allocation and are reserved strictly for true calculation failures
//!   (such as overflow, non-convergence, or complete loss of significance).
//!
//! ### Derivatives
//!
//! Derivatives of Bessel and Hankel functions with respect to the argument $z$ are provided in the
//! [`derivatives`] module. For each function family, both the first derivative (e.g. [`derivatives::bessel_j_p`])
//! and arbitrary $k$-th order derivatives (e.g. [`derivatives::bessel_j_derivative`]) are available:
//!
//! ```rust
//! use amos_bessel_rs::{Scaling, derivatives::{bessel_j_p, bessel_j_derivative}};
//!
//! // First derivative J_0'(1.0):
//! let dj = bessel_j_p(0.0, 1.0).unwrap();
//!
//! // Second derivative (d/dz)^2 J_0(1.0):
//! let d2j = bessel_j_derivative(0.0, 1.0, 2, Scaling::Unscaled).unwrap();
//! ```
//!
//! ## Note on accuracy
//!
//! When the magnitude of `z` or (effective) `order` is extremely large, losses of significance
//! by argument reduction occur in the underlying computations.
//!
//! If either one exceeds `u1 = (0.5/eps).sqrt()` (approx `1.3e8` for `f64`), losses exceeding half
//! of machine precision are likely and `SequenceInfo::partial_loss_of_significance` is set to `true`.
//! If either `z` or `order` is larger than `u2 = 0.5/eps` (approx `1.8e16` for `f64`), then all
//! significance is lost and [BesselError::LossOfSignificance] is returned.
//!
//! For a full mathematical breakdown of the relative error and phase angle accuracy
//! based on the original Amos documentation, please see the
//! [Performance & Accuracy Guide](https://etfrogers.github.io/bessel-rs/).
//!
//! ## Branch cuts and signed zero (`-0.0`)
//!
//! The functions $Y_\nu$, $K_\nu$, $H_\nu^{(1)}$, and $H_\nu^{(2)}$ possess a
//! branch cut along the negative real axis $(-\infty, 0]$.
//!
//! ### Phase Convention
//! Following standard DLMF conventions (DLMF 10.11 and 10.25), the principal branch is defined by:
//!
//! $$-\pi < \arg(z) \le \pi$$
//!
//! The branch cut adheres to the upper half-plane: values on the negative real axis have phase $\arg(z) = +\pi$.
//!
//! ### Signed Zero in Floating-Point Arithmetic
//! In IEEE 754 arithmetic, zero has a distinct sign (`+0.0` vs `-0.0`):
//! - By standard complex analysis conventions, `Complex::new(-x, 0.0)` corresponds to $\arg(z) = +\pi$ (upper edge of the cut).
//! - Conversely, `Complex::new(-x, -0.0)` represents an approach from the lower half-plane ($\arg(z) = -\pi$).
//!
//! In the underlying Amos algorithm, checks are formulated using `z.im < 0.0`. Under IEEE 754 floating-point
//! rules, `-0.0 < 0.0` evaluates to `false`. Consequently, both `Complex::new(-x, 0.0)` and `Complex::new(-x, -0.0)`
//! evaluate consistently on the **upper edge** of the cut ($\arg(z) = +\pi$).
//!
//! If values along the lower edge ($\arg(z) \to -\pi$) are required, the standard DLMF cross-cut
//! continuation relations (DLMF 10.11 and 10.34) should be applied:
//! - $Y_\nu(z e^{-i\pi}) = e^{i\nu\pi} Y_\nu(z) - 2i\cos(\nu\pi) J_\nu(z)$
//! - $K_\nu(z e^{-i\pi}) = e^{i\nu\pi} K_\nu(z) + i\pi I_\nu(z)$
//! - $H_\nu^{(1)}(z e^{-i\pi}) = 2\cos(\nu\pi) H_\nu^{(1)}(z) + e^{-i\nu\pi} H_\nu^{(2)}(z)$
//!
#![cfg_attr(not(feature = "std"), no_std)]

#[cfg(feature = "alloc")]
#[macro_use]
extern crate alloc;

use num::Complex;

/// Container for the complex_\[func\] version of the Bessel and Airy functions
/// for finer control of the calculation and results
pub mod amos;

/// Functions for computing derivatives of Bessel functions with respect to the argument $z$.
pub mod derivatives;
pub(crate) mod reflections;
mod types;

pub use amos::{HankelKind, Scaling};
use amos::{
    complex_airy, complex_airy_b, complex_bessel_i_into, complex_bessel_j_into,
    complex_bessel_k_into, complex_bessel_y_into, complex_hankel1_into, complex_hankel2_into,
};
use types::simple_bessel_wrapper;
pub use types::{BesselError, BesselFloat, BesselInput, SequenceInfo};

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
        .and_then(ZT::back_from)
}

/// Computes the derivative of the Airy function Ai'(z).
pub fn airyp<FT: BesselFloat, ZT: BesselInput<FT>>(z: ZT) -> Result<ZT, BesselError<FT>> {
    complex_airy(z.into(), true, Scaling::Unscaled)
        .map(|x| x.0)
        .and_then(ZT::back_from)
}

/// Computes the Airy function of the second kind Bi(z).
pub fn airy_b<FT: BesselFloat, ZT: BesselInput<FT>>(z: ZT) -> Result<ZT, BesselError<FT>> {
    complex_airy_b(z.into(), false, Scaling::Unscaled)
        .map(|x| x.0)
        .and_then(ZT::back_from)
}

/// Computes the derivative of the Airy function of the second kind Bi'(z).
pub fn airy_bp<FT: BesselFloat, ZT: BesselInput<FT>>(z: ZT) -> Result<ZT, BesselError<FT>> {
    complex_airy_b(z.into(), true, Scaling::Unscaled)
        .map(|x| x.0)
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
