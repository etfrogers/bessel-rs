//! # Bessel Zeros
//!
//! A crate for finding the zeros of the Bessel functions.
//!
//! This crate uses the routine described in:
//! "An Algorithm with ALGOL 60 Program for the Computation of the
//! zeros of the Ordinary Bessel Functions and those of their
//! Derivatives".
//! N. M. Temme
//! Journal of Computational Physics, 32, 270-279 (1979)
//!
//! This crate is inspired by (via several intermediate steps)
//! Adam Wyatt's Matlab version of this code.
//!
//! ## Backends
//!
//! Two backends are available:
//!
//! | Backend | Module / functions | Order type | Speed |
//! |---|---|---|---|
//! | AMOS (`amos-bessel-rs`) | crate root | any `f64`-compatible | slower |
//! | real-bessel | [`fast`] sub-module | `i32` only | faster |
//!
//! ## Primary API (AMOS backend — any real order)
//!
//! - [`bessel_zeros_j`]
//! - [`bessel_zeros_y`]
//! - [`bessel_zeros_jp`]
//! - [`bessel_zeros_yp`]
//! - [`bessel_zeros`] — lower-level, accepts a [`BesselFunType`] and custom precision
//!
//! ## Fast API (real-bessel backend — integer orders only)
//!
//! Use the [`fast`] sub-module for the same functions with integer orders:
//!
//! ```rust
//! use bessel_zeros::fast;
//!
//! let zeros = fast::bessel_zeros_j(0, 5);
//! assert_eq!(zeros.len(), 5);
//! assert!((zeros[0] - 2.40482555769577).abs() < 1e-10);
//! ```
//!
//! ## Example (AMOS backend, non-integer order)
//!
//! ```rust
//! use bessel_zeros::bessel_zeros_j;
//!
//! let zeros = bessel_zeros_j(0.5_f64, 5);
//! assert_eq!(zeros.len(), 5);
//! ```

#![warn(missing_docs)]

pub(crate) mod backend;
pub(crate) mod algorithm;

use algorithm::bessel_zeros_impl;
use backend::amos::AmosBackend;

const DEFAULT_PRECISION: f64 = 1e-14;

/// The type of Bessel function of which to find the zeros.
#[derive(Debug, PartialEq)]
pub enum BesselFunType {
    /// Bessel function of the first kind, Jv(x)
    J,
    /// Bessel function of the second kind, Yv(x)
    Y,
    /// Derivative of the Bessel function of the first kind, Jv'(x)
    JP,
    /// Derivative of the Bessel function of the second kind, Yv'(x)
    YP,
}

impl BesselFunType {
    pub(crate) fn is_non_derivative(&self) -> bool {
        *self == BesselFunType::J || *self == BesselFunType::Y
    }
}

/// Finds the first `n_zeros` zeros of the Bessel function of the first kind Jv(x).
///
/// Uses the AMOS backend; `order` can be any type that converts to `f64`
/// (e.g. `f64`, `f32`, `i32`, `u32`).
pub fn bessel_zeros_j<OT: Into<f64>>(order: OT, n_zeros: usize) -> Vec<f64> {
    bessel_zeros(&BesselFunType::J, order, n_zeros, DEFAULT_PRECISION)
}

/// Finds the first `n_zeros` zeros of the Bessel function of the second kind Yv(x).
///
/// Uses the AMOS backend; `order` can be any type that converts to `f64`.
pub fn bessel_zeros_y<OT: Into<f64>>(order: OT, n_zeros: usize) -> Vec<f64> {
    bessel_zeros(&BesselFunType::Y, order, n_zeros, DEFAULT_PRECISION)
}

/// Finds the first `n_zeros` zeros of the derivative of Jv(x).
///
/// Uses the AMOS backend; `order` can be any type that converts to `f64`.
pub fn bessel_zeros_jp<OT: Into<f64>>(order: OT, n_zeros: usize) -> Vec<f64> {
    bessel_zeros(&BesselFunType::JP, order, n_zeros, DEFAULT_PRECISION)
}

/// Finds the first `n_zeros` zeros of the derivative of Yv(x).
///
/// Uses the AMOS backend; `order` can be any type that converts to `f64`.
pub fn bessel_zeros_yp<OT: Into<f64>>(order: OT, n_zeros: usize) -> Vec<f64> {
    bessel_zeros(&BesselFunType::YP, order, n_zeros, DEFAULT_PRECISION)
}

/// Finds the first `n_zeros` zeros of the specified Bessel function.
///
/// Uses the AMOS backend, so `order` may be any real value (not just integer).
///
/// # Arguments
///
/// * `func_type` - Which Bessel function (J, Y, J', Y').
/// * `order`     - The order (must be ≥ 0); any type convertible to `f64`.
/// * `n_zeros`   - Number of zeros to return.
/// * `precision` - Relative error tolerance on each root.
pub fn bessel_zeros<OT: Into<f64>>(
    func_type: &BesselFunType,
    order: OT,
    n_zeros: usize,
    precision: f64,
) -> Vec<f64> {
    bessel_zeros_impl::<AmosBackend>(func_type, order.into(), n_zeros, precision)
}

/// Fast Bessel zero computation using the `real-bessel` backend.
///
/// This module provides the same API as the crate root but uses the
/// [`real_bessel`] crate internally. It is faster than the AMOS backend
/// but is restricted to **integer orders** (`i32`).
///
/// # Example
///
/// ```rust
/// use bessel_zeros::fast;
///
/// let zeros = fast::bessel_zeros_j(0, 5);
/// assert!((zeros[0] - 2.40482555769577).abs() < 1e-10);
/// ```
pub mod fast {
    use crate::{BesselFunType, DEFAULT_PRECISION, algorithm::bessel_zeros_impl};
    use crate::backend::real::RealBackend;

    /// Finds the first `n_zeros` zeros of Jn(x) using the fast real-bessel backend.
    /// Order must be an integer (`i32`).
    pub fn bessel_zeros_j(order: i32, n_zeros: usize) -> Vec<f64> {
        bessel_zeros(&BesselFunType::J, order, n_zeros, DEFAULT_PRECISION)
    }

    /// Finds the first `n_zeros` zeros of Yn(x) using the fast real-bessel backend.
    /// Order must be an integer (`i32`).
    pub fn bessel_zeros_y(order: i32, n_zeros: usize) -> Vec<f64> {
        bessel_zeros(&BesselFunType::Y, order, n_zeros, DEFAULT_PRECISION)
    }

    /// Finds the first `n_zeros` zeros of Jn'(x) using the fast real-bessel backend.
    /// Order must be an integer (`i32`).
    pub fn bessel_zeros_jp(order: i32, n_zeros: usize) -> Vec<f64> {
        bessel_zeros(&BesselFunType::JP, order, n_zeros, DEFAULT_PRECISION)
    }

    /// Finds the first `n_zeros` zeros of Yn'(x) using the fast real-bessel backend.
    /// Order must be an integer (`i32`).
    pub fn bessel_zeros_yp(order: i32, n_zeros: usize) -> Vec<f64> {
        bessel_zeros(&BesselFunType::YP, order, n_zeros, DEFAULT_PRECISION)
    }

    /// Finds the first `n_zeros` zeros of the specified Bessel function using
    /// the fast real-bessel backend.
    ///
    /// # Arguments
    ///
    /// * `func_type` - Which Bessel function (J, Y, J', Y').
    /// * `order`     - The order (must be a non-negative integer).
    /// * `n_zeros`   - Number of zeros to return.
    /// * `precision` - Relative error tolerance on each root.
    pub fn bessel_zeros(
        func_type: &BesselFunType,
        order: i32,
        n_zeros: usize,
        precision: f64,
    ) -> Vec<f64> {
        bessel_zeros_impl::<RealBackend>(func_type, order as f64, n_zeros, precision)
    }
}
