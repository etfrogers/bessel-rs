#![warn(missing_docs)]
//! Real-valued Bessel functions of the first and second kinds.
//!
//! This crate provides `j0`, `j1`, `jn`, `y0`, `y1` and `yn` for integer
//! orders and real (`f64`) arguments.
//!
//! ## Return types
//!
//! The **J functions** (`j0`, `j1`, `jn`) are real-valued for all real x and
//! return `f64` directly.
//!
//! The **Y functions** (`y0`, `y1`, `yn`) are real-valued only for positive x.
//! For `x ≤ 0` they return `Err(`[`BesselError::NegativeInputForYFunction`]`)`
//! rather than a complex result.
//!
//! ## Special cases shared by all functions
//!
//! | Input | J functions | Y functions |
//! |-------|------------|-------------|
//! | `NaN` | `NaN` | `Ok(NaN)` |
//! | `+∞`  | `0.0` | `Ok(0.0)` |
//! | `x ≤ 0` | real result | `Err(…)` |
//! | `x = 0` | `j0 → 1`, `j1/jn → 0` | `Ok(−∞)` |
//!
//! If you need complex arguments, non-integer orders, or other Bessel
//! varieties (I, K, H), see the
//! [`amos-bessel-rs`](https://crates.io/crates/amos-bessel-rs) crate.
//!
//! ## Example
//!
//! ```
//! use real_bessel::j0;
//!
//! let result = j0(1.0);
//! assert!((result - 0.7651976865579612).abs() < 1e-10);
//! ```
mod types;
pub use types::BesselError;

mod j0;
mod j1;
mod jn;
pub use j0::j0;
pub use j0::y0;
pub use j1::j1;
pub use j1::y1;
pub use jn::jn;
pub use jn::yn;

pub(crate) const TWO_M13: f64 = 1.0 / (1 << 13) as f64; // 2**-13 0x3f20000000000000
pub(crate) const TWO_M27: f64 = 1.0 / (1 << 27) as f64; // 2**-27 0x3e40000000000000
pub(crate) const TWO_M29: f64 = 1.0 / (1 << 29) as f64; // 2**-29
pub(crate) const TWO_M54: f64 = 1.0 / (1_u64 << 54) as f64; // 2**-54 0x3c90000000000000
pub(crate) const TWO_129: f64 = 6.805_647_338_418_77e38; // 2**129 0x4800000000000000
pub(crate) const TWO_302: f64 = 8.148_143_905_337_944e90; // 2**302 0x52D0000000000000
