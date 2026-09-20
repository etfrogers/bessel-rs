# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0]

### Added
- Zero-allocation sequence evaluation APIs (`complex_bessel_j_into`, `complex_bessel_y_into`, `complex_bessel_i_into`, `complex_bessel_k_into`, `complex_hankel1_into`, `complex_hankel2_into`) that evaluate directly into caller-provided mutable slices without heap allocations.
- Full `#![no_std]` compatibility, including pure zero-allocation (`no-alloc`) support on bare-metal embedded targets (e.g. `thumbv7em-none-eabihf`) via `default-features = false`, as well as dynamic allocation via the `alloc` feature.
- Stack-allocated Small Buffer Optimization (`ScratchBuffer`) with a default capacity of 32 elements (`DEFAULT_SBO_CAP = 32`), eliminating heap allocations across internal scratchpads in `complex_bessel_y`, reflection routines, and derivative calculations for $N \le 32$.
- Exact $\pi$-scaled trigonometric suite: `sin_cos_pi`, `sin_pi`, `cos_pi`, `cis_pi`, and `from_polar_pi` in `amos::utils`. These perform exact argument reduction modulo 2, eliminating floating-point rounding errors at integer and half-integer axes, avoiding runtime divisions by $\pi$, and preventing integer overflow on large orders.
- Computation metadata via [`SequenceInfo`] (`n_zeros` and `partial_loss_of_significance`) returned from all sequence functions.
- Arbitrary order derivatives via the `derivatives` module using DLMF binomial identities, capped at order 15 to guarantee stack-only execution.

### Changed
- Replaced the fatal `BesselError::PartialLossOfSignificance` error variant with metadata in `SequenceInfo`. Converged values with degraded precision are now returned safely as `Ok((values, seq_info))` or `Ok(seq_info)` rather than triggering an error path.
- In `reflections.rs`, replaced separate `sinpi` and `cospi` calls with simultaneous `sin_cos_pi` evaluations, halving transcendental evaluations in $J$, $Y$, and Hankel reflections, and simplified $H$ reflection rotations using `cis_pi`.
- Streamlined rotation and phase factors across `algorithms.rs`, `analytic_continuation.rs`, `large_z.rs`, and `uniform.rs` using `cis_pi` and `from_polar_pi`.
- Extensive performance optimizations across core numerical routines:
  - Miller backward recurrence loops refactored into register loops with dedicated fast paths for real arguments.
  - Vectorized Horner polynomial evaluation for uniform asymptotic parameters.
  - Replaced `abs().powi(2)` with `norm_sqr()` in Wronskian and ratio evaluations, eliminating costly `hypot` operations.
  - Recurrence state in power series maintained in registers instead of memory stores.
  - Eliminated `EitherIter` dynamic dispatch in recurrence iteration.
  - Reused internal workspace buffers in-place in Wronskian algorithms.
- Refined underflow handling:
  - Corrected underflow zero counts and trailing term clearing in `scale_k_recurrence`.
  - Added early power series underflow cutoff in `i_power_series`.
  - Restored refined logarithmic underflow check in `check_underflow_uniform_asymp_params`.
  - Preserved phase multipliers in `safe_multiply` underflow branches.
- Converted error detail strings to `&'static str` for pure `no_std` zero-allocation error handling.

### Removed
- Removed deprecated `i_pow_n` helper from `amos::mod`.
- Removed superseded quadrant and integer parity helpers (`integer_part_is_odd`, `integer_half_part_is_odd`, `order_quadrant`).
- Removed `sinpi` and `cospi` wrapper functions from `reflections.rs`.

## [0.4.0] - 2024-05-15
### Added
- Native support for `f32` in all main entry points via generic type parameters, allowing seamless computation in both single and double precision.
- Converted internal math constants (`MachineConsts`) and recursive algorithms to be fully generic over floating-point types.
- Expanded the unit and grid testing suites to comprehensively cover the new `f32` implementations and compare them against the Fortran reference.

### Changed
- Exposed previously internal traits and types that are now required for the new generic `f32`/`f64` interface.
- Switch from `values` to `for` loops in various internal structures for compilation performance and readability.

## [0.3.0] - 2024-04-10
### Changed
- Completed the transition to a pure-Rust implementation! The original Fortran AMOS algorithms (`assypmtotic_i`, `overflow_check`, `i_power_series`, etc.) have been fully translated into safe Rust.
- Cleaned up the public API by hiding internal helper macros, structs, and `gamma_ln` functions.
- Expanded crate-level documentation and function-level documentation (translated from the original AMOS Fortran comments).
- Improved error handling by directly exposing `BesselError` logic and improving `PartialLossOfSignificance` detection.

## [0.2.0] - 2024-03-20
### Added
- Introduced the `BesselInput` trait and generic input/output structures to simplify function signatures.
- Added extensive integration tests and `fortran_limits` feature gating.

### Changed
- Re-architected project from the original monolithic `bessel-rs` into the current `amos-bessel-rs` and workspace structure.

## [0.1.0] - 2024-01-01
### Added
- Initial release porting AMOS Fortran library to Rust using FFI.
