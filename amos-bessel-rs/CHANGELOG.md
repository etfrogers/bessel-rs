# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
