# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2026-07-20
### Added
- Introduced the `fast` module backed by the new `real-bessel` crate for computing zeros of integer-order Bessel functions. This provides a significant speedup over the AMOS backend.
- Added a `criterion`-based benchmark suite to compare the speed of the general AMOS backend vs the new `fast` backend.
- Added robust support for negative integer orders in the `fast` backend (using the mathematical symmetry $J_{-n}$ and $Y_{-n}$).

### Changed
- Refactored the core zero-finding algorithm (Temme's algorithm) to be generic over a `BesselBackend` trait, allowing it to seamlessly switch between AMOS and `real-bessel`.
- Tightened the AMOS-backed `bessel_zeros` API to explicitly require `order >= 0.0` (as the AMOS algorithm is undefined for non-integer negative orders), with clear panics and documentation directing users to the `fast` module for negative integers.

## [0.1.0] - 2026-05-30
### Added
- Initial release.
- Computation of zeros for Bessel functions $J_v(x)$, $Y_v(x)$, and their derivatives $J_v'(x)$, $Y_v'(x)$ using Temme's algorithm.
- Backed by the `amos-bessel-rs` crate (AMOS Fortran library).
