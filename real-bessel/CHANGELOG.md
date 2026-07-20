# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-07-20
### Added
- Initial release of the `real-bessel` crate.
- Idiomatic, pure-Rust implementation for calculating integer-order Bessel functions.
- Supports Bessel functions of the first kind ($J_n$), second kind ($Y_n$), and their derivatives ($J_n'$, $Y_n'$).
- Uses Go/FreeBSD-derived algorithms (forward recurrence for $n \le x$, backward recurrence via continued fractions for $n > x$).
- Does not rely on the Fortran AMOS library, making it significantly faster for integer orders and easier to compile.
- Comprehensive handling of edge cases (NaN, $\pm\infty$, 0, and negative inputs) matching the mathematical definitions.
