Amos Bessel functions in idiomatic Rust
======================================

[![Crates.io](https://img.shields.io/crates/v/amos-bessel-rs.svg)](https://crates.io/crates/amos-bessel-rs)
[![Build Status](https://github.com/etfrogers/bessel-rs/actions/workflows/rust.yml/badge.svg)](https://github.com/etfrogers/bessel-rs/actions/workflows/rust.yml)
[![codecov](https://codecov.io/gh/etfrogers/bessel-rs/branch/main/graph/badge.svg)](https://codecov.io/gh/etfrogers/bessel-rs)
[![Performance & Accuracy Guide](https://img.shields.io/badge/docs-Performance_%26_Accuracy_Guide-blue)](https://etfrogers.github.io/bessel-rs/)

A crate implementing idiomatic, pure Rust translations of [Amos' complex
Bessel function algorithms](https://www.netlib.org/amos/)

Background
----------

The development of this crate originally started as a reaction to finding that:

1) no pure Rust implementations of Bessel functions existed (no longer true!)
2) all the standard library functions (in all languages) seem to use a wrapper
    around Amos' original Fortran.

For example,
- [Python - Scipy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.jv.html#scipy.special.jv)
- [Julia - OpenSpecFunc](https://github.com/JuliaMath/openspecfun)
- [Rust (the best way of calculating Bessel functions at the time)](https://crates.io/crates/complex-bessel-rs/1.2.1)

There are other implementations in some cases for integer order and real argument,
but the general case was all Amos.

Since development of this crate started, [another translation](https://github.com/elgar328/complex-bessel) has been made available. Differences and similarities to this are discussed below.

The aim of this crate is to translate the Amos Fortran code
into idiomatic Rust, with a Rust-style API, while retaining full compatibility with Amos' code if required. Both simplified "just works" API and 
full version, as per Amos, are available.

**📖 Documentation:** Read the [Performance & Accuracy Guide](https://etfrogers.github.io/bessel-rs/) for deep dives into how `amos-bessel-rs` exceeds Fortran performance, mitigates floating-point overflow, and eliminates legacy `IERR` array vulnerabilities.

Usage
-----

```toml
[dependencies]
amos-bessel-rs = "0.4"
```

### `no_std` and `no-alloc` Support

`amos-bessel-rs` fully supports `#![no_std]` environments, both with and without dynamic memory allocation (`alloc`):

- **Pure zero-allocation `no_std` (bare-metal embedded)**:
  By disabling default features (`default-features = false`), the crate operates completely without heap allocations, suitable for bare-metal targets (e.g. `thumbv7em-none-eabihf`). You can use all slice-filling APIs (`_into`), single-value entry points (`bessel_j`, etc.), and derivatives up to order 15, powered by stack-allocated small-buffer optimization (`ScratchBuffer`).
  ```toml
  [dependencies]
  amos-bessel-rs = { version = "1.0", default-features = false }
  ```

- **`no_std` with `alloc`**:
  If a heap allocator is available in your `no_std` environment, enable the `alloc` feature to use allocating sequence functions (returning `Vec<Complex<T>>`) and dynamic buffer expansion for sequences $N > 32$:
  ```toml
  [dependencies]
  amos-bessel-rs = { version = "1.0", default-features = false, features = ["alloc"] }
  ```

When `std` is disabled, the crate relies on pure-Rust software floating-point routines via `libm` and precomputed IEEE-754 machine constants.

Alternatives
------------

To calculate Bessel functions in Rust there are now several alternatives:

- [This crate](https://docs.rs/amos-bessel-rs/latest/amos_bessel_rs/) - A modern, idiomatic pure-Rust translation of Amos' algorithms. Features zero-allocation buffer-passing APIs (`_into`), register-resident recurrence loops, and SIMD-friendly Horner polynomial evaluation. Benchmarks show it consistently matches or outpaces both Fortran AMOS and other translations while offering full `#![no_std]` support.
- [Complex Bessel](http://docs.rs/complex-bessel/latest/complex_bessel/) - A line-by-line translation of Amos code with a very good [comparison tool](https://github.com/elgar328/complex-bessel-test) to confirm both accuracy and computational speed. Optimised for accuracy and speed using FMA tools to aid the compiler.
- [Real Bessel](https://crates.io/crates/real-bessel) - A dedicated pure `core` zero-allocation crate that calculates real-only Bessel functions *J* and *Y* for integer order. Faster for these simple cases.
- [Complex Bessel rs](https://crates.io/crates/complex-bessel-rs/) - A wrapper around the Amos' Fortran functions with a Rust API. Good if you want guarantees that answers will be identical to Fortran, but requires a Fortran compiler in your toolchain to compile.
