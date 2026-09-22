Bessel Rust
===========

[![Crates.io](https://img.shields.io/crates/v/amos-bessel-rs.svg)](https://crates.io/crates/amos-bessel-rs)
[![Build Status](https://github.com/etfrogers/bessel-rs/actions/workflows/rust.yml/badge.svg)](https://github.com/etfrogers/bessel-rs/actions/workflows/rust.yml)
[![codecov](https://codecov.io/gh/etfrogers/bessel-rs/branch/main/graph/badge.svg)](https://codecov.io/gh/etfrogers/bessel-rs)
[![Performance & Accuracy Guide](https://img.shields.io/badge/docs-Performance_%26_Accuracy_Guide-blue)](https://etfrogers.github.io/bessel-rs/)

A workspace for Bessel function implementations in Rust.

Crates
------

- **[amos-bessel-rs](./amos-bessel-rs)**: Idiomatic Rust translation of Amos' Bessel function algorithms (`no_std` and zero-allocation / no-alloc supported).
- **[real-bessel](./real-bessel)**: Faster, real-only Bessel function implementations (pure `no_std`, zero-alloc / no heap needed).
- **[bessel-zeros](./bessel-zeros)**: Tools for finding zeros of Bessel functions (`no_std` + `alloc` supported).
- **[fortran-amos-testing](./fortran-amos-testing)**: Raw Fortran bindings for comparison and testing.

All three pure-Rust crates support `#![no_std]` environments by disabling default features (`real-bessel` and `amos-bessel-rs` support pure zero-allocation execution without a heap allocator, while `bessel-zeros` and allocating APIs in `amos-bessel-rs` use `alloc`).

See the [amos-bessel-rs README](./amos-bessel-rs/README.md) for more details on the main implementation.

📖 **Documentation:** See the [Performance & Accuracy Guide](https://etfrogers.github.io/bessel-rs/) for detailed benchmarks comparing `f32` vs `f64` precision and Pure Rust vs Legacy Fortran performance.
