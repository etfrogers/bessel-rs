Amos Bessel functions in idomatic Rust
======================================

[![Crates.io](https://img.shields.io/crates/v/amos-bessel-rs.svg)](https://crates.io/crates/amos-bessel-rs)
[![Build Status](https://github.com/etfrogers/bessel-rs/actions/workflows/rust.yml/badge.svg)](https://github.com/etfrogers/bessel-rs/actions/workflows/rust.yml)
[![codecov](https://codecov.io/gh/etfrogers/bessel-rs/branch/main/graph/badge.svg)](https://codecov.io/gh/etfrogers/bessel-rs)

A crate implementing idiomatic, pure Rust translations of [Amos' complex
Bessel function algorthims](https://www.netlib.org/amos/)

Background
----------

The development of this crate originally started as a a reaction to finding that:

1) no pure Rust implementations of Bessel functions existed (no longer true!)
2) all the standard library functions (in all languages) seem to use a wrapper
    around Amos' original fortran.

For example,
- [Python - Scipy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.jv.html#scipy.special.jv)
- [Julia - OpenSpecFunc](https://github.com/JuliaMath/openspecfun)
- [Rust (the best way of calculating Bessel functions at the time)](https://crates.io/crates/complex-bessel-rs/1.2.1)

There are other implentations in some cases for integer order and real argument,
but the general case was all Amos.

Since development of this crate started, [another translation](https://github.com/elgar328/complex-bessel) has been made available. Differences and similarites to this are discussed below.

The aim of this crate is to translate the Amos Fortran code
into idomatic Rust, with a Rust-style API, while retaining full compatibility with Amos' code if required. Both simplified "just works" API and 
full version, as per Amos, are available.

Alternatives
------------

To calculate Bessel functions in Rust there are now several alternatives

- [Complex Bessel rs](https://crates.io/crates/complex-bessel-rs/) - A wrapper around the Amos' Fortran functions with a Rust API.
  Good if you want guarantees that answers will be the same as Fortran, but requires a Fortran compliler in your toolchain to compile.
- [Complex Bessel](http://docs.rs/complex-bessel/latest/complex_bessel/) - A line-by-line translation of Amos code with a very good 
  [comparison tool](https://github.com/elgar328/complex-bessel-test) to confirm both accuracy and comnputational speed. Carefully optimised 
  for accuracy and speed using detailed tools (e.g. implementation of FMA) to aid the compiler. 
- [This crate](https://docs.rs/amos-bessel-rs/latest/amos_bessel_rs/) - A more idomatic translation of the Fortran code: using Rust
  tools. Relies on the compiler to optimise as best it can. A fork of the elgar328's [comparison tool](https://github.com/etfrogers/complex-bessel-test) shows similar accuracy and 
  execution speed.
- **Real Bessel** - WIP (soon to be released) crate that calculates real-only Bessel function's *J*, and *Y* for integer order. *J* takes 
- real inputs, *Y* is restircted to positive inputs (to give real answers). This implementation is faster for these simple cases.
