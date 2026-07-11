Real Bessel functions in pure Rust
===================================

[![Crates.io](https://img.shields.io/crates/v/real-bessel.svg)](https://crates.io/crates/real-bessel)
[![Build Status](https://github.com/etfrogers/bessel-rs/actions/workflows/rust.yml/badge.svg)](https://github.com/etfrogers/bessel-rs/actions/workflows/rust.yml)
[![docs.rs](https://docs.rs/real-bessel/badge.svg)](https://docs.rs/real-bessel)

A pure Rust implementation of real-valued Bessel functions of the first and
second kinds, for integer orders and real arguments.

Functions provided
------------------

| Function | Description |
|----------|-------------|
| `j0(x)` | J₀(x) — order-zero, first kind |
| `j1(x)` | J₁(x) — order-one, first kind |
| `jn(n, x)` | Jₙ(x) — integer order, first kind |
| `y0(x)` | Y₀(x) — order-zero, second kind (x > 0) |
| `y1(x)` | Y₁(x) — order-one, second kind (x > 0) |
| `yn(n, x)` | Yₙ(x) — integer order, second kind (x > 0) |

The J functions accept any real `x` (including negative values and ±∞). The Y
functions are real-valued only for positive `x`; a negative argument returns a
[`BesselError`](https://docs.rs/real-bessel/latest/real_bessel/enum.BesselError.html)
rather than a complex result.

Usage
-----

```toml
[dependencies]
real-bessel = "0.1"
```

```rust
use real_bessel::{j0, j1, jn, y0, y1, yn};

// First kind — any real argument
let v = j0(1.0);                         // ≈  0.7652
let v = j1(1.0);                         // ≈  0.4401
let v = jn(5, 2.4);                      // arbitrary integer order

// Second kind — positive argument only; returns Result
let v = y0(1.0).unwrap();                // ≈  0.0883
let v = y1(1.0).unwrap();                // ≈ -0.7812
let v = yn(3, 2.0).unwrap();             // arbitrary integer order

// Negative argument for a Y function returns an error
assert!(y0(-1.0).is_err());
```

Provenance
----------

This is a Rust translation of the Bessel functions from the
[Go standard library](https://cs.opensource.google/go/go/+/master:src/math/j0.go),
which itself is a simplified translation of the C code from
[FreeBSD's libm](https://svnweb.freebsd.org/base/head/lib/msun/src/)
(`e_j0.c`, `e_j1.c`, `e_jn.c`), originally developed at SunPro /
Sun Microsystems, Inc.

The underlying polynomial approximations are well-established and have been
validated across decades of use in production systems.

Alternatives
------------

This crate deliberately covers a narrow scope. If you need something beyond
integer-order real Bessel functions, one of these is likely a better fit:

- **[amos-bessel-rs](https://crates.io/crates/amos-bessel-rs)** — complex
  arguments, non-integer orders, and all Bessel function varieties (I, K, H).
  A pure-Rust idiomatic translation of Amos' algorithms. Use this if you need
  the general case.
- **[complex-bessel-rs](https://crates.io/crates/complex-bessel-rs)** — a
  wrapper around the original Amos Fortran code with a Rust API. Requires a
  Fortran compiler in the toolchain.
- **[complex-bessel](https://docs.rs/complex-bessel)** — a line-by-line Rust
  translation of Amos, carefully optimised for accuracy and speed.

If you only need J and Y for real arguments and integer orders, this crate is
faster than the general-purpose alternatives because it uses the classical
polynomial approximations rather than the Amos recurrence machinery.

The benchmarks to compare this crate against [amos-bessel-rs](https://crates.io/crates/amos-bessel-rs)
are available in the [benches](benches) directory.

License
-------

Licensed under either of

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
- MIT license ([LICENSE-MIT](LICENSE-MIT))

at your option.

The underlying algorithms carry their own upstream notices. The full text of
the Go Authors BSD-3-clause licence and the Sun Microsystems SunPro notice
are reproduced in [NOTICES](NOTICES), as required by the BSD-3-clause licence.
