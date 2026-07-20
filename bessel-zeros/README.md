Bessel zeros
============

A crate for finding the zeros (roots) of the Bessel functions and their derivatives.

## Quick start

```rust
use bessel_zeros::bessel_zeros_j;

// First 5 zeros of J_0(x)
let zeros = bessel_zeros_j(0.0_f64, 5);
assert!((zeros[0] - 2.40482555769577).abs() < 1e-10);
```

## Backends

Two backends are available:

| Backend | Module | Order type | Notes |
|---|---|---|---|
| AMOS (`amos-bessel-rs`) | crate root | any `f64`-compatible | supports non-integer orders |
| real-bessel | `bessel_zeros::fast` | `i32` only | faster; integer orders only |

Use the **crate-root functions** when you need non-integer orders (e.g. `order = 1.5`).
Use the **`fast` module** for integer orders when performance matters.

```rust
use bessel_zeros::{bessel_zeros_j, fast};

// AMOS backend — any real order
let zeros = bessel_zeros_j(1.5_f64, 10);

// fast backend — integer orders only, but faster
let zeros = fast::bessel_zeros_j(1, 10);
```

## Available functions

Bessel functions of the first and second kind (Jν(x) and Yν(x)) and their derivatives 
(Jν'(x) and Yν'(x)) are supported:

| Function | Crate root | `fast` module |
|---|---|---|
| Jν(x) zeros | `bessel_zeros_j(order, n)` | `fast::bessel_zeros_j(order, n)` |
| Yν(x) zeros | `bessel_zeros_y(order, n)` | `fast::bessel_zeros_y(order, n)` |
| Jν'(x) zeros | `bessel_zeros_jp(order, n)` | `fast::bessel_zeros_jp(order, n)` |
| Yν'(x) zeros | `bessel_zeros_yp(order, n)` | `fast::bessel_zeros_yp(order, n)` |

For custom precision, use the lower-level `bessel_zeros` / `fast::bessel_zeros` functions,
which accept a `BesselFunType` and a precision argument.

```rust
use bessel_zeros::{BesselFunType, bessel_zeros, DEFAULT_PRECISION};

let zeros = bessel_zeros(&BesselFunType::J, 0.0_f64, 10, 1e-10);
```

## Algorithm

This crate implements the routine described in:

> N. M. Temme, _"An Algorithm with ALGOL 60 Program for the Computation of the
> zeros of the Ordinary Bessel Functions and those of their Derivatives"_,
> Journal of Computational Physics, **32**, 270–279 (1979)

Inspired by Adam Wyatt's Matlab version of this code.
