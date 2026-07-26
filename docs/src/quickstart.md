# Quick Start

`amos-bessel-rs` is designed to be highly ergonomic and safe, abstracting away the complex multi-return arrays of legacy Fortran into native Rust closures and enums.

## Basic Usage

The crate exposes top-level entry points for standard Bessel functions, which will automatically return the evaluated result or an explicit error type if the underlying numeric algorithm diverges.

```rust
use num::Complex;
use amos_bessel_rs::bessel_j;

fn main() {
    let order = 2.5;
    let z = Complex::new(1.0, 2.0);

    // Compute the Bessel J function
    match bessel_j(order, z) {
        Ok(result) => println!("J_{}(z) = {}", order, result),
        Err(e) => println!("Failed to compute: {}", e),
    }
}
```

## `f32` vs `f64` Generic Precision

All algorithms are generic over the `BesselFloat` trait. This means you can drop in `f32` arguments and inherently gain the `f32` speed boost without modifying your API calls:

```rust
use num::Complex;
use amos_bessel_rs::{bessel_y, BesselError};

// By passing f32 values, the `BesselFloat` trait natively switches 
// the underlying algorithmic precision to single-precision float execution.
let order: f32 = 2.5;
let z: Complex<f32> = Complex::new(1.0, 2.0);

let result: Result<Complex<f32>, BesselError<f32>> = bessel_y(order, z);
```

For workloads like computer graphics, DSP audio loops, or machine learning tensors where 32-bit floats are standard, simply providing `f32` values unlocks substantial performance wins. (See our [Performance Analysis](./performance/f32_vs_f64.md) for more).
