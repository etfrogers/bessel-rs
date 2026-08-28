# Rust vs Fortran Performance

A common question when porting legacy scientific code to a modern language is whether the new implementation sacrifices performance for safety and ergonomics. 

To answer this, we benchmarked `amos-bessel-rs` (our pure Rust implementation) directly against `complex-bessel-rs` (which provides low-level Rust bindings to the original Amos Fortran binaries compiled with `gfortran`).

## Benchmark Results

*Note: All benchmarks were executed on an Apple M2 Pro (ARM64) using Rust 1.95.0. While absolute times will vary depending on your hardware, the comparative speedups should remain relatively consistent across modern CPU architectures.*

The benchmarks compare both implementations running the identical algorithms on two datasets:
- **Scattered**: An array of 46 highly extreme edge cases designed to trigger algorithmic branching.
- **Grid**: A realistic uniform grid of thousands of sequential points, designed to test typical throughput, caching, and branch prediction.

```text
Performance Comparison: Rust vs Fortran
--------------------------------------------------------------------
Function   | Dataset      | Rust (ns)  | Fortran (ns)   | Speedup
--------------------------------------------------------------------
J          | Scattered    | 17246.05   | 19372.35       | 1.12x
J          | Grid         | 10856966.43| 12409125.32    | 1.14x
Y          | Scattered    | 35394.00   | 37235.85       | 1.05x
Y          | Grid         | 21346008.61| 22750567.22    | 1.07x
I          | Scattered    | 14733.58   | 17271.79       | 1.17x
I          | Grid         | 10330456.25| 11532561.09    | 1.12x
K          | Scattered    | 17786.17   | 18692.93       | 1.05x
K          | Grid         | 15016672.92| 15571820.65    | 1.04x
--------------------------------------------------------------------
Overall Geometric Mean Speedup (Rust over Fortran): 1.09x
```

### Analysis

Not only does `amos-bessel-rs` match the performance of the highly-optimized Fortran legacy code, it is actually **9% faster** on average.

This is an incredible result. It proves that you can gain all the benefits of pure Rust—including `f32` generics, memory safety, excellent `cargo` integration, and idiomatic error handling—while maintaining state-of-the-art computational performance.
