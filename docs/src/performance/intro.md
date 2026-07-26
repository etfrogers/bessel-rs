# Performance

This section documents the computational efficiency of `bessel-rs`.

We evaluate performance across two distinct dimensions:
1. **[Rust vs Fortran](./rust_vs_fortran.md)**: A direct 1-to-1 comparison of this pure Rust implementation against the legacy Fortran bindings (via `complex-bessel-rs`).
2. **[Precision: f32 vs f64](./f32_vs_f64.md)**: An analysis of the performance benefits gained when evaluating functions in single precision (`f32`) instead of double precision (`f64`).
