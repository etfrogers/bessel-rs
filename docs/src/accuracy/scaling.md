# Scaling Options

Some mathematical functions natively exhibit extreme exponential growth or decay.

For example, the Modified Bessel Functions of the First and Second Kind—$I_\nu(z)$ and $K_\nu(z)$—grow and decay exponentially as the absolute value of $z$ increases. 

When evaluating these functions at large arguments, the values easily exceed the maximum representable limit of an `f64` float ($\approx 10^{308}$), returning `Infinity` and destroying the calculation.

## `Scaling::Scaled`

To resolve this, the internal Amos algorithms can compute a mathematically scaled version of the function that strips out the dominant exponential factor, keeping the raw floating-point calculations comfortably within machine bounds.

The low-level `complex_bessel_*` APIs require you to pass a `Scaling` enum:

```rust
pub enum Scaling {
    Unscaled,
    Scaled,
}
```

If you select `Scaling::Scaled`, the algorithm will compute the Bessel function multiplied by a scaling factor to counteract the exponential growth:

| Function | Returned Scaled Value |
| :--- | :--- |
| $J_\nu(z)$ | $e^{-|\text{Im}(z)|} J_\nu(z)$ |
| $Y_\nu(z)$ | $e^{-|\text{Im}(z)|} Y_\nu(z)$ |
| $I_\nu(z)$ | $e^{-|\text{Re}(z)|} I_\nu(z)$ |
| $K_\nu(z)$ | $e^{z} K_\nu(z)$ |
| $H^{(1)}_\nu(z)$ | $e^{-\text{Im}(z)} H^{(1)}_\nu(z)$ |
| $H^{(2)}_\nu(z)$ | $e^{\text{Im}(z)} H^{(2)}_\nu(z)$ |

By retrieving this scaled value, your application can safely operate on the mantissa, or defer the exponentiation to a log-space domain where `f64` overflow is no longer a risk.

*(Note: The simplified top-level functions like `bessel_i()` automatically use `Scaling::Unscaled` for standard behavior. You must use `complex_bessel_i()` to request scaled execution).*
