# Loss of Significance

When translating legacy scientific code into modern memory-safe languages, fixing numeric vulnerabilities is just as important as fixing memory leaks. 

In the original Amos Fortran library, mathematical issues during runtime (such as underflows, overflows, or diverging iterations) were handled by silently writing integer error flags to an `IERR` array parameter (e.g., `IERR=3` or `IERR=4`). Because Fortran lacks native sum types, it was incredibly common for users to completely ignore the `IERR` argument and accidentally consume the garbage data left in the uninitialized result arrays.

## Native Rust Error Handling

`amos-bessel-rs` eliminates this vulnerability. It wraps all outputs in a `Result<T, BesselError<T>>`, forcing you to explicitly acknowledge when a mathematical algorithm diverges.

For large arguments `z` or very large orders `nu`, trigonometric phase reductions and continued fraction expansions begin to experience severe precision loss. 

When evaluating these boundary values, `amos-bessel-rs` catches the numerical decay and returns one of two critical errors:

### `BesselError::PartialLossOfSignificance`
This occurs when the calculation succeeded, but the argument reduction caused the internal algorithm to drop below half of machine accuracy. 

Unlike a fatal error, the `PartialLossOfSignificance` enum payload actually **contains the computed value**! 

```rust
// Internally in `types.rs`:
pub enum BesselError<T: BesselFloat = f64> {
    PartialLossOfSignificance {
        y: Vec<Complex<T>>, // The lossy computed values
        n_zeros: usize,          // Elements explicitly zeroed
    },
    // ...
}
```

This acts as a "strict warning". It acknowledges that you *can* use the value returned in `y`, but mathematically advises you that the lower bits of the float are essentially noise.

*(Note: The simplified top-level wrappers like `bessel_j()` implicitly unwrap `PartialLossOfSignificance` to return the data directly as an `Ok(T)` to match standard user expectations. To explicitly catch this warning, use the low-level `complex_bessel_j()` functions).*

### `BesselError::LossOfSignificance`
This occurs when argument reduction resulted in a *complete* loss of accuracy. The output is mathematically meaningless noise. 

In Fortran, this returned `IERR=4` and left garbage in the array. In `amos-bessel-rs`, it returns a strict `Err(LossOfSignificance)`, preventing you from accessing undefined values and halting pipeline corruption.

## Original Amos Notes on Accuracy

In most complex variable computation, one must evaluate elementary functions. When the magnitude of `z` or (effective) `order` is large, losses of significance by argument reduction occur. 

Consequently, if either one exceeds `u1 = (0.5/eps).sqrt()`, then losses exceeding half of machine precision are likely, and an error flag `PartialLossOfSignificance` is triggered, where `eps` is the machine precision (`f64::EPSILON` for double precision). If either `z` or `order` is larger than `u2 = 0.5/eps`, then all significance is lost and `LossOfSignificance` is returned.

In order to use the int function, arguments must be further restricted not to exceed the largest machine integer, `u3 = (i32::MAX as f64) * 0.5`. Thus, the magnitude of `z` and (effective) `order` is restricted by `u2.min(u3)`. With 64-bit precision, `u1`, `u2`, and `u3` are approximately `1.3e8`, `1.8e16`, and `2.1e9`.

The approximate relative error in the magnitude of a complex bessel function can be expressed by `eps * 10.0.pow(s)` where `eps` is the nominal precision and `10.0.pow(s)` represents the increase in error due to argument reduction in the elementary functions. Here, `s = 1.0.max(z.abs().log10().abs().max(order.log10().abs())` approximately. 

However, the phase angle may have only absolute accuracy. This is most likely to occur when one component (in absolute value) is larger than the other by several orders of magnitude. If one component is `10.0.pow(k)` larger than the other, then one can expect only `p.log10().abs() - k).max(0)` significant digits. Stated another way, when `k` exceeds the exponent of `p`, no significant digits remain in the smaller component. However, the phase angle retains absolute accuracy because, in complex arithmetic with precision `p`, the smaller component will not (as a rule) decrease below `p` times the magnitude of the larger component. In these extreme cases, the principal phase angle is on the order of `+p`, `-p`, `pi/2-p`, or `-pi/2+p`.
