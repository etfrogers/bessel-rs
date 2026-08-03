//! Large-order uniform asymptotic expansions for modified Bessel functions $I_\nu(z)$ and $K_\nu(z)$.
//!
//! Originally Amos routines `ZBUNI` (for $I_\nu$) and `ZBUNK` (for $K_\nu$).
//!
//! ### Method Overview
//! Uniform asymptotic expansions (Olver 1954, NIST DLMF 10.20 / 10.41) are accurate when the order
//! $\nu \ge \text{asymptotic\_order\_limit} = 50.0$.
//!
//! - **Domain Splitting**:
//!   - **Standard Sector** ($|\arg(z)| \le \pi/3$): Direct Debye/Airy uniform asymptotics for $I_\nu$ and $K_\nu$
//!     via `i_uniform_asymp1` / `k_uniform_asymp1`.
//!   - **Imaginary-Dominant Sector** ($\pi/3 < |\arg(z)| \le \pi/2$): Uses analytic continuation via rotated
//!     $J_\nu(z e^{\pm i\pi/2})$ and $H^{(2)}_\nu(z e^{\pm i\pi/2})$ via `i_uniform_asymp2` / `k_uniform_asymp2`
//!     to eliminate loss of significance near the imaginary axis.
//!
//! - **Order Thresholding and Backward Recurrence for $I_\nu$**:
//!   Because $I_\nu(z)$ decays with increasing order $\nu$, backward recurrence is numerically stable (Miller's method).
//!   If $\nu_{\max} = \nu + n - 1 < 50$, the order is inflated to $\nu_{\text{adj}} \ge 50$ to evaluate the
//!   asymptotic seed values, then stepped backward to the target orders with scale control.

use num::{Complex, complex::ComplexFloat};

use crate::{
    BesselError, Scaling,
    amos::{
        RotationDirection,
        asymptotics::{i_uniform_asymp1, i_uniform_asymp2, k_uniform_asymp1, k_uniform_asymp2},
        limits::OverflowState,
        recurrence::scale_controlled_recurrence,
        utils::imaginary_dominant,
    },
    types::{BesselFloat, BesselResult},
};

/// Computes the $I$ Bessel sequence $[I_\nu(z), \dots, I_{\nu+n-1}(z)]$ for large argument/order.
///
/// If $\nu_{\max} = \nu + n - 1 < 50$, the order is temporarily increased to `max_order_adjusted` $\ge 50$
/// to evaluate the uniform asymptotic expansions, and then stepped backward to the target orders
/// via scale-controlled recurrence.
///
/// Uses $J_\nu(z e^{\pm i\pi/2})$ when $z$ is imaginary dominant ($\pi/3 < |\arg(z)| \le \pi/2$),
/// and $I_\nu(z)$ directly when $|\arg(z)| \le \pi/3$.
///
/// Originally ZBUNI.
///
/// # Arguments
/// * `z` - Complex argument $\operatorname{Re}(z) \ge 0$.
/// * `order` - Starting order $\nu \ge 0$.
/// * `scaling` - Exponential scaling mode.
/// * `n` - Number of terms to compute in the sequence $[I_\nu(z), \dots, I_{\nu+n-1}(z)]$.
/// * `y` - Output slice of length $\ge n$.
///
/// # Returns
/// `Ok((n_zeros, n_calculated))` where:
/// - `n_zeros` is the number of leading components set to zero due to underflow.
/// - `n_calculated` is the number of elements successfully computed.
pub(crate) fn i_asymp_large_order<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    n: usize,
    y: &mut [Complex<T>],
) -> Result<(usize, usize), BesselError<T>> {
    // To evaluate the uniform asymptotic expansions accurately, the order must be >= asymptotic_order_limit.
    // If it isn't, we increment it up to that limit, evaluate the seeds, then recur backward to find
    // the values at the requested orders.
    let max_order = order + T::from_usize(n - 1);

    let steps_to_asymptotic_limit =
        ((T::MACHINE_CONSTANTS.asymptotic_order_limit - max_order).trunc() + T::one())
            .max(T::zero())
            .to_usize()
            .unwrap();

    let imaginary_dominant = imaginary_dominant(z);

    let evaluate_asymp = |target_order: T, count: usize, out: &mut [Complex<T>]| {
        if imaginary_dominant {
            // Asymptotic expansion for J_nu(z * e^{m*pi/2}) for large nu (pi/3 < |arg(z)| <= pi/2)
            i_uniform_asymp2(z, target_order, scaling, count, out)
        } else {
            // Asymptotic expansion for I_nu(z) for large nu (|arg(z)| <= pi/3)
            i_uniform_asymp1(z, target_order, scaling, count, out)
        }
    };
    if steps_to_asymptotic_limit == 0 {
        return evaluate_asymp(order, n, y);
    }

    let steps_float = T::from_usize(steps_to_asymptotic_limit);
    let max_order_adjusted = max_order + steps_float;
    let mut seeds = [T::C_ZERO; 2];

    // Evaluate asymptotic seeds: seeds[0] = I_{ν_adj}(z), seeds[1] = I_{ν_adj + 1}(z)
    let (n_zeros, n_last) = evaluate_asymp(max_order_adjusted, 2, &mut seeds)?;
    if n_zeros != 0 {
        return Ok((0, n));
    }

    // Pre-scale seeds if they are near underflow/overflow thresholds to protect the recurrence
    let overflow_state = if seeds[0].abs() <= T::MACHINE_CONSTANTS.absolute_approximation_limit {
        OverflowState::NearUnder
    } else if seeds[0].abs() >= T::ONE / T::MACHINE_CONSTANTS.absolute_approximation_limit {
        OverflowState::NearOver
    } else {
        OverflowState::None
    };

    let scaling_factor = overflow_state.scaling_factor::<T>();
    // In backward recurrence: s1 corresponds to order + 1, s2 corresponds to order
    let s1 = seeds[1] * scaling_factor;
    let s2 = seeds[0] * scaling_factor;

    // Stage 1: Step backward from ν_adj down to ν_max = ν + n - 1 (unbuffered)
    let (s1, s2, overflow_state) = scale_controlled_recurrence(
        false,
        max_order,
        z,
        None,
        steps_to_asymptotic_limit,
        steps_to_asymptotic_limit,
        s1,
        s2,
        overflow_state,
    );

    y[n - 1] = s2 * overflow_state.reciprocal_scaling_factor::<T>();
    if n == 1 {
        return Ok((0, n_last));
    }

    // Stage 2: Step backward from ν_max down to ν (buffered into y)
    scale_controlled_recurrence(false, order, z, Some(y), n - 1, n, s1, s2, overflow_state);

    Ok((0, n_last))
}

/// Computes the $K$ Bessel sequence $[K_\nu(z), \dots, K_{\nu+n-1}(z)]$ for $\text{order} > \text{asymptotic\_order\_limit}$ (50.0).
///
/// Unlike $I_\nu(z)$, $K_\nu(z)$ grows with increasing order $\nu$, so forward evaluation across all $n$
/// directly from the uniform asymptotic expansions is stable and does not require order inflation or backward recurrence.
///
/// Uses $H^{(2)}_\nu(z e^{\pm i\pi/2})$ when $z$ is imaginary dominant ($\pi/3 < |\arg(z)| \le \pi/2$),
/// and $K_\nu(z)$ directly when $|\arg(z)| \le \pi/3$.
///
/// Originally ZBUNK.
///
/// # Arguments
/// * `z` - Complex argument $\operatorname{Re}(z) \ge 0$.
/// * `order` - Starting order $\nu \ge 0$.
/// * `scaling` - Exponential scaling mode.
/// * `rotation` - Direction of phase rotation used for analytic continuation into the left half-plane.
/// * `n` - Number of terms to compute in the sequence.
pub(crate) fn k_asymp_large_order<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    rotation: RotationDirection,
    n: usize,
) -> BesselResult<T> {
    if imaginary_dominant(z) {
        // Asymptotic expansion for H^{(2)}_nu(z * e^{m*pi/2}) for large nu (pi/3 < |arg(z)| <= pi/2)
        k_uniform_asymp2(z, order, scaling, rotation, n)
    } else {
        // Asymptotic expansion for K_nu(z) for large nu (|arg(z)| <= pi/3)
        k_uniform_asymp1(z, order, scaling, rotation, n)
    }
}
