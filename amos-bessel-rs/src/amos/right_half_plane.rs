use num::{Complex, Zero, complex::ComplexFloat};

use crate::{
    BesselError, Scaling,
    amos::{
        IKType,
        asymptotics::{i_asymp_large_order, i_asymptotic},
        gamma_ln,
        limits::{OverflowState, check_underflow_uniform_asymp_params},
        power_series::i_power_series,
        recurrence::{i_miller, scale_controlled_recurrence, scale_k_recurrence},
        utils::{two_over_z_safe, will_underflow},
        wronskian::i_wronskian,
    },
    types::{BesselFloat, BesselResult, BesselValues},
};

/// i_right_half_plane computes the i function in the right half z plane
///
/// i_right_half_plane  acts as a master dispatcher for computing the I_v (z)
/// sequence in the right half of the complex plane (Re (z) ≥ 0).
///
/// Because the Modified Bessel Function
///   I_v (z)
/// has wildly different behavioral regimes (growing exponentially for large z, decaying exponentially for large ν),
/// there is no single algorithm that can compute it accurately everywhere. This function checks the magnitude of
/// z and ν and routes the computation to the most numerically stable algorithm.
///
/// ### 1. The Power Series Regime (Small z)
///
/// If |z| is very small (|z| ≤ 2) or ν is relatively large compared to z, it attempts to use a standard Taylor/Power Series ( [i_power_series] ).
///
/// - If the series converges successfully, it returns.
/// - If the highest orders underflow to  0 , it subtracts them from  remaining_n  and passes the rest of the array down to the next algorithms.
///
/// ### 2. Large Argument Asymptotics (Large z, Small ν)
///
/// If |z| is massive (above the machine's asymptotic limit) and ν is relatively small (2|z| ≥ ν²),
/// it uses Hankel's asymptotic expansion for large arguments ( [i_asymptotic] ).
///
/// ### 3. Underflow Truncation
///
/// Before running any intermediate recurrences, if νₘₐₓ > 1, it runs a quick check ( [check_underflow_uniform_asymp_params] ).
/// Because I_v (z) decays extremely fast as ν → ∞, this calculates exactly how many of the highest-order requested values
/// will mathematically underflow to `0.0`.
/// It fills those with `0.0` , updates `remaining_n`, and saves the system from doing pointless, unstable math.
///
/// ### 4. Large Order Asymptotics (Extremely Large ν or z)
///
/// If the requested order or argument exceeds the machine's strict asymptotic limit, standard recurrences break down.
/// It delegates to [i_asymp_large_order], which uses Debye polynomials (Uniform Asymptotic Expansions) to bootstrap
/// the recurrence from a very high order.
///
/// ### 5. Miller's Algorithm (Intermediate Regime)
///
/// If it reaches this point, we are in the "messy middle" where asymptotics fail. It uses Miller's Algorithm,
/// which runs a backward recurrence from a dynamically calculated starting bound.
/// But because Miller's algorithm only gives relative (unscaled) values, they must be normalized.
/// It chooses between two normalization methods:
///
/// - Series Normalization ( [i_miller] ): If |z| is small enough, it normalizes the sequence by
///   plugging it into the known Neumann series identity ∑Iₙ(z) = eᶻ.
/// - Wronskian Normalization ( [i_wronksian] ): For larger |z|, the Neumann series
///   requires too many terms to converge. Instead, it computes K_v (z)
///   independently, and normalizes the I sequence using the Wronskian cross-product identity:
///   I_v K_v+1 + I_v+1 K_v = 1/z
///
/// Originally ZBINU
pub(crate) fn i_right_half_plane<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    n: usize,
) -> BesselResult<T, usize> {
    let mut n_zeros = 0;
    let abs_z = z.abs();
    let mut remaining_n: usize = n;
    let mut max_order = order + T::from_usize(n - 1);
    let mut y = T::c_zeros(n);
    if abs_z <= T::two() || abs_z.powi(2) * T::from_f64(0.25) <= max_order + T::ONE {
        // Power series for small z
        let n_zeros_inner;
        // i_power_series return *signed* n_zeros. As per the docs
        // n_zeros > 0 means that the last n_zeros components were set to zero
        // due to underflow. (As is the normal convention)
        // n_zeros < 0 means underflow occurred, but the
        // condition z.abs() <= 2*(order+1).sqrt() was violated and the
        // computation must be completed in another routine with n=n-abs(n_zeros).
        (y, n_zeros_inner) = i_power_series(z, order, scaling, remaining_n)?;
        let calculation_finished = n_zeros_inner >= 0;
        let n_to_complete: usize = n_zeros_inner.unsigned_abs();
        n_zeros += n_to_complete;
        remaining_n -= n_to_complete;
        if remaining_n == 0 || calculation_finished {
            return Ok((y, n_zeros));
        }
        max_order = order + (T::from_usize(remaining_n) - T::ONE);
    }

    if (abs_z >= T::MACHINE_CONSTANTS.asymptotic_z_limit)
        && ((max_order <= T::ONE) || (max_order.powi(2) <= abs_z + abs_z))
    {
        // Large Argument Asymptotics (Large z, Small order)
        let (cy, n_zeros_asymptotic) = i_asymptotic(z, order, scaling, remaining_n)?;
        debug_assert!(n_zeros_asymptotic == n_zeros);
        return Ok((cy, n_zeros));
    }

    if max_order > T::ONE {
        // Overflow and underflow test on I sequence for Miller algorithm
        let n_zeros_underflow = check_underflow_uniform_asymp_params(
            z,
            order,
            scaling,
            IKType::I,
            remaining_n,
            &mut y,
        )?;
        n_zeros += n_zeros_underflow;
        remaining_n -= n_zeros_underflow;
        if remaining_n == 0 {
            return Ok((y, n_zeros));
        }
        max_order = order + T::from_usize(remaining_n - 1);
    }

    if (max_order > T::MACHINE_CONSTANTS.asymptotic_order_limit)
        || (abs_z > T::MACHINE_CONSTANTS.asymptotic_order_limit)
    {
        let (n_zeros_asymp_lo, remaining_n) =
            i_asymp_large_order(z, order, scaling, remaining_n, &mut y)?;
        n_zeros += n_zeros_asymp_lo;
        if remaining_n == 0 {
            return Ok((y, n_zeros));
        }
    }

    if max_order <= T::ONE && abs_z <= T::MACHINE_CONSTANTS.asymptotic_z_limit {
        // Miller algorithm with series normalization
        let y = i_miller(z, order, scaling, remaining_n)?;
        return Ok((y, n_zeros));
    }

    // Miller algorithm normalized by the Wronskian
    let n_zeros_wr = i_wronskian(z, order, scaling, remaining_n, &mut y)?;
    Ok((y, n_zeros + n_zeros_wr))
}

// Could be moved to utils if it's used outside this file
#[inline]
fn one_over_sinc<T: BesselFloat>(x: T) -> T {
    if x == T::ZERO {
        T::ONE
    } else {
        let pi_x = x * T::PI();
        pi_x / pi_x.sin()
    }
}

/// Computes a sequence of Modified Bessel functions of the second kind, $K_{\nu+k}(z)$
/// for $k = 0, 1, \dots, n-1$, in the right-half plane ($\text{Re}(z) \ge 0$).
/// Originally Amos routine `ZBKNU`.
pub fn k_right_half_plane<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    n: usize,
) -> Result<BesselValues<T, usize>, BesselError<T>> {
    let sqrt_pi_over_2: T = T::from_f64(1.253_314_137_315_500_3);

    let abs_z = z.abs();
    let mut n_zeros = 0;
    let mut underflow_occurred = false;
    let mut overflow_state;
    let two_over_z = two_over_z_safe(z);
    let mut integer_order = (order.round()).to_isize().unwrap(); // round to nearest int
    let small_order_n_eq_1 = integer_order == 0 && n == 1;

    let signed_fractional_order = order - T::from_f64(integer_order as f64); // signed fractional part (-0.5 <= x <= 0.5)
    let frac_order_sqr = if signed_fractional_order.abs() > T::MACHINE_CONSTANTS.abs_error_tolerance
    {
        signed_fractional_order.powi(2)
    } else {
        T::ZERO
    };

    let (mut k_v_minus_1, mut k_v) = if (signed_fractional_order.abs() != T::half())
        && (abs_z <= T::two())
    {
        // Small |z| <= 2.0 (and non-half-integer order): compute seed values K_nu and K_{nu+1}
        // using Temme's power series expansion.
        let (mut s1, mut s2, shinc_mu) =
            compute_small_z_power_series(z, frac_order_sqr, signed_fractional_order);

        if small_order_n_eq_1 {
            // Fast exit: order is small (integer_order == 0) and only n = 1 value was requested,
            // so we can return K_nu directly without running forward recurrence.
            let mut y = s1;
            if scaling == Scaling::Scaled {
                y *= z.exp();
            }
            return Ok((vec![y], 0));
        }

        overflow_state =
            if (order + T::ONE) * shinc_mu.re.abs() > T::MACHINE_CONSTANTS.approximation_limit {
                OverflowState::NearOver
            } else {
                OverflowState::None
            };
        s2 *= overflow_state.scaling_factor::<T>() * two_over_z;
        s1 *= overflow_state.scaling_factor::<T>();
        if scaling == Scaling::Scaled {
            let z_exp = z.exp();
            s1 *= z_exp;
            s2 *= z_exp;
        }
        (s1, s2)
    } else {
        // Large |z| > 2.0 or half-integer orders: compute starting seeds via Miller's algorithm
        // or the exact asymptotic coefficient.
        // If Re(z) is large, enable scaled computation to prevent premature underflow.
        let mut coeff = Complex::<T>::new(sqrt_pi_over_2, T::ZERO) / z.sqrt();
        overflow_state = OverflowState::None;
        if scaling == Scaling::Unscaled {
            if z.re > T::MACHINE_CONSTANTS.approximation_limit {
                underflow_occurred = true;
                overflow_state = OverflowState::NearUnder;
            } else {
                coeff *= overflow_state.scaling_factor::<T>() * (-z).exp();
            }
        }
        let order_rotation = (signed_fractional_order * T::PI()).cos().abs();
        let quarter_minus_nu_sqr = (T::from_f64(0.25) - frac_order_sqr).abs();

        if signed_fractional_order.abs() == T::half()
            || order_rotation == T::ZERO
            || quarter_minus_nu_sqr == T::ZERO
        {
            (coeff, coeff)
        } else {
            compute_large_z_miller_seeds(
                z,
                signed_fractional_order,
                frac_order_sqr,
                order_rotation,
                quarter_minus_nu_sqr,
                coeff,
                small_order_n_eq_1,
            )?
        }
    };

    // Starting seeds k_v_minus_1 (K_nu) and k_v (K_{nu+1}) are ready; proceed to recurrence.
    if n == 1 {
        integer_order -= 1
    };

    if !small_order_n_eq_1 {
        if integer_order > 0 {
            let mut next_offset = 2;
            if underflow_occurred {
                // Scaled forward recurrence to climb out of the underflow region.
                // Once two consecutive non-underflowing values are found and stored in `recovery_buffer`,
                // standard recurrence can safely resume.
                underflow_occurred = false;
                let mut recovery_buffer = [T::C_ZERO; 2];
                let half_exponent_limit = T::half() * T::MACHINE_CONSTANTS.exponent_limit;

                let abs_limit = (-T::MACHINE_CONSTANTS.exponent_limit).exp();

                let mut z_shift = z;
                // Tracks the index of the previous non-underflowing value to ensure two in a row
                let mut last_good_iteration: isize = -1;
                // Alternates between 0 and 1 (via buf_idx = 1 - buf_idx) to buffer two
                // consecutive values in recovery_buffer
                let mut buf_idx = 1;

                // This is the mathematical offset we are about to calculate.
                // We already have offset 0 (k_v_minus_1) and offset 1 (k_v), so we start at 2!
                for offset in 2..=integer_order + 1 {
                    // Record that if we break, the next loop should pick up after us
                    next_offset = offset + 1;

                    let recurrence_factor =
                        (signed_fractional_order + T::from_isize(offset - 1)) * two_over_z;
                    (k_v_minus_1, k_v) = (k_v, k_v * recurrence_factor + k_v_minus_1);
                    let ln_abs_k_v = k_v.abs().ln();
                    if -z_shift.re + ln_abs_k_v >= -T::MACHINE_CONSTANTS.exponent_limit {
                        let trial_k_v =
                            (-z_shift + k_v.ln()).exp() / T::MACHINE_CONSTANTS.abs_error_tolerance;
                        if !will_underflow(trial_k_v) {
                            buf_idx = 1 - buf_idx;
                            recovery_buffer[buf_idx] = trial_k_v;
                            // Two consecutive non-underflowing values found
                            if last_good_iteration == offset - 2 {
                                break;
                            } else {
                                last_good_iteration = offset - 1;
                                continue;
                            }
                        }
                        if ln_abs_k_v < half_exponent_limit {
                            continue;
                        }
                        z_shift.re -= T::MACHINE_CONSTANTS.exponent_limit;
                        k_v_minus_1 *= abs_limit;
                        k_v *= abs_limit;
                    }
                }
                overflow_state = OverflowState::NearUnder;

                k_v = recovery_buffer[buf_idx];
                buf_idx = 1 - buf_idx;
                k_v_minus_1 = recovery_buffer[buf_idx];
            }

            // Resume the standard forward recurrence: K_{v+1} = (2v/z) * K_{v} + K_{v-1}
            // Because K_v grows extremely rapidly with order v, we must dynamically
            // check and scale the values downwards at each step to prevent floating-point overflow.
            let mut reciprocal_scaling_factor = overflow_state.reciprocal_scaling_factor::<T>();
            let mut boundary = overflow_state.boundary();
            for offset in next_offset..=integer_order + 1 {
                let recurrence_factor =
                    (signed_fractional_order + T::from_isize(offset - 1)) * two_over_z;
                (k_v_minus_1, k_v) = (k_v, recurrence_factor * k_v + k_v_minus_1);
                let scaled_k_v = k_v * reciprocal_scaling_factor;

                overflow_state.scale_recurrence(
                    &mut k_v_minus_1,
                    &mut k_v,
                    scaled_k_v,
                    &mut boundary,
                    &mut reciprocal_scaling_factor,
                );
            }
        }
        if n == 1 {
            k_v_minus_1 = k_v;
        }
    }

    let mut y = T::c_zeros(n);

    let n_completed = if underflow_occurred {
        // Seed output array with the starting values
        y[0] = k_v_minus_1;
        if n > 1 {
            y[1] = k_v;
        }
        // Step up through orders until we find two values that don't underflow.
        // scale_k_recurrence places them into y scaled by abs_error_tolerance, which we unscale below.
        scale_k_recurrence(z, order, n, &mut y, &mut n_zeros, two_over_z);
        let n_non_zero = (n - n_zeros) as isize;
        if n_non_zero <= 0 {
            return Ok((y, n_zeros));
        }

        // Unscale the first two valid values by multiplying by abs_error_tolerance
        let mut working_index = n_zeros;
        k_v_minus_1 = y[working_index];
        y[working_index] *= T::MACHINE_CONSTANTS.abs_error_tolerance;
        if n_non_zero > 1 {
            working_index += 1;
            k_v = y[working_index];
            y[working_index] *= T::MACHINE_CONSTANTS.abs_error_tolerance;
        }
        if n_non_zero > 2 {
            // If some values underflowed, the first non-zero values are near the underflow boundary
            overflow_state = OverflowState::NearUnder;
        }
        if n <= 2 {
            return Ok((y, n_zeros));
        }

        working_index + 1
    } else {
        // No underflow occurred: unscale and fill output array
        y[0] = k_v_minus_1 * overflow_state.reciprocal_scaling_factor::<T>();
        if n == 1 {
            return Ok((y, n_zeros));
        }
        y[1] = k_v * overflow_state.reciprocal_scaling_factor::<T>();
        if n == 2 {
            return Ok((y, n_zeros));
        }
        2
    };

    // Fill in remaining values via forward recurrence with scaling control
    scale_controlled_recurrence(
        true,
        order,
        z,
        &mut y,
        n_completed,
        k_v_minus_1,
        k_v,
        overflow_state,
    );
    Ok((y, n_zeros))
}

/// Computes seed values $K_\nu(z)$ and $K_{\nu+1}(z)$ (unscaled) for small $|z| \le 2$
/// using Nico Temme's (1975) power series algorithm.
fn compute_small_z_power_series<T: BesselFloat>(
    z: Complex<T>,
    frac_order_sqr: T,
    signed_fractional_order: T,
) -> (Complex<T>, Complex<T>, Complex<T>) {
    const MAX_ITERATIONS: usize = 1000;
    let gamma_difference_taylor_coeffs: [T; 8] = [
        T::from_f64(5.772_156_649_015_329e-1),
        T::from_f64(-4.200_263_503_409_524e-2),
        T::from_f64(-4.219_773_455_554_433e-2),
        T::from_f64(7.218_943_246_663_1e-3),
        T::from_f64(-2.152_416_741_149_509_8e-4),
        T::from_f64(-2.013_485_478_078_824e-5),
        T::from_f64(1.133_027_231_981_696e-6),
        T::from_f64(6.116_095_104_481_416e-9),
    ];

    let two_over_z = two_over_z_safe(z);

    // mu = nu * ln(2/z), so exp(mu) = (2/z)^nu and exp(-mu) = (z/2)^nu
    let mu = two_over_z.ln() * signed_fractional_order;
    let cosh_mu = mu.cosh();
    let one_over_sinc_nu = one_over_sinc(signed_fractional_order);
    // shinc_mu = sinh(mu) / nu = ((2/z)^nu - (z/2)^nu) / (2*nu)
    // When nu -> 0, lim sinh(mu)/nu = ln(2/z)
    let shinc_mu = if signed_fractional_order == T::ZERO {
        two_over_z.ln()
    } else {
        mu.sinh() / signed_fractional_order
    };
    // Compute 1/Gamma(1+nu) and 1/Gamma(1-nu) using Euler's reflection formula:
    // Gamma(1-nu) * Gamma(1+nu) = pi*nu / sin(pi*nu)
    let recip_gamma_one_plus_nu = (-gamma_ln(T::ONE + signed_fractional_order).unwrap()).exp();
    let recip_gamma_one_minus_nu = T::ONE / (recip_gamma_one_plus_nu * one_over_sinc_nu);

    // Compute (1/Gamma(1-nu) - 1/Gamma(1+nu)) / (2*nu).
    // When |nu| <= 0.1, use Taylor series expansion to avoid 0/0 numerical indeterminacy.
    let gamma_diff_over_2nu = if signed_fractional_order.abs() <= T::from_f64(0.1) {
        let mut sum = T::ZERO;
        for (i, cc) in gamma_difference_taylor_coeffs.iter().enumerate() {
            let term = *cc * frac_order_sqr.powi(i as i32);
            sum += term;
            if term.abs() < T::MACHINE_CONSTANTS.abs_error_tolerance {
                break;
            }
        }
        -sum
    } else {
        (recip_gamma_one_minus_nu - recip_gamma_one_plus_nu) / (T::two() * signed_fractional_order)
    };
    let mean_gamma = (recip_gamma_one_minus_nu + recip_gamma_one_plus_nu) * T::half();

    // Initial terms of Temme's series at k=0:
    // temme_coeff = f_0 = cancellation-free combination of (z/2)^-nu and (z/2)^+nu
    // neg_order_term = p_0 = 0.5 * (z/2)^-nu / Gamma(1+nu)
    // pos_order_term = q_0 = 0.5 * (z/2)^+nu / Gamma(1-nu)
    let mut temme_coeff =
        one_over_sinc_nu * (gamma_diff_over_2nu * cosh_mu + mean_gamma * shinc_mu);
    let mut neg_order_term = T::half() * mu.exp() / recip_gamma_one_plus_nu;
    let mut pos_order_term = (T::half() / mu.exp()) / recip_gamma_one_minus_nu;

    let mut taylor_factor = T::C_ONE;
    let mut term_magnitude = T::ONE;
    let z_sqr_over_4 = T::from_f64(0.25) * z.powu(2);
    let abs_z = z.abs();
    let abs_z_sqr_over_4 = T::from_f64(0.25) * abs_z * abs_z;

    let mut sum_k_nu = temme_coeff;
    let mut sum_k_nu_plus_1 = neg_order_term;
    if abs_z >= T::MACHINE_CONSTANTS.abs_error_tolerance {
        for step in 1..MAX_ITERATIONS {
            let k = T::from_usize(step);
            let k_sqr_minus_nu_sqr = k.powi(2) - frac_order_sqr;
            temme_coeff = (temme_coeff * k + neg_order_term + pos_order_term) / k_sqr_minus_nu_sqr;
            neg_order_term /= k - signed_fractional_order;
            pos_order_term /= k + signed_fractional_order;
            taylor_factor *= z_sqr_over_4 / k;
            sum_k_nu += taylor_factor * temme_coeff;
            sum_k_nu_plus_1 += taylor_factor * (neg_order_term - k * temme_coeff);
            term_magnitude *= abs_z_sqr_over_4 / k;

            if term_magnitude <= T::MACHINE_CONSTANTS.abs_error_tolerance {
                break;
            }
        }
    }
    (sum_k_nu, sum_k_nu_plus_1, shinc_mu)
}

/// Computes seed values $K_\nu(z)$ and $K_{\nu+1}(z)$ for $|z| > 2$
/// using Miller's backward recurrence algorithm normalized by the asymptotic series identity.
fn compute_large_z_miller_seeds<T: BesselFloat>(
    z: Complex<T>,
    signed_fractional_order: T,
    frac_order_sqr: T,
    order_rotation: T,
    quarter_minus_nu_sqr: T,
    coeff: Complex<T>,
    small_order_n_eq_1: bool,
) -> Result<(Complex<T>, Complex<T>), BesselError<T>> {
    let starting_k = determine_miller_starting_k(z, frac_order_sqr, order_rotation)?;
    // Now we have starting_k, run the backward recurrence loop
    // to determine the normalization factor and find K_nu, K_{nu+1}
    let mut unnormalized_k_plus_1 = Complex::<T>::zero();
    let mut unnormalized_k = Complex::<T>::new(T::MACHINE_CONSTANTS.abs_error_tolerance, T::ZERO);
    let mut normalization_sum = unnormalized_k;
    for k_int in (1..=starting_k).rev() {
        let k = T::from_usize(k_int);
        let k_sqr = k.powi(2);
        let backward_recurrence_factor = (z + k) * T::two() / (k + T::ONE);
        (unnormalized_k_plus_1, unnormalized_k) = (
            unnormalized_k,
            (unnormalized_k * backward_recurrence_factor - unnormalized_k_plus_1) * (k_sqr + k)
                / (k_sqr - k + quarter_minus_nu_sqr),
        );
        normalization_sum += unnormalized_k;
    }
    // Normalize the unscaled K_nu using the accumulated sum: K_nu = (P_0 / sum) * coeff
    let mut k_nu = unnormalized_k / normalization_sum.abs();

    normalization_sum = normalization_sum.conj() / normalization_sum.abs();
    k_nu *= coeff * normalization_sum;
    let k_nu_plus_1 = if small_order_n_eq_1 {
        T::C_ZERO
    } else {
        // Numerically stable ratio (P_1 / P_0)
        unnormalized_k_plus_1 /= unnormalized_k.abs();
        unnormalized_k = unnormalized_k.conj() / unnormalized_k.abs();
        (((-(unnormalized_k_plus_1 * unnormalized_k) + signed_fractional_order + T::half()) / z)
            + T::ONE)
            * k_nu
    };
    Ok((k_nu, k_nu_plus_1))
}

fn determine_miller_starting_k<T: BesselFloat>(
    z: Complex<T>,
    frac_order_sqr: T,
    order_rotation: T,
) -> Result<usize, BesselError<T>> {
    let abs_z = z.abs();
    const K_MAX: usize = 30;
    let miller_truncation_heuristic_1: T = T::from_f64(1.909_859_317_102_744);
    let miller_truncation_heuristic_2: T = T::from_f64(1.897_699_993_315_177_5);

    // Compute recurrence_threshold = f(E).
    // If |z| >= R2, use forward recurrence to determine the backward truncation index K.
    // The `recurrence_threshold` is a linear function over the mantissa bits E (12 <= E <= 60).
    let bits = (T::MANTISSA_DIGITS - 1).clamp(12, 60) as f64;
    let recurrence_threshold = T::from_f64((2.0 / 3.0) * bits - 6.0);
    let arg_z = z.arg();

    // Both blocks below are answering the question:
    // How large does our starting index K need to be to achieve K_(nu+K) ≈ 0
    // (that is, equals zero to within machine precision)
    let starting_k = if abs_z > recurrence_threshold {
        // Forward recurrence loop to determine starting_k when z.abs() >= recurrence_threshold
        let convergence_test =
            order_rotation / (T::PI() * abs_z * T::MACHINE_CONSTANTS.abs_error_tolerance);
        if convergence_test <= T::ONE {
            // skip the loop and just return trial index as 1
            1
        } else {
            let mut trial_index = 1;
            let mut p_prev = T::ZERO;
            let mut p_curr = T::ONE;
            let mut converged = false;
            for i in 1..=K_MAX {
                let i_plus_1 = T::from_usize(i + 1);
                let i_float = T::from_usize(i);
                let p_coeff_a = ((i_float - T::half()).powi(2) - frac_order_sqr)
                    / (i_float * (i_float + T::ONE));
                let p_coeff_b = (T::two() * abs_z + T::from_usize(2 * i)) / i_plus_1;
                (p_prev, p_curr) = (p_curr, p_coeff_b * p_curr - p_coeff_a * p_prev);
                trial_index = i + 1;
                if convergence_test < p_curr.abs() * i_plus_1 {
                    converged = true;
                    break;
                }
            }
            if !converged {
                return Err(BesselError::DidNotConverge);
            }
            let raw_k = T::from_usize(trial_index)
                + miller_truncation_heuristic_1 * arg_z * (recurrence_threshold / abs_z).sqrt();
            raw_k.to_usize().unwrap()
        }
    } else {
        // For small z.abs() (< recurrence threshold), we don't bother running the loop above;
        // instead, we use a heuristic equation to calculate the K value directly.
        let precision_factor = order_rotation * miller_truncation_heuristic_2
            / (T::MACHINE_CONSTANTS.abs_error_tolerance * abs_z.sqrt().sqrt());
        let angle_correction_a = T::from_f64(3.0) * arg_z / (T::ONE + abs_z);
        let angle_correction_b = T::from_f64(14.7) * arg_z / (T::from_f64(28.0) + abs_z);
        let heuristic_curve_factor = (precision_factor.ln()
            + abs_z * angle_correction_a.cos() / (T::ONE + T::from_f64(0.008) * abs_z))
            / angle_correction_b.cos();
        let raw_k =
            T::from_f64(0.12125) * heuristic_curve_factor.powi(2) / abs_z + T::from_f64(1.5);
        raw_k.to_usize().unwrap()
    };
    Ok(starting_k)
}
