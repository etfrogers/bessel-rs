use std::cmp::min;

use itertools::Either;
use num::{Complex, complex::ComplexFloat};

use crate::{
    BesselError::{self, DidNotConverge},
    Scaling,
    amos::{
        gamma_ln,
        limits::OverflowState,
        utils::{two_over_z_safe, will_underflow},
    },
    types::BesselFloat,
};

/// i_miller computes the i bessel function for re(z) >= 0.0 by the
/// Miller algorithm normalized by a Neumann series.
/// The Miller algorithm relies on a brilliant trick: you start at some arbitrarily high index  N , assume
///
///    I (z) = 1
///     N
///
///  and
///
///    I   (z) = 0
///     N+1
///
/// , and then iterate backwards down to I₀(z). Because the backward recurrence is numerically stable,
/// you get the correct relative sequence. Then, you sum
/// the sequence using a known normalization identity (like the Neumann series) to find out what the
/// true scaling factor should have been, and scale all the
/// answers up to the truth.
///
/// Originally ZMLRI
pub(crate) fn i_miller<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    n: usize,
) -> Result<Vec<Complex<T>>, BesselError<T>> {
    let scale: T = T::two() * T::MIN_POSITIVE / T::MACHINE_CONSTANTS.abs_error_tolerance;
    let abs_z = z.abs();
    let int_abs_z = abs_z.to_usize().unwrap();
    let int_order = order.to_usize().unwrap();
    let modified_int_order = int_order + n - 1;
    let abs_z_plus_one = T::from_usize(int_abs_z + 1);
    let reciprocal_abs_z = T::one() / abs_z;
    let two_over_z = two_over_z_safe(z);
    let mut fwd_k_minus_1 = T::C_ZERO;
    let mut fwd_k = T::C_ONE;
    let mut abs_recurrence_factor = (abs_z_plus_one + T::one()) * reciprocal_abs_z;
    let rho =
        abs_recurrence_factor + (abs_recurrence_factor * abs_recurrence_factor - T::one()).sqrt();
    let rho_sq = rho * rho;
    let mut convergence_test = (rho_sq + rho_sq) / ((rho_sq - T::one()) * (rho - T::one()));
    convergence_test /= T::MACHINE_CONSTANTS.abs_error_tolerance;
    // Phase 1: Forward Sequence Truncation Bound
    // Run the recurrence forward. The sequence diverges, representing the rapidly growing K_nu(z).
    // We run it until the sequence exceeds a convergence threshold, which tells us how high
    // an index we need to start the backward recurrence from to ensure the truncation error
    // doesn't pollute the final values at our target index.
    let mut converged = false;
    let mut series_trunctation_index = 0;
    for i in 0..80 {
        series_trunctation_index = i + 2;
        let current_index_magnitude = abs_z_plus_one + T::from_usize(i);
        let recurrence_factor = two_over_z * ((abs_z_plus_one + T::from_usize(i * 2)) / T::two());
        (fwd_k_minus_1, fwd_k) = (fwd_k, fwd_k_minus_1 - recurrence_factor * fwd_k);
        if fwd_k.abs() > convergence_test * current_index_magnitude * current_index_magnitude {
            converged = true;
            break;
        }
    }
    if !converged {
        return Err(DidNotConverge);
    }
    let mut ratio_truncation_index = 0;
    if modified_int_order >= int_abs_z {
        // Phase 2: Forward Ratio Truncation Bound
        // If the order is very large compared to |z|, we run a secondary forward
        // recurrence to calculate an even higher truncation index specifically
        // for the Neumann normalisation sum (which requires more terms to converge).
        fwd_k_minus_1 = T::C_ZERO;
        fwd_k = T::C_ONE;
        let starting_order = T::from_f64(modified_int_order as f64) + T::one();
        convergence_test =
            (starting_order * reciprocal_abs_z / T::MACHINE_CONSTANTS.abs_error_tolerance).sqrt();
        let mut hit_loop_end = false;
        converged = false;
        for k in 0..80 {
            ratio_truncation_index = k + 1;
            let recurrence_factor =
                two_over_z * ((starting_order + T::from_usize(k * 2)) / T::two());
            (fwd_k_minus_1, fwd_k) = (fwd_k, fwd_k_minus_1 - recurrence_factor * fwd_k);
            let abs_fwd_k = fwd_k.abs();
            if abs_fwd_k < convergence_test {
                continue;
            }
            if hit_loop_end {
                converged = true;
                break;
            }
            abs_recurrence_factor = recurrence_factor.abs();
            let lambda = abs_recurrence_factor
                + (abs_recurrence_factor * abs_recurrence_factor - T::one()).sqrt();
            let kappa = abs_fwd_k / fwd_k_minus_1.abs();
            let rho = lambda.min(kappa);
            convergence_test *= (rho / (rho * rho - T::one())).sqrt();
            hit_loop_end = true;
        }
        if !converged {
            return Err(DidNotConverge);
        }
    }
    // Phase 3: Backward Recurrence and Neumann Normalisation
    // Run the backward recurrence from the truncation bound down to zero.
    // Simultaneously, accumulate the Neumann series normalisation sum, which mathematically equals e^z.
    // Dividing our unscaled sequence by this sum exactly normalises the whole array.
    let start_index =
        (series_trunctation_index + int_abs_z).max(ratio_truncation_index + modified_int_order);
    let mut kk_float = T::from_f64(start_index as f64);
    let mut val_k_plus_one = T::C_ZERO;
    // Initialize the recurrence starting values and scale them to avoid underflow
    let mut val_k = Complex::<T>::new(scale, T::ZERO);
    let fractional_order = order.fract();
    let twice_fractional_order = fractional_order + fractional_order;
    let mut binomial_coeff = (gamma_ln(kk_float + twice_fractional_order + T::one()).unwrap()
        - gamma_ln(kk_float + T::one()).unwrap()
        - gamma_ln(twice_fractional_order + T::one()).unwrap())
    .exp();
    let mut normalisation_sum = T::C_ZERO;
    // Neumann normalisation loop
    for _ in 0..(start_index - modified_int_order) {
        let pt = val_k;
        val_k = val_k_plus_one + (kk_float + fractional_order) * (two_over_z * val_k);
        val_k_plus_one = pt;
        let binomial_ratio =
            T::one() - twice_fractional_order / (kk_float + twice_fractional_order);
        let next_binomial_coeff = binomial_coeff * binomial_ratio;
        normalisation_sum += (next_binomial_coeff + binomial_coeff) * val_k_plus_one;
        binomial_coeff = next_binomial_coeff;
        kk_float -= T::one();
    }
    let mut y = T::c_zeros(n);
    y[n - 1] = val_k;
    if n != 1 {
        for i in 1..n {
            let pt = val_k;
            val_k = val_k_plus_one + (kk_float + fractional_order) * (two_over_z * pt);
            val_k_plus_one = pt;
            let binomial_ratio =
                T::one() - twice_fractional_order / (kk_float + twice_fractional_order);
            let next_binomial_coeff = binomial_coeff * binomial_ratio;
            normalisation_sum += (next_binomial_coeff + binomial_coeff) * val_k_plus_one;
            binomial_coeff = next_binomial_coeff;
            kk_float -= T::one();
            y[n - (i + 1)] = val_k;
        }
    }
    if int_order > 0 {
        for _i in 0..int_order {
            (val_k_plus_one, val_k) = (
                val_k,
                val_k_plus_one + (kk_float + fractional_order) * (two_over_z * val_k),
            );
            let binomial_ratio =
                T::one() - twice_fractional_order / (kk_float + twice_fractional_order);
            let next_binomial_coeff = binomial_coeff * binomial_ratio;
            normalisation_sum += (next_binomial_coeff + binomial_coeff) * val_k_plus_one;
            binomial_coeff = next_binomial_coeff;
            kk_float -= T::one();
        }
    }

    let mut scaled_z = z;
    if scaling == Scaling::Scaled {
        scaled_z.re = T::ZERO;
    }
    let mut ln_leading_term = -fractional_order * two_over_z.ln() + scaled_z;
    let gamma_term = gamma_ln(T::one() + fractional_order).unwrap();
    ln_leading_term -= gamma_term;
    // Calculate the final normalisation constant.
    // The complex division exp(ln_leading_term) / (normalisation_sum + val_k) is performed
    // by dividing by the magnitude twice, to avoid intermediate overflow from squaring
    // large quantities when computing the complex denominator.
    val_k += normalisation_sum;
    let sum_magnitude = val_k.abs();
    let normalization_constant =
        (ln_leading_term.exp() / sum_magnitude) * val_k.conj() / sum_magnitude;
    for element in y.iter_mut() {
        *element *= normalization_constant;
    }
    Ok(y)
}

/// i_ratios computes ratios of I bessel functions by backward
/// recurrence. The starting index is determined by forward
/// recurrence as described in J. Res. of Nat. Bur. of Standards-B,
/// Mathematical Sciences, vol 77b, p111-114, September, 1973,
/// Bessel functions I and J of complex argument and integer order,
/// by D. J. Sookne.
///
/// Originally ZRATI
pub(crate) fn i_ratios<T: BesselFloat>(z: Complex<T>, order: T, n: usize) -> Vec<Complex<T>> {
    let abs_z = z.abs();
    let integer_order = order.to_isize().unwrap();
    let modified_int_order = integer_order + n as isize - 1;
    let int_abs_z = abs_z.to_isize().unwrap();
    // starting_index is the safe starting order for the forward truncation test loop.
    // To guarantee that the terms are decaying, the truncation test must start looking at an index of at least |z| + 1.
    // However, it also must start at least as
    // high as the maximum order we actually need to calculate in our output array (which is  modified_int_order , or νₘₐₓ).
    // So starting_index is simply max (|z| + 1,νₘₐₓ).
    let starting_index = (int_abs_z + 1).max(modified_int_order);
    // After the forward loop runs for K steps starting from sarting_index,
    // we need to run the backward loop from  staring_index + n_steps  all the way down to our target order νₘₐₓ.
    // How many steps is that? The distance is (FNUP + K) - νₘₐₓ.
    //
    //  • If νₘₐₓ ≥ |z| + 1, then  FNUP  was just νₘₐₓ. The distance is exactly K. In this case,  index_difference  is clamped to 0.
    //  • If νₘₐₓ < |z| + 1, then  FNUP  was |z| + 1. The distance is K + (|z| + 1 - νₘₐₓ).
    // Notice that (|z| + 1 - νₘₐₓ) is exactly  -index_difference !
    let index_difference = modified_int_order - int_abs_z - 1;
    let index_difference = if index_difference > 0 {
        0
    } else {
        index_difference
    };

    let two_over_z = two_over_z_safe(z);
    let mut n_steps = 1;
    let mut abs_fwd_k;
    {
        // First recurr forward to find the placed to start:
        // The sequence divereges, but we want to figure out how high an index
        // K we needs to start from before the truncation error is
        // smaller than machine tolerance.
        let mut fwd_k = -two_over_z * T::from_isize(starting_index);
        let mut fwd_k_minus_1 = T::C_ONE;

        abs_fwd_k = fwd_k.abs();
        let mut abs_fwd_k_minus_1 = fwd_k_minus_1.abs();
        // Scale base_convergence_test and all subsequent fwd_k values by
        // abs_fwd_k_minus_1 to ensure that an overflow does not occur prematurely
        let initial_test_arg = (abs_fwd_k + abs_fwd_k)
            / (abs_fwd_k_minus_1 * T::MACHINE_CONSTANTS.abs_error_tolerance);
        let base_convergence_test = initial_test_arg.sqrt();
        let mut convergence_test = base_convergence_test;
        fwd_k_minus_1 /= abs_fwd_k_minus_1;
        fwd_k /= abs_fwd_k_minus_1;
        abs_fwd_k /= abs_fwd_k_minus_1;
        let mut rough_check = true;

        // we expect to break before the end (i.e. never get to i == 1000)
        // in fortran this was an infinite loop, but here I want the loop index
        for i in 1..1000 {
            // first loop roughly checking that we are in a high-growth region
            n_steps += 1;
            abs_fwd_k_minus_1 = abs_fwd_k;
            let recurrence_factor = two_over_z * T::from_isize(starting_index + i);
            (fwd_k_minus_1, fwd_k) = (fwd_k, fwd_k_minus_1 - (recurrence_factor * fwd_k));

            abs_fwd_k = fwd_k.abs();
            if abs_fwd_k_minus_1 <= convergence_test {
                continue;
            }
            // if we get here, we have reached the high growth region, and move into
            // doing a more refined check. Note that the convergence_test is modified below
            // if we then reach this point with the modified convergence_test, then break
            if !rough_check {
                break;
            }
            rough_check = false;

            let abs_next_recurrence_factor = (recurrence_factor + two_over_z).abs() / T::two();
            let lambda =
                abs_next_recurrence_factor + (abs_next_recurrence_factor.powi(2) - T::one()).sqrt();
            let rho = abs_fwd_k / abs_fwd_k_minus_1.min(lambda);
            convergence_test = base_convergence_test * (rho / (rho.powi(2) - T::one())).sqrt();
        }
    }

    let mut val_k = Complex::<T>::new(T::one() / abs_fwd_k, T::ZERO);
    let mut val_k_plus_1 = T::C_ZERO;

    {
        // Phase 2: Calculate the unscaled top ratio
        // We run a standard Miller backward recurrence starting from the truncation bound,
        // without bothering to calculate the Neumann normalisation sum.
        // We don't care about absolute values, only the ratio of the top two terms.
        let n_backward_steps = n_steps + 1 - index_difference;
        let modified_order = order + T::from_usize(n - 1);
        for k in (1..=n_backward_steps as usize).rev() {
            (val_k, val_k_plus_1) = (
                val_k * (two_over_z * (modified_order + T::from_usize(k))) + val_k_plus_1,
                val_k,
            );
        }
        if val_k.re == T::ZERO && val_k.im == T::ZERO {
            val_k = Complex::<T>::new(
                T::MACHINE_CONSTANTS.abs_error_tolerance,
                T::MACHINE_CONSTANTS.abs_error_tolerance,
            );
        }
    }
    let mut ratios = T::c_zeros(n);
    ratios[n - 1] = val_k_plus_1 / val_k;
    if n > 1 {
        // Phase 3: Evaluate the continued fraction downwards
        // Since R_{k-1} = 1 / (2(\nu+k)/z + R_k), we can simply step downwards
        // using the anchored top ratio to evaluate the continued fraction for the entire array.
        let base_order_term = order * two_over_z;
        for k in (1..n).rev() {
            let mut fraction_denominator =
                base_order_term + T::from_usize(k) * two_over_z + ratios[k];
            let mut abs_pt = fraction_denominator.abs();
            if abs_pt == T::ZERO {
                fraction_denominator = Complex::<T>::new(
                    T::MACHINE_CONSTANTS.abs_error_tolerance,
                    T::MACHINE_CONSTANTS.abs_error_tolerance,
                );
                abs_pt = fraction_denominator.abs();
            }
            ratios[k - 1] = fraction_denominator.conj() / abs_pt.powi(2);
        }
    }
    ratios
}

/// Iterate through k functions (the first number of which may be zeros), set
/// them to zero on underflow, continuing recurrence
/// on scaled functions until the first two members of the sequence *don't* underflow
/// come on scale, then return with min(n_zeros+2,n) values scaled by 1/tol.
///
/// Originally ZKSCL
pub(crate) fn scale_k_recurrence<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    n: usize,
    // input y values are scaled
    y: &mut [Complex<T>],
    n_zeros: &mut usize,
    two_over_z: Complex<T>,
) {
    *n_zeros = 0;
    let mut i_completed = 0;

    // Copy the values by value before we start mutating y
    let original_scaled_0 = y[0];
    let original_scaled_1 = if n > 1 { y[1] } else { T::C_ZERO };

    // repeats twice, unless n < 2
    // This tests the first two values in the loop, and can be the only bit
    // that runs
    for (i, yi) in y.iter_mut().enumerate().take(min(2, n)) {
        let current_val = *yi;
        // Assumption: the value is too small (will underflow)
        *n_zeros += 1;
        *yi = T::C_ZERO;
        if -z.re + current_val.abs().ln() < -T::MACHINE_CONSTANTS.exponent_limit {
            // if the scaling would put the (negative) exponent below the (negative)
            // limit, the the value was too small (assumption true)
            continue;
        }
        // note: this is unscaled by the standard scaling, but is still a factor of
        // abs_error_tolerance smaller than the final answer
        let unscaled_value =
            (current_val.ln() - z).exp() / T::MACHINE_CONSTANTS.abs_error_tolerance;
        if will_underflow(unscaled_value) {
            // if the scaled value would underflow when (later) put back on scale by multiplying
            // by abs_error_tolerance, then the value is still too small, and would underflow
            // so assumption above is true
            continue;
        }
        // Here we know the assumption is false, so set the value properly and
        // decrement n_zeros to undo the increment above
        *yi = unscaled_value;
        i_completed = i;
        *n_zeros -= 1;
    }
    if n <= 2 || *n_zeros == 0 {
        // If there are less than two values requested, we've tested them all, so also
        // return.
        // n_zeros == 0 means that both the first two value were on scale, and
        // we can return.
        return;
    }

    let mut scaled_k_minus_1 = original_scaled_0;
    let mut scaled_k = original_scaled_1;
    let half_exponent_limit = T::half() * T::MACHINE_CONSTANTS.exponent_limit;
    let internal_scaling_factor = (-T::MACHINE_CONSTANTS.exponent_limit).exp();
    let mut effective_z = z;
    // Run the recurrence forward on the scaled values until two consecutive values come
    // on scale. If the scaled sequence grows too large (exceeding e^(limit/2)), we dynamically
    // scale it down to avoid intermediate overflow.
    let mut found_two_good_values = false;
    let mut n_tested = 0;
    for (i, yi) in y.iter_mut().enumerate().skip(2) {
        n_tested = i;
        let recurrence_factor = (order + T::from_usize(i - 1)) * two_over_z;
        (scaled_k_minus_1, scaled_k) = (scaled_k, scaled_k * recurrence_factor + scaled_k_minus_1);

        // Assumption: the value is too small (will underflow)
        *n_zeros += 1;
        *yi = T::C_ZERO;
        if -effective_z.re + scaled_k.abs().ln() >= -T::MACHINE_CONSTANTS.exponent_limit {
            // note: the value below is unscaled by the standard scaling, but is still a factor of
            // abs_error_tolerance smaller than the final answer
            let unscaled_value =
                (scaled_k.ln() - effective_z).exp() / T::MACHINE_CONSTANTS.abs_error_tolerance;
            if !will_underflow(unscaled_value) {
                // by this point, we know the assumption was false, so set the value,
                // and undo the increment we did above
                *yi = unscaled_value;
                *n_zeros -= 1;

                // the if below means:
                // "If we got to this line twice in a row on two iterations of the loop"
                if i_completed == i - 1 {
                    found_two_good_values = true;
                    break;
                }
                i_completed = i;
                continue;
            }
        }

        if scaled_k.abs().ln() > half_exponent_limit {
            effective_z -= T::MACHINE_CONSTANTS.exponent_limit;
            scaled_k_minus_1 *= internal_scaling_factor;
            scaled_k *= internal_scaling_factor;
        }
    }
    if found_two_good_values {
        *n_zeros = n_tested - 2;
    } else {
        *n_zeros = n;
        if i_completed == n {
            // this means we found one good value, on the last iteration
            *n_zeros = n - 1
        }
    }
}

/// Runs a standard Bessel recurrence relation while applying dynamic scaling
/// to prevent floating-point overflow.
///
/// This function executes the standard recurrence `(s1, s2) = (s2, s1 + ck * s2)`
/// across a slice of the output array `y`. On each iteration, it writes the scaled
/// `s2` into the array and updates the provided `OverflowState` to ensure intermediate
/// arithmetic does not exceed exponent limits.
///
/// # Arguments
/// * `forward` - If `true`, iterates forwards (increasing order) and writes to `y[n_offset..]`.
///   If `false`, iterates backwards (decreasing order) and writes to `y[..n_offset]`.
/// * `order` - The base order $\nu$ corresponding to `y[0]`.
/// * `z` - The complex argument $z$.
/// * `y` - The output array where the computed values will be written.
/// * `n_offset` - The index partitioning the array. Defines where the iteration begins
///   depending on the `forward` direction.
/// * `s1`, `s2` - The starting values for the recurrence. For forward recurrence, these
///   typically correspond to $Z_{\nu-1}$ and $Z_{\nu}$ (or related scaled terms).
/// * `overflow_state` - The state tracker used to manage dynamic reciprocal scaling factors.
///
/// # Mathematical Details
/// The recurrence multiplier `recurrence_factor` is computed dynamically from the absolute array index `i`,
/// rendering the loop stateless. Depending on the `forward` flag, `ck` evaluates to exactly
/// $\frac{2}{z}(\nu + i \pm 1)$, correctly mirroring the specific step in the sequence.
#[allow(clippy::too_many_arguments)]
#[inline]
pub(crate) fn scale_controlled_recurrence<T: BesselFloat>(
    forward: bool,
    order: T,
    z: Complex<T>,
    y: &mut [Complex<T>],
    n_offset: usize,
    mut s1: Complex<T>,
    mut s2: Complex<T>,
    mut overflow_state: OverflowState,
) {
    let two_over_z = two_over_z_safe(z);

    let base_iterator = y.iter_mut().enumerate();
    let iterator = if forward {
        Either::Right(base_iterator.skip(n_offset))
    } else {
        Either::Left(base_iterator.take(n_offset).rev())
    };
    let index_adjustment = if forward { -T::one() } else { T::one() };

    let mut recip_scale_factor = overflow_state.reciprocal_scaling_factor::<T>();
    let mut boundary = overflow_state.boundary::<T>();

    for (i, yi) in iterator {
        let recurrence_factor = two_over_z * (order + T::from_usize(i) + index_adjustment);
        (s1, s2) = (s2, s1 + recurrence_factor * s2);
        *yi = s2 * recip_scale_factor;
        overflow_state.scale_recurrence(
            &mut s1,
            &mut s2,
            *yi,
            &mut boundary,
            &mut recip_scale_factor,
        );
    }
}
