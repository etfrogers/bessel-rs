use super::Scaling;
use crate::{
    BesselError::{self, DidNotConverge},
    BesselFloat,
    amos::{
        IKType, MachineConsts, gamma_ln, limits::check_underflow_uniform_asymp_params,
        right_half_plane::k_right_half_plane, utils::two_over_z_safe,
    },
};

use num::{Complex, complex::ComplexFloat};

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
    out: &mut [Complex<T>],
) -> Result<(), BesselError<T>> {
    let mc: &MachineConsts<T> = T::MACHINE_CONSTANTS;
    let n = out.len();
    let scale: T = T::TWO * T::MIN_POSITIVE / mc.abs_error_tolerance;
    let abs_z = z.abs();
    let int_abs_z = abs_z.to_usize().unwrap();
    let int_order = order.to_usize().unwrap();
    let modified_int_order = int_order + n - 1;
    let abs_z_plus_one = T::from_usize(int_abs_z + 1);
    let reciprocal_abs_z = T::ONE / abs_z;
    let two_over_z = two_over_z_safe(z);
    let mut fwd_k_minus_1 = T::C_ZERO;
    let mut fwd_k = T::C_ONE;
    let mut abs_recurrence_factor = (abs_z_plus_one + T::ONE) * reciprocal_abs_z;
    let rho =
        abs_recurrence_factor + (abs_recurrence_factor * abs_recurrence_factor - T::ONE).sqrt();
    let rho_sq = rho * rho;
    let mut convergence_test = (rho_sq + rho_sq) / ((rho_sq - T::ONE) * (rho - T::ONE));
    convergence_test /= mc.abs_error_tolerance;
    // Phase 1: Forward Sequence Truncation Bound
    // Run the recurrence forward. The sequence diverges, representing the rapidly growing K_nu(z).
    // We run it until the sequence exceeds a convergence threshold, which tells us how high
    // an index we need to start the backward recurrence from to ensure the truncation error
    // doesn't pollute the final values at our target index.
    let mut converged = false;
    let mut series_trunctation_index = 0;
    let mut recurrence_factor = two_over_z * (abs_z_plus_one * T::HALF);
    let mut current_index_magnitude = abs_z_plus_one;
    for i in 0..80 {
        series_trunctation_index = i + 2;
        // below is the un-optimised line that is now replaced by incremental addition to reduce computational load
        // un-optimised code is retained for explanation
        // let current_index_magnitude = abs_z_plus_one + T::from_usize(i);
        // let recurrence_factor = two_over_z * ((abs_z_plus_one + T::from_usize(i * 2)) / T::TWO);
        (fwd_k_minus_1, fwd_k) = (fwd_k, fwd_k_minus_1 - recurrence_factor * fwd_k);
        let threshold = convergence_test * current_index_magnitude * current_index_magnitude;
        if fwd_k.norm_sqr() > threshold * threshold {
            converged = true;
            break;
        }
        recurrence_factor += two_over_z;
        current_index_magnitude += T::ONE;
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
        let starting_order = T::from_usize(modified_int_order + 1);
        let mut convergence_test_sqr = starting_order * reciprocal_abs_z / mc.abs_error_tolerance;
        let mut hit_loop_end = false;
        converged = false;
        let mut recurrence_factor = two_over_z * (starting_order * T::HALF);
        for k in 0..80 {
            ratio_truncation_index = k + 1;
            // below is the un-optimised line that is now replaced by incremental addition to reduce computational load
            // un-optimised code is retained for explanation
            // let recurrence_factor = two_over_z * ((starting_order + T::from_usize(i * 2)) / T::TWO);
            (fwd_k_minus_1, fwd_k) = (fwd_k, fwd_k_minus_1 - recurrence_factor * fwd_k);
            let fwd_k_sqr = fwd_k.norm_sqr();

            if fwd_k_sqr < convergence_test_sqr {
                recurrence_factor += two_over_z;
                continue;
            }
            if hit_loop_end {
                converged = true;
                break;
            }
            abs_recurrence_factor = recurrence_factor.abs();

            let lambda = abs_recurrence_factor
                + (abs_recurrence_factor * abs_recurrence_factor - T::ONE).sqrt();
            let abs_fwd_k = fwd_k_sqr.sqrt();
            let kappa = abs_fwd_k / fwd_k_minus_1.abs();
            let rho = lambda.min(kappa);
            convergence_test_sqr *= rho / (rho * rho - T::ONE);
            recurrence_factor += two_over_z;
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
    let mut binomial_coeff = (gamma_ln(kk_float + twice_fractional_order + T::ONE, mc).unwrap()
        - gamma_ln(kk_float + T::ONE, mc).unwrap()
        - gamma_ln(twice_fractional_order + T::ONE, mc).unwrap())
    .exp();
    let mut normalisation_sum = T::C_ZERO;
    // Neumann normalisation loop
    for _ in 0..(start_index - modified_int_order) {
        let pt = val_k;
        val_k = val_k_plus_one + (kk_float + fractional_order) * (two_over_z * val_k);
        val_k_plus_one = pt;
        let binomial_ratio = T::ONE - twice_fractional_order / (kk_float + twice_fractional_order);
        let next_binomial_coeff = binomial_coeff * binomial_ratio;
        normalisation_sum += (next_binomial_coeff + binomial_coeff) * val_k_plus_one;
        binomial_coeff = next_binomial_coeff;
        kk_float -= T::ONE;
    }

    out[n - 1] = val_k;
    if n != 1 {
        for i in 1..n {
            let pt = val_k;
            val_k = val_k_plus_one + (kk_float + fractional_order) * (two_over_z * pt);
            val_k_plus_one = pt;
            let binomial_ratio =
                T::ONE - twice_fractional_order / (kk_float + twice_fractional_order);
            let next_binomial_coeff = binomial_coeff * binomial_ratio;
            normalisation_sum += (next_binomial_coeff + binomial_coeff) * val_k_plus_one;
            binomial_coeff = next_binomial_coeff;
            kk_float -= T::ONE;
            out[n - (i + 1)] = val_k;
        }
    }
    if int_order > 0 {
        for _i in 0..int_order {
            (val_k_plus_one, val_k) = (
                val_k,
                val_k_plus_one + (kk_float + fractional_order) * (two_over_z * val_k),
            );
            let binomial_ratio =
                T::ONE - twice_fractional_order / (kk_float + twice_fractional_order);
            let next_binomial_coeff = binomial_coeff * binomial_ratio;
            normalisation_sum += (next_binomial_coeff + binomial_coeff) * val_k_plus_one;
            binomial_coeff = next_binomial_coeff;
            kk_float -= T::ONE;
        }
    }

    let mut scaled_z = z;
    if scaling == Scaling::Scaled {
        scaled_z.re = T::ZERO;
    }
    let mut ln_leading_term = -fractional_order * two_over_z.ln() + scaled_z;
    let gamma_term = gamma_ln(T::ONE + fractional_order, mc).unwrap();
    ln_leading_term -= gamma_term;
    // Calculate the final normalisation constant.
    // The complex division exp(ln_leading_term) / (normalisation_sum + val_k) is performed
    // by dividing by the magnitude twice, to avoid intermediate overflow from squaring
    // large quantities when computing the complex denominator.
    val_k += normalisation_sum;
    let sum_magnitude = val_k.abs();
    let normalization_constant =
        (ln_leading_term.exp() / sum_magnitude) * val_k.conj() / sum_magnitude;
    for element in out.iter_mut() {
        *element *= normalization_constant;
    }
    Ok(())
}

/// i_ratios computes ratios of I bessel functions by backward
/// recurrence. The starting index is determined by forward
/// recurrence as described in J. Res. of Nat. Bur. of Standards-B,
/// Mathematical Sciences, vol 77b, p111-114, September, 1973,
/// Bessel functions I and J of complex argument and integer order,
/// by D. J. Sookne.
///
/// Originally ZRATI
pub(crate) fn i_ratios<T: BesselFloat>(z: Complex<T>, order: T, out: &mut [Complex<T>]) {
    let mc: &MachineConsts<T> = T::MACHINE_CONSTANTS;
    let n = out.len();
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
        let abs_fwd_k_minus_1 = fwd_k_minus_1.abs();
        // Scale base_convergence_test and all subsequent fwd_k values by
        // abs_fwd_k_minus_1 to ensure that an overflow does not occur prematurely
        let initial_test_arg =
            (abs_fwd_k + abs_fwd_k) / (abs_fwd_k_minus_1 * mc.abs_error_tolerance);
        let base_convergence_test = initial_test_arg;
        let mut convergence_test = base_convergence_test;
        fwd_k_minus_1 /= abs_fwd_k_minus_1;
        fwd_k /= abs_fwd_k_minus_1;
        // abs_fwd_k /= abs_fwd_k_minus_1;
        let mut rough_check = true;

        let mut abs_fwd_k_sqr = abs_fwd_k * abs_fwd_k; //
        // we expect to break before the end (i.e. never get to i == 1000)
        // in fortran this was an infinite loop, but here I want the loop index
        for i in 1..1000 {
            // first loop roughly checking that we are in a high-growth region
            n_steps += 1;
            let abs_fwd_k_minus_1_sqr = abs_fwd_k_sqr;
            let recurrence_factor = two_over_z * T::from_isize(starting_index + i);
            (fwd_k_minus_1, fwd_k) = (fwd_k, fwd_k_minus_1 - (recurrence_factor * fwd_k));

            abs_fwd_k_sqr = fwd_k.norm_sqr();
            if abs_fwd_k_minus_1_sqr <= convergence_test {
                continue;
            }
            // if we get here, we have reached the high growth region, and move into
            // doing a more refined check. Note that the convergence_test is modified below
            // if we then reach this point with the modified convergence_test, then break
            if !rough_check {
                break;
            }
            rough_check = false;

            let abs_next_recurrence_factor = (recurrence_factor + two_over_z).abs() / T::TWO;
            let lambda =
                abs_next_recurrence_factor + (abs_next_recurrence_factor.powi(2) - T::ONE).sqrt();
            abs_fwd_k = fwd_k.abs();
            let rho = abs_fwd_k / abs_fwd_k_minus_1_sqr.sqrt().min(lambda);
            convergence_test = base_convergence_test * (rho / (rho.powi(2) - T::ONE));
        }
        abs_fwd_k = abs_fwd_k_sqr.sqrt();
    }

    let mut val_k = Complex::<T>::new(T::ONE / abs_fwd_k, T::ZERO);
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
            val_k = Complex::<T>::new(mc.abs_error_tolerance, mc.abs_error_tolerance);
        }
    }

    out[n - 1] = val_k_plus_1 / val_k;
    if n > 1 {
        // Phase 3: Evaluate the continued fraction downwards
        // Since R_{k-1} = 1 / (2(\nu+k)/z + R_k), we can simply step downwards
        // using the anchored top ratio to evaluate the continued fraction for the entire array.
        let base_order_term = order * two_over_z;
        for k in (1..n).rev() {
            let mut fraction_denominator = base_order_term + T::from_usize(k) * two_over_z + out[k];
            let mut abs_frac_denom_sqr = fraction_denominator.norm_sqr();
            if abs_frac_denom_sqr == T::ZERO {
                fraction_denominator =
                    Complex::<T>::new(mc.abs_error_tolerance, mc.abs_error_tolerance);
                abs_frac_denom_sqr = fraction_denominator.norm_sqr();
            }
            out[k - 1] = fraction_denominator.conj() / abs_frac_denom_sqr;
        }
    }
}

/// Computes the $I$ Bessel sequence for $\text{Re}(z) \ge 0$ by
/// normalizing the ratios from [i_ratios] using the Wronskian identity with $K_\nu$ and $K_{\nu+1}$.
///
/// Originally ZWRSK
pub(crate) fn i_wronskian<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    out: &mut [Complex<T>],
) -> Result<usize, BesselError<T>> {
    let n = out.len();
    let mc: &MachineConsts<T> = T::MACHINE_CONSTANTS;
    match check_underflow_uniform_asymp_params(
        z,
        order,
        scaling,
        IKType::K,
        2,
        &mut [T::C_ONE; 2],
        mc,
    ) {
        Ok(n_underflow) => {
            if n_underflow > 0 {
                return Err(BesselError::Overflow);
            }
        }
        Err(_) => {
            out.fill(T::C_ZERO);
            return Ok(n);
        }
    }

    // 1. Compute K_nu and K_{nu+1} to serve as Wronskian anchors
    let mut k_values = [T::C_ZERO; 2];
    let _ = k_right_half_plane(z, order, scaling, &mut k_values)?;
    // 2. Compute backward recurrence ratios r_{nu+j} = I_{nu+j+1} / I_{nu+j}
    let y_ratios = out;
    i_ratios(z, order, y_ratios);

    // Initial phase factor for scaled computation (e^{i * Im(z)})
    let mut current_i = if scaling == Scaling::Scaled {
        Complex::<T>::cis(z.im)
    } else {
        T::C_ONE
    };

    // On low-exponent machines, K values can be close to under/overflow limits.
    // Scale K_nu and K_{nu+1} to keep intermediate products well on scale.
    let abs_k_nu_plus_1 = k_values[1].abs();
    let k_scale_factor = if abs_k_nu_plus_1 <= mc.absolute_approximation_limit {
        T::ONE / mc.abs_error_tolerance
    } else if abs_k_nu_plus_1 >= T::ONE / mc.absolute_approximation_limit {
        mc.abs_error_tolerance
    } else {
        T::ONE
    };
    let scaled_k_nu = k_values[0] * k_scale_factor;
    let scaled_k_nu_plus_1 = k_values[1] * k_scale_factor;

    // Evaluate Wronskian denominator: denom = z * (r_nu * K_nu + K_{nu+1}) = 1 / I_nu
    // Performing division as (current_i / |denom|) * (conj(denom) / |denom|) avoids
    // squaring |denom| which could overflow/underflow prematurely.
    let mut wronskian_denom = z * (y_ratios[0] * scaled_k_nu + scaled_k_nu_plus_1);
    let abs_denom = wronskian_denom.abs();
    wronskian_denom = wronskian_denom.conj() / abs_denom;
    current_i = (current_i / abs_denom) * wronskian_denom;

    let out = y_ratios;
    // Multiply by k_scale_factor to restore true scale
    let mut next_ratio = out[0]; // Save r_0 before overwriting out[0]
    out[0] = current_i * k_scale_factor;

    // Step forward: I_{nu+j} = I_{nu+j-1} * r_{nu+j-1}
    for out_i in out.iter_mut().skip(1) {
        let ratio = next_ratio;
        next_ratio = *out_i; // Read ratio r_{i} before overwriting out[i]
        current_i *= ratio;
        *out_i = current_i * k_scale_factor;
    }
    Ok(0)
}
