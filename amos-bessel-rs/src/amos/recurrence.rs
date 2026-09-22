use crate::{
    amos::{
        MachineConsts,
        limits::OverflowState,
        utils::{two_over_z_safe, will_underflow},
    },
    types::BesselFloat,
};
use core::cmp::min;
use num::Complex;

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
    let mc: &MachineConsts<T> = T::MACHINE_CONSTANTS;
    *n_zeros = 0;
    let mut prev_on_scale = false;

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
        if -z.re + (T::HALF * current_val.norm_sqr().ln()) < -mc.exponent_limit {
            // if the scaling would put the (negative) exponent below the (negative)
            // limit, the the value was too small (assumption true)
            continue;
        }
        // note: this is unscaled by the standard scaling, but is still a factor of
        // abs_error_tolerance smaller than the final answer
        let unscaled_value = (current_val.ln() - z).exp() / mc.abs_error_tolerance;
        if will_underflow(unscaled_value, mc) {
            // if the scaled value would underflow when (later) put back on scale by multiplying
            // by abs_error_tolerance, then the value is still too small, and would underflow
            // so assumption above is true
            continue;
        }
        // Here we know the assumption is false, so set the value properly and
        // decrement n_zeros to undo the increment above
        *yi = unscaled_value;
        if i == 1 {
            prev_on_scale = true;
        }
        *n_zeros -= 1;
    }
    if n == 1 {
        return;
    }
    if !prev_on_scale {
        y[0] = T::C_ZERO;
        *n_zeros = 2;
    }
    if n == 2 || *n_zeros == 0 {
        return;
    }

    let mut scaled_k_minus_1 = original_scaled_0;
    let mut scaled_k = original_scaled_1;
    let half_exponent_limit = T::HALF * mc.exponent_limit;
    let internal_scaling_factor = (-mc.exponent_limit).exp();
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

        let ln_abs_k = T::HALF * scaled_k.norm_sqr().ln();

        // Assumption: the value is too small (will underflow)
        *n_zeros += 1;
        *yi = T::C_ZERO;
        if -effective_z.re + ln_abs_k >= -mc.exponent_limit {
            // note: the value below is unscaled by the standard scaling, but is still a factor of
            // abs_error_tolerance smaller than the final answer
            let unscaled_value = (scaled_k.ln() - effective_z).exp() / mc.abs_error_tolerance;
            if !will_underflow(unscaled_value, mc) {
                // by this point, we know the assumption was false, so set the value,
                // and undo the increment we did above
                *yi = unscaled_value;
                *n_zeros -= 1;

                // the if below means:
                // "If we got to this line twice in a row on two iterations of the loop"
                if prev_on_scale {
                    found_two_good_values = true;
                    break;
                }
                prev_on_scale = true;
                continue;
            }
        }

        prev_on_scale = false;
        if ln_abs_k > half_exponent_limit {
            effective_z -= mc.exponent_limit;
            scaled_k_minus_1 *= internal_scaling_factor;
            scaled_k *= internal_scaling_factor;
        }
    }
    if found_two_good_values {
        *n_zeros = n_tested - 1;
    } else {
        *n_zeros = n;
        if prev_on_scale {
            // this means we found one good value, on the last iteration
            *n_zeros = n - 1;
        }
    }
    for yi in &mut y[..*n_zeros] {
        *yi = T::C_ZERO;
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
#[inline]
#[allow(clippy::too_many_arguments)]
pub(crate) fn scale_controlled_recurrence<T: BesselFloat>(
    forward: bool,
    order: T,
    z: Complex<T>,
    y: Option<&mut [Complex<T>]>,
    n_offset: usize,
    n: usize,
    s1: Complex<T>,
    s2: Complex<T>,
    overflow_state: OverflowState,
    mc: &MachineConsts<T>,
) -> (Complex<T>, Complex<T>, OverflowState) {
    let two_over_z = two_over_z_safe(z);

    if forward {
        let base_order = order - T::ONE;
        let iter = n_offset..n;
        match y {
            Some(out) => run_recurrence_core(
                iter,
                base_order,
                two_over_z,
                |i, yi| out[i] = yi,
                s1,
                s2,
                overflow_state,
                mc,
            ),
            None => run_recurrence_core(
                iter,
                base_order,
                two_over_z,
                |_, _| (),
                s1,
                s2,
                overflow_state,
                mc,
            ),
        }
    } else {
        let base_order = order + T::ONE;
        let iter = (0..n_offset).rev();
        match y {
            Some(out) => run_recurrence_core(
                iter,
                base_order,
                two_over_z,
                |i, yi| out[i] = yi,
                s1,
                s2,
                overflow_state,
                mc,
            ),
            None => run_recurrence_core(
                iter,
                base_order,
                two_over_z,
                |_, _| (),
                s1,
                s2,
                overflow_state,
                mc,
            ),
        }
    }
}

#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn run_recurrence_core<T: BesselFloat, I: Iterator<Item = usize>, F: FnMut(usize, Complex<T>)>(
    iter: I,
    base_order: T,
    two_over_z: Complex<T>,
    mut on_yield: F,
    mut s1: Complex<T>,
    mut s2: Complex<T>,
    mut overflow_state: OverflowState,
    mc: &MachineConsts<T>,
) -> (Complex<T>, Complex<T>, OverflowState) {
    let mut recip_scale_factor = overflow_state.reciprocal_scaling_factor::<T>(mc);
    let mut boundary = overflow_state.boundary::<T>(mc);

    for i in iter {
        let recurrence_factor = two_over_z * (base_order + T::from_usize(i));
        (s1, s2) = (s2, s1 + recurrence_factor * s2);
        let yi = s2 * recip_scale_factor;
        on_yield(i, yi);
        overflow_state.scale_recurrence(
            &mut s1,
            &mut s2,
            yi,
            &mut boundary,
            &mut recip_scale_factor,
            mc,
        );
    }
    (s1, s2, overflow_state)
}
