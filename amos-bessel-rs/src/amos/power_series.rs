use num::complex::{Complex, ComplexFloat};

use crate::{
    BesselError, BesselFloat, Scaling,
    amos::{
        MachineConsts, gamma_ln,
        utils::{two_over_z_safe, will_underflow},
    },
};

/// z_power_series computes the I bessel function for `real(z) >= 0.0` by
/// means of the power series for large `z.abs()` in the
/// region `z.abs() <= 2*sqrt(fnu+1)`. n_zeros=0 is a normal return.
/// n_zeros > 0 means that the last n_zeros components were set to zero
/// due to underflow. n_zeros < 0 means underflow occurred, but the
/// condition cabs(z) <= 2*sqrt(fnu+1) was violated and the
/// computation must be completed in another routine with n=n-abs(n_zeros).
///
/// Originally ZSERI
pub fn i_power_series<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    out: &mut [Complex<T>],
) -> Result<isize, BesselError<T>> {
    let mc: &MachineConsts<T> = T::MACHINE_CONSTANTS;
    let mut n_zeros = 0;
    let abs_z = z.abs();
    let n = out.len();

    if abs_z < mc.underflow_limit {
        // If z is zero or very small, can return straight away.
        // If it's zero, then n_zeros = 0 (as y==0), but if its very small but non_zero, then
        // we underflowed, so set n_zeros = n. This is then adjusted for order = 0,
        // as we can set y[0] to one, and return one less n_zeros.
        if order == T::ZERO {
            out[0] = T::C_ONE;
        }
        if abs_z != T::ZERO {
            n_zeros = n.try_into().unwrap();
            if order == T::ZERO {
                n_zeros -= 1;
            }
        }
        return Ok(n_zeros);
    }

    let mut scale_factor = T::one();
    let mut near_underflow = false;
    let mut num_seeded = 0;
    let half_z = z * T::HALF;
    let half_z_sq = if abs_z > mc.underflow_limit.sqrt() {
        half_z.powi(2)
    } else {
        T::C_ZERO
    };
    let two_over_z = two_over_z_safe(z);

    let half_abs_z = abs_z * T::HALF;
    let abs_half_z_sq = half_abs_z * half_abs_z;
    let ln_half_z = half_z.ln();

    let [mut y_k_plus_2, mut y_k_plus_1] = [T::C_ZERO; 2];
    for k in (0..n).rev() {
        if num_seeded < 2 {
            let current_order = order + T::from_usize(k);

            // UNDERFLOW TEST
            // Recur down (setting y to zero) from N until underflow no longer found,
            // then move on to more set last two elements (though still being careful of
            // potential underflow)
            let mut ln_leading_term = ln_half_z * current_order;
            ln_leading_term.re -= gamma_ln(current_order + T::one(), mc).unwrap();
            if scaling == Scaling::Scaled {
                ln_leading_term.re -= z.re;
            }
            if ln_leading_term.re <= -mc.exponent_limit {
                n_zeros += 1;
                out[k] = T::C_ZERO;
                if abs_half_z_sq > current_order {
                    break;
                }
                continue;
            }

            // Now do a more refined underflow test.
            // Note that near_undeflow latches: it does not reset to false on
            // a second pass through this block, only later is it explicitly reset
            if ln_leading_term.re <= (-mc.approximation_limit) {
                near_underflow = true;
                scale_factor = mc.abs_error_tolerance;
            }

            let mut coeff = ln_leading_term.exp();
            if near_underflow {
                coeff *= mc.rtol
            };
            let s1 = single_n_iteration(current_order, half_z_sq, abs_half_z_sq, mc);
            let s2 = s1 * coeff;
            if near_underflow && will_underflow(s2, mc) {
                n_zeros += 1;
                out[k] = T::C_ZERO;
                continue;
            }
            if num_seeded == 0 {
                y_k_plus_2 = s2;
            } else {
                y_k_plus_1 = s2;
            }
            num_seeded += 1;
            out[k] = s2 * scale_factor;
        } else {
            // Continue recurring backward. If underflow was close previously, use scaled values,
            // but the first time that we get out of the underflow region, we can switch
            // to using the unscaled values
            let modified_order = T::from_usize(k + 1) + order;
            let next_val = (two_over_z * y_k_plus_1) * modified_order + y_k_plus_2;
            y_k_plus_2 = y_k_plus_1;
            y_k_plus_1 = next_val;
            if near_underflow {
                out[k] = y_k_plus_1 * scale_factor;
                if out[k].abs() > mc.absolute_approximation_limit {
                    near_underflow = false;
                    // Switch to unscaled values in registers
                    y_k_plus_1 = out[k];
                    y_k_plus_2 = out[k + 1];
                }
            } else {
                out[k] = y_k_plus_1;
            }
        }
    }
    Ok(n_zeros)
}

#[inline]
fn single_n_iteration<T: BesselFloat>(
    current_order: T,
    half_z_sq: Complex<T>,
    abs_half_z_sq: T,
    mc: &MachineConsts<T>,
) -> Complex<T> {
    let order_plus_one = current_order + T::one();
    let tolerance = mc.abs_error_tolerance * abs_half_z_sq / order_plus_one;

    if abs_half_z_sq < mc.abs_error_tolerance * order_plus_one {
        return T::C_ONE;
    }
    // I tried adding zero-division fast paths for the most common orders 0 and 1,
    // but overall performance decreased. The general code below was faster
    let mut series_sum = T::C_ONE;
    let mut current_term = T::C_ONE;
    let mut denominator_step = order_plus_one + T::TWO;
    let mut term_denominator = order_plus_one;
    let mut error_estimate = T::TWO;
    while error_estimate > tolerance {
        let recip_denom = T::ONE / term_denominator;
        current_term *= half_z_sq * recip_denom;
        series_sum += current_term;
        term_denominator += denominator_step;
        denominator_step += T::TWO;
        error_estimate *= abs_half_z_sq * recip_denom;
    }

    series_sum
}
