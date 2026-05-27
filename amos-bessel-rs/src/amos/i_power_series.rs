use num::complex::{Complex, ComplexFloat};

use crate::types::{BesselResult, BesselValues};
use crate::{amos::utils::calc_rz, types::BesselFloat};

use super::{Scaling, gamma_ln, utils::will_underflow};

/// z_power_series computes the I bessel function for `real(z) >= 0.0` by
/// means of the power series for large `z.abs()` in the
/// region `z.abs() <= 2*sqrt(fnu+1)`. nz=0 is a normal return.
/// nz > 0 means that the last nz components were set to zero
/// due to underflow. nz < 0 means underflow occurred, but the
/// condition cabs(z) <= 2*sqrt(fnu+1) was violated and the
/// computation must be completed in another routine with n=n-abs(nz).
///
/// Originally ZSERI
pub fn i_power_series<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    kode: Scaling,
    n: usize,
) -> BesselResult<BesselValues<T, isize>> {
    let mut nz = 0;
    let abs_z = z.abs();
    let mut y = T::c_zeros(n);

    if abs_z < T::MACHINE_CONSTANTS.underflow_limit {
        // If z is zero or very small, can return straight away.
        // If its zero, then nz = 0 (as y==0), but if its very small but nonzero, then
        // we underflowed, so set nz = n. This is then adjusted for order = 0,
        // as we can set y[0] to one, and return one less nz.
        if order == T::zero() {
            y[0] = T::C_ONE;
        }
        if abs_z != T::zero() {
            nz = n.try_into().unwrap();
            if order == T::zero() {
                nz -= 1;
            }
        }
        return Ok((y, nz));
    }

    let mut scale_factor = T::one();
    let mut near_underflow = false;
    let mut first_entries_scaled: Vec<Complex<T>> = vec![];
    let half_z = z * T::half();
    let cz = if abs_z > T::MACHINE_CONSTANTS.underflow_limit.sqrt() {
        half_z.powi(2)
    } else {
        T::C_ZERO
    };
    let rz = calc_rz(z);

    let abs_cz = cz.abs();
    let ln_half_z = half_z.ln();

    let [mut y_k_plus_2, mut y_k_plus_1] = [T::C_ZERO; 2];
    for k in (0..n).rev() {
        if first_entries_scaled.len() < 2 {
            let modified_order = order + T::from_f64(k as f64);

            // UNDERFLOW TEST
            // Recur down (setting y to zero) from N until underflow no longer found,
            // then move on to more set last two elements (though still being careful of
            // potential underflow)
            let mut ak1 = ln_half_z * modified_order;
            ak1.re -= gamma_ln(modified_order + T::one()).unwrap();
            if kode == Scaling::Scaled {
                ak1.re -= z.re;
            }
            if ak1.re <= -T::MACHINE_CONSTANTS.exponent_limit {
                nz += 1;
                y[k] = T::C_ZERO;
                if abs_cz > modified_order {
                    break;
                }
                continue;
            }

            // Now do a more refined underflow test.
            // Note that near_undeflow latches: it does not reset to false on
            // a second pass through this block, only later is it explicitly reset
            if ak1.re <= (-T::MACHINE_CONSTANTS.approximation_limit) {
                near_underflow = true;
                scale_factor = T::MACHINE_CONSTANTS.abs_error_tolerance;
            }

            let mut coeff = ak1.exp();
            if near_underflow {
                coeff *= T::MACHINE_CONSTANTS.rtol
            };
            let s1 = single_n_iteration(modified_order, cz);
            let s2 = s1 * coeff;
            first_entries_scaled.push(s2);
            if near_underflow
                && will_underflow(
                    s2,
                    T::MACHINE_CONSTANTS.absolute_approximation_limit,
                    T::MACHINE_CONSTANTS.abs_error_tolerance,
                )
            {
                nz += 1;
                y[k] = T::C_ZERO;
                continue;
            }
            y[k] = s2 * scale_factor;
        } else {
            // Continue recurring backward. If underflow was close previously, use scaled values,
            // but the first time that we get out of the underflow region, we can switch
            // to using the unscaled values
            if first_entries_scaled.len() == 2 {
                y_k_plus_2 = first_entries_scaled[0];
                y_k_plus_1 = first_entries_scaled[1];
                // the line below makes len == 3, so this block is not called again.
                first_entries_scaled.push(T::C_ZERO);
            }

            let modified_order = T::from_f64((k + 1) as f64) + order;
            if near_underflow {
                // ... using scaled values
                (y_k_plus_2, y_k_plus_1) =
                    (y_k_plus_1, (rz * y_k_plus_1) * modified_order + y_k_plus_2);
                y[k] = y_k_plus_1 * scale_factor;
                if y[k].abs() > T::MACHINE_CONSTANTS.absolute_approximation_limit {
                    near_underflow = false;
                }
            } else {
                // .. using unscaled values
                y[k] = (rz * y[k + 1]) * modified_order + y[k + 2];
            }
        }
    }
    Ok((y, nz))
}

fn single_n_iteration<T: BesselFloat>(modified_order: T, cz: Complex<T>) -> Complex<T> {
    let fnup = modified_order + T::one();
    let abs_cz = cz.abs();
    let atol = T::MACHINE_CONSTANTS.abs_error_tolerance * abs_cz / fnup;

    let mut s1 = T::C_ONE;
    if abs_cz >= T::MACHINE_CONSTANTS.abs_error_tolerance * fnup {
        let mut ak2 = T::C_ONE;
        let mut ak = fnup + T::two();
        let mut s = fnup;
        let mut aa = T::two();
        while aa > atol {
            ak2 *= cz / s;
            s1 += ak2;
            s += ak;
            ak += T::two();
            aa *= abs_cz / s;
        }
    }
    s1
}
