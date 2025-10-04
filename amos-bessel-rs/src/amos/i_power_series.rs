use num::{
    complex::{Complex64, ComplexFloat},
    traits::Pow,
};

use super::{
    BesselResult, BesselValues, MACHINE_CONSTANTS, Scaling, c_one, c_zero, c_zeros, gamma_ln,
    utils::will_underflow,
};

/// z_power_series computes the i bessel function for `real(z) >= 0.0`` by
/// means of the power series for large `z.abs()` in the
/// region `z.abs() <= 2*sqrt(fnu+1)`. nz=0 is a normal return.
/// nz > 0 means that the last nz components were set to zero
/// due to underflow. nz < 0 means underflow occurred, but the
/// condition cabs(z) <= 2*sqrt(fnu+1) was violated and the
/// computation must be completed in another routine with n=n-abs(nz).
///
/// Originally ZSERI
pub fn i_power_series(
    z: Complex64,
    order: f64,
    kode: Scaling,
    n: usize,
) -> BesselResult<BesselValues<isize>> {
    let mut nz = 0;
    let abs_z = z.abs();
    let mut y = c_zeros(n);

    if abs_z < MACHINE_CONSTANTS.underflow_limit {
        // If z is zero or very small, can return straight away.
        // If its zero, then nz = 0 (as y==0), but if its very small but nonzero, then
        // we underflowed, so set nz = n. This is then adjusted for order = 0,
        // as we can set y[0] to one, and return one less nz.
        if order == 0.0 {
            y[0] = c_one();
        }
        if abs_z != 0.0 {
            nz = n.try_into().unwrap();
            if order == 0.0 {
                nz -= 1;
            }
        }
        return Ok((y, nz));
    }

    let mut scale_factor = 1.0;
    let mut near_underflow = false;
    let mut first_entries_scaled: Vec<Complex64> = vec![];
    let half_z = 0.5 * z;
    let cz = if abs_z > MACHINE_CONSTANTS.underflow_limit.sqrt() {
        half_z.pow(2.0)
    } else {
        c_zero()
    };
    let rz = 2.0 * z.conj() / (abs_z.pow(2));

    let abs_cz = cz.abs();
    let ln_half_z = half_z.ln();

    let [mut y_k_plus_2, mut y_k_plus_1] = [c_zero(); 2];
    for k in (0..n).rev() {
        if first_entries_scaled.len() < 2 {
            let modified_order = order + (k as f64);

            // UNDERFLOW TEST
            // Recur down (setting y to zero) from N until underflow no longer found,
            // then move on to more set last two elements (though still being careful of
            // potential underflow)
            let mut ak1 = ln_half_z * modified_order;
            ak1.re -= gamma_ln(modified_order + 1.0).unwrap();
            if kode == Scaling::Scaled {
                ak1.re -= z.re;
            }
            if ak1.re <= -MACHINE_CONSTANTS.exponent_limit {
                nz += 1;
                y[k] = c_zero();
                if abs_cz > modified_order {
                    break;
                }
                continue;
            }

            // Now do a more refined underflow test.
            // Note that near_undeflow latches: it does not reset to false on
            // a second pass through this block, only later is it explicitly reset
            if ak1.re <= (-MACHINE_CONSTANTS.approximation_limit) {
                near_underflow = true;
                scale_factor = MACHINE_CONSTANTS.abs_error_tolerance;
            }

            let mut coeff = ak1.exp();
            if near_underflow {
                coeff *= MACHINE_CONSTANTS.rtol
            };
            let s1 = single_n_iteration(modified_order, cz);
            let s2 = s1 * coeff;
            first_entries_scaled.push(s2);
            if near_underflow
                && will_underflow(
                    s2,
                    MACHINE_CONSTANTS.absolute_approximation_limit,
                    MACHINE_CONSTANTS.abs_error_tolerance,
                )
            {
                nz += 1;
                y[k] = c_zero();
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
                first_entries_scaled.push(c_zero());
            }

            let modified_order = ((k + 1) as f64) + order;
            if near_underflow {
                // ... using scaled values
                (y_k_plus_2, y_k_plus_1) =
                    (y_k_plus_1, modified_order * (rz * y_k_plus_1) + y_k_plus_2);
                y[k] = y_k_plus_1 * scale_factor;
                if y[k].abs() > MACHINE_CONSTANTS.absolute_approximation_limit {
                    near_underflow = false;
                }
            } else {
                // .. using unscaled values
                y[k] = modified_order * (rz * y[k + 1]) + y[k + 2];
            }
        }
    }
    Ok((y, nz))
}

fn single_n_iteration(modified_order: f64, cz: Complex64) -> Complex64 {
    let fnup = modified_order + 1.0;
    let abs_cz = cz.abs();
    let atol = MACHINE_CONSTANTS.abs_error_tolerance * abs_cz / fnup;

    let fnup = modified_order + 1.0;
    let mut s1 = c_one();
    if abs_cz >= MACHINE_CONSTANTS.abs_error_tolerance * fnup {
        let mut ak2 = c_one();
        let mut ak = fnup + 2.0;
        let mut s = fnup;
        let mut aa = 2.0;
        while aa > atol {
            ak2 *= cz / s;
            s1 += ak2;
            s += ak;
            ak += 2.0;
            aa *= abs_cz / s;
        }
    }
    s1
}
