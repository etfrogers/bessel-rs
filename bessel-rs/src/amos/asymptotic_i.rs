use std::f64::consts::PI;

use num::{
    Integer,
    complex::{Complex64, ComplexFloat},
};

use crate::amos::utils::calc_rz;

use super::{
    BesselError::*, BesselResult, MACHINE_CONSTANTS, Scaling, c_one, c_zero, c_zeros, utils::RTPI,
};

/// asymptotic_i computes the I bessel function for real(z) >= 0.0 by
/// means of the asymptotic expansion for large z.abs() in the
/// region z.abs() > rl.max(order.pow(2.0)/2).
///
/// Originally ZASYI
pub fn asymptotic_i(z: Complex64, order: f64, scaling: Scaling, n: usize) -> BesselResult {
    let nz = 0;
    let mut y = c_zeros(n);
    let abs_z = z.abs();
    //-----------------------------------------------------------------------;
    //     OVERFLOW TEST;
    //-----------------------------------------------------------------------;
    let recip_abs_z = 1.0 / abs_z;

    let cz = match scaling {
        Scaling::Unscaled => z,
        Scaling::Scaled => Complex64::new(0.0, z.im),
    };

    if cz.re.abs() > MACHINE_CONSTANTS.exponent_limit {
        return Err(Overflow);
    }
    let scaled_calculations = cz.re.abs() > MACHINE_CONSTANTS.approximation_limit;
    let mut coeff = (RTPI * z.conj() * recip_abs_z.powi(2)).sqrt();
    if !scaled_calculations {
        coeff *= cz.exp();
    }
    let ez = z * 8.0;
    //-----------------------------------------------------------------------;
    // When z is imaginary, the error test must be made relative to the
    // first reciprocal power since this is the leading term of the
    // expansion for the imaginary part.
    //-----------------------------------------------------------------------;
    let abs_ez = 8.0 * abs_z;
    let s = MACHINE_CONSTANTS.abs_error_tolerance / abs_ez;
    let max_iterations = (MACHINE_CONSTANTS.asymptotic_z_limit * 2.0) as i64 + 2;
    let mut p1 = if z.im == 0.0 {
        c_zero()
    } else {
        //-----------------------------------------------------------------------;
        // Calculate (pi*(0.5+fnu+n-il)*i).exp() to minimize losses of;
        // significance when fnu or n is large;
        //-----------------------------------------------------------------------;
        let arg = order.fract() * PI;
        let ak = -arg.sin();
        let mut bk = arg.cos();
        if z.im < 0.0 {
            bk = -bk;
        };
        let p1 = Complex64::new(ak, bk);
        if (order as usize + n).is_even() {
            -p1
        } else {
            p1
        }
    };
    for (k, elem) in y.iter_mut().enumerate().rev().take(2.min(n)) {
        let (mut s1, s2) = {
            // this block is just to contain the large number of mutable variables in a small space
            let modified_order = order + (k as f64);
            let modified_order_sqr = (2.0 * modified_order).powf(2.0);
            let mut sqk = modified_order_sqr - 1.0;
            let atol = s * sqk.abs();
            let mut sign = 1.0;
            let mut cs1 = c_one();
            let mut cs2 = c_one();
            let mut ck = c_one();
            let mut ak = 0.0;
            let mut aa = 1.0;
            let mut bb = abs_ez;
            let mut dk = ez;
            let mut converged = false;
            'convergence: for _ in 0..max_iterations {
                ck *= sqk / dk;
                cs2 += ck;
                sign = -sign;
                cs1 += ck * sign;
                dk += ez;
                aa *= sqk.abs() / bb;
                bb += abs_ez;
                ak += 8.0;
                sqk -= ak;
                if aa <= atol {
                    converged = true;
                    break 'convergence;
                }
            }
            if !converged {
                return Err(DidNotConverge);
            }
            (cs1, cs2)
        };
        if z.re * 2.0 < MACHINE_CONSTANTS.exponent_limit {
            s1 += (-z * 2.0).exp() * p1 * s2;
        }
        p1 = -p1;
        *elem = s1 * coeff;
    }
    if n > 2 {
        let rz = calc_rz(z);
        // recur downward from the last two elements
        for k in (0..n - 2).rev() {
            y[k] = (((k + 1) as f64) + order) * (rz * y[k + 1]) + y[k + 2];
        }
    }
    if scaled_calculations {
        let exp_cz = cz.exp();
        for yi in y.iter_mut() {
            *yi *= exp_cz;
        }
    }
    Ok((y, nz))
}
