use num::{
    Integer,
    complex::{Complex, ComplexFloat},
};

use super::{Scaling, utils::RTPI};
use crate::types::{BesselError::*, BesselResult};
use crate::{amos::utils::calc_rz, types::BesselFloat};

/// asymptotic_i computes the I bessel function for real(z) >= 0.0 by
/// means of the asymptotic expansion for large z.abs() in the
/// region z.abs() > rl.max(order.pow(2.0)/2).
///
/// Originally ZASYI
pub fn asymptotic_i<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    n: usize,
) -> BesselResult<T, usize> {
    let nz = 0;
    let mut y = T::c_zeros(n);
    let abs_z = z.abs();
    //-----------------------------------------------------------------------;
    //     OVERFLOW TEST;
    //-----------------------------------------------------------------------;
    let recip_abs_z = T::one() / abs_z;

    let cz = match scaling {
        Scaling::Unscaled => z,
        Scaling::Scaled => Complex::<T>::new(T::zero(), z.im),
    };

    if cz.re.abs() > T::MACHINE_CONSTANTS.exponent_limit {
        return Err(Overflow);
    }
    let scaled_calculations = cz.re.abs() > T::MACHINE_CONSTANTS.approximation_limit;
    let mut coeff = (T::from_f64(RTPI) * z.conj() * recip_abs_z.powi(2)).sqrt();
    if !scaled_calculations {
        coeff *= cz.exp();
    }
    let ez = z * T::from_f64(8.0);
    //-----------------------------------------------------------------------;
    // When z is imaginary, the error test must be made relative to the
    // first reciprocal power since this is the leading term of the
    // expansion for the imaginary part.
    //-----------------------------------------------------------------------;
    let abs_ez = T::from_f64(8.0) * abs_z;
    let s = T::MACHINE_CONSTANTS.abs_error_tolerance / abs_ez;
    let max_iterations = (T::MACHINE_CONSTANTS.asymptotic_z_limit * T::two())
        .to_i64()
        .unwrap()
        + 2;
    let mut p1 = if z.im == T::zero() {
        T::C_ZERO
    } else {
        //-----------------------------------------------------------------------;
        // Calculate (pi*(0.5+fnu+n-il)*i).exp() to minimize losses of;
        // significance when fnu or n is large;
        //-----------------------------------------------------------------------;
        let arg = order.fract() * T::PI();
        let ak = -arg.sin();
        let mut bk = arg.cos();
        if z.im < T::zero() {
            bk = -bk;
        };
        let p1 = Complex::<T>::new(ak, bk);
        if (order.to_usize().unwrap() + n).is_even() {
            -p1
        } else {
            p1
        }
    };
    for (k, elem) in y.iter_mut().enumerate().rev().take(2.min(n)) {
        let (mut s1, s2) = {
            // this block is just to contain the large number of mutable variables in a small space
            let modified_order = order + T::from_usize(k);
            let modified_order_sqr = (T::two() * modified_order).powf(T::two());
            let mut sqk = modified_order_sqr - T::one();
            let atol = s * sqk.abs();
            let mut sign = T::one();
            let mut cs1 = T::C_ONE;
            let mut cs2 = T::C_ONE;
            let mut ck = T::C_ONE;
            let mut ak = T::zero();
            let mut aa = T::one();
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
                ak += T::from_f64(8.0);
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
        if z.re * T::two() < T::MACHINE_CONSTANTS.exponent_limit {
            s1 += (-z * T::two()).exp() * p1 * s2;
        }
        p1 = -p1;
        *elem = s1 * coeff;
    }
    if n > 2 {
        let rz = calc_rz(z);
        // recur downward from the last two elements
        for k in (0..n - 2).rev() {
            y[k] = (rz * y[k + 1]) * (T::from_usize(k + 1) + order) + y[k + 2];
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
