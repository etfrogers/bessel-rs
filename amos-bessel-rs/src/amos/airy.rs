#![allow(non_snake_case)]
use num::{Complex, complex::ComplexFloat};

use crate::{
    BesselError::{LossOfSignificance, PartialLossOfSignificance},
    Scaling,
    amos::complex_airy,
    types::BesselFloat,
};

pub fn airy_power_series<T: BesselFloat>(
    z: Complex<T>,
    return_derivative: bool,
    coeffs: (f64, f64),
) -> Complex<T> {
    let c1 = T::from_f64(coeffs.0);
    let c2 = T::from_f64(coeffs.1);
    let float_is_derivative = if return_derivative { T::ONE } else { T::ZERO };

    let abs_z = z.abs();
    let z_floor = if abs_z < T::MACHINE_CONSTANTS.underflow_limit {
        T::C_ZERO
    } else {
        z
    };
    let (s1, s2) = if abs_z < T::MACHINE_CONSTANTS.abs_error_tolerance {
        (T::C_ONE, T::C_ONE)
    } else {
        let abs_z_sq = abs_z * abs_z;
        let mut s1 = T::C_ONE;
        let mut s2 = T::C_ONE;

        if abs_z_sq >= T::MACHINE_CONSTANTS.abs_error_tolerance / abs_z {
            let mut term1 = T::C_ONE;
            let mut term2 = T::C_ONE;
            let mut a_term = T::ONE;
            let z3 = z.powf(T::from_f64(3.0));
            let AZ3 = abs_z * abs_z_sq;
            let (AK, BK, CK, DK) = (
                T::from_f64(2.0) + float_is_derivative,
                T::from_f64(3.0) - T::two() * float_is_derivative,
                T::from_f64(4.0) - float_is_derivative,
                T::from_f64(3.0) + T::two() * float_is_derivative,
            );
            let mut D1: T = AK * DK;
            let mut D2 = BK * CK;
            let mut AD = D1.min(D2);
            let mut AK = T::from_f64(24.0) + T::from_f64(9.0) * float_is_derivative;
            let mut BK = T::from_f64(30.0) - T::from_f64(9.0) * float_is_derivative;
            for _ in 0..25 {
                term1 = term1 * z3 / D1;
                s1 += term1;
                term2 = term2 * z3 / D2;
                s2 += term2;
                a_term = a_term * AZ3 / AD;
                D1 += AK;
                D2 += BK;
                AD = D1.min(D2);
                if a_term < T::MACHINE_CONSTANTS.abs_error_tolerance * AD {
                    break;
                }
                AK += T::from_f64(18.0);
                BK += T::from_f64(18.0);
            }
        }
        (s1, s2)
    };

    if return_derivative {
        z_floor.powf(T::two()) * s1 * (c1 / T::two()) - s2 * c2
    } else {
        s1 * c1 - c2 * z_floor * s2
    }
}

pub fn airy_pair<T: BesselFloat>(z: Complex<T>) -> (Complex<T>, Complex<T>) {
    //note that ZAIRY calls in fortran code ignore IERR (using IDUM)
    let airy = match complex_airy(z, false, Scaling::Scaled) {
        Ok((y, _)) => y,
        Err(PartialLossOfSignificance { y, nz: _ }) => y[0],
        // If loss of significance, Fortran code would continue with un-initialised y,
        // which is usually ~=0. Also long as it is << d_airy, the logic below means
        // it will not matter what the precise value is
        Err(LossOfSignificance) => T::C_ZERO,
        Err(err) => {
            panic!(
                "An error {:?} was generated, which is not handled by the Amos code",
                err
            )
        }
    };
    let d_airy = match complex_airy(z, true, Scaling::Scaled) {
        Ok((y, _)) => y,
        Err(PartialLossOfSignificance { y, nz: _ }) => y[0],
        // If loss of significance, Fortran code would continue with un-initialised y,
        // which is usually ~=0. Also long as it is << a_airy, the logic below means
        // it will not matter what the precise value is
        Err(LossOfSignificance) => T::C_ZERO,
        Err(err) => {
            panic!(
                "An error {:?} was generated, which is not handled by the Amos code",
                err
            )
        }
    };
    (airy, d_airy)
}
