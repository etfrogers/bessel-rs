use num::{Complex, complex::ComplexFloat};

use crate::{
    BesselError::{LossOfSignificance, PartialLossOfSignificance},
    Scaling,
    amos::complex_airy,
    types::BesselFloat,
};

const MAX_AIRY_ITERATIONS: usize = 25;
pub fn airy_power_series<T: BesselFloat>(
    z: Complex<T>,
    return_derivative: bool,
    coeffs: (f64, f64),
) -> Complex<T> {
    let c1 = T::from_f64(coeffs.0);
    let c2 = T::from_f64(coeffs.1);

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
            let z_cubed = z.powi(3);
            let abs_z_3_over_2 = abs_z * abs_z_sq;
            let (ak, bk, ck, dk) = if return_derivative {
                (T::from_f64(3.0), T::ONE, T::from_f64(3.0), T::from_f64(5.0))
            } else {
                (
                    T::from_f64(2.0),
                    T::from_f64(3.0),
                    T::from_f64(4.0),
                    T::from_f64(3.0),
                )
            };

            let mut d1: T = ak * dk;
            let mut d2 = bk * ck;
            let mut min_d = d1.min(d2);
            let mut ak = if return_derivative {
                T::from_f64(33.0)
            } else {
                T::from_f64(24.0)
            };
            let mut bk = if return_derivative {
                T::from_f64(21.0)
            } else {
                T::from_f64(30.0)
            };
            for _ in 0..MAX_AIRY_ITERATIONS {
                term1 = term1 * z_cubed / d1;
                s1 += term1;
                term2 = term2 * z_cubed / d2;
                s2 += term2;
                a_term = a_term * abs_z_3_over_2 / min_d;
                d1 += ak;
                d2 += bk;
                min_d = d1.min(d2);
                if a_term < T::MACHINE_CONSTANTS.abs_error_tolerance * min_d {
                    break;
                }
                ak += T::from_f64(18.0);
                bk += T::from_f64(18.0);
            }
        }
        (s1, s2)
    };

    if return_derivative {
        z_floor.powi(2) * s1 * (c1 / T::two()) - s2 * c2
    } else {
        s1 * c1 - c2 * z_floor * s2
    }
}

pub fn airy_pair<T: BesselFloat>(z: Complex<T>) -> (Complex<T>, Complex<T>) {
    //note that ZAIRY calls in fortran code ignore IERR (using IDUM)
    let evaluate_airy_and_unwrap =
        |is_derivative| match complex_airy(z, is_derivative, Scaling::Scaled) {
            Ok((y, _)) => y,
            Err(PartialLossOfSignificance { y, nz: _ }) => y[0],
            // If loss of significance, Fortran code would continue with un-initialised y,
            // which is usually ~=0. As long as it is << d_airy, the logic below means
            // it will not matter what the precise value is
            Err(LossOfSignificance) => T::C_ZERO,
            Err(err) => panic!(
                "An error {:?} was generated, which is not handled by the Amos code",
                err
            ),
        };
    (
        evaluate_airy_and_unwrap(false),
        evaluate_airy_and_unwrap(true),
    )
}
