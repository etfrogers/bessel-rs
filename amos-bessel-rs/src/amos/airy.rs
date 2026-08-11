use num::{Complex, complex::ComplexFloat};

use crate::{
    BesselError, Scaling,
    amos::{MachineConsts, complex_airy},
    types::BesselFloat,
};

const MAX_AIRY_ITERATIONS: usize = 25;

/// Evaluates the Airy function (or its derivative) via its Maclaurin power series.
///
/// The Airy function Ai(z) is given by:
///     Ai(z)  = c1 * f(z) - c2 * g(z)
///     Ai'(z) = c1 * f'(z) - c2 * g'(z)
///
/// Where:
///     f(z) = 1 + 1/3! z^3 + (1*4)/6! z^6 + ...
///     g(z) = z * (1 + 2/4! z^3 + (2*5)/7! z^6 + ...)
///
/// When `return_derivative` is true, it evaluates the series for f'(z) and g'(z).
/// The loop computes the terms for the two series simultaneously (`sum1` and `sum2`).
pub(crate) fn airy_power_series<T: BesselFloat>(
    z: Complex<T>,
    return_derivative: bool,
    coeffs: (f64, f64),
) -> Complex<T> {
    let mc: &MachineConsts<T> = T::MACHINE_CONSTANTS;
    let c1 = T::from_f64(coeffs.0);
    let c2 = T::from_f64(coeffs.1);

    let abs_z = z.abs();

    let (sum1, sum2) = if abs_z < mc.abs_error_tolerance {
        (T::C_ONE, T::C_ONE)
    } else {
        let abs_z_sq = abs_z.powi(2);
        let mut sum1 = T::C_ONE;
        let mut sum2 = T::C_ONE;

        if abs_z_sq >= mc.abs_error_tolerance / abs_z {
            let mut term1 = T::C_ONE;
            let mut term2 = T::C_ONE;
            let mut convergence_term = T::ONE;
            let z_cubed = z.powi(3);
            let abs_z_3_over_2 = abs_z * abs_z_sq;

            // The initial denominator multipliers for the power series terms (k=1).
            // For f(z), the first multiplier is 3 * 2 = 6.0
            // For g(z), the first multiplier is 4 * 3 = 12.0
            // For f'(z) and g'(z), the multipliers shift to 5 * 3 = 15.0 and 3 * 1 = 3.0
            let mut denom1 = T::from_f64(if return_derivative { 15.0 } else { 6.0 });
            let mut denom2 = T::from_f64(if return_derivative { 3.0 } else { 12.0 });
            let mut min_denom = denom1.min(denom2);

            // `step1` and `step2` are the amounts added to the denominators on each iteration.
            let mut step1 = T::from_f64(if return_derivative { 33.0 } else { 24.0 });
            let mut step2 = T::from_f64(if return_derivative { 21.0 } else { 30.0 });

            for _ in 0..MAX_AIRY_ITERATIONS {
                term1 = term1 * z_cubed / denom1;
                sum1 += term1;
                term2 = term2 * z_cubed / denom2;
                sum2 += term2;
                convergence_term = convergence_term * abs_z_3_over_2 / min_denom;
                denom1 += step1;
                denom2 += step2;
                min_denom = denom1.min(denom2);
                if convergence_term < mc.abs_error_tolerance * min_denom {
                    break;
                }
                step1 += T::from_f64(18.0);
                step2 += T::from_f64(18.0);
            }
        }
        (sum1, sum2)
    };

    let z_floor = if abs_z < mc.underflow_limit {
        T::C_ZERO
    } else {
        z
    };

    if return_derivative {
        z_floor.powi(2) * sum1 * (c1 / T::two()) - sum2 * c2
    } else {
        sum1 * c1 - c2 * z_floor * sum2
    }
}

pub(crate) fn airy_pair<T: BesselFloat>(z: Complex<T>) -> (Complex<T>, Complex<T>) {
    //note that ZAIRY calls in fortran code ignore IERR (using IDUM)
    let evaluate_airy_and_unwrap =
        |is_derivative| match complex_airy(z, is_derivative, Scaling::Scaled) {
            Ok((y, _)) => y,
            Err(BesselError::PartialLossOfSignificance { y, n_zeros: _ }) => y[0],
            // If loss of significance, Fortran code would continue with un-initialised y,
            // which is usually ~=0. As long as it is << d_airy, the logic below means
            // it will not matter what the precise value is
            Err(BesselError::LossOfSignificance) => T::C_ZERO,
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
