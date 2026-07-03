#![allow(clippy::excessive_precision)]
use num::complex::{Complex, ComplexFloat};

use crate::types::{BesselError, BesselFloat};

pub const RTPI: f64 = 0.159154943091895336;
pub const RT_THREE: f64 = 1.73205080757;
pub const AIC: f64 = 1.265512123484645396; // == gamma_ln(-0.5).re

/// This slightly odd form of calculation avoids overflow/underflow
/// when z is large/small respectively.
pub(crate) fn calc_rz<T: BesselFloat>(z: Complex<T>) -> Complex<T> {
    //2.0 * z.conj() / abs_z.powi(2)
    let r_abs_z = T::one() / z.abs();
    let intermediate = z.conj() * r_abs_z;
    (intermediate + intermediate) * r_abs_z
}

pub(crate) fn imaginary_dominant<T: BesselFloat>(z: Complex<T>) -> bool {
    z.im.abs() > z.re.abs() * T::from_f64(RT_THREE)
}

/// `y` enters as a scaled quantity whose magnitude is greater than
/// (all names below are from [MachineConsts])
/// `(-approximation_limit).exp() = absolute_approximation_limit
///  = 2.0 * f64::MIN_POSITIVE/abs_error_tolerance`.
/// The test is made to see
/// if the magnitude of the real or imaginary part would underflow
/// when y is scaled (by tol) to its proper value. `y` is accepted
/// if the underflow is at least one precision below the magnitude
/// of the largest component; otherwise the phase angle does not have
/// absolute accuracy and an underflow is assumed
pub(crate) fn will_underflow<T: BesselFloat>(y: Complex<T>, ascle: T, tol: T) -> bool {
    let re_abs = y.re.abs();
    let im_abs = y.im.abs();
    let min_abs_component = re_abs.min(im_abs);
    if min_abs_component > ascle {
        false
    } else {
        let max_abs_component = re_abs.max(im_abs);
        max_abs_component < min_abs_component / tol
    }
}

pub fn is_significance_lost<T: BesselFloat>(
    z_abs: T,
    modified_order: T,
    modify_threshold: bool,
) -> Result<bool, BesselError<T>> {
    let f64_precision_limit = T::half() / T::MACHINE_CONSTANTS.abs_error_tolerance;
    // TODO the below is limited to i32: could push to 64 later, but would change compare to fortran
    let integer_size_limit = T::from_f64((i32::MAX as f64) * 0.5);
    let mut upper_size_limit = f64_precision_limit.min(integer_size_limit);
    if modify_threshold {
        upper_size_limit = upper_size_limit.powf(T::TWO_THIRDS);
    }
    if z_abs > upper_size_limit || modified_order > upper_size_limit {
        return Err(BesselError::LossOfSignificance);
    }
    let scaling_limit = upper_size_limit.sqrt();
    Ok((z_abs > scaling_limit) || (modified_order > scaling_limit))
}

pub(crate) fn sanitise_inputs<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    n: usize,
    check_z_zero: bool,
) -> Result<(), BesselError<T>> {
    let mut err = None;
    if check_z_zero && z.re == T::zero() && z.im == T::zero() {
        err = Some("z must not be zero");
    }
    if order < T::zero() {
        err = Some("order must be positive");
    };
    if n < 1 {
        err = Some("N must be >= 1");
    };
    if let Some(details) = err {
        Err(BesselError::InvalidInput {
            details: details.to_owned(),
        })
    } else {
        Ok(())
    }
}
