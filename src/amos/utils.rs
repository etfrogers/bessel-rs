#![allow(clippy::excessive_precision)]
use num::{complex::Complex64, pow::Pow};

use super::{BesselError, BesselResult, MACHINE_CONSTANTS};

pub const RTPI: f64 = 0.159154943091895336;
pub const TWO_THIRDS: f64 = 6.66666666666666666e-01;
pub const RT_THREE: f64 = 1.73205080757;
pub const AIC: f64 = 1.265512123484645396; // == gamma_ln(-0.5).re

pub(crate) fn imaginary_dominant(z: Complex64) -> bool {
    z.im.abs() > z.re.abs() * RT_THREE
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
pub(crate) fn will_underflow(y: Complex64, ascle: f64, tol: f64) -> bool {
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

pub fn is_sigificance_lost(
    z_abs: f64,
    modified_order: f64,
    modify_threshold: bool,
) -> BesselResult<bool> {
    let f64_precision_limit = 0.5 / MACHINE_CONSTANTS.abs_error_tolerance;
    // TODO the below is limited to i32: could push to 64 later, but would change compare to fortran
    let integer_size_limit = (i32::MAX as f64) * 0.5;
    let mut upper_size_limit = f64_precision_limit.min(integer_size_limit);
    if modify_threshold {
        upper_size_limit = upper_size_limit.pow(TWO_THIRDS);
    }
    if z_abs > upper_size_limit || modified_order > upper_size_limit {
        return Err(BesselError::LossOfSignificance);
    }
    let scaling_limit = upper_size_limit.sqrt();
    Ok((z_abs > scaling_limit) || (modified_order > scaling_limit))
}

pub(crate) fn sanitise_inputs(
    z: Complex64,
    order: f64,
    n: usize,
    check_z_zero: bool,
) -> Result<(), BesselError> {
    let mut err = None;
    if check_z_zero && z.re == 0.0 && z.im == 0.0 {
        err = Some("z must not be zero");
    }
    if order < 0.0_f64 {
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
