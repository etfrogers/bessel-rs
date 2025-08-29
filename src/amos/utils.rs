#![allow(clippy::excessive_precision)]
use num::{complex::Complex64, pow::Pow};

use super::{BesselError, BesselResult, MACHINE_CONSTANTS};

pub const RTPI: f64 = 0.159154943091895336;
pub const TWO_THIRDS: f64 = 6.66666666666666666e-01;

pub fn will_z_underflow(
    y: Complex64, //YR, YI,
    // Outputs: NZ
    ascle: f64,
    tol: f64,
) -> bool {
    // ***BEGIN PROLOGUE  ZUunderflowCHK
    // ***REFER TO z_power_series,ZUOIK,ZUNK1,ZUNK2,ZUNI1,ZUNI2,ZKSCL
    //
    //      Y ENTERS AS A SCALED QUANTITY WHOSE MAGNITUDE IS GREATER THAN
    //      EXP(-ALIM)=ASCLE=1.0E+3*d1mach(1)/TOL. THE TEST IS MADE TO SEE
    //      if THE MAGNITUDE OF THE REAL OR IMAGINARY PART WOULD UNDERFLOW
    //      WHEN Y IS SCALED (BY TOL) TO ITS PROPER VALUE. Y IS ACCEPTED
    //      if THE UNDERFLOW IS AT LEAST ONE PRECISION BELOW THE MAGNITUDE
    //      OF THE LARGEST COMPONENT; OTHERWISE THE PHASE ANGLE DOES NOT HAVE
    //      ABSOLUTE ACCURACY AND AN UNDERFLOW IS ASSUMED.
    //
    // ***ROUTINES CALLED  (NONE)
    // ***END PROLOGUE  ZUunderflowCHK
    //
    //     COMPLEX Y
    // DOUBLE PRECISION ASCLE, SS, ST, TOL, WR, WI, YR, YI
    // INTEGER NZ
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
