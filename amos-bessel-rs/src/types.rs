use approx::relative_eq;
use num::{
    Complex,
    complex::{Complex64, ComplexFloat},
    pow::Pow,
};
use thiserror::Error;

use crate::amos::MACHINE_CONSTANTS;

pub(crate) type BesselValues<T = usize> = (Vec<Complex64>, T);
pub(crate) type BesselResult<T = BesselValues> = Result<T, BesselError>;

pub trait BackTo<T> {
    fn back_to(&self) -> BesselResult<T>;
}

impl BackTo<Complex64> for Complex64 {
    fn back_to(&self) -> BesselResult<Complex64> {
        Ok(*self)
    }
}

impl BackTo<f64> for f64 {
    fn back_to(&self) -> BesselResult<f64> {
        Ok(*self)
    }
}

impl BackTo<f64> for Complex64 {
    fn back_to(&self) -> BesselResult<f64> {
        const MARGIN: f64 = 1000.0;
        let tol = MARGIN * MACHINE_CONSTANTS.abs_error_tolerance;
        if self.im().abs() < tol || self.im().abs() < self.re().abs() * tol {
            Ok(self.re())
        } else {
            Err(BesselError::ComplexOutputForRealInput { output: *self })
        }
    }
}

/// Original Docs:
/// IERR   - ERROR FLAG
///         IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
///         IERR=1, INPUT ERROR   - NO COMPUTATION
///        IERR=2, OVERFLOW      - NO COMPUTATION, order TOO
///                 LARGE OR CABS(Z) TOO SMALL OR BOTH
///         IERR=3, CABS(Z) OR order+N-1 LARGE - COMPUTATION DONE
///                BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
///                 REDUCTION PRODUCE LESS THAN HALF OF MACHINE
///                 ACCURACY
///         IERR=4, CABS(Z) OR order+N-1 TOO LARGE - NO COMPUTA-
///                 TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
///                 CANCE BY ARGUMENT REDUCTION
///         IERR=5, ERROR              - NO COMPUTATION,
///                 ALGORITHM TERMINATION CONDITION NOT MET
#[derive(Error, Debug, PartialEq)]
#[repr(i32)]
pub enum BesselError {
    // This correpsonds to IERR
    // 0 = no error
    #[error("Invalid input: {details}")]
    InvalidInput { details: String } = 1,
    #[error("Overflow: order TOO LARGE OR CABS(Z) TOO SMALL OR BOTH")]
    Overflow = 2, //{ too_large: bool },
    #[error("Partial loss of significance in output. Losssy values returned.")]
    PartialLossOfSignificance { y: Vec<Complex64>, nz: usize } = 3,
    #[error("Loss of too much significance in output")]
    LossOfSignificance = 4,
    #[error("Algorithm failed to terminate")]
    DidNotConverge = 5,
    #[error("Real input returned complex output. Output value {output}")]
    ComplexOutputForRealInput { output: Complex<f64> } = 6,
}

impl BesselError {
    pub fn error_code(&self) -> i32 {
        match self {
            BesselError::InvalidInput { .. } => 1,
            BesselError::Overflow => 2,
            BesselError::PartialLossOfSignificance { .. } => 3,
            BesselError::LossOfSignificance => 4,
            BesselError::DidNotConverge => 5,
            BesselError::ComplexOutputForRealInput { .. } => 6,
        }
    }
}

#[macro_export]
macro_rules! simple_bessel_wrapper {
    (
        $(#[$meta:meta])*
        $base_func:ident // We only pass the base function name now!
    ) => {
        // The paste! macro allows us to create new identifiers
        paste! {
            $(#[$meta])*
            // [<simple_ $base_func>] concatenates into simple_bessel_j
            #[inline]
            fn [<$base_func _single>](order: f64, z: Complex64) -> Result<Complex64, BesselError> {
                let (result_vec, _nz) = [<complex_$base_func>](z, order, Scaling::Unscaled, 1)?;
                Ok(result_vec[0])
            }
        }
    };
}

// Dead code is used in the tests, but not in the main library, so we allow it here.
#[allow(dead_code)]
pub(crate) struct Tolerances {
    pub max_relative: f64,
    pub magnitude_diff: f64,
    pub correction: f64,
    pub tol_re: f64,
    pub tol_im: f64,
    pub used_ref_re: bool,
    pub used_ref_im: bool,
    pub margin: f64,
}

impl Tolerances {
    pub(crate) fn new(
        actual: &Complex64,
        expected: &Complex64,
        reference: Option<&Complex64>,
        margin: f64,
    ) -> Self {
        let max_relative = MACHINE_CONSTANTS.abs_error_tolerance;

        let max_im = actual.im.abs().max(expected.im.abs());
        let max_re = actual.re.abs().max(expected.re.abs());
        let magnitude_diff = (max_re.log10() - max_im.log10()).abs();
        let limit = 1.0 / max_relative;
        let correction = 10.0.pow(magnitude_diff).min(limit);

        let (mut tol_re, mut tol_im) = if max_re > max_im {
            (max_relative, max_relative * correction)
        } else {
            (max_relative * correction, max_relative)
        };

        let mut used_ref_re = false;
        let mut used_ref_im = false;
        if let Some(ref_val) = reference {
            let (_, ref_diffs) = abs_rel_errors_cmplx(expected, ref_val);
            let re_rel_diff = ref_diffs.re;
            let im_rel_diff = ref_diffs.im;

            used_ref_re = re_rel_diff > tol_re;
            if used_ref_re {
                tol_re = re_rel_diff;
            }
            used_ref_im = im_rel_diff > tol_im;
            if used_ref_im {
                tol_im = im_rel_diff;
            }
        }
        Self {
            max_relative,
            magnitude_diff,
            correction,
            tol_re,
            tol_im,
            used_ref_re,
            used_ref_im,
            margin,
        }
    }
}

pub(crate) fn complex_relative_eq(a: &Complex64, b: &Complex64, tolerances: &Tolerances) -> bool {
    if relative_eq!(
        a,
        b,
        max_relative = tolerances.margin * tolerances.max_relative
    ) {
        return true;
    }
    let (_, rel_e_re) = abs_rel_errors(a.re, b.re);
    let (_, rel_e_im) = abs_rel_errors(a.im, b.im);

    rel_e_re < tolerances.margin * tolerances.tol_re
        && rel_e_im < tolerances.margin * tolerances.tol_im
}

fn abs_rel_errors(a: f64, b: f64) -> (f64, f64) {
    let abs_e = (a - b).abs();
    let rel_e = abs_e / a.abs().max(b.abs());
    (abs_e, rel_e)
}

pub(crate) fn abs_rel_errors_cmplx(a: &Complex64, b: &Complex64) -> (Complex64, Complex64) {
    let (abs_e_r, rel_e_r) = abs_rel_errors(a.re, b.re);
    let (abs_e_i, rel_e_i) = abs_rel_errors(a.im, b.im);
    (
        Complex64::new(abs_e_r, abs_e_i),
        Complex64::new(rel_e_r, rel_e_i),
    )
}
