use num::{
    Complex,
    complex::{Complex64, ComplexFloat},
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
        if self.im() > 1000.0 * MACHINE_CONSTANTS.abs_error_tolerance {
            Err(BesselError::ComplexOutputForRealInput { output: *self })
        } else {
            Ok(self.re())
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
