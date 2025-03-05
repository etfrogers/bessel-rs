pub use gamma_ln::{GammaError, gamma_ln};
use thiserror::Error;
pub use translator::zbesj;
pub use z_power_series::z_power_series;
pub mod bindings;
mod gamma_ln;
mod machine;
mod translator;
mod utils;
mod z_power_series;
// mod zbesh;

#[derive(Error, Debug, PartialEq, Eq)]
#[repr(i32)]
pub enum BesselError {
    // This correpsonds IERRR
    // 0 = no error
    #[error("Invalid input: {details}")]
    InvalidInput { details: String } = 1,
    #[error("Overflow: order TOO LARGE OR CABS(Z) TOO SMALL OR BOTH")]
    Overflow = 2, //{ too_large: bool },
    // IERR = 3 is a warning (and hence some return value, I think) and needs handling elsewhere
    #[error("Loss of too much significance in output")]
    LossOfSignificance = 4,
    #[error("Algorithm failed to terminate")]
    DidNotConverge = 5,
    // Original Docs:
    // IERR   - ERROR FLAG
    //         IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
    //         IERR=1, INPUT ERROR   - NO COMPUTATION
    //         IERR=2, OVERFLOW      - NO COMPUTATION, order TOO
    //                 LARGE OR CABS(Z) TOO SMALL OR BOTH
    //         IERR=3, CABS(Z) OR order+N-1 LARGE - COMPUTATION DONE
    //                 BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
    //                 REDUCTION PRODUCE LESS THAN HALF OF MACHINE
    //                 ACCURACY
    //         IERR=4, CABS(Z) OR order+N-1 TOO LARGE - NO COMPUTA-
    //                 TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
    //                 CANCE BY ARGUMENT REDUCTION
    //         IERR=5, ERROR              - NO COMPUTATION,
    //                 ALGORITHM TERMINATION CONDITION NOT MET
    #[error("not yet implemented")]
    NotYetImplemented = 100,
}

impl BesselError {
    pub fn error_code(&self) -> i32 {
        match self {
            BesselError::InvalidInput { .. } => 1,
            BesselError::Overflow => 2,
            BesselError::LossOfSignificance => 4,
            BesselError::DidNotConverge => 5,
            BesselError::NotYetImplemented => 100,
        }
    }
}

pub enum HankelKind {
    First,
    Second,
}

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
#[repr(i32)]
pub enum Scaling {
    Unscaled = 1,
    Scaled = 2,
}
