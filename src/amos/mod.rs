pub use gamma_ln::{gamma_ln, GammaError};
use thiserror::Error;
pub use translator::zbesj;
mod gamma_ln;
mod machine;
mod translator;
// mod zbesh;

const ZEROR: f64 = 0.0;
const ZEROI: f64 = 0.0;
const CONER: f64 = 1.0;
const CONEI: f64 = 0.0;

// const COMPLEX_ONE: Complex64 = Complex64::one();

#[derive(Error, Debug)]
#[repr(u16)]
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
}

pub enum HankelKind {
    First,
    Second,
}

#[derive(Debug, PartialEq, Eq)]
pub enum Scaling {
    Unscaled = 1,
    Scaled = 2,
}
