use thiserror::Error;

pub(crate) type BesselResult<T> = Result<T, BesselError>;

#[derive(Error, Debug, PartialEq)]
#[repr(i32)]
pub enum BesselError {
    #[error("Invalid input: {details}")]
    InvalidInput { details: String } = 1,
    #[error("Overflow: order TOO LARGE OR CABS(Z) TOO SMALL OR BOTH")]
    Overflow = 2,
    #[error("Partial loss of significance in output. Losssy values returned.")]
    PartialLossOfSignificance = 3,
    #[error("Loss of too much significance in output")]
    LossOfSignificance = 4,
    #[error("Algorithm failed to terminate")]
    DidNotConverge = 5,
}
