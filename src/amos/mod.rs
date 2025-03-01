use thiserror::Error;

mod machine;
mod translator;
mod zbesh;

#[derive(Error, Debug)]
pub enum BesselError {
    #[error("Invalid input: {details}")]
    InvalidInput { details: String },
    #[error("Overflow: order TOO LARGE OR CABS(Z) TOO SMALL OR BOTH")]
    Overflow, //{ too_large: bool },
    #[error("Loss of too much significance in output")]
    LossOfSignificance,
    #[error("Algorithm failed to terminate")]
    DidNotConverge,
}

pub enum HankelKind {
    First,
    Second,
}

pub enum Scaling {
    Unscaled,
    Scaled,
}
