use thiserror::Error;

/// Errors that can occur during Bessel function calculations.
///
/// At present the only variant is [`BesselError::NegativeInputForYFunction`],
/// returned by [`y0`](crate::y0), [`y1`](crate::y1) and [`yn`](crate::yn)
/// when called with `x < 0`, because the Y functions are complex for
/// negative arguments and this crate returns only real values.
#[derive(Error, Debug, Clone, PartialEq)]
pub enum BesselError {
    /// Returned by `y0`, `y1` and `yn` when `x < 0`.
    #[error(
        "{function} is complex for x < 0, and this function returns only real values. x = {input}"
    )]
    NegativeInputForYFunction {
        /// The name of the function that caused the error.
        function: String,
        /// The input value that caused the error.
        input: f64,
    },
}
