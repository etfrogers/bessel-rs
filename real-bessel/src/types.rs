use thiserror::Error;

/// BesselError represents errors that can occur during Bessel function calculations.
#[derive(Error, Debug, PartialEq)]
#[repr(i32)]
pub enum BesselError {
    /// Error for negative input values, in the `y0`, `y1` and `yn` functions.
    #[error(
        "{function} is complex for x < 0, and this function returns only real values. x = {input}"
    )]
    NegativeInputForYFunction {
        /// The name of the function that caused the error.
        function: String,
        /// The input value that caused the error.
        input: f64,
    } = 1,
}
