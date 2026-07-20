/// Internal trait abstracting over Bessel function implementations.
///
/// `j` and `y` take the order as `f64` (as used internally by the algorithm),
/// and each backend is responsible for converting to its own order type.
pub(crate) trait BesselBackend {
    fn j(order: f64, x: f64) -> f64;
    fn y(order: f64, x: f64) -> f64;
}

pub(crate) mod amos;
pub(crate) mod real;
