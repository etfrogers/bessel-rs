use super::BesselBackend;
use real_bessel::{jn, yn};

/// Backend using the `real-bessel` crate.
/// Faster than AMOS but restricted to integer orders.
/// The order is received as `f64` from the algorithm and cast to `i32`.
/// The calling function (bessr) is set up so that the input should already be integer values.
pub(crate) struct RealBackend;

impl BesselBackend for RealBackend {
    fn j(order: f64, x: f64) -> f64 {
        // Note: real_bessel::jn returns f64 directly (no Result).
        jn(order as i32, x)
    }

    fn y(order: f64, x: f64) -> f64 {
        // Note: real_bessel::yn returns Result<f64, BesselError> (unlike jn),
        // because yn(x) is undefined for x <= 0 and returns an error there.
        yn(order as i32, x).expect("real_bessel yn computation failed")
    }
}
