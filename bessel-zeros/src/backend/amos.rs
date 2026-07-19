use amos_bessel_rs::{bessel_j, bessel_y};
use super::BesselBackend;

/// Backend using the AMOS algorithm via `amos-bessel-rs`.
/// Supports any real (possibly non-integer) order.
pub(crate) struct AmosBackend;

impl BesselBackend for AmosBackend {
    fn j(order: f64, x: f64) -> f64 {
        bessel_j(order, x).expect("amos bessel_j computation failed")
    }

    fn y(order: f64, x: f64) -> f64 {
        bessel_y(order, x).expect("amos bessel_y computation failed")
    }
}
