use super::BesselBackend;
use amos_bessel_rs::{
    Scaling,
    amos::{complex_bessel_j, complex_bessel_y},
};

/// Backend using the AMOS algorithm via `amos-bessel-rs`.
/// Supports any real (possibly non-integer) order.
pub(crate) struct AmosBackend;

impl BesselBackend for AmosBackend {
    fn j_pair(order: f64, x: f64) -> (f64, f64) {
        let (y, _n_zeros) = complex_bessel_j(x.into(), order, Scaling::Unscaled, 2)
            .expect("amos bessel_j computation failed");
        (y[0].re, y[1].re)
    }

    fn y_pair(order: f64, x: f64) -> (f64, f64) {
        let (y, _n_zeros) = complex_bessel_y(x.into(), order, Scaling::Unscaled, 2)
            .expect("amos bessel_y computation failed");
        (y[0].re, y[1].re)
    }
}
