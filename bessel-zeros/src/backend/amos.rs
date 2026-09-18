use super::BesselBackend;
use amos_bessel_rs::{
    Scaling,
    amos::{complex_bessel_j_into, complex_bessel_y_into},
};
use num::Complex;

/// Backend using the AMOS algorithm via `amos-bessel-rs`.
/// Supports any real (possibly non-integer) order.
pub(crate) struct AmosBackend;

impl BesselBackend for AmosBackend {
    fn j_pair(order: f64, x: f64) -> (f64, f64) {
        let mut buf = [Complex::new(0.0, 0.0); 2];
        complex_bessel_j_into(x.into(), order, Scaling::Unscaled, &mut buf)
            .expect("amos bessel_j computation failed");
        (buf[0].re, buf[1].re)
    }

    fn y_pair(order: f64, x: f64) -> (f64, f64) {
        let mut buf = [Complex::new(0.0, 0.0); 2];
        complex_bessel_y_into(x.into(), order, Scaling::Unscaled, &mut buf)
            .expect("amos bessel_y computation failed");
        (buf[0].re, buf[1].re)
    }
}
