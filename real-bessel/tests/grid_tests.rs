use amos_bessel_rs::{bessel_j, bessel_y};
use approx::assert_relative_eq;
use real_bessel::{bessel_jn, bessel_yn};
use rstest::rstest;

#[macro_use]
mod common;

const Z_PARTS: [f64; 37] = [
    -50.0, -40.0, -30.0, -25.0, -20.0, -15.0, -12.0, -10.0, -8.0, -6.0, -4.0, -3.0, -2.0, -1.0,
    -0.5, -0.1, -0.001, -1e-6, 0.0, 1e-6, 0.001, 0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0,
    12.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0,
];

#[rstest]
fn test_jn_grid(#[values(0, 1, 2, 5, 10, 25, 50, 100, 200, 500, 1000, -2)] order: i32) {
    for &zr in &Z_PARTS {
        let expected = unwrap_real_bessel!(bessel_j, order as f64, zr);
        let actual = bessel_jn(order, zr);
        if expected.is_nan() {
            assert!(actual.is_nan() || actual.is_infinite());
        } else {
            assert_relative_eq!(actual, expected, epsilon = 1e-10, max_relative = 1e-10);
        }
    }
}

#[rstest]
fn test_yn_grid(#[values(0, 1, 2, 5, 10, 25, 50, 100, 200, 500, 1000, -2)] order: i32) {
    for &zr in &Z_PARTS {
        let mut zr_test = zr;
        let yn_local = |z| bessel_yn(order, z);
        let actual = get_real_y_bessel!(yn_local, zr_test, continue);
        let expected = unwrap_real_bessel!(bessel_y, order as f64, zr_test);
        if expected.is_nan() {
            assert!(actual.is_nan() || actual.is_infinite());
        } else {
            assert_relative_eq!(actual, expected, epsilon = 1e-10, max_relative = 1e-10);
        }
    }
}
