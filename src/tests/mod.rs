use approx::assert_relative_eq;
use num::complex::Complex64;
use rstest::rstest;

use crate::amos::bindings::zbesj_wrap;
use crate::amos::{gamma_ln, zbesj};
use crate::{BesselError, GammaError, Scaling, bessel_j};
use complex_bessel_rs::bessel_j::bessel_j as bessel_j_ref;

#[test]
fn test_gamma_ln_hard_coded() {
    for i in 1..=100 {
        let f = i as f64;
        let actual = gamma_ln(f).unwrap();
        let expected = f.gamma().ln();
        assert_relative_eq!(actual, expected)
    }
}

#[test]
fn test_gamma_ln() {
    for f in [0.00000000001, 2.2, 0.5, 2.0_f64.sqrt(), 101.0] {
        // large values cause the "expected" calculation to overflow: the fortran version seems to work!
        let actual = gamma_ln(f).unwrap();
        let expected = f.gamma().ln();
        assert_relative_eq!(actual, expected, epsilon = 1e-10)
    }
}

#[rstest]
fn test_gamma_ln_negative() {
    for f in [0.00000000001, 2.2, 0.5, 2.0_f64.sqrt(), 101.0] {
        // large values cause the "expected" calculation to overflow: the fortran version seems to work!
        let actual = gamma_ln(-f);
        assert_eq!(actual, Err(GammaError::ZLessThanZero))
    }
}

#[rstest]
#[trace]
#[case(4.0, 2.1)] // z_power_series
#[case(5.0, 2.0001)] // z_power_series
#[case(340.0, 35.0001)] // z_power_series, iflag = true
fn test_bessel_j(#[case] order: f64, #[case] z: f64) {
    let actual = bessel_j(order, z).unwrap();
    let expected = bessel_j_ref(order, z.into()).unwrap();
    assert_relative_eq!(actual, expected, epsilon = 1e-10)
}

#[rstest]
fn test_bessel_j_random() {
    for _ in 0..1000000 {
        let order = rand::random_range(std::f64::EPSILON..1000.0);
        let zr = rand::random_range(-1000.0..1000.0);
        let zi = rand::random_range(-1000.0..1000.0);
        let z = Complex64::new(zr, zi);
        let ans = bessel_j(order, z);
        if let Ok(actual) = ans {
            assert_relative_eq!(actual, bessel_j_ref(order, z).unwrap(), epsilon = 1e-10)
        } else {
            assert_eq!(ans, Err(BesselError::NotYetImplemented));
        }
    }
}

#[rstest]
fn test_bessel_j_random_negative() {
    todo!("Add negative values")
}

#[rstest]
#[trace]
#[case(4.0, 2.1)] // z_power_series
#[case(5.0, 2.0001)] // z_power_series
#[case(340.0, 35.0001)] // z_power_series, iflag = true
fn test_bessel_j_large_n_real(
    #[case] order: f64,
    #[case] z: f64,
    #[values(3, 4, 9, 100)] n: usize,
    #[values(Scaling::Unscaled, Scaling::Scaled)] scaling: Scaling,
) {
    let actual = zbesj(z.into(), order, scaling, n);

    let mut cyr: Vec<f64> = Vec::with_capacity(n);
    let mut cyi: Vec<f64> = Vec::with_capacity(n);
    let mut nz = 0;
    let mut ierr = 0;

    let r_uninit = cyr.spare_capacity_mut();
    let i_uninit = cyi.spare_capacity_mut();
    unsafe {
        zbesj_wrap(
            z,
            0.0,
            order,
            scaling as i32,
            n.try_into().unwrap(),
            r_uninit.as_mut_ptr().cast(),
            i_uninit.as_mut_ptr().cast(),
            &mut nz,
            &mut ierr,
        );
        cyr.set_len(n);
        cyi.set_len(n);
    }
    let actual = actual.unwrap();
    assert_eq!(nz, actual.1.try_into().unwrap());
    assert_eq!(ierr, 0);
    for ((re, im), zi) in cyr.into_iter().zip(cyi).zip(actual.0) {
        assert_relative_eq!(re, zi.re);
        assert_relative_eq!(im, zi.im);
    }
}
