use approx::assert_relative_eq;
use rstest::rstest;

use crate::amos::bindings::zbesi_wrap;
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
#[case(4.0, 2.1)]
#[case(5.0, 2.0001)]
fn test_bessel_j(#[case] order: f64, #[case] z: f64) {
    let actual = bessel_j(order, z).unwrap();
    let expected = bessel_j_ref(order, z.into()).unwrap();
    assert_relative_eq!(actual, expected, epsilon = 1e-10)
}

#[rstest]
fn test_bessel_j_random() {
    for _ in 0..1000000 {
        let order = rand::random_range(-1000.0..1000.0);
        let z = rand::random_range(-1000.0..1000.0);
        let ans = bessel_j(order, z);
        if let Ok(actual) = ans {
            assert_relative_eq!(
                actual,
                bessel_j_ref(order, z.into()).unwrap(),
                epsilon = 1e-10
            )
        } else {
            assert_eq!(ans, Err(BesselError::NotYetImplemented));
        }
    }
}

#[rstest]
#[case(4.0, 2.1)]
#[trace]
// #[case(5.0, 2.0001)]
fn test_bessel_j_large_n_real(
    #[case] order: f64,
    #[case] z: f64,
    #[values(3/* , 4, 9, 100*/)] n: usize,
    #[values(Scaling::Unscaled, /*Scaling::Scaled*/)] scaling: Scaling,
) {
    let actual = zbesj(z.into(), order, scaling, n);

    let mut cyr = vec![0.0; n];
    let mut cyi = vec![0.0; n];
    let mut nz = 0;
    let mut ierr = 0;

    unsafe {
        zbesi_wrap(
            z,
            0.0,
            order,
            scaling as i32,
            n.try_into().unwrap(),
            cyr.as_mut_ptr().cast(),
            cyi.as_mut_ptr().cast(),
            &mut nz,
            &mut ierr,
        );
    }
    assert_eq!(nz, 0);
    assert_eq!(ierr, 0);
    let actual = actual.unwrap();
    let actual_real: Vec<_> = actual.0.iter().map(|z| z.re).collect();
    let actual_im: Vec<_> = actual.0.iter().map(|z| z.im).collect();
    assert_eq!(cyr, actual_real);
    assert_eq!(cyi, actual_im);
}
