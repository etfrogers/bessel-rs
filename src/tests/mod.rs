use approx::assert_relative_eq;

use crate::amos::gamma_ln;
use crate::bessel_j;
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

#[test]
fn test_bessel_j() {
    let (order, f) = (4.0, 2.1);
    assert_relative_eq!(
        bessel_j(order, f).unwrap(),
        bessel_j_ref(order, f.into()).unwrap(),
        epsilon = 1e-10
    )
}
