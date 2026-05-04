use approx::assert_relative_eq;
use rstest::rstest;
use rstest_reuse::{apply, template};

use crate::{GammaError, amos::gamma_ln};

// Below is the copy of code from std::sys::cmath to expose the tgamma function
// This avoids using the unstable float_gamma feature, but gives the same functionality
// to allow testing.
unsafe extern "C" {
    pub safe fn tgamma(n: f64) -> f64;
}

#[test]
fn test_gamma_ln_hard_coded() {
    for i in 1..=100 {
        let f = i as f64;
        let actual = gamma_ln(f).unwrap();
        let expected = tgamma(f).ln();
        assert_relative_eq!(actual, expected)
    }
}

#[template]
#[rstest]
fn f_values(
    #[values(0.00000000001, 1.5, 2.2, 0.5, 2.0_f64.sqrt(), 35.45764, 94.3567, 101.0)] f: f64,
) {
}

#[apply(f_values)]
fn test_gamma_ln(f: f64) {
    // large values cause the "expected" calculation to overflow: the fortran version seems to work!
    let actual = gamma_ln(f).unwrap();
    let expected = tgamma(f).ln();
    assert_relative_eq!(actual, expected, max_relative = 1e-10)
}

#[apply(f_values)]
fn test_gamma_ln_negative(f: f64) {
    // large values cause the "expected" calculation to overflow: the fortran version seems to work!
    let actual = gamma_ln(-f);
    assert_eq!(actual, Err(GammaError::ZLessThanZero))
}
