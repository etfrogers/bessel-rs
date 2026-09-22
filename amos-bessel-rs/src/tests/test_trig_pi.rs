use approx::assert_relative_eq;
use num::Complex;

use crate::amos::utils::{cis_pi, cos_pi, from_polar_pi, sin_cos_pi, sin_pi};

#[test]
fn test_sin_cos_pi_exact_axes() {
    let exact_cases: &[(f64, f64, f64)] = &[
        // (x, expected_sin, expected_cos)
        (0.0, 0.0, 1.0),
        (0.5, 1.0, 0.0),
        (1.0, 0.0, -1.0),
        (1.5, -1.0, 0.0),
        (2.0, 0.0, 1.0),
        (2.5, 1.0, 0.0),
        (3.0, 0.0, -1.0),
        (3.5, -1.0, 0.0),
        (4.0, 0.0, 1.0),
        (-0.5, -1.0, 0.0),
        (-1.0, 0.0, -1.0),
        (-1.5, 1.0, 0.0),
        (-2.0, 0.0, 1.0),
        (100.0, 0.0, 1.0),
        (101.0, 0.0, -1.0),
        (100.5, 1.0, 0.0),
        (101.5, -1.0, 0.0),
    ];

    for &(x, exp_sin, exp_cos) in exact_cases {
        let (s, c) = sin_cos_pi(x);
        assert_eq!(s, exp_sin, "sin_pi exactness failed for {x}");
        assert_eq!(c, exp_cos, "cos_pi exactness failed for {x}");
        assert_eq!(sin_pi(x), exp_sin, "sin_pi failed for {x}");
        assert_eq!(cos_pi(x), exp_cos, "cos_pi failed for {x}");
        assert_eq!(
            cis_pi(x),
            Complex::new(exp_cos, exp_sin),
            "cis_pi failed for {x}"
        );
    }
}

#[test]
fn test_sin_cos_pi_intermediate_values() {
    let sqrt2_over_2 = std::f64::consts::FRAC_1_SQRT_2;
    let intermediate: &[(f64, f64, f64)] = &[
        (0.25, sqrt2_over_2, sqrt2_over_2),
        (0.75, sqrt2_over_2, -sqrt2_over_2),
        (1.25, -sqrt2_over_2, -sqrt2_over_2),
        (1.75, -sqrt2_over_2, sqrt2_over_2),
        (-0.25, -sqrt2_over_2, sqrt2_over_2),
        (-0.75, -sqrt2_over_2, -sqrt2_over_2),
    ];

    for &(x, exp_sin, exp_cos) in intermediate {
        let (s, c) = sin_cos_pi(x);
        assert_relative_eq!(s, exp_sin, epsilon = 1e-15);
        assert_relative_eq!(c, exp_cos, epsilon = 1e-15);
        let cis = cis_pi(x);
        assert_relative_eq!(cis.re, exp_cos, epsilon = 1e-15);
        assert_relative_eq!(cis.im, exp_sin, epsilon = 1e-15);
    }
}

#[test]
fn test_from_polar_pi() {
    let z = from_polar_pi(2.5_f64, 0.5);
    assert_eq!(z, Complex::new(0.0, 2.5));

    let z_neg = from_polar_pi(3.0_f64, 1.0);
    assert_eq!(z_neg, Complex::new(-3.0, 0.0));

    let z_diag = from_polar_pi(2.0_f64, 0.25);
    assert_relative_eq!(z_diag.re, 2.0 * std::f64::consts::FRAC_1_SQRT_2, epsilon = 1e-15);
    assert_relative_eq!(z_diag.im, 2.0 * std::f64::consts::FRAC_1_SQRT_2, epsilon = 1e-15);
}

