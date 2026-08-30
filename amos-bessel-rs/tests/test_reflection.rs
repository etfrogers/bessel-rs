use amos_bessel_rs::{
    BesselError, HankelKind, Scaling,
    amos::{
        complex_bessel_i, complex_bessel_j, complex_bessel_k, complex_bessel_y, complex_hankel1,
        complex_hankel2,
    },
    bessel_i, bessel_j, bessel_k, bessel_y, hankel,
};
use num::Complex;
use rstest::rstest;
use std::{assert_matches, f64::consts::PI};

mod common;
use common::assert_complex_arrays_equal;

use crate::common::BesselSig;

// const ORDER: [f64; 4] = [0.0, 5.0, -100.0, -10.0];

const Z_PARTS: [f64; 37] = [
    // -1e-6,
    // 0.0,
    -50.0, -40.0, -30.0, -25.0, -20.0, -15.0, -12.0, -10.0, -8.0, -6.0, -4.0, -3.0, -2.0, -1.0,
    -0.5, -0.1, -0.001, -1e-6, 0.0, 1e-6, 0.001, 0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0,
    12.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0,
];

#[rstest]
fn test_reflection_n_vs_loop(
    #[values(Scaling::Unscaled, Scaling::Scaled)] scaling: Scaling,
    #[values(0.0, 5.0, -3.0, -49.0, -51.0, -100.0, -10.763, -51.1276, 34.23)] order: f64,
    #[values(
        complex_bessel_j as BesselSig,
        // complex_bessel_i as BesselSig,
        // complex_hankel1 as BesselSig,
        // complex_hankel2 as BesselSig,
        // complex_bessel_k as BesselSig,
        // complex_bessel_y as BesselSig,
    )]
    fun: BesselSig,
) {
    let n = 50;
    for re in Z_PARTS {
        for im in Z_PARTS {
            let z = Complex::new(re, im);

            let result = fun(z, order, scaling, n);
            if matches!(result, Err(BesselError::Overflow)) {
                continue;
            }
            if z == Complex::ZERO && order.fract() != 0.0 && order < 0.0 {
                assert_matches!(result, Err(BesselError::InvalidInput { .. }));
                continue;
            }
            let (y, n_zeros) = result.unwrap();
            let mut sum_looped_nz = 0;

            for (i, yi) in y.iter().enumerate() {
                let current_order = order + i as f64;
                let (looped_y, looped_nz) = fun(z, current_order, scaling, 1).unwrap();
                assert_complex_arrays_equal(yi, &looped_y[0], &vec![], 1e6);
                sum_looped_nz += looped_nz;
            }

            // 1. Total underflow count agreement with looped evaluation (allowing +/- 1 boundary tolerance)
            let diff = (n_zeros as isize - sum_looped_nz as isize).abs();
            assert!(
                diff <= 1,
                "n_zeros mismatch: seq n_zeros={n_zeros}, sum looped={sum_looped_nz}, diff={diff} for order={order}, z={z}"
            );

            // 3. Count of exact zeros upper bound
            let count_exact_zeros = y.iter().filter(|&&v| v == Complex::ZERO).count();
            if n_zeros > 0 {
                assert!(
                    count_exact_zeros + 1 >= n_zeros,
                    "Vector has fewer exact zeros ({count_exact_zeros}) than n_zeros ({n_zeros}) for order={order}, z={z}"
                );
            }

            // 4. Total underflow state check
            if n_zeros == n {
                assert!(
                    y.iter().all(|&v| v == Complex::ZERO),
                    "n_zeros == n, but not all elements in y are zero: {:?}",
                    y
                );
            }
        }
    }
}

/// 1. Tests exact spherical/half-integer closed forms for negative orders.
/// These formulas evaluate directly with trigonometric and hyperbolic functions,
/// completely independent of any reflection formulas.
#[test]
fn test_half_integer_closed_forms() {
    let test_points = [
        Complex::new(1.5, 2.0),
        Complex::new(3.0, -1.5),
        Complex::new(-2.0, 3.0),
        Complex::new(-1.0, -2.5),
        Complex::new(5.0, 0.2),
        Complex::new(0.8, 0.5),
        Complex::new(-4.0, 0.1),
    ];

    for &z in &test_points {
        let factor = (Complex::new(2.0 / PI, 0.0) / z).sqrt();
        let k_factor = (Complex::new(PI / 2.0, 0.0) / z).sqrt();

        // J_{-1/2}(z) = sqrt(2 / (pi * z)) * cos(z)
        let actual_j_neg_half: Complex<f64> = bessel_j(-0.5, z).unwrap();
        let expected_j_neg_half = factor * z.cos();
        assert_complex_arrays_equal(&actual_j_neg_half, &expected_j_neg_half, &vec![], 1e6);

        // J_{-3/2}(z) = -sqrt(2 / (pi * z)) * (sin(z) + cos(z)/z)
        let actual_j_neg_three_halves: Complex<f64> = bessel_j(-1.5, z).unwrap();
        let expected_j_neg_three_halves = -factor * (z.sin() + z.cos() / z);
        assert_complex_arrays_equal(
            &actual_j_neg_three_halves,
            &expected_j_neg_three_halves,
            &vec![],
            1e6,
        );

        // Y_{-1/2}(z) = sqrt(2 / (pi * z)) * sin(z)
        let actual_y_neg_half: Complex<f64> = bessel_y(-0.5, z).unwrap();
        let expected_y_neg_half = factor * z.sin();
        assert_complex_arrays_equal(&actual_y_neg_half, &expected_y_neg_half, &vec![], 1e6);

        // Y_{-3/2}(z) = sqrt(2 / (pi * z)) * (cos(z) - sin(z)/z)
        let actual_y_neg_three_halves: Complex<f64> = bessel_y(-1.5, z).unwrap();
        let expected_y_neg_three_halves = factor * (z.cos() - z.sin() / z);
        assert_complex_arrays_equal(
            &actual_y_neg_three_halves,
            &expected_y_neg_three_halves,
            &vec![],
            1e6,
        );

        // I_{-1/2}(z) = sqrt(2 / (pi * z)) * cosh(z)
        let actual_i_neg_half: Complex<f64> = bessel_i(-0.5, z).unwrap();
        let expected_i_neg_half = factor * z.cosh();
        assert_complex_arrays_equal(&actual_i_neg_half, &expected_i_neg_half, &vec![], 1e6);

        // I_{-3/2}(z) = sqrt(2 / (pi * z)) * (sinh(z) - cosh(z)/z)
        let actual_i_neg_three_halves: Complex<f64> = bessel_i(-1.5, z).unwrap();
        let expected_i_neg_three_halves = factor * (z.sinh() - z.cosh() / z);
        assert_complex_arrays_equal(
            &actual_i_neg_three_halves,
            &expected_i_neg_three_halves,
            &vec![],
            1e6,
        );

        // K_{-1/2}(z) = K_{1/2}(z) = sqrt(pi / (2 * z)) * exp(-z)
        let actual_k_neg_half: Complex<f64> = bessel_k(-0.5, z).unwrap();
        let expected_k_neg_half = k_factor * (-z).exp();
        assert_complex_arrays_equal(&actual_k_neg_half, &expected_k_neg_half, &vec![], 1e6);

        // K_{-3/2}(z) = K_{3/2}(z) = sqrt(pi / (2 * z)) * exp(-z) * (1 + 1/z)
        let actual_k_neg_three_halves: Complex<f64> = bessel_k(-1.5, z).unwrap();
        let expected_k_neg_three_halves =
            k_factor * (-z).exp() * (Complex::new(1.0, 0.0) + Complex::new(1.0, 0.0) / z);
        assert_complex_arrays_equal(
            &actual_k_neg_three_halves,
            &expected_k_neg_three_halves,
            &vec![],
            1e6,
        );

        // H_{-1/2}^{(1)}(z) = sqrt(2 / (pi * z)) * exp(i * z)
        let actual_h1_neg_half: Complex<f64> = hankel(-0.5, z, HankelKind::First).unwrap();
        let expected_h1_neg_half = factor * (Complex::<f64>::i() * z).exp();
        assert_complex_arrays_equal(&actual_h1_neg_half, &expected_h1_neg_half, &vec![], 1e6);

        // H_{-1/2}^{(2)}(z) = sqrt(2 / (pi * z)) * exp(-i * z)
        let actual_h2_neg_half: Complex<f64> = hankel(-0.5, z, HankelKind::Second).unwrap();
        let expected_h2_neg_half = factor * (-Complex::<f64>::i() * z).exp();
        assert_complex_arrays_equal(&actual_h2_neg_half, &expected_h2_neg_half, &vec![], 1e6);
    }
}

/// 2. Tests the three-term recurrence relations across negative orders.
/// Recurrence relations hold for all real nu, providing an invariant test
/// that does not depend on reflection formulas.
#[test]
fn test_three_term_recurrence() {
    let test_points = [
        Complex::new(2.5, 1.5),
        Complex::new(4.0, -2.0),
        Complex::new(-3.0, 2.5),
        Complex::new(-2.0, -3.0),
        Complex::new(6.0, 0.0),
    ];

    let test_orders = [-0.5, -1.3, -2.7, -4.0, -5.2, -7.5];

    for &nu in &test_orders {
        for &z in &test_points {
            let two_nu_over_z = Complex::new(2.0 * nu, 0.0) / z;

            // J: J_{nu-1}(z) + J_{nu+1}(z) = (2nu / z) * J_nu(z)
            let j_minus: Complex<f64> = bessel_j(nu - 1.0, z).unwrap();
            let j_mid: Complex<f64> = bessel_j(nu, z).unwrap();
            let j_plus: Complex<f64> = bessel_j(nu + 1.0, z).unwrap();
            let lhs_j = j_minus + j_plus;
            let rhs_j = two_nu_over_z * j_mid;
            assert_complex_arrays_equal(&lhs_j, &rhs_j, &vec![], 1e6);

            // Y: Y_{nu-1}(z) + Y_{nu+1}(z) = (2nu / z) * Y_nu(z)
            let y_minus: Complex<f64> = bessel_y(nu - 1.0, z).unwrap();
            let y_mid: Complex<f64> = bessel_y(nu, z).unwrap();
            let y_plus: Complex<f64> = bessel_y(nu + 1.0, z).unwrap();
            let lhs_y = y_minus + y_plus;
            let rhs_y = two_nu_over_z * y_mid;
            assert_complex_arrays_equal(&lhs_y, &rhs_y, &vec![], 1e6);

            // I: I_{nu-1}(z) - I_{nu+1}(z) = (2nu / z) * I_nu(z)
            let i_minus: Complex<f64> = bessel_i(nu - 1.0, z).unwrap();
            let i_mid: Complex<f64> = bessel_i(nu, z).unwrap();
            let i_plus: Complex<f64> = bessel_i(nu + 1.0, z).unwrap();
            let lhs_i = i_minus - i_plus;
            let rhs_i = two_nu_over_z * i_mid;
            assert_complex_arrays_equal(&lhs_i, &rhs_i, &vec![], 1e6);

            // K: K_{nu-1}(z) - K_{nu+1}(z) = -(2nu / z) * K_nu(z)
            let k_minus: Complex<f64> = bessel_k(nu - 1.0, z).unwrap();
            let k_mid: Complex<f64> = bessel_k(nu, z).unwrap();
            let k_plus: Complex<f64> = bessel_k(nu + 1.0, z).unwrap();
            let lhs_k = k_minus - k_plus;
            let rhs_k = -two_nu_over_z * k_mid;
            assert_complex_arrays_equal(&lhs_k, &rhs_k, &vec![], 1e6);

            // Hankel 1 & 2: H_{nu-1}(z) + H_{nu+1}(z) = (2nu / z) * H_nu(z)
            for &kind in &[HankelKind::First, HankelKind::Second] {
                let h_minus: Complex<f64> = hankel(nu - 1.0, z, kind).unwrap();
                let h_mid: Complex<f64> = hankel(nu, z, kind).unwrap();
                let h_plus: Complex<f64> = hankel(nu + 1.0, z, kind).unwrap();
                let lhs_h = h_minus + h_plus;
                let rhs_h = two_nu_over_z * h_mid;
                assert_complex_arrays_equal(&lhs_h, &rhs_h, &vec![], 1e6);
            }
        }
    }
}

/// 3. Tests integer limits and continuity near negative integers.
/// Compares the negative integer identities with positive evaluations and
/// checks smooth continuity for nu = -n +/- eps.
#[test]
fn test_integer_limit_continuity() {
    let test_points = [
        Complex::new(2.0, 1.0),
        Complex::new(4.5, -2.5),
        Complex::new(-3.0, 2.0),
        Complex::new(-1.5, -3.5),
        Complex::new(5.0, 0.0),
    ];

    let integer_orders = [1, 2, 3, 5, 8];
    let eps = 1e-7;

    for &n in &integer_orders {
        let n_f64 = n as f64;
        let sign = if n % 2 == 0 { 1.0 } else { -1.0 };

        for &z in &test_points {
            // J_{-n}(z) = (-1)^n * J_n(z)
            let j_neg: Complex<f64> = bessel_j(-n_f64, z).unwrap();
            let j_pos: Complex<f64> = bessel_j(n_f64, z).unwrap();
            assert_complex_arrays_equal(&j_neg, &(j_pos * sign), &vec![], 1e6);

            // Continuity check: J_{-n + eps}(z) should be close to J_{-n}(z)
            let j_near: Complex<f64> = bessel_j(-n_f64 + eps, z).unwrap();
            let diff_j = (j_near - j_neg).norm();
            let scale_j = 1.0 + j_neg.norm() + j_near.norm();
            assert!(
                diff_j / scale_j < 1e-2,
                "J continuity failed at n={n}, z={z}: rel_diff={}",
                diff_j / scale_j
            );

            // Y_{-n}(z) = (-1)^n * Y_n(z)
            let y_neg: Complex<f64> = bessel_y(-n_f64, z).unwrap();
            let y_pos: Complex<f64> = bessel_y(n_f64, z).unwrap();
            assert_complex_arrays_equal(&y_neg, &(y_pos * sign), &vec![], 1e6);

            // Continuity check for Y
            let y_near: Complex<f64> = bessel_y(-n_f64 + eps, z).unwrap();
            let diff_y = (y_near - y_neg).norm();
            let scale_y = 1.0 + y_neg.norm() + y_near.norm();
            assert!(
                diff_y / scale_y < 1e-2,
                "Y continuity failed at n={n}, z={z}: rel_diff={}",
                diff_y / scale_y
            );

            // I_{-n}(z) = I_n(z)
            let i_neg: Complex<f64> = bessel_i(-n_f64, z).unwrap();
            let i_pos: Complex<f64> = bessel_i(n_f64, z).unwrap();
            assert_complex_arrays_equal(&i_neg, &i_pos, &vec![], 1e6);

            // Continuity check for I
            let i_near: Complex<f64> = bessel_i(-n_f64 + eps, z).unwrap();
            let diff_i = (i_near - i_neg).norm();
            let scale_i = 1.0 + i_neg.norm() + i_near.norm();
            assert!(
                diff_i / scale_i < 1e-2,
                "I continuity failed at n={n}, z={z}: rel_diff={}",
                diff_i / scale_i
            );

            // K_{-n}(z) = K_n(z)
            let k_neg: Complex<f64> = bessel_k(-n_f64, z).unwrap();
            let k_pos: Complex<f64> = bessel_k(n_f64, z).unwrap();
            assert_complex_arrays_equal(&k_neg, &k_pos, &vec![], 1e6);

            // Continuity check for K
            let k_near: Complex<f64> = bessel_k(-n_f64 + eps, z).unwrap();
            let diff_k = (k_near - k_neg).norm();
            let scale_k = 1.0 + k_neg.norm() + k_near.norm();
            assert!(
                diff_k / scale_k < 1e-2,
                "K continuity failed at n={n}, z={z}: rel_diff={}",
                diff_k / scale_k
            );

            // H_{-n}^{(m)}(z) = (-1)^n * H_n^{(m)}(z)
            for &kind in &[HankelKind::First, HankelKind::Second] {
                let h_neg: Complex<f64> = hankel(-n_f64, z, kind).unwrap();
                let h_pos: Complex<f64> = hankel(n_f64, z, kind).unwrap();
                assert_complex_arrays_equal(&h_neg, &(h_pos * sign), &vec![], 1e6);

                // Continuity check for H
                let h_near: Complex<f64> = hankel(-n_f64 + eps, z, kind).unwrap();
                let diff_h = (h_near - h_neg).norm();
                let scale_h = 1.0 + h_neg.norm() + h_near.norm();
                assert!(
                    diff_h / scale_h < 1e-2,
                    "H continuity failed at n={n}, z={z}: rel_diff={}",
                    diff_h / scale_h
                );
            }
        }
    }
}
