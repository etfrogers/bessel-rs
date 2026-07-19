use amos_bessel_rs::{bessel_j, bessel_y};
use approx::assert_relative_eq;
use bessel_zeros::{fast, *};
use rstest::rstest;

mod common;
use common::*;

// ---- helpers ----------------------------------------------------------------

/// Run a hard-coded zeros check with any callable that takes `(i32, n) -> Vec<f64>`.
fn check_zeros_against_hard_coded(
    expected_zeros: Vec<Vec<f64>>,
    get_zeros: impl Fn(i32, usize) -> Vec<f64>,
) {
    for (order, expected) in (0_i32..).zip(expected_zeros) {
        let actual = get_zeros(order, expected.len());
        expected
            .into_iter()
            .zip(actual)
            .for_each(|(ev, av)| assert_relative_eq!(ev, av, epsilon = 5e-4));
    }
}

/// Verify that every returned zero is actually a zero of Jν, using the AMOS evaluator.
fn check_j_zeros_evaluate_to_zero(zeros: Vec<f64>, order: impl Into<f64> + Copy) {
    for v in zeros {
        assert_relative_eq!(0.0, bessel_j(order, v).unwrap(), epsilon = 1e-6);
    }
}

/// Verify that every returned zero is actually a zero of Yν, using the AMOS evaluator.
fn check_y_zeros_evaluate_to_zero(zeros: Vec<f64>, order: impl Into<f64> + Copy) {
    for v in zeros {
        assert_relative_eq!(0.0, bessel_y(order, v).unwrap(), epsilon = 1e-6);
    }
}

// ---- basic sanity -----------------------------------------------------------

#[rstest]
fn test_specialized_apis_basic() {
    let n = 1;
    let order = 0;

    // AMOS backend
    let j = bessel_zeros_j(order, n);
    let y = bessel_zeros_y(order, n);
    let jp = bessel_zeros_jp(order, n);
    let yp = bessel_zeros_yp(order, n);
    assert_eq!(j.len(), n);
    assert_eq!(y.len(), n);
    assert_eq!(jp.len(), n);
    assert_eq!(yp.len(), n);
    assert!((j[0] - 2.40482555769577).abs() < 1e-10);

    // fast backend
    let j = fast::bessel_zeros_j(order, n);
    let y = fast::bessel_zeros_y(order, n);
    let jp = fast::bessel_zeros_jp(order, n);
    let yp = fast::bessel_zeros_yp(order, n);
    assert_eq!(j.len(), n);
    assert_eq!(y.len(), n);
    assert_eq!(jp.len(), n);
    assert_eq!(yp.len(), n);
    assert!((j[0] - 2.40482555769577).abs() < 1e-10);
}

// ---- hard-coded zeros: both backends ----------------------------------------

#[rstest]
fn test_j_zeros_against_hard_coded() {
    let zeros = parse_zeros_wa(J_ZEROS);
    check_zeros_against_hard_coded(zeros.clone(), |order, n| bessel_zeros_j(order, n));
    check_zeros_against_hard_coded(zeros, fast::bessel_zeros_j);
}

#[rstest]
fn test_y_zeros_against_hard_coded() {
    let zeros = parse_zeros_python(Y_ZEROS);
    check_zeros_against_hard_coded(zeros.clone(), |order, n| bessel_zeros_y(order, n));
    check_zeros_against_hard_coded(zeros, fast::bessel_zeros_y);
}

#[rstest]
fn test_jp_zeros_against_hard_coded() {
    let zeros = parse_zeros_wa(JP_ZEROS);
    check_zeros_against_hard_coded(zeros.clone(), |order, n| bessel_zeros_jp(order, n));
    check_zeros_against_hard_coded(zeros, fast::bessel_zeros_jp);
}

#[rstest]
fn test_yp_zeros_against_hard_coded() {
    let zeros = parse_zeros_python(YP_ZEROS);
    check_zeros_against_hard_coded(zeros.clone(), |order, n| bessel_zeros_yp(order, n));
    check_zeros_against_hard_coded(zeros, fast::bessel_zeros_yp);
}

// ---- evaluate at returned zeros: integer orders, both backends --------------

#[rstest]
fn test_evaluation_at_zero_j_int() {
    for order in 0..20 {
        check_j_zeros_evaluate_to_zero(bessel_zeros_j(order, 50), order);
        check_j_zeros_evaluate_to_zero(fast::bessel_zeros_j(order, 50), order);
    }
}

#[rstest]
fn test_evaluation_at_zero_y_int() {
    for order in 0..20 {
        check_y_zeros_evaluate_to_zero(bessel_zeros_y(order, 50), order);
        check_y_zeros_evaluate_to_zero(fast::bessel_zeros_y(order, 50), order);
    }
}

// ---- evaluate at returned zeros: float orders, AMOS backend only ------------

#[rstest]
fn test_evaluation_at_zero_j_float() {
    for order_int in 0..200 {
        let order = order_int as f64 / 50.0;
        let zeros = bessel_zeros_j(order, 50);
        check_j_zeros_evaluate_to_zero(zeros, order);
    }
}

#[rstest]
fn test_evaluation_at_zero_y_float() {
    for order_int in 0..200 {
        let order = order_int as f64 / 50.0;
        let zeros = bessel_zeros_y(order, 50);
        check_y_zeros_evaluate_to_zero(zeros, order);
    }
}
