use amos_bessel_rs::{bessel_j, bessel_y};
use bessel_zeros::{fast, *};
use rstest::rstest;

mod common;
use common::*;

// ---- helpers ----------------------------------------------------------------

/// Run a hard-coded zeros check with any callable that takes `(i32, n) -> Vec<f64>`.
fn check_zeros_against_hard_coded(
    expected_zeros: Vec<Vec<f64>>,
    backend: &str,
    get_zeros: impl Fn(i32, usize) -> Vec<f64>,
) {
    for (order, expected) in (0_i32..).zip(expected_zeros) {
        let actual = get_zeros(order, expected.len());
        expected
            .into_iter()
            .zip(actual)
            .for_each(|(ev, av)| {
                assert!(
                    approx::relative_eq!(ev, av, epsilon = 5e-4),
                    "{backend} backend: expected {ev}, got {av} (order={order})"
                )
            });
    }
}

/// Verify returned J zeros actually satisfy J_ν(x) ≈ 0.
fn check_j_zeros_evaluate_to_zero(zeros: Vec<f64>, order: impl Into<f64> + Copy, backend: &str) {
    let order_f64 = order.into();
    for v in zeros {
        let val = bessel_j(order, v).unwrap();
        assert!(
            approx::relative_eq!(0.0, val, epsilon = 1e-6),
            "{backend} backend: J_{order_f64}({v}) = {val}, expected ≈ 0"
        );
    }
}

/// Verify returned Y zeros actually satisfy Y_ν(x) ≈ 0.
fn check_y_zeros_evaluate_to_zero(zeros: Vec<f64>, order: impl Into<f64> + Copy, backend: &str) {
    let order_f64 = order.into();
    for v in zeros {
        let val = bessel_y(order, v).unwrap();
        assert!(
            approx::relative_eq!(0.0, val, epsilon = 1e-6),
            "{backend} backend: Y_{order_f64}({v}) = {val}, expected ≈ 0"
        );
    }
}

/// J_n'(x) via the recurrence  J_n' = n/x * J_n - J_{n+1}  (x ≠ 0).
fn bessel_j_prime(order: i32, x: f64) -> f64 {
    let n = order as f64;
    n / x * bessel_j(order, x).unwrap() - bessel_j(order + 1, x).unwrap()
}

/// Y_n'(x) via the recurrence  Y_n' = n/x * Y_n - Y_{n+1}  (x > 0).
fn bessel_y_prime(order: i32, x: f64) -> f64 {
    let n = order as f64;
    n / x * bessel_y(order, x).unwrap() - bessel_y(order + 1, x).unwrap()
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
    check_zeros_against_hard_coded(zeros.clone(), "amos", |order, n| bessel_zeros_j(order, n));
    check_zeros_against_hard_coded(zeros, "fast", fast::bessel_zeros_j);
}

#[rstest]
fn test_y_zeros_against_hard_coded() {
    let zeros = parse_zeros_python(Y_ZEROS);
    check_zeros_against_hard_coded(zeros.clone(), "amos", |order, n| bessel_zeros_y(order, n));
    check_zeros_against_hard_coded(zeros, "fast", fast::bessel_zeros_y);
}

#[rstest]
fn test_jp_zeros_against_hard_coded() {
    let zeros = parse_zeros_wa(JP_ZEROS);
    check_zeros_against_hard_coded(zeros.clone(), "amos", |order, n| bessel_zeros_jp(order, n));
    check_zeros_against_hard_coded(zeros, "fast", fast::bessel_zeros_jp);
}

#[rstest]
fn test_yp_zeros_against_hard_coded() {
    let zeros = parse_zeros_python(YP_ZEROS);
    check_zeros_against_hard_coded(zeros.clone(), "amos", |order, n| bessel_zeros_yp(order, n));
    check_zeros_against_hard_coded(zeros, "fast", fast::bessel_zeros_yp);
}

// ---- evaluate at returned zeros: J and Y, integer orders, both backends -----

#[rstest]
fn test_evaluation_at_zero_j_int() {
    for order in 0..20 {
        check_j_zeros_evaluate_to_zero(bessel_zeros_j(order, 50), order, "amos");
        check_j_zeros_evaluate_to_zero(fast::bessel_zeros_j(order, 50), order, "fast");
    }
}

#[rstest]
fn test_evaluation_at_zero_y_int() {
    for order in 0..20 {
        check_y_zeros_evaluate_to_zero(bessel_zeros_y(order, 50), order, "amos");
        check_y_zeros_evaluate_to_zero(fast::bessel_zeros_y(order, 50), order, "fast");
    }
}

// ---- evaluate at returned zeros: JP and YP, integer orders, both backends ---

#[rstest]
fn test_evaluation_at_zero_jp_int() {
    for order in 0..20 {
        for (zeros, backend) in [
            (bessel_zeros_jp(order, 50), "amos"),
            (fast::bessel_zeros_jp(order, 50), "fast"),
        ] {
            for v in zeros {
                if v == 0.0 { continue; } // J_0'(0) = 0 exactly; skip to avoid n/x term
                let val = bessel_j_prime(order, v);
                assert!(
                    approx::relative_eq!(0.0, val, epsilon = 1e-6),
                    "{backend} backend: J_{order}'({v}) = {val}, expected ≈ 0"
                );
            }
        }
    }
}

#[rstest]
fn test_evaluation_at_zero_yp_int() {
    for order in 0..20 {
        for (zeros, backend) in [
            (bessel_zeros_yp(order, 50), "amos"),
            (fast::bessel_zeros_yp(order, 50), "fast"),
        ] {
            for v in zeros {
                let val = bessel_y_prime(order, v);
                assert!(
                    approx::relative_eq!(0.0, val, epsilon = 1e-6),
                    "{backend} backend: Y_{order}'({v}) = {val}, expected ≈ 0"
                );
            }
        }
    }
}

// ---- evaluate at returned zeros: float orders, AMOS backend only ------------

#[rstest]
fn test_evaluation_at_zero_j_float() {
    for order_int in 0..200 {
        let order = order_int as f64 / 50.0;
        let zeros = bessel_zeros_j(order, 50);
        check_j_zeros_evaluate_to_zero(zeros, order, "amos");
    }
}

#[rstest]
fn test_evaluation_at_zero_y_float() {
    for order_int in 0..200 {
        let order = order_int as f64 / 50.0;
        let zeros = bessel_zeros_y(order, 50);
        check_y_zeros_evaluate_to_zero(zeros, order, "amos");
    }
}

// ---- edge cases -------------------------------------------------------------

#[rstest]
fn test_n_zeros_zero_returns_empty() {
    assert_eq!(bessel_zeros_j(0, 0), vec![]);
    assert_eq!(bessel_zeros_y(0, 0), vec![]);
    assert_eq!(bessel_zeros_jp(0, 0), vec![]);
    assert_eq!(bessel_zeros_yp(0, 0), vec![]);

    assert_eq!(fast::bessel_zeros_j(0, 0), vec![]);
    assert_eq!(fast::bessel_zeros_y(0, 0), vec![]);
    assert_eq!(fast::bessel_zeros_jp(0, 0), vec![]);
    assert_eq!(fast::bessel_zeros_yp(0, 0), vec![]);
}
