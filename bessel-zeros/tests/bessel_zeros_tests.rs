use amos_bessel_rs::{bessel_j, bessel_y};
use approx::assert_relative_eq;
use bessel_zeros::{fast, BesselFunType, bessel_zeros};
use rstest::rstest;

mod common;
use common::*;

// ---- helpers ----------------------------------------------------------------

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

// ---- hard-coded zeros: both backends ----------------------------------------

#[rstest]
#[case(BesselFunType::J, parse_zeros_wa(J_ZEROS))]
#[case(BesselFunType::JP, parse_zeros_wa(JP_ZEROS))]
#[case(BesselFunType::Y, parse_zeros_python(Y_ZEROS))]
#[case(BesselFunType::YP, parse_zeros_python(YP_ZEROS))]
fn test_against_hard_coded_zeros(#[case] fun_type: BesselFunType, #[case] zeros: Vec<Vec<f64>>) {
    for (order, expected) in (0_i32..).zip(zeros) {
        // AMOS backend
        let actual = bessel_zeros(fun_type, order as f64, expected.len(), 1e-6);
        expected
            .iter()
            .zip(&actual)
            .for_each(|(ev, av)| assert_relative_eq!(ev, av, epsilon = 5e-4));

        // fast (real-bessel) backend
        let actual_fast = fast::bessel_zeros(fun_type, order, expected.len(), 1e-6);
        expected
            .into_iter()
            .zip(actual_fast)
            .for_each(|(ev, av)| assert_relative_eq!(ev, av, epsilon = 5e-4));
    }
}

// ---- evaluate at returned zeros: J and Y, integer orders, both backends -----

#[rstest]
fn test_evaluation_at_zero_j_integer() {
    for order in 0..20 {
        for (zeros, backend) in [
            (bessel_zeros(BesselFunType::J, order as f64, 100, 0.1e-6), "amos"),
            (fast::bessel_zeros(BesselFunType::J, order, 100, 0.1e-6), "fast"),
        ] {
            for v in zeros {
                let val = bessel_j(order, v).unwrap();
                assert!(
                    approx::relative_eq!(0.0, val, epsilon = 1e-6),
                    "{backend} backend: J_{order}({v}) = {val}, expected ≈ 0"
                );
            }
        }
    }
}

#[rstest]
fn test_evaluation_at_zero_y_integer() {
    for order in 0..20 {
        for (zeros, backend) in [
            (bessel_zeros(BesselFunType::Y, order as f64, 100, 0.1e-6), "amos"),
            (fast::bessel_zeros(BesselFunType::Y, order, 100, 0.1e-6), "fast"),
        ] {
            for v in zeros {
                let val = bessel_y(order, v).unwrap();
                assert!(
                    approx::relative_eq!(0.0, val, epsilon = 1e-6),
                    "{backend} backend: Y_{order}({v}) = {val}, expected ≈ 0"
                );
            }
        }
    }
}

// ---- evaluate at returned zeros: JP and YP, integer orders, both backends ---

#[rstest]
fn test_evaluation_at_zero_jp_integer() {
    for order in 0..20 {
        for (zeros, backend) in [
            (bessel_zeros(BesselFunType::JP, order as f64, 100, 0.1e-6), "amos"),
            (fast::bessel_zeros(BesselFunType::JP, order, 100, 0.1e-6), "fast"),
        ] {
            for v in zeros {
                // J_0'(0) = 0 exactly; skip x=0 to avoid the n/x term.
                if v == 0.0 {
                    continue;
                }
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
fn test_evaluation_at_zero_yp_integer() {
    for order in 0..20 {
        for (zeros, backend) in [
            (bessel_zeros(BesselFunType::YP, order as f64, 100, 0.1e-6), "amos"),
            (fast::bessel_zeros(BesselFunType::YP, order, 100, 0.1e-6), "fast"),
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
        let zeros = bessel_zeros(BesselFunType::J, order, 100, 0.1e-6);
        for v in zeros {
            assert_relative_eq!(0.0, bessel_j(order, v).unwrap(), epsilon = 1e-6)
        }
    }
}

#[rstest]
fn test_evaluation_at_zero_y_float() {
    for order_int in 0..200 {
        let order = order_int as f64 / 50.0;
        let zeros = bessel_zeros(BesselFunType::Y, order, 100, 0.1e-6);
        for v in zeros {
            assert_relative_eq!(0.0, bessel_y(order, v).unwrap(), epsilon = 1e-6)
        }
    }
}

// ---- edge cases -------------------------------------------------------------

#[rstest]
fn test_n_zeros_zero_returns_empty() {
    for fun_type in [BesselFunType::J, BesselFunType::Y, BesselFunType::JP, BesselFunType::YP] {
        assert_eq!(bessel_zeros(fun_type, 0.0, 0, 1e-14), vec![]);
        assert_eq!(fast::bessel_zeros(fun_type, 0, 0, 1e-14), vec![]);
    }
}
