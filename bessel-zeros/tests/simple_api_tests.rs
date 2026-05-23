use amos_bessel_rs::{bessel_j, bessel_y};
use approx::assert_relative_eq;
use bessel_zeros::*;
use rstest::rstest;

mod common;
use common::*;

#[rstest]
fn test_specialized_apis_basic() {
    let n = 1;
    let order = 0.0;

    let j = bessel_zeros_j(order, n);
    let y = bessel_zeros_y(order, n);
    let jp = bessel_zeros_jp(order, n);
    let yp = bessel_zeros_yp(order, n);

    assert_eq!(j.len(), n);
    assert_eq!(y.len(), n);
    assert_eq!(jp.len(), n);
    assert_eq!(yp.len(), n);

    assert!((j[0] - 2.40482555769577).abs() < 1e-10);
}

#[rstest]
fn test_j_zeros_against_hard_coded() {
    let zeros = parse_zeros_wa(J_ZEROS);
    for (order, expected) in (0_i32..).zip(zeros) {
        let actual = bessel_zeros_j(order, expected.len());
        expected
            .into_iter()
            .zip(actual)
            .for_each(|(ev, av)| assert_relative_eq!(ev, av, epsilon = 5e-4));
    }
}

#[rstest]
fn test_y_zeros_against_hard_coded() {
    let zeros = parse_zeros_python(Y_ZEROS);
    for (order, expected) in (0_i32..).zip(zeros) {
        let actual = bessel_zeros_y(order, expected.len());
        expected
            .into_iter()
            .zip(actual)
            .for_each(|(ev, av)| assert_relative_eq!(ev, av, epsilon = 5e-4));
    }
}

#[rstest]
fn test_jp_zeros_against_hard_coded() {
    let zeros = parse_zeros_wa(JP_ZEROS);
    for (order, expected) in (0_i32..).zip(zeros) {
        let actual = bessel_zeros_jp(order, expected.len());
        expected
            .into_iter()
            .zip(actual)
            .for_each(|(ev, av)| assert_relative_eq!(ev, av, epsilon = 5e-4));
    }
}

#[rstest]
fn test_yp_zeros_against_hard_coded() {
    let zeros = parse_zeros_python(YP_ZEROS);
    for (order, expected) in (0_i32..).zip(zeros) {
        let actual = bessel_zeros_yp(order, expected.len());
        expected
            .into_iter()
            .zip(actual)
            .for_each(|(ev, av)| assert_relative_eq!(ev, av, epsilon = 5e-4));
    }
}

#[rstest]
fn test_evaluation_at_zero_j_int() {
    for order in 0..20 {
        let zeros = bessel_zeros_j(order, 50);
        for v in zeros {
            assert_relative_eq!(0.0, bessel_j(order, v).unwrap(), epsilon = 1e-6)
        }
    }
}

#[rstest]
fn test_evaluation_at_zero_y_int() {
    for order in 0..20 {
        let zeros = bessel_zeros_y(order, 50);
        for v in zeros {
            assert_relative_eq!(0.0, bessel_y(order, v).unwrap(), epsilon = 1e-6);
        }
    }
}

#[rstest]
fn test_evaluation_at_zero_j_float() {
    for order_int in 0..200 {
        let order = order_int as f64 / 50.0;
        let zeros = bessel_zeros_j(order, 50);
        for v in zeros {
            assert_relative_eq!(0.0, bessel_j(order, v).unwrap(), epsilon = 1e-6)
        }
    }
}

#[rstest]
fn test_evaluation_at_zero_y_float() {
    for order_int in 0..200 {
        let order = order_int as f64 / 50.0;
        let zeros = bessel_zeros_y(order, 50);
        for v in zeros {
            assert_relative_eq!(0.0, bessel_y(order, v).unwrap(), epsilon = 1e-6);
        }
    }
}
