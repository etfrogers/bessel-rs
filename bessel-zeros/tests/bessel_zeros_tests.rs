use amos_bessel_rs::{bessel_j, bessel_y};
use approx::assert_relative_eq;
use bessel_zeros::{BesselFunType, bessel_zeros};
use rstest::rstest;

mod common;
use common::*;

#[rstest]
#[case(BesselFunType::J, parse_zeros_wa(J_ZEROS))]
#[case(BesselFunType::JP, parse_zeros_wa(JP_ZEROS))]
#[case(BesselFunType::Y, parse_zeros_python(Y_ZEROS))]
#[case(BesselFunType::YP, parse_zeros_python(YP_ZEROS))]
fn test_against_hard_coded_zeros(#[case] fun_type: BesselFunType, #[case] zeros: Vec<Vec<f64>>) {
    for (order, expected) in (0_i32..).zip(zeros) {
        let actual = bessel_zeros(&fun_type, order, expected.len(), 1e-6);
        expected
            .into_iter()
            .zip(actual)
            .for_each(|(ev, av)| assert_relative_eq!(ev, av, epsilon = 5e-4));
    }
}

#[rstest]
fn test_evaluation_at_zero_j_integer() {
    for order in 0..20 {
        let zeros = bessel_zeros(&BesselFunType::J, order, 100, 0.1e-6);
        for v in zeros {
            assert_relative_eq!(0.0, bessel_j(order, v).unwrap(), epsilon = 1e-6)
        }
    }
}

#[rstest]
fn test_evaluation_at_zero_y_integer() {
    for order in 0..20 {
        let zeros = bessel_zeros(&BesselFunType::Y, order, 100, 0.1e-6);
        for v in zeros {
            assert_relative_eq!(0.0, bessel_y(order, v).unwrap(), epsilon = 1e-6);
        }
    }
}

#[rstest]
fn test_evaluation_at_zero_j_float() {
    for order_int in 0..200 {
        let order = order_int as f64 / 50.0;
        let zeros = bessel_zeros(&BesselFunType::J, order, 100, 0.1e-6);
        for v in zeros {
            assert_relative_eq!(0.0, bessel_j(order, v).unwrap(), epsilon = 1e-6)
        }
    }
}

#[rstest]
fn test_evaluation_at_zero_y_float() {
    for order_int in 0..200 {
        let order = order_int as f64 / 50.0;
        let zeros = bessel_zeros(&BesselFunType::Y, order, 100, 0.1e-6);
        for v in zeros {
            assert_relative_eq!(0.0, bessel_y(order, v).unwrap(), epsilon = 1e-6);
        }
    }
}
