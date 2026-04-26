mod parametrisation;
use parametrisation::bessel_cases;

use approx::assert_relative_eq;
use bessel_rs::{BesselError, bessel_j, bessel_y};
use real_bessel::{bessel_j0, bessel_j1, bessel_jn, bessel_y0, bessel_y1};
use rstest::rstest;
use rstest_reuse::apply;

#[rstest]
#[case(0.0, 1.0)]
#[case(1.0, 0.7651976865579666)]
#[case(2.0, 0.22389077914123562)]
#[case(3.0, -0.2600519549019335)]
#[case(4.0, -0.39714980986384737)]
#[case(5.0, -0.1775967713143383)]
fn test_bessel_j0(#[case] input: f64, #[case] expected: f64) {
    let result = bessel_j0(input);
    assert_relative_eq!(result, expected, epsilon = 1e-10,);
}

#[apply(bessel_cases)]
fn test_bessel_j0_against_amos(#[case] _order: f64, #[case] zr: f64, #[case] _zi: f64) {
    let result = bessel_j(0.0, zr);
    let expected = match result {
        Ok(v) => v,
        Err(bessel_rs::BesselError::PartialLossOfSignificance { y, nz }) => {
            assert_eq!(nz, 0);
            y[0].re
        }
        Err(e) => panic!("Unexpected error from bessel_j: {:?}", e),
    };
    let actual = bessel_j0(zr);
    assert_relative_eq!(actual, expected, epsilon = 1e-10,);
}

#[apply(bessel_cases)]
fn test_bessel_y0_against_amos(#[case] _order: f64, #[case] mut zr: f64, #[case] _zi: f64) {
    if zr == 0.0 {
        // bessel_y0 returns -inf for 0, while bessel_y returns an error.
        return;
    }

    let actual_result = bessel_y0(zr);
    // bessel_y0 does not support negative inputs, but bessel_y does.
    // therefore if we get an invalid input error, we retry with the absolute value.
    // we could ignore this case, but it's better to test something.
    let actual = match actual_result {
        Ok(v) => v,
        Err(details) => {
            if details.contains("z < 0") {
                zr = zr.abs();
                bessel_y0(zr).unwrap()
            } else {
                panic!("Unexpected error from bessel_y0: {:?}", details);
            }
        }
    };

    let result = bessel_y(0.0, zr);
    let expected = match result {
        Ok(v) => v,
        Err(BesselError::PartialLossOfSignificance { y, nz }) => {
            assert_eq!(nz, 0);
            y[0].re
        }
        Err(e) => panic!("Unexpected error from bessel_y: {:?}", e),
    };

    assert_relative_eq!(actual, expected, epsilon = 1e-10,);
}

#[apply(bessel_cases)]
fn test_bessel_j1_against_amos(#[case] _order: f64, #[case] zr: f64, #[case] _zi: f64) {
    let result = bessel_j(1.0, zr);
    let expected = match result {
        Ok(v) => v,
        Err(BesselError::PartialLossOfSignificance { y, nz }) => {
            assert_eq!(nz, 0);
            y[0].re
        }
        Err(e) => panic!("Unexpected error from bessel_j: {:?}", e),
    };
    let actual = bessel_j1(zr);
    assert_relative_eq!(actual, expected, epsilon = 1e-10,);
}

#[apply(bessel_cases)]
fn test_bessel_y1_against_amos(#[case] _order: f64, #[case] mut zr: f64, #[case] _zi: f64) {
    if zr == 0.0 {
        // bessel_y1 returns -inf for 0, while bessel_y returns an error.
        return;
    }

    let actual_result = bessel_y1(zr);
    // bessel_y1 does not support negative inputs, but bessel_y does.
    // therefore if we get an invalid input error, we retry with the absolute value.
    // we could ignore this case, but it's better to test something.
    let actual = match actual_result {
        Ok(v) => v,
        Err(details) => {
            if details.contains("z < 0") {
                zr = zr.abs();
                bessel_y1(zr).unwrap()
            } else {
                panic!("Unexpected error from bessel_y1: {:?}", details);
            }
        }
    };

    let result = bessel_y(1.0, zr);
    let expected = match result {
        Ok(v) => v,
        Err(bessel_rs::BesselError::PartialLossOfSignificance { y, nz }) => {
            assert_eq!(nz, 0);
            y[0].re
        }
        Err(e) => panic!("Unexpected error from bessel_y: {:?}", e),
    };

    assert_relative_eq!(actual, expected, epsilon = 1e-10,);
}

#[apply(bessel_cases)]
#[trace]
fn test_bessel_jn_against_amos(
    #[case] order: f64,
    #[case] zr: f64,
    #[case] _zi: f64,
    // #[values(2, 3, 4, 5, 6, 7, 8, 9, 16, 125, 120923)] n: i32,
) {
    let integer_order = order as i32;
    let result = bessel_j(integer_order as f64, zr);
    let expected = match result {
        Ok(v) => v,
        Err(bessel_rs::BesselError::PartialLossOfSignificance { y, nz: _ }) => {
            // assert_eq!(nz, 0);
            y[0].re
        }
        Err(e) => panic!("Unexpected error from bessel_j: {:?}", e),
    };
    let actual = bessel_jn(integer_order, zr);
    assert_relative_eq!(actual, expected, epsilon = 1e-10,);
}
