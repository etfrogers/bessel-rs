use crate::BesselError;
use crate::go::{bessel_j0_real, bessel_y0_real};
use crate::tests::bessel_cases;
use approx::assert_relative_eq;
use rstest::rstest;
use rstest_reuse::apply;

#[rstest]
#[case(0.0, 1.0)]
#[case(1.0, 0.7651976865579666)]
#[case(2.0, 0.22389077914123562)]
#[case(3.0, -0.2600519549019335)]
#[case(4.0, -0.39714980986384737)]
#[case(5.0, -0.1775967713143383)]
fn test_bessel_j0_real(#[case] input: f64, #[case] expected: f64) {
    let result = bessel_j0_real(input).unwrap();
    assert_relative_eq!(result, expected, epsilon = 1e-10,);
}

#[apply(bessel_cases)]
fn test_bessel_j0_against_amos(#[case] _order: f64, #[case] zr: f64, #[case] _zi: f64) {
    let result = crate::bessel_j(0.0, zr);
    let expected = match result {
        Ok(v) => v,
        Err(BesselError::PartialLossOfSignificance { y, nz }) => {
            assert_eq!(nz, 0);
            y[0].re
        }
        Err(e) => panic!("Unexpected error from bessel_j: {:?}", e),
    };
    let actual = bessel_j0_real(zr).unwrap();
    assert_relative_eq!(actual, expected, epsilon = 1e-10,);
}

#[apply(bessel_cases)]
fn test_bessel_y0_against_amos(#[case] _order: f64, #[case] mut zr: f64, #[case] _zi: f64) {
    if zr == 0.0 {
        // bessel_y0_real returns -inf for 0, while bessel_y returns an error.
        return;
    }

    let actual_result = bessel_y0_real(zr);
    // bessel_y0_real does not support negative inputs, but bessel_y does.
    // therefore if we get an invalid input error, we retry with the absolute value.
    // we could ignore this case, but it's better to test something.
    let actual = match actual_result {
        Ok(v) => v,
        Err(BesselError::InvalidInput { details }) => {
            if details.contains("z < 0") {
                zr = zr.abs();
                bessel_y0_real(zr).unwrap()
            } else {
                panic!("Unexpected error from bessel_y0_real: {:?}", details);
            }
        }
        Err(e) => panic!("Unexpected error from bessel_y0_real: {:?}", e),
    };

    let result = crate::bessel_y(0.0, zr);
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
    let result = crate::bessel_j(1.0, zr);
    let expected = match result {
        Ok(v) => v,
        Err(BesselError::PartialLossOfSignificance { y, nz }) => {
            assert_eq!(nz, 0);
            y[0].re
        }
        Err(e) => panic!("Unexpected error from bessel_j: {:?}", e),
    };
    let actual = crate::go::bessel_j1_real(zr).unwrap();
    assert_relative_eq!(actual, expected, epsilon = 1e-10,);
}

#[apply(bessel_cases)]
fn test_bessel_y1_against_amos(#[case] _order: f64, #[case] mut zr: f64, #[case] _zi: f64) {
    if zr == 0.0 {
        // bessel_y1_real returns -inf for 0, while bessel_y returns an error.
        return;
    }

    let actual_result = crate::go::bessel_y1_real(zr);
    // bessel_y1_real does not support negative inputs, but bessel_y does.
    // therefore if we get an invalid input error, we retry with the absolute value.
    // we could ignore this case, but it's better to test something.
    let actual = match actual_result {
        Ok(v) => v,
        Err(BesselError::InvalidInput { details }) => {
            if details.contains("z < 0") {
                zr = zr.abs();
                crate::go::bessel_y1_real(zr).unwrap()
            } else {
                panic!("Unexpected error from bessel_y1_real: {:?}", details);
            }
        }
        Err(e) => panic!("Unexpected error from bessel_y1_real: {:?}", e),
    };

    let result = crate::bessel_y(1.0, zr);
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
