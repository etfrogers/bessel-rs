use crate::BesselError;
use crate::go::bessel_j0_real;
use crate::tests::bessel_cases;
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
    let result = bessel_j0_real(input);
    let diff = (result - expected).abs();
    assert!(
        diff < 1e-10,
        "bessel_j0_real({}) = {}, expected {}, diff {}",
        input,
        result,
        expected,
        diff
    );
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
    let result = bessel_j0_real(zr);
    let diff = (result - expected).abs();
    assert!(
        diff < 1e-10,
        "bessel_j0_real({}) = {}, expected {}, diff {}",
        zr,
        result,
        expected,
        diff
    );
}
