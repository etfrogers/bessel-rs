use super::bessel_cases;

use complex_bessel_rs::bessel_i::bessel_i as bessel_i_ref;
use complex_bessel_rs::bessel_j::bessel_j as bessel_j_ref;
use complex_bessel_rs::bessel_k::bessel_k as bessel_k_ref;
use complex_bessel_rs::bessel_y::bessel_y as bessel_y_ref;

use crate::{
    BesselError, HankelKind, Scaling, airy, airy_b, airy_bp, airyp,
    amos::{complex_bessel_i, complex_bessel_j, complex_bessel_k, complex_bessel_y},
    bessel_i, bessel_j, bessel_k, bessel_y, complex_hankel1, complex_hankel2, hankel,
    test_utils::{
        BesselSig, ComplexConversions, DiagnosticBesselFloat, airy_ref, assert_results_are_equal,
        assert_results_are_equal_floats, bessel_h_ref, biry_ref, sig_airy, sig_airyp, sig_biry,
        sig_biryp, zeros_are_not_equivalent,
    },
    types::{BesselFloat, BesselValues},
};
use num::{Complex, complex::Complex64};
use rstest::rstest;
use rstest_reuse::apply;

fn single_to_bessel_values<T: BesselFloat>(val: Complex<T>) -> BesselValues<T> {
    (vec![val], 0)
}

fn assert_results_are_equal_fotran<T: DiagnosticBesselFloat>(
    actual: Result<Complex<T>, BesselError<T>>,
    expected: Result<Complex<f64>, i32>,
    margin: f64,
) {
    let actual = actual.map(single_to_bessel_values);
    let expected = expected
        .map_err(|val| BesselError::from_i32(val).unwrap())
        .map(single_to_bessel_values);
    if actual.is_ok()
        && matches!(
            expected,
            Err(BesselError::PartialLossOfSignificance { y: _, nz: _ })
        )
    {
        // the single output functions unwrap partial loss of significance, where
        // complex_bessel_rs returns an error code.
        return;
    }
    assert_results_are_equal_floats(&actual, &expected, margin);
}

#[apply(bessel_cases)]
#[rstest]
fn test_bessel_j(#[case] order: f64, #[case] zr: f64, #[case] zi: f64) {
    let z64 = Complex64::new(zr, zi);
    let z32 = z64.to_c32();
    if zeros_are_not_equivalent(z64, z32, order) {
        return;
    }
    let actual = bessel_j(order as f32, z32);

    let expected = bessel_j_ref(order, z64.into());
    assert_results_are_equal_fotran(actual, expected, 1e6)
}

#[apply(bessel_cases)]
#[rstest]
fn test_bessel_i(#[case] order: f64, #[case] zr: f64, #[case] zi: f64) {
    let z64 = Complex64::new(zr, zi);
    let z32 = z64.to_c32();
    if zeros_are_not_equivalent(z64, z32, order) {
        return;
    }
    let actual = bessel_i(order as f32, z32);

    let expected = bessel_i_ref(order, z64.into());
    assert_results_are_equal_fotran(actual, expected, 1e6)
}

#[apply(bessel_cases)]
#[trace]
fn test_bessel_k(#[case] order: f64, #[case] zr: f64, #[case] zi: f64) {
    let z = Complex64::new(zr, zi);
    let actual = bessel_k(order, z);

    let expected = bessel_k_ref(order, z.into());
    assert_results_are_equal(&actual, &expected, &Vec::new(), 1e6);
}

#[apply(bessel_cases)]
#[trace]
fn test_bessel_y(#[case] order: f64, #[case] zr: f64, #[case] zi: f64) {
    let z = Complex64::new(zr, zi);
    let actual = bessel_y(order, z);

    let expected = bessel_y_ref(order, z.into());
    assert_results_are_equal(&actual, &expected, &Vec::new(), 1e6);
}

#[apply(bessel_cases)]
#[trace]
fn test_bessel_h(
    #[case] order: f64,
    #[case] zr: f64,
    #[case] zi: f64,
    #[values(HankelKind::First, HankelKind::Second)] kind: HankelKind,
) {
    let z = Complex64::new(zr, zi);
    let actual = hankel(order, z, kind);
    // dbg!(&actual);
    let expected = bessel_h_ref(order, z.into(), kind);
    // dbg!(&expected);
    assert_results_are_equal(&actual, &expected, &Vec::new(), 1e6);
}

#[apply(bessel_cases)]
#[trace]
fn test_airy(
    #[case] _order: f64,
    #[case] zr: f64,
    #[case] zi: f64,
    #[values(false, true)] is_derivative: bool,
) {
    let z = Complex64::new(zr, zi);

    let actual = if is_derivative { airyp(z) } else { airy(z) };
    let expected = airy_ref(z, is_derivative);
    assert_results_are_equal(&actual, &expected, &Vec::new(), 1e6);
}

#[apply(bessel_cases)]
#[trace]
fn test_biry(
    #[case] _order: f64,
    #[case] zr: f64,
    #[case] zi: f64,
    #[values(false, true)] is_derivative: bool,
) {
    let z = Complex64::new(zr, zi);

    let actual = if is_derivative { airy_bp(z) } else { airy_b(z) };
    let expected = biry_ref(z, is_derivative);
    assert_results_are_equal(&actual, &expected, &Vec::new(), 1e6);
}

#[apply(bessel_cases)]
#[trace]
#[rstest]
fn test_bessel_funcs_f32(
    #[case] order: f64,
    #[case] zr: f64,
    #[case] zi: f64,
    #[values(Scaling::Unscaled, Scaling::Scaled)] scaling: Scaling,
    #[values(
        (complex_bessel_j as BesselSig<f64>, complex_bessel_j as BesselSig<f32>),
        (complex_bessel_i as BesselSig<f64>, complex_bessel_i as BesselSig<f32>),
        (complex_hankel1 as BesselSig<f64>, complex_hankel1 as BesselSig<f32>),
        (complex_hankel2 as BesselSig<f64>, complex_hankel2 as BesselSig<f32>),
        (complex_bessel_k as BesselSig<f64>, complex_bessel_k as BesselSig<f32>),
        (complex_bessel_y as BesselSig<f64>, complex_bessel_y as BesselSig<f32>),
        (sig_airy::<f64> as BesselSig<f64>, sig_airy::<f32> as BesselSig<f32>),
        (sig_airyp::<f64> as BesselSig<f64>, sig_airyp::<f32> as BesselSig<f32>),
        (sig_biry::<f64> as BesselSig<f64>, sig_biry::<f32> as BesselSig<f32>),
        (sig_biryp::<f64> as BesselSig<f64>, sig_biryp::<f32> as BesselSig<f32>),
    )]
    (fn_64, fn_32): (BesselSig<f64>, BesselSig<f32>),
) {
    let n = 1;
    let z64 = Complex64::new(zr, zi);
    let z32 = z64.to_c32();
    if zeros_are_not_equivalent(z64, z32, order) {
        return;
    }

    let f64_vals = fn_64(z64, order, scaling, n);
    let f32_vals = fn_32(z32, order as f32, scaling, n);

    assert_results_are_equal_floats(&f32_vals, &f64_vals, 1e4);
}
