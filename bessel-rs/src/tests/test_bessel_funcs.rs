use super::{
    BesselFortranSig, BesselSig, airy_ref, bessel_cases, bessel_h_ref, biry_ref,
    check_against_fortran, sig_airy, sig_airy_fortran, sig_airyp, sig_airyp_fortran, sig_biry,
    sig_biry_fortran, sig_biryp, sig_biryp_fortran, utils::assert_results_are_equal, zbesh_first,
    zbesh_fortran_first, zbesh_fortran_second, zbesh_second, zbesi_fortran, zbesj_fortran,
    zbesk_fortran, zbesy_fortran,
};

use complex_bessel_rs::bessel_i::bessel_i as bessel_i_ref;
use complex_bessel_rs::bessel_j::bessel_j as bessel_j_ref;
use complex_bessel_rs::bessel_k::bessel_k as bessel_k_ref;
use complex_bessel_rs::bessel_y::bessel_y as bessel_y_ref;

use crate::{
    HankelKind, Scaling, airy, airy_b, airy_bp, airyp,
    amos::{complex_bessel_i, complex_bessel_j, complex_bessel_k, complex_bessel_y},
    bessel_i, bessel_j, bessel_k, bessel_y, hankel,
};
use num::complex::Complex64;
use rstest::rstest;
use rstest_reuse::apply;

#[apply(bessel_cases)]
fn test_bessel_j(#[case] order: f64, #[case] zr: f64, #[case] zi: f64) {
    let z = Complex64::new(zr, zi);
    let actual = bessel_j(order, z);

    let expected = bessel_j_ref(order, z.into());
    assert_results_are_equal(&actual, &expected, &Vec::new());
}

#[apply(bessel_cases)]
#[trace]
fn test_bessel_i(#[case] order: f64, #[case] zr: f64, #[case] zi: f64) {
    let z = Complex64::new(zr, zi);
    let actual = bessel_i(order, z);
    // dbg!(&actual);

    let expected = bessel_i_ref(order, z.into());
    assert_results_are_equal(&actual, &expected, &Vec::new());
}

#[apply(bessel_cases)]
#[trace]
fn test_bessel_k(#[case] order: f64, #[case] zr: f64, #[case] zi: f64) {
    let z = Complex64::new(zr, zi);
    let actual = bessel_k(order, z);

    let expected = bessel_k_ref(order, z.into());
    assert_results_are_equal(&actual, &expected, &Vec::new());
}

#[apply(bessel_cases)]
#[trace]
fn test_bessel_y(#[case] order: f64, #[case] zr: f64, #[case] zi: f64) {
    let z = Complex64::new(zr, zi);
    let actual = bessel_y(order, z);

    let expected = bessel_y_ref(order, z.into());
    assert_results_are_equal(&actual, &expected, &Vec::new());
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
    assert_results_are_equal(&actual, &expected, &Vec::new());
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
    assert_results_are_equal(&actual, &expected, &Vec::new());
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
    assert_results_are_equal(&actual, &expected, &Vec::new());
}

#[rstest]
fn test_bessel_extremes(
    #[values(
        (complex_bessel_j as BesselSig, zbesj_fortran as BesselFortranSig),
        (complex_bessel_i as BesselSig, zbesi_fortran as BesselFortranSig),
        (zbesh_first as BesselSig , zbesh_fortran_first as BesselFortranSig),
        (zbesh_second as BesselSig , zbesh_fortran_second as BesselFortranSig),
        (complex_bessel_k as BesselSig, zbesk_fortran as BesselFortranSig),
        (complex_bessel_y as BesselSig, zbesy_fortran as BesselFortranSig),
        (sig_airy as BesselSig, sig_airy_fortran as BesselFortranSig),
        (sig_airyp as BesselSig, sig_airyp_fortran as BesselFortranSig),
        (sig_biry as BesselSig, sig_biry_fortran as BesselFortranSig),
        (sig_biryp as BesselSig, sig_biryp_fortran as BesselFortranSig),
    )]
    (rust_fn, fortran_fn): (BesselSig, BesselFortranSig),

    #[values(0.0, f64::EPSILON, f64::MIN_POSITIVE, 1.0, f64::MAX)] order: f64,
    #[values(
        f64::MIN,
        -f64::EPSILON,
        -f64::MIN_POSITIVE,
        0.0,
        f64::MIN_POSITIVE,
        f64::EPSILON,
        1.0,
        f64::MAX
    )]
    zr: f64,
    #[values(
        f64::MIN,
        -f64::EPSILON,
        -f64::MIN_POSITIVE,
        0.0,
        f64::MIN_POSITIVE,
        f64::EPSILON,
        1.0,
        f64::MAX
        )]
    zi: f64,
    #[values(Scaling::Unscaled, Scaling::Scaled)] scaling: Scaling,
) {
    let n = 9;
    let z = Complex64::new(zr, zi);
    check_against_fortran(order, z, scaling, n, rust_fn, fortran_fn);
}
