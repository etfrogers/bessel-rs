use super::{bessel_cases, zbesi_fortran, zbesj_fortran, zbesk_fortran, zbesy_fortran};

use complex_bessel_rs::bessel_i::bessel_i as bessel_i_ref;
use complex_bessel_rs::bessel_j::bessel_j as bessel_j_ref;
use complex_bessel_rs::bessel_k::bessel_k as bessel_k_ref;
use complex_bessel_rs::bessel_y::bessel_y as bessel_y_ref;

use crate::{
    HankelKind, Scaling, airy, airy_b, airy_bp, airyp,
    amos::{complex_bessel_i, complex_bessel_j, complex_bessel_k, complex_bessel_y},
    bessel_i, bessel_j, bessel_k, bessel_y, complex_hankel1, complex_hankel2, hankel,
    test_utils::{
        BesselFortranSig, BesselSig, ComplexConversions, airy_ref, assert_results_are_equal,
        assert_results_are_equal_floats, bessel_h_ref, biry_ref, check_against_fortran, sig_airy,
        sig_airy_fortran, sig_airyp, sig_airyp_fortran, sig_biry, sig_biry_fortran, sig_biryp,
        sig_biryp_fortran, zbesh_fortran_first, zbesh_fortran_second, zeros_are_not_equivalent,
    },
};
use num::complex::Complex64;
use rstest::rstest;
use rstest_reuse::apply;

#[apply(bessel_cases)]
fn test_bessel_j(#[case] order: f64, #[case] zr: f64, #[case] zi: f64) {
    let z = Complex64::new(zr, zi);
    let actual = bessel_j(order, z);

    let expected = bessel_j_ref(order, z.into());
    assert_results_are_equal(&actual, &expected, &Vec::new(), 1e6);
}

#[apply(bessel_cases)]
#[trace]
fn test_bessel_i(#[case] order: f64, #[case] zr: f64, #[case] zi: f64) {
    let z = Complex64::new(zr, zi);
    let actual = bessel_i(order, z);
    // dbg!(&actual);

    let expected = bessel_i_ref(order, z.into());
    assert_results_are_equal(&actual, &expected, &Vec::new(), 1e6);
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

#[rstest]
fn test_bessel_extremes(
    #[values(
        (complex_bessel_j as BesselSig, zbesj_fortran as BesselFortranSig),
        (complex_bessel_i as BesselSig, zbesi_fortran as BesselFortranSig),
        (complex_hankel1 as BesselSig , zbesh_fortran_first as BesselFortranSig),
        (complex_hankel2 as BesselSig , zbesh_fortran_second as BesselFortranSig),
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
    check_against_fortran(order, z, scaling, n, rust_fn, fortran_fn, 1e6);
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

// #[rstest]
// fn test_bessel_extremes_f32(
//     #[values(
//         (complex_bessel_j as BesselSig<f32>, zbesj_fortran as BesselFortranSig),
//         (complex_bessel_i as BesselSig<f32>, zbesi_fortran as BesselFortranSig),
//         (complex_hankel1 as BesselSig<f32> , zbesh_fortran_first as BesselFortranSig),
//         (complex_hankel2 as BesselSig<f32> , zbesh_fortran_second as BesselFortranSig),
//         (complex_bessel_k as BesselSig<f32>, zbesk_fortran as BesselFortranSig),
//         (complex_bessel_y as BesselSig<f32>, zbesy_fortran as BesselFortranSig),
//         (sig_airy as BesselSig<f32>, sig_airy_fortran as BesselFortranSig),
//         (sig_airyp as BesselSig<f32>, sig_airyp_fortran as BesselFortranSig),
//         (sig_biry as BesselSig<f32>, sig_biry_fortran as BesselFortranSig),
//         (sig_biryp as BesselSig<f32>, sig_biryp_fortran as BesselFortranSig),
//     )]
//     (rust_fn, fortran_fn): (BesselSig<f32>, BesselFortranSig),

//     #[values(0.0, f32::EPSILON, f32::MIN_POSITIVE, 1.0, f32::MAX)] order: f32,
//     #[values(
//         f32::MIN,
//         -f32::EPSILON,
//         // -f32::MIN_POSITIVE,
//         0.0,
//         // f32::MIN_POSITIVE, MIN_POSITIVE causes an underflow, which is not caused by fortran called with an f64.
//         f32::EPSILON,
//         1.0,
//         f32::MAX
//     )]
//     zr: f32,
//     #[values(
//         f32::MIN,
//         -f32::EPSILON,
//         -f32::MIN_POSITIVE,
//         0.0,
//         f32::MIN_POSITIVE,
//         f32::EPSILON,
//         1.0,
//         f32::MAX
//     )]
//     zi: f32,
//     #[values(Scaling::Unscaled, Scaling::Scaled)] scaling: Scaling,
// ) {
//     // if zr.abs() > 25.0 || zi.abs() > 25.0 || order > 25.0 {
//     //     return;
//     // }
//     // if zr != 0.0 && (zr as f32) == 0.0 {
//     //     return;
//     // }
//     // if zi != 0.0 && (zi as f32) == 0.0 {
//     //     return;
//     // }
//     // if order != 0.0 && (order as f32) == 0.0 {
//     //     return;
//     // }

//     let n = 9;
//     let z = Complex::new(zr, zi);

//     if z.abs() < f32::MACHINE_CONSTANTS.underflow_limit {
//         // MIN_POSITIVE causes an underflow, which is not caused by fortran called with an f64.
//         return;
//     }

//     check_against_fortran(order, z, scaling, n, rust_fn, fortran_fn, 1e4);
//     // let n = 9;
//     // let z = Complex::<f32>::new(zr, zi);
//     // let f32_result = rust_fn(z, order, scaling, n);
//     // let z64 = Complex::<f64>::new(zr as f64, zi as f64);
//     // let fortran_vals = fortran_fn(order as f64, z64, scaling as i32, n);

//     // let fortran_vals_f32 = (
//     //     fortran_vals
//     //         .0
//     //         .iter()
//     //         .map(|x| Complex::<f64>::new(x.re as f64, x.im as f64))
//     //         .collect(),
//     //     fortran_vals.1,
//     //     fortran_vals.2,
//     // );
//     // assert_results_are_equal(&f32_result, &fortran_vals_f32, &vec![], 1e4_f32);
// }
