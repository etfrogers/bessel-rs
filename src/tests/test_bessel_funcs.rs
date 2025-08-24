use super::{
    BesselFortranSig, BesselSig, bessel_cases, bessel_h_ref, check_against_fortran,
    check_complex_arrays_equal, zbesh_first, zbesh_fortran_first, zbesh_fortran_second,
    zbesh_second, zbesi_fortran, zbesj_fortran,
};
use complex_bessel_rs::bessel_i::bessel_i as bessel_i_ref;
use complex_bessel_rs::bessel_j::bessel_j as bessel_j_ref;

use crate::{
    BesselError, HankelKind, Scaling,
    amos::{zbesi, zbesj},
    bessel_i, bessel_j, hankel,
};
use num::complex::Complex64;
use rstest::rstest;
use rstest_reuse::apply;

#[apply(bessel_cases)]
fn test_bessel_j(#[case] order: f64, #[case] zr: f64, #[case] zi: f64) {
    let z = Complex64::new(zr, zi);
    let actual = bessel_j(order, z);

    if actual == Err(BesselError::NotYetImplemented) {
        todo!()
    }

    let expected = bessel_j_ref(order, z.into());
    if let Ok(actual) = actual {
        check_complex_arrays_equal(&actual, &expected.unwrap(), &Vec::new());
    } else {
        assert_eq!(actual.unwrap_err().error_code(), expected.unwrap_err())
    }
}

#[apply(bessel_cases)]
#[trace]
fn test_bessel_i(#[case] order: f64, #[case] zr: f64, #[case] zi: f64) {
    let z = Complex64::new(zr, zi);
    let actual = bessel_i(order, z);
    dbg!(&actual);
    if actual == Err(BesselError::NotYetImplemented) {
        todo!()
    }

    let expected = bessel_i_ref(order, z.into());
    if let Ok(actual) = actual {
        check_complex_arrays_equal(&actual, &expected.unwrap(), &Vec::new());
    } else {
        assert_eq!(actual.unwrap_err().error_code(), expected.unwrap_err())
    }
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
    if actual == Err(BesselError::NotYetImplemented) {
        return;
        //todo!()
    }

    let expected = bessel_h_ref(order, z.into(), kind);
    // dbg!(&expected);
    if let Ok(actual) = actual {
        check_complex_arrays_equal(&actual, &expected.unwrap(), &Vec::new());
    } else {
        assert_eq!(actual.unwrap_err().error_code(), expected.unwrap_err())
    }
}

#[rstest]
fn test_bessel_extremes(
    #[values(
        (zbesj as BesselSig, zbesj_fortran as BesselFortranSig),
        (zbesi as BesselSig, zbesi_fortran as BesselFortranSig),
        (zbesh_first as BesselSig , zbesh_fortran_first as BesselFortranSig),
        (zbesh_second as BesselSig , zbesh_fortran_second as BesselFortranSig),
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
