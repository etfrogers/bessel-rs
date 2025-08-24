use f64;
use rstest_reuse::apply;

use num::complex::Complex64;
use rstest::rstest;

use crate::amos::bindings::{zbesh_wrap, zbesi_wrap, zbesj_wrap};
use crate::amos::{BesselResult, HankelKind, zbesh, zbesi, zbesj};
use crate::{BesselError, Scaling, bessel_i, bessel_j, hankel};

use complex_bessel_rs::bessel_i::bessel_i as bessel_i_ref;
use complex_bessel_rs::bessel_j::bessel_j as bessel_j_ref;
use parametrisation::bessel_j_cases;
pub use utils::{check_against_fortran, check_complex_arrays_equal};

#[cfg(feature = "random_tests")]
pub use utils::fortran_bess_loop;

mod parametrisation;
#[cfg(feature = "random_tests")]
mod random_tests;
mod test_gamma_ln;
mod test_machine_consts;
mod utils;

const TOLERANCE_MARGIN: f64 = 1e8;

type BesselSig = fn(Complex64, f64, Scaling, usize) -> BesselResult;
type BesselFortranSig = fn(f64, Complex64, Scaling, usize) -> (Vec<Complex64>, usize, i32);

#[apply(bessel_j_cases)]
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

#[apply(bessel_j_cases)]
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

#[apply(bessel_j_cases)]
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

#[apply(bessel_j_cases)]
#[trace]
fn test_bessel_large_n_real(
    #[case] order: f64,
    #[case] zr: f64,
    #[case] _zi: f64,
    #[values(3, 4, 9, 100)] n: usize,
    #[values(Scaling::Unscaled, Scaling::Scaled)] scaling: Scaling,
    #[values(
        (zbesj as BesselSig, zbesj_fortran as BesselFortranSig),
        (zbesi as BesselSig, zbesi_fortran as BesselFortranSig),
        (zbesh_first as BesselSig , zbesh_fortran_first as BesselFortranSig),
        (zbesh_second as BesselSig , zbesh_fortran_second as BesselFortranSig),
    )]
    (rust_fn, fortran_fn): (BesselSig, BesselFortranSig),
    // #[values(4)] n: usize,
    // #[values(Scaling::Unscaled)] scaling: Scaling,
) {
    // ignores the zi input
    check_against_fortran(order, zr.into(), scaling, n, rust_fn, fortran_fn);
}

#[apply(bessel_j_cases)]
#[trace]
fn test_bessel_large_n_complex(
    #[case] order: f64,
    #[case] zr: f64,
    #[case] zi: f64,
    #[values(3, 4, 9, 100)] n: usize,
    #[values(Scaling::Unscaled, Scaling::Scaled)] scaling: Scaling,
    #[values(
        (zbesj as BesselSig, zbesj_fortran as BesselFortranSig),
        (zbesi as BesselSig, zbesi_fortran as BesselFortranSig),
        (zbesh_first as BesselSig , zbesh_fortran_first as BesselFortranSig),
        (zbesh_second as BesselSig , zbesh_fortran_second as BesselFortranSig),
    )]
    (rust_fn, fortran_fn): (BesselSig, BesselFortranSig),
    // #[values(9)] n: usize,
    // #[values(Scaling::Unscaled)] scaling: Scaling,
) {
    let z = Complex64::new(zr, zi);
    check_against_fortran(order, z, scaling, n, rust_fn, fortran_fn);
}

fn zbesh_first(z: Complex64, order: f64, scaling: Scaling, n: usize) -> BesselResult {
    zbesh(z, order, scaling, HankelKind::First, n)
}

fn zbesh_fortran_first(
    order: f64,
    z: Complex64,
    scaling: Scaling,
    n: usize,
) -> (Vec<Complex64>, usize, i32) {
    zbesh_fortran(order, z, scaling, HankelKind::First, n)
}

fn zbesh_second(z: Complex64, order: f64, scaling: Scaling, n: usize) -> BesselResult {
    zbesh(z, order, scaling, HankelKind::Second, n)
}

fn zbesh_fortran_second(
    order: f64,
    z: Complex64,
    scaling: Scaling,
    n: usize,
) -> (Vec<Complex64>, usize, i32) {
    zbesh_fortran(order, z, scaling, HankelKind::Second, n)
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

fn zbesj_fortran(
    order: f64,
    z: Complex64,
    scaling: Scaling,
    n: usize,
) -> (Vec<Complex64>, usize, i32) {
    let mut cyr: Vec<f64> = Vec::with_capacity(n);
    let mut cyi: Vec<f64> = Vec::with_capacity(n);
    let mut nz = 0;
    let mut ierr = 0;

    let r_uninit = cyr.spare_capacity_mut();
    let i_uninit = cyi.spare_capacity_mut();
    unsafe {
        zbesj_wrap(
            z.re,
            z.im,
            order,
            scaling as i32,
            n.try_into().unwrap(),
            r_uninit.as_mut_ptr().cast(),
            i_uninit.as_mut_ptr().cast(),
            &mut nz,
            &mut ierr,
        );
        cyr.set_len(n);
        cyi.set_len(n);
    }
    let cy = cyr
        .into_iter()
        .zip(cyi)
        .map(|(r, i)| Complex64::new(r, i))
        .collect();
    (cy, nz.try_into().unwrap(), ierr)
}

fn zbesi_fortran(
    order: f64,
    z: Complex64,
    scaling: Scaling,
    n: usize,
) -> (Vec<Complex64>, usize, i32) {
    let mut cyr: Vec<f64> = Vec::with_capacity(n);
    let mut cyi: Vec<f64> = Vec::with_capacity(n);
    let mut nz = 0;
    let mut ierr = 0;

    let r_uninit = cyr.spare_capacity_mut();
    let i_uninit = cyi.spare_capacity_mut();
    unsafe {
        zbesi_wrap(
            z.re,
            z.im,
            order,
            scaling as i32,
            n.try_into().unwrap(),
            r_uninit.as_mut_ptr().cast(),
            i_uninit.as_mut_ptr().cast(),
            &mut nz,
            &mut ierr,
        );
        cyr.set_len(n);
        cyi.set_len(n);
    }
    let cy = cyr
        .into_iter()
        .zip(cyi)
        .map(|(r, i)| Complex64::new(r, i))
        .collect();
    (cy, nz.try_into().unwrap(), ierr)
}

fn bessel_h_ref(order: f64, z: Complex64, kind: HankelKind) -> Result<Complex64, i32> {
    let (y, _, ierr) = zbesh_fortran(order, z, Scaling::Unscaled, kind, 1);
    if ierr != 0 { Err(ierr) } else { Ok(y[0]) }
}

fn zbesh_fortran(
    order: f64,
    z: Complex64,
    scaling: Scaling,
    kind: HankelKind,
    n: usize,
) -> (Vec<Complex64>, usize, i32) {
    let mut cyr: Vec<f64> = Vec::with_capacity(n);
    let mut cyi: Vec<f64> = Vec::with_capacity(n);
    let mut nz = 0;
    let mut ierr = 0;

    let r_uninit = cyr.spare_capacity_mut();
    let i_uninit = cyi.spare_capacity_mut();
    unsafe {
        zbesh_wrap(
            z.re,
            z.im,
            order,
            scaling as i32,
            kind.into(),
            n.try_into().unwrap(),
            r_uninit.as_mut_ptr().cast(),
            i_uninit.as_mut_ptr().cast(),
            &mut nz,
            &mut ierr,
        );
        cyr.set_len(n);
        cyi.set_len(n);
    }
    let cy = cyr
        .into_iter()
        .zip(cyi)
        .map(|(r, i)| Complex64::new(r, i))
        .collect();
    (cy, nz.try_into().unwrap(), ierr)
}
