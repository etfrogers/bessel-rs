use f64;
use rstest_reuse::{apply, template};

use num::complex::Complex64;
use rstest::rstest;

use crate::amos::bindings::{zbesh_wrap, zbesi_wrap, zbesj_wrap};
use crate::amos::{BesselResult, HankelKind, zbesh, zbesi, zbesj};
use crate::{BesselError, Scaling, bessel_i, bessel_j, hankel};

use complex_bessel_rs::bessel_i::bessel_i as bessel_i_ref;
use complex_bessel_rs::bessel_j::bessel_j as bessel_j_ref;
pub use utils::{check_against_fortran, check_complex_arrays_equal};

#[cfg(feature = "random_tests")]
pub use utils::fortran_bess_loop;

#[cfg(feature = "random_tests")]
mod random_tests;
mod test_gamma_ln;
mod test_machine_consts;
mod utils;

const TOLERANCE_MARGIN: f64 = 10_000.0;

type BesselSig = fn(Complex64, f64, Scaling, usize) -> BesselResult;
type BesselFortranSig = fn(f64, Complex64, Scaling, usize) -> (Vec<Complex64>, usize, i32);

#[template]
#[rstest]
#[case(4.0, 2.1, 0.0)]
#[case(5.0, 2.0001, 0.0)]
#[case(340.0, 35.0001, 0.0)]
#[case(407.3,-325.1, 635.2)]
#[case(465.0,-867.0, -448.0)] // 5
#[case(10.711871220659752, -6.89931119845653, -9.408182256887017)]
#[case(8.544950848779838, -8.645033163963603, 18.976439189605003,)]
#[case(21.04, 53.19, -40.77)]
#[case(4.0, 2.1, 0.0)]
#[case(5.0, 2.0001, 0.0)] // 10
#[case(340.0, 35.0001, 0.0)]
#[case(899.6,-35.7,317.8)]
#[case(531.0,-106.7,-16.0)]
#[case(531.0,-106.0,-16.0)]
#[case(433.0,16.874,-38.17)] //15
#[case(433.0,16.8,-38.17)]
#[case(311.2078694557452,-10.990141170168044,-25.70154097357056,)]
#[case(8.544950848779838, -8.645033163963603, 18.976439189605003,)]
#[case(17.5, 70.3, 37.4)]
#[case(13.337522865795481, -29.8266399174247, 17.66323218839807)] //20
#[case(5423.24, -7915.11, -3113.95)]
#[case(2213.0, -1813.0, -1033.0)]
#[case(5514.86274463943, -9489.650336481069, 4951.6909981261)]
#[case(2.74e-288, 6.33e-166, 7.53e-275)]
#[case(1.51e-150, -3.07e-118, 3.51e-42)] //25
#[case(2.637e-27, -4.01e-50, 0.0)]
#[case(4.0e-132, 0e0, 445.0)]
#[case(8714.0, 8904.0, -10.0)]
#[case(60.9, 246.2, -982.5)]
#[case(40.5, 1673.3, -4.0)] // 30
#[case(2634.5, -2634.5, 14.1)]
#[case(5.007e-14, 4.401331657952316e-5, -3.6e-6)]
#[case(1719.3, 920.1, 0.0)]
#[case(3.5695132850479827e3, -2.2313404290100934e3, 8.646324128723001e3)]
#[case(0.28008208034835413, -2435.84398720043, -9106.813568430613)] // 35
#[case(35.42423142304685, 2689.1019240048972, -688.7899868054337)]
#[case(1.0111752223029848, 7037.518427975952, -685.0803465010631)]
#[case(9491.159287083694, -2404.8869667701747, -6391.664651975572)]
#[case(3.468367867017804e0, -1.8067397106295227e-254, -3.0255676077184667e-21)]
#[case(6.946702885186345e-149, 0e0, -6.691424259254966e2)] // 40
#[case(3.684122892548987e3, -5.107972475729046e3, 5.916387337090975e3)]
#[case(7.107636998006379e3, -1.867258055869096e3, 4.865284129480511e3)]
fn bessel_j_cases(#[case] order: f64, #[case] zr: f64, #[case] zi: f64) {}

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
