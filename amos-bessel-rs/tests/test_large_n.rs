extern crate fortran_amos_testing;
mod common;
use common::parametrisation::bessel_cases;
use rstest::rstest;
use rstest_reuse::apply;

use amos_bessel_rs::{
    Scaling, complex_bessel_i, complex_bessel_j, complex_bessel_k, complex_bessel_y,
    complex_hankel1, complex_hankel2,
    test_utils::{
        BesselFortranSig, BesselSig, check_against_fortran, zbesh_fortran_first,
        zbesh_fortran_second,
    },
};
use fortran_amos_testing::{zbesi_fortran, zbesj_fortran, zbesk_fortran, zbesy_fortran};
use num::complex::Complex64;

#[apply(bessel_cases)]
#[trace]
fn test_bessel_large_n_real(
    #[case] order: f64,
    #[case] zr: f64,
    #[case] _zi: f64,
    #[values(3, 4, 9, 100)] n: usize,
    #[values(Scaling::Unscaled, Scaling::Scaled)] scaling: Scaling,
    #[values(
        (complex_bessel_j as BesselSig, zbesj_fortran as BesselFortranSig),
        (complex_bessel_i as BesselSig, zbesi_fortran as BesselFortranSig),
        (complex_hankel1 as BesselSig , zbesh_fortran_first as BesselFortranSig),
        (complex_hankel2 as BesselSig , zbesh_fortran_second as BesselFortranSig),
        (complex_bessel_k as BesselSig, zbesk_fortran as BesselFortranSig),
        (complex_bessel_y as BesselSig, zbesy_fortran as BesselFortranSig)
    )]
    (rust_fn, fortran_fn): (BesselSig, BesselFortranSig),
) {
    // ignores the zi input
    // let z = Complex64::new(zr, zi);
    // check_against_fortran(order, z, scaling, n, rust_fn, fortran_fn);
    check_against_fortran(order, zr.into(), scaling, n, rust_fn, fortran_fn, 1e6);
}

#[apply(bessel_cases)]
#[trace]
fn test_bessel_large_n_complex(
    #[case] order: f64,
    #[case] zr: f64,
    #[case] zi: f64,
    #[values(3, 4, 9, 100)] n: usize,
    #[values(Scaling::Unscaled, Scaling::Scaled)] scaling: Scaling,
    #[values(
        (complex_bessel_j as BesselSig, zbesj_fortran as BesselFortranSig),
        (complex_bessel_i as BesselSig, zbesi_fortran as BesselFortranSig),
        (complex_hankel1 as BesselSig , zbesh_fortran_first as BesselFortranSig),
        (complex_hankel2 as BesselSig , zbesh_fortran_second as BesselFortranSig),
        (complex_bessel_k as BesselSig, zbesk_fortran as BesselFortranSig),
        (complex_bessel_y as BesselSig, zbesy_fortran as BesselFortranSig)
    )]
    (rust_fn, fortran_fn): (BesselSig, BesselFortranSig),
) {
    let z = Complex64::new(zr, zi);
    check_against_fortran(order, z, scaling, n, rust_fn, fortran_fn, 1e6);
}
