use rstest::rstest;
use rstest_reuse::apply;

use super::{
    BesselFortranSig, BesselSig, bessel_cases, check_against_fortran, zbesh_first,
    zbesh_fortran_first, zbesh_fortran_second, zbesh_second, zbesi_fortran, zbesj_fortran,
    zbesk_fortran, zbesy_fortran,
};
use crate::{
    Scaling,
    amos::{complex_bessel_k, complex_bessel_y, complex_bessel_j, complex_bessel_i},
};
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
        (zbesh_first as BesselSig , zbesh_fortran_first as BesselFortranSig),
        (zbesh_second as BesselSig , zbesh_fortran_second as BesselFortranSig),
        (complex_bessel_k as BesselSig, zbesk_fortran as BesselFortranSig),
        (complex_bessel_y as BesselSig, zbesy_fortran as BesselFortranSig)
    )]
    (rust_fn, fortran_fn): (BesselSig, BesselFortranSig),
    // #[values(4)] n: usize,
    // #[values(Scaling::Unscaled)] scaling: Scaling,
) {
    // ignores the zi input
    check_against_fortran(order, zr.into(), scaling, n, rust_fn, fortran_fn);
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
        (zbesh_first as BesselSig , zbesh_fortran_first as BesselFortranSig),
        (zbesh_second as BesselSig , zbesh_fortran_second as BesselFortranSig),
        (complex_bessel_k as BesselSig, zbesk_fortran as BesselFortranSig),
        (complex_bessel_y as BesselSig, zbesy_fortran as BesselFortranSig)
    )]
    (rust_fn, fortran_fn): (BesselSig, BesselFortranSig),
    // #[values(9)] n: usize,
    // #[values(Scaling::Unscaled)] scaling: Scaling,
) {
    let z = Complex64::new(zr, zi);
    check_against_fortran(order, z, scaling, n, rust_fn, fortran_fn);
}
