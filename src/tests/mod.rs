use f64;
pub use rstest_reuse::apply;

use num::complex::Complex64;
use rstest::rstest;

use crate::Scaling;
use crate::amos::{BesselResult, HankelKind, zbesi, zbesj};

pub use bessel_h_wrappers::*;
pub use fortran_calls::*;
use parametrisation::bessel_cases;

pub use utils::{check_against_fortran, check_complex_arrays_equal};

#[cfg(feature = "random_tests")]
pub use utils::fortran_bess_loop;

mod bessel_h_wrappers;
mod fortran_calls;
mod parametrisation;
#[cfg(feature = "random_tests")]
mod random_tests;
mod test_bessel_funcs;
mod test_gamma_ln;
mod test_machine_consts;
mod utils;

const TOLERANCE_MARGIN: f64 = 1e8;

type BesselSig = fn(Complex64, f64, Scaling, usize) -> BesselResult;
type BesselFortranSig = fn(f64, Complex64, Scaling, usize) -> (Vec<Complex64>, usize, i32);

#[apply(bessel_cases)]
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

#[apply(bessel_cases)]
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

// This function needed as complex-bessel-rs (which is used for the other *_ref functions) does not
// provide a bessel_h function.
fn bessel_h_ref(order: f64, z: Complex64, kind: HankelKind) -> Result<Complex64, i32> {
    let (y, _, ierr) = zbesh_fortran(order, z, Scaling::Unscaled, kind, 1);
    if ierr != 0 { Err(ierr) } else { Ok(y[0]) }
}
