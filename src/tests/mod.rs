use f64;

use num::complex::Complex64;

use crate::Scaling;
use crate::amos::{BesselResult, HankelKind};

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
mod test_large_n;
mod test_machine_consts;
mod utils;

const TOLERANCE_MARGIN: f64 = 1e8;

type BesselSig = fn(Complex64, f64, Scaling, usize) -> BesselResult;
type BesselFortranSig = fn(f64, Complex64, Scaling, usize) -> (Vec<Complex64>, usize, i32);

// This function needed as complex-bessel-rs (which is used for the other *_ref functions) does not
// provide a bessel_h function.
fn bessel_h_ref(order: f64, z: Complex64, kind: HankelKind) -> Result<Complex64, i32> {
    let (y, _, ierr) = zbesh_fortran(order, z, Scaling::Unscaled, kind, 1);
    if ierr != 0 { Err(ierr) } else { Ok(y[0]) }
}

// TODO test bad inputs
