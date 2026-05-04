use crate::HankelKind;
use num::complex::Complex64;

use crate::tests::zbesh_fortran;

pub fn zbesh_fortran_first(
    order: f64,
    z: Complex64,
    scaling: i32,
    n: usize,
) -> (Vec<Complex64>, usize, i32) {
    zbesh_fortran(order, z, scaling, HankelKind::First as i32, n)
}

pub fn zbesh_fortran_second(
    order: f64,
    z: Complex64,
    scaling: i32,
    n: usize,
) -> (Vec<Complex64>, usize, i32) {
    zbesh_fortran(order, z, scaling, HankelKind::Second as i32, n)
}
