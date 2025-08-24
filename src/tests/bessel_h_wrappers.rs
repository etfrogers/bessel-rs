use num::complex::Complex64;

use crate::{
    Scaling,
    amos::{BesselResult, HankelKind, zbesh},
    tests::zbesh_fortran,
};

pub fn zbesh_first(z: Complex64, order: f64, scaling: Scaling, n: usize) -> BesselResult {
    zbesh(z, order, scaling, HankelKind::First, n)
}

pub fn zbesh_fortran_first(
    order: f64,
    z: Complex64,
    scaling: Scaling,
    n: usize,
) -> (Vec<Complex64>, usize, i32) {
    zbesh_fortran(order, z, scaling, HankelKind::First, n)
}

pub fn zbesh_second(z: Complex64, order: f64, scaling: Scaling, n: usize) -> BesselResult {
    zbesh(z, order, scaling, HankelKind::Second, n)
}

pub fn zbesh_fortran_second(
    order: f64,
    z: Complex64,
    scaling: Scaling,
    n: usize,
) -> (Vec<Complex64>, usize, i32) {
    zbesh_fortran(order, z, scaling, HankelKind::Second, n)
}
