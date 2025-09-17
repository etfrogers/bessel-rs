use crate::amos::{
    HankelKind,
    bindings::{zbesh_wrap, zbesi_wrap, zbesj_wrap, zbesk_wrap},
};
use num::complex::Complex64;

use crate::Scaling;

pub fn zbesj_fortran(
    order: f64,
    z: Complex64,
    scaling: Scaling,
    n: usize,
) -> (Vec<Complex64>, usize, i32) {
    let mut cyr = vec![0.0; n];
    let mut cyi = vec![0.0; n];
    let mut nz = 0;
    let mut ierr = 0;

    unsafe {
        zbesj_wrap(
            z.re,
            z.im,
            order,
            scaling as i32,
            n.try_into().unwrap(),
            cyr.as_mut_ptr(),
            cyi.as_mut_ptr(),
            &mut nz,
            &mut ierr,
        );
    }
    let cy = cyr
        .into_iter()
        .zip(cyi)
        .map(|(r, i)| Complex64::new(r, i))
        .collect();
    (cy, nz.try_into().unwrap(), ierr)
}

pub fn zbesi_fortran(
    order: f64,
    z: Complex64,
    scaling: Scaling,
    n: usize,
) -> (Vec<Complex64>, usize, i32) {
    let mut cyr = vec![0.0; n];
    let mut cyi = vec![0.0; n];
    let mut nz = 0;
    let mut ierr = 0;

    unsafe {
        zbesi_wrap(
            z.re,
            z.im,
            order,
            scaling as i32,
            n.try_into().unwrap(),
            cyr.as_mut_ptr(),
            cyi.as_mut_ptr(),
            &mut nz,
            &mut ierr,
        );
    }
    let cy = cyr
        .into_iter()
        .zip(cyi)
        .map(|(r, i)| Complex64::new(r, i))
        .collect();
    (cy, nz.try_into().unwrap(), ierr)
}

pub fn zbesk_fortran(
    order: f64,
    z: Complex64,
    scaling: Scaling,
    n: usize,
) -> (Vec<Complex64>, usize, i32) {
    let mut cyr = vec![0.0; n];
    let mut cyi = vec![0.0; n];
    let mut nz = 0;
    let mut ierr = 0;

    unsafe {
        zbesk_wrap(
            z.re,
            z.im,
            order,
            scaling as i32,
            n.try_into().unwrap(),
            cyr.as_mut_ptr(),
            cyi.as_mut_ptr(),
            &mut nz,
            &mut ierr,
        );
    }
    let cy = cyr
        .into_iter()
        .zip(cyi)
        .map(|(r, i)| Complex64::new(r, i))
        .collect();
    (cy, nz.try_into().unwrap(), ierr)
}

pub fn zbesh_fortran(
    order: f64,
    z: Complex64,
    scaling: Scaling,
    kind: HankelKind,
    n: usize,
) -> (Vec<Complex64>, usize, i32) {
    let mut cyr = vec![0.0; n];
    let mut cyi = vec![0.0; n];
    let mut nz = 0;
    let mut ierr = 0;

    unsafe {
        zbesh_wrap(
            z.re,
            z.im,
            order,
            scaling as i32,
            kind as i32,
            n.try_into().unwrap(),
            cyr.as_mut_ptr(),
            cyi.as_mut_ptr(),
            &mut nz,
            &mut ierr,
        );
    }
    let cy = cyr
        .into_iter()
        .zip(cyi)
        .map(|(r, i)| Complex64::new(r, i))
        .collect();
    (cy, nz.try_into().unwrap(), ierr)
}
