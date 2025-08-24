use crate::amos::{
    HankelKind,
    bindings::{zbesh_wrap, zbesi_wrap, zbesj_wrap},
};
use num::complex::Complex64;

use crate::Scaling;

pub fn zbesj_fortran(
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

pub fn zbesi_fortran(
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

pub fn zbesh_fortran(
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
