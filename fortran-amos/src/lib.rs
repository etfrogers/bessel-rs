use num::complex::Complex64;
use std::ffi::{c_double, c_int};

#[link(name = "amos_testing", kind = "static")]
#[link(name = "gfortran")]
unsafe extern "C" {
    pub fn zbesy_wrap_testing(
        zr: c_double,
        zi: c_double,
        nu: c_double,
        kode: c_int,
        N: c_int,
        cyr: *mut c_double,
        cyi: *mut c_double,
        nz: *mut c_int,
        cwrkr: *mut c_double,
        cwrki: *mut c_double,
        ierr: *mut c_int,
    );

    pub fn zbesj_wrap_testing(
        zr: c_double,
        zi: c_double,
        nu: c_double,
        kode: c_int,
        n: c_int,
        cyr: *mut c_double,
        cyi: *mut c_double,
        nz: *mut c_int,
        ierr: *mut c_int,
    );

    pub fn zbesk_wrap_testing(
        zr: c_double,
        zi: c_double,
        nu: c_double,
        kode: c_int,
        n: c_int,
        cyr: *mut c_double,
        cyi: *mut c_double,
        nz: *mut c_int,
        ierr: *mut c_int,
    );

    pub fn zbesi_wrap_testing(
        zr: c_double,
        zi: c_double,
        nu: c_double,
        kode: c_int,
        n: c_int,
        cyr: *mut c_double,
        cyi: *mut c_double,
        nz: *mut c_int,
        ierr: *mut c_int,
    );

    pub fn zbesh_wrap_testing(
        zr: c_double,
        zi: c_double,
        nu: c_double,
        kode: c_int,
        M: c_int,
        N: c_int,
        cyr: *mut c_double,
        cyi: *mut c_double,
        nz: *mut c_int,
        ierr: *mut c_int,
    );

    pub unsafe fn zairy_wrap_testing(
        zr: c_double,
        zi: c_double,
        id: c_int,
        kode: c_int,
        air: *mut c_double,
        aii: *mut c_double,
        nz: *mut c_int,
        ierr: *mut c_int,
    );

    pub unsafe fn zbiry_wrap_testing(
        zr: c_double,
        zi: c_double,
        id: c_int,
        kode: c_int,
        bir: *mut c_double,
        bii: *mut c_double,
        ierr: *mut c_int,
    );
}

pub fn zbesj_fortran(
    order: f64,
    z: Complex64,
    scaling: i32,
    n: usize,
) -> (Vec<Complex64>, usize, i32) {
    let mut cyr = vec![0.0; n];
    let mut cyi = vec![0.0; n];
    let mut nz = 0;
    let mut ierr = 0;

    unsafe {
        zbesj_wrap_testing(
            z.re,
            z.im,
            order,
            scaling,
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
    scaling: i32,
    n: usize,
) -> (Vec<Complex64>, usize, i32) {
    let mut cyr = vec![0.0; n];
    let mut cyi = vec![0.0; n];
    let mut nz = 0;
    let mut ierr = 0;

    unsafe {
        zbesi_wrap_testing(
            z.re,
            z.im,
            order,
            scaling,
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
    scaling: i32,
    n: usize,
) -> (Vec<Complex64>, usize, i32) {
    let mut cyr = vec![0.0; n];
    let mut cyi = vec![0.0; n];
    let mut nz = 0;
    let mut ierr = 0;

    unsafe {
        zbesk_wrap_testing(
            z.re,
            z.im,
            order,
            scaling,
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

pub fn zbesy_fortran(
    order: f64,
    z: Complex64,
    scaling: i32,
    n: usize,
) -> (Vec<Complex64>, usize, i32) {
    let mut cyr = vec![0.0; n];
    let mut cyi = vec![0.0; n];
    let mut cwrkr = vec![0.0; n];
    let mut cwrki = vec![0.0; n];
    let mut nz = 0;
    let mut ierr = 0;

    unsafe {
        zbesy_wrap_testing(
            z.re,
            z.im,
            order,
            scaling,
            n.try_into().unwrap(),
            cyr.as_mut_ptr(),
            cyi.as_mut_ptr(),
            &mut nz,
            cwrkr.as_mut_ptr(),
            cwrki.as_mut_ptr(),
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
    scaling: i32,
    kind: i32,
    n: usize,
) -> (Vec<Complex64>, usize, i32) {
    let mut cyr = vec![0.0; n];
    let mut cyi = vec![0.0; n];
    let mut nz = 0;
    let mut ierr = 0;

    unsafe {
        zbesh_wrap_testing(
            z.re,
            z.im,
            order,
            scaling,
            kind,
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

pub fn zairy_fortran(z: Complex64, is_derivative: bool, scaling: i32) -> (Complex64, usize, i32) {
    let mut yr = 0.0;
    let mut yi = 0.0;
    let mut nz = 0;
    let mut ierr = 0;

    unsafe {
        zairy_wrap_testing(
            z.re,
            z.im,
            is_derivative as i32,
            scaling,
            &mut yr,
            &mut yi,
            &mut nz,
            &mut ierr,
        );
    }
    (Complex64::new(yr, yi), nz.try_into().unwrap(), ierr)
}

pub fn zbiry_fortran(z: Complex64, is_derivative: bool, scaling: i32) -> (Complex64, usize, i32) {
    let mut yr = 0.0;
    let mut yi = 0.0;
    let mut ierr = 0;

    unsafe {
        zbiry_wrap_testing(
            z.re,
            z.im,
            is_derivative as i32,
            scaling,
            &mut yr,
            &mut yi,
            &mut ierr,
        );
    }
    (Complex64::new(yr, yi), 0, ierr)
}
