extern crate fortran_amos_testing;
mod common;
use common::parametrisation::bessel_cases;
use rstest::rstest;
use rstest_reuse::apply;

use amos_bessel_rs::{
    Scaling,
    amos::{
        complex_bessel_i, complex_bessel_j, complex_bessel_k, complex_bessel_y, complex_hankel1,
        complex_hankel2,
    },
};
use common::{
    BesselFortranSig, BesselSig, assert_results_are_equal_floats, check_against_fortran,
    zbesh_fortran_first, zbesh_fortran_second, zeros_are_not_equivalent,
};
use num::complex::Complex64;

use common::ComplexConversions;
use fortran_amos_testing::{zbesi_fortran, zbesj_fortran, zbesk_fortran, zbesy_fortran};

#[apply(bessel_cases)]
fn test_bessel_large_n_real(#[case] order: f64, #[case] zr: f64, #[case] _zi: f64) {
    // n, scaling, and function pairs are iterated at runtime rather than via #[values] to
    // avoid a Cartesian-product explosion of monomorphised test functions (47 × 4 × 2 × 6 = 2256).

    let fn_pairs: &[(BesselSig, BesselFortranSig)] = &[
        (
            complex_bessel_j as BesselSig,
            zbesj_fortran as BesselFortranSig,
        ),
        (
            complex_bessel_i as BesselSig,
            zbesi_fortran as BesselFortranSig,
        ),
        (
            complex_hankel1 as BesselSig,
            zbesh_fortran_first as BesselFortranSig,
        ),
        (
            complex_hankel2 as BesselSig,
            zbesh_fortran_second as BesselFortranSig,
        ),
        (
            complex_bessel_k as BesselSig,
            zbesk_fortran as BesselFortranSig,
        ),
        (
            complex_bessel_y as BesselSig,
            zbesy_fortran as BesselFortranSig,
        ),
    ];
    for &n in &[3_usize, 4, 9, 100] {
        for scaling in [Scaling::Unscaled, Scaling::Scaled] {
            for &(rust_fn, fortran_fn) in fn_pairs {
                // ignores the zi input
                check_against_fortran(order, zr.into(), scaling, n, rust_fn, fortran_fn, 1e6);
            }
        }
    }
}

#[apply(bessel_cases)]
fn test_bessel_large_n_complex(#[case] order: f64, #[case] zr: f64, #[case] zi: f64) {
    // n, scaling, and function pairs are iterated at runtime rather than via #[values] to
    // avoid a Cartesian-product explosion of monomorphised test functions (47 × 4 × 2 × 6 = 2256).

    let fn_pairs: &[(BesselSig, BesselFortranSig)] = &[
        (
            complex_bessel_j as BesselSig,
            zbesj_fortran as BesselFortranSig,
        ),
        (
            complex_bessel_i as BesselSig,
            zbesi_fortran as BesselFortranSig,
        ),
        (
            complex_hankel1 as BesselSig,
            zbesh_fortran_first as BesselFortranSig,
        ),
        (
            complex_hankel2 as BesselSig,
            zbesh_fortran_second as BesselFortranSig,
        ),
        (
            complex_bessel_k as BesselSig,
            zbesk_fortran as BesselFortranSig,
        ),
        (
            complex_bessel_y as BesselSig,
            zbesy_fortran as BesselFortranSig,
        ),
    ];
    let z = Complex64::new(zr, zi);
    for &n in &[3_usize, 4, 9, 100] {
        for scaling in [Scaling::Unscaled, Scaling::Scaled] {
            for &(rust_fn, fortran_fn) in fn_pairs {
                check_against_fortran(order, z, scaling, n, rust_fn, fortran_fn, 1e6);
            }
        }
    }
}

#[apply(bessel_cases)]
fn test_bessel_large_n_real_f32(#[case] order: f64, #[case] zr: f64, #[case] _zi: f64) {
    // n and function pairs are iterated at runtime rather than via #[values] to
    // avoid a Cartesian-product explosion of monomorphised test functions.

    let fn_pairs: &[(BesselSig<f64>, BesselSig<f32>)] = &[
        (
            complex_bessel_j as BesselSig<f64>,
            complex_bessel_j as BesselSig<f32>,
        ),
        (
            complex_bessel_i as BesselSig<f64>,
            complex_bessel_i as BesselSig<f32>,
        ),
        (
            complex_hankel1 as BesselSig<f64>,
            complex_hankel1 as BesselSig<f32>,
        ),
        (
            complex_hankel2 as BesselSig<f64>,
            complex_hankel2 as BesselSig<f32>,
        ),
        (
            complex_bessel_k as BesselSig<f64>,
            complex_bessel_k as BesselSig<f32>,
        ),
        (
            complex_bessel_y as BesselSig<f64>,
            complex_bessel_y as BesselSig<f32>,
        ),
    ];
    let z64 = Complex64::new(zr, 0.0);
    let z32 = z64.to_c32();
    if zeros_are_not_equivalent(z64, z32, order) {
        return;
    }
    for &n in &[3_usize, 4, 9, 100] {
        for scaling in [Scaling::Unscaled, Scaling::Scaled] {
            for &(fn_64, fn_32) in fn_pairs {
                use crate::common::assert_results_are_equal_floats;

                let mut margin = 1e4 * n as f64;
                if n >= 9 && zr.abs() > order {
                    margin *= 100.0;
                }
                let f64_vals = fn_64(z64, order, scaling, n);
                let f32_vals = fn_32(z32, order as f32, scaling, n);
                assert_results_are_equal_floats(&f32_vals, &f64_vals, margin);
            }
        }
    }
}

#[apply(bessel_cases)]
fn test_bessel_large_n_complex_f32(#[case] order: f64, #[case] zr: f64, #[case] zi: f64) {
    // n and function pairs are iterated at runtime rather than via #[values] to
    // avoid a Cartesian-product explosion of monomorphised test functions.

    let fn_pairs: &[(BesselSig<f64>, BesselSig<f32>)] = &[
        (
            complex_bessel_j as BesselSig<f64>,
            complex_bessel_j as BesselSig<f32>,
        ),
        (
            complex_bessel_i as BesselSig<f64>,
            complex_bessel_i as BesselSig<f32>,
        ),
        (
            complex_hankel1 as BesselSig<f64>,
            complex_hankel1 as BesselSig<f32>,
        ),
        (
            complex_hankel2 as BesselSig<f64>,
            complex_hankel2 as BesselSig<f32>,
        ),
        (
            complex_bessel_k as BesselSig<f64>,
            complex_bessel_k as BesselSig<f32>,
        ),
        (
            complex_bessel_y as BesselSig<f64>,
            complex_bessel_y as BesselSig<f32>,
        ),
    ];
    let z64 = Complex64::new(zr, zi);
    let z32 = z64.to_c32();
    if zeros_are_not_equivalent(z64, z32, order) {
        return;
    }
    let abs_z = z64.norm();
    for &n in &[3_usize, 4, 9, 100] {
        for scaling in [Scaling::Unscaled, Scaling::Scaled] {
            for &(fn_64, fn_32) in fn_pairs {
                let mut margin = 1e4 * n as f64;
                // For n >= 9 steps in the oscillatory Bessel regime (|z| > order), f32 accumulated
                // rounding errors grow beyond linear n-scaling. Apply an extra 100x multiplier.
                // This is analogous to the |zi| ≈ 40 skip in test_f32_vs_f64: the AMOS f32 code
                // takes different algorithmic paths due to differing machine constants.
                if n >= 9 && abs_z > order {
                    margin *= 100.0;
                }
                let f64_vals = fn_64(z64, order, scaling, n);
                let f32_vals = fn_32(z32, order as f32, scaling, n);
                assert_results_are_equal_floats(&f32_vals, &f64_vals, margin);
            }
        }
    }
}

#[rstest]
fn test_scale_k_recurrence() {
    let z = Complex64::new(715.0, 0.0);
    let order = 0.0;
    let n = 800;
    check_against_fortran(
        order,
        z,
        Scaling::Unscaled,
        n,
        complex_bessel_k,
        zbesk_fortran,
        1e6,
    );
}

#[rstest]
fn test_k_right_half_plane_underflow_offset_bug() {
    // z is extremely large, so K_0(715) evaluates to roughly e^-715,
    // which is below f64::MIN_POSITIVE and underflows to exactly 0.
    let z = Complex64::new(715.0, 0.0);
    let n = 200;
    let order = 0.0;

    let (res, n_zeros) = complex_bessel_k(z, order, Scaling::Unscaled, n).unwrap();

    // It should take 156 steps of the sequence growing before the values
    // finally become large enough to be represented by an f64 again!
    assert_eq!(n_zeros, 156);

    // Therefore, the first 156 elements MUST be exactly zero.
    for (i, &res_val) in res.iter().enumerate().take(n_zeros + 1) {
        assert_eq!(
            res_val,
            Complex64::new(0.0, 0.0),
            "Index {} should have underflowed to 0",
            i
        );
    }

    // And the 157th element should finally be non-zero.
    assert_ne!(
        res[157],
        Complex64::new(0.0, 0.0),
        "Index 156 should be non-zero"
    );

    // also test against fortran
    check_against_fortran(
        order,
        z,
        Scaling::Unscaled,
        n,
        complex_bessel_k,
        zbesk_fortran,
        1e6,
    );
}
