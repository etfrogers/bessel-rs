use crate::tests::utils::check_complex_arrays_equal;

use crate::tests::{
    BesselFortranSig, BesselSig, check_against_fortran, sig_airy, sig_airy_fortran, sig_airyp,
    sig_airyp_fortran, sig_biry, sig_biry_fortran, sig_biryp, sig_biryp_fortran,
    zbesh_fortran_first, zbesh_fortran_second, zbesi_fortran, zbesj_fortran, zbesk_fortran,
    zbesy_fortran,
};
use crate::{
    HankelKind, Scaling, bessel_i, bessel_j, bessel_k, bessel_y, complex_bessel_i,
    complex_bessel_j, complex_bessel_k, complex_bessel_y, complex_hankel1, complex_hankel2, hankel,
    types::BesselError,
};

use complex_bessel::besseli as bessel_i_ref;
use complex_bessel::besselj as bessel_j_ref;
use complex_bessel::besselk as bessel_k_ref;
use complex_bessel::bessely as bessel_y_ref;
use complex_bessel::hankel1 as hankel1_ref;
use complex_bessel::hankel2 as hankel2_ref;

use num::Complex;
#[cfg(test)]
use rstest::rstest;

const ORDERS: [f64; 21] = [
    0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 5.0, 10.0, 25.0, 50.0, 75.0, 85.0, 90.0, 100.0, 150.0, 200.0,
    500.0, 1000.0, -0.5, -1.5, -2.0,
];

const Z_PARTS: [f64; 37] = [
    -50.0, -40.0, -30.0, -25.0, -20.0, -15.0, -12.0, -10.0, -8.0, -6.0, -4.0, -3.0, -2.0, -1.0,
    -0.5, -0.1, -0.001, -1e-6, 0.0, 1e-6, 0.001, 0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0,
    12.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0,
];

// Too many values create macro invocation errors
#[rstest]
fn test_bessel_grid_fortran(
    #[values(Scaling::Unscaled, Scaling::Scaled)] scaling: Scaling,
    #[values(
        (complex_bessel_j as BesselSig, zbesj_fortran as BesselFortranSig),
        (complex_bessel_i as BesselSig, zbesi_fortran as BesselFortranSig),
        (complex_hankel1 as BesselSig , zbesh_fortran_first as BesselFortranSig),
        (complex_hankel2 as BesselSig , zbesh_fortran_second as BesselFortranSig),
        (complex_bessel_k as BesselSig, zbesk_fortran as BesselFortranSig),
        (complex_bessel_y as BesselSig, zbesy_fortran as BesselFortranSig),

    )]
    (rust_fn, fortran_fn): (BesselSig, BesselFortranSig),
) {
    let n = 1;
    for order in ORDERS {
        for re in Z_PARTS {
            for im in Z_PARTS {
                let z = Complex::new(re, im);
                check_against_fortran(order, z, scaling, n, rust_fn, fortran_fn);
            }
        }
    }
}

use complex_bessel::Error as RefError;
type Sig<T> = fn(f64, Complex<f64>) -> Result<Complex<f64>, T>;

fn hankel1(order: f64, z: Complex<f64>) -> Result<Complex<f64>, BesselError> {
    hankel(order, z, HankelKind::First)
}

fn hankel2(order: f64, z: Complex<f64>) -> Result<Complex<f64>, BesselError> {
    hankel(order, z, HankelKind::Second)
}

///Compute relative error between computed and reference complex values.
/// Returns None when both values are near zero (comparison meaningless).
/// math.hypot is used to avoid overflow on large intermediate values.
fn complex_bessel_test_relative_error(
    computed: Complex<f64>,
    ref_val: Complex<f64>,
) -> Option<f64> {
    let diff_re = computed.re - ref_val.re;
    let diff_im = computed.im - ref_val.im;
    let diff_mag = pymath::math::hypot(&[diff_re, diff_im]);

    let ref_mag = pymath::math::hypot(&[ref_val.re, ref_val.im]);

    // Both values near zero: comparison meaningless (e.g. analytic zeros,
    // f64 underflow at extreme orders). Skip these points entirely.
    if ref_mag < 1e-14 && diff_mag < 1e-14 {
        return None;
    }

    if ref_mag < 1e-300 {
        // Reference essentially zero but computed is not: real problem
        return Some(diff_mag);
    }

    Some(diff_mag / ref_mag)
}

/// Tests for a grid of values, comparing against the complex-bessel crate.
/// This is not intended to be exhaustive, but to catch any major issues with
/// the implementation across a wide range of inputs. The test_bessel_extremes
/// test is more exhaustive for edge cases.
/// Complex-bessel is used as a reference implementation here, as it implements
/// negative orders. It is tested (externally) against scipy and mpmath
/// for a range of values.
/// Complex-bessel-rs does implement negative orders, but does not appear
/// to do so correctly, or test them.
/// Airy tests are not needed here, as the fortran tests are exhaustive (and negative orders
/// are not relevant for Airy functions).
#[rstest]
fn test_bessel_grid_complex_besssel(
    #[values(
        (bessel_j as Sig<BesselError>, bessel_j_ref as Sig<RefError>),
        (bessel_i as Sig<BesselError>, bessel_i_ref as Sig<RefError>),
        (bessel_k as Sig<BesselError>, bessel_k_ref as Sig<RefError>),
        (bessel_y as Sig<BesselError>, bessel_y_ref as Sig<RefError>),
        (hankel1 as Sig<BesselError>, hankel1_ref as Sig<RefError>),
        (hankel2 as Sig<BesselError>, hankel2_ref as Sig<RefError>),
    )]
    (rust_fn, ref_fn): (Sig<BesselError>, Sig<RefError>),
) {
    for order in ORDERS {
        for re in Z_PARTS {
            for im in Z_PARTS {
                let z = Complex::new(re, im);
                let actual = rust_fn(order, z);
                let expected = ref_fn(order, z);
                if let Err(BesselError::InvalidInput { details: _ }) = actual {
                    assert!(
                        matches!(expected, Err(RefError::InvalidInput { .. })),
                        "Expected an InvalidInput error for order {order} and z {z}, but got {expected:?}"
                    );
                    return;
                }
                let actual = actual.unwrap();
                let expected = expected.unwrap();
                let rel_err = complex_bessel_test_relative_error(actual, expected);
                print!(
                    "\norder: {order}\nz: {z}\nactual: {actual:?}\nexpected: {expected:?}\nRelative Error: {rel_err:?}\n"
                );

                if let Some(msg) = check_complex_arrays_equal(&actual, &expected, &Vec::new()) {
                    panic!(
                        "Grid test failed\norder: {order}\nz: {z}
                    {msg}"
                    )
                }
                // below is the measure which is used by complex_bessel-test (though it only
                // prints the relative error: it doesn't assert.)
                // It's a very different measure to the one used in the fortran tests,
                // which is based on the number of matching significant digits.
                assert!(
                    rel_err.is_none() || rel_err.unwrap() < 1e-10,
                    "Relative error {rel_err:?} exceeds threshold for order {order} and z {z}",
                );
            }
        }
    }
}

#[rstest]
fn test_airy_grid_fortran(
    #[values(Scaling::Unscaled, Scaling::Scaled)] scaling: Scaling,
    #[values((sig_airy as BesselSig, sig_airy_fortran as BesselFortranSig),
        (sig_airyp as BesselSig, sig_airyp_fortran as BesselFortranSig),
        (sig_biry as BesselSig, sig_biry_fortran as BesselFortranSig),
        (sig_biryp as BesselSig, sig_biryp_fortran as BesselFortranSig))]
    (rust_fn, fortran_fn): (BesselSig, BesselFortranSig),
) {
    let dummy_order = 0.0;
    let n = 1;
    for re in Z_PARTS {
        for im in Z_PARTS {
            let z = Complex::new(re, im);
            check_against_fortran(dummy_order, z, scaling, n, rust_fn, fortran_fn);
        }
    }
}
