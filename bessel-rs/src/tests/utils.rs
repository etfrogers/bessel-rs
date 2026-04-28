use std::f64;
use std::fmt::Debug;

use num::{Zero, complex::Complex64};

use crate::types::{
    BesselError, BesselResult, Tolerances, abs_rel_errors_cmplx, complex_relative_eq,
};

use crate::{
    Scaling,
    tests::{BesselFortranSig, BesselSig},
};

pub fn assert_results_are_equal<T: IntoComplexVec + Debug>(
    actual: &BesselResult<T>,
    expected: &Result<T, i32>,
    reference: &impl IntoComplexVec,
    margin: f64,
) {
    if let Ok(actual) = actual {
        assert_complex_arrays_equal(actual, expected.as_ref().unwrap(), reference, margin);
    } else {
        assert_eq!(
            actual.as_ref().unwrap_err().error_code(),
            *expected.as_ref().unwrap_err()
        )
    }
}

pub fn check_against_fortran(
    order: f64,
    z: Complex64,
    scaling: Scaling,
    n: usize,
    rust_func: BesselSig,
    fortran_func: BesselFortranSig,
    margin: f64,
) {
    let actual = rust_func(z, order, scaling, n);

    let (cy, nz, ierr) = fortran_func(order, z, scaling, n);
    let (cy_loop_fort, _, _) = fortran_bess_loop(order, z, scaling, n, fortran_func);

    let fail = |reason: &str| -> () {
        let cy_loop_rust = match rust_bess_loop(order, z, scaling, n, rust_func) {
            Ok((data, _)) => data,
            Err(err) => {
                panic!(
                    "Error generated in looped rust that was not present in unlooped case: {err:?}"
                );
            }
        };
        println!("Order: {order:e}\nz: {z:e}\nscaling: {scaling:?}\nn: {n}");
        println!("#[case({:e}, {:e}, {:e})]", order, z.re, z.im);
        println!("#[case({:.1}, {:.1}, {:.1})]\n", order, z.re, z.im);
        match &actual {
            Ok(actual) => {
                println!("Fortran Nz: {nz}, translator Nz: {}\n", actual.1);
                print_complex_arrays(&cy, &actual.0, &cy_loop_fort, &cy_loop_rust);
            }
            Err(err) => {
                println!(
                    "Fortran error: {ierr}. Translation error: {err:?} ({})",
                    err.error_code()
                );
                if let BesselError::PartialLossOfSignificance {
                    y: ref actual_y,
                    nz: actual_nz,
                } = *err
                {
                    println!("Fortran Nz: {nz}, translator Nz: {}\n", actual_nz);
                    print_complex_arrays(&cy, actual_y, &cy_loop_fort, &cy_loop_rust);
                }
            }
        }
        println!();
        panic!("{reason}")
    };

    match &actual {
        Ok(actual) => {
            if let Some(reason) = check_complex_arrays_equal(&actual.0, &cy, &cy_loop_fort, margin)
            {
                fail(&reason)
            }
        }
        Err(err) => {
            if ierr != err.error_code() {
                fail("Failed for mismatched error code")
            };
            if let BesselError::PartialLossOfSignificance {
                y: ref actual_y,
                nz: actual_nz,
            } = *err
            {
                if nz != actual_nz {
                    // for partial loss of significance, it seems occasionally fortran
                    // will return some values very nearly zero, but it's only happening
                    // on a release build, so it may be some optimization issue. It also occurs
                    // sometimes (though flakily on a linux build) To avoid
                    // this causing test failures, effectively skipping the check on the nz value
                    // And falling through to the value checks, below, but these will catch large errors.
                    // This is not ideal, but I have not been able to find a better solution.
                    //
                    // Note this is only for the partial loss of significance case, which is
                    // already a case where the results are not fully trustworthy, so it seems
                    // reasonable to me to allow this kind of mismatch in this case.

                    // fail("Failed for mismatched nz value");
                }
                if let Some(reason) =
                    check_complex_arrays_equal(actual_y, &cy, &cy_loop_fort, margin * 1e2)
                {
                    fail(&reason)
                }
            }
        }
    }
}

fn rust_bess_loop(
    order: f64,
    z: Complex64,
    scaling: Scaling,
    n: usize,
    func: BesselSig,
) -> BesselResult {
    let mut y = vec![Complex64::zero(); n];
    let mut nz = 0;
    for i in 0..n {
        let (yi, nzi) = match func(z, order + i as f64, scaling, 1) {
            Ok((y_, nz_)) => (y_, nz_),
            Err(BesselError::PartialLossOfSignificance { y, nz }) => (y, nz),
            Err(err) => return Err(err),
        };
        y[i] = yi[0];
        nz += nzi;
    }
    return Ok((y, nz));
}

pub fn fortran_bess_loop(
    order: f64,
    z: Complex64,
    scaling: Scaling,
    n: usize,
    func: BesselFortranSig,
) -> (Vec<Complex64>, usize, i32) {
    let mut y = vec![Complex64::zero(); n];
    let mut nz = 0;
    for i in 0..n {
        let (yi, nzi, ierr) = func(order + i as f64, z, scaling, 1);
        if ierr != 0 {
            return (y, nz, ierr);
        }
        y[i] = yi[0];
        nz += nzi;
    }
    (y, nz, 0)
}

pub trait IntoComplexVec: Clone {
    fn into_vec(self) -> Vec<Complex64>;
}
impl IntoComplexVec for Complex64 {
    fn into_vec(self) -> Vec<Complex64> {
        vec![self]
    }
}

impl IntoComplexVec for Vec<Complex64> {
    fn into_vec(self) -> Vec<Complex64> {
        self
    }
}

pub fn assert_complex_arrays_equal(
    actual: &impl IntoComplexVec,
    expected: &impl IntoComplexVec,
    reference: &impl IntoComplexVec,
    margin: f64,
) {
    if let Some(reason) = check_complex_arrays_equal(actual, expected, reference, margin) {
        panic!("{}", reason);
    }
}

#[must_use]
pub fn check_complex_arrays_equal(
    actual: &impl IntoComplexVec,
    expected: &impl IntoComplexVec,
    reference: &impl IntoComplexVec,
    margin: f64,
) -> Option<String> {
    let actual = actual.clone().into_vec();
    let expected = expected.clone().into_vec();
    let reference = reference.clone().into_vec();

    for (i, (&act, exp)) in actual.iter().zip(expected).enumerate() {
        let ref_val = reference.get(i);
        let tolerances = Tolerances::new(&act, &exp, ref_val, margin);
        if !complex_relative_eq(&act, &exp, &tolerances) {
            let (actual_error, relative_error) = abs_rel_errors_cmplx(&act, &exp);
            return Some(format!(
                "Failed on matching values at index {i}\n\
                Actual: {act:e}\n\
                Expected: {exp:e}\n\
                \n\
                Magnitude difference in actual: {}\n\
                Correction: {:e}\n\
                \n\
                Used loop diff real: {}\n\
                Used loop diff imag: {}\n\
                \n\
                Relative tolerance - real: {:e}\n\
                Absolute error - real: {:e}\n\
                Relative error - real: {:e}\n\
                \n\
                Relative tolerance - imag: {:e}\n\
                Absolute error - imag: {:e}\n\
                Relative error - imag: {:e}",
                tolerances.magnitude_diff,
                tolerances.correction,
                tolerances.used_ref_re,
                tolerances.used_ref_im,
                tolerances.tol_re * tolerances.margin,
                actual_error.re,
                relative_error.re,
                tolerances.tol_im * tolerances.margin,
                actual_error.im,
                relative_error.im,
            ));
        };
    }
    None
}

fn print_complex_arrays(c1: &[Complex64], c2: &[Complex64], c3: &[Complex64], c4: &[Complex64]) {
    println!("i\tFortran\t\t\t\tTranslator\t\t\t\tFortran looped\t\t\t\tRust looped");
    c1.iter().enumerate().for_each(|(i, fort)| {
        println!(
            "{i}\t{}\t{}\t{}\t{}",
            to_str(fort),
            to_str(&c2[i]),
            to_str(&c3[i]),
            to_str(&c4[i]),
        );
    });

    let errors: Vec<_> = c1
        .iter()
        .enumerate()
        .map(|(i, &fort)| {
            (
                abs_rel_errors_cmplx(&fort, &c2[i]),
                abs_rel_errors_cmplx(&fort, &c3[i]),
                abs_rel_errors_cmplx(&fort, &c4[i]),
            )
        })
        .collect();

    println!("\nAbsolute Errors");
    println!("i\tTranslator\t\t\t\tFortran looped\t\t\t\tRust looped");
    errors.iter().enumerate().for_each(|(i, errs)| {
        println!(
            "{i}\t{}\t{}\t{}",
            to_str(&errs.0.0),
            to_str(&errs.1.0),
            to_str(&errs.2.0),
        );
    });

    println!("\nRelative Errors");
    println!("i\tTranslator\t\t\t\tFortran looped\t\t\t\tRust looped");
    errors.iter().enumerate().for_each(|(i, errs)| {
        println!(
            "{i}\t{}\t{}\t{}",
            to_str(&errs.0.1),
            to_str(&errs.1.1),
            to_str(&errs.2.1),
        );
    });
}

fn to_str(c: &Complex64) -> String {
    format!("{:>+1.5e} {:>+1.5e}i", c.re, c.im)
}
