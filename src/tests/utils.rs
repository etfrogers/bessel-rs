use approx::relative_eq;
use num::{complex::Complex64, pow::Pow};

use crate::{
    BesselError, Scaling,
    amos::{MACHINE_CONSTANTS, zbesj},
    tests::{TOLERANCE_MARGIN, zbesj_fortran_loop, zbesj_loop},
};

pub fn check_against_fortran(order: f64, z: Complex64, scaling: Scaling, n: usize) {
    let actual = zbesj(z, order, scaling, n);
    if let Err(ref err) = actual {
        if *err == BesselError::NotYetImplemented {
            // return;
            panic!("NotYetImplemented should not occur")
        }
    }

    let (cy, nz, ierr) = super::zbesj_fortran(order, z, scaling, n);
    let (cy_loop_fort, _, _) = zbesj_fortran_loop(order, z, scaling, n);

    let fail = |reason: &str| -> () {
        let (cy_loop_rust, _) = zbesj_loop(order, z, scaling, n).unwrap();
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
            if let Some(reason) = check_complex_arrays_equal(&actual.0, &cy, &cy_loop_fort) {
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
                if nz != actual_nz.try_into().unwrap() {
                    fail("Failed for mismatched nz value");
                }
                if let Some(reason) = check_complex_arrays_equal(actual_y, &cy, &cy_loop_fort) {
                    fail(&reason)
                }
            }
        }
    }
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

pub fn check_complex_arrays_equal(
    actual: &impl IntoComplexVec,
    expected: &impl IntoComplexVec,
    reference: &impl IntoComplexVec,
) -> Option<String> {
    let actual = actual.clone().into_vec();
    let expected = expected.clone().into_vec();
    let reference = reference.clone().into_vec();
    let exp_error = MACHINE_CONSTANTS.abs_error_tolerance;

    for (i, (&act, exp)) in actual.iter().zip(expected).enumerate() {
        let ref_val = reference.get(i);
        let tolerances = Tolerances::new(act, exp, ref_val, exp_error);
        if !complex_relative_eq(act, exp, &tolerances) {
            let (actual_error, relative_error) = abs_rel_errors_cmplx(act, exp);
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
                tolerances.tol_re * TOLERANCE_MARGIN,
                actual_error.re,
                relative_error.re,
                tolerances.tol_im * TOLERANCE_MARGIN,
                actual_error.im,
                relative_error.im,
            ));
        };
    }
    None
}

struct Tolerances {
    max_relative: f64,
    magnitude_diff: f64,
    correction: f64,
    tol_re: f64,
    tol_im: f64,
    used_ref_re: bool,
    used_ref_im: bool,
}

impl Tolerances {
    fn new(
        actual: Complex64,
        expected: Complex64,
        reference: Option<&Complex64>,
        max_relative: f64,
    ) -> Self {
        let max_im = actual.im.abs().max(expected.im.abs());
        let max_re = actual.re.abs().max(expected.re.abs());
        let magnitude_diff = (max_re.log10() - max_im.log10()).abs();
        let limit = 1.0 / max_relative;
        let correction = 10.0.pow(magnitude_diff).min(limit);

        let (mut tol_re, mut tol_im) = if max_re > max_im {
            (max_relative, max_relative * correction)
        } else {
            (max_relative * correction, max_relative)
        };

        let mut used_ref_re = false;
        let mut used_ref_im = false;
        if let Some(ref_val) = reference {
            let (_, ref_diffs) = abs_rel_errors_cmplx(expected, *ref_val);
            let re_rel_diff = ref_diffs.re;
            let im_rel_diff = ref_diffs.im;

            used_ref_re = re_rel_diff > tol_re;
            if used_ref_re {
                tol_re = re_rel_diff;
            }
            used_ref_im = im_rel_diff > tol_im;
            if used_ref_im {
                tol_im = im_rel_diff;
            }
        }
        Self {
            max_relative,
            magnitude_diff,
            correction,
            tol_re,
            tol_im,
            used_ref_re,
            used_ref_im,
        }
    }
}

fn complex_relative_eq(a: Complex64, b: Complex64, tolerances: &Tolerances) -> bool {
    if relative_eq!(
        a,
        b,
        max_relative = TOLERANCE_MARGIN * tolerances.max_relative
    ) {
        return true;
    }
    let (_, rel_e_re) = abs_rel_errors(a.re, b.re);
    let (_, rel_e_im) = abs_rel_errors(a.im, b.im);

    rel_e_re < TOLERANCE_MARGIN * tolerances.tol_re
        && rel_e_im < TOLERANCE_MARGIN * tolerances.tol_im
}

fn abs_rel_errors(a: f64, b: f64) -> (f64, f64) {
    let abs_e = (a - b).abs();
    let rel_e = abs_e / a.abs().max(b.abs());
    (abs_e, rel_e)
}

fn abs_rel_errors_cmplx(a: Complex64, b: Complex64) -> (Complex64, Complex64) {
    let (abs_e_r, rel_e_r) = abs_rel_errors(a.re, b.re);
    let (abs_e_i, rel_e_i) = abs_rel_errors(a.im, b.im);
    (
        Complex64::new(abs_e_r, abs_e_i),
        Complex64::new(rel_e_r, rel_e_i),
    )
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
                abs_rel_errors_cmplx(fort, c2[i]),
                abs_rel_errors_cmplx(fort, c3[i]),
                abs_rel_errors_cmplx(fort, c4[i]),
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
