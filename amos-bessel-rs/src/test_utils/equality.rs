use approx::relative_eq;
use num::{Complex, pow::Pow};
use std::fmt::Debug;

use crate::{amos::MACHINE_CONSTANTS, types::BesselResult};

pub(crate) struct Tolerances {
    pub max_relative: f64,
    pub magnitude_diff: f64,
    pub correction: f64,
    pub tol_re: f64,
    pub tol_im: f64,
    pub used_ref_re: bool,
    pub used_ref_im: bool,
    pub margin: f64,
}

impl Tolerances {
    pub(crate) fn new(
        actual: &Complex<f64>,
        expected: &Complex<f64>,
        reference: Option<&Complex<f64>>,
        margin: f64,
    ) -> Self {
        let max_relative = MACHINE_CONSTANTS.abs_error_tolerance;

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
            let (_, ref_diffs) = abs_rel_errors_cmplx(expected, ref_val);
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
            margin,
        }
    }
}

pub(crate) fn complex_relative_eq(
    a: &Complex<f64>,
    b: &Complex<f64>,
    tolerances: &Tolerances,
) -> bool {
    let abs_floor = MACHINE_CONSTANTS.abs_error_tolerance * 1e2;

    if relative_eq!(
        a,
        b,
        epsilon = abs_floor,
        max_relative = tolerances.margin * tolerances.max_relative
    ) {
        return true;
    }

    // AMOS Component-Wise Check
    // This evaluates the real and imaginary parts using your dynamically
    // generated tolerances, but safely ignores differences below the noise floor.
    let re_matches = relative_eq!(
        a.re,
        b.re,
        epsilon = abs_floor,
        max_relative = tolerances.margin * tolerances.tol_re
    );

    let im_matches = relative_eq!(
        a.im,
        b.im,
        epsilon = abs_floor,
        max_relative = tolerances.margin * tolerances.tol_im
    );

    re_matches && im_matches
}

fn abs_rel_errors(a: f64, b: f64) -> (f64, f64) {
    let abs_e = (a - b).abs();
    let rel_e = abs_e / a.abs().max(b.abs());
    (abs_e, rel_e)
}

pub(crate) fn abs_rel_errors_cmplx(
    a: &Complex<f64>,
    b: &Complex<f64>,
) -> (Complex<f64>, Complex<f64>) {
    let (abs_e_r, rel_e_r) = abs_rel_errors(a.re, b.re);
    let (abs_e_i, rel_e_i) = abs_rel_errors(a.im, b.im);
    (
        Complex::<f64>::new(abs_e_r, abs_e_i),
        Complex::<f64>::new(rel_e_r, rel_e_i),
    )
}

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

pub trait IntoComplexVec: Clone {
    fn into_vec(self) -> Vec<Complex<f64>>;
}
impl IntoComplexVec for Complex<f64> {
    fn into_vec(self) -> Vec<Complex<f64>> {
        vec![self]
    }
}

impl IntoComplexVec for Vec<Complex<f64>> {
    fn into_vec(self) -> Vec<Complex<f64>> {
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

pub fn print_complex_arrays(
    c1: &[Complex<f64>],
    c2: &[Complex<f64>],
    c3: &[Complex<f64>],
    c4: &[Complex<f64>],
) {
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

fn to_str(c: &Complex<f64>) -> String {
    format!("{:>+1.5e} {:>+1.5e}i", c.re, c.im)
}
