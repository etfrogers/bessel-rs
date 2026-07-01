use approx::relative_eq;
use num::{Complex, Float};
use std::fmt::{Debug, LowerExp};

use super::DiagnosticBesselFloat;
use crate::{
    BesselError,
    types::{BesselFloat, BesselValues},
};

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
    pub(crate) fn new<T: BesselFloat>(
        actual: &Complex<T>,
        expected: &Complex<f64>,
        reference: Option<&Complex<f64>>,
        margin: f64,
    ) -> Self {
        let max_relative = T::MACHINE_CONSTANTS.abs_error_tolerance.to_f64().unwrap();

        let max_im = actual.im.abs().to_f64().unwrap().max(expected.im.abs());
        let max_re = actual.re.abs().to_f64().unwrap().max(expected.re.abs());
        let magnitude_diff = (max_re.log10() - max_im.log10()).abs();
        let limit = 1.0 / max_relative;
        let correction = 10.0.powf(magnitude_diff).min(limit);

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

fn errors_eq<T: DiagnosticBesselFloat>(
    lhs: &BesselError<T>,
    rhs: &BesselError<f64>,
    margin: f64,
) -> bool {
    match (lhs, rhs) {
        (
            BesselError::InvalidInput { details: l_details },
            BesselError::InvalidInput { details: r_details },
        ) => l_details == r_details,
        (
            BesselError::PartialLossOfSignificance { y: l_y, nz: l_nz },
            BesselError::PartialLossOfSignificance { y: r_y, nz: r_nz },
        ) => check_complex_arrays_equal(l_y, r_y, &vec![], margin).is_none() && l_nz == r_nz,
        (
            BesselError::ComplexOutputForRealInput { output: l_output },
            BesselError::ComplexOutputForRealInput { output: r_output },
        ) => {
            let tolerances = &Tolerances::new(l_output, r_output, None, margin);
            complex_relative_eq(l_output, r_output, tolerances)
        }

        _ => core::mem::discriminant(&lhs.to_f32()) == core::mem::discriminant(&rhs.to_f32()),
    }
}

pub(crate) fn complex_relative_eq<T: DiagnosticBesselFloat>(
    a: &Complex<T>,
    b: &Complex<f64>,
    tolerances: &Tolerances,
) -> bool {
    let abs_floor = f64::MACHINE_CONSTANTS.abs_error_tolerance * 1e2;

    if relative_eq!(
        &a.to_c64(),
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
        a.re.to_f64().unwrap(),
        b.re,
        epsilon = abs_floor,
        max_relative = tolerances.margin * tolerances.tol_re
    );

    let im_matches = relative_eq!(
        a.im.to_f64().unwrap(),
        b.im,
        epsilon = abs_floor,
        max_relative = tolerances.margin * tolerances.tol_im
    );

    re_matches && im_matches
}

fn abs_rel_errors<T: BesselFloat>(a: T, b: f64) -> (T, T) {
    let b = T::from_f64(b);

    let abs_e = (a - b).abs();
    let rel_e = abs_e / a.abs().max(b.abs());
    (abs_e, rel_e)
}

pub(crate) fn abs_rel_errors_cmplx<T: BesselFloat>(
    a: &Complex<T>,
    b: &Complex<f64>,
) -> (Complex<T>, Complex<T>) {
    let (abs_e_r, rel_e_r) = abs_rel_errors(a.re, b.re);
    let (abs_e_i, rel_e_i) = abs_rel_errors(a.im, b.im);
    (
        Complex::new(abs_e_r, abs_e_i),
        Complex::new(rel_e_r, rel_e_i),
    )
}

pub trait ComplexConversions {
    fn to_c32(&self) -> Complex<f32>;
    fn to_c64(&self) -> Complex<f64>;
}

impl<T: BesselFloat> ComplexConversions for Complex<T> {
    fn to_c32(&self) -> Complex<f32> {
        Complex::new(self.re.to_f32().unwrap(), self.im.to_f32().unwrap())
    }

    fn to_c64(&self) -> Complex<f64> {
        Complex::new(self.re.to_f64().unwrap(), self.im.to_f64().unwrap())
    }
}

pub fn assert_results_are_equal_floats<T: DiagnosticBesselFloat>(
    actual: &Result<BesselValues<T>, BesselError<T>>,
    expected: &Result<BesselValues<f64>, BesselError<f64>>,
    margin: f64,
) {
    match (actual, expected) {
        (Ok(actual_vals), Ok(expected_vals)) => {
            // if let (Ok(actual_vals), Ok(expected_vals)) = (actual, expected) {
            if actual_vals.1 > 0 || expected_vals.1 > 0 {
                // If either calculation experienced an underflow (nz > 0),
                // f32 and f64 will completely diverge. Skip comparison.
                return;
            }

            let actual_vec = &actual_vals.0;
            let expected_vec = &expected_vals.0;
            assert_complex_arrays_equal(actual_vec, expected_vec, &vec![], margin);
        }
        (Err(BesselError::Overflow), _) => {
            // Overflow for f32 does not imply overflow for f64
            return;
        }
        (Err(BesselError::LossOfSignificance), Err(BesselError::Overflow)) => {
            // As the loss of significance check is early in the code path, it may prevent a later
            // Overflow form occurring.
            return;
        }

        (
            Err(BesselError::LossOfSignificance),
            Err(BesselError::PartialLossOfSignificance { y: _, nz: _ }),
        ) => {
            // Possible for f32 to lose all siginifcance, and f64 to retain some. That's OK.
            return;
        }
        (
            Err(BesselError::PartialLossOfSignificance { y: actual_y, nz: _ }),
            Err(BesselError::PartialLossOfSignificance {
                y: expected_y,
                nz: _,
            }),
        ) => {
            // If they both lose significance, it is unlikley that the values in there will be the same, but that's ok.
            // Just check that they are the same order of magnitude
            // margin is used as err < margin * abs_error_tolerance
            // therefore margin to make an order-of-magnitude check
            println!("Both lost significance: \n{:?}\n {:?}", actual, expected);
            let oom_margin = 1.0 / T::MACHINE_CONSTANTS.abs_error_tolerance.to_f64().unwrap();
            assert_complex_arrays_equal(actual_y, expected_y, &vec![], oom_margin);
        }
        (Err(BesselError::PartialLossOfSignificance { y: actual_y, nz: _ }), Ok(expected_vals)) => {
            // In this case f32 has lost significance, but f64 hasn't. Again, check that answers are within
            // an order of magnitude
            println!("One lost significance: \n{:?}\n {:?}", actual, expected);
            let oom_margin = 1.0 / T::MACHINE_CONSTANTS.abs_error_tolerance.to_f64().unwrap();
            assert_complex_arrays_equal(actual_y, &expected_vals.0, &vec![], oom_margin);
        }

        (Err(actual_err), Err(expected_err)) => {
            assert!(
                errors_eq(actual_err, expected_err, margin),
                "{actual:?} != {expected:?}",
            )
        }
        _ => {
            panic!(
                "Mismatched result types: actual={:?}, expected={:?}",
                actual,   //.as_ref()//.map(|_| "Ok"),
                expected  //.as_ref()//.map(|_| "Ok")
            );
        }
    }
}

pub fn assert_results_are_equal<T: DiagnosticBesselFloat>(
    actual: &Result<impl IntoComplexVec<T> + Debug, BesselError<T>>,
    expected: &Result<impl IntoComplexVec<f64> + Debug, i32>,
    reference: &impl IntoComplexVec<f64>,
    margin: f64,
)
where
// Complex<T>: IntoComplexVec<T>,
{
    if let Ok(actual) = actual {
        // The single-output function unwrap PartialLossOfSignificance errors, and return Ok(),
        // This corresponds to the ref function returning ierr=3, we need to allow expected
        // to be error code 3.
        // Not that the reference functions don't return any values to check against if
        // ierr != 0
        // Tests against Fortran should allow test of values.
        if expected.as_ref().is_err_and(|err| *err == 3) {
            return;
        }
        assert_complex_arrays_equal(actual, expected.as_ref().unwrap(), reference, margin);
    } else {
        assert_eq!(
            actual.as_ref().unwrap_err().error_code(),
            *expected.as_ref().unwrap_err()
        )
    }
}

pub trait IntoComplexVec<T>: Clone {
    fn into_vec(self) -> Vec<Complex<T>>;
}

impl<T: BesselFloat> IntoComplexVec<T> for Complex<T> {
    fn into_vec(self) -> Vec<Complex<T>> {
        vec![self]
    }
}

impl<T: BesselFloat> IntoComplexVec<T> for Vec<Complex<T>> {
    fn into_vec(self) -> Vec<Complex<T>> {
        self
    }
}

pub fn assert_complex_arrays_equal<T: DiagnosticBesselFloat>(
    actual: &impl IntoComplexVec<T>,
    expected: &impl IntoComplexVec<f64>,
    reference: &impl IntoComplexVec<f64>,
    margin: f64,
) {
    if let Some(reason) = check_complex_arrays_equal(actual, expected, reference, margin) {
        panic!("{}", reason);
    }
}

#[must_use]
pub fn check_complex_arrays_equal<T1: DiagnosticBesselFloat>(
    actual: &impl IntoComplexVec<T1>,
    expected: &impl IntoComplexVec<f64>,
    reference: &impl IntoComplexVec<f64>,
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
                Relative error - imag: {:e}\n\
                \n\
                Margin: {:e}",
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
                tolerances.margin,
            ));
        };
    }
    None
}

// Note that fortran values are always f64, but rust can vary.
pub fn print_complex_arrays<T: DiagnosticBesselFloat>(
    exp: &[Complex<f64>],
    actual: &[Complex<T>],
    exp_looped: &[Complex<f64>],
    actual_looped: &[Complex<T>],
) {
    println!("i\tFortran\t\t\t\tTranslator\t\t\t\tFortran looped\t\t\t\tRust looped");
    exp.iter().enumerate().for_each(|(i, fort)| {
        println!(
            "{i}\t{}\t{}\t{}\t{}",
            to_str(fort),
            to_str(&actual[i]),
            to_str(&exp_looped[i]),
            to_str(&actual_looped[i]),
        );
    });

    let errors: Vec<_> = exp
        .iter()
        .enumerate()
        .map(|(i, &fort)| {
            (
                abs_rel_errors_cmplx(&actual[i], &fort),
                abs_rel_errors_cmplx(&exp_looped[i], &fort),
                abs_rel_errors_cmplx(&actual_looped[i], &fort),
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

fn to_str<T: Float + LowerExp>(c: &Complex<T>) -> String {
    format!("{:>+1.5e} {:>+1.5e}i", c.re, c.im)
}
