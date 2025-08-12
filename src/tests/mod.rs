use f64;
use num::Zero;
use rstest_reuse::{apply, template};

use approx::relative_eq;
use num::complex::Complex64;
use num::pow::Pow;
use rstest::rstest;

use crate::amos::bindings::zbesj_wrap;
use crate::amos::{BesselResult, MACHINE_CONSTANTS, MachineConsts, zbesj};
use crate::{BesselError, Scaling, bessel_j};
use complex_bessel_rs::bessel_i::bessel_i as bessel_i_ref;
use complex_bessel_rs::bessel_j::bessel_j as bessel_j_ref;

#[cfg(feature = "random_tests")]
mod random_tests;
mod test_gamma_ln;
mod test_machine_consts;

const TOLERANCE_MARGIN: f64 = 10_000.0;

#[template]
#[rstest]
#[case(4.0, 2.1, 0.0)]
#[case(5.0, 2.0001, 0.0)]
#[case(340.0, 35.0001, 0.0)]
#[case(407.3,-325.1, 635.2)]
#[case(465.0,-867.0, -448.0)] // 5
#[case(10.711871220659752, -6.89931119845653, -9.408182256887017)]
#[case(8.544950848779838, -8.645033163963603, 18.976439189605003,)]
#[case(21.04, 53.19, -40.77)]
#[case(4.0, 2.1, 0.0)]
#[case(5.0, 2.0001, 0.0)] // 10
#[case(340.0, 35.0001, 0.0)]
#[case(899.6,-35.7,317.8)]
#[case(531.0,-106.7,-16.0)]
#[case(531.0,-106.0,-16.0)]
#[case(433.0,16.874,-38.17)] //15
#[case(433.0,16.8,-38.17)]
#[case(311.2078694557452,-10.990141170168044,-25.70154097357056,)]
#[case(8.544950848779838, -8.645033163963603, 18.976439189605003,)]
#[case(17.5, 70.3, 37.4)]
#[case(13.337522865795481, -29.8266399174247, 17.66323218839807)] //20
#[case(5423.24, -7915.11, -3113.95)]
#[case(2213.0, -1813.0, -1033.0)]
#[case(5514.86274463943, -9489.650336481069, 4951.6909981261)]
#[case(2.74e-288, 6.33e-166, 7.53e-275)]
#[case(1.51e-150, -3.07e-118, 3.51e-42)] //25
#[case(2.637e-27, -4.01e-50, 0.0)]
#[case(4.0e-132, 0e0, 445.0)]
#[case(8714.0, 8904.0, -10.0)]
#[case(60.9, 246.2, -982.5)]
#[case(40.5, 1673.3, -4.0)] // 30
#[case(2634.5, -2634.5, 14.1)]
#[case(5.007e-14, 4.401331657952316e-5, -3.6e-6)]
#[case(1719.3, 920.1, 0.0)]
#[case(3.5695132850479827e3, -2.2313404290100934e3, 8.646324128723001e3)]
fn bessel_j_cases(#[case] order: f64, #[case] zr: f64, #[case] zi: f64) {}

#[apply(bessel_j_cases)]
fn test_bessel_j(#[case] order: f64, #[case] zr: f64, #[case] zi: f64) {
    let z = Complex64::new(zr, zi);
    let actual = bessel_j(order, z);

    if actual == Err(BesselError::NotYetImplemented) {
        todo!()
    }

    let expected = bessel_j_ref(order, z.into());
    if let Ok(actual) = actual {
        check_complex_arrays_equal(&actual, &expected.unwrap(), &Vec::new());
    } else {
        assert_eq!(actual.unwrap_err().error_code(), expected.unwrap_err())
    }
}

#[apply(bessel_j_cases)]
#[trace]
fn test_bessel_i(#[case] order: f64, #[case] zr: f64, #[case] zi: f64) {
    use crate::bessel_i;

    let z = Complex64::new(zr, zi);
    let actual = bessel_i(order, z);
    dbg!(&actual);
    if actual == Err(BesselError::NotYetImplemented) {
        todo!()
    }

    let expected = bessel_i_ref(order, z.into());
    if let Ok(actual) = actual {
        check_complex_arrays_equal(&actual, &expected.unwrap(), &Vec::new());
    } else {
        assert_eq!(actual.unwrap_err().error_code(), expected.unwrap_err())
    }
}

#[apply(bessel_j_cases)]
#[trace]
fn test_bessel_j_large_n_real(
    #[case] order: f64,
    #[case] zr: f64,
    #[case] _zi: f64,
    #[values(3, 4, 9, 100)] n: usize,
    #[values(Scaling::Unscaled, Scaling::Scaled)] scaling: Scaling,
    // #[values(3)] n: usize,
    // #[values(Scaling::Scaled)] scaling: Scaling,
) {
    // ignores the zi input
    check_against_fortran(order, zr.into(), scaling, n);
}

#[apply(bessel_j_cases)]
fn test_bessel_j_large_n_complex(
    #[case] order: f64,
    #[case] zr: f64,
    #[case] zi: f64,
    #[values(3, 4, 9, 100)] n: usize,
    #[values(Scaling::Unscaled, Scaling::Scaled)] scaling: Scaling,
    // #[values(9)] n: usize,
    // #[values(Scaling::Unscaled)] scaling: Scaling,
) {
    let z = Complex64::new(zr, zi);
    check_against_fortran(order, z, scaling, n);
}

#[rstest]
fn test_bessel_j_extremes(
    #[values(0.0, f64::EPSILON, f64::MIN_POSITIVE, 1.0, f64::MAX)] order: f64,
    #[values(
        f64::MIN,
        -f64::EPSILON,
        -f64::MIN_POSITIVE,
        0.0,
        f64::MIN_POSITIVE,
        f64::EPSILON,
        1.0,
        f64::MAX
    )]
    zr: f64,
    #[values(
        f64::MIN,
        -f64::EPSILON,
        -f64::MIN_POSITIVE,
        0.0,
        f64::MIN_POSITIVE,
        f64::EPSILON,
        1.0,
        f64::MAX
        )]
    zi: f64,
    #[values(Scaling::Unscaled, Scaling::Scaled)] scaling: Scaling,
) {
    let n = 9;
    let z = Complex64::new(zr, zi);
    check_against_fortran(order, z, scaling, n);
}

fn check_against_fortran(order: f64, z: Complex64, scaling: Scaling, n: usize) {
    let actual = zbesj(z, order, scaling, n);
    if let Err(ref err) = actual {
        if *err == BesselError::NotYetImplemented {
            // return;
            panic!("NotYetImplemented should not occur")
        }
    }

    let (cy, nz, ierr) = zbesj_fortran(order, z, scaling, n);
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
trait IntoComplexVec: Clone {
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

fn check_complex_arrays_equal(
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

fn zbesj_loop(order: f64, z: Complex64, scaling: Scaling, n: usize) -> BesselResult {
    let mut y = vec![Complex64::zero(); n];
    let mut nz = 0;
    for i in 0..n {
        let (yi, nzi) = zbesj(z, order + i as f64, scaling, 1)?;
        y[i] = yi[0];
        nz += nzi;
    }
    return Ok((y, nz));
}

fn zbesj_fortran_loop(
    order: f64,
    z: Complex64,
    scaling: Scaling,
    n: usize,
) -> (Vec<Complex64>, usize, i32) {
    let mut y = vec![Complex64::zero(); n];
    let mut nz = 0;
    for i in 0..n {
        let (yi, nzi, ierr) = zbesj_fortran(order + i as f64, z, scaling, 1);
        if ierr != 0 {
            return (y, nz, ierr);
        }
        y[i] = yi[0];
        nz += nzi;
    }
    (y, nz, 0)
}

fn zbesj_fortran(
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
