use approx::assert_relative_eq;
use num::complex::Complex64;
use rstest::rstest;

use crate::amos::bindings::zbesj_wrap;
use crate::amos::{gamma_ln, zbesj};
use crate::{BesselError, GammaError, Scaling, bessel_j};
use complex_bessel_rs::bessel_j::bessel_j as bessel_j_ref;

const RANDOM_LIMIT: f64 = 10_000.0;

#[test]
fn test_gamma_ln_hard_coded() {
    for i in 1..=100 {
        let f = i as f64;
        let actual = gamma_ln(f).unwrap();
        let expected = f.gamma().ln();
        assert_relative_eq!(actual, expected)
    }
}

#[test]
fn test_gamma_ln() {
    for f in [0.00000000001, 2.2, 0.5, 2.0_f64.sqrt(), 101.0] {
        // large values cause the "expected" calculation to overflow: the fortran version seems to work!
        let actual = gamma_ln(f).unwrap();
        let expected = f.gamma().ln();
        assert_relative_eq!(actual, expected, epsilon = 1e-10)
    }
}

#[rstest]
fn test_gamma_ln_negative() {
    for f in [0.00000000001, 2.2, 0.5, 2.0_f64.sqrt(), 101.0] {
        // large values cause the "expected" calculation to overflow: the fortran version seems to work!
        let actual = gamma_ln(-f);
        assert_eq!(actual, Err(GammaError::ZLessThanZero))
    }
}

#[rstest]
#[trace]
#[case(4.0, 2.1, 0.0)] // z_power_series
#[case(5.0, 2.0001, 0.0)] // z_power_series
#[case(340.0, 35.0001, 0.0)]
// z_power_series, iflag = true
// #[case(407.332478234955,-325.1302407029058, 635.2191950523381)] //not yet implemented
// #[case(465.45726205167904,-867.9390777060459, -448.2782267806473)] //not yet implemented
#[case(10.711871220659752, -6.89931119845653, -9.408182256887017)]
#[case(8.544950848779838, -8.645033163963603, 18.976439189605003,)]
#[case(21.044832034798453, 53.19358867279834, -40.77080989865806)]
fn test_bessel_j(#[case] order: f64, #[case] zr: f64, #[case] zi: f64) {
    let z = Complex64::new(zr, zi);
    let actual = bessel_j(order, z);

    if actual == Err(BesselError::NotYetImplemented) {
        todo!()
    }

    let expected = bessel_j_ref(order, z.into());
    if let Ok(actual) = actual {
        assert_relative_eq!(
            actual,
            expected.unwrap(),
            epsilon = 1e-10,
            max_relative = 1e-10
        )
    } else {
        assert_eq!(actual.unwrap_err().error_code(), expected.unwrap_err())
    }
}

#[rstest]
fn test_bessel_j_random() {
    for _ in 0..1000000 {
        let order = rand::random_range(std::f64::EPSILON..RANDOM_LIMIT);
        let zr = random_val();
        let zi = random_val();
        let z = Complex64::new(zr, zi);
        // dbg!(order, &z);
        // println!("#[case({}, {}, {})]", order, z.re, z.im);
        let ans = bessel_j(order, z);
        let expected = bessel_j_ref(order, z);
        if let Ok(actual) = ans {
            assert_relative_eq!(
                actual,
                expected.unwrap(),
                epsilon = 1e-10,
                max_relative = 1e-10
            )
        } else {
            let err = ans.unwrap_err();
            if err == BesselError::NotYetImplemented {
                continue;
            }
            // dbg!(&err, &expected);
            assert_eq!(err.error_code(), expected.unwrap_err());
        }
    }
}

#[rstest]
fn test_bessel_j_random_negative() {
    todo!("Add negative values")
}

#[rstest]
#[trace]
#[case(4.0, 2.1)] // z_power_series
#[case(5.0, 2.0001)] // z_power_series
#[case(340.0, 35.0001)] // z_power_series, iflag = true
fn test_bessel_j_large_n_real(
    #[case] order: f64,
    #[case] z: f64,
    #[values(3, 4, 9, 100)] n: usize,
    #[values(Scaling::Unscaled, Scaling::Scaled)] scaling: Scaling,
) {
    check_against_fortran(order, z.into(), scaling, n);
}

#[rstest]
#[case(899.6186771978487,-35.73821920707451,317.85710422430134)]
#[case(531.0,-106.7,-16.0)]
#[case(531.0,-106.0,-16.0)]
#[case(433.0,16.874,-38.17)]
#[case(433.0,16.8,-38.17)]
#[case(311.2078694557452,-10.990141170168044,-25.70154097357056,)]
#[case(8.544950848779838, -8.645033163963603, 18.976439189605003,)]
#[case(17.556977911963312, 70.34021294440504, 37.416997183283456)]
#[case(13.337522865795481, -29.8266399174247, 17.66323218839807)]
#[case(5423.246927434604, -7915.1124370237285, -3113.950242590895)]
#[case(2213.61988214781, -1813.3484572476455, -1033.3403805479065)]
#[trace]
fn test_bessel_j_large_n_complex(
    #[case] order: f64,
    #[case] zr: f64,
    #[case] zi: f64,
    #[values(3, 4, 9, 100)] n: usize,
    #[values(Scaling::Unscaled, Scaling::Scaled)] scaling: Scaling,
) {
    let z = Complex64::new(zr, zi);
    check_against_fortran(order, z, scaling, n);
}

enum NumType {
    Real,
    Imaginary,
    Complex,
}

fn random_val() -> f64 {
    rand::random_range(-RANDOM_LIMIT..RANDOM_LIMIT)
}

#[rstest]
fn test_bessel_j_large_n_random(
    #[values(Scaling::Unscaled, Scaling::Scaled)] scaling: Scaling,
    #[values(NumType::Real, NumType::Imaginary, NumType::Complex)] num_type: NumType,
) {
    let n = 9;
    for _ in 0..100000 {
        let order = rand::random_range(std::f64::EPSILON..RANDOM_LIMIT);
        let (zr, zi) = match num_type {
            NumType::Real => (random_val(), 0.0),
            NumType::Imaginary => (0.0, random_val()),
            NumType::Complex => (random_val(), random_val()),
        };
        let z = Complex64::new(zr, zi);
        // dbg!(order, &z);
        // println!("#[case({}, {}, {})]", order, z.re, z.im);
        check_against_fortran(order, z, scaling, n);
    }
}

fn check_against_fortran(order: f64, z: Complex64, scaling: Scaling, n: usize) {
    let actual = zbesj(z, order, scaling, n);
    if let Err(ref err) = actual {
        if *err == BesselError::NotYetImplemented {
            return;
        }
    }

    let (cy, nz, ierr) = zbesj_fortran(order, z, scaling, n);

    if let Err(err) = actual {
        assert_eq!(ierr, err.error_code());
    } else {
        let actual = actual.unwrap();
        assert_eq!(nz, actual.1.try_into().unwrap());

        for (czi, zi) in cy.into_iter().zip(actual.0) {
            assert_relative_eq!(czi, zi, epsilon = 1e-10, max_relative = 1e-10);
            // assert_relative_eq!(im, zi.im, epsilon = 1e-10, max_relative = 1e-10);
        }
    }
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
