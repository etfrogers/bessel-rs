use complex_bessel_rs::bessel_j::bessel_j as bessel_j_ref;
use num::complex::Complex64;
use num::pow::Pow;
use rand::seq::IndexedRandom;
use rand::{Rng, SeedableRng, rngs::SmallRng};
use rstest::{fixture, rstest};
use std::f64::consts::{FRAC_PI_2, PI};

use approx::assert_relative_eq;

use super::{check_against_fortran, check_complex_arrays_equal, zbesj_fortran_loop};
use crate::{BesselError, Scaling, bessel_j};

const RANDOM_LIMIT: f64 = 10_000.0;

enum NumType {
    Real,
    Imaginary,
    Complex,
}

fn random_val_rng(rng: &mut SmallRng) -> f64 {
    rng.random_range(-RANDOM_LIMIT..RANDOM_LIMIT)
}

#[rstest]
fn test_bessel_j_random(mut rng: SmallRng) {
    for _ in 0..1000000 {
        let order = rng.random_range(f64::EPSILON..RANDOM_LIMIT);
        let zr = random_val_rng(&mut rng);
        let zi = random_val_rng(&mut rng);
        let z = Complex64::new(zr, zi);
        // dbg!(order, &z);
        // println!("#[case({}, {}, {})]", order, z.re, z.im);
        let ans = bessel_j(order, z);
        let expected = bessel_j_ref(order, z);
        if let Ok(actual) = ans {
            let (cy_loop_fort, _, _) = zbesj_fortran_loop(order, z, Scaling::Unscaled, 1);
            check_complex_arrays_equal(&actual, &expected.unwrap(), &cy_loop_fort);
        } else {
            let err = ans.unwrap_err();
            if err == BesselError::NotYetImplemented {
                // continue;
                panic!("NotYetImplemented should not occur")
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
fn test_bessel_j_large_n_random(
    #[values(Scaling::Unscaled, Scaling::Scaled)] scaling: Scaling,
    #[values(NumType::Real, NumType::Imaginary, NumType::Complex)] num_type: NumType,
    mut rng: SmallRng,
) {
    let n = 9;
    for _ in 0..100000 {
        let order = rng.random_range(f64::EPSILON..RANDOM_LIMIT);
        let (zr, zi) = match num_type {
            NumType::Real => (random_val_rng(&mut rng), 0.0),
            NumType::Imaginary => (0.0, random_val_rng(&mut rng)),
            NumType::Complex => (random_val_rng(&mut rng), random_val_rng(&mut rng)),
        };
        let z = Complex64::new(zr, zi);
        // dbg!(order, &z);
        // println!("#[case({}, {}, {})]", order, z.re, z.im);
        check_against_fortran(order, z, scaling, n);
    }
}

#[fixture]
fn rng() -> SmallRng {
    SmallRng::seed_from_u64(42)
}

#[rstest]
fn test_bessel_j_random_logspace(
    #[values(Scaling::Unscaled, Scaling::Scaled)] scaling: Scaling,
    #[values(NumType::Real, NumType::Imaginary, NumType::Complex)] num_type: NumType,
    mut rng: SmallRng,
) {
    let n = 9;
    let mut random_val = |pos: bool| -> f64 {
        let ll = f64::MIN_POSITIVE;
        let ul = f64::MAX;
        let r = rng.random_range(ll.log10()..ul.log10());
        let sign = if pos {
            1.0
        } else {
            *([-1.0, 1.0].choose(&mut rng).unwrap())
        };
        sign * 10.0.pow(r)
    };
    for _ in 0..500000 {
        let order = random_val(true);
        let (zr, zi) = match num_type {
            NumType::Real => (random_val(false), 0.0),
            NumType::Imaginary => (0.0, random_val(false)),
            NumType::Complex => (random_val(false), random_val(false)),
        };
        let z = Complex64::new(zr, zi);
        check_against_fortran(order, z, scaling, n);
    }
}

#[rstest]
fn test_fortran_ang(mut rng: SmallRng) {
    let mut random_val = || random_val_rng(&mut rng);
    const THREE_PI_BY_2: f64 = 4.71238898038468986e+00;

    let fortran_ang = |zth: Complex64| -> f64 {
        let mut ang = THREE_PI_BY_2;
        if !(zth.re >= 0.0 && zth.im < 0.0) {
            ang = FRAC_PI_2;
            if zth.re != 0.0 {
                // ang = zth.arg();
                ang = (zth.im / zth.re).atan();
            }
            if zth.re < 0.0 {
                ang += PI;
            }
        }
        ang
    };

    let shift_arg = |zth: Complex64| -> f64 {
        let mut ang = zth.arg();
        if ang < 0.0 {
            ang = (PI * 2.0) + ang;
        }
        ang.clamp(0.0, THREE_PI_BY_2)
    };

    for _ in 0..1000000 {
        let z = Complex64::new(random_val(), random_val());
        let z_re = Complex64::new(z.re, 0.0);
        let z_im = Complex64::new(0.0, z.im);

        assert_relative_eq!(fortran_ang(z), shift_arg(z));
        assert_relative_eq!(fortran_ang(z_re), shift_arg(z_re));
        assert_relative_eq!(fortran_ang(z_im), shift_arg(z_im));
    }
}
