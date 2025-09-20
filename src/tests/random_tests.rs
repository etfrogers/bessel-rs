use complex_bessel_rs::{
    bessel_i::bessel_i as bessel_i_ref, bessel_j::bessel_j as bessel_j_ref,
    bessel_k::bessel_k as bessel_k_ref, bessel_y::bessel_y as bessel_y_ref,
};
use num::complex::Complex64;
use num::pow::Pow;
use rand::seq::IndexedRandom;
use rand::{Rng, SeedableRng, rngs::SmallRng};
use rstest::{fixture, rstest};
use std::f64::consts::{FRAC_PI_2, PI};

use approx::assert_relative_eq;

use super::{
    BesselFortranSig, BesselSig, TOLERANCE_MARGIN, airy_ref, bessel_h_ref, check_against_fortran,
    check_complex_arrays_equal, fortran_bess_loop, zbesh_first, zbesh_fortran_first,
    zbesh_fortran_second, zbesh_second, zbesi_fortran, zbesj_fortran, zbesk_fortran, zbesy_fortran,
};

use crate::{
    BesselError, HankelKind, Scaling, airy, airyp,
    amos::{
        MACHINE_CONSTANTS, complex_bessel_i, complex_bessel_j, complex_bessel_k, complex_bessel_y,
    },
    bessel_i, bessel_j, bessel_k, bessel_y, hankel,
};

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
            let (cy_loop_fort, _, _) =
                fortran_bess_loop(order, z, Scaling::Unscaled, 1, zbesj_fortran);
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
fn test_bessel_i_random(mut rng: SmallRng) {
    for _ in 0..1000000 {
        let order = rng.random_range(f64::EPSILON..RANDOM_LIMIT);
        let zr = random_val_rng(&mut rng);
        let zi = random_val_rng(&mut rng);
        let z = Complex64::new(zr, zi);
        // dbg!(order, &z);
        // println!("#[case({}, {}, {})]", order, z.re, z.im);
        let ans = bessel_i(order, z);
        let expected = bessel_i_ref(order, z);
        if let Ok(actual) = ans {
            let (cy_loop_fort, _, _) =
                fortran_bess_loop(order, z, Scaling::Unscaled, 1, zbesj_fortran);
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
fn test_bessel_k_random(mut rng: SmallRng) {
    for _ in 0..1000000 {
        let order = rng.random_range(f64::EPSILON..RANDOM_LIMIT);
        let zr = random_val_rng(&mut rng);
        let zi = random_val_rng(&mut rng);
        let z = Complex64::new(zr, zi);
        // dbg!(order, &z);
        // println!("#[case({}, {}, {})]", order, z.re, z.im);
        let ans = bessel_k(order, z);
        let expected = bessel_k_ref(order, z);
        if let Ok(actual) = ans {
            let (cy_loop_fort, _, _) =
                fortran_bess_loop(order, z, Scaling::Unscaled, 1, zbesk_fortran);
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
fn test_bessel_y_random(mut rng: SmallRng) {
    for _ in 0..1000000 {
        let order = rng.random_range(f64::EPSILON..RANDOM_LIMIT);
        let zr = random_val_rng(&mut rng);
        let zi = random_val_rng(&mut rng);
        let z = Complex64::new(zr, zi);
        // dbg!(order, &z);
        // println!("#[case({}, {}, {})]", order, z.re, z.im);
        let ans = bessel_y(order, z);
        let expected = bessel_y_ref(order, z);
        if let Ok(actual) = ans {
            let (cy_loop_fort, _, _) =
                fortran_bess_loop(order, z, Scaling::Unscaled, 1, zbesk_fortran);
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
fn test_bessel_h_random(
    mut rng: SmallRng,
    #[values(HankelKind::First, HankelKind::Second)] kind: HankelKind,
) {
    for _ in 0..1000000 {
        let order = rng.random_range(f64::EPSILON..RANDOM_LIMIT);
        let zr = random_val_rng(&mut rng);
        let zi = random_val_rng(&mut rng);
        let z = Complex64::new(zr, zi);
        // dbg!(order, &z);
        // println!("#[case({}, {}, {})]", order, z.re, z.im);
        let ans = hankel(order, z, kind);
        let expected = bessel_h_ref(order, z, kind);
        if let Ok(actual) = ans {
            let (cy_loop_fort, _, _) =
                fortran_bess_loop(order, z, Scaling::Unscaled, 1, zbesj_fortran);
            check_complex_arrays_equal(&actual, &expected.unwrap(), &cy_loop_fort);
        } else {
            let err = ans.unwrap_err();
            if err == BesselError::NotYetImplemented {
                continue;
                // panic!("NotYetImplemented should not occur")
            }
            // dbg!(&err, &expected);
            assert_eq!(err.error_code(), expected.unwrap_err());
        }
    }
}

#[rstest]
fn test_airy_random(mut rng: SmallRng, #[values(true, false)] is_derivative: bool) {
    for _ in 0..1000000 {
        let zr = random_val_rng(&mut rng);
        let zi = random_val_rng(&mut rng);

        let z = Complex64::new(zr, zi);
        // dbg!(&z);
        // println!("#[case({}, {}, {})]", 1, z.re, z.im);

        let actual = if is_derivative { airyp(z) } else { airy(z) };
        let expected = airy_ref(z, is_derivative);
        if let Ok(actual) = actual {
            assert_relative_eq!(
                actual,
                &expected.unwrap(),
                max_relative = MACHINE_CONSTANTS.abs_error_tolerance * TOLERANCE_MARGIN
            );
        } else {
            assert_eq!(actual.unwrap_err().error_code(), expected.unwrap_err())
        }
    }
}

#[rstest]
fn test_bessel_j_random_negative() {
    todo!("Add negative values")
}

#[rstest]
fn test_bessel_large_n_random(
    #[values(Scaling::Unscaled, Scaling::Scaled)] scaling: Scaling,
    #[values(NumType::Real, NumType::Imaginary, NumType::Complex)] num_type: NumType,
    #[values(
        (complex_bessel_j as BesselSig, zbesj_fortran as BesselFortranSig),
        (complex_bessel_i as BesselSig, zbesi_fortran as BesselFortranSig),
        (zbesh_first as BesselSig , zbesh_fortran_first as BesselFortranSig),
        (zbesh_second as BesselSig , zbesh_fortran_second as BesselFortranSig),
        (complex_bessel_k as BesselSig, zbesk_fortran as BesselFortranSig),
        (complex_bessel_y as BesselSig, zbesy_fortran as BesselFortranSig)
    )]
    (rust_fn, fortran_fn): (BesselSig, BesselFortranSig),
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
        check_against_fortran(order, z, scaling, n, rust_fn, fortran_fn);
    }
}

#[fixture]
fn rng() -> SmallRng {
    SmallRng::seed_from_u64(42)
}

#[rstest]
fn test_bessel_random_logspace(
    #[values(Scaling::Unscaled, Scaling::Scaled)] scaling: Scaling,
    #[values(NumType::Real, NumType::Imaginary, NumType::Complex)] num_type: NumType,
    mut rng: SmallRng,
    #[values(
        (complex_bessel_j as BesselSig, zbesj_fortran as BesselFortranSig),
        (complex_bessel_i as BesselSig, zbesi_fortran as BesselFortranSig),
        (zbesh_first as BesselSig , zbesh_fortran_first as BesselFortranSig),
        (zbesh_second as BesselSig , zbesh_fortran_second as BesselFortranSig),
        (complex_bessel_k as BesselSig, zbesk_fortran as BesselFortranSig),
        (complex_bessel_y as BesselSig, zbesy_fortran as BesselFortranSig)
    )]
    (rust_fn, fortran_fn): (BesselSig, BesselFortranSig),
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
        check_against_fortran(order, z, scaling, n, rust_fn, fortran_fn);
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
