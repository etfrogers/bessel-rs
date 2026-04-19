use crate::tests::{
    BesselFortranSig, BesselSig, check_against_fortran, sig_airy, sig_airy_fortran, sig_airyp,
    sig_airyp_fortran, sig_biry, sig_biry_fortran, sig_biryp, sig_biryp_fortran, zbesh_first,
    zbesh_fortran_first, zbesh_fortran_second, zbesh_second, zbesi_fortran, zbesj_fortran,
    zbesk_fortran, zbesy_fortran,
};
use crate::{Scaling, complex_bessel_i, complex_bessel_j, complex_bessel_k, complex_bessel_y};
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
fn test_bessel_grid(
    #[values(Scaling::Unscaled, Scaling::Scaled)] scaling: Scaling,
    #[values(
        (complex_bessel_j as BesselSig, zbesj_fortran as BesselFortranSig),
        (complex_bessel_i as BesselSig, zbesi_fortran as BesselFortranSig),
        (zbesh_first as BesselSig , zbesh_fortran_first as BesselFortranSig),
        (zbesh_second as BesselSig , zbesh_fortran_second as BesselFortranSig),
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

#[rstest]
fn test_airy_grid(
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
