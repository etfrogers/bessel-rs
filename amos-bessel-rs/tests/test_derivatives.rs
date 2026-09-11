use amos_bessel_rs::{
    BesselError, bessel_j,
    derivatives::{bessel_j_derivative, bessel_j_p},
};
use approx::assert_relative_eq;
use num::Complex;
use rstest::rstest;
use std::f64::consts::PI;

mod common;
use common::{ORDERS, Z_PARTS};

#[rstest]
#[case(0, -0.44005)]
#[case(1, 0.32515)]
#[case(2, 0.21033)]
fn test_bessel_j_hardcoded(#[case] order: u32, #[case] expected: f64) {
    let dz = bessel_j_derivative(order, 1.0, 1).unwrap();
    assert_relative_eq!(dz, expected, epsilon = 1e-4);
}

// -----------------------------------------------------------------------------
// Test 1: Bessel's Differential Equation on the grid
// z^2 · w''(z) + z · w'(z) + (z^2 - ν^2) · w(z) = 0
// -----------------------------------------------------------------------------
#[rstest]
fn test_bessel_differential_equation_grid() {
    for order in ORDERS {
        for zr in Z_PARTS {
            for zi in Z_PARTS {
                let z = Complex::new(zr, zi);
                if z == Complex::ZERO {
                    continue;
                }
                let j0 = match bessel_j(order, z) {
                    Ok(val) => val,
                    Err(_) => continue,
                };
                let j1 = match bessel_j_derivative(order, z, 1) {
                    Ok(val) => val,
                    Err(_) => continue,
                };
                let j2 = match bessel_j_derivative(order, z, 2) {
                    Ok(val) => val,
                    Err(_) => continue,
                };

                let term2 = z * z * j2;
                let term1 = z * j1;
                let term0 = (z * z - Complex::new(order * order, 0.0)) * j0;
                let ode = term2 + term1 + term0;

                let scale = term2.norm() + term1.norm() + term0.norm();
                if scale < 1e-50 || j0.norm() < 1e-50 || j2.norm() < 1e-50 {
                    continue;
                }
                let rel_residual = ode.norm() / scale;
                assert!(
                    rel_residual < 1e-7,
                    "ODE residual failed for order={}, z={}: rel_residual={:e}",
                    order,
                    z,
                    rel_residual
                );
            }
        }
    }
}

// -----------------------------------------------------------------------------
// Test 2: Standard 1st Derivative Recurrence Relations on the grid
// 1) J_ν'(z) = 1/2 · [J_{ν-1}(z) - J_{ν+1}(z)]
// 2) J_ν'(z) = J_{ν-1}(z) - (ν/z) · J_ν(z)
// 3) J_0'(z) = -J_1(z)
// 4) bessel_j_p(ν, z) == bessel_j_derivative(ν, z, 1)
// -----------------------------------------------------------------------------
#[rstest]
fn test_first_derivative_recurrence_grid() {
    for order in ORDERS {
        for zr in Z_PARTS {
            for zi in Z_PARTS {
                let z = Complex::new(zr, zi);
                if z == Complex::ZERO {
                    continue;
                }
                let dz = match bessel_j_derivative(order, z, 1) {
                    Ok(val) => val,
                    Err(_) => continue,
                };
                let dz_p = match bessel_j_p(order, z) {
                    Ok(val) => val,
                    Err(_) => continue,
                };
                assert_relative_eq!(dz, dz_p, max_relative = 1e-14);

                // Identity 1: J_ν'(z) = 0.5 * (J_{ν-1}(z) - J_{ν+1}(z))
                if let (Ok(j_prev), Ok(j_next)) =
                    (bessel_j(order - 1.0, z), bessel_j(order + 1.0, z))
                {
                    let expected_sym = (j_prev - j_next) * 0.5;
                    let scale = j_prev.norm() + j_next.norm();
                    if scale > 1e-100 {
                        let diff = (dz - expected_sym).norm() / scale;
                        assert!(
                            diff < 1e-9,
                            "Symmetric recurrence failed for order={}, z={}: diff={:e}",
                            order,
                            z,
                            diff
                        );
                    }
                }

                // Identity 2: J_ν'(z) = J_{ν-1}(z) - (ν / z) * J_ν(z)
                if let (Ok(j_prev), Ok(j_curr)) = (bessel_j(order - 1.0, z), bessel_j(order, z)) {
                    let expected_asym = j_prev - (j_curr * (order / z));
                    let scale = j_prev.norm() + (j_curr * (order / z)).norm();
                    if scale > 1e-100 {
                        let diff = (dz - expected_asym).norm() / scale;
                        assert!(
                            diff < 1e-8,
                            "Asymmetric recurrence failed for order={}, z={}: diff={:e}",
                            order,
                            z,
                            diff
                        );
                    }
                }

                // Identity 3: Special case for ν = 0 -> J_0'(z) = -J_1(z)
                if order == 0.0 {
                    if let Ok(j1) = bessel_j(1.0, z) {
                        let scale = dz.norm() + j1.norm();
                        if scale > 1e-100 {
                            let diff = (dz - (-j1)).norm() / scale;
                            assert!(
                                diff < 1e-10,
                                "J_0'(z) = -J_1(z) failed for z={}: diff={:e}",
                                z,
                                diff
                            );
                        }
                    }
                }
            }
        }
    }
}

// -----------------------------------------------------------------------------
// Test 3: 0-th Derivative Consistency on the grid
// J_ν^(0)(z) ≡ J_ν(z)
// -----------------------------------------------------------------------------
#[rstest]
fn test_zero_derivative_order_grid() {
    for order in ORDERS {
        for zr in Z_PARTS {
            for zi in Z_PARTS {
                let z = Complex::new(zr, zi);
                if let (Ok(d0), Ok(direct)) = (
                    bessel_j_derivative(order, z, 0),
                    bessel_j(order, z),
                ) {
                    assert_relative_eq!(d0, direct, max_relative = 1e-14);
                }
            }
        }
    }
}

// -----------------------------------------------------------------------------
// Test 4: Integer Reflection Parity on the grid
// (d/dz)^k J_{-n}(z) = (-1)^n · (d/dz)^k J_n(z)
// -----------------------------------------------------------------------------
#[rstest]
fn test_integer_reflection_parity_grid(#[values(1, 2, 3)] k: u32) {
    for order in ORDERS {
        if order.fract() != 0.0 {
            continue; // only integer orders
        }
        let n = order.abs() as i64;
        let parity = if n % 2 == 0 { 1.0 } else { -1.0 };
        for zr in Z_PARTS {
            for zi in Z_PARTS {
                let z = Complex::new(zr, zi);
                if let (Ok(d_pos), Ok(d_neg)) = (
                    bessel_j_derivative(n as f64, z, k),
                    bessel_j_derivative(-n as f64, z, k),
                ) {
                    let scale = d_pos.norm() + d_neg.norm();
                    if scale > 1e-100 {
                        let diff = (d_neg - d_pos * parity).norm() / scale;
                        assert!(
                            diff < 1e-10,
                            "Parity failed for n={}, k={}, z={}: diff={:e}",
                            n,
                            k,
                            z,
                            diff
                        );
                    }
                }
            }
        }
    }
}

// -----------------------------------------------------------------------------
// Test 5: Elementary Closed-Form Values for Half-Integer Orders on the grid
// J_{1/2}'(z)  = √(2 / (π·z)) · [cos(z) - sin(z) / (2z)]
// J_{-1/2}'(z) = √(2 / (π·z)) · [-sin(z) - cos(z) / (2z)]
// -----------------------------------------------------------------------------
#[rstest]
fn test_half_integer_closed_forms_grid() {
    for zr in Z_PARTS {
        for zi in Z_PARTS {
            let z = Complex::new(zr, zi);
            if z.norm() < 1e-3 {
                continue; // Avoid singularity at z = 0
            }
            let prefactor = (Complex::new(2.0 / PI, 0.0) / z).sqrt();

            // ν = 1/2
            if let Ok(computed) = bessel_j_derivative(0.5, z, 1) {
                let expected = prefactor * (z.cos() - z.sin() / (z * 2.0));
                let scale = computed.norm() + expected.norm();
                if scale > 1e-100 {
                    let diff = (computed - expected).norm() / scale;
                    assert!(
                        diff < 1e-9,
                        "Closed form J_{{1/2}}' failed for z={}: diff={:e}",
                        z,
                        diff
                    );
                }
            }

            // ν = -1/2
            if let Ok(computed) = bessel_j_derivative(-0.5, z, 1) {
                let expected = prefactor * (-z.sin() - z.cos() / (z * 2.0));
                let scale = computed.norm() + expected.norm();
                if scale > 1e-100 {
                    let diff = (computed - expected).norm() / scale;
                    assert!(
                        diff < 1e-9,
                        "Closed form J_{{-1/2}}' failed for z={}: diff={:e}",
                        z,
                        diff
                    );
                }
            }
        }
    }
}

// -----------------------------------------------------------------------------
// Test 6: High-Precision Ground Truth via Cauchy's Integral Formula on the grid
// (d/dz)^k f(z) = (k! / (M · r^k)) ∑_{m=0}^{M-1} f(z + r·e^{iθ_m}) · e^{-i k θ_m}
// -----------------------------------------------------------------------------
fn cauchy_derivative(order: f64, z: Complex<f64>, k: u32) -> Result<Complex<f64>, BesselError> {
    let m = 32;
    // Adapt radius to ensure circle does not enclose 0
    let r = (0.2 * z.norm()).clamp(0.01, 0.4);
    let mut sum = Complex::new(0.0, 0.0);

    for j in 0..m {
        let theta = 2.0 * PI * (j as f64) / (m as f64);
        let w = z + Complex::from_polar(r, theta);
        let f_w = bessel_j(order, w)?;
        let phase = Complex::from_polar(1.0, -(k as f64) * theta);
        sum += f_w * phase;
    }

    let mut fact = 1.0;
    for i in 1..=k {
        fact *= i as f64;
    }
    let prefactor = fact / ((m as f64) * r.powi(k as i32));
    Ok(sum * prefactor)
}

#[rstest]
fn test_against_cauchy_derivative_grid(
    #[values(1, 2, 3, 4, 5)] derivative_order: u32,
    #[values(0.0, 1.0, 2.0, 5.0, -2.0)] order: f64,
) {
    // Select representative grid coordinates from Z_PARTS to maintain fast test runtimes
    for zr in Z_PARTS {
        for zi in Z_PARTS {
            let z = Complex::new(zr, zi);
            if z.norm() < 0.5 {
                continue; // Avoid circle enclosing z = 0
            }

            let dz = match bessel_j_derivative(order, z, derivative_order) {
                Ok(val) => val,
                Err(_) => continue,
            };
            let dz_cauchy = match cauchy_derivative(order, z, derivative_order) {
                Ok(val) => val,
                Err(_) => continue,
            };

            let scale = dz.norm() + dz_cauchy.norm();
            if scale < 1e-100 {
                continue;
            }
            let diff = (dz - dz_cauchy).norm() / scale;
            let tol = match derivative_order {
                1..=2 => 1e-9,
                3 => 1e-8,
                _ => 1e-6,
            };
            assert!(
                diff < tol,
                "Cauchy derivative mismatch for order={}, k={}, z={}: diff={:e}",
                order,
                derivative_order,
                z,
                diff
            );
        }
    }
}
