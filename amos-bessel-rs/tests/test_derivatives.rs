use amos_bessel_rs::{
    BesselError, BesselFloat, HankelKind, bessel_i, bessel_j, bessel_k, bessel_y, hankel,
    derivatives::{
        bessel_i_derivative, bessel_i_p, bessel_j_derivative, bessel_j_p,
        bessel_k_derivative, bessel_k_p,
        bessel_y_derivative, bessel_y_p,
        hankel_derivative, hankel_p, hankel1_derivative, hankel1_p,
        hankel2_derivative, hankel2_p,
    },
};
use approx::assert_relative_eq;
use num::Complex;
use rstest::rstest;
use std::f64::consts::PI;

mod common;
use common::{ORDERS, Z_PARTS};

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum RecurrenceKind {
    Cylinder,
    ModifiedI,
    ModifiedK,
}

#[allow(type_alias_bounds)]
pub type BesselSimpleSig<T: BesselFloat = f64> =
    fn(T, Complex<T>) -> Result<Complex<T>, BesselError<T>>;

#[allow(type_alias_bounds)]
pub type DerivSig<T: BesselFloat = f64> =
    fn(T, Complex<T>, u32) -> Result<Complex<T>, BesselError<T>>;

fn hankel1(order: f64, z: Complex<f64>) -> Result<Complex<f64>, BesselError> {
    hankel(order, z, HankelKind::First)
}

fn hankel2(order: f64, z: Complex<f64>) -> Result<Complex<f64>, BesselError> {
    hankel(order, z, HankelKind::Second)
}

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
#[case(bessel_j, bessel_j_derivative, false)]
#[case(bessel_y, bessel_y_derivative, false)]
#[case(bessel_i, bessel_i_derivative, true)]
#[case(bessel_k, bessel_k_derivative, true)]
#[case(hankel1, hankel1_derivative, false)]
#[case(hankel2, hankel2_derivative, false)]
fn test_bessel_differential_equation_grid(
    #[case] func: BesselSimpleSig,
    #[case] d_func: DerivSig,
    #[case] is_modified: bool,
) {
    for order in ORDERS {
        for zr in Z_PARTS {
            for zi in Z_PARTS {
                let z = Complex::new(zr, zi);
                if z == Complex::ZERO {
                    continue;
                }
                let j0 = match func(order, z) {
                    Ok(val) => val,
                    Err(_) => continue,
                };
                let j1 = match d_func(order, z, 1) {
                    Ok(val) => val,
                    Err(_) => continue,
                };
                let j2 = match d_func(order, z, 2) {
                    Ok(val) => val,
                    Err(_) => continue,
                };

                let term2 = z * z * j2;
                let term1 = z * j1;
                let term0 = if is_modified {
                    -(z * z + Complex::new(order * order, 0.0)) * j0
                } else {
                    (z * z - Complex::new(order * order, 0.0)) * j0
                };
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
// Cylinder:
// 1) C_ν'(z) = 1/2 · [C_{ν-1}(z) - C_{ν+1}(z)]
// 2) C_ν'(z) = C_{ν-1}(z) - (ν/z) · C_ν(z)
// 3) C_0'(z) = -C_1(z)
// Modified I:
// 1) I_ν'(z) = 1/2 · [I_{ν-1}(z) + I_{ν+1}(z)]
// 2) I_ν'(z) = I_{ν-1}(z) - (ν/z) · I_ν(z)
// 3) I_0'(z) = +I_1(z)
// Modified K:
// 1) K_ν'(z) = -1/2 · [K_{ν-1}(z) + K_{ν+1}(z)]
// 2) K_ν'(z) = -K_{ν-1}(z) - (ν/z) · K_ν(z)
// 3) K_0'(z) = -K_1(z)
// All: func_p(ν, z) == d_func(ν, z, 1)
// -----------------------------------------------------------------------------
#[rstest]
#[case(bessel_j, bessel_j_derivative, bessel_j_p, RecurrenceKind::Cylinder)]
#[case(bessel_y, bessel_y_derivative, bessel_y_p, RecurrenceKind::Cylinder)]
#[case(bessel_i, bessel_i_derivative, bessel_i_p, RecurrenceKind::ModifiedI)]
#[case(bessel_k, bessel_k_derivative, bessel_k_p, RecurrenceKind::ModifiedK)]
#[case(hankel1, hankel1_derivative, hankel1_p, RecurrenceKind::Cylinder)]
#[case(hankel2, hankel2_derivative, hankel2_p, RecurrenceKind::Cylinder)]
fn test_first_derivative_recurrence_grid(
    #[case] func: BesselSimpleSig,
    #[case] d_func: DerivSig,
    #[case] p_func: BesselSimpleSig,
    #[case] kind: RecurrenceKind,
) {
    for order in ORDERS {
        for zr in Z_PARTS {
            for zi in Z_PARTS {
                let z = Complex::new(zr, zi);
                if z == Complex::ZERO {
                    continue;
                }
                let dz = match d_func(order, z, 1) {
                    Ok(val) => val,
                    Err(_) => continue,
                };
                let dz_p = match p_func(order, z) {
                    Ok(val) => val,
                    Err(_) => continue,
                };
                assert_relative_eq!(dz, dz_p, max_relative = 1e-14);

                // Identity 1: Symmetric recurrence
                if let (Ok(f_prev), Ok(f_next)) = (func(order - 1.0, z), func(order + 1.0, z)) {
                    let expected_sym = match kind {
                        RecurrenceKind::Cylinder => (f_prev - f_next) * 0.5,
                        RecurrenceKind::ModifiedI => (f_prev + f_next) * 0.5,
                        RecurrenceKind::ModifiedK => -(f_prev + f_next) * 0.5,
                    };
                    let scale = f_prev.norm() + f_next.norm();
                    if scale > 1e-100 {
                        let diff = (dz - expected_sym).norm() / scale;
                        assert!(
                            diff < 1e-8,
                            "Symmetric recurrence failed for order={}, z={}: diff={:e}",
                            order,
                            z,
                            diff
                        );
                    }
                }

                // Identity 2: Asymmetric recurrence
                if let (Ok(f_prev), Ok(f_curr)) = (func(order - 1.0, z), func(order, z)) {
                    let expected_asym = match kind {
                        RecurrenceKind::Cylinder | RecurrenceKind::ModifiedI => {
                            f_prev - (f_curr * (order / z))
                        }
                        RecurrenceKind::ModifiedK => -f_prev - (f_curr * (order / z)),
                    };
                    let scale = f_prev.norm() + (f_curr * (order / z)).norm();
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

                // Identity 3: Special case for ν = 0
                if order == 0.0
                    && let Ok(f1) = func(1.0, z)
                {
                    let expected_0 = match kind {
                        RecurrenceKind::Cylinder | RecurrenceKind::ModifiedK => -f1,
                        RecurrenceKind::ModifiedI => f1,
                    };
                    let scale = dz.norm() + f1.norm();
                    if scale > 1e-100 {
                        let diff = (dz - expected_0).norm() / scale;
                        assert!(
                            diff < 1e-10,
                            "Order 0 recurrence failed for z={}: diff={:e}",
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
// Test 3: 0-th Derivative Consistency on the grid
// f_ν^(0)(z) ≡ f_ν(z)
// -----------------------------------------------------------------------------
#[rstest]
#[case(bessel_j, bessel_j_derivative)]
#[case(bessel_y, bessel_y_derivative)]
#[case(bessel_i, bessel_i_derivative)]
#[case(bessel_k, bessel_k_derivative)]
#[case(hankel1, hankel1_derivative)]
#[case(hankel2, hankel2_derivative)]
fn test_zero_derivative_order_grid(
    #[case] func: BesselSimpleSig,
    #[case] d_func: DerivSig,
) {
    for order in ORDERS {
        for zr in Z_PARTS {
            for zi in Z_PARTS {
                let z = Complex::new(zr, zi);
                if let (Ok(d0), Ok(direct)) = (d_func(order, z, 0), func(order, z)) {
                    assert_relative_eq!(d0, direct, max_relative = 1e-14);
                }
            }
        }
    }
}

// -----------------------------------------------------------------------------
// Test 4: Integer Reflection Parity on the grid
// (d/dz)^k C_{-n}(z) = (-1)^n · (d/dz)^k C_n(z)
// (d/dz)^k I_{-n}(z) = (+1)   · (d/dz)^k I_n(z)
// (d/dz)^k K_{-n}(z) = (+1)   · (d/dz)^k K_n(z)
// -----------------------------------------------------------------------------
#[rstest]
#[case(bessel_j_derivative, false)]
#[case(bessel_y_derivative, false)]
#[case(bessel_i_derivative, true)]
#[case(bessel_k_derivative, true)]
#[case(hankel1_derivative, false)]
#[case(hankel2_derivative, false)]
fn test_integer_reflection_parity_grid(
    #[case] d_func: DerivSig,
    #[case] is_modified: bool,
    #[values(1, 2, 3)] k: u32,
) {
    for order in ORDERS {
        if order.fract() != 0.0 {
            continue; // only integer orders
        }
        let n = order.abs() as i64;
        let parity = if is_modified || n % 2 == 0 { 1.0 } else { -1.0 };
        for zr in Z_PARTS {
            for zi in Z_PARTS {
                let z = Complex::new(zr, zi);
                if let (Ok(d_pos), Ok(d_neg)) = (
                    d_func(n as f64, z, k),
                    d_func(-n as f64, z, k),
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

// -----------------------------------------------------------------------------
// Test 7: Hankel Derivative Relations on the grid
// 1) H_ν^(1)'(z) = J_ν'(z) + i · Y_ν'(z)
// 2) H_ν^(2)'(z) = J_ν'(z) - i · Y_ν'(z)
// 3) hankel_derivative(order, z, kind, k) matches hankel1/2_derivative
// -----------------------------------------------------------------------------
#[rstest]
fn test_hankel_derivative_relations_grid(
    #[values(1, 2)] k: u32,
    #[values(0.0, 1.0, 2.5, -0.5, -2.0)] order: f64,
) {
    let i_unit = Complex::new(0.0, 1.0);
    for zr in Z_PARTS {
        for zi in Z_PARTS {
            let z = Complex::new(zr, zi);
            if z.norm() < 1e-4 {
                continue;
            }

            // Test dispatch: hankel_p and hankel_derivative
            if let Ok(h1_d) = hankel1_derivative(order, z, k) {
                let h_d = hankel_derivative(order, z, HankelKind::First, k).unwrap();
                assert_relative_eq!(h1_d, h_d, max_relative = 1e-14);
                if k == 1 {
                    let h1_p_val = hankel1_p(order, z).unwrap();
                    let h_p_val = hankel_p(order, z, HankelKind::First).unwrap();
                    assert_relative_eq!(h1_d, h1_p_val, max_relative = 1e-14);
                    assert_relative_eq!(h1_d, h_p_val, max_relative = 1e-14);
                }
            }

            if let Ok(h2_d) = hankel2_derivative(order, z, k) {
                let h_d = hankel_derivative(order, z, HankelKind::Second, k).unwrap();
                assert_relative_eq!(h2_d, h_d, max_relative = 1e-14);
                if k == 1 {
                    let h2_p_val = hankel2_p(order, z).unwrap();
                    let h_p_val = hankel_p(order, z, HankelKind::Second).unwrap();
                    assert_relative_eq!(h2_d, h2_p_val, max_relative = 1e-14);
                    assert_relative_eq!(h2_d, h_p_val, max_relative = 1e-14);
                }
            }

            // Test linear relationship with J' and Y':
            // H^(1)(k) = J(k) + i * Y(k)
            // H^(2)(k) = J(k) - i * Y(k)
            // Note: For large |Im(z)|, H is exponentially decaying while J and Y are exponentially growing,
            // so evaluating J +/- i*Y in f64 loses all precision to subtractive cancellation.
            // Therefore, compare linear relation for moderate |Im(z)| <= 2.0.
            if z.im.abs() <= 2.0
                && let (Ok(dj), Ok(dy), Ok(dh1), Ok(dh2)) = (
                    bessel_j_derivative(order, z, k),
                    bessel_y_derivative(order, z, k),
                    hankel1_derivative(order, z, k),
                    hankel2_derivative(order, z, k),
                )
            {
                let expected_h1 = dj + i_unit * dy;
                let abs_diff1 = (dh1 - expected_h1).norm();
                if abs_diff1 > 1e-13 {
                    let scale1 = dh1.norm() + expected_h1.norm();
                    let diff1 = abs_diff1 / scale1;
                    assert!(
                        diff1 < 1e-9,
                        "H1 derivative linear relation failed for order={}, k={}, z={}: diff={:e}",
                        order,
                        k,
                        z,
                        diff1
                    );
                }

                let expected_h2 = dj - i_unit * dy;
                let abs_diff2 = (dh2 - expected_h2).norm();
                if abs_diff2 > 1e-13 {
                    let scale2 = dh2.norm() + expected_h2.norm();
                    let diff2 = abs_diff2 / scale2;
                    assert!(
                        diff2 < 1e-9,
                        "H2 derivative linear relation failed for order={}, k={}, z={}: diff={:e}",
                        order,
                        k,
                        z,
                        diff2
                    );
                }
            }
        }
    }
}

// -----------------------------------------------------------------------------
// Test 8: Wronskian Invariant Relations on the grid
// 1) W{J_ν, Y_ν}(z) = J_ν(z) · Y_ν'(z) - J_ν'(z) · Y_ν(z) = 2 / (π · z)
// 2) W{H^(1)_ν, H^(2)_ν}(z) = H^(1)_ν(z) · H^(2)_ν'(z) - H^(1)_ν'(z) · H^(2)_ν(z) = -4i / (π · z)
// -----------------------------------------------------------------------------
#[rstest]
fn test_wronskian_j_y_grid() {
    for order in ORDERS {
        for zr in Z_PARTS {
            for zi in Z_PARTS {
                let z = Complex::new(zr, zi);
                if z == Complex::ZERO {
                    continue;
                }
                let (j, dj, y, dy) = match (
                    bessel_j(order, z),
                    bessel_j_derivative(order, z, 1),
                    bessel_y(order, z),
                    bessel_y_derivative(order, z, 1),
                ) {
                    (Ok(j), Ok(dj), Ok(y), Ok(dy)) => (j, dj, y, dy),
                    _ => continue,
                };

                // Avoid underflow regime where J or Y cannot be represented with full significance
                if j.norm() < 1e-50 || dj.norm() < 1e-50 || y.norm() < 1e-50 || dy.norm() < 1e-50 {
                    continue;
                }

                let w = j * dy - dj * y;
                let w_exact = Complex::new(2.0 / PI, 0.0) / z;

                let term_scale = (j * dy).norm() + (dj * y).norm();
                if !term_scale.is_finite() || !w.re.is_finite() || !w.im.is_finite() {
                    continue;
                }
                if term_scale < 1e-100 {
                    continue;
                }
                let diff = (w - w_exact).norm() / term_scale;
                if !diff.is_finite() {
                    continue;
                }
                assert!(
                    diff < 1e-6,
                    "Wronskian W{{J, Y}} failed for order={}, z={}: diff={:e}",
                    order,
                    z,
                    diff
                );
            }
        }
    }
}

#[rstest]
fn test_wronskian_hankel_grid() {
    let i_unit = Complex::new(0.0, 1.0);
    for order in ORDERS {
        for zr in Z_PARTS {
            for zi in Z_PARTS {
                let z = Complex::new(zr, zi);
                if z == Complex::ZERO {
                    continue;
                }
                let (h1, dh1, h2, dh2) = match (
                    hankel1(order, z),
                    hankel1_derivative(order, z, 1),
                    hankel2(order, z),
                    hankel2_derivative(order, z, 1),
                ) {
                    (Ok(h1), Ok(dh1), Ok(h2), Ok(dh2)) => (h1, dh1, h2, dh2),
                    _ => continue,
                };

                if h1.norm() < 1e-50 || dh1.norm() < 1e-50 || h2.norm() < 1e-50 || dh2.norm() < 1e-50 {
                    continue;
                }

                let w = h1 * dh2 - dh1 * h2;
                let w_exact = -i_unit * (4.0 / PI) / z;

                let term_scale = (h1 * dh2).norm() + (dh1 * h2).norm();
                if !term_scale.is_finite() || !w.re.is_finite() || !w.im.is_finite() {
                    continue;
                }
                if term_scale < 1e-100 {
                    continue;
                }
                let diff = (w - w_exact).norm() / term_scale;
                if !diff.is_finite() {
                    continue;
                }
                assert!(
                    diff < 1e-6,
                    "Wronskian W{{H1, H2}} failed for order={}, z={}: diff={:e}",
                    order,
                    z,
                    diff
                );
            }
        }
    }
}

// -----------------------------------------------------------------------------
// Test 9: Wronskian Invariant for Modified Bessel Functions on the grid
// W{I_ν, K_ν}(z) = I_ν(z) · K_ν'(z) - I_ν'(z) · K_ν(z) = -1 / z
// -----------------------------------------------------------------------------
#[rstest]
fn test_wronskian_i_k_grid() {
    for order in ORDERS {
        for zr in Z_PARTS {
            for zi in Z_PARTS {
                let z = Complex::new(zr, zi);
                if z == Complex::ZERO {
                    continue;
                }
                let (i, di, k, dk) = match (
                    bessel_i(order, z),
                    bessel_i_derivative(order, z, 1),
                    bessel_k(order, z),
                    bessel_k_derivative(order, z, 1),
                ) {
                    (Ok(i), Ok(di), Ok(k), Ok(dk)) => (i, di, k, dk),
                    _ => continue,
                };

                if i.norm() < 1e-50 || di.norm() < 1e-50 || k.norm() < 1e-50 || dk.norm() < 1e-50 {
                    continue;
                }

                let w = i * dk - di * k;
                let w_exact = -Complex::new(1.0, 0.0) / z;

                let term_scale = (i * dk).norm() + (di * k).norm();
                if !term_scale.is_finite() || !w.re.is_finite() || !w.im.is_finite() {
                    continue;
                }
                if term_scale < 1e-100 {
                    continue;
                }
                let diff = (w - w_exact).norm() / term_scale;
                if !diff.is_finite() {
                    continue;
                }
                assert!(
                    diff < 1e-6,
                    "Wronskian W{{I, K}} failed for order={}, z={}: diff={:e}",
                    order,
                    z,
                    diff
                );
            }
        }
    }
}

// -----------------------------------------------------------------------------
// Test 10: Elementary Closed-Form Values for K Half-Integer Orders on the grid
// K_{1/2}'(z)  = -√(π / (2z)) · e^{-z} · [1 + 1 / (2z)]
// K_{-1/2}'(z) = K_{1/2}'(z)
// -----------------------------------------------------------------------------
#[rstest]
fn test_k_half_integer_closed_forms_grid() {
    for zr in Z_PARTS {
        for zi in Z_PARTS {
            let z = Complex::new(zr, zi);
            if z.re <= 0.0 || z.norm() < 1e-3 {
                continue; // Avoid branch cut on non-positive real half plane and singularity at z = 0
            }
            let prefactor = (Complex::new(PI / 2.0, 0.0) / z).sqrt() * (-z).exp();
            let expected_k_prime =
                -prefactor * (Complex::new(1.0, 0.0) + Complex::new(0.5, 0.0) / z);

            // ν = 1/2
            if let Ok(computed) = bessel_k_derivative(0.5, z, 1) {
                let scale = computed.norm() + expected_k_prime.norm();
                if scale > 1e-100 {
                    let diff = (computed - expected_k_prime).norm() / scale;
                    assert!(
                        diff < 1e-9,
                        "Closed form K_{{1/2}}' failed for z={}: diff={:e}",
                        z,
                        diff
                    );
                }
            }

            // ν = -1/2
            if let Ok(computed) = bessel_k_derivative(-0.5, z, 1) {
                let scale = computed.norm() + expected_k_prime.norm();
                if scale > 1e-100 {
                    let diff = (computed - expected_k_prime).norm() / scale;
                    assert!(
                        diff < 1e-9,
                        "Closed form K_{{-1/2}}' failed for z={}: diff={:e}",
                        z,
                        diff
                    );
                }
            }
        }
    }
}

#[test]
fn test_derivatives_f32() {
    let dz = bessel_j_derivative(0.0f32, 1.0f32, 1).unwrap();
    assert_relative_eq!(dz, -0.4400506f32, max_relative = 1e-5);

    let z_cpx = Complex::new(1.0f32, 0.5f32);
    let dz_cpx = bessel_j_derivative(1.0f32, z_cpx, 2).unwrap();
    assert!(dz_cpx.re.is_finite() && dz_cpx.im.is_finite());
}

#[test]
fn test_derivative_order_too_large() {
    let err = bessel_j_derivative(0.0, 1.0, 61).unwrap_err();
    match err {
        BesselError::InvalidInput { details } => {
            assert!(details.contains("too large"));
        }
        _ => panic!("Expected InvalidInput error, got {:?}", err),
    }
}

#[test]
fn test_derivative_real_input_and_complex_error() {
    let val: f64 = bessel_j_derivative(0.0, 2.0, 1).unwrap();
    assert_relative_eq!(val, -0.5767248077568734, max_relative = 1e-12);

    let err = bessel_y_derivative(0.0, -3.0, 1).unwrap_err();
    match err {
        BesselError::ComplexOutputForRealInput { output } => {
            assert!(output.norm() > 0.0);
        }
        _ => panic!("Expected ComplexOutputForRealInput, got {:?}", err),
    }
}

