use crate::{
    BesselError, BesselFloat, HankelKind, Scaling, airy, airy_b, airy_bp, airyp, bessel_i,
    bessel_j, bessel_k, bessel_y,
    derivatives::{
        bessel_i_derivative, bessel_j_derivative, bessel_k_derivative, bessel_y_derivative,
        hankel1_derivative, hankel2_derivative,
    },
    hankel,
};
use approx::assert_relative_eq;
use num::Complex;
use rstest::rstest;
use std::f64::consts::PI;

const ORDER_SMOOTHNESS_TOLERANCE: f64 = 1e-3; // Order derivative difference is O(delta) - we don't have analytic order derivative
const Z_SMOOTHNESS_TOLERANCE: f64 = 1e-7; // z Taylor-shooting tolerance
const DELTA: f64 = 1e-4;

/// Amos parameter RL: lower boundary of asymptotic expansion
/// for large z (≈ 21.784_271_729_432_426).
fn asymptotic_z_limit() -> f64 {
    f64::MACHINE_CONSTANTS.asymptotic_z_limit
}

/// Amos parameter FNUL: lower boundary of uniform asymptotic series
/// for large order (≈ 85.921_358_647_162_12).
fn asymptotic_order_limit() -> f64 {
    f64::MACHINE_CONSTANTS.asymptotic_order_limit
}

/// Amos parameter R2: boundary separating heuristic Miller truncation
/// from forward recurrence loop (≈ 28.666_666_666_666_668).
fn recurrence_threshold() -> f64 {
    let bits = (f64::MANTISSA_DIGITS - 1).clamp(12, 60) as f64;
    (2.0 / 3.0) * bits - 6.0
}

fn hankel1(order: f64, z: Complex<f64>) -> Result<Complex<f64>, BesselError> {
    hankel(order, z, HankelKind::First)
}

fn hankel2(order: f64, z: Complex<f64>) -> Result<Complex<f64>, BesselError> {
    hankel(order, z, HankelKind::Second)
}

// -----------------------------------------------------------------------------
// Generic Smoothness Assertion Helpers
// -----------------------------------------------------------------------------

/// Verifies that a function is C¹ smooth across a radial boundary circle |z| = r_boundary
/// at angle theta. Steps radially by ±delta along u = exp(i·theta) and compares the
/// finite-difference slope across the boundary against the exact analytical derivative.
fn assert_radial_smoothness<F, DF>(func: F, deriv: DF, r_boundary: f64, context: &str)
where
    F: Fn(Complex<f64>) -> Result<Complex<f64>, BesselError>,
    DF: Fn(Complex<f64>, u32) -> Result<Complex<f64>, BesselError>,
{
    const TEST_ANGLES: &[f64] = &[
        0.0,
        PI / 6.0,
        PI / 4.0,
        PI / 2.0,
        2.0 * PI / 3.0,
        3.0 * PI / 4.0,
        -PI / 4.0,
        -PI / 2.0,
    ];

    for &theta in TEST_ANGLES {
        let direction = Complex::from_polar(1.0, theta);
        let z_center = Complex::from_polar(r_boundary, theta);
        let theta_context = format!("{context} (θ = {theta:.3} rad, r = {r_boundary})");

        assert_directional_smoothness(&func, &deriv, z_center, direction, &theta_context);
    }
}

/// Verifies that a function is C¹ smooth across a Cartesian boundary line (e.g. Re(z)=0 or Im(z)=0).
/// Steps by ±delta along `direction` and compares the directional slope against the analytical derivative.
fn assert_directional_smoothness<F, DF>(
    func: F,
    deriv: DF,
    z_center: Complex<f64>,
    direction: Complex<f64>,

    context: &str,
) where
    F: Fn(Complex<f64>) -> Result<Complex<f64>, BesselError>,
    DF: Fn(Complex<f64>, u32) -> Result<Complex<f64>, BesselError>,
{
    assert_relative_eq!(direction.norm(), 1.0);
    let dz = DELTA * direction;

    let z_minus = z_center - dz;
    let z_plus = z_center + dz;

    let f_minus = func(z_minus)
        .unwrap_or_else(|e| panic!("{context}: Left evaluation failed at z={z_minus}: {e:?}"));
    let f_plus = func(z_plus)
        .unwrap_or_else(|e| panic!("{context}: Right evaluation failed at z={z_plus}: {e:?}"));
    let f_center = func(z_center)
        .unwrap_or_else(|e| panic!("{context}: Center evaluation failed at z={z_center}: {e:?}"));

    let df_minus = deriv(z_minus, 1).unwrap_or_else(|e| {
        panic!("{context}: Analytical derivative failed at z={z_minus}: {e:?}")
    });
    let df_plus = deriv(z_plus, 1)
        .unwrap_or_else(|e| panic!("{context}: Analytical derivative failed at z={z_plus}: {e:?}"));

    let d2f_minus = deriv(z_minus, 2).unwrap_or_else(|e| {
        panic!("{context}: Analytical derivative failed at z={z_minus}: {e:?}")
    });
    let d2f_plus = deriv(z_plus, 2)
        .unwrap_or_else(|e| panic!("{context}: Analytical derivative failed at z={z_plus}: {e:?}"));

    let from_minus = f_minus + (df_minus * dz) + 0.5 * d2f_minus * dz.powi(2);
    let from_plus = f_plus - (df_plus * dz) + 0.5 * d2f_plus * dz.powi(2);
    // let df_numerical = (f_plus - f_minus) / (2.0 * delta * direction);
    let assert_eq = |a: Complex<f64>, b: Complex<f64>| {
        let scale = a.norm().max(b.norm()).max(1e-15);
        let diff = (a - b).norm();
        let rel_err = diff / scale;

        assert!(
            rel_err < Z_SMOOTHNESS_TOLERANCE,
            "{context} at z={z_center}: rel_err={rel_err:e} > tol={Z_SMOOTHNESS_TOLERANCE:e} \
         (a={a}, at_center={b})"
        );
    };

    assert_eq(from_minus, f_center);
    assert_eq(from_plus, f_center);
}

/// Verifies that a function is continuous and smooth with respect to the order parameter ν
/// across an order boundary nu_boundary. Compares left-side and right-side one-sided finite differences.
fn assert_order_smoothness<F>(func: F, nu_boundary: f64, z: Complex<f64>, context: &str)
where
    F: Fn(f64, Complex<f64>) -> Result<Complex<f64>, BesselError>,
{
    let f_minus = func(nu_boundary - DELTA, z).unwrap_or_else(|e| {
        panic!(
            "{context}: Left order eval failed at nu={}: {e:?}",
            nu_boundary - DELTA
        )
    });
    let f_center = func(nu_boundary, z).unwrap_or_else(|e| {
        panic!("{context}: Center order eval failed at nu={nu_boundary}: {e:?}")
    });
    let f_plus = func(nu_boundary + DELTA, z).unwrap_or_else(|e| {
        panic!(
            "{context}: Right order eval failed at nu={}: {e:?}",
            nu_boundary + DELTA
        )
    });

    let slope_left = (f_center - f_minus) / DELTA;
    let slope_right = (f_plus - f_center) / DELTA;

    let diff = (slope_right - slope_left).norm();
    let scale = slope_left
        .norm()
        .max(slope_right.norm())
        .max(f_center.norm())
        .max(1e-15);
    let rel_err = diff / scale;

    assert!(
        rel_err < ORDER_SMOOTHNESS_TOLERANCE,
        "{context} at nu={nu_boundary}, z={z}: rel_err={rel_err:e} > tol={ORDER_SMOOTHNESS_TOLERANCE:e} \
         (slope_L={slope_left}, slope_R={slope_right})"
    );
}

// -----------------------------------------------------------------------------
// 1. Radial Threshold Tests (|z| boundaries)
// -----------------------------------------------------------------------------

/// Amos Airy Diagram: Boundary |z| = 1.0 separating Maclaurin power series (|z| <= 1.0)
/// from Bessel K_{1/3}, K_{2/3} formulation (|z| > 1.0).
#[rstest]
fn test_airy_radial_smoothness_at_z_1() {
    // Ai: derivative is Ai'
    assert_radial_smoothness(
        airy::<f64, _>,
        |z, order| match order {
            1 => airyp(z),
            2 => airy(z).map(|ai| z * ai),
            _ => unreachable!(),
        },
        1.0,
        "Airy Ai at |z|=1.0",
    );

    // Ai': derivative is Ai'' = z · Ai(z) (from Airy ODE w'' - z w = 0)
    assert_radial_smoothness(
        airyp::<f64, _>,
        |z, order| match order {
            1 => airy(z).map(|ai| z * ai),
            2 => Ok(airy(z)? + z * airyp(z)?),
            _ => unreachable!(),
        },
        1.0,
        "Airy Ai' at |z|=1.0",
    );

    // Bi: derivative is Bi'
    assert_radial_smoothness(
        airy_b::<f64, _>,
        |z, order| match order {
            1 => airy_bp(z),
            2 => airy_b(z).map(|ai| z * ai),
            _ => unreachable!(),
        },
        1.0,
        "Airy Bi at |z|=1.0",
    );

    // Bi': derivative is Bi'' = z · Bi(z)
    assert_radial_smoothness(
        airy_bp::<f64, _>,
        |z, order| match order {
            1 => airy_b(z).map(|ai| z * ai),
            2 => Ok(airy_b(z)? + z * airy_bp(z)?),
            _ => unreachable!(),
        },
        1.0,
        "Airy Bi' at |z|=1.0",
    );
}

/// Amos Figure 2 (K_ν): Boundary |z| = 2.0 in the right-half plane separating
/// Temme's power series expansion (|z| <= 2.0) from Miller's backward recurrence (|z| > 2.0).
#[rstest]
#[case(0.1)]
#[case(0.4)]
#[case(0.7)]
#[case(1.2)]
#[case(2.5)]
#[case(85.0)] // just less than asymptotic_order_limit
fn test_bessel_k_radial_smoothness_at_z_2(#[case] nu: f64) {
    assert_radial_smoothness(
        |z| bessel_k(nu, z),
        |z, order| bessel_k_derivative(nu, z, order, Scaling::Unscaled),
        2.0,
        &format!("Bessel K_{nu} at |z|=2.0"),
    );
}

/// Amos Figure 1 (I_ν): Parabolic boundary |z| = 2√(ν + 1) (which equals |z| = 2.0 at ν = 0)
/// separating Taylor power series from Miller's backward recurrence.
#[rstest]
#[case(0.0, 2.0)]
#[case(0.5, 2.0 * 1.5_f64.sqrt())]
#[case(1.0, 2.0 * 2.0_f64.sqrt())]
#[case(3.0, 4.0)]
#[case(8.0, 6.0)]
fn test_bessel_i_radial_smoothness_power_series_boundary(#[case] nu: f64, #[case] r_boundary: f64) {
    assert_radial_smoothness(
        |z| bessel_i(nu, z),
        |z, order| bessel_i_derivative(nu, z, order, Scaling::Unscaled),
        r_boundary,
        &format!("Bessel I_{nu} at parabolic boundary |z|={r_boundary:.3}"),
    );
}

/// Amos Figure 1 (I_ν): Parabolic boundary |z| = ν² / 2 separating Hankel large-argument
/// asymptotics (Domain II) from Miller with Wronskian (Domain V, for ν <= √(2·FNUL))
/// or uniform Debye asymptotics with recurrence (Domain IV, for ν > √(2·FNUL)).
#[rstest]
#[case(7.0, 7.0 * 7.0 / 2.0)]
#[case(8.0, 8.0 * 8.0 / 2.0)]
#[case(10.0, 10.0 * 10.0 / 2.0)]
#[case(12.0, 12.0 * 12.0 / 2.0)]
#[case(15.0, 15.0 * 15.0 / 2.0)]
fn test_bessel_i_radial_smoothness_parabolic_asymptotic_boundary(
    #[case] nu: f64,
    #[case] r_boundary: f64,
) {
    assert_radial_smoothness(
        |z| bessel_i(nu, z),
        |z, order| bessel_i_derivative(nu, z, order, Scaling::Unscaled),
        r_boundary,
        &format!("Bessel I_{nu} at parabolic asymptotic boundary |z|={r_boundary:.3}"),
    );
}

/// Amos Figure 1 (I_ν): Boundary |z| = FNUL ≈ 85.921 separating Miller's algorithm
/// with Wronskian normalization (Domain V) from uniform Debye asymptotics with
/// backward recurrence (Domain IV), active for √(2·FNUL) < ν <= FNUL (≈ 13.11 < ν <= 85.921).
#[rstest]
#[case(15.0)]
#[case(20.0)]
#[case(45.0)]
#[case(80.0)]
#[case(85.0)]
fn test_bessel_i_radial_smoothness_at_asymptotic_order_limit(#[case] nu: f64) {
    assert_radial_smoothness(
        |z| bessel_i(nu, z),
        |z, order| bessel_i_derivative(nu, z, order, Scaling::Unscaled),
        asymptotic_order_limit(),
        &format!("Bessel I_{nu} at |z|=FNUL"),
    );
}

/// Amos Figure 1 (I_ν): Boundary |z| = RL ≈ 21.784 separating Miller's algorithm
/// from Hankel's large-argument asymptotic series for small orders (ν² <= 2|z|)
/// and recurrence + Wronskian for larger orders.
#[rstest]
#[case(0.0)]
#[case(0.5)]
#[case(1.0)]
#[case(3.0)]
#[case(10.0)]
#[case(45.23)]
#[case(85.0)]
fn test_bessel_i_radial_smoothness_at_asymptotic_z_limit(#[case] nu: f64) {
    assert_radial_smoothness(
        |z| bessel_i(nu, z),
        |z, order| bessel_i_derivative(nu, z, order, Scaling::Unscaled),
        asymptotic_z_limit(),
        &format!("Bessel I_{nu} at |z|=RL"),
    );
}

/// Amos internal heuristic: Boundary |z| = R2 ≈ 28.666 in `determine_miller_starting_k`
/// separating heuristic curve factor from forward recurrence loop.
#[rstest]
#[case(0.3)]
#[case(1.5)]
#[case(4.0)]
fn test_bessel_k_radial_smoothness_at_recurrence_threshold(#[case] nu: f64) {
    assert_radial_smoothness(
        |z| bessel_k(nu, z),
        |z, order| bessel_k_derivative(nu, z, order, Scaling::Unscaled),
        recurrence_threshold(),
        &format!("Bessel K_{nu} at Miller recurrence threshold |z|=R2"),
    );
}

/// Amos uniform asymptotics angular boundary |arg(z)| = π/3 (60°):
/// Separates direct uniform Airy expansions (ZUNI1 / ZUNK1) for |arg(z)| <= π/3
/// from rotated uniform expansions (ZUNI2 / ZUNK2) for π/3 < |arg(z)| <= π/2.
#[rstest]
#[case(100.0, 30.0)]
#[case(100.0, 80.0)]
#[case(120.0, 50.0)]
fn test_uniform_asymptotics_angular_boundary_pi_over_3(#[case] nu: f64, #[case] r: f64) {
    let theta = PI / 3.0;
    let z_center = Complex::from_polar(r, theta);
    // Normal direction pointing across the ray in the direction of increasing argument theta:
    // d/dθ (r e^{iθ}) / r = i e^{iθ} = exp(i(θ + π/2))
    let direction = Complex::from_polar(1.0, theta + PI / 2.0);

    assert_directional_smoothness(
        |z| bessel_k(nu, z),
        |z, order| bessel_k_derivative(nu, z, order, Scaling::Unscaled),
        z_center,
        direction,
        &format!("Bessel K_{nu} across uniform asymptotics angular boundary θ=π/3 (r={r})"),
    );

    assert_directional_smoothness(
        |z| bessel_i(nu, z),
        |z, order| bessel_i_derivative(nu, z, order, Scaling::Unscaled),
        z_center,
        direction,
        &format!("Bessel I_{nu} across uniform asymptotics angular boundary θ=π/3 (r={r})"),
    );
}

// -----------------------------------------------------------------------------
// 2. Order Threshold Tests (ν boundaries)
// -----------------------------------------------------------------------------

/// Amos Figure 1 (I_ν): Boundary ν = 1.0 separating Miller with Neumann series normalization (ν <= 1.0)
/// from Miller with Wronskian normalization (ν > 1.0).
#[rstest]
#[case(Complex::new(3.0, 0.0))]
#[case(Complex::new(5.0, 2.0))]
#[case(Complex::new(10.0, -3.0))]
#[case(Complex::new(15.0, 5.0))]
fn test_bessel_i_order_smoothness_at_nu_1(#[case] z: Complex<f64>) {
    assert_order_smoothness(
        bessel_i,
        1.0,
        z,
        "Bessel I_ν at ν=1.0 (Neumann vs Wronskian)",
    );
}

/// Amos Figure 1 & 2: Boundary ν = FNUL ≈ 85.921 separating recurrence / series from
/// Debye uniform asymptotic expansions for large orders.
#[rstest]
#[case(1.5)] // power-series -> assymptotic expansion for k
#[case(15.0)]
#[case(30.0)]
#[case(60.0)]
#[case(90.0)] // asymptotic -> backward recursion regime for i
fn test_bessel_order_smoothness_at_asymptotic_order_limit(#[case] zr: f64) {
    let z = Complex::new(zr, 0.0);

    assert_order_smoothness(
        bessel_k,
        asymptotic_order_limit(),
        z,
        "Bessel K_ν at ν=FNUL (Recurrence vs Uniform Asymptotics)",
    );

    assert_order_smoothness(
        bessel_i,
        asymptotic_order_limit(),
        z,
        "Bessel I_ν at ν=FNUL (Recurrence vs Uniform Asymptotics)",
    );

    assert_order_smoothness(
        hankel1,
        asymptotic_order_limit(),
        z,
        "Hankel H^(1)_ν at ν=FNUL (Recurrence vs Uniform Asymptotics)",
    );
}

/// K_ν order boundary transitions across low orders:
/// - ν = 0.0: Negative-order reflection K_{-ν} = K_ν and ln(2/z) limit in Temme series
/// - ν = 0.5: Exact closed form K_{1/2}(z) = √(π/2z) e^{-z} (0 recurrence steps)
/// - ν = 1.0: Integer order limit and fractional part sign-flip
/// - ν = 1.5: Exact closed form + boundary where forward recurrence begins (1 recurrence step)
#[rstest]
#[case(0.0)]
#[case(0.1)] // Temme series Taylor-vs-Gamma threshold |nu - round(nu)| = 0.1
#[case(0.5)]
#[case(1.0)]
#[case(1.5)]
fn test_bessel_k_order_smoothness_low_orders(#[case] nu: f64) {
    // 1. Temme power series regime (|z| <= 2.0)
    let z_small = Complex::new(1.0, 0.5);
    assert_order_smoothness(
        bessel_k,
        nu,
        z_small,
        &format!("Bessel K_ν at ν={nu} in Temme series regime (|z| <= 2)"),
    );

    // 2. Miller backward recurrence regime (|z| > 2.0)
    let z_large = Complex::new(4.0, 1.0);
    assert_order_smoothness(
        bessel_k,
        nu,
        z_large,
        &format!("Bessel K_ν at ν={nu} in Miller regime (|z| > 2)"),
    );
}

// -----------------------------------------------------------------------------
// 3. Half-Plane & Axis Continuation Tests
// -----------------------------------------------------------------------------

/// Imaginary Axis Re(z) = 0: Separates right-half plane direct calculation from
/// left-half plane analytic continuation for I_ν and K_ν.
#[rstest]
#[case(0.3, 2.0)]
#[case(0.3, 5.0)]
#[case(1.0, 3.0)]
#[case(2.5, 4.0)]
#[case(2.5, 10.0)]
fn test_imaginary_axis_smoothness_i_and_k(#[case] nu: f64, #[case] y: f64) {
    let z_center = Complex::new(0.0, y);
    let horizontal = Complex::new(1.0, 0.0);

    // I_ν: Re(z) >= 0 direct vs Re(z) < 0 continuation
    assert_directional_smoothness(
        |z| bessel_i(nu, z),
        |z, order| bessel_i_derivative(nu, z, order, Scaling::Unscaled),
        z_center,
        horizontal,
        &format!("Bessel I_{nu} crossing Re(z)=0 at y={y}"),
    );

    // K_ν: Re(z) >= 0 direct vs Re(z) < 0 continuation
    assert_directional_smoothness(
        |z| bessel_k(nu, z),
        |z, order| bessel_k_derivative(nu, z, order, Scaling::Unscaled),
        z_center,
        horizontal,
        &format!("Bessel K_{nu} crossing Re(z)=0 at y={y}"),
    );
}

/// Positive Real Axis Im(z) = 0 (x > 0): Separates upper-half plane calculation
/// from lower-half plane conjugate symmetry reflection for J_ν and Y_ν.
#[rstest]
#[case(0.0, 1.5)]
#[case(0.0, 6.0)]
#[case(0.7, 3.0)]
#[case(1.5, 8.0)]
#[case(2.0, 12.0)]
fn test_real_axis_smoothness_j_and_y(#[case] nu: f64, #[case] x: f64) {
    let z_center = Complex::new(x, 0.0);
    let vertical = Complex::new(0.0, 1.0); // Direction along imaginary axis

    // J_ν: Im(z) >= 0 rotation vs Im(z) < 0 conjugate reflection
    assert_directional_smoothness(
        |z| bessel_j(nu, z),
        |z, order| bessel_j_derivative(nu, z, order, Scaling::Unscaled),
        z_center,
        vertical,
        &format!("Bessel J_{nu} crossing Im(z)=0 at x={x}"),
    );

    // Y_ν: Im(z) >= 0 rotation vs Im(z) < 0 conjugate reflection
    assert_directional_smoothness(
        |z| bessel_y(nu, z),
        |z, order| bessel_y_derivative(nu, z, order, Scaling::Unscaled),
        z_center,
        vertical,
        &format!("Bessel Y_{nu} crossing Im(z)=0 at x={x}"),
    );
}

/// Positive Real Axis Im(z) = 0 (x > 0): For Hankel functions H^(1) and H^(2),
/// tests analytic continuation across the rotated cut boundary.
#[rstest]
#[case(0.0, 2.0)]
#[case(1.0, 5.0)]
#[case(2.5, 8.0)]
fn test_real_axis_smoothness_hankel(#[case] nu: f64, #[case] x: f64) {
    let z_center = Complex::new(x, 0.0);
    let vertical = Complex::new(0.0, 1.0);

    assert_directional_smoothness(
        |z| hankel1(nu, z),
        |z, order| hankel1_derivative(nu, z, order, Scaling::Unscaled),
        z_center,
        vertical,
        &format!("Hankel H^(1)_{nu} crossing Im(z)=0 at x={x}"),
    );

    assert_directional_smoothness(
        |z| hankel2(nu, z),
        |z, order| hankel2_derivative(nu, z, order, Scaling::Unscaled),
        z_center,
        vertical,
        &format!("Hankel H^(2)_{nu} crossing Im(z)=0 at x={x}"),
    );
}
