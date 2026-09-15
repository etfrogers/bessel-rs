use amos_bessel_rs::{
    BesselError, HankelKind, Scaling, airy, airy_b, airy_bp, airyp,
    amos::{
        complex_bessel_i, complex_bessel_j, complex_bessel_k, complex_bessel_y, complex_hankel1,
        complex_hankel2,
    },
    bessel_i, bessel_j, bessel_k, bessel_y,
    derivatives::{
        bessel_i_derivative, bessel_j_derivative, bessel_k_derivative, bessel_y_derivative,
        hankel_derivative, hankel1_derivative, hankel2_derivative,
    },
    hankel,
};
use num::Complex;
use rand::{Rng, SeedableRng, rngs::SmallRng};
use std::f64::consts::PI;
use std::panic::{AssertUnwindSafe, catch_unwind};

const RNG_SEED: u64 = 0x5EED_BEEF_C0DE_1234;

fn sample_complex_disc(rng: &mut SmallRng, min_r: f64, max_r: f64) -> Complex<f64> {
    let r = rng.random_range(min_r..max_r);
    let theta = rng.random_range(-PI..PI);
    Complex::from_polar(r, theta)
}

fn hankel1(order: f64, z: Complex<f64>) -> Result<Complex<f64>, BesselError> {
    hankel(order, z, HankelKind::First)
}

fn hankel2(order: f64, z: Complex<f64>) -> Result<Complex<f64>, BesselError> {
    hankel(order, z, HankelKind::Second)
}

// -----------------------------------------------------------------------------
// Suite 1: Three-Term Recurrences (Unscaled and Scaled)
// Cylinder: C_{ν-1}(z) + C_{ν+1}(z) = (2ν/z) C_ν(z)
// Modified I: I_{ν-1}(z) - I_{ν+1}(z) = (2ν/z) I_ν(z)
// Modified K: K_{ν+1}(z) - K_{ν-1}(z) = (2ν/z) K_ν(z)
// -----------------------------------------------------------------------------

#[test]
fn fuzz_three_term_recurrence_j_y() {
    let mut rng = SmallRng::seed_from_u64(RNG_SEED);

    for _ in 0..250 {
        let nu = rng.random_range(-25.0..25.0);
        let z = sample_complex_disc(&mut rng, 0.5, 60.0);

        // --- Unscaled J ---
        if let (Ok(j_prev), Ok(j_curr), Ok(j_next)) = (
            bessel_j(nu - 1.0, z),
            bessel_j(nu, z),
            bessel_j(nu + 1.0, z),
        ) {
            let lhs = j_prev + j_next;
            let rhs = (2.0 * nu / z) * j_curr;
            let scale = j_prev.norm() + j_next.norm() + rhs.norm();
            if scale > 1e-50 && scale.is_finite() {
                let rel_err = (lhs - rhs).norm() / scale;
                assert!(
                    rel_err < 1e-7,
                    "J recurrence failed for nu={nu}, z={z}: rel_err={rel_err:e}"
                );
            }
        }

        // --- Scaled J ---
        if let (Ok((yj_prev, _)), Ok((yj_curr, _)), Ok((yj_next, _))) = (
            complex_bessel_j(z, nu - 1.0, Scaling::Scaled, 1),
            complex_bessel_j(z, nu, Scaling::Scaled, 1),
            complex_bessel_j(z, nu + 1.0, Scaling::Scaled, 1),
        ) {
            let lhs = yj_prev[0] + yj_next[0];
            let rhs = (2.0 * nu / z) * yj_curr[0];
            let scale = yj_prev[0].norm() + yj_next[0].norm() + rhs.norm();
            if scale > 1e-50 && scale.is_finite() {
                let rel_err = (lhs - rhs).norm() / scale;
                assert!(
                    rel_err < 1e-7,
                    "Scaled J recurrence failed for nu={nu}, z={z}: rel_err={rel_err:e}"
                );
            }
        }

        // --- Unscaled Y ---
        if let (Ok(y_prev), Ok(y_curr), Ok(y_next)) = (
            bessel_y(nu - 1.0, z),
            bessel_y(nu, z),
            bessel_y(nu + 1.0, z),
        ) {
            let lhs = y_prev + y_next;
            let rhs = (2.0 * nu / z) * y_curr;
            let scale = y_prev.norm() + y_next.norm() + rhs.norm();
            if scale > 1e-50 && scale.is_finite() {
                let rel_err = (lhs - rhs).norm() / scale;
                assert!(
                    rel_err < 1e-7,
                    "Y recurrence failed for nu={nu}, z={z}: rel_err={rel_err:e}"
                );
            }
        }

        // --- Scaled Y ---
        if let (Ok((yy_prev, _)), Ok((yy_curr, _)), Ok((yy_next, _))) = (
            complex_bessel_y(z, nu - 1.0, Scaling::Scaled, 1),
            complex_bessel_y(z, nu, Scaling::Scaled, 1),
            complex_bessel_y(z, nu + 1.0, Scaling::Scaled, 1),
        ) {
            let lhs = yy_prev[0] + yy_next[0];
            let rhs = (2.0 * nu / z) * yy_curr[0];
            let scale = yy_prev[0].norm() + yy_next[0].norm() + rhs.norm();
            if scale > 1e-50 && scale.is_finite() {
                let rel_err = (lhs - rhs).norm() / scale;
                assert!(
                    rel_err < 1e-7,
                    "Scaled Y recurrence failed for nu={nu}, z={z}: rel_err={rel_err:e}"
                );
            }
        }
    }
}

#[test]
fn fuzz_three_term_recurrence_modified_i_k() {
    let mut rng = SmallRng::seed_from_u64(RNG_SEED ^ 0x1111);

    for _ in 0..250 {
        let nu = rng.random_range(-25.0..25.0);
        let z = sample_complex_disc(&mut rng, 0.5, 60.0);

        // --- Unscaled I: I_{ν-1}(z) - I_{ν+1}(z) = (2ν/z) I_ν(z) ---
        if let (Ok(i_prev), Ok(i_curr), Ok(i_next)) = (
            bessel_i(nu - 1.0, z),
            bessel_i(nu, z),
            bessel_i(nu + 1.0, z),
        ) {
            let lhs = i_prev - i_next;
            let rhs = (2.0 * nu / z) * i_curr;
            let scale = i_prev.norm() + i_next.norm() + rhs.norm();
            if scale > 1e-50 && scale.is_finite() {
                let rel_err = (lhs - rhs).norm() / scale;
                assert!(
                    rel_err < 1e-7,
                    "I recurrence failed for nu={nu}, z={z}: rel_err={rel_err:e}"
                );
            }
        }

        // --- Scaled I ---
        if let (Ok((yi_prev, _)), Ok((yi_curr, _)), Ok((yi_next, _))) = (
            complex_bessel_i(z, nu - 1.0, Scaling::Scaled, 1),
            complex_bessel_i(z, nu, Scaling::Scaled, 1),
            complex_bessel_i(z, nu + 1.0, Scaling::Scaled, 1),
        ) {
            let lhs = yi_prev[0] - yi_next[0];
            let rhs = (2.0 * nu / z) * yi_curr[0];
            let scale = yi_prev[0].norm() + yi_next[0].norm() + rhs.norm();
            if scale > 1e-50 && scale.is_finite() {
                let rel_err = (lhs - rhs).norm() / scale;
                assert!(
                    rel_err < 1e-7,
                    "Scaled I recurrence failed for nu={nu}, z={z}: rel_err={rel_err:e}"
                );
            }
        }

        // --- Unscaled K: K_{ν+1}(z) - K_{ν-1}(z) = (2ν/z) K_ν(z) ---
        if let (Ok(k_prev), Ok(k_curr), Ok(k_next)) = (
            bessel_k(nu - 1.0, z),
            bessel_k(nu, z),
            bessel_k(nu + 1.0, z),
        ) {
            let lhs = k_next - k_prev;
            let rhs = (2.0 * nu / z) * k_curr;
            let scale = k_prev.norm() + k_next.norm() + rhs.norm();
            if scale > 1e-50 && scale.is_finite() {
                let rel_err = (lhs - rhs).norm() / scale;
                assert!(
                    rel_err < 1e-7,
                    "K recurrence failed for nu={nu}, z={z}: rel_err={rel_err:e}"
                );
            }
        }

        // --- Scaled K ---
        if let (Ok((yk_prev, _)), Ok((yk_curr, _)), Ok((yk_next, _))) = (
            complex_bessel_k(z, nu - 1.0, Scaling::Scaled, 1),
            complex_bessel_k(z, nu, Scaling::Scaled, 1),
            complex_bessel_k(z, nu + 1.0, Scaling::Scaled, 1),
        ) {
            let lhs = yk_next[0] - yk_prev[0];
            let rhs = (2.0 * nu / z) * yk_curr[0];
            let scale = yk_prev[0].norm() + yk_next[0].norm() + rhs.norm();
            if scale > 1e-50 && scale.is_finite() {
                let rel_err = (lhs - rhs).norm() / scale;
                assert!(
                    rel_err < 1e-7,
                    "Scaled K recurrence failed for nu={nu}, z={z}: rel_err={rel_err:e}"
                );
            }
        }
    }
}

#[test]
fn fuzz_three_term_recurrence_hankel() {
    let mut rng = SmallRng::seed_from_u64(RNG_SEED ^ 0x2222);

    for _ in 0..200 {
        let nu = rng.random_range(-25.0..25.0);
        let z = sample_complex_disc(&mut rng, 0.5, 50.0);

        // H^(1)
        if let (Ok(h1_prev), Ok(h1_curr), Ok(h1_next)) =
            (hankel1(nu - 1.0, z), hankel1(nu, z), hankel1(nu + 1.0, z))
        {
            let lhs = h1_prev + h1_next;
            let rhs = (2.0 * nu / z) * h1_curr;
            let scale = h1_prev.norm() + h1_next.norm() + rhs.norm();
            if scale > 1e-50 && scale.is_finite() {
                let rel_err = (lhs - rhs).norm() / scale;
                assert!(
                    rel_err < 1e-7,
                    "H1 recurrence failed for nu={nu}, z={z}: rel_err={rel_err:e}"
                );
            }
        }

        // Scaled H^(1)
        if let (Ok((yh1_prev, _)), Ok((yh1_curr, _)), Ok((yh1_next, _))) = (
            complex_hankel1(z, nu - 1.0, Scaling::Scaled, 1),
            complex_hankel1(z, nu, Scaling::Scaled, 1),
            complex_hankel1(z, nu + 1.0, Scaling::Scaled, 1),
        ) {
            let lhs = yh1_prev[0] + yh1_next[0];
            let rhs = (2.0 * nu / z) * yh1_curr[0];
            let scale = yh1_prev[0].norm() + yh1_next[0].norm() + rhs.norm();
            if scale > 1e-50 && scale.is_finite() {
                let rel_err = (lhs - rhs).norm() / scale;
                assert!(
                    rel_err < 1e-7,
                    "Scaled H1 recurrence failed for nu={nu}, z={z}: rel_err={rel_err:e}"
                );
            }
        }

        // H^(2)
        if let (Ok(h2_prev), Ok(h2_curr), Ok(h2_next)) =
            (hankel2(nu - 1.0, z), hankel2(nu, z), hankel2(nu + 1.0, z))
        {
            let lhs = h2_prev + h2_next;
            let rhs = (2.0 * nu / z) * h2_curr;
            let scale = h2_prev.norm() + h2_next.norm() + rhs.norm();
            if scale > 1e-50 && scale.is_finite() {
                let rel_err = (lhs - rhs).norm() / scale;
                assert!(
                    rel_err < 1e-7,
                    "H2 recurrence failed for nu={nu}, z={z}: rel_err={rel_err:e}"
                );
            }
        }

        // Scaled H^(2)
        if let (Ok((yh2_prev, _)), Ok((yh2_curr, _)), Ok((yh2_next, _))) = (
            complex_hankel2(z, nu - 1.0, Scaling::Scaled, 1),
            complex_hankel2(z, nu, Scaling::Scaled, 1),
            complex_hankel2(z, nu + 1.0, Scaling::Scaled, 1),
        ) {
            let lhs = yh2_prev[0] + yh2_next[0];
            let rhs = (2.0 * nu / z) * yh2_curr[0];
            let scale = yh2_prev[0].norm() + yh2_next[0].norm() + rhs.norm();
            if scale > 1e-50 && scale.is_finite() {
                let rel_err = (lhs - rhs).norm() / scale;
                assert!(
                    rel_err < 1e-7,
                    "Scaled H2 recurrence failed for nu={nu}, z={z}: rel_err={rel_err:e}"
                );
            }
        }
    }
}

// -----------------------------------------------------------------------------
// Suite 2: Wronskian Invariants
// W{J_ν, Y_ν} = 2 / (π z)
// W{I_ν, K_ν} = -1 / z
// W{H^(1)_ν, H^(2)_ν} = -4i / (π z)
// W{J_ν, J_{-ν}} = -2 sin(νπ) / (π z)
// W{Ai, Bi} = 1 / π
// -----------------------------------------------------------------------------

#[test]
fn fuzz_wronskian_invariants() {
    let mut rng = SmallRng::seed_from_u64(RNG_SEED ^ 0x3333);

    for _ in 0..200 {
        let nu = rng.random_range(-20.0..20.0);
        // Restrict |Im(z)| <= 2.0 to avoid catastrophic subtractive cancellation in unscaled products
        let zr = rng.random_range(0.3..25.0);
        let zi = rng.random_range(-2.0..2.0);
        let z = Complex::new(zr, zi);

        // 1) W{J_ν, Y_ν} = 2 / (π z)
        if let (Ok(j), Ok(dj), Ok(y), Ok(dy)) = (
            bessel_j(nu, z),
            bessel_j_derivative(nu, z, 1, Scaling::Unscaled),
            bessel_y(nu, z),
            bessel_y_derivative(nu, z, 1, Scaling::Unscaled),
        ) {
            let w = j * dy - dj * y;
            let w_exact = Complex::new(2.0 / PI, 0.0) / z;
            let cancellation_scale = (j * dy).norm() + (dj * y).norm();
            if cancellation_scale > 1e-50 && cancellation_scale.is_finite() {
                let rel_err = (w - w_exact).norm() / cancellation_scale;
                assert!(
                    rel_err < 1e-7,
                    "W{{J, Y}} failed for nu={nu}, z={z}: rel_err={rel_err:e}"
                );
            }
        }

        // 2) W{I_ν, K_ν} = -1 / z (for Re(z) > 0)
        if let (Ok(i), Ok(di), Ok(k), Ok(dk)) = (
            bessel_i(nu, z),
            bessel_i_derivative(nu, z, 1, Scaling::Unscaled),
            bessel_k(nu, z),
            bessel_k_derivative(nu, z, 1, Scaling::Unscaled),
        ) {
            let w = i * dk - di * k;
            let w_exact = -Complex::new(1.0, 0.0) / z;
            let cancellation_scale = (i * dk).norm() + (di * k).norm();
            if cancellation_scale > 1e-50 && cancellation_scale.is_finite() {
                let rel_err = (w - w_exact).norm() / cancellation_scale;
                assert!(
                    rel_err < 1e-7,
                    "W{{I, K}} failed for nu={nu}, z={z}: rel_err={rel_err:e}"
                );
            }
        }

        // 3) W{H^(1)_ν, H^(2)_ν} = -4i / (π z)
        if let (Ok(h1), Ok(dh1), Ok(h2), Ok(dh2)) = (
            hankel1(nu, z),
            hankel1_derivative(nu, z, 1, Scaling::Unscaled),
            hankel2(nu, z),
            hankel2_derivative(nu, z, 1, Scaling::Unscaled),
        ) {
            let w = h1 * dh2 - dh1 * h2;
            let w_exact = Complex::new(0.0, -4.0 / PI) / z;
            let cancellation_scale = (h1 * dh2).norm() + (dh1 * h2).norm();
            if cancellation_scale > 1e-50 && cancellation_scale.is_finite() {
                let rel_err = (w - w_exact).norm() / cancellation_scale;
                assert!(
                    rel_err < 1e-7,
                    "W{{H1, H2}} failed for nu={nu}, z={z}: rel_err={rel_err:e}"
                );
            }
        }

        // 4) W{J_ν, J_{-ν}} = -2 sin(νπ) / (π z) for non-integer ν
        if (nu.fract().abs() > 0.05)
            && ((1.0 - nu.fract().abs()) > 0.05)
            && let (Ok(j_pos), Ok(dj_pos), Ok(j_neg), Ok(dj_neg)) = (
                bessel_j(nu, z),
                bessel_j_derivative(nu, z, 1, Scaling::Unscaled),
                bessel_j(-nu, z),
                bessel_j_derivative(-nu, z, 1, Scaling::Unscaled),
            )
        {
            let w = j_pos * dj_neg - dj_pos * j_neg;
            let w_exact = Complex::new(-2.0 * (nu * PI).sin() / PI, 0.0) / z;
            let cancellation_scale = (j_pos * dj_neg).norm() + (dj_pos * j_neg).norm();
            if cancellation_scale > 1e-50 && cancellation_scale.is_finite() {
                let rel_err = (w - w_exact).norm() / cancellation_scale;
                assert!(
                    rel_err < 1e-7,
                    "W{{J_ν, J_{{-ν}}}} failed for nu={nu}, z={z}: rel_err={rel_err:e}"
                );
            }
        }
    }

    // 5) Airy Wronskian: W{Ai(z), Bi(z)} = 1 / π across the complex plane
    for _ in 0..100 {
        let z = sample_complex_disc(&mut rng, 0.1, 8.0);
        if let (Ok(ai), Ok(aip), Ok(bi), Ok(bip)) = (
            airy::<f64, _>(z),
            airyp::<f64, _>(z),
            airy_b::<f64, _>(z),
            airy_bp::<f64, _>(z),
        ) {
            let w = ai * bip - aip * bi;
            let w_exact = Complex::new(1.0 / PI, 0.0);
            let cancellation_scale = (ai * bip).norm() + (aip * bi).norm();
            if cancellation_scale > 1e-50 && cancellation_scale.is_finite() {
                let rel_err = (w - w_exact).norm() / cancellation_scale;
                assert!(
                    rel_err < 1e-7,
                    "W{{Ai, Bi}} failed for z={z}: rel_err={rel_err:e}"
                );
            }
        }
    }
}

// -----------------------------------------------------------------------------
// Suite 3: Negative Order Reflection Invariants
// DLMF 10.4.1:  J_{-ν}(z) = cos(νπ) J_ν(z) - sin(νπ) Y_ν(z)
// DLMF 10.4.1:  Y_{-ν}(z) = sin(νπ) J_ν(z) + cos(νπ) Y_ν(z)
// DLMF 10.27.2: I_{-ν}(z) = I_ν(z) + (2/π) sin(νπ) K_ν(z)
// DLMF 10.27.3: K_{-ν}(z) = K_ν(z)
// DLMF 10.4.6:  H_{-ν}^{(1)}(z) = e^{iνπ} H_ν^{(1)}(z)
// DLMF 10.4.7:  H_{-ν}^{(2)}(z) = e^{-iνπ} H_ν^{(2)}(z)
// -----------------------------------------------------------------------------

#[test]
fn fuzz_negative_order_reflection_identities() {
    let mut rng = SmallRng::seed_from_u64(RNG_SEED ^ 0x4444);

    for _ in 0..200 {
        let nu = rng.random_range(0.1..20.0);
        let z = sample_complex_disc(&mut rng, 0.5, 30.0);
        let cos_pi_nu = (nu * PI).cos();
        let sin_pi_nu = (nu * PI).sin();

        // 1) J reflection
        if let (Ok(j_pos), Ok(y_pos), Ok(j_neg)) =
            (bessel_j(nu, z), bessel_y(nu, z), bessel_j(-nu, z))
        {
            let expected_j_neg = cos_pi_nu * j_pos - sin_pi_nu * y_pos;
            let scale = j_neg.norm() + expected_j_neg.norm();
            if scale > 1e-50 && scale.is_finite() {
                let rel_err = (j_neg - expected_j_neg).norm() / scale;
                assert!(
                    rel_err < 1e-7,
                    "J reflection failed for nu={nu}, z={z}: rel_err={rel_err:e}"
                );
            }
        }

        // 2) Y reflection
        if let (Ok(j_pos), Ok(y_pos), Ok(y_neg)) =
            (bessel_j(nu, z), bessel_y(nu, z), bessel_y(-nu, z))
        {
            let expected_y_neg = sin_pi_nu * j_pos + cos_pi_nu * y_pos;
            let scale = y_neg.norm() + expected_y_neg.norm();
            if scale > 1e-50 && scale.is_finite() {
                let rel_err = (y_neg - expected_y_neg).norm() / scale;
                assert!(
                    rel_err < 1e-7,
                    "Y reflection failed for nu={nu}, z={z}: rel_err={rel_err:e}"
                );
            }
        }

        // 3) I reflection
        if let (Ok(i_pos), Ok(k_pos), Ok(i_neg)) =
            (bessel_i(nu, z), bessel_k(nu, z), bessel_i(-nu, z))
        {
            let expected_i_neg = i_pos + Complex::new(2.0 * sin_pi_nu / PI, 0.0) * k_pos;
            let scale = i_neg.norm() + expected_i_neg.norm();
            if scale > 1e-50 && scale.is_finite() {
                let rel_err = (i_neg - expected_i_neg).norm() / scale;
                assert!(
                    rel_err < 1e-7,
                    "I reflection failed for nu={nu}, z={z}: rel_err={rel_err:e}"
                );
            }
        }

        // 4) K reflection: K_{-ν}(z) == K_ν(z)
        if let (Ok(k_pos), Ok(k_neg)) = (bessel_k(nu, z), bessel_k(-nu, z)) {
            let scale = k_pos.norm() + k_neg.norm();
            if scale > 1e-50 && scale.is_finite() {
                let rel_err = (k_neg - k_pos).norm() / scale;
                assert!(
                    rel_err < 1e-7,
                    "K reflection failed for nu={nu}, z={z}: rel_err={rel_err:e}"
                );
            }
        }

        // 5) H1 reflection: H_{-ν}^{(1)}(z) = e^{iνπ} H_ν^{(1)}(z)
        if let (Ok(h1_pos), Ok(h1_neg)) = (hankel1(nu, z), hankel1(-nu, z)) {
            let phase = Complex::from_polar(1.0, nu * PI);
            let expected = phase * h1_pos;
            let scale = h1_neg.norm() + expected.norm();
            if scale > 1e-50 && scale.is_finite() {
                let rel_err = (h1_neg - expected).norm() / scale;
                assert!(
                    rel_err < 1e-7,
                    "H1 reflection failed for nu={nu}, z={z}: rel_err={rel_err:e}"
                );
            }
        }

        // 6) H2 reflection: H_{-ν}^{(2)}(z) = e^{-iνπ} H_ν^{(2)}(z)
        if let (Ok(h2_pos), Ok(h2_neg)) = (hankel2(nu, z), hankel2(-nu, z)) {
            let phase = Complex::from_polar(1.0, -nu * PI);
            let expected = phase * h2_pos;
            let scale = h2_neg.norm() + expected.norm();
            if scale > 1e-50 && scale.is_finite() {
                let rel_err = (h2_neg - expected).norm() / scale;
                assert!(
                    rel_err < 1e-7,
                    "H2 reflection failed for nu={nu}, z={z}: rel_err={rel_err:e}"
                );
            }
        }
    }
}

// -----------------------------------------------------------------------------
// Suite 4: Bessel Differential Equation Invariants
// Cylinder: z^2 w''(z) + z w'(z) + (z^2 - ν^2) w(z) = 0
// Modified: z^2 w''(z) + z w'(z) - (z^2 + ν^2) w(z) = 0
// -----------------------------------------------------------------------------

#[test]
fn fuzz_bessel_differential_equation_derivatives() {
    let mut rng = SmallRng::seed_from_u64(RNG_SEED ^ 0x5555);

    for _ in 0..150 {
        let nu = rng.random_range(-20.0..20.0);
        let z = sample_complex_disc(&mut rng, 0.5, 40.0);

        // Cylinder functions: J, Y, H1, H2
        let cylinder_cases = [
            (
                "J",
                bessel_j(nu, z),
                bessel_j_derivative(nu, z, 1, Scaling::Unscaled),
                bessel_j_derivative(nu, z, 2, Scaling::Unscaled),
            ),
            (
                "Y",
                bessel_y(nu, z),
                bessel_y_derivative(nu, z, 1, Scaling::Unscaled),
                bessel_y_derivative(nu, z, 2, Scaling::Unscaled),
            ),
            (
                "H1",
                hankel1(nu, z),
                hankel1_derivative(nu, z, 1, Scaling::Unscaled),
                hankel1_derivative(nu, z, 2, Scaling::Unscaled),
            ),
            (
                "H2",
                hankel2(nu, z),
                hankel2_derivative(nu, z, 1, Scaling::Unscaled),
                hankel2_derivative(nu, z, 2, Scaling::Unscaled),
            ),
        ];

        for (name, f0_res, f1_res, f2_res) in cylinder_cases {
            if let (Ok(w), Ok(dw), Ok(d2w)) = (f0_res, f1_res, f2_res) {
                let t2 = z * z * d2w;
                let t1 = z * dw;
                let t0 = (z * z - Complex::new(nu * nu, 0.0)) * w;
                let ode = t2 + t1 + t0;
                let scale = t2.norm() + t1.norm() + t0.norm();
                if scale > 1e-50 && scale.is_finite() {
                    let rel_err = ode.norm() / scale;
                    assert!(
                        rel_err < 1e-7,
                        "Cylinder ODE for {name} failed for nu={nu}, z={z}: rel_err={rel_err:e}"
                    );
                }
            }
        }

        // Modified functions: I, K
        let modified_cases = [
            (
                "I",
                bessel_i(nu, z),
                bessel_i_derivative(nu, z, 1, Scaling::Unscaled),
                bessel_i_derivative(nu, z, 2, Scaling::Unscaled),
            ),
            (
                "K",
                bessel_k(nu, z),
                bessel_k_derivative(nu, z, 1, Scaling::Unscaled),
                bessel_k_derivative(nu, z, 2, Scaling::Unscaled),
            ),
        ];

        for (name, f0_res, f1_res, f2_res) in modified_cases {
            if let (Ok(w), Ok(dw), Ok(d2w)) = (f0_res, f1_res, f2_res) {
                let t2 = z * z * d2w;
                let t1 = z * dw;
                let t0 = -(z * z + Complex::new(nu * nu, 0.0)) * w;
                let ode = t2 + t1 + t0;
                let scale = t2.norm() + t1.norm() + t0.norm();
                if scale > 1e-50 && scale.is_finite() {
                    let rel_err = ode.norm() / scale;
                    assert!(
                        rel_err < 1e-7,
                        "Modified ODE for {name} failed for nu={nu}, z={z}: rel_err={rel_err:e}"
                    );
                }
            }
        }
    }
}

// -----------------------------------------------------------------------------
// Suite 5: Crash-Free Robustness on Extreme / Random Floats
// Verifies that arbitrary float inputs (NaN, +/-Inf, subnormals, extreme magnitudes)
// never cause unhandled panics across all library entry points.
// -----------------------------------------------------------------------------

#[test]
fn fuzz_crash_free_robustness_arbitrary_inputs() {
    let mut rng = SmallRng::seed_from_u64(RNG_SEED ^ 0x6666);

    let edge_cases = [
        0.0,
        -0.0,
        f64::MIN_POSITIVE,
        f64::MIN_POSITIVE * 0.1, // subnormal
        f64::EPSILON,
        1.0,
        -1.0,
        1e-15,
        -1e-15,
        1e-300,
        1e300,
        1e150,
        -1e150,
        f64::MAX,
        f64::MIN,
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
    ];

    // Combine explicit edge cases and random bit patterns
    let mut float_pool = Vec::from(edge_cases);
    for _ in 0..100 {
        float_pool.push(f64::from_bits(rng.random()));
    }

    for &order in &float_pool[..30] {
        for &zr in &float_pool[..30] {
            let zi = 0.0;
            let z_complex = Complex::new(zr, zi);

            // Test catch_unwind on complex entry points
            let res = catch_unwind(AssertUnwindSafe(|| {
                let _ = bessel_j(order, z_complex);
                let _ = bessel_y(order, z_complex);
                let _ = bessel_i(order, z_complex);
                let _ = bessel_k(order, z_complex);
                let _ = hankel(order, z_complex, HankelKind::First);
                let _ = hankel(order, z_complex, HankelKind::Second);
                let _ = bessel_j_derivative(order, z_complex, 1, Scaling::Unscaled);
                let _ = bessel_y_derivative(order, z_complex, 1, Scaling::Unscaled);
                let _ = bessel_i_derivative(order, z_complex, 1, Scaling::Unscaled);
                let _ = bessel_k_derivative(order, z_complex, 1, Scaling::Unscaled);
                let _ =
                    hankel_derivative(order, z_complex, HankelKind::First, 1, Scaling::Unscaled);
                let _ = airy::<f64, _>(z_complex);
                let _ = airyp::<f64, _>(z_complex);
                let _ = airy_b::<f64, _>(z_complex);
                let _ = airy_bp::<f64, _>(z_complex);
            }));

            assert!(
                res.is_ok(),
                "Panic encountered for order={order:?}, z={z_complex:?}"
            );
        }
    }
}
