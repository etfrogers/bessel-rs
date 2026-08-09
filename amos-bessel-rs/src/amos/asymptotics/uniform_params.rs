//! Parameter calculation for uniform asymptotic expansions of Bessel functions.
//!
//! This module provides geometric and asymptotic parameters for evaluating Bessel functions
//! of large orders $\nu$:
//! - **Debye Uniform Expansions** for modified Bessel functions $I_\nu(z)$ and $K_\nu(z)$
//!   ([`DebyeGeometry`], [`DebyeParams`]), originally implemented in Fortran `ZUNIK`.
//! - **Airy Uniform Expansions** for Bessel and Hankel functions $J_\nu(z), Y_\nu(z), H^{(1)}_\nu(z), H^{(2)}_\nu(z)$
//!   ([`AiryGeometry`], [`AiryState`], [`AiryParams`]), originally implemented in Fortran `ZUNHJ`.
//!
//! # References
//! - Abramowitz & Stegun, *Handbook of Mathematical Functions* (1965), Chapter 9 (§9.3.35–9.3.46, §9.7.7–9.7.10).
//! - F. W. J. Olver, *Asymptotics and Special Functions*, Academic Press (1974), pp. 376–382, 420.
//! - NIST Digital Library of Mathematical Functions (DLMF), §10.20 and §10.41.

use std::f64::consts::FRAC_PI_2;

use num::{Complex, complex::ComplexFloat};

use crate::{
    amos::{
        ComplexExt, MachineConsts,
        asymptotics::consts::{
            AIRY_ASYMP_COEFFS_A, AIRY_ASYMP_COEFFS_B, AIRY_HJ_POLYNOMIAL_COEFFS,
            DEBYE_IK_POLYNOMIAL_COEFFS, IK_NORMALIZATION_FACTORS, TRANSITION_AIRY_A_COEFFS,
            TRANSITION_AIRY_B_COEFFS, TRANSITION_AIRY_B0_COEFFS, TURNING_POINT_ZETA_COEFFS,
        },
    },
    types::BesselFloat,
};

/// Determines whether $|z|$ is sufficiently small relative to $\nu$ to trigger an underflow condition
/// in the geometric parameters.
fn geom_underflow<T: BesselFloat>(z: Complex<T>, order: T) -> bool {
    let mc: &MachineConsts<T> = T::MACHINE_CONSTANTS;

    let test = order * mc.underflow_limit;
    z.re.abs() <= test && z.im.abs() <= test
}

/// Generates fallback $(\zeta_1, \zeta_2)$ scaling parameters when an underflow condition occurs.
fn underflow_zetas<T: BesselFloat>(order: T) -> (Complex<T>, Complex<T>) {
    let mc: &MachineConsts<T> = T::MACHINE_CONSTANTS;

    let zeta1 = Complex::<T>::new(T::two() * mc.underflow_limit.ln().abs() + order, T::zero());
    let zeta2 = Complex::<T>::new(order, T::zero());
    (zeta1, zeta2)
}

/// Geometric and scaling parameters for the Debye uniform asymptotic expansion of $I_\nu(z)$ and $K_\nu(z)$.
///
/// Computes the intermediate transformation variables:
/// - $t = z / \nu$
/// - $s = 1 + t^2$
/// - $\zeta_1 = \nu \ln\left(\frac{1 + \sqrt{1 + t^2}}{t}\right)$
/// - $\zeta_2 = \nu \sqrt{1 + t^2}$
/// - $\Phi_I = \frac{(1 + t^2)^{-1/4}}{\sqrt{2\pi \nu}}, \quad \Phi_K = \sqrt{\frac{\pi}{2\nu}} (1 + t^2)^{-1/4}$
///
/// Ref: Abramowitz & Stegun (9.7.7–9.7.10), NIST DLMF (10.41). Originally part of Fortran `ZUNIK` (`IPMTR = 1`).
#[derive(Clone, Copy, Debug)]
pub(crate) struct DebyeGeometry<T: BesselFloat> {
    /// Leading prefactor $\Phi_I = \frac{(1 + (z/\nu)^2)^{-1/4}}{\sqrt{2\pi\nu}}$ for $I_\nu(z)$.
    pub phi_i: Complex<T>,
    /// Leading prefactor $\Phi_K = \sqrt{\frac{\pi}{2\nu}} (1 + (z/\nu)^2)^{-1/4}$ for $K_\nu(z)$.
    pub phi_k: Complex<T>,
    /// Logarithmic scaling exponent $\zeta_1 = \nu \ln\left(\frac{1 + \sqrt{1 + (z/\nu)^2}}{z/\nu}\right)$.
    pub zeta1: Complex<T>,
    /// Radical scaling exponent $\zeta_2 = \nu \sqrt{1 + (z/\nu)^2}$.
    pub zeta2: Complex<T>,

    /// Intermediate square $s = 1 + (z/\nu)^2$.
    s: Complex<T>,
    /// Scaled radical $s_r = \frac{1}{\nu \sqrt{1 + (z/\nu)^2}}$.
    expansion_step: Complex<T>,
    /// Reciprocal order $1 / \nu$.
    reciprocal_order: T,
    /// Indicates whether underflow occurred during parameter evaluation.
    is_underflow: bool,
}

impl<T: BesselFloat> DebyeGeometry<T> {
    /// Computes the Debye geometric and scaling parameters for a given complex argument $z$ and order $\nu$.
    pub fn compute(z: Complex<T>, order: T) -> Self {
        if geom_underflow(z, order) {
            let (zeta1, zeta2) = underflow_zetas(order);

            return Self {
                phi_i: T::C_ONE,
                phi_k: T::C_ONE,
                zeta1,
                zeta2,
                s: T::C_ONE,
                expansion_step: T::C_ZERO,
                reciprocal_order: T::ZERO,
                is_underflow: true,
            };
        }

        let reciprocal_order = T::one() / order;
        let z_over_order = z * reciprocal_order; // t = z / nu
        let s = T::C_ONE + z_over_order.powi(2); // s = 1 + t^2
        let sqrt_s = s.sqrt(); // sqrt(1 + t^2)
        let log_arg = (T::C_ONE + sqrt_s) / z_over_order; // (1 + sqrt(1 + t^2)) / t
        let zeta1 = log_arg.ln() * order; // zeta1 = nu * ln((1 + sqrt(s)) / t)
        let zeta2 = sqrt_s * order; // zeta2 = nu * sqrt(s)
        let expansion_step = Complex::<T>::from(reciprocal_order) / sqrt_s; // 1 / (nu * sqrt(1 + t^2))
        let phi_base = expansion_step.sqrt(); // unnormalized (1 + t^2)^(-1/4) / sqrt(nu)
        let phi_i = phi_base * T::from_f64(IK_NORMALIZATION_FACTORS[0]);
        let phi_k = phi_base * T::from_f64(IK_NORMALIZATION_FACTORS[1]);

        Self {
            phi_i,
            phi_k,
            zeta1,
            zeta2,
            s,
            expansion_step,
            reciprocal_order,
            is_underflow: false,
        }
    }
}

/// Full parameters for the Debye uniform asymptotic expansion of modified Bessel functions $I_\nu(z)$ and $K_\nu(z)$.
///
/// Encapsulates both the geometric scaling parameters ($\Phi_I, \Phi_K, \zeta_1, \zeta_2$) and the Olver
/// asymptotic polynomial series sums:
/// - $I_\nu(z) \sim \Phi_I e^{\zeta_2 - \zeta_1} \sum_I$
/// - $K_\nu(z) \sim \Phi_K e^{\zeta_1 - \zeta_2} \sum_K$
///
/// where $\sum_I = \sum_{k=0}^{14} \frac{u_k(t)}{\nu^k}$ and $\sum_K = \sum_{k=0}^{14} (-1)^k \frac{u_k(t)}{\nu^k}$.
///
/// Ref: Abramowitz & Stegun (9.7.7–9.7.10), NIST DLMF (10.41). Originally Fortran `ZUNIK` (`IPMTR = 0`).
#[derive(Clone, Copy, Debug)]
pub(crate) struct DebyeParams<T: BesselFloat> {
    /// Leading prefactor $\Phi_I$ for $I_\nu(z)$.
    pub phi_i: Complex<T>,
    /// Leading prefactor $\Phi_K$ for $K_\nu(z)$.
    pub phi_k: Complex<T>,
    /// Logarithmic scaling exponent $\zeta_1 = \nu \ln\left(\frac{1 + \sqrt{1 + (z/\nu)^2}}{z/\nu}\right)$.
    pub zeta1: Complex<T>,
    /// Radical scaling exponent $\zeta_2 = \nu \sqrt{1 + (z/\nu)^2}$.
    pub zeta2: Complex<T>,
    /// Olver asymptotic series sum for $I_\nu(z)$ ($\sum_{k=0}^{14} u_k(t) / \nu^k$).
    pub sum_i: Complex<T>,
    /// Olver asymptotic series sum for $K_\nu(z)$ ($\sum_{k=0}^{14} (-1)^k u_k(t) / \nu^k$).
    pub sum_k: Complex<T>,
}

impl<T: BesselFloat> DebyeParams<T> {
    /// Computes the complete Debye asymptotic parameters and Olver series sums for $z$ and order $\nu$.
    pub fn compute(z: Complex<T>, order: T) -> Self {
        let geom = DebyeGeometry::compute(z, order);
        if geom.is_underflow {
            Self {
                phi_i: geom.phi_i,
                phi_k: geom.phi_k,
                zeta1: geom.zeta1,
                zeta2: geom.zeta2,
                sum_i: T::C_ZERO,
                sum_k: T::C_ZERO,
            }
        } else {
            let mc: &MachineConsts<T> = T::MACHINE_CONSTANTS;

            let poly_arg = T::C_ONE / geom.s; // 1 / (1 + t^2), argument for Olver polynomials
            let mut sum_i = T::C_ONE;
            let mut sum_k = T::C_ONE;
            let mut term_sign = T::one();
            let mut sr_power = T::C_ONE;
            let mut recip_order_pow = T::one();

            for coeffs in DEBYE_IK_POLYNOMIAL_COEFFS.iter().skip(1) {
                // Evaluate the k-th Olver polynomial p_k(1/s) using Horner's method
                let poly_val = evaluate_horner(coeffs, poly_arg);
                sr_power *= geom.expansion_step;
                let term = sr_power * poly_val;
                sum_i += term;
                term_sign = -term_sign; // (-1)^k alternating sign for K_nu
                sum_k += term * term_sign;

                recip_order_pow *= geom.reciprocal_order;
                let test = term.l1_norm();
                // Early exit if term magnitude and (1/nu)^k drop below machine tolerance
                if recip_order_pow < mc.abs_error_tolerance && test < mc.abs_error_tolerance {
                    break;
                }
            }

            Self {
                phi_i: geom.phi_i,
                phi_k: geom.phi_k,
                zeta1: geom.zeta1,
                zeta2: geom.zeta2,
                sum_i,
                sum_k,
            }
        }
    }
}

/// Geometric and scaling parameters for Airy-type uniform asymptotic expansions of $J_\nu, Y_\nu, H^{(1)}_\nu, H^{(2)}_\nu$.
///
/// Computes the Airy argument $\arg = \zeta \nu^{2/3}$, scaling exponents $\zeta_1, \zeta_2$, and prefactor $\Phi$
/// where $\frac{2}{3} \nu \zeta^{3/2} = \zeta_1 - \zeta_2$.
///
/// Captures regime-specific scratch state in [`AiryState`] to allow efficient evaluation of the Olver $A(\zeta)$ and
/// $B(\zeta)$ series in [`AiryParams::compute`] without redundant recomputation.
///
/// Ref: Abramowitz & Stegun (9.3.35–9.3.46), NIST DLMF (10.20). Originally part of Fortran `ZUNHJ` (`IPMTR = 1`).
#[derive(Debug, Clone, Copy)]
pub(crate) struct AiryGeometry<T: BesselFloat> {
    /// Leading prefactor $\Phi = \left(\frac{4\zeta}{1 - (z/\nu)^2}\right)^{1/4} \nu^{-1/3}$.
    pub phi: Complex<T>,
    /// Scaled argument to the Airy function: $\arg = \zeta \nu^{2/3}$.
    pub arg: Complex<T>,
    /// Logarithmic scaling exponent $\zeta_1 = \nu \ln\left(\frac{1 + \sqrt{1 - (z/\nu)^2}}{z/\nu}\right)$.
    pub zeta1: Complex<T>,
    /// Radical scaling exponent $\zeta_2 = \nu \sqrt{1 - (z/\nu)^2}$.
    pub zeta2: Complex<T>,

    /// Intermediate regime state passed to [`AiryParams`].
    state: AiryState<T>,
}

impl<T: BesselFloat> AiryGeometry<T> {
    const THREE_PI_BY_2: f64 = 3.0 * FRAC_PI_2;

    /// Computes the Airy geometry and scaling parameters for a given complex argument $z$ and order $\nu$.
    pub fn compute(z: Complex<T>, order: T) -> Self {
        if geom_underflow(z, order) {
            let (zeta1, zeta2) = underflow_zetas(order);
            let phi = T::C_ONE;
            let arg = T::C_ONE;
            return Self {
                phi,
                arg,
                zeta1,
                zeta2,
                state: AiryState::Underflow,
            };
        }
        let mc: &MachineConsts<T> = T::MACHINE_CONSTANTS;

        let reciprocal_order = T::ONE / order;
        let z_over_order = z * reciprocal_order; // t = z / nu
        let recip_order_sqr = reciprocal_order * reciprocal_order;

        let order_one_third = order.powf(T::from_f64(1.0 / 3.0)); // nu^(1/3)
        let order_two_thirds = order_one_third * order_one_third; // nu^(2/3)
        let recip_order_one_third = T::one() / order_one_third; // nu^(-1/3)

        let w_sqr = T::C_ONE - z_over_order.powi(2); // w^2 = 1 - (z/nu)^2
        let abs_w_sqr = w_sqr.abs();

        let power_series = abs_w_sqr <= T::from_f64(0.25);
        if power_series {
            // Power series in w^2 near turning point (|w^2| <= 0.25, z ≈ nu)
            let mut k_max = 1;
            let mut p = [T::C_ZERO; 30];
            let mut abs_p = [T::zero(); 30];
            p[0] = T::C_ONE;
            let mut zeta_sum =
                Complex::<T>::new(T::from_f64(TURNING_POINT_ZETA_COEFFS[0]), T::zero());
            abs_p[0] = T::one();
            if abs_w_sqr >= mc.abs_error_tolerance {
                for k in 1..30 {
                    k_max = k + 1;
                    p[k] = p[k - 1] * w_sqr; // p[k] = (w^2)^k
                    zeta_sum += p[k] * T::from_f64(TURNING_POINT_ZETA_COEFFS[k]); // zeta / w^2 = sum c_k (w^2)^k
                    abs_p[k] = abs_p[k - 1] * abs_w_sqr;
                    if abs_p[k] < mc.abs_error_tolerance {
                        break;
                    }
                }
            }
            let zeta = w_sqr * zeta_sum; // zeta = w^2 * (zeta / w^2)
            let arg = zeta * order_two_thirds; // Airy argument: arg = zeta * nu^(2/3)
            let sqrt_zeta_over_w = zeta_sum.sqrt(); // sqrt(zeta / w^2) = sqrt(zeta) / w
            let w = w_sqr.sqrt();
            let zeta2 = w * order; // zeta2 = nu * w
            let zeta1 = (T::C_ONE + T::TWO_THIRDS * zeta * sqrt_zeta_over_w) * zeta2; // zeta1 = zeta2 + (2/3)*nu*zeta^(3/2)
            let phi_base = (T::two() * sqrt_zeta_over_w).sqrt(); // phi_base = (4*zeta / (1 - t^2))^(1/4)
            let phi = phi_base * recip_order_one_third; // Phi = phi_base * nu^(-1/3)
            Self {
                phi,
                arg,
                zeta1,
                zeta2,
                state: AiryState::Transition {
                    p,
                    abs_p,
                    k_max,
                    recip_order_sqr,
                    reciprocal_order,
                    recip_order_one_third,
                },
            }
        } else {
            // Asymptotic expansion away from turning point (|w^2| > 0.25)
            let mut w = w_sqr.sqrt(); // w = sqrt(1 - (z/nu)^2)
            if w.re < T::zero() {
                w.re = T::zero()
            };
            if w.im < T::zero() {
                w.im = T::zero()
            };

            let log_arg = (T::C_ONE + w) / z_over_order; // (1 + w) / (z/nu)
            let mut log_term = log_arg.ln(); // ln((1 + w)/t) = zeta1 / nu
            log_term.im = log_term.im.clamp(T::zero(), T::from_f64(FRAC_PI_2));
            if log_term.re < T::zero() {
                log_term.re = T::zero()
            };
            let zeta_3_2 = (log_term - w) * T::from_f64(1.5); // zeta^(3/2) = (3/2) * (zeta1 - zeta2)/nu
            let zeta1 = log_term * order; // zeta1 = nu * ln((1 + w)/t)
            let zeta2 = w * order; // zeta2 = nu * w
            let abs_zeta_3_2 = zeta_3_2.abs();
            let mut ang = zeta_3_2.parg();
            ang = ang.clamp(T::zero(), T::from_f64(Self::THREE_PI_BY_2));

            // Recover zeta = (zeta^(3/2))^(2/3) via polar form:
            let abs_zeta = abs_zeta_3_2.powf(T::TWO_THIRDS); // |zeta| = |zeta^(3/2)|^(2/3)
            ang *= T::TWO_THIRDS; // arg(zeta) = (2/3) * arg(zeta^(3/2))
            let mut zeta = Complex::<T>::cis(ang) * abs_zeta;
            if zeta.im < T::zero() {
                zeta.im = T::zero()
            };
            let arg = zeta * order_two_thirds; // Airy argument: arg = zeta * nu^(2/3)
            let sqrt_zeta = zeta_3_2 / zeta; // sqrt(zeta) = zeta^(3/2) / zeta
            let sqrt_zeta_over_w = sqrt_zeta / w; // sqrt(zeta) / w
            let phi_base = (T::two() * sqrt_zeta_over_w).sqrt(); // phi_base = (4*zeta / (1 - t^2))^(1/4)
            let phi = phi_base * recip_order_one_third; // Phi = phi_base * nu^(-1/3)

            Self {
                phi,
                arg,
                zeta1,
                zeta2,
                state: AiryState::Asymptotic {
                    w_sqr,
                    abs_w_sqr,
                    w,
                    zeta_3_2,
                    abs_zeta_3_2,
                    sqrt_zeta,
                    reciprocal_order,
                    recip_order_sqr,
                    recip_order_one_third,
                },
            }
        }
    }
}

/// Regime-specific intermediate scratch state preserved across Airy geometry and parameter evaluation.
///
/// Distinguishes between:
/// - [`AiryState::Transition`]: Near the turning point ($|1 - (z/\nu)^2| \le 0.25$), where power series in $w^2 = 1 - (z/\nu)^2$ are used.
/// - [`AiryState::Asymptotic`]: Away from the turning point ($|1 - (z/\nu)^2| > 0.25$), where asymptotic expansions of the Airy parameters are used.
/// - [`AiryState::Underflow`]: Extreme underflow conditions.
#[derive(Debug, Clone, Copy)]
pub(crate) enum AiryState<T: BesselFloat> {
    /// Underflow fallback condition.
    Underflow,
    /// Turning-point transition regime ($|w^2| \le 0.25$).
    Transition {
        /// Powers of $w^2$: $p[k] = w^{2k}$.
        p: [Complex<T>; 30],
        /// Magnitudes $|p[k]| = |w^2|^k$ for convergence monitoring.
        abs_p: [T; 30],
        /// Number of valid populated terms in `p` and `abs_p`.
        k_max: usize,
        /// Reciprocal order $1 / \nu$.
        reciprocal_order: T,
        /// Squared reciprocal order $1 / \nu^2$.
        recip_order_sqr: T,
        /// Order power $\nu^{-1/3}$.
        recip_order_one_third: T,
    },
    /// Asymptotic regime away from the turning point ($|w^2| > 0.25$).
    Asymptotic {
        /// Distance from turning point $w^2 = 1 - (z/\nu)^2$.
        w_sqr: Complex<T>,
        /// Magnitude $|w^2|$.
        abs_w_sqr: T,
        /// Radical $w = \sqrt{1 - (z/\nu)^2}$.
        w: Complex<T>,
        /// Scaled exponent $\zeta^{3/2} = \frac{3}{2}(\zeta_1 - \zeta_2) / \nu$.
        zeta_3_2: Complex<T>,
        /// Magnitude $|\zeta^{3/2}|$.
        abs_zeta_3_2: T,
        /// Radical $\sqrt{\zeta} = \zeta^{3/2} / \zeta$.
        sqrt_zeta: Complex<T>,
        /// Reciprocal order $1 / \nu$.
        reciprocal_order: T,
        /// Squared reciprocal order $1 / \nu^2$.
        recip_order_sqr: T,
        /// Order power $\nu^{-1/3}$.
        recip_order_one_third: T,
    },
}

/// Full parameters for the Airy-type uniform asymptotic expansion of $J_\nu, Y_\nu, H^{(1)}_\nu, H^{(2)}_\nu$.
///
/// Encapsulates the geometry ($\Phi, \arg, \zeta_1, \zeta_2$) and the Olver asymptotic series sums $A(\zeta)$ and $B(\zeta)$:
/// $$C_\nu(z) \sim C_1 \Phi \left[ A(\zeta) \text{Ai}(\text{arg}) + \frac{C_2}{\nu^{2/3}} B(\zeta) \text{Ai}'(\text{arg}) \right]$$
///
/// where $\text{arg} = \zeta \nu^{2/3}$, $A(\zeta) = \sum_{s=0}^\infty \frac{A_s(\zeta)}{\nu^{2s}}$, and $B(\zeta) = \sum_{s=0}^\infty \frac{B_s(\zeta)}{\nu^{2s}}$.
///
/// Ref: Abramowitz & Stegun (9.3.35–9.3.46), NIST DLMF (10.20). Originally Fortran `ZUNHJ` (`IPMTR = 0`).
#[derive(Debug, Clone, Copy)]
pub(crate) struct AiryParams<T: BesselFloat> {
    /// Leading prefactor $\Phi = \left(\frac{4\zeta}{1 - (z/\nu)^2}\right)^{1/4} \nu^{-1/3}$.
    pub phi: Complex<T>,
    /// Scaled argument to the Airy function: $\arg = \zeta \nu^{2/3}$.
    pub arg: Complex<T>,
    /// Logarithmic scaling exponent $\zeta_1 = \nu \ln\left(\frac{1 + \sqrt{1 - (z/\nu)^2}}{z/\nu}\right)$.
    pub zeta1: Complex<T>,
    /// Radical scaling exponent $\zeta_2 = \nu \sqrt{1 - (z/\nu)^2}$.
    pub zeta2: Complex<T>,
    /// Olver asymptotic series sum $A(\zeta)$ associated with the Airy function $\text{Ai}(\text{arg})$.
    pub asum: Complex<T>,
    /// Olver asymptotic series sum $B(\zeta)$ associated with the Airy derivative $\text{Ai}'(\text{arg})$.
    pub bsum: Complex<T>,
}

impl<T: BesselFloat> AiryParams<T> {
    /// Computes the complete Airy uniform asymptotic parameters and Olver series sums for $z$ and order $\nu$.
    pub fn compute(z: Complex<T>, order: T) -> Self {
        let geom = AiryGeometry::compute(z, order);
        let phi = geom.phi;
        let arg = geom.arg;
        let zeta1 = geom.zeta1;
        let zeta2 = geom.zeta2;
        let mc: &MachineConsts<T> = T::MACHINE_CONSTANTS;

        match geom.state {
            AiryState::Underflow => {
                let (zeta1, zeta2) = underflow_zetas(order);
                let phi = T::C_ONE;
                let arg = T::C_ONE;
                Self {
                    phi,
                    arg,
                    zeta1,
                    zeta2,
                    asum: T::C_ZERO,
                    bsum: T::C_ZERO,
                }
            }
            AiryState::Transition {
                p,
                abs_p,
                k_max,
                recip_order_sqr,
                reciprocal_order,
                recip_order_one_third,
            } => {
                // Transition regime (|1 - (z/nu)^2| <= 0.25): Power series in w^2 (NIST DLMF §10.20.7)
                // A(zeta) = 1 + sum_{s=1}^6 A_s(w^2) / nu^(2s)
                // B(zeta) = B_0(w^2) + sum_{s=1}^6 B_s(w^2) / nu^(2s)
                let b_first_term: Complex<T> = p[..k_max]
                    .iter()
                    .zip(TRANSITION_AIRY_B0_COEFFS)
                    .map(|(p, b)| p * T::from_f64(b))
                    .sum();
                let a_first_term = T::C_ONE;

                let mut asum = a_first_term;
                let mut bsum = b_first_term;

                if recip_order_sqr >= mc.abs_error_tolerance {
                    let btol = mc.abs_error_tolerance * (b_first_term.l1_norm());
                    let mut poly_tol = mc.abs_error_tolerance;

                    let mut a_converged = false;
                    let mut b_converged = false;
                    let mut recip_order_power = T::one();

                    for (a_block, b_block) in TRANSITION_AIRY_A_COEFFS
                        .iter()
                        .zip(&TRANSITION_AIRY_B_COEFFS)
                    {
                        poly_tol /= recip_order_sqr;
                        recip_order_power *= recip_order_sqr;
                        if !a_converged {
                            let a_poly =
                                evaluate_transition_poly(&p, &abs_p, a_block, k_max, poly_tol);
                            asum += a_poly * recip_order_power;
                            if recip_order_power < mc.abs_error_tolerance {
                                a_converged = true
                            };
                        }
                        if !b_converged {
                            let b_poly =
                                evaluate_transition_poly(&p, &abs_p, b_block, k_max, poly_tol);
                            bsum += b_poly * recip_order_power;
                            if recip_order_power < btol {
                                b_converged = true;
                            }
                        }
                        if a_converged && b_converged {
                            break;
                        }
                    }
                }
                bsum *= reciprocal_order * recip_order_one_third;

                Self {
                    phi,
                    arg,
                    zeta1,
                    zeta2,
                    asum,
                    bsum,
                }
            }
            AiryState::Asymptotic {
                w_sqr,
                abs_w_sqr,
                w,
                zeta_3_2,
                abs_zeta_3_2,
                sqrt_zeta,
                reciprocal_order,
                recip_order_sqr,
                recip_order_one_third,
            } => {
                // Asymptotic regime (|1 - (z/nu)^2| > 0.25): Olver Debye-Airy expansion (NIST DLMF §10.20(iii))
                // A(zeta) = sum_{s=0}^6 A_s / nu^(2s),  B(zeta) = sum_{s=0}^6 B_s / nu^(2s)
                // where A_s and B_s are discrete convolutions of Debye polynomials u_k(p) and Airy coefficients.

                let recip_abs_w_sqr = T::one() / abs_w_sqr;
                let recip_w_sqr = w_sqr.conj() * recip_abs_w_sqr * recip_abs_w_sqr;
                let recip_nu_w = w.conj() * recip_abs_w_sqr * reciprocal_order;

                let recip_abs_zeta_3_2 = T::one() / abs_zeta_3_2;
                let recip_nu_zeta_3_2 =
                    zeta_3_2.conj() * recip_abs_zeta_3_2 * recip_abs_zeta_3_2 * reciprocal_order;

                let airy_a1_term = recip_nu_zeta_3_2 * T::from_f64(AIRY_ASYMP_COEFFS_A[1]);
                let mut u_polys = [T::C_ZERO; 14];
                u_polys[1] =
                    evaluate_horner(AIRY_HJ_POLYNOMIAL_COEFFS[1], recip_w_sqr) * recip_nu_w;

                let first_b_term = u_polys[1] + airy_a1_term;
                let first_a_term = T::C_ONE;

                let mut asum = first_a_term;
                let mut bsum = first_b_term;

                if reciprocal_order >= mc.abs_error_tolerance {
                    let mut recip_nu_zeta_3_2_power = recip_nu_zeta_3_2;
                    let mut recip_nu_w_power = recip_nu_w;
                    u_polys[0] = T::C_ONE;
                    let mut recip_order_power = T::one();
                    let btol = mc.abs_error_tolerance * (bsum.l1_norm());

                    let mut a_converged = false;
                    let mut b_converged = false;
                    let mut a_coeffs = [T::C_ZERO; 12];
                    let mut b_coeffs = [T::C_ZERO; 12];
                    for lr in (2..=12).step_by(2) {
                        // Step s = lr / 2: compute 2 new polynomial terms u_{2s}, u_{2s+1} and Airy scaling factors
                        for k in lr..=lr + 1 {
                            // Horner's method on polynomial u_k(1/w^2)
                            let poly_val =
                                evaluate_horner(AIRY_HJ_POLYNOMIAL_COEFFS[k], recip_w_sqr);
                            recip_nu_w_power *= recip_nu_w;
                            u_polys[k] = recip_nu_w_power * poly_val;

                            // Airy derivative coefficient v_{k-1} / (nu * zeta^(3/2))^(k-1)
                            a_coeffs[k - 2] =
                                recip_nu_zeta_3_2_power * T::from_f64(AIRY_ASYMP_COEFFS_B[k - 1]);
                            recip_nu_zeta_3_2_power *= recip_nu_zeta_3_2;
                            // Airy function coefficient u_k / (nu * zeta^(3/2))^k
                            b_coeffs[k - 2] =
                                recip_nu_zeta_3_2_power * T::from_f64(AIRY_ASYMP_COEFFS_A[k]);
                        }
                        recip_order_power *= recip_order_sqr;

                        // A_s = u_{2s} + sum_{j=1}^{2s} v_j * u_{2s - j}
                        if !a_converged {
                            let a_term =
                                u_polys[lr] + convolve_asymptotic_series(&a_coeffs, &u_polys, lr);
                            asum += a_term;
                            if recip_order_power < mc.abs_error_tolerance
                                && a_term.l1_norm() < mc.abs_error_tolerance
                            {
                                a_converged = true
                            }
                        }

                        // B_s = u_{2s+1} + u_{2s} * a_1 + sum_{j=1}^{2s} u_{j+1} * u_{2s - j}
                        if !b_converged {
                            let b_term = u_polys[lr + 1]
                                + u_polys[lr] * airy_a1_term
                                + convolve_asymptotic_series(&b_coeffs, &u_polys, lr);
                            bsum += b_term;
                            if recip_order_power < btol && b_term.l1_norm() < btol {
                                b_converged = true;
                            }
                        }
                        if a_converged && b_converged {
                            break;
                        }
                    }
                }

                bsum = (-bsum * recip_order_one_third) / sqrt_zeta;

                Self {
                    phi,
                    arg,
                    zeta1,
                    zeta2,
                    asum,
                    bsum,
                }
            }
        }
    }
}

/// Evaluates a transition polynomial $P(w^2) = \sum_{k=0}^{k_{\max}-1} c_k (w^2)^k$
/// using precomputed powers $p[k] = (w^2)^k$, with early exit on convergence.
///
/// Ref: NIST DLMF §10.20.7.
#[inline]
fn evaluate_transition_poly<T: BesselFloat>(
    p: &[Complex<T>; 30],
    abs_p: &[T; 30],
    coeffs: &[f64; 30],
    k_max: usize,
    tolerance: T,
) -> Complex<T> {
    let mut value = T::C_ZERO;
    for k in 0..k_max {
        value += p[k] * T::from_f64(coeffs[k]);
        if abs_p[k] < tolerance {
            break;
        }
    }
    value
}

/// Evaluates a polynomial $P(x) = \sum_{j=0}^n c_j x^{n-j}$ with real coefficients
/// using Horner's method: $P(x) = (\dots((c_0 x + c_1) x + c_2) \dots) x + c_n$.
#[inline]
fn evaluate_horner<T: BesselFloat>(coeffs: &[f64], x: Complex<T>) -> Complex<T> {
    let mut val = Complex::<T>::new(T::from_f64(coeffs[0]), T::zero());
    for &c in &coeffs[1..] {
        val = val * x + T::from_f64(c);
    }
    val
}

/// Computes the discrete convolution between Airy asymptotic scaling coefficients $c_j$
/// and the Olver polynomials $u_k(p)$:
/// $$\sum_{j=1}^{\text{len}} c_j \cdot u_{\text{len} - j}(p)$$
///
/// Ref: NIST DLMF §10.20(iii).
#[inline]
fn convolve_asymptotic_series<T: BesselFloat>(
    coeffs: &[Complex<T>],
    u_polys: &[Complex<T>],
    len: usize,
) -> Complex<T> {
    coeffs[..len]
        .iter()
        .zip(u_polys[..len].iter().rev())
        .map(|(&c, &u)| c * u)
        .sum()
}
