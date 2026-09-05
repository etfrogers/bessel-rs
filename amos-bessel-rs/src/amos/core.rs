use num::{Complex, Integer, complex::ComplexFloat};

use crate::{
    BesselError::{self, *},
    BesselFloat, Scaling,
    amos::{
        ComplexExt, HankelKind, IKType, MachineConsts, RotationDirection,
        airy::airy_power_series,
        analytic_continuation::{airy_analytic_continuation, analytic_continuation},
        asymptotics::k_asymp_large_order,
        i_pow_n,
        limits::check_underflow_uniform_asymp_params,
        right_half_plane::{i_right_half_plane, k_right_half_plane},
        utils::{is_significance_lost, sanitise_inputs},
    },
    types::BesselResult,
};

/// Core Amos implementation for Hankel functions $H_\nu^{(1)}(z)$ and $H_\nu^{(2)}(z)$ ($\nu \ge 0$).
///
/// Corresponds to Amos routine `ZBESK` / `CBESH`.
///
/// # Algorithm & Mechanics
/// - Computes $H_\nu^{(m)}(z)$ via $K_\nu$:
///   $$H_\nu^{(m)}(z) = -\text{fmm} \cdot \frac{i}{\pi/2} z_t^\nu K_\nu(-z \cdot z_t)$$
///   where $z_t = \exp(-i \cdot \text{fmm} \cdot \pi / 2)$ and $\text{fmm} = 3 - 2m$.
/// - **Right half plane** ($\text{Re}(z_t \cdot z) > 0$): computes $K_\nu$ directly via `k_right_half_plane`.
/// - **Left half plane**: applies analytic continuation via `analytic_continuation`.
/// - **Large orders** ($\nu > \text{asymptotic\_order\_limit}$): uses uniform asymptotic expansion `k_asymp_large_order`.
///
/// Precondition: `order >= 0.0`, `n >= 1`, $z \ne 0$.
pub(crate) fn complex_bessel_h<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    hankel_kind: HankelKind,
    n: usize,
) -> BesselResult<T> {
    sanitise_inputs(z, order, n, true)?;
    let mc: &MachineConsts<T> = T::MACHINE_CONSTANTS;
    let mut n_zeros = 0;

    let max_order = order + T::from_usize(n - 1);

    let rotation = hankel_kind.get_rotation();
    let rotation_factor = -T::I * rotation.to_float::<T>();
    let mut z_rotated = rotation_factor * z;
    let abs_z = z.abs();
    let partial_loss_of_significance = is_significance_lost(abs_z, max_order, false, mc)?;

    // Test for overflow on the maximum order
    if abs_z < mc.underflow_limit {
        return Err(Overflow);
    }
    let (mut y, n_zeros) = if order < mc.asymptotic_order_limit {
        if max_order > T::ONE {
            if max_order > T::TWO {
                let mut y = T::c_zeros(n);
                let n_underflow = check_underflow_uniform_asymp_params(
                    z_rotated,
                    order,
                    scaling,
                    IKType::K,
                    n,
                    &mut y,
                    mc,
                )?;

                n_zeros += n_underflow;

                if n == n_underflow {
                    return if z_rotated.re < T::ZERO {
                        Err(Overflow)
                    } else if partial_loss_of_significance {
                        Err(PartialLossOfSignificance { y, n_zeros })
                    } else {
                        Ok((y, n_zeros))
                    };
                }
            }
            if abs_z <= mc.abs_error_tolerance
                && -max_order * (T::HALF * abs_z).ln() > mc.exponent_limit
            {
                return Err(Overflow);
            }
        }
        // z_rotated is in the right half plane (or on the positive imaginary axis for H1)
        if z_rotated.re > T::ZERO
            || (z_rotated.re == T::ZERO
                && z_rotated.im > T::ZERO
                && hankel_kind == HankelKind::First)
        {
            // Right half plane: compute K directly
            k_right_half_plane(z_rotated, order, scaling, n)?
        } else {
            // Left half plane: use analytic continuation
            analytic_continuation(z_rotated, order, scaling, -rotation, n)?
        }
    } else {
        // Large order: use uniform asymptotic expansion.
        // If z_rotated is in the left half plane, negate it and set a rotation
        // so the asymptotic expansion can work in the right half plane.
        let mut asymptotic_rotation = RotationDirection::None;
        if (z_rotated.re < T::ZERO)
            || (z_rotated.re == T::ZERO
                && z_rotated.im < T::ZERO
                && hankel_kind == HankelKind::Second)
        {
            asymptotic_rotation = -rotation;
            if !(z_rotated.re != T::ZERO || z_rotated.im >= T::ZERO) {
                z_rotated = -z_rotated;
            }
        }
        let (y, n_zeros_k) =
            k_asymp_large_order(z_rotated, order, scaling, asymptotic_rotation, n)?;
        n_zeros += n_zeros_k;
        (y, n_zeros)
    };

    // Convert K results to H via: H_m(ν,z) = -fmm·(i/(π/2))·zₜᵛ·K(ν, -z·zₜ)
    // where zₜ = exp(-i·fmm·π/2) = -i·fmm, fmm = 3 - 2m
    let sign = -T::FRAC_PI_2() * T::from_f64(rotation.signum());
    // Compute exp(i·ν·π/2) via order mod 2 to avoid significance loss for large orders
    let arg = (order % T::TWO) * sign;
    let mut phase_multiplier = (T::ONE / sign) * T::I * Complex::<T>::cis(arg);
    if (order.to_i64().unwrap() / 2).is_odd() {
        phase_multiplier = -phase_multiplier;
    }

    for yi in y.iter_mut().take(n - n_zeros) {
        *yi = safe_multiply(*yi, phase_multiplier, mc);
        phase_multiplier *= rotation_factor;
    }
    if partial_loss_of_significance {
        Err(PartialLossOfSignificance { y, n_zeros })
    } else {
        Ok((y, n_zeros))
    }
}

/// Core Amos implementation for modified Bessel functions $I_\nu(z)$ ($\nu \ge 0$).
///
/// Corresponds to Amos routine `ZBESI` / `CBESI`.
///
/// # Algorithm & Mechanics
/// - Evaluates in the right half-plane ($\text{Re}(z) \ge 0$) using `i_right_half_plane`:
///   - Power series for small $|z| \le 2\sqrt{\nu + 1}$.
///   - Asymptotic expansion for large $|z|$.
///   - Miller algorithm normalized by Wronskian / Neumann series for intermediate $|z|$.
///   - Uniform asymptotic expansion for large $\nu$.
/// - **Left half plane** ($\text{Re}(z) < 0$): continued via $I_\nu(z) = \exp(\pm i\pi\nu) I_\nu(-z)$.
///
/// Precondition: `order >= 0.0`, `n >= 1`.
pub(crate) fn complex_bessel_i<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    n: usize,
) -> BesselResult<T, usize> {
    sanitise_inputs(z, order, n, false)?;
    let mc: &MachineConsts<T> = T::MACHINE_CONSTANTS;

    let abs_z = z.abs();
    let max_order = order + T::from_usize(n - 1);
    let partial_significance_loss = is_significance_lost(abs_z, max_order, false, mc)?;

    let (z_right_half_plane, mut continuation_phase) = if z.re >= T::ZERO {
        (z, T::C_ONE)
    } else {
        // Compute exp(i·ν·π) via fractional part to avoid significance loss for large orders
        let integer_order = order.to_usize().unwrap();
        let arg = order.fract() * T::PI() * if z.im < T::ZERO { -T::ONE } else { T::ONE };
        let mut continuation_phase = Complex::<T>::cis(arg);
        if !integer_order.is_even() {
            continuation_phase = -continuation_phase;
        }
        (-z, continuation_phase)
    };
    let (mut y, n_zeros) = i_right_half_plane(z_right_half_plane, order, scaling, n)?;
    let remaining_n = n - n_zeros;
    if z.re < T::ZERO && remaining_n > 0 {
        // Left half plane: apply continuation I(ν,z) = exp(±iπν)·I(ν,-z)
        for yi in y.iter_mut().take(remaining_n) {
            *yi = safe_multiply(*yi, continuation_phase, mc);
            continuation_phase = -continuation_phase;
        }
    }

    if partial_significance_loss {
        Err(PartialLossOfSignificance { y, n_zeros })
    } else {
        Ok((y, n_zeros))
    }
}

/// Core Amos implementation for Bessel functions of the first kind $J_\nu(z)$ ($\nu \ge 0$).
///
/// Corresponds to Amos routine `ZBESJ` / `CBESJ`.
///
/// # Algorithm & Mechanics
/// - Transforms $J_\nu(z)$ to $I_\nu$:
///   $$J_\nu(z) = \exp(i\nu\pi/2) I_\nu(-iz) \quad (\text{Im}(z) \ge 0)$$
///   and uses conjugate symmetry for $\text{Im}(z) < 0$.
/// - Evaluates $I_\nu$ in the right half-plane via `i_right_half_plane`.
/// - Multiplies by phase $\exp(i\nu\pi/2)$ using order mod 2 arithmetic to prevent loss of significance.
///
/// Precondition: `order >= 0.0`, `n >= 1`.
pub(crate) fn complex_bessel_j<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    n: usize,
) -> BesselResult<T> {
    sanitise_inputs(z, order, n, false)?;
    let mc: &MachineConsts<T> = T::MACHINE_CONSTANTS;

    let partial_significance_loss =
        is_significance_lost(z.abs(), order + T::from_usize(n - 1), false, mc)?;
    // Compute exp(i·ν·π/2) via order mod 2 to avoid significance loss for large orders
    let arg = (order % T::TWO) * T::FRAC_PI_2();
    let mut phase_multiplier = Complex::<T>::cis(arg);
    if (order.to_i64().unwrap() / 2).is_odd() {
        phase_multiplier = -phase_multiplier;
    }
    // J(ν,z) = exp(iνπ/2)·I(ν,-iz) for Im(z) ≥ 0; conjugate symmetry handles Im(z) < 0
    let mut sign_selector = T::ONE;
    let mut z_rotated = -T::I * z;
    if z.im < T::ZERO {
        z_rotated = -z_rotated;
        phase_multiplier.im = -phase_multiplier.im;
        sign_selector = -sign_selector;
    }
    let (mut y, n_zeros) = i_right_half_plane(z_rotated, order, scaling, n)?;
    for yi in y.iter_mut().take(n - n_zeros) {
        *yi = safe_multiply(*yi, phase_multiplier, mc);
        phase_multiplier *= T::I * sign_selector;
    }
    if partial_significance_loss {
        Err(PartialLossOfSignificance { y, n_zeros })
    } else {
        Ok((y, n_zeros))
    }
}

/// Core Amos implementation for modified Bessel functions of the second kind $K_\nu(z)$ ($\nu \ge 0$).
///
/// Corresponds to Amos routine `ZBESK` / `CBESK`.
///
/// # Algorithm & Mechanics
/// - **Small orders**: computes $K_\nu$ and $K_{\nu+1}$ in the right half-plane via `k_right_half_plane`
///   and uses forward recurrence for higher orders.
/// - **Left half plane** ($\text{Re}(z) < 0$): continued via
///   $$K_\nu(z e^{\pm i\pi}) = e^{\mp i\pi\nu} K_\nu(z) \mp i\pi I_\nu(z)$$
/// - **Large orders** ($\nu > \text{asymptotic\_order\_limit}$): uniform asymptotic expansions via `k_asymp_large_order`.
///
/// Precondition: `order >= 0.0`, `n >= 1`, $z \ne 0$.
pub(crate) fn complex_bessel_k<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    n: usize,
) -> BesselResult<T> {
    sanitise_inputs(z, order, n, true)?;
    let mc: &MachineConsts<T> = T::MACHINE_CONSTANTS;
    let abs_z = z.abs();
    let max_order = order + T::from_usize(n - 1);
    let partial_significance_loss = is_significance_lost(abs_z, max_order, false, mc)?;

    // Overflow: K diverges as z → 0
    if abs_z < mc.underflow_limit {
        return Err(Overflow);
    }

    let mut n_zeros = 0;
    if order > mc.asymptotic_order_limit {
        // Large order: use uniform asymptotic expansion
        let rotation = if z.re >= T::ZERO {
            RotationDirection::None
        } else if z.im < T::ZERO {
            RotationDirection::Left
        } else {
            RotationDirection::Right
        };

        let (y, n_zeros) = k_asymp_large_order(z, order, scaling, rotation, n)?;
        return if partial_significance_loss {
            Err(PartialLossOfSignificance { y, n_zeros })
        } else {
            Ok((y, n_zeros))
        };
    }

    if max_order > T::TWO {
        let mut y = T::c_zeros(n);
        let n_underflow =
            check_underflow_uniform_asymp_params(z, order, scaling, IKType::K, n, &mut y, mc)?;
        n_zeros += n_underflow;

        if n_underflow == n {
            return if z.re < T::ZERO {
                Err(Overflow)
            } else if partial_significance_loss {
                Err(PartialLossOfSignificance { y, n_zeros })
            } else {
                Ok((y, n_zeros))
            };
        }
    }
    // For very small |z| and large order, K grows as (z/2)^{-ν}, check this doesn't overflow
    if (max_order > T::ONE) && abs_z <= mc.abs_error_tolerance {
        let half_abs_z = T::HALF * abs_z;
        if -max_order * half_abs_z.ln() > mc.exponent_limit {
            return Err(Overflow);
        }
    }
    let (y, n_zeros) = if z.re >= T::ZERO {
        // Right half plane
        k_right_half_plane(z, order, scaling, n)?
    } else {
        // Left half plane: use analytic continuation
        // If any orders already underflowed, the continuation will overflow
        if n_zeros != 0 {
            return Err(Overflow);
        }
        let rotation = if z.im < T::ZERO {
            RotationDirection::Left
        } else {
            RotationDirection::Right
        };
        analytic_continuation(z, order, scaling, rotation, n)?
    };
    if partial_significance_loss {
        Err(PartialLossOfSignificance { y, n_zeros })
    } else {
        Ok((y, n_zeros))
    }
}

/// Core Amos implementation for Bessel functions of the second kind $Y_\nu(z)$ ($\nu \ge 0$).
///
/// Corresponds to Amos routine `ZBESY` / `CBESY`.
///
/// # Algorithm & Mechanics
/// - Computes $Y_\nu(z)$ in the right half-plane via $I_\nu$ and $K_\nu$:
///   $$Y_\nu(z) = i e^{i\pi\nu/2} I_\nu(-iz) - \frac{2}{\pi} e^{-i\pi\nu/2} K_\nu(-iz) \quad (\text{Im}(z) \ge 0)$$
/// - For $\text{Im}(z) < 0$, uses conjugate symmetry $Y_\nu(z) = \overline{Y_\nu(\bar{z})}$.
/// - Handles exponential scaling correction when `Scaling::Scaled` is requested.
///
/// Precondition: `order >= 0.0`, `n >= 1`, $z \ne 0$.
pub(crate) fn complex_bessel_y<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    n: usize,
) -> BesselResult<T> {
    sanitise_inputs(z, order, n, true)?;
    let mc: &MachineConsts<T> = T::MACHINE_CONSTANTS;
    // Use conjugate symmetry: Y(ν,z) = conj(Y(ν,conj(z))) for Im(z) < 0
    let z_upper_half_plane = if z.im < T::ZERO { z.conj() } else { z };
    let z_rotated = -T::I * z_upper_half_plane;
    let mut partial_loss_of_significance = false;

    let mut unwrap_psl = |result: BesselResult<T>| match result {
        Ok((y_, n_zeros_)) => Ok((y_, n_zeros_)),
        Err(PartialLossOfSignificance {
            y: y_,
            n_zeros: n_zeros_,
        }) => {
            partial_loss_of_significance = true;
            Ok((y_, n_zeros_))
        }
        err => err,
    };

    let (bess_i, n_zeros_i) = unwrap_psl(complex_bessel_i(z_rotated, order, scaling, n))?;
    let (bess_k, n_zeros_k) = unwrap_psl(complex_bessel_k(z_rotated, order, scaling, n))?;

    let mut n_zeros = n_zeros_i.min(n_zeros_k);
    let frac_order = order.fract();
    let integer_order = order.to_usize().unwrap();
    let mut i_coeff = Complex::<T>::cis(T::FRAC_PI_2() * frac_order);
    i_coeff *= i_pow_n(integer_order);
    let mut k_coeff = i_coeff.conj() * T::FRAC_2_PI();
    i_coeff *= T::I;

    let mut exponential_correction = T::ONE;
    if scaling == Scaling::Scaled {
        let phase_correction = Complex::<T>::cis(z.re);
        let two_abs_z = T::TWO * z.im.abs();
        exponential_correction = if two_abs_z < mc.exponent_limit {
            (-two_abs_z).exp()
        } else {
            T::ZERO
        };
        k_coeff *= phase_correction * exponential_correction;
        n_zeros = 0;
    }
    let mut y: Vec<Complex<T>> = bess_i
        .iter()
        .zip(bess_k)
        .map(|(&z_i, z_k)| {
            let z_k = scaled_multiply(z_k, k_coeff, scaling, mc);
            let z_i = scaled_multiply(z_i, i_coeff, scaling, mc);
            let val = z_i - z_k;
            if scaling == Scaling::Scaled && val == T::C_ZERO && exponential_correction == T::ZERO {
                n_zeros += 1;
            }
            i_coeff *= T::I;
            k_coeff *= -T::I;
            val
        })
        .collect();

    if z.im < T::ZERO {
        y.iter_mut().for_each(|v| *v = v.conj());
    }
    if partial_loss_of_significance {
        Err(PartialLossOfSignificance { y, n_zeros })
    } else {
        Ok((y, n_zeros))
    }
}

#[inline]
fn scaled_multiply<T: BesselFloat>(
    z: Complex<T>,
    coeff: Complex<T>,
    scaling: Scaling,
    mc: &MachineConsts<T>,
) -> Complex<T> {
    match scaling {
        Scaling::Unscaled => z * coeff,
        Scaling::Scaled => safe_multiply(z, coeff, mc),
    }
}

#[inline]
fn safe_multiply<T: BesselFloat>(
    z: Complex<T>,
    coeff: Complex<T>,
    mc: &MachineConsts<T>,
) -> Complex<T> {
    if z.linf_norm() <= mc.absolute_approximation_limit {
        (z * mc.rtol) * mc.abs_error_tolerance
    } else {
        z * coeff
    }
}

/// Computes the Airy function Ai(z) or its derivative dAi(z)/dz for a complex argument.
///
/// This function computes the complex Airy function Ai(z) or its derivative dAi(z)/dz.
/// A scaling option is provided to remove the exponential decay in `-PI/3 < z.arg() < PI/3`
/// and the exponential growth in `PI/3 < z.arg().abs() < PI`.
///
/// While the Airy functions Ai(z) and dAi(z)/dz are analytic in the whole z-plane,
/// the corresponding scaled functions have a cut along the negative real axis.
///
/// Ai(z) AND dAi(z)/dz are computed for `z.abs() > 1.0` from the K Bessel functions
/// by the formulae:
///
/// `Ai(z) = c * z.sqrt() * K(1/3, zeta)`, and
/// `dAi(z)/dz = -c * z * K(2/3, zeta)`
///
/// where `c = 1.0 / (PI * (3.0).sqrt())`
/// and `zeta = (2/3) * z.powf(3/2)`
///
/// and with the power series for `z.abs() <= 1.0`.
///
/// # Arguments
///
/// * `z` - Complex argument `z`.
/// * `return_derivative` - A boolean indicating whether to compute the derivative.
///     * `false`: computes `Ai(z)`.
///     * `true`: computes `dAi(z)/dz`.
/// * `scaling` - A parameter to indicate the scaling option.
///     * `Scaling::Unscaled`: returns `Ai(z)` or `dAi(z)/dz`.
///     * `Scaling::Scaled`: returns `zeta.exp() * Ai(z)` or `zeta.exp() * dAi(z)/dz`,
///       where `zeta = (2/3) * z * z.sqrt()`.
///
/// # Returns
///
/// A tuple containing:
/// * The complex result of the Airy function computation.
/// * An underflow indicator (`0` for normal return, `1` for underflow).
pub fn complex_airy<T: BesselFloat>(
    z: Complex<T>,
    return_derivative: bool,
    scaling: Scaling,
) -> Result<(Complex<T>, usize), BesselError<T>> {
    const POWER_SERIES_COEFFS: (f64, f64) = (3.550_280_538_878_172e-1, 2.588_194_037_928_068e-1);
    const FRAC_1_PI_SQRT_3: f64 = 1.837_762_984_739_306_8e-1;
    let mc: &MachineConsts<T> = T::MACHINE_CONSTANTS;

    let abs_z = z.abs();
    // significance loss only tested against z, not order, so 0.0 is used to never cause significance loss
    let partial_loss_of_significance = is_significance_lost(abs_z, T::ZERO, true, mc)?;

    let return_values = if abs_z <= T::ONE {
        // Power series for small |z|
        let ai = airy_power_series(z, return_derivative, POWER_SERIES_COEFFS);
        (
            match scaling {
                Scaling::Scaled => ai * (T::TWO_THIRDS * z * z.sqrt()).exp(),
                Scaling::Unscaled => ai,
            },
            0,
        )
    } else {
        // Large |z|: use K Bessel functions
        let order = if return_derivative {
            T::TWO_THIRDS
        } else {
            T::ONE_THIRD
        };
        let ln_abs_z = abs_z.ln();

        let sqrt_z = z.sqrt();
        let mut zeta = T::TWO_THIRDS * z * sqrt_z;
        // Ensure Re(ζ) ≤ 0 when Re(z) < 0 (especially for small Im(z))
        let mut scale_factor = T::ONE;
        if z.re < T::ZERO {
            zeta.re = -zeta.re.abs();
        }
        if z.im == T::ZERO && z.re <= T::ZERO {
            zeta.re = T::ZERO;
        }
        let re_zeta = zeta.re;
        let (y, n_zeros) = if re_zeta < T::ZERO || z.re <= T::ZERO {
            // Overflow test for unscaled mode
            if scaling == Scaling::Unscaled && re_zeta <= -mc.approximation_limit {
                scale_factor = mc.abs_error_tolerance;
                if (-re_zeta + T::from_f64(0.25) * ln_abs_z) > mc.exponent_limit {
                    return Err(Overflow);
                }
            }
            // In scaled mode, k_right_half_plane and analytic_continuation return exp(ζ)·K(ν,ζ)
            let rotation = if z.im < T::ZERO {
                RotationDirection::Left
            } else {
                RotationDirection::Right
            };
            airy_analytic_continuation(zeta, order, scaling, rotation)?
        } else {
            // Underflow test for unscaled mode
            let mut retval = None;
            if scaling == Scaling::Unscaled && re_zeta > mc.approximation_limit {
                scale_factor = T::ONE / mc.abs_error_tolerance;
                if (-re_zeta - T::from_f64(0.25) * ln_abs_z) < -mc.exponent_limit {
                    retval = Some(Ok((T::c_zeros(1), 1)));
                }
            }
            retval.unwrap_or_else(|| k_right_half_plane(zeta, order, scaling, 1))?
        };

        let mut y = y[0] * T::from_f64(FRAC_1_PI_SQRT_3) * scale_factor;
        y *= if return_derivative { -z } else { sqrt_z };
        (y / scale_factor, n_zeros)
    };

    if partial_loss_of_significance {
        Err(PartialLossOfSignificance {
            y: vec![return_values.0],
            n_zeros: return_values.1,
        })
    } else {
        Ok(return_values)
    }
}

/// Computes the Airy function Bi(z) or its derivative dBi(z)/dz for a complex argument.
///
/// This function computes the complex Airy function Bi(z) or its derivative dBi(z)/dz.
/// A scaling option is provided to remove the exponential behavior in both the left
/// and right half-planes.
///
/// Bi and dBi are computed for `z.abs() > 1.0` from the I Bessel functions by
///
/// Bi(z) = c* z.sqrt() * ( I(-1/3, zeta) + I(1/3, zeta) )
/// dBi(z) = c *  z  * ( I(-2/3, zeta) + I(2/3, zeta) )
///
/// where `c = 1.0 / (3.0).sqrt()` and `zeta = (2/3) * z.powf(3/2)`
///
/// and with the power series for `z.abs() <= 1.0`.
///
///
/// # Arguments
///
/// * `z` - Complex argument `z`.
/// * `return_derivative` - A boolean indicating whether to compute the derivative.
///     * `false`: computes `Bi(z)`.
///     * `true`: computes `dBi(z)/dz`.
/// * `scaling` - A parameter to indicate the scaling option.
///     * `Scaling::Unscaled`: returns `Bi(z)` or `dBi(z)/dz`.
///     * `Scaling::Scaled`: returns `(-zeta.re.abs()).exp() * Bi(z)` or `(-zeta.re.abs()).exp() * dBi(z)/dz`,
///       where `zeta = (2/3) * z.powf(3/2)`.
///
/// # Returns
///
/// The complex result of the Airy function computation.
pub fn complex_airy_b<T: BesselFloat>(
    z: Complex<T>,
    return_derivative: bool,
    scaling: Scaling,
) -> Result<Complex<T>, BesselError<T>> {
    const POWER_SERIES_COEFFS: (f64, f64) = (6.149_266_274_460_007e-1, -4.482_883_573_538_264e-1);
    const FRAC_1_SQRT_3: f64 = 5.773_502_691_896_257e-1;
    let mc: &MachineConsts<T> = T::MACHINE_CONSTANTS;

    let (order1, order2) = if return_derivative {
        (T::TWO_THIRDS, T::ONE_THIRD)
    } else {
        (T::ONE_THIRD, T::TWO_THIRDS)
    };

    let abs_z = z.abs();
    let mut partial_loss_of_significance = false;

    let y = if abs_z <= T::ONE {
        // Power series for small |z|
        let y = airy_power_series(z, return_derivative, POWER_SERIES_COEFFS);
        match scaling {
            Scaling::Scaled => {
                //TODO zeta used many places with similar definition
                let zeta = T::TWO_THIRDS * (z * z.sqrt());
                y * (-(zeta.re.abs())).exp()
            }
            Scaling::Unscaled => y,
        }
    } else {
        // Large |z|: use I Bessel functions
        // significance loss only tested against z, not order, so 0.0 is used to never cause significance loss
        partial_loss_of_significance = is_significance_lost(abs_z, T::ZERO, true, mc)?;
        let mut scale_factor = T::ONE;
        let mut zeta = T::TWO_THIRDS * (z * z.sqrt());

        // Ensure Re(ζ) ≤ 0 when Re(z) < 0 (especially for small Im(z))
        if z.re < T::ZERO {
            zeta.re = -zeta.re.abs();
        }
        if z.im == T::ZERO && z.re < T::ZERO {
            zeta.re = T::ZERO;
        }
        if scaling == Scaling::Unscaled {
            // Overflow test for unscaled mode
            let re_zeta = zeta.re.abs();
            if re_zeta > mc.approximation_limit {
                scale_factor = mc.abs_error_tolerance;
                if re_zeta + T::from_f64(0.25) * abs_z.ln() > mc.exponent_limit {
                    return Err(Overflow);
                }
            }
        }
        let mut rotation_angle = T::ZERO;
        if zeta.re < T::ZERO || z.re <= T::ZERO {
            rotation_angle = T::PI();
            if z.im < T::ZERO {
                rotation_angle = -T::PI();
            }
            zeta *= -T::ONE;
        }
        // Compute I(ν₁,ζ) and I(ν₂,ζ); in scaled mode these return exp(-|Re(ζ)|)·I(ν,ζ)
        // rotation_angle provides the analytic continuation factor for left half plane
        let (i1, _) = i_right_half_plane(zeta, order1, scaling, 1)?;
        let i_pos_term = Complex::<T>::cis(rotation_angle * order1) * i1[0] * scale_factor;
        let (mut i2, _) = i_right_half_plane(zeta, order2, scaling, 2)?;
        i2[0] *= scale_factor;
        i2[1] *= scale_factor;

        // Backward recurrence one step for negative order: I(-ν,ζ) = (2ν/ζ)·I(ν,ζ) + I(ν+1,ζ)
        let i_neg_term = (T::TWO * order2) * (i2[0] / zeta) + i2[1];
        let bi_unscaled = T::from_f64(FRAC_1_SQRT_3)
            * (i_pos_term + i_neg_term * Complex::<T>::cis(rotation_angle * (order2 - T::ONE)));
        let z_factor = if return_derivative { z } else { z.sqrt() };
        bi_unscaled * z_factor / scale_factor
    };

    if partial_loss_of_significance {
        Err(PartialLossOfSignificance {
            y: vec![y],
            n_zeros: 0,
        })
    } else {
        Ok(y)
    }
}
