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

/// Computes the H-Bessel functions (Hankel functions) of a complex argument.
///
/// This function computes a sequence of complex Hankel (Bessel) functions
/// `y[j] = H(order + j - 1, z)` real, non-negative
/// orders `order + j - 1` (`j = 1, ..., n`), and a complex argument `z` which is
/// not equal to `(0.0, 0.0)`. The computation is valid in the cut plane
/// `-PI < z.arg() <= PI`.
///
/// The kind of the Hankel function is specified by the hankel_kind parameter,
/// which can take values [HankelKind::First] or [HankelKind::Second]
///
/// When `scaling` is `Scaling::Scaled`, this function returns the scaled Hankel
/// functions, which remove the exponential behavior in both the upper and
/// lower half-planes.
/// `y(j) = (-(3 - 2 * m)*z*i).exp() * H(order + j - 1, z)` where `m` depends
/// on the kind of Hankel function (1 for First, 2 for second).
///
/// # Arguments
///
/// * `z` - Complex argument `z`, `z != (0.0, 0.0)`, `-PI < z.arg() <= PI`.
/// * `order` - Order of the initial H function, `order >= 0.0`.
/// * `scaling` - A parameter to indicate the scaling option.
///     * `Scaling::Unscaled`: returns `y(j) = H(order + j - 1, z)`.
///     * `Scaling::Scaled`: returns `y(j) = H(m, order + j - 1, z) * (-i * z * (3 - 2*m)).exp()`
///       where `m` is determined by the kind of Hankel function (1 for First, 2 for second).
/// * `hankel_kind` - Kind of Hankel function.
/// * `n` - Number of members in the sequence, `n >= 1`.
///
/// # Returns
///
/// A tuple containing:
/// * `y`: A vector of complex numbers containing the values of the Hankel
///   functions for orders `[order, order + 1, ..., order + n - 1]`.
/// * `n_zeros`: The number of components in `y` set to zero due to underflow.
pub fn complex_bessel_h<T: BesselFloat>(
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

    for element in y.iter_mut() {
        let scaling = if element.linf_norm() < mc.absolute_approximation_limit {
            *element *= mc.rtol;
            mc.abs_error_tolerance
        } else {
            T::ONE
        };
        *element *= phase_multiplier * scaling;
        phase_multiplier *= rotation_factor;
    }
    if partial_loss_of_significance {
        Err(PartialLossOfSignificance { y, n_zeros })
    } else {
        Ok((y, n_zeros))
    }
}

/// Computes the Hankel function of the first kind H1v(z) for a complex argument.
///
/// An alternative interface to [`complex_bessel_h`]
#[inline]
pub fn complex_hankel1<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    n: usize,
) -> BesselResult<T> {
    complex_bessel_h(z, order, scaling, HankelKind::First, n)
}

/// Computes the Hankel function of the second kind H2v(z) for a complex argument.
///
/// An alternative interface to [`complex_bessel_h`]
#[inline]
pub fn complex_hankel2<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    n: usize,
) -> BesselResult<T> {
    complex_bessel_h(z, order, scaling, HankelKind::Second, n)
}

/// Computes the I-Bessel function of a complex argument.
///
/// This function computes a sequence of complex Bessel functions `y(j) = I(order + j - 1, z)`
/// for real, non-negative orders `order + j - 1` (`j = 1, ..., n`) and a complex argument `z`
/// in the cut plane `-PI < z.arg() <= PI`.
///
/// When `scaling` is `Scaling::Scaled`, this function returns the scaled functions
/// `y(j) = (-(z.re.abs())).exp() * I(order + j - 1, z)` which remove the
/// exponential growth in both the left and right half-planes for `z` to infinity.
///
/// The computation is carried out by the power series for small `z.abs()`,
/// the asymptotic expansion for large `z.abs()`,
/// the Miller algorithm normalized by the Wronskian and a Neumann
/// series for intermediate magnitudes, and the
/// uniform asymptotic expansions for I(order, z) and J(order, z)
/// for large orders. Backward recurrence is used to generate
/// sequences or reduce orders when necessary.
///
/// The calculations above are done in the right half-plane and
/// continued into the left half-plane by the formula
/// `I(order, z * (m * PI).exp()) = (m * PI * order).exp() * I(order, z),   z.re > 0.0`
/// with `m = +i OR -i`,  (`i` is the imaginary unit).
//
///
/// # Arguments
///
/// * `z` - Complex argument `z`, `-PI < z.arg() <= PI`.
/// * `order` - Order of the initial I function, `order >= 0.0`.
/// * `scaling` - A parameter to indicate the scaling option.
///     * `Scaling::Unscaled`: returns `y(j) = I(order + j - 1, z)`.
///     * `Scaling::Scaled`: returns `y(j) = I(order + j - 1, z) * (-z.re().abs()).exp()`.
/// * `n` - Number of members of the sequence, `n >= 1`.
///
/// # Returns
///
/// A tuple containing:
/// * `y`: A vector of complex numbers containing the values of the Bessel
///   functions for orders `[order, order + 1, ..., order + n - 1]`.
/// * `n_zeros`: The number of components in `y` set to zero due to underflow.
pub fn complex_bessel_i<T: BesselFloat>(
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
            let correction = if yi.linf_norm() <= mc.absolute_approximation_limit {
                *yi *= mc.rtol;
                mc.abs_error_tolerance
            } else {
                T::ONE
            };
            *yi *= continuation_phase;
            *yi *= correction;
            continuation_phase = -continuation_phase;
        }
    }

    if partial_significance_loss {
        Err(PartialLossOfSignificance { y, n_zeros })
    } else {
        Ok((y, n_zeros))
    }
}

/// Computes the J-Bessel function of a complex argument.
///
/// This function computes a sequence of complex Bessel functions `y(j) = J(order + j - 1, z)`
/// for real, non-negative orders `order + j - 1` (`j = 1, ..., n`) and a complex argument `z`
/// in the cut plane `-PI < z.arg() <= PI`.
///
/// When `scaling` is `Scaling::Scaled`, this function returns the scaled functions
/// `y(j) = (-(z.im.abs())).exp() * J(order + j - 1, z)`, which removes the
/// exponential growth in both the upper and lower half-planes for `z` to infinity.
///
/// The computation is carried out by the formula
///
/// `J(order, Z) = ( order * PI * i / 2).exp() * I(order, -i*z)`    if `z.im >= 0.0`
///
/// `J(order, Z) = (-order * PI * i / 2).exp() * I(order, i*z)`    if `z.im < 0.0`
///
/// where `i` is the imaginary unit and `I(order, z)` is the I Bessel function.
///
/// # Arguments
///
/// * `z` - Complex argument `z`, `-PI < z.arg() <= PI`.
/// * `order` - Order of the initial J function, `order >= 0.0`.
/// * `scaling` - A parameter to indicate the scaling option.
///     * `Scaling::Unscaled`: returns `y(j) = J(order + j - 1, z)`.
///     * `Scaling::Scaled`: returns `y(j) = J(order + j - 1, z) * (-(z.im.abs())).exp()`.
/// * `n` - Number of members of the sequence, `n >= 1`.
///
/// # Returns
///
/// A tuple containing:
/// * `y`: A vector of complex numbers containing the values of the Bessel
///   functions for orders `[order, order + 1, ..., order + n - 1]`.
/// * `n_zeros`: The number of components in `y` set to zero due to underflow.
pub fn complex_bessel_j<T: BesselFloat>(
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
        let mut scaling = T::ONE;
        // TODO is the below a pattern?
        if yi.linf_norm() <= mc.absolute_approximation_limit {
            *yi *= mc.rtol;
            scaling = mc.abs_error_tolerance;
        }
        *yi *= phase_multiplier * scaling;
        phase_multiplier *= T::I * sign_selector;
    }
    if partial_significance_loss {
        Err(PartialLossOfSignificance { y, n_zeros })
    } else {
        Ok((y, n_zeros))
    }
}

/// Computes the K-Bessel function of a complex argument.
///
/// This function computes a sequence of complex Bessel functions `y(j) = K(order + j - 1, z)`
/// for real, non-negative orders `order + j - 1` (`j = 1, ..., n`) and a complex argument `z`
/// which is not equal to `(0.0, 0.0)`. The computation is valid in the cut plane
/// `-PI < z.arg() <= PI`.
///
/// When `scaling` is `Scaling::Scaled`, this function returns the scaled K functions,
/// `y(j) = z.exp() * K(order + j - 1, z)`, which remove the exponential behavior in both
/// the left and right half-planes for `z` to infinity.
///
/// EQUATIONS ARE IMPLEMENTED FOR SMALL ORDERS
/// order AND order + 1.0 IN THE RIGHT HALF PLANE X >= 0.0. FORWARD
/// RECURRENCE GENERATES HIGHER ORDERS. K IS CONTINUED TO THE LEFT
/// HALF PLANE BY THE RELATION
///
/// `K(order, z * mp.exp()) = (-mp * order).exp() * K(order, z) - mp * I(order, z)`
///
/// where `mp = mr * PI * i`, `mr = +1 OR -1`, `z.re > 0`, `i` is the imaginary unit
/// and `I(order, Z)` is the I Bessel function.
///
/// For large order, `order > MACHINE_CONSTANTS.asymptotic_order_limit`, the K function is computed
/// by means of its uniform asymptotic expansions.
///
/// # Arguments
///
/// * `z` - Complex argument `z`, `z != (0.0, 0.0)`, `-PI < z.arg() <= PI`.
/// * `order` - Order of the initial K function, `order >= 0.0`.
/// * `scaling` - A parameter to indicate the scaling option.
///     * `Scaling::Unscaled`: returns `y(j) = K(order + j - 1, z)`.
///     * `Scaling::Scaled`: returns `y(j) = K(order + j - 1, z) * z.exp()`.
/// * `n` - Number of members of the sequence, `n >= 1`.
///
/// # Returns
///
/// A tuple containing:
/// * `y`: A vector of complex numbers containing the values of the Bessel
///   functions for orders `[order, order + 1, ..., order + n - 1]`.
/// * `n_zeros`: The number of components in `y` set to zero due to underflow.
pub fn complex_bessel_k<T: BesselFloat>(
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

/// Computes the Y-Bessel function of a complex argument.
///
/// This function computes a sequence of complex Bessel functions `cy(j) = Y(order + j - 1, z)`
/// for real, non-negative orders `order + j - 1` (`j = 1, ..., n`) and a complex argument `z`
/// which is not equal to `(0.0, 0.0)`. The computation is valid in the cut plane
/// `-PI < z.arg() <= PI`.
///
/// When `scaling` is `Scaling::Scaled`, this function returns the scaled functions
/// `cy(j) = (-(z.im.abs())).exp() * Y(order + j - 1, z)`, which remove the
/// exponential growth in both the upper and lower half-planes for `z` to infinity.
///
/// The computation is carried out in terms of the I(order, z) and
/// K(order, z) Bessel functions in the right half-plane by
///
/// `Y(order, z) = i * cc * I(order, arg) - (2/PI) * cc.conj() * K(order, arg)` if `z.im >= 0`
///
/// `Y(order, z) = Y(order, z.conj()).conj()` if `z.im < 0`
///
/// where
/// `cc = (i* PI * order / 2).exp()`, `arg = z * (-i * PI / 2).exp()` and `i` is the imaginary unit.
///
/// # Arguments
///
/// * `z` - Complex argument `z`, `z != (0.0, 0.0)`, `-PI < z.arg() <= PI`.
/// * `order` - Order of the initial Y function, `order >= 0.0`.
/// * `scaling` - A parameter to indicate the scaling option.
///     * `Scaling::Unscaled`: returns `cy(j) = Y(order + j - 1, z)`.
///     * `Scaling::Scaled`: returns `cy(j) = Y(order + j - 1, z) * (-(z.im.abs())).exp()`.
/// * `n` - Number of members of the sequence, `n >= 1`.
///
/// # Returns
///
/// A tuple containing:
/// * `cy`: A vector of complex numbers containing the values of the Bessel
///   functions for orders `[order, order + 1, ..., order + n - 1]`.
/// * `n_zeros`: The number of components in `cy` set to zero due to underflow.
pub fn complex_bessel_y<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    n: usize,
) -> BesselResult<T> {
    sanitise_inputs(z, order, n, true)?;
    let mc: &MachineConsts<T> = T::MACHINE_CONSTANTS;
    let zz = if z.im < T::zero() { z.conj() } else { z };
    let zn = -T::I * zz;
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

    let (bess_i, n_zeros_i) = unwrap_psl(complex_bessel_i(zn, order, scaling, n))?;
    let (bess_k, n_zeros_k) = unwrap_psl(complex_bessel_k(zn, order, scaling, n))?;

    let mut n_zeros = n_zeros_i.min(n_zeros_k);
    let frac_order = order.fract();
    let integer_order = order.to_usize().unwrap();
    let mut i_coeff = Complex::<T>::cis(T::FRAC_PI_2() * frac_order);
    i_coeff *= i_pow_n(integer_order);
    let mut k_coeff = i_coeff.conj() * T::FRAC_2_PI();
    i_coeff *= T::I;

    let mut ey = T::one();
    if scaling == Scaling::Scaled {
        let ex = Complex::<T>::cis(z.re);
        let two_abs_z = T::two() * z.im.abs();
        ey = if two_abs_z < mc.exponent_limit {
            (-two_abs_z).exp()
        } else {
            T::zero()
        };
        k_coeff *= ex * ey;
        n_zeros = 0;
    }
    let mut y: Vec<Complex<T>> = bess_i
        .iter()
        .zip(bess_k)
        .map(|(&z_i, z_k)| {
            //----------------------------------------------------------------------;
            //       cy(I) = CSGN*cy(I)-CSPN*CWRK(I): PRODUCTS ARE COMPUTED IN;
            //       SCALED MODE if cy(I) OR CWRK(I) ARE CLOSE TO UNDERFLOW TO;
            //       PREVENT UNDERFLOW IN AN INTERMEDIATE COMPUTATION.;
            //----------------------------------------------------------------------;
            let z_k = scaled_multiply(z_k, k_coeff, scaling, mc);
            let z_i = scaled_multiply(z_i, i_coeff, scaling, mc);
            let val = z_i - z_k;
            if scaling == Scaling::Scaled && val == T::C_ZERO && ey == T::zero() {
                n_zeros += 1;
            }
            i_coeff *= T::I;
            k_coeff *= -T::I;
            val
        })
        .collect();

    if z.im < T::zero() {
        y.iter_mut().for_each(|v| *v = v.conj());
    }
    if partial_loss_of_significance {
        Err(PartialLossOfSignificance { y, n_zeros })
    } else {
        Ok((y, n_zeros))
    }
}

fn scaled_multiply<T: BesselFloat>(
    mut z: Complex<T>,
    coeff: Complex<T>,
    scaling: Scaling,
    mc: &MachineConsts<T>,
) -> Complex<T> {
    match scaling {
        Scaling::Unscaled => z * coeff,
        Scaling::Scaled => {
            let atol = if z.linf_norm() <= mc.absolute_approximation_limit {
                z *= mc.rtol;
                mc.abs_error_tolerance
            } else {
                T::one()
            };
            (z * coeff) * atol
        }
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
    const COEFF: f64 = 1.837_762_984_739_306_8e-1;
    let mc: &MachineConsts<T> = T::MACHINE_CONSTANTS;

    let abs_z = z.abs();
    //--------------------------------------------------------------------------
    //     TEST FOR PROPER RANGE
    //-----------------------------------------------------------------------
    // significance loss only tested against z, not order, so 0.0 is used to never cause significance loss
    let partial_loss_of_significance = is_significance_lost(abs_z, T::zero(), true, mc)?;

    let retval = if abs_z <= T::one() {
        //-----------------------------------------------------------------------
        //     POWER SERIES FOR z.abs() <= 1.
        //-----------------------------------------------------------------------
        let ai = airy_power_series(z, return_derivative, POWER_SERIES_COEFFS);
        (
            match scaling {
                Scaling::Scaled => ai * (T::TWO_THIRDS * z * z.sqrt()).exp(),
                Scaling::Unscaled => ai,
            },
            0,
        )
    } else {
        //-----------------------------------------------------------------------
        //     CASE FOR CABS(z) > 1.0
        //-----------------------------------------------------------------------
        let order = (if return_derivative {
            T::two()
        } else {
            T::one()
        }) / T::from_f64(3.0);
        let ln_abs_z = abs_z.ln();

        let sqrt_z = z.sqrt();
        let mut zeta = T::TWO_THIRDS * z * sqrt_z;
        //-----------------------------------------------------------------------
        //     RE(zeta) <= 0 WHEN RE(z) < 0, ESPECIALLY WHEN IM(z) IS SMALL
        //-----------------------------------------------------------------------
        let mut scale_factor = T::one();
        if z.re < T::zero() {
            zeta.re = -zeta.re.abs();
        }
        if z.im == T::zero() && z.re <= T::zero() {
            zeta.re = T::zero();
        }
        let re_zeta = zeta.re;
        let (cy, n_zeros) = if re_zeta < T::zero() || z.re <= T::zero() {
            //-----------------------------------------------------------------------
            //     OVERFLOW TEST
            //-----------------------------------------------------------------------
            if scaling == Scaling::Unscaled && re_zeta <= -mc.approximation_limit {
                scale_factor = mc.abs_error_tolerance;
                if (-re_zeta + T::from_f64(0.25) * ln_abs_z) > mc.exponent_limit {
                    return Err(Overflow);
                }
            }
            //-----------------------------------------------------------------------
            //     CBKNU AND CACON RETURN EXP(zeta)*K(order,zeta) ON KODE=2
            //-----------------------------------------------------------------------
            let rotation = if z.im < T::zero() {
                RotationDirection::Left
            } else {
                RotationDirection::Right
            };
            airy_analytic_continuation(zeta, order, scaling, rotation)?
        } else {
            //-----------------------------------------------------------------------
            //     UNDERFLOW TEST
            //-----------------------------------------------------------------------
            let mut retval = None;
            if scaling == Scaling::Unscaled && re_zeta > mc.approximation_limit {
                scale_factor = T::one() / mc.abs_error_tolerance;
                if (-re_zeta - T::from_f64(0.25) * ln_abs_z) < -mc.exponent_limit {
                    retval = Some(Ok((T::c_zeros(1), 1)));
                }
            }
            retval.unwrap_or_else(|| k_right_half_plane(zeta, order, scaling, 1))?
        };

        let mut s1 = cy[0] * T::from_f64(COEFF) * scale_factor;
        s1 *= if return_derivative { -z } else { sqrt_z };
        (s1 / scale_factor, n_zeros)
    };
    if partial_loss_of_significance {
        Err(PartialLossOfSignificance {
            y: vec![retval.0],
            n_zeros: retval.1,
        })
    } else {
        Ok(retval)
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
    const COEF: f64 = 5.773_502_691_896_257e-1;
    let mc: &MachineConsts<T> = T::MACHINE_CONSTANTS;

    let (order1, order2) = if return_derivative {
        (T::two() / T::from_f64(3.0), T::one() / T::from_f64(3.0))
    } else {
        (T::one() / T::from_f64(3.0), T::two() / T::from_f64(3.0))
    };

    let abs_z = z.abs();
    let mut partial_loss_of_significance = false;

    let bi = if abs_z <= T::one() {
        //-----------------------------------------------------------------------
        //     POWER SERIES FOR CABS(z) <= 1.
        //-----------------------------------------------------------------------
        let bi = airy_power_series(z, return_derivative, POWER_SERIES_COEFFS);
        match scaling {
            Scaling::Scaled => {
                //TODO zeta used many places with similar definition
                let zeta = T::TWO_THIRDS * (z * z.sqrt());
                bi * (-(zeta.re.abs())).exp()
            }
            Scaling::Unscaled => bi,
        }
    } else {
        //-----------------------------------------------------------------------;
        //     CASE FOR CABS(z) > 1.0;
        //-----------------------------------------------------------------------;
        //-----------------------------------------------------------------------;
        //     TEST FOR RANGE;
        //-----------------------------------------------------------------------;
        // significance loss only tested against z, not order, so 0.0 is used to never cause significance loss
        partial_loss_of_significance = is_significance_lost(abs_z, T::zero(), true, mc)?;
        let mut scale_factor = T::one();
        let mut zeta = T::TWO_THIRDS * (z * z.sqrt());

        //-----------------------------------------------------------------------;
        //     RE(zeta) <= 0 WHEN RE(z) < 0, ESPECIALLY WHEN IM(z) IS SMALL;
        //-----------------------------------------------------------------------;
        if z.re < T::zero() {
            zeta.re = -zeta.re.abs();
        }
        if z.im == T::zero() && z.re < T::zero() {
            zeta.re = T::zero();
        }
        if scaling == Scaling::Unscaled {
            //-----------------------------------------------------------------------;
            //     OVERFLOW TEST;
            //-----------------------------------------------------------------------;
            let re_zeta = zeta.re.abs();
            if re_zeta > mc.approximation_limit {
                scale_factor = mc.abs_error_tolerance;
                if re_zeta + T::from_f64(0.25) * abs_z.ln() > mc.exponent_limit {
                    return Err(Overflow);
                }
            }
        }
        let mut rotation_angle = T::zero();
        if zeta.re < T::zero() || z.re <= T::zero() {
            rotation_angle = T::PI();
            if z.im < T::zero() {
                rotation_angle = -T::PI();
            }
            zeta *= -T::one();
        }
        //-----------------------------------------------------------------------;
        //     AA=FACTOR FOR ANALYTIC CONTINUATION OF I(order,zeta);
        //     KODE=2 RETURNS EXP(-ABS(Xzeta))*I(order,zeta) FROM ZBESI;
        //-----------------------------------------------------------------------;
        let (cy, _) = i_right_half_plane(zeta, order1, scaling, 1)?;
        let mut s1 = Complex::<T>::cis(rotation_angle * order1) * cy[0] * scale_factor;
        let (mut cy, _) = i_right_half_plane(zeta, order2, scaling, 2)?;
        cy[0] *= scale_factor;
        cy[1] *= scale_factor;

        //-----------------------------------------------------------------------;
        //     BACKWARD RECUR ONE STEP FOR ORDERS -1/3 OR -2/3;
        //-----------------------------------------------------------------------;
        let s2 = (T::two() * order2) * (cy[0] / zeta) + cy[1];
        s1 =
            T::from_f64(COEF) * (s1 + s2 * Complex::<T>::cis(rotation_angle * (order2 - T::one())));
        let z_factor = if return_derivative { z } else { z.sqrt() };
        s1 * z_factor / scale_factor
    };
    if partial_loss_of_significance {
        Err(PartialLossOfSignificance {
            y: vec![bi],
            n_zeros: 0,
        })
    } else {
        Ok(bi)
    }
}
