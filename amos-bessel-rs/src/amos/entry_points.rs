use std::f64::consts::{FRAC_2_PI, FRAC_PI_2, PI};

use num::complex::{Complex64, ComplexFloat};

use crate::{
    Scaling,
    amos::{
        CIP, HankelKind, IKType, MACHINE_CONSTANTS, RotationDirection, c_one, c_zero, c_zeros,
        max_abs_component,
        overflow_checks::check_underflow_uniform_asymp_params,
        translator::{
            ZACAI, ZBUNK, airy_power_series, analytic_continuation, i_right_half_plane,
            k_right_half_plane,
        },
        utils::{TWO_THIRDS, is_significance_lost, sanitise_inputs},
    },
    types::{BesselError::*, BesselResult},
};

/// Computes the H-Bessel functions (Hankel functions) of a complex argument.
///
/// This function computes a sequence of complex Hankel (Bessel) functions
/// `cy[j] = H(order + j - 1, z)` real, non-negative
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
/// `cy(j) = (-(3 - 2 * m)*z*i).exp() * H(order + j - 1, z)` where `m` depends
/// on the kind of Hankel function (1 for First, 2 for second).
///
/// # Arguments
///
/// * `z` - Complex argument `z`, `z != (0.0, 0.0)`, `-PI < z.arg() <= PI`.
/// * `order` - Order of the initial H function, `order >= 0.0`.
/// * `scaling` - A parameter to indicate the scaling option.
///     * `Scaling::Unscaled`: returns `cy(j) = H(order + j - 1, z)`.
///     * `Scaling::Scaled`: returns `cy(j) = H(m, order + j - 1, z) * (-i * z * (3 - 2*m)).exp()`
///       where `m` is determined by the kind of Hankel function (1 for First, 2 for second).
/// * `hankel_kind` - Kind of Hankel function.
/// * `n` - Number of members in the sequence, `n >= 1`.
///
/// # Returns
///
/// A tuple containing:
/// * `cy`: A vector of complex numbers containing the values of the Hankel
///   functions for orders `[order, order + 1, ..., order + n - 1]`.
/// * `nz`: The number of components in `cy` set to zero due to underflow.
pub fn complex_bessel_h(
    z: Complex64,
    order: f64,
    scaling: Scaling,
    hankel_kind: HankelKind,
    n: usize,
) -> BesselResult {
    sanitise_inputs(z, order, n, true)?;
    let mut nz = 0;

    let modified_order = order + ((n - 1) as f64);

    let rotation = hankel_kind.get_rotation();
    let rotation_f64: f64 = rotation.into();
    let mut zn = -Complex64::I * rotation_f64 * z;
    //-----------------------------------------------------------------------
    //     TEST FOR PROPER RANGE
    //-----------------------------------------------------------------------
    let abs_z = z.abs();
    let partial_loss_of_significance = is_significance_lost(abs_z, modified_order, false)?;
    //-----------------------------------------------------------------------
    //     OVERFLOW TEST ON THE LAST MEMBER OF THE SEQUENCE
    //-----------------------------------------------------------------------
    if abs_z < MACHINE_CONSTANTS.underflow_limit {
        return Err(Overflow);
    }
    let (mut cy, nz) = if order < MACHINE_CONSTANTS.asymptotic_order_limit {
        if modified_order > 1.0 {
            if modified_order > 2.0 {
                let mut cy = c_zeros(n);
                let n_underflow = check_underflow_uniform_asymp_params(
                    zn,
                    order,
                    scaling,
                    IKType::K,
                    n,
                    &mut cy,
                )?;

                nz += n_underflow;

                // Here nn=n or nn=0 since n_underflow=(0 or nn) on return from
                // check_underflow_uniform_asymp_params (for ik_type = k)
                //
                // if nuf=nn, then cy[i]=c_zero() for all i
                if n == n_underflow {
                    return if zn.re < 0.0 {
                        Err(Overflow)
                    } else if partial_loss_of_significance {
                        Err(PartialLossOfSignificance { y: cy, nz })
                    } else {
                        Ok((cy, nz))
                    };
                }
            }
            if abs_z <= MACHINE_CONSTANTS.abs_error_tolerance
                && -modified_order * (0.5 * abs_z).ln() > MACHINE_CONSTANTS.exponent_limit
            {
                return Err(Overflow);
            }
        }
        if !((zn.re < 0.0) || (zn.re == 0.0 && zn.im < 0.0 && hankel_kind == HankelKind::Second)) {
            //-----------------------------------------------------------------------
            //     RIGHT HALF PLANE COMPUTATION, XN >= 0. && (XN != 0. ||
            //     YN >= 0. || M=1)
            //-----------------------------------------------------------------------
            k_right_half_plane(zn, order, scaling, n)?
        } else {
            //-----------------------------------------------------------------------
            //     LEFT HALF PLANE COMPUTATION
            //-----------------------------------------------------------------------
            analytic_continuation(zn, order, scaling, -rotation, n)?
        }
    } else {
        //-----------------------------------------------------------------------
        //     UNIFORM ASYMPTOTIC EXPANSIONS FOR order > asymptotic_order_limit
        //-----------------------------------------------------------------------
        let mut asymptotic_rotation = RotationDirection::None;
        if !((zn.re >= 0.0) && (zn.re != 0.0 || zn.im >= 0.0 || hankel_kind != HankelKind::Second))
        {
            asymptotic_rotation = -rotation;
            if !(zn.re != 0.0 || zn.im >= 0.0) {
                zn = -zn;
            }
        }
        let (cy, nw) = ZBUNK(zn, order, scaling, asymptotic_rotation, n)?;
        nz += nw;
        (cy, nz)
    };
    //-----------------------------------------------------------------------
    //     H(M,order,z) = -FMM*(I/FRAC_PI_2)*(ZT**order)*K(order,-z*ZT)
    //
    //     ZT=(-FMM*FRAC_PI_2*I).exp() = CMPLX(0.0,-FMM), FMM=3-2*M, M=1,2
    //-----------------------------------------------------------------------
    let sign = -FRAC_PI_2 * rotation.signum();
    //-----------------------------------------------------------------------
    //     CALCULATE (order*FRAC_PI_2*I).exp() TO MINIMIZE LOSSES OF SIGNIFICANCE
    //     WHEN order IS LARGE
    //-----------------------------------------------------------------------
    let int_order = order as i64;
    let half_int_order = int_order / 2;
    let int_remain = int_order - 2 * half_int_order;
    let arg = (order - ((int_order - int_remain) as f64)) * sign;
    let mut csgn = (1.0 / sign) * Complex64::I * Complex64::cis(arg);
    if half_int_order % 2 != 0 {
        csgn = -csgn;
    }
    for element in cy.iter_mut() {
        let scaling =
            if max_abs_component(*element) < MACHINE_CONSTANTS.absolute_approximation_limit {
                *element *= MACHINE_CONSTANTS.rtol;
                MACHINE_CONSTANTS.abs_error_tolerance
            } else {
                1.0
            };
        *element *= csgn * scaling;
        csgn *= Complex64::I * -rotation_f64;
    }
    if partial_loss_of_significance {
        Err(PartialLossOfSignificance { y: cy, nz })
    } else {
        Ok((cy, nz))
    }
}

/// Computes the Hankel function of the first kind H1v(z) for a complex argument.
///
/// An alternative interface to [`complex_bessel_h`]
#[inline]
pub fn complex_hankel1(z: Complex64, order: f64, scaling: Scaling, n: usize) -> BesselResult {
    complex_bessel_h(z, order, scaling, HankelKind::First, n)
}

/// Computes the Hankel function of the second kind H2v(z) for a complex argument.
///
/// An alternative interface to [`complex_bessel_h`]
#[inline]
pub fn complex_hankel2(z: Complex64, order: f64, scaling: Scaling, n: usize) -> BesselResult {
    complex_bessel_h(z, order, scaling, HankelKind::Second, n)
}

/// Computes the I-Bessel function of a complex argument.
///
/// This function computes a sequence of complex Bessel functions `cy(j) = I(order + j - 1, z)`
/// for real, non-negative orders `order + j - 1` (`j = 1, ..., n`) and a complex argument `z`
/// in the cut plane `-PI < z.arg() <= PI`.
///
/// When `scaling` is `Scaling::Scaled`, this function returns the scaled functions
/// `cy(j) = (-(z.re.abs())).exp() * I(order + j - 1, z)` which remove the
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
///     * `Scaling::Unscaled`: returns `cy(j) = I(order + j - 1, z)`.
///     * `Scaling::Scaled`: returns `cy(j) = I(order + j - 1, z) * (-z.re().abs()).exp()`.
/// * `n` - Number of members of the sequence, `n >= 1`.
///
/// # Returns
///
/// A tuple containing:
/// * `cy`: A vector of complex numbers containing the values of the Bessel
///   functions for orders `[order, order + 1, ..., order + n - 1]`.
/// * `nz`: The number of components in `cy` set to zero due to underflow.
pub fn complex_bessel_i(z: Complex64, order: f64, scaling: Scaling, n: usize) -> BesselResult {
    sanitise_inputs(z, order, n, false)?;

    let abs_z = z.abs();
    let modified_order = order + ((n - 1) as f64);
    let partial_significance_loss = is_significance_lost(abs_z, modified_order, false)?;

    let (zn, mut csgn) = if z.re >= 0.0 {
        (z, c_one())
    } else {
        //-----------------------------------------------------------------------
        //     CALCULATE CSGN=(order*PI*I).exp() TO MINIMIZE LOSSES OF SIGNIFICANCE
        //     WHEN order IS LARGE
        //-----------------------------------------------------------------------
        let integer_order = order as usize;
        let arg = order.fract() * PI * if z.im < 0.0 { -1.0 } else { 1.0 };
        let mut csgn = Complex64::cis(arg);
        if !integer_order.is_multiple_of(2) {
            csgn = -csgn;
        }
        (-z, csgn)
    };
    let (mut y, nz) = i_right_half_plane(zn, order, scaling, n)?;
    let remaining_n = n - nz;
    if z.re < 0.0 && remaining_n > 0 {
        //-----------------------------------------------------------------------
        //     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE
        //-----------------------------------------------------------------------
        for yi in y.iter_mut().take(remaining_n) {
            let correction =
                if max_abs_component(*yi) <= MACHINE_CONSTANTS.absolute_approximation_limit {
                    *yi *= MACHINE_CONSTANTS.rtol;
                    MACHINE_CONSTANTS.abs_error_tolerance
                } else {
                    1.0
                };
            *yi *= csgn;
            *yi *= correction;
            csgn = -csgn;
        }
    }

    if partial_significance_loss {
        Err(PartialLossOfSignificance { y, nz })
    } else {
        Ok((y, nz))
    }
}

/// Computes the J-Bessel function of a complex argument.
///
/// This function computes a sequence of complex Bessel functions `cy(j) = J(order + j - 1, z)`
/// for real, non-negative orders `order + j - 1` (`j = 1, ..., n`) and a complex argument `z`
/// in the cut plane `-PI < z.arg() <= PI`.
///
/// When `scaling` is `Scaling::Scaled`, this function returns the scaled functions
/// `cy(j) = (-(z.im.abs())).exp() * J(order + j - 1, z)`, which removes the
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
///     * `Scaling::Unscaled`: returns `cy(j) = J(order + j - 1, z)`.
///     * `Scaling::Scaled`: returns `cy(j) = J(order + j - 1, z) * (-(z.im.abs())).exp()`.
/// * `n` - Number of members of the sequence, `n >= 1`.
///
/// # Returns
///
/// A tuple containing:
/// * `cy`: A vector of complex numbers containing the values of the Bessel
///   functions for orders `[order, order + 1, ..., order + n - 1]`.
/// * `nz`: The number of components in `cy` set to zero due to underflow.
pub fn complex_bessel_j(z: Complex64, order: f64, scaling: Scaling, n: usize) -> BesselResult {
    sanitise_inputs(z, order, n, false)?;

    let partial_significance_loss = is_significance_lost(z.abs(), order + ((n - 1) as f64), false)?;
    //-----------------------------------------------------------------------
    //     CALCULATE CSGN=EXP(order*FRAC_PI_2*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
    //     WHEN order IS LARGE
    //-----------------------------------------------------------------------
    let order_int = order as i64;
    let half_order_int = order_int / 2;
    let order_rounded_down_to_even = 2 * half_order_int;
    let arg = (order - (order_rounded_down_to_even as f64)) * FRAC_PI_2;
    let mut csgn = Complex64::cis(arg);
    if (half_order_int % 2) != 0 {
        csgn = -csgn;
    }
    //-----------------------------------------------------------------------
    //     ZN IS IN THE RIGHT HALF PLANE
    //-----------------------------------------------------------------------
    let mut sign_selector = 1.0;
    let mut zn = -Complex64::I * z;
    if z.im < 0.0 {
        zn = -zn;
        csgn.im = -csgn.im;
        sign_selector = -sign_selector;
    }
    let (mut cy, nz) = i_right_half_plane(zn, order, scaling, n)?;
    for cyi in cy.iter_mut().take(n - nz) {
        let mut scaling = 1.0;
        // TODO is the below a pattern?
        if (max_abs_component(*cyi)) <= MACHINE_CONSTANTS.absolute_approximation_limit {
            *cyi *= MACHINE_CONSTANTS.rtol;
            scaling = MACHINE_CONSTANTS.abs_error_tolerance;
        }
        *cyi *= csgn * scaling;
        csgn *= sign_selector * Complex64::I;
    }
    if partial_significance_loss {
        Err(PartialLossOfSignificance { y: cy, nz })
    } else {
        Ok((cy, nz))
    }
}

/// Computes the K-Bessel function of a complex argument.
///
/// This function computes a sequence of complex Bessel functions `cy(j) = K(order + j - 1, z)`
/// for real, non-negative orders `order + j - 1` (`j = 1, ..., n`) and a complex argument `z`
/// which is not equal to `(0.0, 0.0)`. The computation is valid in the cut plane
/// `-PI < z.arg() <= PI`.
///
/// When `scaling` is `Scaling::Scaled`, this function returns the scaled K functions,
/// `cy(j) = z.exp() * K(order + j - 1, z)`, which remove the exponential behavior in both
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
///     * `Scaling::Unscaled`: returns `cy(j) = K(order + j - 1, z)`.
///     * `Scaling::Scaled`: returns `cy(j) = K(order + j - 1, z) * z.exp()`.
/// * `n` - Number of members of the sequence, `n >= 1`.
///
/// # Returns
///
/// A tuple containing:
/// * `cy`: A vector of complex numbers containing the values of the Bessel
///   functions for orders `[order, order + 1, ..., order + n - 1]`.
/// * `nz`: The number of components in `cy` set to zero due to underflow.
pub fn complex_bessel_k(z: Complex64, order: f64, scaling: Scaling, n: usize) -> BesselResult {
    sanitise_inputs(z, order, n, true)?;
    //-----------------------------------------------------------------------------;
    //     TEST FOR PROPER RANGE;
    //-----------------------------------------------------------------------;
    let abs_z = z.abs();
    let modified_order = order + ((n - 1) as f64);
    let partial_significance_loss = is_significance_lost(abs_z, modified_order, false)?;

    //-----------------------------------------------------------------------;
    //     OVERFLOW TEST ON THE LAST MEMBER OF THE SEQUENCE;
    //-----------------------------------------------------------------------;
    if abs_z < MACHINE_CONSTANTS.underflow_limit {
        return Err(Overflow);
    }

    let mut nz = 0;
    if order > MACHINE_CONSTANTS.asymptotic_order_limit {
        //-----------------------------------------------------------------------
        //     UNIFORM ASYMPTOTIC EXPANSIONS FOR order > asymptotic_order_limit
        //-----------------------------------------------------------------------
        let rotation = if z.re >= 0.0 {
            RotationDirection::None
        } else if z.im < 0.0 {
            RotationDirection::Left
        } else {
            RotationDirection::Right
        };

        let (y, nz) = ZBUNK(z, order, scaling, rotation, n)?;
        return if partial_significance_loss {
            Err(PartialLossOfSignificance { y, nz })
        } else {
            Ok((y, nz))
        };
    }

    if modified_order > 2.0 {
        let mut y = c_zeros(n);
        let n_underflow =
            check_underflow_uniform_asymp_params(z, order, scaling, IKType::K, n, &mut y)?;
        nz += n_underflow;

        //-----------------------------------------------------------------------;
        //     HERE NN=n OR NN=0 SINCE NUF=0,NN, OR -1 ON RETURN FROM CUOIK;
        //     if NUF=NN, THEN cy(I)=CZERO FOR ALL I;
        //-----------------------------------------------------------------------;
        if n_underflow == n {
            return if z.re < 0.0 {
                Err(Overflow)
            } else if partial_significance_loss {
                Err(PartialLossOfSignificance { y, nz })
            } else {
                Ok((y, nz))
            };
        }
    }
    if (modified_order > 1.0) && abs_z <= MACHINE_CONSTANTS.abs_error_tolerance {
        let half_abs_z = 0.5 * abs_z;
        if -modified_order * half_abs_z.ln() > MACHINE_CONSTANTS.exponent_limit {
            return Err(Overflow);
        }
    }
    let (y, nz) = if z.re >= 0.0 {
        //-----------------------------------------------------------------------;
        //     RIGHT HALF PLANE COMPUTATION, REAL(z) >= 0.;
        //-----------------------------------------------------------------------;
        k_right_half_plane(z, order, scaling, n)?
    } else {
        //-----------------------------------------------------------------------;
        //     LEFT HALF PLANE COMPUTATION;
        //     PI/2 < z.arg() <= PI AND -PI < z.arg() < -PI/2.;
        //-----------------------------------------------------------------------;
        if nz != 0 {
            return Err(Overflow);
        }
        let rotation = if z.im < 0.0 {
            RotationDirection::Left
        } else {
            RotationDirection::Right
        };
        analytic_continuation(z, order, scaling, rotation, n)?
    };
    if partial_significance_loss {
        Err(PartialLossOfSignificance { y, nz })
    } else {
        Ok((y, nz))
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
/// * `nz`: The number of components in `cy` set to zero due to underflow.
pub fn complex_bessel_y(z: Complex64, order: f64, scaling: Scaling, n: usize) -> BesselResult {
    sanitise_inputs(z, order, n, true)?;
    let zz = if z.im < 0.0 { z.conj() } else { z };
    let zn = -Complex64::I * zz;
    let mut partial_loss_of_significance = false;

    let mut unwrap_psl = |result: BesselResult| match result {
        Ok((y_, nz_)) => Ok((y_, nz_)),
        Err(PartialLossOfSignificance { y: y_, nz: nz_ }) => {
            partial_loss_of_significance = true;
            Ok((y_, nz_))
        }
        err => err,
    };

    let (bess_i, nz_i) = unwrap_psl(complex_bessel_i(zn, order, scaling, n))?;
    let (bess_k, nz_k) = unwrap_psl(complex_bessel_k(zn, order, scaling, n))?;

    let mut nz = nz_i.min(nz_k);
    let frac_order = order.fract();
    let integer_order = order as usize;
    let mut csgn = Complex64::cis(FRAC_PI_2 * frac_order);
    let index = integer_order % 4;
    csgn *= CIP[index];
    let mut cspn = csgn.conj() * FRAC_2_PI;
    csgn *= Complex64::I;

    let mut ey = 1.0;
    if scaling == Scaling::Scaled {
        let ex = Complex64::cis(z.re);
        let two_abs_z = 2.0 * z.im.abs();
        ey = if two_abs_z < MACHINE_CONSTANTS.exponent_limit {
            (-two_abs_z).exp()
        } else {
            0.0
        };
        cspn *= ex * ey;
        nz = 0;
    }
    let mut y: Vec<Complex64> = bess_i
        .iter()
        .zip(bess_k)
        .map(|(&z_i, z_k)| {
            //----------------------------------------------------------------------;
            //       cy(I) = CSGN*cy(I)-CSPN*CWRK(I): PRODUCTS ARE COMPUTED IN;
            //       SCALED MODE if cy(I) OR CWRK(I) ARE CLOSE TO UNDERFLOW TO;
            //       PREVENT UNDERFLOW IN AN INTERMEDIATE COMPUTATION.;
            //----------------------------------------------------------------------;
            let z_k = scaled_multiply(z_k, cspn, scaling);
            let z_i = scaled_multiply(z_i, csgn, scaling);
            let val = z_i - z_k;
            if scaling == Scaling::Scaled && val == c_zero() && ey == 0.0 {
                nz += 1;
            }
            csgn *= Complex64::I;
            cspn *= -Complex64::I;
            val
        })
        .collect();

    if z.im < 0.0 {
        y.iter_mut().for_each(|v| *v = v.conj());
    }
    if partial_loss_of_significance {
        Err(PartialLossOfSignificance { y, nz })
    } else {
        Ok((y, nz))
    }
}

fn scaled_multiply(mut z: Complex64, coeff: Complex64, scaling: Scaling) -> Complex64 {
    match scaling {
        Scaling::Unscaled => z * coeff,
        Scaling::Scaled => {
            let atol = if max_abs_component(z) <= MACHINE_CONSTANTS.absolute_approximation_limit {
                z *= MACHINE_CONSTANTS.rtol;
                MACHINE_CONSTANTS.abs_error_tolerance
            } else {
                1.0
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
pub fn complex_airy(
    z: Complex64,
    return_derivative: bool,
    scaling: Scaling,
) -> BesselResult<(Complex64, usize)> {
    const POWER_SERIES_COEFFS: (f64, f64) = (3.550_280_538_878_172e-1, 2.588_194_037_928_068e-1);
    const COEFF: f64 = 1.837_762_984_739_306_8e-1;

    let abs_z = z.abs();
    let float_is_derivative = if return_derivative { 1.0 } else { 0.0 };
    //--------------------------------------------------------------------------
    //     TEST FOR PROPER RANGE
    //-----------------------------------------------------------------------
    // significance loss only tested against z, not order, so 0.0 is used to never cause significance loss
    let partial_loss_of_significance = is_significance_lost(abs_z, 0.0, true)?;

    let retval = if abs_z <= 1.0 {
        //-----------------------------------------------------------------------
        //     POWER SERIES FOR z.abs() <= 1.
        //-----------------------------------------------------------------------
        let ai = airy_power_series(z, return_derivative, POWER_SERIES_COEFFS);
        (
            match scaling {
                Scaling::Scaled => ai * (TWO_THIRDS * z * z.sqrt()).exp(),
                Scaling::Unscaled => ai,
            },
            0,
        )
    } else {
        //-----------------------------------------------------------------------
        //     CASE FOR CABS(z) > 1.0
        //-----------------------------------------------------------------------
        let order = (1.0 + float_is_derivative) / 3.0;
        let ln_abs_z = abs_z.ln();

        let sqrt_z = z.sqrt();
        let mut zeta = TWO_THIRDS * z * sqrt_z;
        //-----------------------------------------------------------------------
        //     RE(zeta) <= 0 WHEN RE(z) < 0, ESPECIALLY WHEN IM(z) IS SMALL
        //-----------------------------------------------------------------------
        let mut scale_factor = 1.0;
        if z.re < 0.0 {
            zeta.re = -zeta.re.abs();
        }
        if z.im == 0.0 && z.re <= 0.0 {
            zeta.re = 0.0;
        }
        let re_zeta = zeta.re;
        let (cy, nz) = if re_zeta < 0.0 || z.re <= 0.0 {
            //-----------------------------------------------------------------------
            //     OVERFLOW TEST
            //-----------------------------------------------------------------------
            if scaling == Scaling::Unscaled && re_zeta <= -MACHINE_CONSTANTS.approximation_limit {
                scale_factor = MACHINE_CONSTANTS.abs_error_tolerance;
                if (-re_zeta + 0.25 * ln_abs_z) > MACHINE_CONSTANTS.exponent_limit {
                    return Err(Overflow);
                }
            }
            //-----------------------------------------------------------------------
            //     CBKNU AND CACON RETURN EXP(zeta)*K(order,zeta) ON KODE=2
            //-----------------------------------------------------------------------
            let rotation = if z.im < 0.0 {
                RotationDirection::Left
            } else {
                RotationDirection::Right
            };
            ZACAI(zeta, order, scaling, rotation, 1)?
        } else {
            //-----------------------------------------------------------------------
            //     UNDERFLOW TEST
            //-----------------------------------------------------------------------
            let mut retval = None;
            if scaling == Scaling::Unscaled && re_zeta > MACHINE_CONSTANTS.approximation_limit {
                scale_factor = 1.0 / MACHINE_CONSTANTS.abs_error_tolerance;
                if (-re_zeta - 0.25 * ln_abs_z) < -MACHINE_CONSTANTS.exponent_limit {
                    retval = Some(Ok((c_zeros(1), 1)));
                }
            }
            retval.unwrap_or_else(|| k_right_half_plane(zeta, order, scaling, 1))?
        };

        let mut s1 = cy[0] * COEFF * scale_factor;
        s1 *= if return_derivative { -z } else { sqrt_z };
        (s1 / scale_factor, nz)
    };
    if partial_loss_of_significance {
        Err(PartialLossOfSignificance {
            y: vec![retval.0],
            nz: retval.1,
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
pub fn complex_airy_b(
    z: Complex64,
    return_derivative: bool,
    scaling: Scaling,
) -> BesselResult<Complex64> {
    const POWER_SERIES_COEFFS: (f64, f64) = (6.149_266_274_460_007e-1, -4.482_883_573_538_264e-1);
    const COEF: f64 = 5.773_502_691_896_257e-1;

    let abs_z = z.abs();
    let float_is_derivative = if return_derivative { 1.0 } else { 0.0 };
    let mut partial_loss_of_significance = false;

    let bi = if abs_z <= 1.0 {
        //-----------------------------------------------------------------------
        //     POWER SERIES FOR CABS(z) <= 1.
        //-----------------------------------------------------------------------
        let bi = airy_power_series(z, return_derivative, POWER_SERIES_COEFFS);
        match scaling {
            Scaling::Scaled => {
                //TODO zeta used many places with similar definition
                let zeta = TWO_THIRDS * (z * z.sqrt());
                bi * (-(zeta.re.abs())).exp()
            }
            Scaling::Unscaled => bi,
        }
    } else {
        //-----------------------------------------------------------------------;
        //     CASE FOR CABS(z) > 1.0;
        //-----------------------------------------------------------------------;
        let order = (1.0 + float_is_derivative) / 3.0;
        //-----------------------------------------------------------------------;
        //     TEST FOR RANGE;
        //-----------------------------------------------------------------------;
        // significance loss only tested against z, not order, so 0.0 is used to never cause significance loss
        partial_loss_of_significance = is_significance_lost(abs_z, 0.0, true)?;
        let mut scale_factor = 1.0;
        let mut zeta = TWO_THIRDS * (z * z.sqrt());

        //-----------------------------------------------------------------------;
        //     RE(zeta) <= 0 WHEN RE(z) < 0, ESPECIALLY WHEN IM(z) IS SMALL;
        //-----------------------------------------------------------------------;
        if z.re < 0.0 {
            zeta.re = -zeta.re.abs();
        }
        if z.im == 0.0 && z.re < 0.0 {
            zeta.re = 0.0;
        }
        if scaling == Scaling::Unscaled {
            //-----------------------------------------------------------------------;
            //     OVERFLOW TEST;
            //-----------------------------------------------------------------------;
            let re_zeta = zeta.re.abs();
            if re_zeta > MACHINE_CONSTANTS.approximation_limit {
                scale_factor = MACHINE_CONSTANTS.abs_error_tolerance;
                if re_zeta + 0.25 * abs_z.ln() > MACHINE_CONSTANTS.exponent_limit {
                    return Err(Overflow);
                }
            }
        }
        let mut rotation_angle = 0.0;
        if zeta.re < 0.0 || z.re <= 0.0 {
            rotation_angle = PI;
            if z.im < 0.0 {
                rotation_angle = -PI;
            }
            zeta *= -1.0;
        }
        //-----------------------------------------------------------------------;
        //     AA=FACTOR FOR ANALYTIC CONTINUATION OF I(order,zeta);
        //     KODE=2 RETURNS EXP(-ABS(Xzeta))*I(order,zeta) FROM ZBESI;
        //-----------------------------------------------------------------------;
        let (cy, _) = i_right_half_plane(zeta, order, scaling, 1)?;
        let mut s1 = Complex64::cis(rotation_angle * order) * cy[0] * scale_factor;
        let order = (2.0 - float_is_derivative) / 3.0;
        let (mut cy, _) = i_right_half_plane(zeta, order, scaling, 2)?;
        cy[0] *= scale_factor;
        cy[1] *= scale_factor;

        //-----------------------------------------------------------------------;
        //     BACKWARD RECUR ONE STEP FOR ORDERS -1/3 OR -2/3;
        //-----------------------------------------------------------------------;
        let s2 = (2.0 * order) * (cy[0] / zeta) + cy[1];
        s1 = COEF * (s1 + s2 * Complex64::cis(rotation_angle * (order - 1.0)));
        let z_factor = if return_derivative { z } else { z.sqrt() };
        s1 * z_factor / scale_factor
    };
    if partial_loss_of_significance {
        Err(PartialLossOfSignificance { y: vec![bi], nz: 0 })
    } else {
        Ok(bi)
    }
}
