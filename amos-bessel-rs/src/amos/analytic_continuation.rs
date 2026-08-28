use num::{Complex, Integer, complex::ComplexFloat};

use crate::{
    BesselError, BesselFloat, Scaling,
    amos::{
        MachineConsts, RotationDirection,
        asymptotics::i_asymptotic,
        i_computation::i_miller,
        limits::{OverflowState, underflow_add_i_k},
        power_series::i_power_series,
        right_half_plane::{i_right_half_plane, k_right_half_plane},
        utils::two_over_z_safe,
    },
    types::{BesselResult, BesselValues},
};

/// Applies the analytic continuation formula
///     K(fnu, zn*exp(mp))=K(fnu, zn)*exp(-mp*fnu) - mp*I(fnu, zn)
///             mp=pi*mr*cmplx(0.0,1.0)
///
/// to continue the k function from the right half to the left
/// half z plane
///
/// Originally ZACON
pub fn analytic_continuation<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    rotation: RotationDirection,
    n: usize,
) -> Result<BesselValues<T>, BesselError<T>> {
    let mc: &MachineConsts<T> = T::MACHINE_CONSTANTS;
    let mut n_zeros = 0;
    let negative_z = -z;
    let (i_values, _) = i_right_half_plane(negative_z, order, scaling, n)?;
    //-----------------------------------------------------------------------
    // Analytic continuation to the left half plane for the K function
    //-----------------------------------------------------------------------
    let (k_seeds, n_zeros_k) = k_right_half_plane(negative_z, order, scaling, 2.min(n))?;
    if n_zeros_k > 0 {
        return Err(BesselError::Overflow);
        // Amos also handled  n_zeros_inner = -1 or -2  as error cases, but in rust these
        // are handled by k_right_half_plane returning an error,
        // The amos code defaults to an overflow, if n_zeros_inner != 0
    }
    // The base continuation formula is K_v(z) = K_v(-z)*exp(-v*i*pi*m) - i*pi*m*I_v(-z).
    // `rotation_angle` represents `-pi*m`.
    let rotation_angle = -T::PI() * T::from_f64(rotation.signum());

    // This initializes the `-i*pi*m` coefficient for the I_v(-z) term.
    let mut i_continuation_coeff = Complex::<T>::new(T::ZERO, rotation_angle);
    if scaling == Scaling::Scaled {
        // In scaled mode, we must return K_v(z)*exp(z).
        // The I_v(-z) term from the right-half plane solver is scaled by exp(-Re(-z)) = exp(Re(z)).
        // To achieve the full exp(z) = exp(Re(z))*exp(i*Im(z)) scaling, we multiply the
        // coefficient by the missing phase factor exp(i*Im(z)) which equals exp(-i*Im(-z)).
        i_continuation_coeff *= Complex::<T>::cis(-negative_z.im);
    }
    //-----------------------------------------------------------------------
    // Calculate continuation coefficient = exp(order * pi * i) to minimize
    // losses of significance when order is large
    //-----------------------------------------------------------------------
    let mut k_continuation_coeff = Complex::<T>::cis(order.fract() * rotation_angle);
    if order.to_usize().unwrap().is_odd() {
        k_continuation_coeff = -k_continuation_coeff;
    }

    let mut k_component = k_seeds[0];
    let mut i_component = i_values[0];
    if scaling == Scaling::Scaled
        && underflow_add_i_k(negative_z, &mut k_component, &mut i_component, mc)
    {
        n_zeros += 1;
    }

    // re-use the i_values buffer for y, modifying it in place below.
    // Saves memory allocation
    let mut y = i_values;
    y[0] = k_continuation_coeff * k_component + i_continuation_coeff * i_component;
    if n == 1 {
        return Ok((y, n_zeros));
    }
    k_continuation_coeff = -k_continuation_coeff;

    let mut k_prev = k_seeds[0];
    let mut k_curr = k_seeds[1];
    k_component = k_curr;
    i_component = y[1]; // y = i_values here
    let mut scaled_k_component = if scaling == Scaling::Scaled {
        if underflow_add_i_k(negative_z, &mut k_component, &mut i_component, mc) {
            n_zeros += 1;
        }
        Some(k_component)
    } else {
        None
    };

    y[1] = k_continuation_coeff * k_component + i_continuation_coeff * i_component;
    if n == 2 {
        return Ok((y, n_zeros));
    }

    k_continuation_coeff = -k_continuation_coeff;
    let two_over_z = two_over_z_safe(negative_z);

    //-----------------------------------------------------------------------
    // Scale near exponent extremes during recurrence on K functions
    //-----------------------------------------------------------------------
    let abs_k_curr = k_curr.abs();
    let mut overflow_state = if abs_k_curr <= OverflowState::boundary(&OverflowState::NearUnder, mc)
    {
        OverflowState::NearUnder
    } else if abs_k_curr > OverflowState::boundary(&OverflowState::None, mc) {
        OverflowState::NearOver
    } else {
        OverflowState::None
    };
    let mut boundary = overflow_state.boundary::<T>(mc);
    k_prev *= overflow_state.scaling_factor::<T>(mc);
    k_curr *= overflow_state.scaling_factor::<T>(mc);
    let mut recip_scaling_factor = overflow_state.reciprocal_scaling_factor::<T>(mc);
    let mut out_of_underflow_regime = false;
    // y contains the I values. The K values are computed and combined with I,
    // then assigned back to y in place below.
    let mut n_without_underflow = 0;
    for (i, y_val) in y.iter_mut().enumerate().skip(2) {
        let recurrence_factor = (order + T::from_usize(i - 1)) * two_over_z;
        (k_prev, k_curr) = (k_curr, recurrence_factor * k_curr + k_prev);
        k_component = k_curr * recip_scaling_factor;
        let mut unscaled_k_curr = k_component;
        i_component = *y_val;
        if scaling == Scaling::Scaled && !out_of_underflow_regime {
            let full_underflow =
                underflow_add_i_k(negative_z, &mut k_component, &mut i_component, mc);
            let partial_underflow = k_component == T::C_ZERO;

            if full_underflow {
                n_zeros += 1;
            }
            if full_underflow || partial_underflow {
                n_without_underflow = 0; // either full underflow or K was squashed to zero, reset the counter
            } else {
                n_without_underflow += 1; // Only increment if K actually survived scaling!
            }

            let saved_k_component = scaled_k_component.unwrap();
            scaled_k_component = Some(k_component);
            if n_without_underflow == 3 {
                out_of_underflow_regime = true;
                k_prev = saved_k_component * overflow_state.scaling_factor::<T>(mc);
                k_curr = k_component * overflow_state.scaling_factor::<T>(mc);
                unscaled_k_curr = k_component;
            }
        }
        *y_val = k_continuation_coeff * k_component + i_continuation_coeff * i_component;
        k_continuation_coeff = -k_continuation_coeff;
        overflow_state.scale_recurrence(
            &mut k_prev,
            &mut k_curr,
            unscaled_k_curr,
            &mut boundary,
            &mut recip_scaling_factor,
            mc,
        );
    }
    Ok((y, n_zeros))
}

/// Applies the analytic continuation formula
///
/// K(fnu, zn*exp(mp)) = K(fnu, zn) * exp(-mp*fnu) - mp*I(fnu,zn)
/// mp = pi * mr * cmplx(0.0,1.0)
///
/// to continue the k function from the right half to the left
/// half z plane for use with complex_airy where fnu=1/3 or 2/3 and n=1.
///
/// airy_analytic_continuation is the same as analytic_continuation (zacon)
/// with the parts for larger orders and recurrence removed.
///
/// Amos introduced it to prevent a cyclic call graph (which early Fortran compilers forbade)
/// complex_airy  →  analytic_continuation  →  k_right_half_plane  →  complex_airy
/// This call chain never occurs at runtime, but Fortran forbade even the possibility.
///
/// The function is left in place for the Rust version as it provides an optimised option
/// for airy evaluation without the overhead of the general function.
///
/// Originally ZACAI
pub fn airy_analytic_continuation<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    rotation: RotationDirection,
) -> BesselResult<T> {
    let mc: &MachineConsts<T> = T::MACHINE_CONSTANTS;
    let mut n_zeros = 0;
    let negative_z = -z;
    let abs_z = z.abs();

    let i_value = if (abs_z * abs_z * T::from_f64(0.25) <= order + T::one()) || (abs_z <= T::TWO) {
        //-----------------------------------------------------------------------
        // Power series for the I function
        //-----------------------------------------------------------------------
        let (y, n_zeros_inner_signed) = i_power_series(negative_z, order, scaling, 1)?;
        // While some calls to i_power_series can return negative values,
        // the call here should not
        debug_assert!(n_zeros_inner_signed >= 0);
        y[0]
    } else if abs_z >= mc.asymptotic_z_limit {
        //-----------------------------------------------------------------------
        // Asymptotic expansion for large z for the I function
        //-----------------------------------------------------------------------
        let (y, _) = i_asymptotic(negative_z, order, scaling, 1)?;
        y[0]

    //-----------------------------------------------------------------------
    // Miller algorithm normalized by the series for the I function
    //-----------------------------------------------------------------------
    } else {
        let y = i_miller(negative_z, order, scaling, 1)?;
        y[0]
    };
    //-----------------------------------------------------------------------
    // Analytic continuation to the left half plane for the K function
    //-----------------------------------------------------------------------
    let (k_value, n_zeros_k) = k_right_half_plane(negative_z, order, scaling, 1)?;
    if n_zeros_k != 0 {
        return Err(BesselError::Overflow);
    }
    // The base continuation formula is K_v(z) = K_v(-z)*exp(-v*i*pi*m) - i*pi*m*I_v(-z).
    // `rotation_angle` represents `-pi*m`.
    let rotation_angle = -T::PI() * T::from_f64(rotation.signum());

    // This initializes the `-i*pi*m` coefficient for the I_v(-z) term.
    let mut i_coeff = Complex::<T>::new(T::ZERO, rotation_angle);
    if scaling == Scaling::Scaled {
        // In scaled mode, we must return K_v(z)*exp(z).
        // The I_v(-z) term from the right-half plane solver is scaled by exp(-Re(-z)) = exp(Re(z)).
        // To achieve the full exp(z) = exp(Re(z))*exp(i*Im(z)) scaling, we multiply the
        // coefficient by the missing phase factor exp(i*Im(z)) which equals exp(-i*Im(-z)).
        i_coeff *= Complex::<T>::cis(-negative_z.im);
    }
    //-----------------------------------------------------------------------
    // Calculate continuation coefficient = exp(order * pi * i) to minimize
    // losses of significance when order is large
    //-----------------------------------------------------------------------

    let mut k_coeff = Complex::<T>::cis(order.fract() * rotation_angle);
    if order.to_usize().unwrap().is_odd() {
        k_coeff = -k_coeff;
    }
    let mut k_value = k_value[0];
    let mut i_value = i_value;
    if scaling == Scaling::Scaled && underflow_add_i_k(negative_z, &mut k_value, &mut i_value, mc) {
        n_zeros += 1;
    }
    let y = vec![k_coeff * k_value + i_coeff * i_value];
    Ok((y, n_zeros))
}
