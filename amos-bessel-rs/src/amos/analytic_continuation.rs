use num::{Complex, Integer, complex::ComplexFloat};

use crate::{
    BesselError, Scaling,
    amos::{
        RotationDirection,
        asymptotics::i_asymptotic,
        i_power_series,
        limits::{OverflowState, underflow_add_i_k},
        max_abs_component,
        recurrence::i_miller,
        right_half_plane::{i_right_half_plane, k_right_half_plane},
        utils::calc_rz,
    },
    types::{BesselFloat, BesselResult, BesselValues},
};

/// Applies the analytic continuation formula
///     K(fnu, zn*exp(mp))=K(fnu, zn)*exp(-mp*fnu) - mp*I(fnu, zn)
///             mp=pi*mr*cmplx(0.0,1.0)
///
/// to continue the k function from the right half to the left
/// half z plane
/// Originally ZACON
pub fn analytic_continuation<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    rotation: RotationDirection,
    n: usize,
) -> Result<BesselValues<T>, BesselError<T>> {
    let mut n_zeros = 0;
    let negative_z = -z;
    let (i_values, _) = i_right_half_plane(negative_z, order, scaling, n)?;
    //-----------------------------------------------------------------------
    //     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
    //-----------------------------------------------------------------------
    let (k_values, n_zeros_inner) = k_right_half_plane(negative_z, order, scaling, 2.min(n))?;
    if n_zeros_inner > 0 {
        return Err(BesselError::Overflow);
        // Amos also handled  n_zeros_inner = -1 or -2  as error cases, but in rust these
        // are handled by k_right_half_plane returning an error,
        // The amos code defaults to an overflow, if n_zeros_inner != 0
    }
    let rotation_angle = -T::PI() * T::from_f64(rotation.signum());
    let mut i_continuation_coeff = Complex::<T>::new(T::ZERO, rotation_angle);
    if scaling == Scaling::Scaled {
        i_continuation_coeff *= Complex::<T>::cis(-negative_z.im);
    }
    //-----------------------------------------------------------------------
    //     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
    //     WHEN FNU IS LARGE
    //-----------------------------------------------------------------------
    let mut k_continuation_coeff = Complex::<T>::cis(order.fract() * rotation_angle);
    if order.to_usize().unwrap().is_odd() {
        k_continuation_coeff = -k_continuation_coeff;
    }
    let mut n_good = 0;
    let mut k_prev = k_values[0];
    let mut k_component = k_prev;
    let mut i_component = i_values[0];
    if scaling == Scaling::Scaled
        && underflow_add_i_k(negative_z, &mut k_component, &mut i_component, &mut n_good)
    {
        n_zeros += 1;
    }

    let mut y = T::c_zeros(n);
    y[0] = k_continuation_coeff * k_component + i_continuation_coeff * i_component;
    if n == 1 {
        return Ok((y, n_zeros));
    }

    k_continuation_coeff = -k_continuation_coeff;
    let mut k_curr = k_values[1];
    k_component = k_curr;
    i_component = i_values[1];
    let mut scaled_k_component = if scaling == Scaling::Scaled {
        if underflow_add_i_k(negative_z, &mut k_component, &mut i_component, &mut n_good) {
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
    let reciprocal_z = calc_rz(negative_z);
    let mut recurrence_factor = (order + T::one()) * reciprocal_z;
    //-----------------------------------------------------------------------
    //     SCALE NEAR EXPONENT EXTREMES DURING RECURRENCE ON K FUNCTIONS
    //-----------------------------------------------------------------------
    let abs_s2 = k_curr.abs();
    let mut overflow_state = if abs_s2 <= OverflowState::boundary(&OverflowState::NearUnder) {
        OverflowState::NearUnder
    } else if abs_s2 > OverflowState::boundary(&OverflowState::None) {
        OverflowState::NearOver
    } else {
        OverflowState::None
    };
    let mut boundary = overflow_state.boundary::<T>();
    k_prev *= overflow_state.scaling_factor::<T>();
    k_curr *= overflow_state.scaling_factor::<T>();
    let mut recip_scaling_factor = overflow_state.reciprocal_scaling_factor::<T>();
    for (yi, ii) in y.iter_mut().zip(i_values).skip(2) {
        //TODO common pattern below
        (k_prev, k_curr) = (k_curr, recurrence_factor * k_curr + k_prev);
        k_component = k_curr * recip_scaling_factor;
        let mut unscaled_k_curr = k_component;
        i_component = ii;
        if scaling == Scaling::Scaled && n_good >= 0 {
            if underflow_add_i_k(negative_z, &mut k_component, &mut i_component, &mut n_good) {
                n_zeros += 1;
            }
            let saved_k_component = scaled_k_component.unwrap();
            scaled_k_component = Some(k_component);
            if n_good == 3 {
                n_good = -4;
                k_prev = saved_k_component * overflow_state.scaling_factor::<T>();
                k_curr = k_component * overflow_state.scaling_factor::<T>();
                unscaled_k_curr = k_component;
            }
        }
        *yi = k_continuation_coeff * k_component + i_continuation_coeff * i_component;
        recurrence_factor += reciprocal_z;
        k_continuation_coeff = -k_continuation_coeff;
        if overflow_state != OverflowState::NearOver && max_abs_component(k_component) < boundary {
            overflow_state.increment();
            boundary = overflow_state.boundary::<T>();
            k_prev *= recip_scaling_factor;
            k_curr = unscaled_k_curr;
            k_prev *= overflow_state.scaling_factor::<T>();
            k_curr *= overflow_state.scaling_factor::<T>();
            recip_scaling_factor = overflow_state.reciprocal_scaling_factor::<T>();
        }
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
/// zacai is the same as analytic_continuation (zacon) with the parts for larger orders and
/// recurrence removed. A recursive call to zacon can result if zacon
/// is called from complex_airy.
///
/// Originally ZACAI
pub fn airy_analytic_continuation<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    rotation: RotationDirection,
) -> BesselResult<T> {
    let mut n_zeros = 0;
    let negative_z = -z;
    let abs_z = z.abs();

    let (i_value, _) =
        if (abs_z * abs_z * T::from_f64(0.25) <= order + T::one()) || (abs_z <= T::two()) {
            //-----------------------------------------------------------------------
            //     POWER SERIES FOR THE I FUNCTION
            //-----------------------------------------------------------------------
            let (y, n_zeros_inner_signed) = i_power_series(negative_z, order, scaling, 1)?;
            // While some calls to i_power_series can return negative values,
            // the call here should not
            debug_assert!(n_zeros_inner_signed >= 0);
            (y, n_zeros_inner_signed.unsigned_abs())
        } else if abs_z >= T::MACHINE_CONSTANTS.asymptotic_z_limit {
            //-----------------------------------------------------------------------
            //     ASYMPTOTIC EXPANSION FOR LARGE Z FOR THE I FUNCTION
            //-----------------------------------------------------------------------
            i_asymptotic(negative_z, order, scaling, 1)?
        //-----------------------------------------------------------------------
        //     MILLER ALGORITHM NORMALIZED BY THE SERIES FOR THE I FUNCTION
        //-----------------------------------------------------------------------
        } else {
            i_miller(negative_z, order, scaling, 1)?
        };
    //-----------------------------------------------------------------------
    //     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
    //-----------------------------------------------------------------------s
    let (k_value, n_zeros_k) = k_right_half_plane(negative_z, order, scaling, 1)?;
    if n_zeros_k != 0 {
        return Err(BesselError::Overflow);
    }
    let rotation_angle = -T::PI() * T::from_f64(rotation.signum());
    let mut i_coeff = Complex::<T>::new(T::ZERO, rotation_angle);
    if scaling == Scaling::Scaled {
        i_coeff *= Complex::<T>::cis(-negative_z.im);
    }
    //-----------------------------------------------------------------------
    //     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
    //     WHEN FNU IS LARGE
    //-----------------------------------------------------------------------

    let mut k_coeff = Complex::<T>::cis(order.fract() * rotation_angle);
    if !order.to_usize().unwrap().is_even() {
        k_coeff = -k_coeff;
    }
    let mut k_value = k_value[0];
    let mut i_value = i_value[0];
    if scaling == Scaling::Scaled {
        let mut n_good_dummy = 0;
        if underflow_add_i_k(negative_z, &mut k_value, &mut i_value, &mut n_good_dummy) {
            n_zeros += 1;
        }
    }
    let y = vec![k_coeff * k_value + i_coeff * i_value];
    Ok((y, n_zeros))
}
