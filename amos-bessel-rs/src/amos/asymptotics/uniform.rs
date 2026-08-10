use num::{Complex, Integer, complex::ComplexFloat};

use crate::{
    BesselError, Scaling,
    amos::{
        IKType, MachineConsts, RotationDirection,
        airy::airy_pair,
        asymptotics::{
            AiryGeometry, DebyeGeometry,
            uniform_params::{AiryParams, DebyeParams},
        },
        i_pow,
        limits::{OverflowState, check_underflow_uniform_asymp_params, underflow_add_i_k},
        recurrence::scale_controlled_recurrence,
        utils::{AIC, two_over_z_safe, will_underflow},
    },
    types::{BesselFloat, BesselResult},
};

/// i_uniform_asymp1 computes I(fnu,z)  by means of the uniform asymptotic
/// expansion for I(fnu,z) in -pi/3 <= arg(z) <= pi/3.
///
/// asymptotic_order_limit is the smallest order permitted for the asymptotic
/// expansion.
///
/// nlast=0 means all of the y values were set.
/// nlast != 0 is the number left to be computed by another formula for orders `order..order+nlast-1`
/// because `fnu+nlast-1 < asymptotic_order_limit`.
///
/// y[i] = zero for i = nlast+1,n
///
/// Originally ZUNI1
pub(crate) fn i_uniform_asymp1<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    n: usize,
    y: &mut [Complex<T>],
) -> Result<(usize, usize), BesselError<T>> {
    let mc: &MachineConsts<T> = T::MACHINE_CONSTANTS;
    let mut n_zeros = 0;
    let mut n_remaining = n;

    // First check for complete underflow and overflow on the first member (n=1)
    let limited_order = order.max(T::one());
    let DebyeGeometry { zeta1, zeta2, .. } = DebyeGeometry::compute(z, limited_order);
    let exponent = scaling.scale_zetas(z, limited_order, zeta1, zeta2);

    // phi = 1 is chosen here for refined tests to equal the original tests.
    // However, the was_refined flag is never checked, so the value used for
    // refinement has no effect anyway
    match OverflowState::check(exponent.re, T::C_ONE, T::ZERO, mc) {
        OverflowState::Over { was_refined: _ } => return Err(BesselError::Overflow),
        OverflowState::Under { was_refined: _ } => return Ok((n, 0)),
        _ => (),
    }
    let mut overflow_state = OverflowState::None; // this value should never be used
    let mut recurrence_seeds = [T::C_ZERO; 2];

    // If handle_underflow returns true, either the search has failed (all elements zero)
    // or the order has dropped below the assymptotic limit.
    // Either way the outer function should return.
    // If the closure returns false, the search can continue.
    let mut handle_underflow =
        |n_remaining: &mut usize, y: &mut [Complex<T>]| -> Result<bool, BesselError<T>> {
            // Set the y value to zero, and inc/decrement counters.
            y[*n_remaining - 1] = T::C_ZERO;
            n_zeros += 1;
            *n_remaining -= 1;

            // If we have no more values to test, tell the outer function to return
            if *n_remaining == 0 {
                return Ok(true);
            }

            // The line below lets us test for underflow on many elements in one go.
            let n_underflow = check_underflow_uniform_asymp_params(
                z,
                order,
                scaling,
                IKType::I,
                *n_remaining,
                y,
                mc,
            )?;

            // Again, inc/decrement and return if we ran out of values to test
            *n_remaining -= n_underflow;
            n_zeros += n_underflow;
            if *n_remaining == 0 {
                return Ok(true);
            }

            // Now check whether the decremented effective order has dropped
            // below the assymptotic limit
            let effective_order = order + T::from_usize(*n_remaining - 1);
            if effective_order < mc.asymptotic_order_limit {
                return Ok(true);
            }
            Ok(false)
        };

    // As the high-order I values can underflow to zero easily (and break the backward recurrence)
    // This outer loop is retrying the inner loop until it succeds in finding two non-zero
    // seeds for recurrence
    'retry_find_seeds: loop {
        for i in 0..2.min(n_remaining) {
            let current_order = order + T::from_usize(n_remaining - (i + 1));
            let DebyeParams {
                phi_i: phi,
                zeta1,
                zeta2,
                sum_i: sum,
                ..
            } = DebyeParams::compute(z, current_order);

            let mut exponent = scaling.scale_zetas(z, current_order, zeta1, zeta2);
            if scaling == Scaling::Scaled {
                exponent += Complex::<T>::new(T::ZERO, z.im);
            }

            let overflow = OverflowState::check(exponent.re, phi, T::ZERO, mc);
            if i == 0 {
                overflow_state = overflow;
            }
            match overflow {
                OverflowState::Over { .. } => return Err(BesselError::Overflow),
                OverflowState::Under { .. } => {
                    if handle_underflow(&mut n_remaining, y)? {
                        return Ok((n_zeros, n_remaining));
                    }
                    continue 'retry_find_seeds;
                }
                _ => (),
            }

            let amplitude = phi * sum;
            let exp_factor = overflow_state.scaling_factor::<T>(mc) * exponent.exp();
            let bessel_value = amplitude * exp_factor;
            if overflow_state == OverflowState::NearUnder && will_underflow(bessel_value, mc) {
                if handle_underflow(&mut n_remaining, y)? {
                    return Ok((n_zeros, n_remaining));
                }
                continue 'retry_find_seeds;
            }
            recurrence_seeds[i] = bessel_value;
            y[n_remaining - i - 1] =
                bessel_value * overflow_state.reciprocal_scaling_factor::<T>(mc);
        }
        // if the loop above completed sucessfully, we found two non-zero seeds
        break 'retry_find_seeds;
    }

    if n_remaining > 2 {
        let [rec_prev, rec_curr] = recurrence_seeds;
        scale_controlled_recurrence(
            false,
            order,
            z,
            Some(y),
            n_remaining - 2,
            n,
            rec_prev,
            rec_curr,
            overflow_state,
            mc,
        );
    }

    Ok((n_zeros, 0))
}

/// i_uniform_asymp2 computes I(fnu, z) in the right half plane by means of
/// uniform asymptotic expansion for J(fnu, zn) where zn is z*i
/// or -z*i and zn is in the right half plane also.
///
/// asymptotic_order_limit is the smallest order permitted for the asymptotic
/// expansion.
///
/// nlast=0 means all of the y values were set.
/// nlast != 0 is the number left to be computed by another
/// formula for orders `order..order+nlast-1` because
/// `order+nlast-1 < asymptotic_order_limit`.
///
/// y[i] = czero for i in nlast+1..n
///
/// The logic is very similar to i_uniform_asymp2 and the flow control comments from that
/// function apply here too.
///
/// Originally ZUNI2
pub(crate) fn i_uniform_asymp2<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    n: usize,
    y: &mut [Complex<T>],
) -> Result<(usize, usize), BesselError<T>> {
    let mc: &MachineConsts<T> = T::MACHINE_CONSTANTS;
    let mut n_zeros = 0;
    let mut n_remaining = n;

    // We force z into the upper half plane to calculate, and then conjugate the answer
    // back down at the very end if necessary.
    // z_uppper is always in the upper half plane (z_upper.im > 0)
    // z_rotated is in the right half plane after rotation by I or -I
    let (sign_of_i, z_was_flipped, z_upper) = if z.im <= T::ZERO {
        (T::ONE, true, z.conj())
    } else {
        (-T::ONE, false, z)
    };
    let z_rotated = -T::I * z_upper;
    let integer_order = order.to_usize().unwrap();

    // Computes i^order by separating integer/fractional parts to avoid precision loss in trig functions
    let calculate_rotation_factor = |effective_n: usize| {
        let r = Complex::<T>::cis(T::FRAC_PI_2() * order.fract())
            * i_pow(integer_order + effective_n - 1);
        if z_was_flipped { r.conj() } else { r }
    };

    let mut rotation_factor = calculate_rotation_factor(n);

    // First check for complete underflow and overflow on the first member (n=1)
    let limited_order = order.max(T::one());
    let AiryGeometry { zeta1, zeta2, .. } = AiryGeometry::compute(z_rotated, limited_order);

    // 1. Calculate the exponent
    let exponent = scaling.scale_zetas(z_upper, limited_order, zeta1, zeta2);

    // phi = 1 is chosen here for refined tests to equal the original tests.
    // However, the was_refined flag is never checked, so the value used for
    // refinement has no effect anyway
    match OverflowState::check(exponent.re, T::C_ONE, T::zero(), mc) {
        OverflowState::Over { .. } => return Err(BesselError::Overflow),
        OverflowState::Under { .. } => return Ok((n, 0)),
        _ => (),
    }

    debug_assert!(limited_order + T::from_usize(n - 1) > mc.asymptotic_order_limit);

    // If handle_underflow returns true, either the search has failed (all elements zero)
    // or the order has dropped below the assymptotic limit.
    // Either way the outer function should return.
    // If the closure returns false, the search can continue.
    //
    // The underflow handling and retry loop below is identical in structure to i_uniform_asymp1.
    // See the comments in that function for a detailed explanation of the recurrence seeding logic.
    let mut handle_underflow = |n_remaining: &mut usize,
                                rotation_factor: &mut Complex<T>,
                                y: &mut [Complex<T>]|
     -> Result<bool, BesselError<T>> {
        y[*n_remaining - 1] = T::C_ZERO;
        n_zeros += 1;
        *n_remaining -= 1;

        if *n_remaining == 0 {
            return Ok(true);
        }
        let n_underflow = check_underflow_uniform_asymp_params(
            z,
            order,
            scaling,
            IKType::I,
            *n_remaining,
            y,
            mc,
        )?;
        *n_remaining -= n_underflow;
        n_zeros += n_underflow;
        if *n_remaining == 0 {
            return Ok(true);
        }
        let modified_order = order + T::from_usize(*n_remaining - 1);
        if modified_order < mc.asymptotic_order_limit {
            return Ok(true);
        }
        *rotation_factor = calculate_rotation_factor(*n_remaining);
        Ok(false)
    };

    let mut overflow_state = OverflowState::NearUnder;
    let mut recurrence_seeds = [T::C_ZERO; 2];
    // This outer loop says: keep looping until you found two non-zero seeds for the recurrence
    // formula
    'retry_find_seeds: loop {
        for i in 0..2.min(n_remaining) {
            let current_order = order + T::from_usize(n_remaining - (i + 1));
            let AiryParams {
                phi,
                arg,
                zeta1,
                zeta2,
                asum,
                bsum,
                ..
            } = AiryParams::compute(z_rotated, current_order);

            let mut exponent = scaling.scale_zetas(z_upper, current_order, zeta1, zeta2);
            if scaling == Scaling::Scaled {
                exponent += T::I * z.im.abs();
            }

            let of = OverflowState::check(
                exponent.re,
                phi,
                T::from_f64(-0.25) * arg.abs().ln() - T::from_f64(AIC),
                mc,
            );
            if i == 0 {
                overflow_state = of;
            }
            match of {
                OverflowState::Over { .. } => return Err(BesselError::Overflow),
                OverflowState::Under { .. } => {
                    if handle_underflow(&mut n_remaining, &mut rotation_factor, y)? {
                        return Ok((n_zeros, n_remaining));
                    }
                    continue 'retry_find_seeds;
                }
                _ => (),
            }
            let (a_airy, d_airy) = airy_pair(arg);

            // 2. Combine the Airy functions, Debye sums, and phi into the amplitude
            let airy_sum = a_airy * asum + d_airy * bsum;
            let amplitude = phi * airy_sum;

            // 3. Exponentiate the exponent
            let exp_factor = overflow_state.scaling_factor::<T>(mc) * exponent.exp();

            // 4. Combine them to get the un-rotated Bessel evaluation
            let mut bessel_value = amplitude * exp_factor;

            if overflow_state == OverflowState::NearUnder && will_underflow(bessel_value, mc) {
                if handle_underflow(&mut n_remaining, &mut rotation_factor, y)? {
                    return Ok((n_zeros, n_remaining));
                }
                continue 'retry_find_seeds;
            }

            // 5. Finally, rotate the value from the J domain back into the I domain
            if z_was_flipped {
                bessel_value = bessel_value.conj();
            }
            bessel_value *= rotation_factor;

            recurrence_seeds[i] = bessel_value;
            y[n_remaining - i - 1] =
                bessel_value * overflow_state.reciprocal_scaling_factor::<T>(mc);

            // Step the rotation factor backwards for the next order (n - 1)
            rotation_factor *= sign_of_i * T::I;
        }
        break 'retry_find_seeds;
    }

    // If we found two seeds, and still have come n to calculate, then do it
    // with backward recurrence
    if n_remaining > 2 {
        let [s1, s2] = recurrence_seeds;
        scale_controlled_recurrence(
            false,
            order,
            z,
            Some(y),
            n_remaining - 2,
            n,
            s1,
            s2,
            overflow_state,
            mc,
        );
    }
    Ok((n_zeros, 0))
}

/// k_uniform_asymp1 computes K(order, z) and its analytic continuation from the
/// right half plane to the left half plane by means of the  uniform asymptotic expansion.
///
/// `rotation` indicates the direction of rotation for analytic continuation.
///
/// Originally ZUNK1
pub(crate) fn k_uniform_asymp1<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    rotation: RotationDirection,
    n: usize,
) -> BesselResult<T> {
    let mc: &MachineConsts<T> = T::MACHINE_CONSTANTS;

    let mut found_one_good_seed = false;
    let mut n_zeros = 0;
    let (z_right_half, z_was_flipped) = if z.re < T::ZERO {
        (-z, true)
    } else {
        (z, false)
    };
    let mut debye_seeds: [Option<DebyeParams<T>>; 2] = [None, None];
    let mut recurrence_seeds = [T::C_ZERO; 2];

    let mut n_elements_set = 0;
    let mut y = T::c_zeros(n);
    let mut k_overflow_state = OverflowState::NearUnder;

    for i in 0..n {
        n_elements_set = i + 1;
        let modified_order = order + T::from_usize(i);

        // Note: use z_right_half so Re(z) >= 0
        let params = DebyeParams::compute(z_right_half, modified_order);
        if i < 2 {
            debye_seeds[i] = Some(params);
        }

        // Use the K fields:
        let phi = params.phi_k;
        let sum = params.sum_k;
        let exponent =
            -scaling.scale_zetas(z_right_half, modified_order, params.zeta1, params.zeta2);
        let overflow = OverflowState::check(exponent.re, phi, T::ZERO, mc);
        if !found_one_good_seed {
            k_overflow_state = overflow;
        }
        match overflow {
            OverflowState::Over { .. } => return Err(BesselError::Overflow),
            OverflowState::Under { .. } => {
                if z_was_flipped {
                    return Err(BesselError::Overflow);
                }
                found_one_good_seed = false;
                y[i] = T::C_ZERO;
                n_zeros += 1;
            }
            OverflowState::None | OverflowState::NearOver | OverflowState::NearUnder => {
                let amplitude = phi * sum;
                let exp_factor = k_overflow_state.scaling_factor::<T>(mc) * exponent.exp();
                let bessel_value = amplitude * exp_factor;

                let will_underflow = will_underflow(bessel_value, mc);
                if k_overflow_state != OverflowState::NearUnder || !will_underflow {
                    recurrence_seeds[found_one_good_seed as usize] = bessel_value;
                    y[i] = bessel_value * k_overflow_state.reciprocal_scaling_factor::<T>(mc);
                    if found_one_good_seed {
                        // if we already found one, we've now found another so break out of the loop
                        break;
                    }
                    found_one_good_seed = true;
                } else if will_underflow {
                    if z_was_flipped {
                        return Err(BesselError::Overflow);
                    }
                    y[i] = T::C_ZERO;
                    n_zeros += 1;
                    if i > 0 && y[i - 1] != T::C_ZERO {
                        y[i - 1] = T::C_ZERO;
                        n_zeros += 1
                    }
                }
            }
        };
    }

    let two_over_z = two_over_z_safe(z_right_half);
    if n_elements_set < n {
        // Test the last member of the requested sequence for overflow/underflow.
        // If it underflows, the rest of the sequence is zeroed out since K_v(z) grows with order.
        let max_order = order + T::from_usize(n - 1);
        let DebyeGeometry {
            phi_k: phi,
            zeta1: zet1d,
            zeta2: zet2d,
            ..
        } = DebyeGeometry::compute(z_right_half, max_order);
        let overflow_test = -scaling.scale_zetas(z_right_half, max_order, zet1d, zet2d);

        match OverflowState::check(overflow_test.re.abs(), phi, T::ZERO, mc) {
            OverflowState::Over { .. } => return Err(BesselError::Overflow),
            OverflowState::Under { .. } => {
                return if z_was_flipped {
                    Err(BesselError::Overflow)
                } else {
                    Ok((vec![T::C_ZERO; n], n))
                };
            }
            _ => (),
        }
        // Forward recurse to populate the remainder of the K sequence.
        let [s1, s2] = recurrence_seeds;
        scale_controlled_recurrence(
            true,
            order,
            z_right_half,
            Some(&mut y),
            n_elements_set,
            n,
            s1,
            s2,
            k_overflow_state,
            mc,
        );
    }
    if rotation == RotationDirection::None {
        return Ok((y, n_zeros));
    }
    // Perform analytic continuation if z was originally in the left half-plane (or if rotation was explicitly requested).
    // The continuation formula is: K_v(z * e^m*pi*i) = e^-m*v*pi*i * K_v(z) - i*pi * I_v(z)
    n_zeros = 0;
    let rotation_angle = -T::PI() * T::from_f64(rotation.signum());

    let integer_order = order.to_i64().unwrap();
    let order_frac = order.fract();
    let modified_int_order = integer_order + (n as i64) - 1;
    let mut continuation_phase = Complex::<T>::cis(order_frac * rotation_angle);
    if (modified_int_order % 2) != 0 {
        continuation_phase = -continuation_phase;
    }
    let mut found_one_good_entry = false;
    let mut i_overflow_state = OverflowState::None;
    let mut remaining_n = n;
    for (i, yi) in y.iter_mut().enumerate().rev() {
        remaining_n = i;
        let current_order = order + T::from_usize(i);

        // Reuse from stack if i is 0 or 1, otherwise compute fresh:
        let params = if i < 2 {
            debye_seeds[i].unwrap_or_else(|| DebyeParams::compute(z_right_half, current_order))
        } else {
            DebyeParams::compute(z_right_half, current_order)
        };

        // Use the I fields:
        let phid = params.phi_i;
        let sumd = params.sum_i;

        let exponent = scaling.scale_zetas(z_right_half, current_order, params.zeta1, params.zeta2);
        let overflow = OverflowState::check(exponent.re, phid, T::ZERO, mc);
        if !found_one_good_entry && !matches!(overflow, OverflowState::Under { .. }) {
            i_overflow_state = overflow;
        }
        let mut i_bessel_value = match overflow {
            OverflowState::Over { .. } => {
                return Err(BesselError::Overflow);
            }
            OverflowState::Under { .. } => T::C_ZERO,
            OverflowState::NearOver | OverflowState::NearUnder | OverflowState::None => {
                // 1. Calculate the I amplitude
                let amplitude = phid * sumd;

                // 2. Multiply by ±iπ for the analytic continuation
                let continuation_amplitude = T::I * amplitude * rotation_angle;

                // 3. Exponentiate the I exponent
                let exp_factor = exponent.exp() * i_overflow_state.scaling_factor::<T>(mc);

                // 4. Combine into the final I value
                let mut i_bessel_value = continuation_amplitude * exp_factor;

                if i_overflow_state == OverflowState::NearUnder
                    && will_underflow(i_bessel_value, mc)
                {
                    i_bessel_value = T::C_ZERO;
                }
                i_bessel_value
            }
        };
        recurrence_seeds[found_one_good_entry as usize] = i_bessel_value;
        let underflowed = i_bessel_value == T::C_ZERO;
        i_bessel_value *= i_overflow_state.reciprocal_scaling_factor::<T>(mc);
        // Handle the underflow of I + K and combine them for the continuation
        let mut k_bessel_value = *yi;
        let mut dummy_n_good = 0;
        if scaling == Scaling::Scaled
            && underflow_add_i_k(
                z_right_half,
                &mut k_bessel_value,
                &mut i_bessel_value,
                &mut dummy_n_good,
                mc,
            )
        {
            n_zeros += 1;
        }
        *yi = k_bessel_value * continuation_phase + i_bessel_value;
        continuation_phase = -continuation_phase;
        if underflowed {
            found_one_good_entry = false;
        } else {
            if found_one_good_entry {
                break;
            }
            found_one_good_entry = true;
        }
    }
    if remaining_n > 0 {
        // Recurse backward to populate the remainder of the I sequence, adding it to the
        // existing K sequence. The recurrence is dynamically scaled to prevent overflow.
        let [mut s1, mut s2] = recurrence_seeds;
        let mut reciprocal_scale_factor = i_overflow_state.reciprocal_scaling_factor::<T>(mc);
        let mut boundary = i_overflow_state.boundary::<T>(mc);
        for (i, yi) in y.iter_mut().enumerate().take(remaining_n).rev() {
            let current_order = order + T::from_usize(i + 1);
            // 1. Calculate the next terms in the sequence
            (s1, s2) = (s2, s1 + current_order * (two_over_z * s2));

            // 2. Prepare the I and K values
            let mut i_bessel_value = s2 * reciprocal_scale_factor;
            let mut k_bessel_value = *yi;

            let mut dummy_n_good = 0;
            if scaling == Scaling::Scaled
                && underflow_add_i_k(
                    z_right_half,
                    &mut k_bessel_value,
                    &mut i_bessel_value,
                    &mut dummy_n_good,
                    mc,
                )
            {
                n_zeros += 1;
            }
            // 3. Combine them using the analytic continuation formula!
            *yi = k_bessel_value * continuation_phase + i_bessel_value;
            continuation_phase = -continuation_phase;
            i_overflow_state.scale_recurrence(
                &mut s1,
                &mut s2,
                i_bessel_value,
                &mut boundary,
                &mut reciprocal_scale_factor,
                mc,
            );
        }
    }
    Ok((y, n_zeros))
}

/// k_uniform_asymp2 computes K(fnu,z) and its analytic continuation from the
/// right half plane to the left half plane by means of the
/// uniform asymptotic expansions for H(kind,fnu,zn) and J(fnu,zn)
/// where zn is in the right half plane, kind=(3-mr)/2, mr=+1 or
/// -1. here zn=zr*i or -zr*i where zr=z if z is in the right
/// half plane or zr=-z if z is in the left half plane.
///
/// `rotation` indicates the direction of rotation for analytic continuation.
///
/// Originally ZUNK2
pub(crate) fn k_uniform_asymp2<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    rotation: RotationDirection,
    n: usize,
) -> BesselResult<T> {
    let mc: &MachineConsts<T> = T::MACHINE_CONSTANTS;
    let cr1: Complex<T> = Complex::<T>::new(T::one(), T::from_f64(1.732_050_807_568_877_2));
    let cr2: Complex<T> = Complex::<T>::new(-T::half(), -T::from_f64(8.660_254_037_844_386e-1));

    let mut found_one_good_entry = false;
    let mut n_zeros = 0;
    let mut y = T::c_zeros(n);
    let zr = if z.re < T::ZERO { -z } else { z };
    let mut zn = -T::I * zr;
    let mut zb = zr;
    let integer_order = order.to_usize().unwrap();
    let order_fract = order.fract();
    let angle = -T::FRAC_PI_2() * order_fract;
    let c2 = -T::I * Complex::<T>::from_polar(T::FRAC_PI_2(), angle);
    let mut cs = cr1 * c2 * i_pow(integer_order).conj();
    if zr.im <= T::ZERO {
        zn.re = -zn.re;
        zb.im = -zb.im;
    }
    //-----------------------------------------------------------------------
    //     K(FNU,Z) IS COMPUTED FROM H(2,FNU,-I*Z) WHERE Z IS IN THE FIRST
    //     QUADRANT. FOURTH QUADRANT VALUES (YY <= 0.0E0) ARE COMPUTED BY
    //     CONJUGATION SINCE THE K FUNCTION IS REAL ON THE POSITIVE REAL AXIS
    //-----------------------------------------------------------------------

    let mut phi = [T::C_ZERO; 2];
    let mut arg = [T::C_ZERO; 2];
    let mut zeta1 = [T::C_ZERO; 2];
    let mut zeta2 = [T::C_ZERO; 2];
    let mut asum = [T::C_ZERO; 2];
    let mut bsum = [T::C_ZERO; 2];
    let mut cy = [T::C_ZERO; 2];
    let mut j = 1;
    let mut k_overflow_state = OverflowState::None;
    let mut n_elements_set = 0;

    for i in 0..n {
        n_elements_set = i + 1;
        // j flip-flops between 0 and 1 using  = 1-j
        j = 1 - j;
        let current_order = order + T::from_usize(i);
        let AiryParams {
            phi: phi_,
            arg: arg_,
            zeta1: zeta1_,
            zeta2: zeta2_,
            asum: asum_,
            bsum: bsum_,
            ..
        } = AiryParams::compute(zn, current_order);
        phi[j] = phi_;
        arg[j] = arg_;
        zeta1[j] = zeta1_;
        zeta2[j] = zeta2_;
        asum[j] = asum_;
        bsum[j] = bsum_;
        // (phi[j], arg[j], zeta1[j], zeta2[j], asum[j], bsum[j]) =
        //     hj_uniform_asymp_params(zn, current_order, false);
        let s1 = -scaling.scale_zetas(zb, current_order, zeta1[j], zeta2[j]);
        let of = OverflowState::check(
            s1.re,
            phi[j],
            T::from_f64(-0.25) * arg[j].abs().ln() - T::from_f64(AIC),
            mc,
        );

        let mut handle_underflow = |of_already: &mut bool, cs_: &mut Complex<T>| {
            //-----------------------------------------------------------------------
            //     FOR ZR < 0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
            //-----------------------------------------------------------------------
            if z.re < T::ZERO {
                return Err(BesselError::Overflow);
            }
            *of_already = false;
            y[i] = T::C_ZERO;
            n_zeros += 1;
            *cs_ *= -T::I;
            if i != 0 && y[i - 1] != T::C_ZERO {
                y[i - 1] = T::C_ZERO;
                n_zeros += 1;
            }
            Ok(())
        };

        if !found_one_good_entry {
            k_overflow_state = of;
        }

        match of {
            OverflowState::Over { .. } => return Err(BesselError::Overflow),

            OverflowState::Under { .. } => handle_underflow(&mut found_one_good_entry, &mut cs)?,
            OverflowState::NearOver | OverflowState::NearUnder | OverflowState::None => {
                //-----------------------------------------------------------------------;
                //     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR;
                //     EXPONENT EXTREMES;
                //-----------------------------------------------------------------------;
                let c2 = cr2 * arg[j];

                let (airy, d_airy) = airy_pair(c2);
                let pt = ((d_airy * bsum[j]) * cr2 + (airy * asum[j])) * phi[j];
                let mut s2 = pt * cs;
                let s1 = s1.exp() * k_overflow_state.scaling_factor::<T>(mc);
                s2 *= s1;
                if k_overflow_state == OverflowState::NearUnder && will_underflow(s2, mc) {
                    handle_underflow(&mut found_one_good_entry, &mut cs)?
                }
                if zr.im <= T::ZERO {
                    s2 = s2.conj();
                }
                cy[found_one_good_entry as usize] = s2;
                y[i] = s2 * k_overflow_state.reciprocal_scaling_factor::<T>(mc);
                cs = -T::I * cs;
                if found_one_good_entry {
                    break;
                }
                found_one_good_entry = true;
            }
        };
    }

    let two_over_z = two_over_z_safe(zr);
    let mut phid = T::C_ZERO;
    let mut argd = T::C_ZERO;
    let mut zeta1d = T::C_ZERO;
    let mut zeta2d = T::C_ZERO;
    let mut asumd = T::C_ZERO;
    let mut bsumd = T::C_ZERO;
    let do_overflow_check = n_elements_set < n;
    if do_overflow_check {
        //-----------------------------------------------------------------------;
        //     TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW. SET SEQUENCE TO ZERO;
        //     ON UNDERFLOW.;
        //-----------------------------------------------------------------------;
        let max_order = order + T::from_usize(n - 1);
        AiryParams {
            phi: phid,
            arg: argd,
            zeta1: zeta1d,
            zeta2: zeta2d,
            asum: asumd,
            bsum: bsumd,
            ..
        } = AiryParams::compute(zn, max_order);
        // (phid, argd, zeta1d, zeta2d, asumd, bsumd) =
        //     hj_uniform_asymp_params(zn, max_order, rotation == RotationDirection::None);
        let s1 = -scaling.scale_zetas(zb, max_order, zeta1d, zeta2d);
        match OverflowState::check(s1.re, phid, T::ZERO, mc) {
            OverflowState::Over { .. } => return Err(BesselError::Overflow),

            OverflowState::Under { .. } => {
                if z.re < T::ZERO {
                    return Err(BesselError::Overflow);
                }
                return Ok((T::c_zeros(n), n_zeros));
            }
            OverflowState::NearOver | OverflowState::None | OverflowState::NearUnder => (),
        }
        let [s1, s2] = cy;
        scale_controlled_recurrence(
            true,
            order,
            zr,
            Some(&mut y),
            n_elements_set,
            n,
            s1,
            s2,
            k_overflow_state,
            mc,
        );
    }
    if rotation == RotationDirection::None {
        return Ok((y, n_zeros));
    }
    //-----------------------------------------------------------------------
    //     ANALYTIC CONTINUATION FOR RE(Z) < 0.0
    //-----------------------------------------------------------------------
    n_zeros = 0;
    let sgn = -T::PI() * T::from_f64(rotation.signum());
    //-----------------------------------------------------------------------
    //     CSPN AND CSGN ARE COEFF OF K AND I FUNCIONS RESP.
    //-----------------------------------------------------------------------
    let csgn = if zr.im <= T::ZERO { -sgn } else { sgn };
    let modified_integer_order = integer_order + n - 1;
    let mut cspn = Complex::<T>::cis(order_fract * sgn);
    if modified_integer_order.is_odd() {
        cspn = -cspn;
    }
    //-----------------------------------------------------------------------
    //     CS=COEFF OF THE J FUNCTION TO GET THE I FUNCTION. I(FNU,Z) IS
    //     COMPUTED FROM EXP(I*FNU*FRAC_PI_2)*J(FNU,-I*Z) WHERE Z IS IN THE FIRST
    //     QUADRANT. FOURTH QUADRANT VALUES (YY <= 0.0E0) ARE COMPUTED BY
    //     CONJUGATION SINCE THE I FUNCTION IS REAL ON THE POSITIVE REAL AXIS
    //-----------------------------------------------------------------------;
    // TODO what's the actual maths below?
    let cos_sin = Complex::<T>::cis(angle);
    // let mut cs = Complex::<T>::I * Complex::<T>::from_polar(CSGNI, ANG);
    let mut cs = csgn * Complex::<T>::new(cos_sin.im, cos_sin.re);
    cs *= i_pow(modified_integer_order);
    let mut iuf = 0;

    found_one_good_entry = false;
    let mut i_overflow_state = OverflowState::None;
    let mut remaining_n = n;
    for (i, yi) in y.iter_mut().enumerate().rev() {
        remaining_n = i;
        let current_order = order + T::from_usize(i);
        //-----------------------------------------------------------------------
        //     LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K
        //     FUNCTION ABOVE
        //-----------------------------------------------------------------------
        // Note that, is the overflow check was done, the ___d are already set, and
        // valid for kk == n-1. Also that kk == n-1 on the first pas through this loop.
        let use_preset_overflow = (i == n - 1) && do_overflow_check;
        // these where the last two kk values where phi etc where recorded in the previous run.
        // Would it be better to store all of them?!
        let in_last_two_set = (i == n_elements_set - 1) || (i == n_elements_set - 2);
        if n <= 2 || (!use_preset_overflow) && in_last_two_set {
            phid = phi[j];
            argd = arg[j];
            zeta1d = zeta1[j];
            zeta2d = zeta2[j];
            asumd = asum[j];
            bsumd = bsum[j];
            j = 1 - j;
        } else if !(use_preset_overflow || in_last_two_set) {
            AiryParams {
                phi: phid,
                arg: argd,
                zeta1: zeta1d,
                zeta2: zeta2d,
                asum: asumd,
                bsum: bsumd,
                ..
            } = AiryParams::compute(zn, current_order);
        } else {
            // Case were overflow check has already set the ___d variables ?
        }
        let mut s1 = scaling.scale_zetas(zb, current_order, zeta1d, zeta2d);

        let of = OverflowState::check(
            s1.re,
            phid,
            T::from_f64(-0.25) * argd.abs().ln() - T::from_f64(AIC),
            mc,
        );
        if !found_one_good_entry {
            i_overflow_state = if matches!(of, OverflowState::Under { .. }) {
                OverflowState::None
            } else {
                of
            };
        }
        let mut s2 = match of {
            OverflowState::Over { .. } => return Err(BesselError::Overflow),
            OverflowState::Under { .. } => T::C_ZERO,
            OverflowState::NearOver | OverflowState::None | OverflowState::NearUnder => {
                let (airy, d_airy) = airy_pair(argd);
                let pt = ((d_airy * bsumd) + (airy * asumd)) * phid;
                let mut s2 = pt * cs;
                s1 = s1.exp() * i_overflow_state.scaling_factor::<T>(mc);
                s2 *= s1;
                if i_overflow_state == OverflowState::NearUnder && will_underflow(s2, mc) {
                    s2 = T::C_ZERO;
                }
                s2
            }
        };
        if zr.im <= T::ZERO {
            s2 = s2.conj();
        }
        cy[found_one_good_entry as usize] = s2;
        let c2 = s2;
        s2 *= i_overflow_state.reciprocal_scaling_factor::<T>(mc);
        //-----------------------------------------------------------------------;
        //     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N;
        //-----------------------------------------------------------------------;
        s1 = *yi;
        if scaling == Scaling::Scaled && underflow_add_i_k(zr, &mut s1, &mut s2, &mut iuf, mc) {
            n_zeros += 1;
        }
        *yi = s1 * cspn + s2;
        cspn = -cspn;
        cs *= -T::I;
        if c2 == T::C_ZERO {
            found_one_good_entry = false;
        } else {
            if found_one_good_entry {
                break;
            }
            found_one_good_entry = true;
        }
    }

    if remaining_n > 0 {
        //-----------------------------------------------------------------------
        //     RECUR BACKWARD FOR REMAINDER OF I SEQUENCE AND ADD IN THE
        //     K FUNCTIONS, SCALING THE I SEQUENCE DURING RECURRENCE TO KEEP
        //     INTERMEDIATE ARITHMETIC ON SCALE NEAR EXPONENT EXTREMES.
        //-----------------------------------------------------------------------
        let [mut s1, mut s2] = cy;

        let mut recip_scale_factor = i_overflow_state.reciprocal_scaling_factor::<T>(mc);
        let mut ascle = i_overflow_state.boundary::<T>(mc);
        let mut ck = (order + T::from_usize(remaining_n)) * two_over_z;
        for yi in y.iter_mut().take(remaining_n).rev() {
            (s1, s2) = (s2, s1 + ck * s2);
            ck -= two_over_z;
            let mut c2 = s2 * recip_scale_factor;
            let old_c2 = c2;
            let mut c1 = *yi;
            if scaling == Scaling::Scaled && underflow_add_i_k(zr, &mut c1, &mut c2, &mut iuf, mc) {
                n_zeros += 1;
            }
            *yi = c1 * cspn + c2;
            cspn = -cspn;
            i_overflow_state.scale_recurrence(
                &mut s1,
                &mut s2,
                old_c2,
                &mut ascle,
                &mut recip_scale_factor,
                mc,
            );
        }
    }
    Ok((y, n_zeros))
}
