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
        i_pow_n,
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
/// y[i] = T::C_ZERO for i in nlast+1..n
///
/// The logic is very similar to i_uniform_asymp1 and the flow control comments from that
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
            * i_pow_n(integer_order + effective_n - 1);
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
    match OverflowState::check(exponent.re, T::C_ONE, T::ZERO, mc) {
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
    let mut param_seeds: [Option<DebyeParams<T>>; 2] = [None, None];
    let mut i_recurrence_seeds = [T::C_ZERO; 2];

    let mut n_elements_set = 0;
    let mut y = T::c_zeros(n);
    let mut k_overflow_state = OverflowState::NearUnder;

    for i in 0..n {
        n_elements_set = i + 1;
        let modified_order = order + T::from_usize(i);

        // Note: use z_right_half so Re(z) >= 0
        let params = DebyeParams::compute(z_right_half, modified_order);
        param_seeds[1] = param_seeds[0];
        param_seeds[0] = Some(params);

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
                    i_recurrence_seeds[found_one_good_seed as usize] = bessel_value;
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
        let [rec_prev, rec_curr] = i_recurrence_seeds;
        scale_controlled_recurrence(
            true,
            order,
            z_right_half,
            Some(&mut y),
            n_elements_set,
            n,
            rec_prev,
            rec_curr,
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

        // Reuse from stack if i is the last or second-to-last item from the first loop:
        let params = if n_elements_set > 0 && i == n_elements_set - 1 {
            param_seeds[0].unwrap()
        } else if n_elements_set > 1 && i == n_elements_set - 2 {
            param_seeds[1].unwrap()
        } else {
            DebyeParams::compute(z_right_half, current_order)
        };

        // Use the I fields:
        let phi = params.phi_i;
        let sum = params.sum_i;

        let exponent = scaling.scale_zetas(z_right_half, current_order, params.zeta1, params.zeta2);
        let overflow = OverflowState::check(exponent.re, phi, T::ZERO, mc);
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
                let amplitude = phi * sum;

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
        i_recurrence_seeds[found_one_good_entry as usize] = i_bessel_value;
        let underflowed = i_bessel_value == T::C_ZERO;
        i_bessel_value *= i_overflow_state.reciprocal_scaling_factor::<T>(mc);
        // Handle the underflow of I + K and combine them for the continuation
        let mut k_bessel_value = *yi;
        if scaling == Scaling::Scaled
            && underflow_add_i_k(z_right_half, &mut k_bessel_value, &mut i_bessel_value, mc)
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

    Ok(i_k_mixing_recurrence(
        order,
        z_right_half,
        scaling,
        two_over_z,
        i_recurrence_seeds,
        y,
        continuation_phase,
        i_overflow_state,
        n_zeros,
        remaining_n,
        mc,
    ))
}

/// k_uniform_asymp2 computes K(order, z) and its analytic continuation from the
/// right half plane to the left half plane by means of the
/// uniform asymptotic expansions for H(kind, order, zn) and J(order, zn)
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
    let h1_airy_rotation: Complex<T> =
        Complex::<T>::new(T::one(), T::from_f64(1.732_050_807_568_877_2));
    let h2_airy_rotation: Complex<T> =
        Complex::<T>::new(-T::half(), -T::from_f64(8.660_254_037_844_386e-1));

    let mut n_zeros = 0;
    let mut y = T::c_zeros(n);

    let integer_order = order.to_usize().unwrap();
    let order_fract = order.fract();
    let angle = -T::FRAC_PI_2() * order_fract;
    let conversion_phase_fract = -T::I * Complex::<T>::from_polar(T::FRAC_PI_2(), angle);
    let mut h2_to_k_factor =
        h1_airy_rotation * conversion_phase_fract * i_pow_n(integer_order).conj();

    let (z_right_half, z_was_flipped_into_right_half) = if z.re < T::ZERO {
        (-z, true)
    } else {
        (z, false)
    };
    let (z_first_quadrant, z_was_flipped_up) = if z_right_half.im <= T::ZERO {
        (z_right_half.conj(), true)
    } else {
        (z_right_half, false)
    };
    let z_rotated = -T::I * z_first_quadrant;

    // To calculate K_v(z), this function uses the analytic continuation from the Hankel function:
    //     K_v(z) = -i * (pi/2) * e^(-i * v * pi/2) * H^{(2)}_v(-i * z)
    //
    // It evaluates the uniform asymptotic expansion for H^{(2)}_v(-i * z) in the first
    // quadrant, and then applies the conversion factors to map it back to K_v(z).
    //
    // - z_first_quadrant: The input z forced into the first quadrant where the expansion is valid.
    // - z_rotated: The value (-i * z_first_quadrant) which is the rotated argument passed into H^{(2)}_v.
    // - h1_airy_rotation: (2 * e^(i * pi/3)), the overall amplitude factor from the H^{(2)}_v expansion.
    // - h2_airy_rotation: (e^(-i * 2*pi/3)), the rotation applied to the inner Airy function arguments.
    // - conversion_phase_fract: The fractional order part of the (-i * (pi/2) * e^(-i * v * pi/2)) term.
    // - h2_to_k_factor: The baked-in proportionality constant that maps the raw Airy sum directly into K_v.
    let mut i_recurrence_seeds = [T::C_ZERO; 2];
    let mut param_seeds: [Option<AiryParams<T>>; 2] = [None, None];
    let mut k_overflow_state = OverflowState::None;
    let mut n_elements_set = 0;

    // Loop 1: Evaluate K_v(z) moving upwards in order (from n=0 to N).
    // The sequence evaluates from the lowest order to quickly identify the boundary
    // of underflow. The loop breaks early as soon as it finds two consecutive valid
    // non-underflowing values, which will serve as seeds for recurrence. It keeps a
    // rolling window of the last two computed AiryParams (`param_seeds`), which
    // caches the exact parameters at the underflow boundary for the second loop.
    let mut found_one_good_entry = false;
    for i in 0..n {
        n_elements_set = i + 1;

        let current_order = order + T::from_usize(i);
        let params = AiryParams::compute(z_rotated, current_order);
        param_seeds[1] = param_seeds[0];
        param_seeds[0] = Some(params);
        let AiryParams {
            phi,
            arg,
            zeta1,
            zeta2,
            asum,
            bsum,
            ..
        } = params;

        let exponent = -scaling.scale_zetas(z_first_quadrant, current_order, zeta1, zeta2);
        let overflow = OverflowState::check(
            exponent.re,
            phi,
            T::from_f64(-0.25) * arg.abs().ln() - T::from_f64(AIC),
            mc,
        );

        let mut handle_underflow = |of_already: &mut bool, h2_to_k_factor_: &mut Complex<T>| {
            // if z.re < 0.0, then the I function will overflow, so return an error
            if z_was_flipped_into_right_half {
                return Err(BesselError::Overflow);
            }
            *of_already = false;
            y[i] = T::C_ZERO;
            n_zeros += 1;
            *h2_to_k_factor_ *= -T::I;
            if i != 0 && y[i - 1] != T::C_ZERO {
                y[i - 1] = T::C_ZERO;
                n_zeros += 1;
            }
            Ok(())
        };

        if !found_one_good_entry {
            k_overflow_state = overflow;
        }

        match overflow {
            OverflowState::Over { .. } => return Err(BesselError::Overflow),
            OverflowState::Under { .. } => {
                handle_underflow(&mut found_one_good_entry, &mut h2_to_k_factor)?
            }
            OverflowState::NearOver | OverflowState::NearUnder | OverflowState::None => {
                let rotated_airy_arg = h2_airy_rotation * arg;

                let (airy, d_airy) = airy_pair(rotated_airy_arg);
                let inner_sum = ((d_airy * bsum) * h2_airy_rotation + (airy * asum)) * phi;
                let mut k_bessel_value = inner_sum * h2_to_k_factor;
                let exp_factor = exponent.exp() * k_overflow_state.scaling_factor::<T>(mc);
                k_bessel_value *= exp_factor;
                if k_overflow_state == OverflowState::NearUnder
                    && will_underflow(k_bessel_value, mc)
                {
                    handle_underflow(&mut found_one_good_entry, &mut h2_to_k_factor)?
                }
                if z_was_flipped_up {
                    k_bessel_value = k_bessel_value.conj();
                }
                i_recurrence_seeds[found_one_good_entry as usize] = k_bessel_value;
                y[i] = k_bessel_value * k_overflow_state.reciprocal_scaling_factor::<T>(mc);
                h2_to_k_factor = -T::I * h2_to_k_factor;
                if found_one_good_entry {
                    break;
                }
                found_one_good_entry = true;
            }
        };
    }

    let two_over_z = two_over_z_safe(z_right_half);
    let do_overflow_check = n_elements_set < n;
    let mut max_order_params: Option<AiryParams<T>> = None;
    if do_overflow_check {
        // Test the last member for overflow/underflow to shortcut the rest.
        let max_order = order + T::from_usize(n - 1);

        let params = AiryParams::compute(z_rotated, max_order);
        max_order_params = Some(params);
        let AiryParams {
            phi, zeta1, zeta2, ..
        } = params;
        let exponent = -scaling.scale_zetas(z_first_quadrant, max_order, zeta1, zeta2);
        match OverflowState::check(exponent.re, phi, T::ZERO, mc) {
            OverflowState::Over { .. } => return Err(BesselError::Overflow),

            OverflowState::Under { .. } => {
                if z_was_flipped_into_right_half {
                    return Err(BesselError::Overflow);
                }
                return Ok((T::c_zeros(n), n_zeros));
            }
            OverflowState::NearOver | OverflowState::None | OverflowState::NearUnder => (),
        }
        let [rec_prev, rec_curr] = i_recurrence_seeds;
        scale_controlled_recurrence(
            true,
            order,
            z_right_half,
            Some(&mut y),
            n_elements_set,
            n,
            rec_prev,
            rec_curr,
            k_overflow_state,
            mc,
        );
    }
    if rotation == RotationDirection::None {
        return Ok((y, n_zeros));
    }

    // When Re(z) < 0.0, the K function is not single-valued and requires an analytic
    // continuation from the right half-plane using the relation:
    //     K_v(z) = e^{\mp i v \pi} K_v(-z) \mp i \pi I_v(-z)
    // The K_v(-z) term is already stored in `y`. The following block computes the
    // second term, ±iπ I_v(-z), and adds it to complete the continuation.
    n_zeros = 0;
    let sgn = -T::PI() * T::from_f64(rotation.signum());
    let signed_pi = if z_was_flipped_up { -sgn } else { sgn };
    let modified_integer_order = integer_order + n - 1;
    let mut continuation_phase = Complex::<T>::cis(order_fract * sgn);
    if modified_integer_order.is_odd() {
        continuation_phase = -continuation_phase;
    }
    // To compute the ±iπ I_v(z) term for the analytic continuation, this block evaluates
    // the asymptotic expansion for J_v(-iz). The variable `j_to_continuation_i_factor`
    // computes the combined transformation coefficient (±iπ * e^{i v π/2}) that maps
    // the un-rotated J_v(-iz) Airy sum directly into the final ±iπ I_v(z) term.
    let cos_sin = Complex::<T>::cis(angle);
    let mut j_to_continuation_i_factor = signed_pi * Complex::<T>::new(cos_sin.im, cos_sin.re);
    j_to_continuation_i_factor *= i_pow_n(modified_integer_order);

    found_one_good_entry = false;
    let mut i_overflow_state = OverflowState::None;
    let mut remaining_n = n;

    // Loop 2: Evaluate the analytic continuation term I_v(z) for Re(z) < 0.
    // Here, we iterate downwards from the maximum order (N). Like Loop 1, this loop
    // breaks early if it finds two consecutive valid non-underflowing I_v(z) values.
    // If it reaches the underflow boundary found by Loop 1, it efficiently reuses the
    // cached `param_seeds` to avoid recomputing Airy parameters. As it evaluates
    // I_v(z), it merges it with the corresponding K_v(z) value.
    for (i, yi) in y.iter_mut().enumerate().rev() {
        remaining_n = i;
        let current_order = order + T::from_usize(i);
        let params = if n_elements_set > 0 && i == n_elements_set - 1 {
            param_seeds[0].unwrap()
        } else if n_elements_set > 1 && i == n_elements_set - 2 {
            param_seeds[1].unwrap()
        } else if i == n - 1 && do_overflow_check {
            max_order_params.unwrap()
        } else {
            AiryParams::compute(z_rotated, current_order)
        };
        let AiryParams {
            phi,
            arg,
            zeta1,
            zeta2,
            asum,
            bsum,
            ..
        } = params;
        let exponent = scaling.scale_zetas(z_first_quadrant, current_order, zeta1, zeta2);

        let overflow = OverflowState::check(
            exponent.re,
            phi,
            T::from_f64(-0.25) * arg.abs().ln() - T::from_f64(AIC),
            mc,
        );
        if !found_one_good_entry {
            i_overflow_state = if matches!(overflow, OverflowState::Under { .. }) {
                OverflowState::None
            } else {
                overflow
            };
        }
        let mut i_bessel_value = match overflow {
            OverflowState::Over { .. } => return Err(BesselError::Overflow),
            OverflowState::Under { .. } => T::C_ZERO,
            OverflowState::NearOver | OverflowState::None | OverflowState::NearUnder => {
                let (airy, d_airy) = airy_pair(arg);
                let inner_sum = ((d_airy * bsum) + (airy * asum)) * phi;
                let mut continuation_i_term = inner_sum * j_to_continuation_i_factor;
                let exp_factor = exponent.exp() * i_overflow_state.scaling_factor::<T>(mc);
                continuation_i_term *= exp_factor;
                if i_overflow_state == OverflowState::NearUnder
                    && will_underflow(continuation_i_term, mc)
                {
                    continuation_i_term = T::C_ZERO;
                }
                continuation_i_term
            }
        };
        if z_was_flipped_up {
            i_bessel_value = i_bessel_value.conj();
        }
        i_recurrence_seeds[found_one_good_entry as usize] = i_bessel_value;

        let underflowed = i_bessel_value == T::C_ZERO;
        i_bessel_value *= i_overflow_state.reciprocal_scaling_factor::<T>(mc);

        let mut k_bessel_value = *yi;
        if scaling == Scaling::Scaled
            && underflow_add_i_k(z_right_half, &mut k_bessel_value, &mut i_bessel_value, mc)
        {
            n_zeros += 1;
        }
        *yi = k_bessel_value * continuation_phase + i_bessel_value;
        continuation_phase = -continuation_phase;
        j_to_continuation_i_factor *= -T::I;
        if underflowed {
            found_one_good_entry = false;
        } else {
            if found_one_good_entry {
                break;
            }
            found_one_good_entry = true;
        }
    }

    Ok(i_k_mixing_recurrence(
        order,
        z_right_half,
        scaling,
        two_over_z,
        i_recurrence_seeds,
        y,
        continuation_phase,
        i_overflow_state,
        n_zeros,
        remaining_n,
        mc,
    ))
}

#[allow(clippy::too_many_arguments)]
fn i_k_mixing_recurrence<T: BesselFloat>(
    order: T,
    z_right_half: Complex<T>,
    scaling: Scaling,
    two_over_z: Complex<T>,
    i_recurrence_seeds: [Complex<T>; 2],
    mut y: Vec<Complex<T>>,
    mut continuation_phase: Complex<T>,
    mut i_overflow_state: OverflowState,
    mut n_zeros: usize,
    remaining_n: usize,
    mc: &MachineConsts<T>,
) -> (Vec<Complex<T>>, usize) {
    if remaining_n == 0 {
        return (y, n_zeros);
    }

    // Loop 3: Fill the rest of the array via backward recurrence.
    // If Loop 2 stopped early because it found two valid non-underflowing seeds,
    // we can compute the remainder of the I_v sequence much faster using the
    // backward recurrence relation: I_{v-1}(z) = I_{v+1}(z) + (2v/z) * I_v(z).
    // The I_v terms are dynamically scaled during recurrence to avoid overflow.
    // As they are computed, they are combined with the K_v(z) values already in `y`.

    // 1. Calculate the next terms in the sequence
    let [mut i_rec_prev, mut i_rec_curr] = i_recurrence_seeds;

    let mut recip_scale_factor = i_overflow_state.reciprocal_scaling_factor::<T>(mc);
    let mut ascle = i_overflow_state.boundary::<T>(mc);
    for (i, yi) in y.iter_mut().enumerate().take(remaining_n).rev() {
        let current_order = order + T::from_usize(i + 1);
        (i_rec_prev, i_rec_curr) = (
            i_rec_curr,
            i_rec_prev + (current_order * two_over_z) * i_rec_curr,
        );

        // 2. Prepare the I and K values
        let mut i_bessel_value = i_rec_curr * recip_scale_factor;
        let i_bessel_value_unscaled = i_bessel_value;
        let mut k_bessel_value = *yi;
        if scaling == Scaling::Scaled
            && underflow_add_i_k(z_right_half, &mut k_bessel_value, &mut i_bessel_value, mc)
        {
            n_zeros += 1;
        }
        // 3. Combine them using the analytic continuation formula!
        *yi = k_bessel_value * continuation_phase + i_bessel_value;
        continuation_phase = -continuation_phase;
        i_overflow_state.scale_recurrence(
            &mut i_rec_prev,
            &mut i_rec_curr,
            i_bessel_value_unscaled,
            &mut ascle,
            &mut recip_scale_factor,
            mc,
        );
    }
    (y, n_zeros)
}
