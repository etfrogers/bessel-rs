#![allow(non_snake_case)]
use num::{Complex, complex::ComplexFloat};

use crate::{
    BesselError, Scaling,
    amos::{
        RotationDirection,
        asymptotics::asymptotic_i,
        i_power_series,
        limits::{Overflow, underflow_add_i_k},
        max_abs_component,
        recurrence::i_miller,
        translator::{i_right_half_plane, k_right_half_plane},
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
    let mut nz = 0;
    let zn = -z;
    let (mut y, _) = i_right_half_plane(zn, order, scaling, n)?;
    //-----------------------------------------------------------------------
    //     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
    //-----------------------------------------------------------------------
    let (cy, NW) = k_right_half_plane(zn, order, scaling, 2.min(n))?;
    if NW > 0 {
        return Err(BesselError::Overflow);
        // the NW = -1 or -2 is handled by ZBNKU returning an error,
        // but the amos code defaults to an overflow, if NW != 0
    }
    let mut s1 = cy[0];
    let SGN = -T::PI() * T::from_f64(rotation.signum());
    let mut csgn = Complex::<T>::new(T::ZERO, SGN);
    if scaling == Scaling::Scaled {
        csgn *= Complex::<T>::cis(-zn.im);
    }
    //-----------------------------------------------------------------------
    //     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
    //     WHEN FNU IS LARGE
    //-----------------------------------------------------------------------
    let mut cspn = Complex::<T>::cis(order.fract() * SGN);
    if order.to_i64().unwrap() % 2 != 0 {
        cspn = -cspn;
    }
    let mut n_good = 0;
    let mut c1 = s1;
    let mut c2 = y[0];
    if scaling == Scaling::Scaled {
        nz += underflow_add_i_k(zn, &mut c1, &mut c2, &mut n_good);
    }
    y[0] = cspn * c1 + csgn * c2;
    if n == 1 {
        return Ok((y, nz));
    }

    cspn = -cspn;
    let mut s2 = cy[1];
    c1 = s2;
    c2 = y[1];
    // this value never used, as initialised and used if scaling is needed
    let mut scaled_c2 = T::C_ZERO * T::NAN;
    if scaling == Scaling::Scaled {
        nz += underflow_add_i_k(zn, &mut c1, &mut c2, &mut n_good);
        scaled_c2 = c1;
    }
    y[1] = cspn * c1 + csgn * c2;
    if n == 2 {
        return Ok((y, nz));
    }

    cspn = -cspn;
    let rz = calc_rz(zn);
    let FN = order + T::one();
    let mut ck = FN * rz;
    //-----------------------------------------------------------------------
    //     SCALE NEAR EXPONENT EXTREMES DURING RECURRENCE ON K FUNCTIONS
    //-----------------------------------------------------------------------
    let abs_s2 = s2.abs();
    let mut overflow_state = if abs_s2 <= T::MACHINE_CONSTANTS.overflow_boundary[0] {
        Overflow::NearUnder
    } else if abs_s2 > T::MACHINE_CONSTANTS.overflow_boundary[1] {
        Overflow::NearOver
    } else {
        Overflow::None
    };
    let mut boundary = T::MACHINE_CONSTANTS.overflow_boundary[overflow_state];
    s1 *= T::MACHINE_CONSTANTS.scaling_factors[overflow_state];
    s2 *= T::MACHINE_CONSTANTS.scaling_factors[overflow_state];
    let mut recip_scaling_factor = T::MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
    for yi in y.iter_mut().skip(2) {
        //TODO common pattern below
        (s1, s2) = (s2, ck * s2 + s1);
        c1 = s2 * recip_scaling_factor;
        let mut st = c1;
        c2 = *yi;
        if scaling == Scaling::Scaled && n_good >= 0 {
            nz += underflow_add_i_k(zn, &mut c1, &mut c2, &mut n_good);
            let saved_c2 = scaled_c2;
            scaled_c2 = c1;
            if n_good == 3 {
                n_good = -4;
                s1 = saved_c2 * T::MACHINE_CONSTANTS.scaling_factors[overflow_state];
                s2 = scaled_c2 * T::MACHINE_CONSTANTS.scaling_factors[overflow_state];
                st = scaled_c2;
            }
        }
        *yi = cspn * c1 + csgn * c2;
        ck += rz;
        cspn = -cspn;
        if overflow_state != Overflow::NearOver && max_abs_component(c1) < boundary {
            overflow_state.increment();
            boundary = T::MACHINE_CONSTANTS.overflow_boundary[overflow_state];
            s1 *= recip_scaling_factor;
            s2 = st;
            s1 *= T::MACHINE_CONSTANTS.scaling_factors[overflow_state];
            s2 *= T::MACHINE_CONSTANTS.scaling_factors[overflow_state];
            recip_scaling_factor = T::MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
        }
    }
    Ok((y, nz))
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
    KODE: Scaling,
    rotation: RotationDirection,
    N: usize,
) -> BesselResult<T> {
    let mut NZ = 0;
    let zn = -z;
    let AZ = z.abs();
    let NN = N;
    let DFNU = order + T::from_usize(N - 1);
    let (mut y, _) = if (AZ * AZ * T::from_f64(0.25) <= DFNU + T::one()) || (AZ <= T::two()) {
        //-----------------------------------------------------------------------
        //     POWER SERIES FOR THE I FUNCTION
        //-----------------------------------------------------------------------
        let (y, NW_signed) = i_power_series(zn, order, KODE, NN)?;
        debug_assert!(NW_signed >= 0);
        (y, NW_signed.unsigned_abs())
    } else if AZ >= T::MACHINE_CONSTANTS.asymptotic_z_limit {
        //-----------------------------------------------------------------------
        //     ASYMPTOTIC EXPANSION FOR LARGE Z FOR THE I FUNCTION
        //-----------------------------------------------------------------------
        asymptotic_i(zn, order, KODE, NN)?
    //-----------------------------------------------------------------------
    //     MILLER ALGORITHM NORMALIZED BY THE SERIES FOR THE I FUNCTION
    //-----------------------------------------------------------------------
    } else {
        i_miller(zn, order, KODE, NN)?
    };
    //-----------------------------------------------------------------------
    //     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
    //-----------------------------------------------------------------------s
    let (cy, nz) = k_right_half_plane(zn, order, KODE, 1)?;
    if nz != 0 {
        return Err(BesselError::Overflow);
    }
    let SGN = -T::PI() * T::from_f64(rotation.signum());
    let mut csgn = Complex::<T>::new(T::ZERO, SGN);
    if KODE == Scaling::Scaled {
        csgn = T::I * csgn.im * Complex::<T>::cis(-zn.im);
    }
    //-----------------------------------------------------------------------
    //     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
    //     WHEN FNU IS LARGE
    //-----------------------------------------------------------------------
    let INU = order.to_usize().unwrap();
    let mut cspn = Complex::<T>::cis(order.fract() * SGN);
    if !INU.is_multiple_of(2) {
        cspn = -cspn;
    }
    let mut c1 = cy[0];
    let mut c2 = y[0];
    if KODE == Scaling::Scaled {
        let mut IUF = 0;
        let NW = underflow_add_i_k(zn, &mut c1, &mut c2, &mut IUF);
        NZ += NW;
    }
    y[0] = cspn * c1 + csgn * c2;
    Ok((y, NZ))
}
