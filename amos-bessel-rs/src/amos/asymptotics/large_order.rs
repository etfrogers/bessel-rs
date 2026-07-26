#![allow(non_snake_case)]
use num::{Complex, complex::ComplexFloat};

use crate::{
    BesselError, Scaling,
    amos::{
        RotationDirection,
        asymptotics::{i_uniform_asymp1, i_uniform_asymp2, k_uniform_asymp1, k_uniform_asymp2},
        limits::Overflow,
        max_abs_component,
        utils::{calc_rz, imaginary_dominant},
    },
    types::{BesselFloat, BesselResult},
};

/// i_asymp_large_order computes the i bessel function for large cabs(z) >
/// asymptotic_order_limit and fnu+n-1 < asymptotic_order_limit.
/// The order is increased from
/// fnu+n-1 greater than asymptotic_order_limit by adding nui and computing
/// according to the uniform asymptotic expansion for J(fnu,z)
///  if z is imaginary_dominant and the expansion for I(fnu,z)
/// if z is _not_ imaginary_dominant.
///
/// Originally ZBUNI
pub(crate) fn i_asymp_large_order<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    KODE: Scaling,
    n: usize,
    NUI: usize,
    y: &mut [Complex<T>],
) -> Result<(usize, usize), BesselError<T>> {
    let imaginary_dominant = imaginary_dominant(z);
    if NUI != 0 {
        let mut FNUI = T::from_usize(NUI);
        let DFNU = order + T::from_usize(n - 1);
        let GNU = DFNU + FNUI;
        let mut cy = [T::C_ZERO; 2];
        let (NW, NLAST) = if imaginary_dominant {
            //-----------------------------------------------------------------------
            //     ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*FRAC_PI_2)) FOR LARGE FNU
            //     APPLIED IN PI/3 < ABS(ARG(Z)) <= PI/2 WHERE M=+I OR -I
            //     AND FRAC_PI_2=PI/2
            //-----------------------------------------------------------------------
            i_uniform_asymp2(z, GNU, KODE, 2, &mut cy)?
        } else {
            //-----------------------------------------------------------------------
            //     ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN
            //     -PI/3 <= ARG(Z) <= PI/3
            //-----------------------------------------------------------------------
            i_uniform_asymp1(z, GNU, KODE, 2, &mut cy)?
        };
        if NW != 0 {
            return Ok((0, n));
        }
        //----------------------------------------------------------------------
        //     SCALE BACKWARD RECURRENCE, BRY(3) IS DEFINED BUT NEVER USED
        //----------------------------------------------------------------------
        let (mut overflow_state, mut ASCLE, mut CSCLR) =
            if cy[0].abs() <= T::MACHINE_CONSTANTS.absolute_approximation_limit {
                (
                    Overflow::NearUnder,
                    T::MACHINE_CONSTANTS.absolute_approximation_limit,
                    T::one() / T::MACHINE_CONSTANTS.abs_error_tolerance,
                )
            } else if cy[0].abs()
                >= T::from_f64(1.0) / T::MACHINE_CONSTANTS.absolute_approximation_limit
            {
                (
                    Overflow::NearOver,
                    T::max_value() / T::from_f64(2.0),
                    T::MACHINE_CONSTANTS.abs_error_tolerance,
                )
            } else {
                (
                    Overflow::None,
                    T::from_f64(1.0) / T::MACHINE_CONSTANTS.absolute_approximation_limit,
                    T::one(),
                )
            };

        let mut CSCRR = T::one() / CSCLR;
        let mut s1 = cy[1] * CSCLR;
        let mut s2 = cy[0] * CSCLR;
        // working out rz in multiple steps seems to give different floating point answer.
        let rz = calc_rz(z);

        for _ in 0..NUI {
            let st = s2;
            s2 = (DFNU + FNUI) * rz * s2 + s1;
            s1 = st;
            FNUI -= T::one();
            if overflow_state == Overflow::NearOver {
                continue;
            }
            let st = s2 * CSCRR;
            if max_abs_component(st) <= ASCLE {
                continue;
            }
            overflow_state.increment();
            ASCLE = overflow_state.boundary::<T>();
            s1 *= CSCRR;
            s2 = st;
            CSCLR *= T::MACHINE_CONSTANTS.abs_error_tolerance;
            CSCRR = T::one() / CSCLR;
            s1 *= CSCLR;
            s2 *= CSCLR;
        }
        y[n - 1] = s2 * CSCRR;
        if n == 1 {
            return Ok((0, NLAST));
        }
        let NL = n - 1;
        FNUI = T::from_usize(NL);
        let mut K = NL;
        for _ in 0..NL {
            let st = s2;
            s2 = (order + FNUI) * (rz * s2) + s1;
            s1 = st;
            y[K - 1] = s2 * CSCRR;
            FNUI -= T::one();
            K -= 1;
            if overflow_state == Overflow::NearOver {
                continue;
            }
            // using K (rather than K-1) below as Amos "saved" the y value before K was decremented
            if max_abs_component(y[K]) <= ASCLE {
                continue;
            }
            overflow_state.increment();
            ASCLE = overflow_state.boundary::<T>();
            s1 *= CSCRR;
            s2 = y[K - 1];
            CSCLR *= T::MACHINE_CONSTANTS.abs_error_tolerance;
            CSCRR = T::one() / CSCLR;
            s1 *= CSCLR;
            s2 *= CSCLR;
        }
        return Ok((0, NLAST));
    }
    let (NW, NLAST) = if imaginary_dominant {
        //-----------------------------------------------------------------------
        //     ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*FRAC_PI_2)) FOR LARGE FNU
        //     APPLIED IN PI/3 < ABS(ARG(Z)) <= PI/2 WHERE M=+I OR -I
        //     AND FRAC_PI_2=PI/2
        //-----------------------------------------------------------------------
        i_uniform_asymp2(z, order, KODE, n, y)?
    } else {
        //-----------------------------------------------------------------------
        //     ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN
        //     -PI/3 <= ARG(Z) <= PI/3
        //-----------------------------------------------------------------------
        i_uniform_asymp1(z, order, KODE, n, y)?
    };

    Ok((NW, NLAST))
}

/// zbunk computes the K Bessel function for order > asymptotic_order_limit.
/// according to the uniform asymptotic expansion for K(order, z)
/// in zunk1 and the expansion for H(2, order, z) in zunk2
pub fn k_asymp_large_order<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    rotation: RotationDirection,
    n: usize,
) -> BesselResult<T> {
    if imaginary_dominant(z) {
        //-----------------------------------------------------------------------
        //     ASYMPTOTIC EXPANSION FOR H(2,FNU,Z*EXP(M*FRAC_PI_2)) FOR LARGE FNU
        //     APPLIED IN PI/3 < ABS(ARG(Z)) <= PI/2 WHERE M=+I OR -I
        //     AND FRAC_PI_2=PI/2
        //-----------------------------------------------------------------------
        k_uniform_asymp2(z, order, scaling, rotation, n)
    } else {
        //-----------------------------------------------------------------------
        //     ASYMPTOTIC EXPANSION FOR K(FNU,Z) FOR LARGE FNU APPLIED IN
        //     -PI/3 <= ARG(Z) <= PI/3
        //-----------------------------------------------------------------------
        k_uniform_asymp1(z, order, scaling, rotation, n)
    }
}
