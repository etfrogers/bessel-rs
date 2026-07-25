#![allow(non_snake_case)]
use std::cmp::min;

use itertools::Either;
use num::{Complex, Zero, complex::ComplexFloat};

use crate::{
    BesselError::DidNotConverge,
    Scaling,
    amos::{
        gamma_ln,
        limits::Overflow,
        max_abs_component,
        utils::{calc_rz, will_underflow},
    },
    types::{BesselFloat, BesselResult},
};

/// i_miller computes the i bessel function for re(z) >= 0.0 by the
/// Miller algorithm normalized by a Neumann series.
/// Originally ZMLRI
pub(crate) fn i_miller<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    n: usize,
) -> BesselResult<T, usize> {
    let scale: T = T::two() * T::MIN_POSITIVE / T::MACHINE_CONSTANTS.abs_error_tolerance;
    let nz = 0;
    let abs_z = z.abs();
    let int_abs_z = abs_z.to_usize().unwrap();
    let int_order = order.to_usize().unwrap();
    let INU = int_order + n - 1;
    let AT = T::from_f64(int_abs_z as f64) + T::one();
    let RAZ = T::one() / abs_z;
    let mut ck = z.conj() * RAZ * RAZ * AT;
    let rz = calc_rz(z);
    let mut p1 = T::C_ZERO;
    let mut p2 = T::C_ONE;
    let mut ACK = (AT + T::one()) * RAZ;
    let RHO = ACK + (ACK * ACK - T::one()).sqrt();
    let RHO2 = RHO * RHO;
    let mut TST = (RHO2 + RHO2) / ((RHO2 - T::one()) * (RHO - T::one()));
    TST /= T::MACHINE_CONSTANTS.abs_error_tolerance;
    //-----------------------------------------------------------------------
    //     COMPUTE RELATIVE TRUNCATION ERROR INDEX FOR SERIES
    //-----------------------------------------------------------------------
    let mut AK = AT;
    let mut converged = false;
    let mut I = 0;
    for i in 0..80 {
        I = i + 1;
        let pt = p2;
        p2 = p1 - ck * p2;
        p1 = pt;
        ck += rz;
        if p2.abs() > TST * AK * AK {
            converged = true;
            break;
        }
        AK += T::one();
    }
    if !converged {
        return Err(DidNotConverge);
    }
    I += 1;
    let mut K = 0;
    if INU >= int_abs_z {
        //-----------------------------------------------------------------------
        //     COMPUTE RELATIVE TRUNCATION ERROR FOR RATIOS
        //-----------------------------------------------------------------------
        p1 = T::C_ZERO;
        p2 = T::C_ONE;
        let AT = T::from_f64(INU as f64) + T::one();
        ck = z.conj() * RAZ * RAZ * AT;
        ACK = AT * RAZ;
        TST = (ACK / T::MACHINE_CONSTANTS.abs_error_tolerance).sqrt();
        let mut hit_loop_end = false;
        converged = false;
        for k in 0..80 {
            K = k;
            let pt = p2;
            p2 = p1 - ck * pt;
            p1 = pt;
            ck += rz;
            let AP = p2.abs();
            if AP < TST {
                continue;
            }
            if hit_loop_end {
                converged = true;
                break;
            }
            ACK = ck.abs();
            let FLAM = ACK + (ACK * ACK - T::one()).sqrt();
            let FKAP = AP / p1.abs();
            let RHO = FLAM.min(FKAP);
            TST *= (RHO / (RHO * RHO - T::one())).sqrt();
            hit_loop_end = true;
        }
        if !converged {
            return Err(DidNotConverge);
        }
    }
    //-----------------------------------------------------------------------
    //     BACKWARD RECURRENCE AND SUM NORMALIZING RELATION
    //-----------------------------------------------------------------------
    K += 1;
    let KK = (I + int_abs_z).max(K + INU);
    let mut kk_float = T::from_f64(KK as f64);
    let mut p1 = T::C_ZERO;
    //-----------------------------------------------------------------------
    //     SCALE P2 AND SUM BY SCLE
    //-----------------------------------------------------------------------
    let mut p2 = Complex::<T>::new(scale, T::ZERO);
    let fractional_order = order.fract();
    let twice_fractional_order = fractional_order + fractional_order;
    let mut BK = (gamma_ln(kk_float + twice_fractional_order + T::one()).unwrap()
        - gamma_ln(kk_float + T::one()).unwrap()
        - gamma_ln(twice_fractional_order + T::one()).unwrap())
    .exp();
    let mut sumr = T::C_ZERO;
    for _ in 0..(KK - INU) {
        let pt = p2;
        p2 = p1 + (kk_float + fractional_order) * (rz * p2);
        p1 = pt;
        AK = T::one() - twice_fractional_order / (kk_float + twice_fractional_order);
        ACK = BK * AK;
        sumr += (ACK + BK) * p1;
        BK = ACK;
        kk_float -= T::one();
    }
    let mut y = T::c_zeros(n);
    y[n - 1] = p2;
    if n != 1 {
        for i in 1..n {
            let pt = p2;
            p2 = p1 + (kk_float + fractional_order) * (rz * pt);
            p1 = pt;
            AK = T::one() - twice_fractional_order / (kk_float + twice_fractional_order);
            ACK = BK * AK;
            sumr += (ACK + BK) * p1;
            BK = ACK;
            kk_float -= T::one();
            y[n - (i + 1)] = p2;
        }
    }
    if int_order > 0 {
        for _i in 0..int_order {
            let pt = p2;
            p2 = p1 + (kk_float + fractional_order) * (rz * pt);
            p1 = pt;
            AK = T::one() - twice_fractional_order / (kk_float + twice_fractional_order);
            ACK = BK * AK;
            sumr += (ACK + BK) * p1;
            BK = ACK;
            kk_float -= T::one();
        }
    }

    let mut pt = z;
    if scaling == Scaling::Scaled {
        pt.re = T::ZERO;
    }
    p1 = -fractional_order * rz.ln() + pt;
    let AP = gamma_ln(T::one() + fractional_order).unwrap();
    p1 -= AP;
    //-----------------------------------------------------------------------
    //     THE DIVISION CEXP(PT)/(SUM+P2) IS ALTERED TO AVOID OVERFLOW
    //     IN THE DENOMINATOR BY SQUARING LARGE QUANTITIES
    //-----------------------------------------------------------------------
    p2 += sumr;
    let AP = p2.abs();
    ck = p1.exp() / AP;
    let cnorm = ck * p2.conj() / AP;
    for element in y.iter_mut() {
        *element *= cnorm;
    }
    Ok((y, nz))
}

/// i_ratios computes ratios of I bessel functions by backward
/// recurrence. The starting index is determined by forward
/// recurrence as described in J. Res. of Nat. Bur. of Standards-B,
/// Mathematical Sciences, vol 77b, p111-114, September, 1973,
/// Bessel functions I and J of complex argument and integer order,
/// by D. J. Sookne.
///
/// Originally ZRATI
pub(crate) fn i_ratios<T: BesselFloat>(z: Complex<T>, order: T, n: usize) -> Vec<Complex<T>> {
    let abs_z = z.abs();
    let integer_order = order.to_usize().unwrap();
    let modified_int_order = integer_order + n - 1;
    let int_abs_z = abs_z.to_isize().unwrap();
    let FNUP = T::from_f64((int_abs_z + 1).max(modified_int_order as isize) as f64);
    let ID_ = modified_int_order as isize - int_abs_z - 1;
    let ID = if ID_ > 0 { 0 } else { ID_ };

    let rz = calc_rz(z);
    let mut K = 1;
    let mut abs_p2;
    {
        let mut t1 = rz * FNUP;
        let mut p2 = -t1;
        let mut p1 = T::C_ONE;
        t1 += rz;

        abs_p2 = p2.abs();
        let mut abs_p1 = p1.abs();
        //-----------------------------------------------------------------------
        //     THE OVERFLOW TEST ON K(FNU+I-1,Z) BEFORE THE CALL TO CBKNU
        //     GUARANTEES THAT P2 IS ON SCALE. SCALE TEST1 AND ALL SUBSEQUENT
        //     P2 VALUES BY AP1 TO ENSURE THAT AN OVERFLOW DOES NOT OCCUR
        //     PREMATURELY.
        //-----------------------------------------------------------------------
        let ARG = (abs_p2 + abs_p2) / (abs_p1 * T::MACHINE_CONSTANTS.abs_error_tolerance);
        let TEST1 = ARG.sqrt();
        let mut TEST = TEST1;
        p1 /= abs_p1;
        p2 /= abs_p1;
        abs_p2 /= abs_p1;
        let mut first_pass = true;
        'l10: loop {
            K += 1;
            abs_p1 = abs_p2;
            (p1, p2) = (p2, p1 - (t1 * p2));
            t1 += rz;
            abs_p2 = p2.abs();
            if abs_p1 <= TEST {
                continue;
            }
            if !first_pass {
                break 'l10;
            }
            {
                let ak = t1.abs() / T::two();
                let flam = ak + (ak.powi(2) - T::one()).sqrt();
                let rho = abs_p2 / abs_p1.min(flam);
                TEST = TEST1 * (rho / (rho.powi(2) - T::one())).sqrt();
            }
            first_pass = false;
        }
    }

    let mut p1 = Complex::<T>::new(T::one() / abs_p2, T::ZERO);
    let mut p2 = T::C_ZERO;

    {
        let kk: usize = (K as isize + 1 - ID).try_into().unwrap();
        let mut t1 = Complex::<T>::from(T::from_usize(kk));
        let modified_order = order + T::from_usize(n - 1);
        for _ in 0..kk {
            (p1, p2) = (p1 * (rz * (modified_order + t1.re)) + p2, p1);
            t1.re -= T::one();
        }
        if p1.re == T::ZERO && p1.im == T::ZERO {
            p1 = Complex::<T>::new(
                T::MACHINE_CONSTANTS.abs_error_tolerance,
                T::MACHINE_CONSTANTS.abs_error_tolerance,
            );
        }
    }
    let mut cy = T::c_zeros(n);
    cy[n - 1] = p2 / p1;
    if n > 1 {
        let mut t1 = Complex::<T>::from(T::from_usize(n - 1));
        let cdfnu = order * rz;
        for k in (1..n).rev() {
            let mut pt = cdfnu + t1 * rz + cy[k];
            let mut abs_pt = pt.abs();
            if abs_pt == T::ZERO {
                pt = Complex::<T>::new(
                    T::MACHINE_CONSTANTS.abs_error_tolerance,
                    T::MACHINE_CONSTANTS.abs_error_tolerance,
                );
                abs_pt = pt.abs();
            }
            cy[k - 1] = pt.conj() / abs_pt.powi(2);
            t1 -= T::one();
        }
    }
    cy
}

/// Set k functions to zero on underflow, continue recurrence
/// on scaled functions until two members come on scale, then
/// return with min(nz+2,n) values scaled by 1/tol.
///
/// Originally ZKSCL
pub(crate) fn scale_k_recurrence<T: BesselFloat>(
    zr: Complex<T>,
    order: T,
    n: usize,
    y: &mut [Complex<T>],
    nz: &mut usize,
    rz: Complex<T>,
    absolute_approximation_limit: T,
) {
    *nz = 0;
    // let NN = min(2, n);
    let mut cy = [T::C_ZERO; 2];
    let mut i_completed = 0;
    // repeats twice, unless n < 2
    for i in 0..min(2, n) {
        let s1 = y[i];
        cy[i] = s1;
        *nz += 1;
        y[i] = T::C_ZERO;
        if -zr.re + s1.abs().ln() < -T::MACHINE_CONSTANTS.exponent_limit {
            continue;
        }

        let cs = (s1.ln() - zr).exp() / T::MACHINE_CONSTANTS.abs_error_tolerance;
        if will_underflow(
            cs,
            absolute_approximation_limit,
            T::MACHINE_CONSTANTS.abs_error_tolerance,
        ) {
            continue;
        }
        y[i] = cs;
        i_completed = i;
        *nz -= 1;
    }
    if n <= 2 || *nz == 0 {
        return;
    }
    // if i_completed < 1 {
    //     y[0] = T::C_ZERO;
    //     *nz = 2;
    // }
    // if n == 2 {
    //     return;
    // }
    // if *nz == 0 {
    //     return;
    // }
    let FN = order + T::one();
    let mut ck = FN * rz;
    let mut s1 = cy[0];
    let mut s2 = cy[1];
    let half_elim = T::half() * T::MACHINE_CONSTANTS.exponent_limit;
    let ELM = (-T::MACHINE_CONSTANTS.exponent_limit).exp();
    let CELMR = ELM;
    let mut zd = zr;
    //     FIND TWO CONSECUTIVE Y VALUES ON SCALE. SCALE RECURRENCE if
    //     S2 GETS LARGER THAN EXP(ELIM/2)
    let mut skip_to_40 = false;
    let mut I = 0;
    for (i, yi) in y.iter_mut().enumerate().skip(2) {
        I = i;
        let mut cs = s2;
        s2 = cs * ck + s1;
        s1 = cs;
        ck += rz;
        let ALAS = s2.abs().ln();
        *nz += 1;
        *yi = Complex::<T>::zero();
        if -zd.re + s2.abs().ln() >= -T::MACHINE_CONSTANTS.exponent_limit {
            cs = s2.ln() - zd;
            cs = cs.exp() / T::MACHINE_CONSTANTS.abs_error_tolerance;
            if !will_underflow(
                cs,
                absolute_approximation_limit,
                T::MACHINE_CONSTANTS.abs_error_tolerance,
            ) {
                *yi = cs;
                *nz -= 1;
                if i_completed == i - 1 {
                    skip_to_40 = true;
                    break;
                }
                i_completed = i;
                continue;
            }
        }

        if ALAS < half_elim {
            continue;
        }
        zd -= T::MACHINE_CONSTANTS.exponent_limit;
        s1 *= CELMR;
        s2 *= CELMR;
    }
    if !skip_to_40 {
        *nz = n;
        if i_completed == n {
            *nz = n - 1
        };
    } else {
        *nz = I - 2;
    }
    for element in y.iter_mut().take(*nz) {
        *element = T::C_ZERO;
    }
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn backward_recurrence<T: BesselFloat>(
    forward: bool,
    order: T,
    z: Complex<T>,
    y: &mut [Complex<T>],
    n_offset: usize,
    mut s1: Complex<T>,
    mut s2: Complex<T>,
    mut overflow_state: Overflow,
) {
    let rz = calc_rz(z);

    let base_iterator = y.iter_mut().enumerate();
    let iterator = if forward {
        Either::Right(base_iterator.skip(n_offset))
    } else {
        Either::Left(base_iterator.take(n_offset).rev())
    };
    let index_adjustment = if forward { -T::one() } else { T::one() };

    let mut recip_scale_factor = T::MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
    let mut boundary = T::MACHINE_CONSTANTS.overflow_boundary[overflow_state];

    for (i, yi) in iterator {
        let modified_order = order + T::from_usize(i) + index_adjustment;
        (s1, s2) = (s2, s1 + modified_order * rz * s2);
        *yi = s2 * recip_scale_factor;
        if overflow_state != Overflow::NearOver && max_abs_component(*yi) > boundary {
            overflow_state.increment();
            boundary = T::MACHINE_CONSTANTS.overflow_boundary[overflow_state];
            s1 *= recip_scale_factor;
            s2 = *yi;
            s1 *= T::MACHINE_CONSTANTS.scaling_factors[overflow_state];
            s2 *= T::MACHINE_CONSTANTS.scaling_factors[overflow_state];
            recip_scale_factor = T::MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
        }
    }
}
