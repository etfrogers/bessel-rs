#![allow(non_snake_case, clippy::excessive_precision)]
use super::{
    IKType, Scaling, gamma_ln, i_power_series,
    overflow_checks::check_underflow_uniform_asymp_params, utils::will_underflow,
};
use crate::{
    amos::asymptotics::i_asymp_large_order,
    types::{
        BesselError::{self, *},
        BesselFloat, BesselResult, BesselValues,
    },
};

use crate::amos::{
    asymptotic_i::asymptotic_i, max_abs_component, overflow_checks::Overflow, utils::calc_rz,
};
use itertools::Either;
use num::{Complex, Zero, complex::ComplexFloat};
use std::cmp::min;

/// zbknu computes the k bessel function in the right half z plane.
/// Originally ZBKNU
pub fn k_right_half_plane<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    n: usize,
) -> Result<BesselValues<T, usize>, BesselError<T>> {
    const KMAX: usize = 30;
    let RTFRAC_PI_2: T = T::from_f64(1.25331413731550025);
    let SPI: T = T::from_f64(1.90985931710274403);
    let FPI: T = T::from_f64(1.89769999331517738);
    let CC: [T; 8] = [
        T::from_f64(5.77215664901532861e-01),
        T::from_f64(-4.20026350340952355e-02),
        T::from_f64(-4.21977345555443367e-02),
        T::from_f64(7.21894324666309954e-03),
        T::from_f64(-2.15241674114950973e-04),
        T::from_f64(-2.01348547807882387e-05),
        T::from_f64(1.13302723198169588e-06),
        T::from_f64(6.11609510448141582e-09),
    ];

    let abs_z = z.abs();
    let mut nz = 0;
    let mut underflow_occurred = false;
    let mut overflow_state;
    let rz = calc_rz(z);
    let mut integer_order = (order.round()).to_isize().unwrap(); // round to nearest int
    let simple_case = integer_order == 0 && n == 1;

    let signed_fractional_order = order - T::from_f64(integer_order as f64); // signed fractional part (-0.5 < DNU < 0.5 )
    let frac_order_sqr = if signed_fractional_order.abs() > T::MACHINE_CONSTANTS.abs_error_tolerance
    {
        signed_fractional_order * signed_fractional_order
    } else {
        T::ZERO
    };

    let (mut s1, mut s2) = if (signed_fractional_order.abs() != T::half()) && (abs_z <= T::two()) {
        // series for (z.abs() <= 2.0) and not half integer order
        let mut fc = T::one();
        let mut smu = rz.ln();
        let fmu = smu * signed_fractional_order;
        let csh = fmu.sinh();
        let cch = fmu.cosh();
        if signed_fractional_order != T::ZERO {
            fc = signed_fractional_order * T::PI();
            fc /= fc.sin();
            smu = csh / signed_fractional_order;
        }
        //-----------------------------------------------------------------------;
        //     GAM(1-Z)*GAM(1+Z)=PI*Z/SIN(PI*Z), T1=1/GAM(1-DNU), T2=1/GAM(1+DNU);
        //-----------------------------------------------------------------------;
        let t2 = (-gamma_ln(T::one() + signed_fractional_order).unwrap()).exp();
        let t1 = T::one() / (t2 * fc);

        let g1 = if signed_fractional_order.abs() <= T::from_f64(0.1) {
            //-----------------------------------------------------------------------;
            //     SERIES FOR F0 TO RESOLVE INDETERMINACY FOR SMALL ABS(DNU);
            //-----------------------------------------------------------------------;
            let mut ak = T::one();
            let mut sum = CC[0];
            for cc in CC[1..].iter() {
                ak *= frac_order_sqr;
                let tm = *cc * ak;
                sum += tm;
                if tm.abs() < T::MACHINE_CONSTANTS.abs_error_tolerance {
                    break;
                }
            }
            -sum
        } else {
            (t1 - t2) / (T::two() * signed_fractional_order)
        };
        let g2 = (t1 + t2) * T::half();
        let f = fc * (g1 * cch + g2 * smu);
        let p = T::half() * fmu.exp() / t2;
        let q = (T::half() / fmu.exp()) / t1;

        let ck = T::C_ONE;
        if simple_case {
            //-----------------------------------------------------------------------;
            //     SPECIAL CASE
            //     GENERATE K(FNU,Z), 0.0  <=  FNU  <  0.5 AND N=1;
            //-----------------------------------------------------------------------;
            let (s1, _) =
                k_right_half_plane_helper(z, frac_order_sqr, signed_fractional_order, f, p, q, ck);

            let mut y = s1;
            if scaling == Scaling::Scaled {
                y *= z.exp();
            }
            return Ok((vec![y], nz));
        }

        //-----------------------------------------------------------------------;
        //     GENERATE K(DNU,Z) AND K(DNU+1,Z) FOR FORWARD RECURRENCE;
        //-----------------------------------------------------------------------;
        let (mut s1, mut s2) =
            k_right_half_plane_helper(z, frac_order_sqr, signed_fractional_order, f, p, q, ck);

        overflow_state =
            if (order + T::one()) * smu.re.abs() > T::MACHINE_CONSTANTS.approximation_limit {
                Overflow::NearOver
            } else {
                Overflow::None
            };
        s2 *= T::MACHINE_CONSTANTS.scaling_factors[overflow_state] * rz;
        s1 *= T::MACHINE_CONSTANTS.scaling_factors[overflow_state];
        if scaling == Scaling::Scaled {
            let z_exp = z.exp();
            s1 *= z_exp;
            s2 *= z_exp;
        }
        (s1, s2)
    } else {
        // alternative series for z.abs() > 2.0 or half integer order
        //-----------------------------------------------------------------------;
        //     underflow_occured=0 MEANS NO UNDERFLOW OCCURRED;
        //     underflow_occured=1 MEANS AN UNDERFLOW OCCURRED- COMPUTATION PROCEEDS WITH;
        //     KODED=2 AND A TEST FOR ON SCALE VALUES IS MADE DURING FORWARD;
        //     RECURSION;
        //-----------------------------------------------------------------------;
        let mut coeff = Complex::<T>::new(RTFRAC_PI_2, T::ZERO) / z.sqrt();
        overflow_state = Overflow::None;
        if scaling == Scaling::Unscaled {
            if z.re > T::MACHINE_CONSTANTS.approximation_limit {
                underflow_occurred = true;
                overflow_state = Overflow::NearUnder;
            } else {
                coeff *= T::MACHINE_CONSTANTS.scaling_factors[overflow_state] * (-z).exp();
            }
        }
        let mut AK = (signed_fractional_order * T::PI()).cos().abs();
        let mut FHS = (T::from_f64(0.25) - frac_order_sqr).abs();

        if signed_fractional_order.abs() == T::half() || AK == T::ZERO || FHS == T::ZERO {
            (coeff, coeff)
        } else {
            //-----------------------------------------------------------------------
            //     MILLER ALGORITHM FOR CABS(Z) > R1;
            //-----------------------------------------------------------------------
            //-----------------------------------------------------------------------
            //     COMPUTE R2=F(E). if CABS(Z) >= R2, USE FORWARD RECURRENCE TO
            //     DETERMINE THE BACKWARD INDEX K. R2=F(E) IS A STRAIGHT LINE ON
            //     12 <= E <= 60. E IS COMPUTED FROM 2**(-E)=B**(1-i1mach(14))=
            //     TOL WHERE B IS THE BASE OF THE ARITHMETIC.
            //-----------------------------------------------------------------------
            // let f64_significant_digits =
            // (f64::MANTISSA_DIGITS - 1) as f64 * (f64::RADIX as f64).log10();
            let bits = T::from_f64((T::MANTISSA_DIGITS - 1) as f64);
            let determiner = bits.clamp(T::from_f64(12.0), T::from_f64(60.0));
            let recurrence_threshold = T::TWO_THIRDS * determiner - T::from_f64(6.0);
            let arg_z = z.arg();

            let (FK, FHS) = if abs_z > recurrence_threshold {
                //-----------------------------------------------------------------------;
                //     FORWARD RECURRENCE LOOP WHEN CABS(Z) >= R2;
                //-----------------------------------------------------------------------;
                let convergence_test =
                    AK / (T::PI() * abs_z * T::MACHINE_CONSTANTS.abs_error_tolerance);
                let mut FK = T::one();
                if convergence_test >= T::one() {
                    let mut FKS = T::two();
                    let mut CKR = abs_z + abs_z + T::two();
                    let mut p1 = T::ZERO;
                    let mut p2 = T::one();
                    let mut converged = false;
                    for _ in 0..KMAX {
                        let AK = FHS / FKS;
                        let CBR = CKR / (FK + T::one());
                        let pt = p2;
                        p2 = CBR * p2 - AK * p1;
                        p1 = pt;
                        CKR += T::two();
                        FKS += FK + FK + T::two();
                        FHS += FK + FK;
                        FK += T::one();
                        if convergence_test < p2.abs() * FK {
                            converged = true;
                            break;
                        }
                    }
                    if !converged {
                        return Err(DidNotConverge);
                    }
                    FK += SPI * arg_z * (recurrence_threshold / abs_z).sqrt();
                    FHS = (T::from_f64(0.25) - frac_order_sqr).abs();
                }
                (FK, FHS)
            } else {
                //-----------------------------------------------------------------------;
                //     COMPUTE BACKWARD INDEX K FOR CABS(Z) < R2;
                //-----------------------------------------------------------------------;
                AK *= FPI / (T::MACHINE_CONSTANTS.abs_error_tolerance * abs_z.sqrt().sqrt());
                let AA = T::from_f64(3.0) * arg_z / (T::one() + abs_z);
                let BB = T::from_f64(14.7) * arg_z / (T::from_f64(28.0) + abs_z);
                AK = (AK.ln() + abs_z * AA.cos() / (T::one() + T::from_f64(0.008) * abs_z))
                    / BB.cos();
                let FK = T::from_f64(0.12125) * AK * AK / abs_z + T::from_f64(1.5);
                (FK, FHS)
            };
            //-----------------------------------------------------------------------;
            //     BACKWARD RECURRENCE LOOP FOR MILLER ALGORITHM;
            //-----------------------------------------------------------------------;
            let K = FK.to_usize().unwrap();
            let mut k_squared = FK.floor().powi(2);
            let mut p1 = Complex::<T>::zero();
            let mut p2 = Complex::<T>::new(T::MACHINE_CONSTANTS.abs_error_tolerance, T::ZERO);
            let mut cs = p2;
            for i in (0..K).rev() {
                let k_float = T::from_usize(i + 1);
                let cb = (z + k_float) * T::two() / (k_float + T::one());
                (p1, p2) = (
                    p2,
                    (p2 * cb - p1) * (k_squared + k_float) / (k_squared - k_float + FHS),
                );
                cs += p2;
                k_squared -= (T::two() * k_float) - T::one();
            }
            //-----------------------------------------------------------------------;
            //     COMPUTE (P2/CS)=(P2/CABS(CS))*(CONJG(CS)/CABS(CS)) FOR BETTER;
            //     SCALING;
            //-----------------------------------------------------------------------;
            let mut s1 = p2 / cs.abs();
            let mut s2 = Complex::<T>::zero();
            cs = cs.conj() / cs.abs();
            s1 *= coeff * cs;
            if !simple_case {
                //-----------------------------------------------------------------------;
                //     COMPUTE P1/P2=(P1/CABS(P2)*CONJG(P2)/CABS(P2) FOR SCALING;
                //-----------------------------------------------------------------------;
                p1 /= p2.abs();
                p2 = p2.conj() / p2.abs();
                s2 = (((-(p1 * p2) + signed_fractional_order + T::half()) / z) + T::one()) * s1;
            }
            (s1, s2)
        }
    };

    // Now s1, s2 set up, we can go to recurrence

    //-----------------------------------------------------------------------
    //     FORWARD RECURSION ON THE THREE TERM RECURSION WITH RELATION WITH
    //     SCALING NEAR EXPONENT EXTREMES ON KFLAG=1 OR KFLAG=3
    //-----------------------------------------------------------------------
    let mut ck = (signed_fractional_order + T::one()) * rz;
    if n == 1 {
        integer_order -= 1
    };

    if !simple_case {
        if integer_order > 0 {
            let mut n_tested = 1;
            if underflow_occurred {
                underflow_occurred = false;
                //-----------------------------------------------------------------------;
                //     underflow_occured=1 CASES, FORWARD RECURRENCE ON SCALED VALUES ON UNDERFLOW;
                //-----------------------------------------------------------------------;
                let mut cy = [T::C_ZERO; 2];
                let half_exponent_limit = T::half() * T::MACHINE_CONSTANTS.exponent_limit;

                let abs_limit = (-T::MACHINE_CONSTANTS.exponent_limit).exp();
                let ASCLE = T::MACHINE_CONSTANTS.overflow_boundary[0];
                let mut zd = z;
                let mut IC: isize = -1;
                let mut J = 1;
                for i in 0..integer_order {
                    n_tested = i + 2;
                    // TODO same calculation as other loops - this one is over different range and sets cy
                    // (so is designed to run until cy is set, and record this in INUB)
                    (s1, s2) = (s2, s2 * ck + s1);
                    ck += rz;
                    let abs_ln_s2 = s2.abs().ln();
                    if -zd.re + abs_ln_s2 >= -T::MACHINE_CONSTANTS.exponent_limit {
                        let p1 = (-zd + s2.ln()).exp() / T::MACHINE_CONSTANTS.abs_error_tolerance;
                        if !will_underflow(p1, ASCLE, T::MACHINE_CONSTANTS.abs_error_tolerance) {
                            J = 1 - J;
                            cy[J] = p1;
                            // below implies we got here twice in a row
                            if IC == i - 1 {
                                // underflow_occurred = true; //implies 270
                                break;
                            } else {
                                IC = i;
                                continue;
                            }
                        }
                        if abs_ln_s2 < half_exponent_limit {
                            continue;
                        }
                        zd.re -= T::MACHINE_CONSTANTS.exponent_limit;
                        s1 *= abs_limit;
                        s2 *= abs_limit;
                    }
                }
                overflow_state = Overflow::NearUnder;

                s2 = cy[J];
                J = 1 - J;
                s1 = cy[J];
            }

            let mut P1R = T::MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
            let mut ASCLE = T::MACHINE_CONSTANTS.overflow_boundary[overflow_state];
            for _ in n_tested..=integer_order {
                // TODO same loop as below?
                // TODO and same recurrence logic used in ZUNKX, ZUNIX?
                (s1, s2) = (s2, ck * s2 + s1);
                ck += rz;
                if overflow_state == Overflow::NearOver {
                    continue;
                }
                let p2 = s2 * P1R;
                if max_abs_component(p2) <= ASCLE {
                    continue;
                }
                overflow_state.increment();
                ASCLE = T::MACHINE_CONSTANTS.overflow_boundary[overflow_state];
                s1 *= P1R;
                s2 = p2;
                s1 *= T::MACHINE_CONSTANTS.scaling_factors[overflow_state];
                s2 *= T::MACHINE_CONSTANTS.scaling_factors[overflow_state];
                P1R = T::MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
            }
        }
        if n == 1 {
            s1 = s2;
        }
    }

    let mut y = T::c_zeros(n);
    let n_completed = if !underflow_occurred {
        // ********* basic setup
        y[0] = s1 * T::MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
        if n > 1 {
            y[1] = s2 * T::MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
            2
        } else {
            1
        }

        // ********* End Basic Setup
    } else {
        // ******** Alternative setup if underflow_occured

        y[0] = s1;
        if n > 1 {
            y[1] = s2;
        }
        ZKSCL(
            z,
            order,
            n,
            &mut y,
            &mut nz,
            rz,
            T::MACHINE_CONSTANTS.absolute_approximation_limit,
        );
        let n_non_zero = (n - nz) as isize;
        if n_non_zero <= 0 {
            return Ok((y, nz));
        }
        let mut working_index = nz;
        s1 = y[working_index];
        y[working_index] *= T::MACHINE_CONSTANTS.reciprocal_scaling_factors[0];
        if n_non_zero > 1 {
            // if n_non_zero == 1 {
            //     return Ok((y, nz));
            // }
            working_index += 1;
            s2 = y[working_index];
            y[working_index] *= T::MACHINE_CONSTANTS.reciprocal_scaling_factors[0];
        }
        if n_non_zero > 2 {
            ck = (order + T::from_usize(working_index)) * rz;
            overflow_state = Overflow::NearUnder;
        }
        working_index + 1
    };
    // End Setup
    if n_completed >= n {
        return Ok((y, nz));
    }
    let mut P1R = T::MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
    let mut ASCLE = T::MACHINE_CONSTANTS.overflow_boundary[overflow_state];
    for y_elem in y.iter_mut().skip(n_completed) {
        // TODO same loops as above
        (s1, s2) = (s2, ck * s2 + s1);
        ck += rz;
        *y_elem = s2 * P1R;
        if overflow_state == Overflow::NearOver {
            continue;
        };
        if max_abs_component(*y_elem) <= ASCLE {
            continue;
        }
        overflow_state.increment();
        ASCLE = T::MACHINE_CONSTANTS.overflow_boundary[overflow_state];
        s1 *= P1R;
        s2 = *y_elem;
        s1 *= T::MACHINE_CONSTANTS.scaling_factors[overflow_state];
        s2 *= T::MACHINE_CONSTANTS.scaling_factors[overflow_state];
        P1R = T::MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
    }
    Ok((y, nz))
}

fn k_right_half_plane_helper<T: BesselFloat>(
    z: Complex<T>,
    frac_order_sqr: T,
    signed_fractional_order: T,
    mut f: Complex<T>,
    mut p: Complex<T>,
    mut q: Complex<T>,
    mut ck: Complex<T>,
) -> (Complex<T>, Complex<T>) {
    let mut a1 = T::one();
    let cz_sqr_over_4 = T::from_f64(0.25) * z.powu(2);
    let abs_z = z.abs();
    let abs_z_sqr_over_4 = T::from_f64(0.25) * abs_z * abs_z;
    let mut ak = T::one();
    let mut bk = T::one() - frac_order_sqr;

    let mut s1 = f;
    let mut s2 = p;
    if abs_z >= T::MACHINE_CONSTANTS.abs_error_tolerance {
        while a1 > T::MACHINE_CONSTANTS.abs_error_tolerance {
            f = (f * ak + p + q) / bk;
            p /= ak - signed_fractional_order;
            q /= ak + signed_fractional_order;
            ck *= cz_sqr_over_4 / ak;
            s1 += ck * f;
            s2 += ck * (p - ak * f);
            a1 *= abs_z_sqr_over_4 / ak;
            bk += (T::two() * ak) + T::one();
            ak += T::one();
        }
    }
    (s1, s2)
}

/// Set k functions to zero on underflow, continue recurrence
/// on scaled functions until two members come on scale, then
/// return with min(nz+2,n) values scaled by 1/tol.
///
/// Originally ZKSCL
fn ZKSCL<T: BesselFloat>(
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

/// ratios_i computes ratios of I bessel functions by backward
/// recurrence. The starting index is determined by forward
/// recurrence as described in J. Res. of Nat. Bur. of Standards-B,
/// Mathematical Sciences, vol 77b, p111-114, September, 1973,
/// Bessel functions I and J of complex argument and integer order,
/// by D. J. Sookne.
///
/// Originally ZRATI
fn ratios_i<T: BesselFloat>(z: Complex<T>, order: T, n: usize) -> Vec<Complex<T>> {
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

// i_wronksian computes the i bessel function for re(z) >= 0.0 by
// normalizing the i function ratios from zrati by the Wronskian
// Originally ZWRSK
fn i_wronksian<T: BesselFloat>(
    zr: Complex<T>,
    order: T,
    scaling: Scaling,
    n: usize,
    y: &mut [Complex<T>],
) -> Result<usize, BesselError<T>> {
    //-----------------------------------------------------------------------
    //     I(FNU+I-1,Z) BY BACKWARD RECURRENCE FOR RATIOS
    //     Y(I)=I(FNU+I,Z)/I(FNU+I-1,Z) FROM CRATI NORMALIZED BY THE
    //     WRONSKIAN WITH K(FNU,Z) AND K(FNU+1,Z) FROM CBKNU.
    //-----------------------------------------------------------------------
    let nz = 0;
    let (cw, _) = k_right_half_plane(zr, order, scaling, 2)?;
    let y_ratios = ratios_i(zr, order, n);
    //-----------------------------------------------------------------------
    //     RECUR FORWARD ON I(FNU+1,Z) = R(FNU,Z)*I(FNU,Z),
    //     R(FNU+J-1,Z)=Y(J),  J=1,...,N
    //-----------------------------------------------------------------------
    let mut cinu = T::C_ONE;
    if scaling == Scaling::Scaled {
        cinu = Complex::<T>::cis(zr.im);
    }
    //-----------------------------------------------------------------------
    //     ON LOW EXPONENT MACHINES THE K FUNCTIONS CAN BE CLOSE TO BOTH
    //     THE UNDER AND OVERFLOW LIMITS AND THE NORMALIZATION MUST BE
    //     SCALED TO PREVENT OVER OR UNDERFLOW. CUOIK HAS DETERMINED THAT
    //     THE RESULT IS ON SCALE.
    //-----------------------------------------------------------------------
    let acw = cw[1].abs();
    let CSCLR = if acw <= T::MACHINE_CONSTANTS.absolute_approximation_limit {
        T::one() / T::MACHINE_CONSTANTS.abs_error_tolerance
    } else if acw >= T::one() / T::MACHINE_CONSTANTS.absolute_approximation_limit {
        T::MACHINE_CONSTANTS.abs_error_tolerance
    } else {
        T::one()
    };
    let c1 = cw[0] * CSCLR;
    let c2 = cw[1] * CSCLR;
    //-----------------------------------------------------------------------
    //     CINU=CINU*(CONJG(CT)/CABS(CT))*(1.0/CABS(CT) PREVENTS
    //     UNDER- OR OVERFLOW PREMATURELY BY SQUARING CABS(CT)
    //-----------------------------------------------------------------------
    let mut ct = zr * (y_ratios[0] * c1 + c2);
    let ct_abs = ct.abs();
    ct = ct.conj() / ct_abs;
    cinu = (cinu / ct_abs) * ct;
    y[0] = cinu * CSCLR;
    for i in 1..n {
        cinu *= y_ratios[i - 1];
        y[i] = cinu * CSCLR;
    }
    Ok(nz)
}

/// i_right_half_plane computes the i function in the right half z plane
/// Originally ZBINU
pub fn i_right_half_plane<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    KODE: Scaling,
    N: usize,
) -> BesselResult<T, usize> {
    let mut NZ = 0;
    let AZ = z.abs();
    let mut NN: usize = N;
    let mut DFNU = order + T::from_usize(N - 1);
    let mut cy = T::c_zeros(N);
    if AZ <= T::two() || AZ * AZ * T::from_f64(0.25) <= DFNU + T::one() {
        //-----------------------------------------------------------------------
        //     POWER SERIES
        //-----------------------------------------------------------------------
        let NW;
        (cy, NW) = i_power_series(z, order, KODE, NN)?;
        let INW: usize = NW.unsigned_abs();
        NZ += INW;
        NN -= INW;
        if NN == 0 || NW >= 0 {
            return Ok((cy, NZ));
        }

        DFNU = order + (T::from_usize(NN) - T::one());
    }

    if (AZ >= T::MACHINE_CONSTANTS.asymptotic_z_limit)
        && ((DFNU <= T::one()) || (AZ + AZ >= DFNU * DFNU))
    {
        //-----------------------------------------------------------------------
        //     ASYMPTOTIC EXPANSION FOR LARGE Z
        //-----------------------------------------------------------------------
        let (cy, nw) = asymptotic_i(z, order, KODE, NN)?;
        debug_assert!(nw == NZ);
        return Ok((cy, NZ));
    }
    let mut skip_az_rl_check = true;
    if DFNU > T::one() {
        skip_az_rl_check = false;
        //-----------------------------------------------------------------------
        //     OVERFLOW AND UNDERFLOW TEST ON I SEQUENCE FOR MILLER ALGORITHM
        //-----------------------------------------------------------------------
        let nw = check_underflow_uniform_asymp_params(z, order, KODE, IKType::I, NN, &mut cy)?;
        NZ += nw;
        NN -= nw;
        if NN == 0 {
            return Ok((cy, NZ));
        }
        DFNU = order + T::from_usize(NN - 1);
    }
    if (DFNU > T::MACHINE_CONSTANTS.asymptotic_order_limit)
        || (AZ > T::MACHINE_CONSTANTS.asymptotic_order_limit)
    {
        //-----------------------------------------------------------------------
        //     INCREMENT FNU+NN-1 UP TO FNUL, COMPUTE AND RECUR BACKWARD
        //-----------------------------------------------------------------------
        let NUI = ((T::MACHINE_CONSTANTS.asymptotic_order_limit - DFNU).trunc() + T::one())
            .max(T::zero())
            .to_usize()
            .unwrap();

        let (NW, NLAST) = i_asymp_large_order(z, order, KODE, NN, NUI, &mut cy)?;
        NZ += NW;
        if NLAST == 0 {
            return Ok((cy, NZ));
        }
        NN = NLAST;
    }
    if !skip_az_rl_check && AZ <= T::MACHINE_CONSTANTS.asymptotic_z_limit {
        //-----------------------------------------------------------------------
        //     MILLER ALGORITHM NORMALIZED BY THE SERIES
        //-----------------------------------------------------------------------
        let (cy, _) = i_miller(z, order, KODE, NN)?;
        return Ok((cy, NZ)); //}
    }
    //-----------------------------------------------------------------------
    //     MILLER ALGORITHM NORMALIZED BY THE WRONSKIAN
    //-----------------------------------------------------------------------
    //-----------------------------------------------------------------------
    //     OVERFLOW TEST ON K FUNCTIONS USED IN WRONSKIAN
    //-----------------------------------------------------------------------
    if let Ok(NW) =
        check_underflow_uniform_asymp_params(z, order, KODE, IKType::K, 2, &mut [T::C_ONE; 2])
    {
        if NW > 0 {
            Err(Overflow)
        } else {
            let nz = i_wronksian(z, order, KODE, NN, &mut cy)?;
            Ok((cy, nz))
        }
    } else {
        Ok((vec![T::C_ONE; NN], NN))
    }
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn recurr<T: BesselFloat>(
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
