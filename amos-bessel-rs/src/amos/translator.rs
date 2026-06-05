#![allow(non_snake_case, clippy::excessive_precision)]
use super::{
    IKType, Scaling, gamma_ln, i_power_series,
    overflow_checks::{check_underflow_uniform_asymp_params, underflow_add_i_k, zunik},
    utils::{AIC, will_underflow},
};
use crate::types::{
    BesselError::{self, *},
    BesselFloat, BesselResult, BesselValues,
};

use crate::amos::{
    CIP, RotationDirection,
    asymptotic_i::asymptotic_i,
    complex_airy, max_abs_component,
    overflow_checks::{Overflow, zunhj},
    utils::{calc_rz, imaginary_dominant},
};
use itertools::Either;
use num::{Complex, Integer, Zero, complex::ComplexFloat};
use std::cmp::min;

pub fn airy_power_series<T: BesselFloat>(
    z: Complex<T>,
    return_derivative: bool,
    coeffs: (f64, f64),
) -> Complex<T> {
    let c1 = T::from_f64(coeffs.0);
    let c2 = T::from_f64(coeffs.1);
    let float_is_derivative = if return_derivative { T::ONE } else { T::ZERO };

    let abs_z = z.abs();
    let z_floor = if abs_z < T::MACHINE_CONSTANTS.underflow_limit {
        T::C_ZERO
    } else {
        z
    };
    let (s1, s2) = if abs_z < T::MACHINE_CONSTANTS.abs_error_tolerance {
        (T::C_ONE, T::C_ONE)
    } else {
        let abs_z_sq = abs_z * abs_z;
        let mut s1 = T::C_ONE;
        let mut s2 = T::C_ONE;

        if abs_z_sq >= T::MACHINE_CONSTANTS.abs_error_tolerance / abs_z {
            let mut term1 = T::C_ONE;
            let mut term2 = T::C_ONE;
            let mut a_term = T::ONE;
            let z3 = z.powf(T::from_f64(3.0));
            let AZ3 = abs_z * abs_z_sq;
            let (AK, BK, CK, DK) = (
                T::from_f64(2.0) + float_is_derivative,
                T::from_f64(3.0) - T::two() * float_is_derivative,
                T::from_f64(4.0) - float_is_derivative,
                T::from_f64(3.0) + T::two() * float_is_derivative,
            );
            let mut D1: T = AK * DK;
            let mut D2 = BK * CK;
            let mut AD = D1.min(D2);
            let mut AK = T::from_f64(24.0) + T::from_f64(9.0) * float_is_derivative;
            let mut BK = T::from_f64(30.0) - T::from_f64(9.0) * float_is_derivative;
            for _ in 0..25 {
                term1 = term1 * z3 / D1;
                s1 += term1;
                term2 = term2 * z3 / D2;
                s2 += term2;
                a_term = a_term * AZ3 / AD;
                D1 += AK;
                D2 += BK;
                AD = D1.min(D2);
                if a_term < T::MACHINE_CONSTANTS.abs_error_tolerance * AD {
                    break;
                }
                AK += T::from_f64(18.0);
                BK += T::from_f64(18.0);
            }
        }
        (s1, s2)
    };

    if return_derivative {
        z_floor.powf(T::two()) * s1 * (c1 / T::two()) - s2 * c2
    } else {
        s1 * c1 - c2 * z_floor * s2
    }
}

pub fn airy_pair<T: BesselFloat>(z: Complex<T>) -> (Complex<T>, Complex<T>) {
    //note that ZAIRY calls in fortran code ignore IERR (using IDUM)
    let airy = match complex_airy(z, false, Scaling::Scaled) {
        Ok((y, _)) => y,
        Err(PartialLossOfSignificance { y, nz: _ }) => y[0],
        // If loss of significance, Fortran code would continue with un-initialised y,
        // which is usually ~=0. Also long as it is << d_airy, the logic below means
        // it will not matter what the precise value is
        Err(LossOfSignificance) => T::C_ZERO,
        Err(err) => {
            panic!(
                "An error {:?} was generated, which is not handled by the Amos code",
                err
            )
        }
    };
    let d_airy = match complex_airy(z, true, Scaling::Scaled) {
        Ok((y, _)) => y,
        Err(PartialLossOfSignificance { y, nz: _ }) => y[0],
        // If loss of significance, Fortran code would continue with un-initialised y,
        // which is usually ~=0. Also long as it is << a_airy, the logic below means
        // it will not matter what the precise value is
        Err(LossOfSignificance) => T::C_ZERO,
        Err(err) => {
            panic!(
                "An error {:?} was generated, which is not handled by the Amos code",
                err
            )
        }
    };
    (airy, d_airy)
}

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

/// zbunk computes the K Bessel function for order > asymptotic_order_limit.
/// according to the uniform asymptotic expansion for K(order, z)
/// in zunk1 and the expansion for H(2, order, z) in zunk2
pub fn ZBUNK<T: BesselFloat>(
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
        ZUNK2(z, order, scaling, rotation, n)
    } else {
        //-----------------------------------------------------------------------
        //     ASYMPTOTIC EXPANSION FOR K(FNU,Z) FOR LARGE FNU APPLIED IN
        //     -PI/3 <= ARG(Z) <= PI/3
        //-----------------------------------------------------------------------
        ZUNK1(z, order, scaling, rotation, n)
    }
}

/// i_miller computes the i bessel function for re(z) >= 0.0 by the
/// Miller algorithm normalized by a Neumann series.
/// Originally ZMLRI
fn i_miller<T: BesselFloat>(
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
        return Err(Overflow);
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

        let (NW, NLAST) = ZBUNI(z, order, KODE, NN, NUI, &mut cy)?;
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
pub fn ZACAI<T: BesselFloat>(
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
        return Err(Overflow);
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

/// zunk1 computes K(fnu,z) and its analytic continuation from the
/// right half plane to the left half plane by means of the
/// uniform asymptotic expansion.
/// `rotation` indicates the direction of rotation for analytic continuation.
///
/// Originally ZUNK1
fn ZUNK1<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    rotation: RotationDirection,
    n: usize,
) -> BesselResult<T> {
    let mut found_one_good_entry = false;
    let mut nz = 0;
    //-----------------------------------------------------------------------
    //     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION GREATER THAN
    //     THE UNDERFLOW LIMIT
    //-----------------------------------------------------------------------
    let zr = if z.re < T::ZERO { -z } else { z };
    let mut phi = [T::C_ZERO; 2];
    let mut zeta1 = [T::C_ZERO; 2];
    let mut zeta2 = [T::C_ZERO; 2];
    let mut sum = [T::C_ZERO; 2];
    let mut cy = [T::C_ZERO; 2];
    let mut n_elements_set = 0;
    let mut y = T::c_zeros(n);
    let mut k_overflow_state = Overflow::NearUnder;

    let mut j = 1;
    for i in 0..n {
        n_elements_set = i + 1;
        // j flip-flops between 0 and 1 using j = 1-j
        j = 1 - j;
        let modified_order = order + T::from_usize(i);

        let sum_opt;
        (phi[j], zeta1[j], zeta2[j], sum_opt) = zunik(zr, modified_order, IKType::K, false);
        sum[j] = sum_opt.unwrap();
        let mut s1 = -scaling.scale_zetas(zr, modified_order, zeta1[j], zeta2[j]);
        let of = Overflow::find_overflow(s1.re, phi[j], T::ZERO);
        if !found_one_good_entry {
            k_overflow_state = of;
        }
        match of {
            Overflow::Over(_) => return Err(Overflow),
            Overflow::Under(_) => {
                if z.re < T::ZERO {
                    return Err(Overflow);
                }
                found_one_good_entry = false;
                y[i] = T::C_ZERO;
                nz += 1;
            }
            Overflow::None | Overflow::NearOver | Overflow::NearUnder => {
                //-----------------------------------------------------------------------
                //     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
                //     EXPONENT EXTREMES
                //-----------------------------------------------------------------------
                let mut s2 = phi[j] * sum[j];
                s1 = T::MACHINE_CONSTANTS.scaling_factors[k_overflow_state] * s1.exp();
                s2 *= s1;
                let will_underflow = will_underflow(
                    s2,
                    T::MACHINE_CONSTANTS.overflow_boundary[0],
                    T::MACHINE_CONSTANTS.abs_error_tolerance,
                );
                if k_overflow_state != Overflow::NearUnder || !will_underflow {
                    cy[found_one_good_entry as usize] = s2;
                    y[i] = s2 * T::MACHINE_CONSTANTS.reciprocal_scaling_factors[k_overflow_state];
                    if found_one_good_entry {
                        break;
                    }
                    found_one_good_entry = true;
                } else if will_underflow {
                    if z.re < T::ZERO {
                        return Err(Overflow);
                    }
                    y[i] = T::C_ZERO;
                    nz += 1;
                    if i > 0 && y[i - 1] != T::C_ZERO {
                        y[i - 1] = T::C_ZERO;
                        nz += 1
                    }
                }
            }
        };
    }

    let rz = calc_rz(zr);
    if n_elements_set < n {
        //-----------------------------------------------------------------------
        //     TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW. SET SEQUENCE TO ZERO
        //     ON UNDERFLOW.
        //-----------------------------------------------------------------------
        let modified_order = order + T::from_usize(n - 1);
        let (phi, zet1d, zet2d, _sumd) = zunik(
            zr,
            modified_order,
            IKType::K,
            rotation == RotationDirection::None,
        );
        let overflow_test = -scaling.scale_zetas(zr, modified_order, zet1d, zet2d);

        match Overflow::find_overflow(overflow_test.re.abs(), phi, T::ZERO) {
            Overflow::Over(_) => return Err(Overflow),
            Overflow::Under(_) => {
                return if z.re < T::ZERO {
                    Err(Overflow)
                } else {
                    Ok((vec![T::C_ZERO; n], n))
                };
            }
            _ => (),
        }
        //---------------------------------------------------------------------------
        //     FORWARD RECUR FOR REMAINDER OF THE SEQUENCE
        //----------------------------------------------------------------------------
        let [s1, s2] = cy;
        recurr(
            true,
            order,
            zr,
            &mut y,
            n_elements_set,
            s1,
            s2,
            k_overflow_state,
        );
    }
    if rotation == RotationDirection::None {
        return Ok((y, nz));
    }
    //-----------------------------------------------------------------------
    //     ANALYTIC CONTINUATION FOR RE(Z) < 0.0
    //-----------------------------------------------------------------------
    nz = 0;
    let rotation_angle = -T::PI() * T::from_f64(rotation.signum());
    //-----------------------------------------------------------------------
    //     CSPN AND CSGN ARE COEFF OF K AND I FUNCTIONS RESP.
    //-----------------------------------------------------------------------

    let integer_order = order.to_i64().unwrap();
    let order_frac = order.fract();
    let modified_int_order = integer_order + (n as i64) - 1;
    let mut cspn = Complex::<T>::cis(order_frac * rotation_angle);
    if (modified_int_order % 2) != 0 {
        cspn = -cspn;
    }
    let mut dummy_n_good = 0;
    let mut found_one_good_entry = false;
    let mut i_overflow_state = Overflow::None;
    let mut remaining_n = n;
    for (i, yi) in y.iter_mut().enumerate().rev() {
        remaining_n = i;
        let modified_order = order + T::from_usize(i);
        //-----------------------------------------------------------------------
        //     LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K
        //     FUNCTION ABOVE
        //-----------------------------------------------------------------------
        // TODO no logic needed! as zunik ccahes the values. Should the other similar functions, too?
        let (phid, zet1d, zet2d, sumd) = zunik(zr, modified_order, IKType::I, false); //, &mut INITD);
        let sumd = sumd.unwrap();
        let mut s1 = scaling.scale_zetas(zr, modified_order, zet1d, zet2d);
        //-----------------------------------------------------------------------
        //     TEST FOR UNDERFLOW AND OVERFLOW
        //-----------------------------------------------------------------------
        let of = Overflow::find_overflow(s1.re, phid, T::ZERO);
        if !found_one_good_entry && !matches!(of, Overflow::Under(_)) {
            i_overflow_state = of;
        }
        let mut s2 = match of {
            Overflow::Over(_) => {
                return Err(Overflow);
            }
            Overflow::Under(_) => T::C_ZERO,
            Overflow::NearOver | Overflow::NearUnder | Overflow::None => {
                let st = phid * sumd;
                let mut s2 = T::I * st * rotation_angle;
                s1 = s1.exp() * T::MACHINE_CONSTANTS.scaling_factors[i_overflow_state];
                s2 *= s1;
                if i_overflow_state == Overflow::NearUnder
                    && will_underflow(
                        s2,
                        T::MACHINE_CONSTANTS.overflow_boundary[0],
                        T::MACHINE_CONSTANTS.abs_error_tolerance,
                    )
                {
                    s2 = T::C_ZERO;
                }
                s2
            }
        };
        cy[found_one_good_entry as usize] = s2;
        let c2 = s2;
        s2 *= T::MACHINE_CONSTANTS.reciprocal_scaling_factors[i_overflow_state];
        //-----------------------------------------------------------------------
        //     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N
        //-----------------------------------------------------------------------
        s1 = *yi;
        if scaling == Scaling::Scaled {
            nz += underflow_add_i_k(zr, &mut s1, &mut s2, &mut dummy_n_good);
        }
        *yi = s1 * cspn + s2;
        cspn = -cspn;
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
        let mut reciprocal_scale_factor =
            T::MACHINE_CONSTANTS.reciprocal_scaling_factors[i_overflow_state];
        let mut absolute_approximation_limit =
            T::MACHINE_CONSTANTS.overflow_boundary[i_overflow_state];
        for (i, yi) in y.iter_mut().enumerate().take(remaining_n).rev() {
            let modified_order = order + T::from_usize(i + 1);
            (s1, s2) = (s2, s1 + modified_order * (rz * s2));
            let mut unscaled_s2 = s2 * reciprocal_scale_factor;
            let ck = unscaled_s2;

            let mut c1 = *yi;
            if scaling == Scaling::Scaled {
                nz += underflow_add_i_k(zr, &mut c1, &mut unscaled_s2, &mut dummy_n_good);
            }
            *yi = c1 * cspn + unscaled_s2;
            cspn = -cspn;
            if i_overflow_state == Overflow::NearOver {
                continue;
            }
            if max_abs_component(unscaled_s2) <= absolute_approximation_limit {
                continue;
            }
            i_overflow_state.increment();
            absolute_approximation_limit = T::MACHINE_CONSTANTS.overflow_boundary[i_overflow_state];
            s1 *= reciprocal_scale_factor;
            s2 = ck; // ck is previously calculated s2 * reciprocal_scale_factor
            s1 *= T::MACHINE_CONSTANTS.scaling_factors[i_overflow_state];
            s2 *= T::MACHINE_CONSTANTS.scaling_factors[i_overflow_state];
            reciprocal_scale_factor =
                T::MACHINE_CONSTANTS.reciprocal_scaling_factors[i_overflow_state];
        }
    }
    Ok((y, nz))
}

/// zunk2 computes K(fnu,z) and its analytic continuation from the
/// right half plane to the left half plane by means of the
/// uniform asymptotic expansions for H(kind,fnu,zn) and J(fnu,zn)
/// where zn is in the right half plane, kind=(3-mr)/2, mr=+1 or
/// -1. here zn=zr*i or -zr*i where zr=z if z is in the right
/// half plane or zr=-z if z is in the left half plane.
///
/// `rotation` indicates the direction of rotation for analytic continuation.
///
/// Originally ZUNK2
fn ZUNK2<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    rotation: RotationDirection,
    n: usize,
) -> BesselResult<T> {
    let CR1: Complex<T> = Complex::<T>::new(T::one(), T::from_f64(1.73205080756887729));
    let CR2: Complex<T> = Complex::<T>::new(-T::half(), -T::from_f64(8.66025403784438647e-01));

    let mut found_one_good_entry = false;
    let mut nz = 0;
    let mut y = T::c_zeros(n);
    let zr = if z.re < T::ZERO { -z } else { z };
    let mut zn = -T::I * zr;
    let mut zb = zr;
    let integer_order = order.to_usize().unwrap();
    let order_fract = order.fract();
    let ANG = -T::FRAC_PI_2() * order_fract;
    let c2 = -T::I * Complex::<T>::from_polar(T::FRAC_PI_2(), ANG);
    let mut cs = CR1 * c2 * T::from_cpx64(CIP[integer_order % 4].conj());
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
    let mut asum = [None; 2];
    let mut bsum = [None; 2];
    let mut cy = [T::C_ZERO; 2];
    let mut j = 1;
    let mut overflow_state_k = Overflow::None;
    let mut n_elements_set = 0;

    for i in 0..n {
        n_elements_set = i + 1;
        // j flip-flops between 0 and 1 using  = 1-j
        j = 1 - j;
        let modified_order = order + T::from_usize(i);
        (phi[j], arg[j], zeta1[j], zeta2[j], asum[j], bsum[j]) = zunhj(zn, modified_order, false);
        let s1 = -scaling.scale_zetas(zb, modified_order, zeta1[j], zeta2[j]);
        let of = Overflow::find_overflow(
            s1.re,
            phi[j],
            T::from_f64(-0.25) * arg[j].abs().ln() - T::from_f64(AIC),
        );

        let mut handle_underflow = |of_already: &mut bool, cs_: &mut Complex<T>| {
            //-----------------------------------------------------------------------
            //     FOR ZR < 0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
            //-----------------------------------------------------------------------
            if z.re < T::ZERO {
                return Err(Overflow);
            }
            *of_already = false;
            y[i] = T::C_ZERO;
            nz += 1;
            *cs_ *= -T::I;
            if i != 0 && y[i - 1] != T::C_ZERO {
                y[i - 1] = T::C_ZERO;
                nz += 1;
            }
            Ok(())
        };

        if !found_one_good_entry {
            overflow_state_k = of;
        }

        match of {
            Overflow::Over(_) => return Err(Overflow),

            Overflow::Under(_) => handle_underflow(&mut found_one_good_entry, &mut cs)?,
            Overflow::NearOver | Overflow::NearUnder | Overflow::None => {
                //-----------------------------------------------------------------------;
                //     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR;
                //     EXPONENT EXTREMES;
                //-----------------------------------------------------------------------;
                let c2 = CR2 * arg[j];

                let (airy, d_airy) = airy_pair(c2);
                let pt = ((d_airy * bsum[j].unwrap()) * CR2 + (airy * asum[j].unwrap())) * phi[j];
                let mut s2 = pt * cs;
                let s1 = s1.exp() * T::MACHINE_CONSTANTS.scaling_factors[overflow_state_k];
                s2 *= s1;
                if overflow_state_k == Overflow::NearUnder
                    && will_underflow(
                        s2,
                        T::MACHINE_CONSTANTS.overflow_boundary[0],
                        T::MACHINE_CONSTANTS.abs_error_tolerance,
                    )
                {
                    handle_underflow(&mut found_one_good_entry, &mut cs)?
                }
                if zr.im <= T::ZERO {
                    s2 = s2.conj();
                }
                cy[found_one_good_entry as usize] = s2;
                y[i] = s2 * T::MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state_k];
                cs = -T::I * cs;
                if found_one_good_entry {
                    break;
                }
                found_one_good_entry = true;
            }
        };
    }

    let rz = calc_rz(zr);
    let mut phid = T::C_ZERO;
    let mut argd = T::C_ZERO;
    let mut zeta1d = T::C_ZERO;
    let mut zeta2d = T::C_ZERO;
    let mut asumd = None;
    let mut bsumd = None;
    let do_overflow_check = n_elements_set < n;
    if do_overflow_check {
        //-----------------------------------------------------------------------;
        //     TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW. SET SEQUENCE TO ZERO;
        //     ON UNDERFLOW.;
        //-----------------------------------------------------------------------;
        let modified_order = order + T::from_usize(n - 1);
        (phid, argd, zeta1d, zeta2d, asumd, bsumd) =
            zunhj(zn, modified_order, rotation == RotationDirection::None);
        let s1 = -scaling.scale_zetas(zb, modified_order, zeta1d, zeta2d);
        match Overflow::find_overflow(s1.re, phid, T::ZERO) {
            Overflow::Over(_) => return Err(Overflow),

            Overflow::Under(_) => {
                if z.re < T::ZERO {
                    return Err(Overflow);
                }
                return Ok((T::c_zeros(n), nz));
            }
            Overflow::NearOver | Overflow::None | Overflow::NearUnder => (),
        }
        let [s1, s2] = cy;
        recurr(
            true,
            order,
            zr,
            &mut y,
            n_elements_set,
            s1,
            s2,
            overflow_state_k,
        );
    }
    if rotation == RotationDirection::None {
        return Ok((y, nz));
    }
    //-----------------------------------------------------------------------
    //     ANALYTIC CONTINUATION FOR RE(Z) < 0.0
    //-----------------------------------------------------------------------
    nz = 0;
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
    let cos_sin = Complex::<T>::cis(ANG);
    // let mut cs = Complex::<T>::I * Complex::<T>::from_polar(CSGNI, ANG);
    let mut cs = csgn * Complex::<T>::new(cos_sin.im, cos_sin.re);
    cs *= T::from_cpx64(CIP[modified_integer_order % 4]);
    let mut IUF = 0;

    found_one_good_entry = false;
    let mut overflow_state_i = Overflow::None;
    let mut remaining_n = n;
    for (i, yi) in y.iter_mut().enumerate().rev() {
        remaining_n = i;
        let modified_order = order + T::from_usize(i);
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
            (phid, argd, zeta1d, zeta2d, asumd, bsumd) = zunhj(zn, modified_order, false);
        } else {
            // Case were overflow check has already set the ___d variables ?
        }
        let mut s1 = scaling.scale_zetas(zb, modified_order, zeta1d, zeta2d);

        let of = Overflow::find_overflow(
            s1.re,
            phid,
            T::from_f64(-0.25) * argd.abs().ln() - T::from_f64(AIC),
        );
        if !found_one_good_entry {
            overflow_state_i = if matches!(of, Overflow::Under(_)) {
                Overflow::None
            } else {
                of
            };
        }
        let mut s2 = match of {
            Overflow::Over(_) => return Err(Overflow),
            Overflow::Under(_) => T::C_ZERO,
            Overflow::NearOver | Overflow::None | Overflow::NearUnder => {
                let (airy, d_airy) = airy_pair(argd);
                let pt = ((d_airy * bsumd.unwrap()) + (airy * asumd.unwrap())) * phid;
                let mut s2 = pt * cs;
                s1 = s1.exp() * T::MACHINE_CONSTANTS.scaling_factors[overflow_state_i];
                s2 *= s1;
                if overflow_state_i == Overflow::NearUnder
                    && will_underflow(
                        s2,
                        T::MACHINE_CONSTANTS.overflow_boundary[0],
                        T::MACHINE_CONSTANTS.abs_error_tolerance,
                    )
                {
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
        s2 *= T::MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state_i];
        //-----------------------------------------------------------------------;
        //     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N;
        //-----------------------------------------------------------------------;
        s1 = *yi;
        if scaling == Scaling::Scaled {
            nz += underflow_add_i_k(zr, &mut s1, &mut s2, &mut IUF);
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

        let mut recip_scale_factor =
            T::MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state_i];
        let mut ascle = T::MACHINE_CONSTANTS.overflow_boundary[overflow_state_i];
        // TODO recurr with assignment fn
        for (i, yi) in y.iter_mut().enumerate().take(remaining_n).rev() {
            let modified_order = order + T::from_usize(i + 1);
            (s1, s2) = (s2, s1 + modified_order * (rz * s2));
            let mut c2 = s2 * recip_scale_factor;
            let old_c2 = c2;
            let mut c1 = *yi;
            if scaling == Scaling::Scaled {
                nz += underflow_add_i_k(zr, &mut c1, &mut c2, &mut IUF);
            }
            *yi = c1 * cspn + c2;
            cspn = -cspn;
            if overflow_state_i != Overflow::NearOver && max_abs_component(c2) > ascle {
                overflow_state_i.increment();
                ascle = T::MACHINE_CONSTANTS.overflow_boundary[overflow_state_i];
                s1 *= recip_scale_factor;
                s2 = old_c2;
                s1 *= T::MACHINE_CONSTANTS.scaling_factors[overflow_state_i];
                s2 *= T::MACHINE_CONSTANTS.scaling_factors[overflow_state_i];
                recip_scale_factor =
                    T::MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state_i];
            }
        }
    }
    Ok((y, nz))
}

/// zbuni computes the i bessel function for large cabs(z) >
/// asymptotic_order_limit and fnu+n-1 < asymptotic_order_limit.
/// The order is increased from
/// fnu+n-1 greater than asymptotic_order_limit by adding nui and computing
/// according to the uniform asymptotic expansion for J(fnu,z)
///  if z is imaginary_dominant and the expansion for I(fnu,z)
/// if z is _not_ imaginary_dominant.
///
/// Originally ZBUNI
fn ZBUNI<T: BesselFloat>(
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
            ZUNI2(z, GNU, KODE, 2, &mut cy)?
        } else {
            //-----------------------------------------------------------------------
            //     ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN
            //     -PI/3 <= ARG(Z) <= PI/3
            //-----------------------------------------------------------------------
            ZUNI1(z, GNU, KODE, 2, &mut cy)?
        };
        if NW != 0 {
            return Ok((0, n));
        }
        //----------------------------------------------------------------------
        //     SCALE BACKWARD RECURRENCE, BRY(3) IS DEFINED BUT NEVER USED
        //----------------------------------------------------------------------
        let (mut overflow_state, mut ASCLE, mut CSCLR) =
            if cy[0].abs() <= T::MACHINE_CONSTANTS.overflow_boundary[0] {
                (
                    Overflow::NearUnder,
                    T::MACHINE_CONSTANTS.overflow_boundary[0],
                    T::one() / T::MACHINE_CONSTANTS.abs_error_tolerance,
                )
            } else if cy[0].abs() >= T::MACHINE_CONSTANTS.overflow_boundary[1] {
                (
                    Overflow::NearOver,
                    T::MACHINE_CONSTANTS.overflow_boundary[2],
                    T::MACHINE_CONSTANTS.abs_error_tolerance,
                )
            } else {
                (
                    Overflow::None,
                    T::MACHINE_CONSTANTS.overflow_boundary[1],
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
            ASCLE = T::MACHINE_CONSTANTS.overflow_boundary[overflow_state];
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
            ASCLE = T::MACHINE_CONSTANTS.overflow_boundary[overflow_state];
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
        ZUNI2(z, order, KODE, n, y)?
    } else {
        //-----------------------------------------------------------------------
        //     ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN
        //     -PI/3 <= ARG(Z) <= PI/3
        //-----------------------------------------------------------------------
        ZUNI1(z, order, KODE, n, y)?
    };

    Ok((NW, NLAST))
}

/// zuni1 computes I(fnu,z)  by means of the uniform asymptotic
/// expansion for I(fnu,z) in -pi/3 <= arg z <= pi/3.
///
/// asymptotic_order_limit is the smallest order permitted for the asymptotic
/// expansion. nlast=0 means all of the y values were set.
/// nlast != 0 is the number left to be computed by another
/// formula for orders fnu to fnu+nlast-1 because
/// fnu+nlast-1 < asymptotic_order_limit.
/// y(i)=czero for i = nlast+1,n
///
/// Originally ZUNI1
fn ZUNI1<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    n: usize,
    y: &mut [Complex<T>],
) -> Result<(usize, usize), BesselError<T>> {
    let mut nz = 0;
    let mut n_remaining = n;
    //-----------------------------------------------------------------------
    //     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
    //-----------------------------------------------------------------------
    let mut modified_order = order.max(T::one());
    let (_, zeta1, zeta2, _) = zunik(z, modified_order, IKType::I, true);
    let s1 = scaling.scale_zetas(z, modified_order, zeta1, zeta2);
    // phi is chosen here for refined tests to equal the original tests
    // which don't test refinement
    match Overflow::find_overflow(s1.re, T::C_ONE, T::ZERO) {
        Overflow::Over(_) => return Err(Overflow),
        Overflow::Under(_) => return Ok((n, 0)),
        _ => (),
    }
    let mut overflow_state = Overflow::None; // this value should never be used
    let mut cy = [T::C_ZERO; 2];
    let mut handle_underflow = |n_remaining: &mut usize,
                                y: &mut [Complex<T>]|
     -> Result<bool, BesselError<T>> {
        //-----------------------------------------------------------------------
        //     SET UNDERFLOW AND UPDATE PARAMETERS
        //-----------------------------------------------------------------------
        y[*n_remaining - 1] = T::C_ZERO;
        nz += 1;
        *n_remaining -= 1;
        if *n_remaining == 0 {
            return Ok(true);
        }
        let n_underflow =
            check_underflow_uniform_asymp_params(z, order, scaling, IKType::I, *n_remaining, y)?;
        *n_remaining -= n_underflow;
        nz += n_underflow;
        if *n_remaining == 0 {
            return Ok(true);
        }
        let modified_order = order + T::from_usize(*n_remaining - 1);
        if modified_order < T::MACHINE_CONSTANTS.asymptotic_order_limit {
            return Ok(true);
        }
        Ok(false)
    };

    'outer: loop {
        for i in 0..2.min(n_remaining) {
            modified_order = order + T::from_usize(n_remaining - (i + 1));
            let (phi, zeta1, zeta2, sum) = zunik(z, modified_order, IKType::I, false);
            let sum = sum.unwrap();
            let mut s1 = scaling.scale_zetas(z, modified_order, zeta1, zeta2);
            if scaling == Scaling::Scaled {
                s1 += Complex::<T>::new(T::ZERO, z.im);
            }

            let of = Overflow::find_overflow(s1.re, phi, T::ZERO);
            if i == 0 {
                overflow_state = of;
            }
            match of {
                Overflow::Over(_) => return Err(Overflow),
                Overflow::Under(_) => {
                    if handle_underflow(&mut n_remaining, y)? {
                        return Ok((nz, n_remaining));
                    }
                    continue 'outer;
                }
                _ => (),
            }
            //-----------------------------------------------------------------------
            //     SCALE S1 if CABS(S1) < ASCLE
            //-----------------------------------------------------------------------
            let mut s2 = phi * sum;
            s1 = T::MACHINE_CONSTANTS.scaling_factors[overflow_state] * s1.exp();
            s2 *= s1;
            if overflow_state == Overflow::NearUnder
                && will_underflow(
                    s2,
                    T::MACHINE_CONSTANTS.overflow_boundary[0],
                    T::MACHINE_CONSTANTS.abs_error_tolerance,
                )
            {
                if handle_underflow(&mut n_remaining, y)? {
                    return Ok((nz, n_remaining));
                }
                continue 'outer;
            }
            cy[i] = s2;
            y[n_remaining - i - 1] =
                s2 * T::MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
        }
        break 'outer;
    }
    if n_remaining > 2 {
        let [s1, s2] = cy;
        recurr(false, order, z, y, n_remaining - 2, s1, s2, overflow_state);
    }
    Ok((nz, 0))
}

/// zuni2 computes I(fnu,z) in the right half plane by means of
/// uniform asymptotic expansion for J(fnu,zn) where zn is z*i
/// or -z*i and zn is in the right half plane also.
///
/// asymptotic_order_limit is the smallest order permitted for the asymptotic
/// expansion. nlast=0 means all of the y values were set.
/// nlast != 0 is the number left to be computed by another
/// formula for orders fnu to fnu+nlast-1 because fnu+nlast-1 < asymptotic_order_limit.
/// y(i)=czero for i=nlast+1,n
fn ZUNI2<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    n: usize,
    y: &mut [Complex<T>],
) -> Result<(usize, usize), BesselError<T>> {
    let mut nz = 0;
    let mut n_remaining = n;

    //-----------------------------------------------------------------------
    //     ZN IS IN THE RIGHT HALF PLANE AFTER ROTATION BY CI OR -CI
    //-----------------------------------------------------------------------
    let mut zn = Complex::<T>::new(z.im, -z.re);
    let mut zb = z;
    let integer_order = order.to_usize().unwrap();

    let build_c2 = |effective_n: usize| {
        let index = (integer_order + effective_n - 1) % 4;
        let mut c2 = Complex::<T>::cis(T::FRAC_PI_2() * order.fract()) * T::from_cpx64(CIP[index]);
        if z.im <= T::zero() {
            c2 = c2.conj();
        }
        c2
    };
    let mut c2 = build_c2(n);
    let sign_of_i = if z.im <= T::zero() {
        zn.re = -zn.re;
        zb.im = -zb.im;
        T::one()
    } else {
        -T::one()
    };
    //-----------------------------------------------------------------------
    //     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
    //-----------------------------------------------------------------------
    let mut modified_order = order.max(T::one());
    let (_, _, zeta1, zeta2, _, _) = zunhj(zn, modified_order, true);

    let s1 = scaling.scale_zetas(zb, modified_order, zeta1, zeta2);

    // phi is chosen here for refined tests to equal the original tests
    // which don't test refinement
    match Overflow::find_overflow(s1.re, T::C_ONE, T::zero()) {
        Overflow::Over(_) => return Err(Overflow),
        Overflow::Under(_) => return Ok((n, 0)),
        _ => (),
    }

    debug_assert!(
        modified_order + T::from_usize(n - 1) > T::MACHINE_CONSTANTS.asymptotic_order_limit
    );

    let mut overflow_state = Overflow::NearUnder;
    let mut cy = [T::C_ZERO; 2];
    let mut handle_underflow = |n_remaining: &mut usize,
                                c2: &mut Complex<T>,
                                y: &mut [Complex<T>]|
     -> Result<bool, BesselError<T>> {
        //-----------------------------------------------------------------------
        //     SET UNDERFLOW AND UPDATE PARAMETERS
        //-----------------------------------------------------------------------
        y[*n_remaining - 1] = T::C_ZERO;
        nz += 1;
        *n_remaining -= 1;
        if *n_remaining == 0 {
            return Ok(true);
        }
        let n_underflow =
            check_underflow_uniform_asymp_params(z, order, scaling, IKType::I, *n_remaining, y)?;
        *n_remaining -= n_underflow;
        nz += n_underflow;
        if *n_remaining == 0 {
            return Ok(true);
        }
        let modified_order = order + T::from_usize(*n_remaining - 1);
        if modified_order < T::MACHINE_CONSTANTS.asymptotic_order_limit {
            return Ok(true);
        }
        *c2 = build_c2(*n_remaining);
        Ok(false)
    };
    'outer: loop {
        for i in 0..2.min(n_remaining) {
            modified_order = order + T::from_usize(n_remaining - (i + 1));
            let (phi, arg, zeta1, zeta2, asum, bsum) = zunhj(zn, modified_order, false);
            let asum = asum.unwrap();
            let bsum = bsum.unwrap();
            let mut s1 = scaling.scale_zetas(zb, modified_order, zeta1, zeta2);
            if scaling == Scaling::Scaled {
                s1 += T::I * z.im.abs();
            }

            //-----------------------------------------------------------------------
            //     TEST FOR UNDERFLOW AND OVERFLOW
            //-----------------------------------------------------------------------
            let of = Overflow::find_overflow(
                s1.re,
                phi,
                T::from_f64(-0.25) * arg.abs().ln() - T::from_f64(AIC),
            );
            if i == 0 {
                overflow_state = of;
            }
            match of {
                Overflow::Over(_) => return Err(Overflow),
                Overflow::Under(_) => {
                    if handle_underflow(&mut n_remaining, &mut c2, y)? {
                        return Ok((nz, n_remaining));
                    }
                    continue 'outer;
                }
                _ => (),
            }
            //-----------------------------------------------------------------------
            //     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
            //     EXPONENT EXTREMES
            //-----------------------------------------------------------------------
            let (a_airy, d_airy) = airy_pair(arg);

            let mut s2 = phi * (d_airy * bsum + a_airy * asum);
            let s1 = T::MACHINE_CONSTANTS.scaling_factors[overflow_state] * s1.exp();
            s2 *= s1;
            if overflow_state == Overflow::NearUnder
                && will_underflow(
                    s2,
                    T::MACHINE_CONSTANTS.overflow_boundary[0],
                    T::MACHINE_CONSTANTS.abs_error_tolerance,
                )
            {
                if handle_underflow(&mut n_remaining, &mut c2, y)? {
                    return Ok((nz, n_remaining));
                }
                continue 'outer;
            }
            if z.im <= T::ZERO {
                s2 = s2.conj();
            }
            s2 *= c2;
            cy[i] = s2;
            y[n_remaining - i - 1] =
                s2 * T::MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
            c2 *= sign_of_i * T::I;
        }
        break 'outer;
    }
    if n_remaining > 2 {
        let [s1, s2] = cy;
        recurr(false, order, z, y, n_remaining - 2, s1, s2, overflow_state);
    }
    Ok((nz, 0))
}

#[allow(clippy::too_many_arguments)]
fn recurr<T: BesselFloat>(
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
