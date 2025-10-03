#![allow(non_snake_case, clippy::excessive_precision)]
use super::{
    BesselError, BesselResult, IKType, Scaling, c_one, c_zero, c_zeros, gamma_ln, i_power_series,
    overflow_checks::{check_underflow_uniform_asymp_params, underflow_add_i_k, zunik},
    utils::{AIC, TWO_THIRDS, will_underflow},
};
use crate::amos::{
    BesselError::*,
    CIP, MACHINE_CONSTANTS, RotationDirection,
    asymptotic_i::asymptotic_i,
    complex_airy, max_abs_component,
    overflow_checks::{Overflow, zunhj},
    utils::imaginary_dominant,
};
use itertools::Either;
use num::{
    Integer, Zero,
    complex::{Complex64, ComplexFloat},
    pow::Pow,
};
use std::{
    cmp::min,
    f64::consts::{FRAC_PI_2, PI},
};

pub fn airy_power_series(z: Complex64, return_derivative: bool, coeffs: (f64, f64)) -> Complex64 {
    let float_is_derivative = if return_derivative { 1.0 } else { 0.0 };

    let abs_z = z.abs();
    let z_floor = if abs_z < MACHINE_CONSTANTS.underflow_limit {
        c_zero()
    } else {
        z
    };
    let (s1, s2) = if abs_z < MACHINE_CONSTANTS.abs_error_tolerance {
        (c_one(), c_one())
    } else {
        let abs_z_sq = abs_z * abs_z;
        let mut s1 = c_one();
        let mut s2 = c_one();

        if abs_z_sq >= MACHINE_CONSTANTS.abs_error_tolerance / abs_z {
            let mut term1 = c_one();
            let mut term2 = c_one();
            let mut a_term = 1.0;
            let z3 = z.pow(3.0);
            let AZ3 = abs_z * abs_z_sq;
            let (AK, BK, CK, DK) = (
                2.0 + float_is_derivative,
                3.0 - 2.0 * float_is_derivative,
                4.0 - float_is_derivative,
                3.0 + 2.0 * float_is_derivative,
            );
            let mut D1: f64 = AK * DK;
            let mut D2 = BK * CK;
            let mut AD = D1.min(D2);
            let mut AK = 24.0 + 9.0 * float_is_derivative;
            let mut BK = 30.0 - 9.0 * float_is_derivative;
            for _ in 0..25 {
                term1 = term1 * z3 / D1;
                s1 += term1;
                term2 = term2 * z3 / D2;
                s2 += term2;
                a_term = a_term * AZ3 / AD;
                D1 += AK;
                D2 += BK;
                AD = D1.min(D2);
                if a_term < MACHINE_CONSTANTS.abs_error_tolerance * AD {
                    break;
                }
                AK += 18.0;
                BK += 18.0;
            }
        }
        (s1, s2)
    };
    let (c1, c2) = coeffs;
    if return_derivative {
        (c1 / 2.0) * z_floor.pow(2.0) * s1 - s2 * c2
    } else {
        s1 * c1 - c2 * z_floor * s2
    }
}

pub fn airy_pair(z: Complex64) -> (Complex64, Complex64) {
    //note that ZAIRY calls in fortran code ignore IERR (using IDUM)
    let airy = match complex_airy(z, false, Scaling::Scaled) {
        Ok((y, _)) => y,
        Err(PartialLossOfSignificance { y, nz: _ }) => y[0],
        // If loss of significance, Fortran code would continue with un-initialised y,
        // which is usually ~=0. Also long as it is << d_airy, the logic below means
        // it will not matter what the precise value is
        Err(LossOfSignificance) => c_zero(),
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
        Err(LossOfSignificance) => c_zero(),
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
pub fn k_right_half_plane(z: Complex64, order: f64, scaling: Scaling, n: usize) -> BesselResult {
    const KMAX: usize = 30;
    const RTFRAC_PI_2: f64 = 1.25331413731550025;
    const SPI: f64 = 1.90985931710274403;
    const FPI: f64 = 1.89769999331517738;
    const CC: [f64; 8] = [
        5.77215664901532861e-01,
        -4.20026350340952355e-02,
        -4.21977345555443367e-02,
        7.21894324666309954e-03,
        -2.15241674114950973e-04,
        -2.01348547807882387e-05,
        1.13302723198169588e-06,
        6.11609510448141582e-09,
    ];

    let abs_z = z.abs();
    let mut nz = 0;
    let mut underflow_occurred = false;
    let mut overflow_state;
    let rz = 2.0 * z.conj() / abs_z.powi(2);
    let mut integer_order = (order + 0.5) as isize; // round to nearest int
    let simple_case = integer_order == 0 && n == 1;

    let signed_fractional_order = order - (integer_order as f64); // signed fractional part (-0.5 < DNU < 0.5 )
    let frac_order_sqr = if signed_fractional_order.abs() > MACHINE_CONSTANTS.abs_error_tolerance {
        signed_fractional_order * signed_fractional_order
    } else {
        0.0
    };

    let (mut s1, mut s2) = if (signed_fractional_order.abs() != 0.5) && (abs_z <= 2.0) {
        // series for (z.abs() <= 2.0) and not half integer order
        let mut fc = 1.0;
        let mut smu = rz.ln();
        let fmu = smu * signed_fractional_order;
        let csh = fmu.sinh();
        let cch = fmu.cosh();
        if signed_fractional_order != 0.0 {
            fc = signed_fractional_order * PI;
            fc /= fc.sin();
            smu = csh / signed_fractional_order;
        }
        //-----------------------------------------------------------------------;
        //     GAM(1-Z)*GAM(1+Z)=PI*Z/SIN(PI*Z), T1=1/GAM(1-DNU), T2=1/GAM(1+DNU);
        //-----------------------------------------------------------------------;
        let t2 = (-gamma_ln(1.0 + signed_fractional_order).unwrap()).exp();
        let t1 = 1.0 / (t2 * fc);

        let g1 = if signed_fractional_order.abs() <= 0.1 {
            //-----------------------------------------------------------------------;
            //     SERIES FOR F0 TO RESOLVE INDETERMINACY FOR SMALL ABS(DNU);
            //-----------------------------------------------------------------------;
            let mut ak = 1.0;
            let mut sum = CC[0];
            for cc in CC[1..].iter() {
                ak *= frac_order_sqr;
                let tm = cc * ak;
                sum += tm;
                if tm.abs() < MACHINE_CONSTANTS.abs_error_tolerance {
                    break;
                }
            }
            -sum
        } else {
            (t1 - t2) / (2.0 * signed_fractional_order)
        };
        let g2 = (t1 + t2) * 0.5;
        let f = fc * (g1 * cch + g2 * smu);
        let p = 0.5 * fmu.exp() / t2;
        let q = (0.5 / fmu.exp()) / t1;

        let ck = c_one();
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

        overflow_state = if (order + 1.0) * smu.re.abs() > MACHINE_CONSTANTS.approximation_limit {
            Overflow::NearOver
        } else {
            Overflow::None
        };
        s2 *= MACHINE_CONSTANTS.scaling_factors[overflow_state] * rz;
        s1 *= MACHINE_CONSTANTS.scaling_factors[overflow_state];
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
        let mut coeff = Complex64::new(RTFRAC_PI_2, 0.0) / z.sqrt();
        overflow_state = Overflow::None;
        if scaling == Scaling::Unscaled {
            if z.re > MACHINE_CONSTANTS.approximation_limit {
                underflow_occurred = true;
                overflow_state = Overflow::NearUnder;
            } else {
                coeff *= MACHINE_CONSTANTS.scaling_factors[overflow_state] * (-z).exp();
            }
        }
        let mut AK = (signed_fractional_order * PI).cos().abs();
        let mut FHS = (0.25 - frac_order_sqr).abs();

        if signed_fractional_order.abs() == 0.5 || AK == 0.0 || FHS == 0.0 {
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
            let f64_significant_digits =
                (f64::MANTISSA_DIGITS - 1) as f64 * (f64::RADIX as f64).log10();
            let determiner = (f64_significant_digits * std::f64::consts::LOG2_10).clamp(12.0, 60.0);
            let recurrence_threshold = TWO_THIRDS * determiner - 6.0;
            let arg_z = z.arg();

            let (FK, FHS) = if abs_z > recurrence_threshold {
                //-----------------------------------------------------------------------;
                //     FORWARD RECURRENCE LOOP WHEN CABS(Z) >= R2;
                //-----------------------------------------------------------------------;
                let convergence_test = AK / (PI * abs_z * MACHINE_CONSTANTS.abs_error_tolerance);
                let mut FK = 1.0;
                if convergence_test >= 1.0 {
                    let mut FKS = 2.0;
                    let mut CKR = abs_z + abs_z + 2.0;
                    let mut p1 = 0.0;
                    let mut p2 = 1.0;
                    let mut converged = false;
                    for _ in 0..KMAX {
                        let AK = FHS / FKS;
                        let CBR = CKR / (FK + 1.0);
                        let pt = p2;
                        p2 = CBR * p2 - AK * p1;
                        p1 = pt;
                        CKR += 2.0;
                        FKS += FK + FK + 2.0;
                        FHS += FK + FK;
                        FK += 1.0;
                        if convergence_test < p2.abs() * FK {
                            converged = true;
                            break;
                        }
                    }
                    if !converged {
                        return Err(DidNotConverge);
                    }
                    FK += SPI * arg_z * (recurrence_threshold / abs_z).sqrt();
                    FHS = (0.25 - frac_order_sqr).abs();
                }
                (FK, FHS)
            } else {
                //-----------------------------------------------------------------------;
                //     COMPUTE BACKWARD INDEX K FOR CABS(Z) < R2;
                //-----------------------------------------------------------------------;
                AK *= FPI / (MACHINE_CONSTANTS.abs_error_tolerance * abs_z.sqrt().sqrt());
                let AA = 3.0 * arg_z / (1.0 + abs_z);
                let BB = 14.7 * arg_z / (28.0 + abs_z);
                AK = (AK.ln() + abs_z * AA.cos() / (1.0 + 0.008 * abs_z)) / BB.cos();
                let FK = 0.12125 * AK * AK / abs_z + 1.5;
                (FK, FHS)
            };
            //-----------------------------------------------------------------------;
            //     BACKWARD RECURRENCE LOOP FOR MILLER ALGORITHM;
            //-----------------------------------------------------------------------;
            let K = FK as usize;
            let mut k_squared = FK.floor().pow(2);
            let mut p1 = Complex64::zero();
            let mut p2 = Complex64::new(MACHINE_CONSTANTS.abs_error_tolerance, 0.0);
            let mut cs = p2;
            for i in (0..K).rev() {
                let k_f64 = (i + 1) as f64;
                let cb = (z + k_f64) * 2.0 / (k_f64 + 1.0);
                (p1, p2) = (
                    p2,
                    (p2 * cb - p1) * (k_squared + k_f64) / (k_squared - k_f64 + FHS),
                );
                cs += p2;
                k_squared -= (2.0 * k_f64) - 1.0;
            }
            //-----------------------------------------------------------------------;
            //     COMPUTE (P2/CS)=(P2/CABS(CS))*(CONJG(CS)/CABS(CS)) FOR BETTER;
            //     SCALING;
            //-----------------------------------------------------------------------;
            let mut s1 = p2 / cs.abs();
            let mut s2 = Complex64::zero();
            cs = cs.conj() / cs.abs();
            s1 *= coeff * cs;
            if !simple_case {
                //-----------------------------------------------------------------------;
                //     COMPUTE P1/P2=(P1/CABS(P2)*CONJG(P2)/CABS(P2) FOR SCALING;
                //-----------------------------------------------------------------------;
                p1 /= p2.abs();
                p2 = p2.conj() / p2.abs();
                s2 = (((signed_fractional_order + 0.5 - (p1 * p2)) / z) + 1.0) * s1;
            }
            (s1, s2)
        }
    };

    // Now s1, s2 set up, we can go to recurrence

    //-----------------------------------------------------------------------
    //     FORWARD RECURSION ON THE THREE TERM RECURSION WITH RELATION WITH
    //     SCALING NEAR EXPONENT EXTREMES ON KFLAG=1 OR KFLAG=3
    //-----------------------------------------------------------------------
    let mut ck = (signed_fractional_order + 1.0) * rz;
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
                let mut cy = [c_zero(); 2];
                let half_exponent_limit = 0.5 * MACHINE_CONSTANTS.exponent_limit;

                let abs_limit = (-MACHINE_CONSTANTS.exponent_limit).exp();
                let ASCLE = MACHINE_CONSTANTS.smallness_threshold[0];
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
                    if -zd.re + abs_ln_s2 >= -MACHINE_CONSTANTS.exponent_limit {
                        let p1 = (-zd + s2.ln()).exp() / MACHINE_CONSTANTS.abs_error_tolerance;
                        if !will_underflow(p1, ASCLE, MACHINE_CONSTANTS.abs_error_tolerance) {
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
                        zd.re -= MACHINE_CONSTANTS.exponent_limit;
                        s1 *= abs_limit;
                        s2 *= abs_limit;
                    }
                }
                overflow_state = Overflow::NearUnder;

                s2 = cy[J];
                J = 1 - J;
                s1 = cy[J];
            }

            let mut P1R = MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
            let mut ASCLE = MACHINE_CONSTANTS.smallness_threshold[overflow_state];
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
                ASCLE = MACHINE_CONSTANTS.smallness_threshold[overflow_state];
                s1 *= P1R;
                s2 = p2;
                s1 *= MACHINE_CONSTANTS.scaling_factors[overflow_state];
                s2 *= MACHINE_CONSTANTS.scaling_factors[overflow_state];
                P1R = MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
            }
        }
        if n == 1 {
            s1 = s2;
        }
    }

    let mut y = c_zeros(n);
    let n_completed = if !underflow_occurred {
        // ********* basic setup
        y[0] = s1 * MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
        if n > 1 {
            y[1] = s2 * MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
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
            MACHINE_CONSTANTS.absolute_approximation_limit,
        );
        let n_non_zero = (n - nz) as isize;
        if n_non_zero <= 0 {
            return Ok((y, nz));
        }
        let mut working_index = nz;
        s1 = y[working_index];
        y[working_index] *= MACHINE_CONSTANTS.reciprocal_scaling_factors[0];
        if n_non_zero > 1 {
            // if n_non_zero == 1 {
            //     return Ok((y, nz));
            // }
            working_index += 1;
            s2 = y[working_index];
            y[working_index] *= MACHINE_CONSTANTS.reciprocal_scaling_factors[0];
        }
        if n_non_zero > 2 {
            ck = (order + (working_index as f64)) * rz;
            overflow_state = Overflow::NearUnder;
        }
        working_index + 1
    };
    // End Setup
    if n_completed >= n {
        return Ok((y, nz));
    }
    let mut P1R = MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
    let mut ASCLE = MACHINE_CONSTANTS.smallness_threshold[overflow_state];
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
        ASCLE = MACHINE_CONSTANTS.smallness_threshold[overflow_state];
        s1 *= P1R;
        s2 = *y_elem;
        s1 *= MACHINE_CONSTANTS.scaling_factors[overflow_state];
        s2 *= MACHINE_CONSTANTS.scaling_factors[overflow_state];
        P1R = MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
    }
    Ok((y, nz))
}

fn k_right_half_plane_helper(
    z: Complex64,
    frac_order_sqr: f64,
    signed_fractional_order: f64,
    mut f: Complex64,
    mut p: Complex64,
    mut q: Complex64,
    mut ck: Complex64,
) -> (Complex64, Complex64) {
    let mut a1 = 1.0;
    let cz_sqr_over_4 = 0.25 * z.powu(2);
    let abs_z = z.abs();
    let abs_z_sqr_over_4 = 0.25 * abs_z * abs_z;
    let mut ak = 1.0;
    let mut bk = 1.0 - frac_order_sqr;

    let mut s1 = f;
    let mut s2 = p;
    if abs_z >= MACHINE_CONSTANTS.abs_error_tolerance {
        while a1 > MACHINE_CONSTANTS.abs_error_tolerance {
            f = (f * ak + p + q) / bk;
            p /= ak - signed_fractional_order;
            q /= ak + signed_fractional_order;
            ck *= cz_sqr_over_4 / ak;
            s1 += ck * f;
            s2 += ck * (p - ak * f);
            a1 *= abs_z_sqr_over_4 / ak;
            bk += (2.0 * ak) + 1.0;
            ak += 1.0;
        }
    }
    (s1, s2)
}

/// Set k functions to zero on underflow, continue recurrence
/// on scaled functions until two members come on scale, then
/// return with min(nz+2,n) values scaled by 1/tol.
fn ZKSCL(
    zr: Complex64,
    order: f64,
    n: usize,
    y: &mut [Complex64],
    nz: &mut usize,
    rz: Complex64,
    ASCLE: f64,
) {
    *nz = 0;
    // let NN = min(2, n);
    let mut cy = [c_zero(); 2];
    let mut i_completed = 0;
    // repeats twice, unless n < 2
    for i in 0..min(2, n) {
        let s1 = y[i];
        cy[i] = s1;
        *nz += 1;
        y[i] = c_zero();
        if -zr.re + s1.abs().ln() < -MACHINE_CONSTANTS.exponent_limit {
            continue;
        }

        let cs = (s1.ln() - zr).exp() / MACHINE_CONSTANTS.abs_error_tolerance;
        if will_underflow(cs, ASCLE, MACHINE_CONSTANTS.abs_error_tolerance) {
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
    //     y[0] = c_zero();
    //     *nz = 2;
    // }
    // if n == 2 {
    //     return;
    // }
    // if *nz == 0 {
    //     return;
    // }
    let FN = order + 1.0;
    let mut ck = FN * rz;
    let mut s1 = cy[0];
    let mut s2 = cy[1];
    let half_elim = 0.5 * MACHINE_CONSTANTS.exponent_limit;
    let ELM = (-MACHINE_CONSTANTS.exponent_limit).exp();
    let CELMR = ELM;
    let mut zd = zr;
    //     FIND TWO CONSECUTIVE Y VALUES ON SCALE. SCALE RECURRENCE if
    //     S2 GETS LARGER THAN EXP(ELIM/2)
    let mut skip_to_40 = false;
    let mut I = 0;
    for i in 2..n {
        I = i;
        let mut cs = s2;
        s2 = cs * ck + s1;
        s1 = cs;
        ck += rz;
        let ALAS = s2.abs().ln();
        *nz += 1;
        y[i] = Complex64::zero();
        if -zd.re + s2.abs().ln() >= -MACHINE_CONSTANTS.exponent_limit {
            cs = s2.ln() - zd;
            cs = cs.exp() / MACHINE_CONSTANTS.abs_error_tolerance;
            if !will_underflow(cs, ASCLE, MACHINE_CONSTANTS.abs_error_tolerance) {
                y[i] = cs;
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
        zd -= MACHINE_CONSTANTS.exponent_limit;
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
        *element = c_zero();
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
fn ratios_i(z: Complex64, order: f64, n: usize) -> Vec<Complex64> {
    let abs_z = z.abs();
    let integer_order = order as usize;
    let modified_int_order = integer_order + n - 1;
    let int_abs_z = abs_z as isize;
    let FNUP = (int_abs_z + 1).max(modified_int_order as isize) as f64;
    let ID_ = modified_int_order as isize - int_abs_z - 1;
    let ID = if ID_ > 0 { 0 } else { ID_ };

    let rz = 2.0 * z.conj() / abs_z.powi(2);
    let mut K = 1;
    let mut abs_p2;
    {
        let mut t1 = rz * FNUP;
        let mut p2 = -t1;
        let mut p1 = c_one();
        t1 += rz;

        abs_p2 = p2.abs();
        let mut abs_p1 = p1.abs();
        //-----------------------------------------------------------------------
        //     THE OVERFLOW TEST ON K(FNU+I-1,Z) BEFORE THE CALL TO CBKNU
        //     GUARANTEES THAT P2 IS ON SCALE. SCALE TEST1 AND ALL SUBSEQUENT
        //     P2 VALUES BY AP1 TO ENSURE THAT AN OVERFLOW DOES NOT OCCUR
        //     PREMATURELY.
        //-----------------------------------------------------------------------
        let ARG = (abs_p2 + abs_p2) / (abs_p1 * MACHINE_CONSTANTS.abs_error_tolerance);
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
                let ak = t1.abs() / 2.0;
                let flam = ak + (ak.powi(2) - 1.0).sqrt();
                let rho = abs_p2 / abs_p1.min(flam);
                TEST = TEST1 * (rho / (rho.powi(2) - 1.0)).sqrt();
            }
            first_pass = false;
        }
    }

    let mut p1 = Complex64::new(1.0 / abs_p2, 0.0);
    let mut p2 = c_zero();

    {
        let kk: usize = (K as isize + 1 - ID).try_into().unwrap();
        let mut t1 = Complex64::new(kk as f64, 0.0);
        let modified_order = order + ((n - 1) as f64);
        for _ in 0..kk {
            (p1, p2) = (p1 * (rz * (modified_order + t1.re)) + p2, p1);
            t1.re -= 1.0;
        }
        if p1.re == 0.0 && p1.im == 0.0 {
            p1 = Complex64::new(
                MACHINE_CONSTANTS.abs_error_tolerance,
                MACHINE_CONSTANTS.abs_error_tolerance,
            );
        }
    }
    let mut cy = c_zeros(n);
    cy[n - 1] = p2 / p1;
    if n > 1 {
        let mut t1 = Complex64::new((n - 1) as f64, 0.0);
        let cdfnu = order * rz;
        for k in (1..n).rev() {
            let mut pt = cdfnu + t1 * rz + cy[k];
            let mut abs_pt = pt.abs();
            if abs_pt == 0.0 {
                pt = Complex64::new(
                    MACHINE_CONSTANTS.abs_error_tolerance,
                    MACHINE_CONSTANTS.abs_error_tolerance,
                );
                abs_pt = pt.abs();
            }
            cy[k - 1] = pt.conj() / abs_pt.powi(2);
            t1 -= 1.0;
        }
    }
    cy
}

pub fn ZBUNK(
    z: Complex64,
    order: f64,
    scaling: Scaling,
    rotation: RotationDirection,
    N: usize,
) -> BesselResult {
    // ***BEGIN PROLOGUE  ZBUNK
    // ***REFER TO  ZBESK,ZBESH
    //
    //     ZBUNK COMPUTES THE K BESSEL FUNCTION FOR FNU > FNUL.
    //     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR K(FNU,Z)
    //     IN ZUNK1 AND THE EXPANSION FOR H(2,FNU,Z) IN ZUNK2
    //
    if imaginary_dominant(z) {
        //-----------------------------------------------------------------------
        //     ASYMPTOTIC EXPANSION FOR H(2,FNU,Z*EXP(M*FRAC_PI_2)) FOR LARGE FNU
        //     APPLIED IN PI/3 < ABS(ARG(Z)) <= PI/2 WHERE M=+I OR -I
        //     AND FRAC_PI_2=PI/2
        //-----------------------------------------------------------------------
        ZUNK2(z, order, scaling, rotation, N)
    } else {
        //-----------------------------------------------------------------------
        //     ASYMPTOTIC EXPANSION FOR K(FNU,Z) FOR LARGE FNU APPLIED IN
        //     -PI/3 <= ARG(Z) <= PI/3
        //-----------------------------------------------------------------------
        ZUNK1(z, order, scaling, rotation, N)
    }
}

fn i_miller(z: Complex64, order: f64, KODE: Scaling, N: usize) -> BesselResult {
    // ***BEGIN PROLOGUE  ZMLRI
    // ***REFER TO  ZBESI,ZBESK
    //
    //     ZMLRI COMPUTES THE I BESSEL FUNCTION FOR RE(Z) >= 0.0 BY THE
    //     MILLER ALGORITHM NORMALIZED BY A NEUMANN SERIES.
    //

    let SCLE: f64 = 2.0 * f64::MIN_POSITIVE / MACHINE_CONSTANTS.abs_error_tolerance;
    let NZ = 0;
    let AZ = z.abs();
    let IAZ = AZ as usize;
    let IFNU = order as usize;
    let INU = IFNU + N - 1;
    let AT = (IAZ as f64) + 1.0;
    let RAZ = 1.0 / AZ;
    let mut ck = z.conj() * RAZ * RAZ * AT;
    let rz = z.conj() * 2.0 * RAZ * RAZ;
    let mut p1 = c_zero();
    let mut p2 = c_one();
    let mut ACK = (AT + 1.0) * RAZ;
    let RHO = ACK + (ACK * ACK - 1.0).sqrt();
    let RHO2 = RHO * RHO;
    let mut TST = (RHO2 + RHO2) / ((RHO2 - 1.0) * (RHO - 1.0));
    TST /= MACHINE_CONSTANTS.abs_error_tolerance;
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
        AK += 1.0;
    }
    if !converged {
        return Err(DidNotConverge);
    }
    I += 1;
    let mut K = 0;
    if INU >= IAZ {
        //-----------------------------------------------------------------------
        //     COMPUTE RELATIVE TRUNCATION ERROR FOR RATIOS
        //-----------------------------------------------------------------------
        p1 = c_zero();
        p2 = c_one();
        let AT = (INU as f64) + 1.0;
        ck = z.conj() * RAZ * RAZ * AT;
        ACK = AT * RAZ;
        TST = (ACK / MACHINE_CONSTANTS.abs_error_tolerance).sqrt();
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
            let FLAM = ACK + (ACK * ACK - 1.0).sqrt();
            let FKAP = AP / p1.abs();
            let RHO = FLAM.min(FKAP);
            TST *= (RHO / (RHO * RHO - 1.0)).sqrt();
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
    let KK = (I + IAZ).max(K + INU);
    let mut FKK = KK as f64;
    let mut p1 = c_zero();
    //-----------------------------------------------------------------------
    //     SCALE P2 AND SUM BY SCLE
    //-----------------------------------------------------------------------
    let mut p2 = Complex64::new(SCLE, 0.0);
    let FNF = order - (IFNU as f64);
    let TFNF = FNF + FNF;
    let mut BK = (gamma_ln(FKK + TFNF + 1.0).unwrap()
        - gamma_ln(FKK + 1.0).unwrap()
        - gamma_ln(TFNF + 1.0).unwrap())
    .exp();
    let mut sumr = c_zero();
    for _ in 0..(KK - INU) {
        let pt = p2;
        p2 = p1 + (FKK + FNF) * (rz * p2);
        p1 = pt;
        AK = 1.0 - TFNF / (FKK + TFNF);
        ACK = BK * AK;
        sumr += (ACK + BK) * p1;
        BK = ACK;
        FKK -= 1.0;
    }
    let mut y = c_zeros(N);
    y[N - 1] = p2;
    if N != 1 {
        for i in 1..N {
            let pt = p2;
            p2 = p1 + (FKK + FNF) * (rz * pt);
            p1 = pt;
            AK = 1.0 - TFNF / (FKK + TFNF);
            ACK = BK * AK;
            sumr += (ACK + BK) * p1;
            BK = ACK;
            FKK -= 1.0;
            y[N - (i + 1)] = p2;
        }
    }
    if IFNU > 0 {
        for _i in 0..IFNU {
            let pt = p2;
            p2 = p1 + (FKK + FNF) * (rz * pt);
            p1 = pt;
            AK = 1.0 - TFNF / (FKK + TFNF);
            ACK = BK * AK;
            sumr += (ACK + BK) * p1;
            BK = ACK;
            FKK -= 1.0;
        }
    }

    let mut pt = z;
    if KODE == Scaling::Scaled {
        pt.re = 0.0;
    }
    p1 = -FNF * rz.ln() + pt;
    let AP = gamma_ln(1.0 + FNF).unwrap();
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
    Ok((y, NZ))
}

fn i_wronksian(
    zr: Complex64,
    order: f64,
    KODE: Scaling,
    N: usize,
    y: &mut [Complex64],
) -> BesselResult<usize> {
    // ***BEGIN PROLOGUE  ZWRSK
    // ***REFER TO  ZBESI,ZBESK
    //
    //     ZWRSK COMPUTES THE I BESSEL FUNCTION FOR RE(Z) >= 0.0 BY
    //     NORMALIZING THE I FUNCTION RATIOS FROM ZRATI BY THE WRONSKIAN
    //

    //-----------------------------------------------------------------------
    //     I(FNU+I-1,Z) BY BACKWARD RECURRENCE FOR RATIOS
    //     Y(I)=I(FNU+I,Z)/I(FNU+I-1,Z) FROM CRATI NORMALIZED BY THE
    //     WRONSKIAN WITH K(FNU,Z) AND K(FNU+1,Z) FROM CBKNU.
    //-----------------------------------------------------------------------
    let NZ = 0;
    let (cw, _) = k_right_half_plane(zr, order, KODE, 2)?;
    let y_ratios = ratios_i(zr, order, N);
    //-----------------------------------------------------------------------
    //     RECUR FORWARD ON I(FNU+1,Z) = R(FNU,Z)*I(FNU,Z),
    //     R(FNU+J-1,Z)=Y(J),  J=1,...,N
    //-----------------------------------------------------------------------
    let mut cinu = c_one();
    if KODE == Scaling::Scaled {
        cinu = Complex64::cis(zr.im);
    }
    //-----------------------------------------------------------------------
    //     ON LOW EXPONENT MACHINES THE K FUNCTIONS CAN BE CLOSE TO BOTH
    //     THE UNDER AND OVERFLOW LIMITS AND THE NORMALIZATION MUST BE
    //     SCALED TO PREVENT OVER OR UNDERFLOW. CUOIK HAS DETERMINED THAT
    //     THE RESULT IS ON SCALE.
    //-----------------------------------------------------------------------
    let acw = cw[1].abs();
    let CSCLR = if acw <= MACHINE_CONSTANTS.absolute_approximation_limit {
        1.0 / MACHINE_CONSTANTS.abs_error_tolerance
    } else if acw >= 1.0 / MACHINE_CONSTANTS.absolute_approximation_limit {
        MACHINE_CONSTANTS.abs_error_tolerance
    } else {
        1.0
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
    for i in 1..N {
        cinu *= y_ratios[i - 1];
        y[i] = cinu * CSCLR;
    }
    Ok(NZ)
}

pub fn analytic_continuation(
    z: Complex64,
    order: f64,
    scaling: Scaling,
    rotation: RotationDirection,
    n: usize,
) -> BesselResult {
    // ***BEGIN PROLOGUE  ZACON
    // ***REFER TO  ZBESK,ZBESH
    //
    //     ZACON APPLIES THE ANALYTIC CONTINUATION FORMULA
    //
    //         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
    //                 MP=PI*MR*CMPLX(0.0,1.0)
    //
    //     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT
    //     HALF Z PLANE
    //

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
    let SGN = -PI * rotation.signum();
    let mut csgn = Complex64::new(0.0, SGN);
    if scaling == Scaling::Scaled {
        csgn *= Complex64::cis(-zn.im);
    }
    //-----------------------------------------------------------------------
    //     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
    //     WHEN FNU IS LARGE
    //-----------------------------------------------------------------------
    let mut cspn = Complex64::cis(order.fract() * SGN);
    if order as i64 % 2 != 0 {
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
    let mut scaled_c2 = c_zero() * f64::NAN;
    if scaling == Scaling::Scaled {
        nz += underflow_add_i_k(zn, &mut c1, &mut c2, &mut n_good);
        scaled_c2 = c1;
    }
    y[1] = cspn * c1 + csgn * c2;
    if n == 2 {
        return Ok((y, nz));
    }

    cspn = -cspn;
    let rz = 2.0 * (zn.conj()) / zn.abs().powi(2);
    let FN = order + 1.0;
    let mut ck = FN * rz;
    //-----------------------------------------------------------------------
    //     SCALE NEAR EXPONENT EXTREMES DURING RECURRENCE ON K FUNCTIONS
    //-----------------------------------------------------------------------
    let abs_s2 = s2.abs();
    let mut overflow_state = if abs_s2 <= MACHINE_CONSTANTS.smallness_threshold[0] {
        Overflow::NearUnder
    } else if abs_s2 > MACHINE_CONSTANTS.smallness_threshold[1] {
        Overflow::NearOver
    } else {
        Overflow::None
    };
    let mut boundary = MACHINE_CONSTANTS.smallness_threshold[overflow_state];
    s1 *= MACHINE_CONSTANTS.scaling_factors[overflow_state];
    s2 *= MACHINE_CONSTANTS.scaling_factors[overflow_state];
    let mut recip_scaling_factor = MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
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
                s1 = saved_c2 * MACHINE_CONSTANTS.scaling_factors[overflow_state];
                s2 = scaled_c2 * MACHINE_CONSTANTS.scaling_factors[overflow_state];
                st = scaled_c2;
            }
        }
        *yi = cspn * c1 + csgn * c2;
        ck += rz;
        cspn = -cspn;
        if overflow_state != Overflow::NearOver && max_abs_component(c1) < boundary {
            overflow_state.increment();
            boundary = MACHINE_CONSTANTS.smallness_threshold[overflow_state];
            s1 *= recip_scaling_factor;
            s2 = st;
            s1 *= MACHINE_CONSTANTS.scaling_factors[overflow_state];
            s2 *= MACHINE_CONSTANTS.scaling_factors[overflow_state];
            recip_scaling_factor = MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
        }
    }
    Ok((y, nz))
}

/// i_right_half_plane computes the i function in the right half z plane
/// Originally ZBINU
pub fn i_right_half_plane(z: Complex64, order: f64, KODE: Scaling, N: usize) -> BesselResult {
    let mut NZ = 0;
    let AZ = z.abs();
    let mut NN: usize = N;
    let mut DFNU = order + ((N - 1) as f64);
    let mut cy = c_zeros(N);
    if AZ <= 2.0 || AZ * AZ * 0.25 <= DFNU + 1.0 {
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

        DFNU = order + ((NN as f64) - 1.0);
    }

    if (AZ >= MACHINE_CONSTANTS.asymptotic_z_limit) && ((DFNU <= 1.0) || (AZ + AZ >= DFNU * DFNU)) {
        //-----------------------------------------------------------------------
        //     ASYMPTOTIC EXPANSION FOR LARGE Z
        //-----------------------------------------------------------------------
        let (cy, nw) = asymptotic_i(z, order, KODE, NN)?;
        debug_assert!(nw == NZ);
        return Ok((cy, NZ));
    }
    let mut skip_az_rl_check = true;
    if DFNU > 1.0 {
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
        DFNU = order + ((NN - 1) as f64);
    }
    if (DFNU > MACHINE_CONSTANTS.asymptotic_order_limit)
        || (AZ > MACHINE_CONSTANTS.asymptotic_order_limit)
    {
        //-----------------------------------------------------------------------
        //     INCREMENT FNU+NN-1 UP TO FNUL, COMPUTE AND RECUR BACKWARD
        //-----------------------------------------------------------------------
        let NUI_isize = (MACHINE_CONSTANTS.asymptotic_order_limit - DFNU) as isize + 1;
        let NUI = NUI_isize.max(0) as usize;
        let (NW, NLAST) = ZBUNI(z, order, KODE, NN, NUI, &mut cy)?;
        NZ += NW;
        if NLAST == 0 {
            return Ok((cy, NZ));
        }
        NN = NLAST;
    }
    if !skip_az_rl_check && AZ <= MACHINE_CONSTANTS.asymptotic_z_limit {
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
        check_underflow_uniform_asymp_params(z, order, KODE, IKType::K, 2, &mut [c_one(); 2])
    {
        if NW > 0 {
            Err(Overflow)
        } else {
            let nz = i_wronksian(z, order, KODE, NN, &mut cy)?;
            Ok((cy, nz))
        }
    } else {
        Ok((vec![c_one(); NN], NN))
    }
}

pub fn ZACAI(
    z: Complex64,
    order: f64,
    KODE: Scaling,
    rotation: RotationDirection,
    N: usize,
) -> BesselResult {
    // ***BEGIN PROLOGUE  ZACAI
    // ***REFER TO  ZAIRY
    //
    //     ZACAI APPLIES THE ANALYTIC CONTINUATION FORMULA
    //
    //         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
    //                 MP=PI*MR*CMPLX(0.0,1.0)
    //
    //     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT
    //     HALF Z PLANE FOR USE WITH ZAIRY WHERE FNU=1/3 OR 2/3 AND N=1.
    //     ZACAI IS THE SAME AS ZACON WITH THE PARTS FOR LARGER ORDERS AND
    //     RECURRENCE REMOVED. A RECURSIVE CALL TO ZACON CAN RESULT if ZACON
    //     IS CALLED FROM ZAIRY.
    //

    let mut NZ = 0;
    let zn = -z;
    let AZ = z.abs();
    let NN = N;
    let DFNU = order + ((N - 1) as f64);
    let (mut y, _) = if (AZ * AZ * 0.25 <= DFNU + 1.0) || (AZ <= 2.0) {
        //-----------------------------------------------------------------------
        //     POWER SERIES FOR THE I FUNCTION
        //-----------------------------------------------------------------------
        let (y, NW_signed) = i_power_series(zn, order, KODE, NN)?;
        debug_assert!(NW_signed >= 0);
        (y, NW_signed.unsigned_abs())
    } else if AZ >= MACHINE_CONSTANTS.asymptotic_z_limit {
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
    let SGN = -PI * rotation.signum();
    let mut csgn = Complex64::new(0.0, SGN);
    if KODE == Scaling::Scaled {
        csgn = Complex64::I * csgn.im * Complex64::cis(-zn.im);
    }
    //-----------------------------------------------------------------------
    //     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
    //     WHEN FNU IS LARGE
    //-----------------------------------------------------------------------
    let INU = order as usize;
    let mut cspn = Complex64::cis(order.fract() * SGN);
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

fn ZUNK1(
    z: Complex64,
    order: f64,
    scaling: Scaling,
    rotation: RotationDirection,
    n: usize,
) -> BesselResult {
    //     ZUNK1 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE
    //     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE
    //     UNIFORM ASYMPTOTIC EXPANSION.
    //     MR INDICATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION.

    let mut found_one_good_entry = false;
    let mut nz = 0;
    //-----------------------------------------------------------------------
    //     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION GREATER THAN
    //     THE UNDERFLOW LIMIT
    //-----------------------------------------------------------------------
    let zr = if z.re < 0.0 { -z } else { z };
    let mut phi = [c_zero(); 2];
    let mut zeta1 = [c_zero(); 2];
    let mut zeta2 = [c_zero(); 2];
    let mut sum = [c_zero(); 2];
    let mut cy = [c_zero(); 2];
    let mut n_elements_set = 0;
    let mut y = c_zeros(n);
    let mut k_overflow_state = Overflow::NearUnder;

    let mut j = 1;
    for i in 0..n {
        n_elements_set = i + 1;
        // j flip-flops between 0 and 1 using j = 1-j
        j = 1 - j;
        let modified_order = order + (i as f64);

        let sum_opt;
        (phi[j], zeta1[j], zeta2[j], sum_opt) = zunik(zr, modified_order, IKType::K, false);
        sum[j] = sum_opt.unwrap();
        let mut s1 = -scaling.scale_zetas(zr, modified_order, zeta1[j], zeta2[j]);
        let of = Overflow::find_overflow(s1.re, phi[j], 0.0);
        if !found_one_good_entry {
            k_overflow_state = of;
        }
        match of {
            Overflow::Over(_) => return Err(Overflow),
            Overflow::Under(_) => {
                if z.re < 0.0 {
                    return Err(Overflow);
                }
                found_one_good_entry = false;
                y[i] = c_zero();
                nz += 1;
            }
            Overflow::None | Overflow::NearOver | Overflow::NearUnder => {
                //-----------------------------------------------------------------------
                //     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
                //     EXPONENT EXTREMES
                //-----------------------------------------------------------------------
                let mut s2 = phi[j] * sum[j];
                s1 = MACHINE_CONSTANTS.scaling_factors[k_overflow_state] * s1.exp();
                s2 *= s1;
                let will_underflow = will_underflow(
                    s2,
                    MACHINE_CONSTANTS.smallness_threshold[0],
                    MACHINE_CONSTANTS.abs_error_tolerance,
                );
                if k_overflow_state != Overflow::NearUnder || !will_underflow {
                    cy[found_one_good_entry as usize] = s2;
                    y[i] = s2 * MACHINE_CONSTANTS.reciprocal_scaling_factors[k_overflow_state];
                    if found_one_good_entry {
                        break;
                    }
                    found_one_good_entry = true;
                } else if will_underflow {
                    if z.re < 0.0 {
                        return Err(Overflow);
                    }
                    y[i] = c_zero();
                    nz += 1;
                    if i > 0 && y[i - 1] != c_zero() {
                        y[i - 1] = c_zero();
                        nz += 1
                    }
                }
            }
        };
    }

    let rz = 2.0 * zr.conj() / zr.abs().powi(2);
    if n_elements_set < n {
        //-----------------------------------------------------------------------
        //     TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW. SET SEQUENCE TO ZERO
        //     ON UNDERFLOW.
        //-----------------------------------------------------------------------
        let modified_order = order + ((n - 1) as f64);
        let (phi, zet1d, zet2d, _sumd) = zunik(
            zr,
            modified_order,
            IKType::K,
            rotation == RotationDirection::None,
        );
        let overflow_test = -scaling.scale_zetas(zr, modified_order, zet1d, zet2d);

        match Overflow::find_overflow(overflow_test.re.abs(), phi, 0.0) {
            Overflow::Over(_) => return Err(Overflow),
            Overflow::Under(_) => {
                return if z.re < 0.0 {
                    Err(Overflow)
                } else {
                    Ok((vec![c_zero(); n], n))
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
        if rotation == RotationDirection::None {
            return Ok((y, nz));
        }
        //-----------------------------------------------------------------------
        //     ANALYTIC CONTINUATION FOR RE(Z) < 0.0
        //-----------------------------------------------------------------------
        nz = 0;
        let rotation_angle = -PI * rotation.signum();
        //-----------------------------------------------------------------------
        //     CSPN AND CSGN ARE COEFF OF K AND I FUNCTIONS RESP.
        //-----------------------------------------------------------------------

        let integer_order = order as i64;
        let order_frac = order.fract();
        let modified_int_order = integer_order + (n as i64) - 1;
        let mut cspn = Complex64::cis(order_frac * rotation_angle);
        if (modified_int_order % 2) != 0 {
            cspn = -cspn;
        }
        let mut dummy_n_good = 0;
        let mut found_one_good_entry = false;
        let mut i_overflow_state = Overflow::None;
        let mut remaining_n = n;
        for (i, yi) in y.iter_mut().enumerate().rev() {
            remaining_n = i;
            let modified_order = order + (i as f64);
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
            let of = Overflow::find_overflow(s1.re, phid, 0.0);
            if !found_one_good_entry && !matches!(of, Overflow::Under(_)) {
                i_overflow_state = of;
            }
            let mut s2 = match of {
                Overflow::Over(_) => {
                    return Err(Overflow);
                }
                Overflow::Under(_) => c_zero(),
                Overflow::NearOver | Overflow::NearUnder | Overflow::None => {
                    let st = phid * sumd;
                    let mut s2 = Complex64::I * st * rotation_angle;
                    s1 = s1.exp() * MACHINE_CONSTANTS.scaling_factors[i_overflow_state];
                    s2 *= s1;
                    if i_overflow_state == Overflow::NearUnder
                        && will_underflow(
                            s2,
                            MACHINE_CONSTANTS.smallness_threshold[0],
                            MACHINE_CONSTANTS.abs_error_tolerance,
                        )
                    {
                        s2 = c_zero();
                    }
                    s2
                }
            };
            cy[found_one_good_entry as usize] = s2;
            let c2 = s2;
            s2 *= MACHINE_CONSTANTS.reciprocal_scaling_factors[i_overflow_state];
            //-----------------------------------------------------------------------
            //     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N
            //-----------------------------------------------------------------------
            s1 = *yi;
            if scaling == Scaling::Scaled {
                nz += underflow_add_i_k(zr, &mut s1, &mut s2, &mut dummy_n_good);
            }
            *yi = s1 * cspn + s2;
            cspn = -cspn;
            if c2 == c_zero() {
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
                MACHINE_CONSTANTS.reciprocal_scaling_factors[i_overflow_state];
            let mut ASCLE = MACHINE_CONSTANTS.smallness_threshold[i_overflow_state];
            for (i, yi) in y.iter_mut().enumerate().take(remaining_n).rev() {
                let modified_order = order + (i + 1) as f64;
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
                if max_abs_component(unscaled_s2) <= ASCLE {
                    continue;
                }
                i_overflow_state.increment();
                ASCLE = MACHINE_CONSTANTS.smallness_threshold[i_overflow_state];
                s1 *= reciprocal_scale_factor;
                s2 = ck; // ck is previously calculated s2 * reciprocal_scale_factor
                s1 *= MACHINE_CONSTANTS.scaling_factors[i_overflow_state];
                s2 *= MACHINE_CONSTANTS.scaling_factors[i_overflow_state];
                reciprocal_scale_factor =
                    MACHINE_CONSTANTS.reciprocal_scaling_factors[i_overflow_state];
            }
        }
    }
    Ok((y, nz))
}

fn ZUNK2(
    z: Complex64,
    order: f64,
    scaling: Scaling,
    rotation: RotationDirection,
    n: usize,
) -> BesselResult {
    //     ZUNK2 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE
    //     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE
    //     UNIFORM ASYMPTOTIC EXPANSIONS FOR H(KIND,FNU,ZN) AND J(FNU,ZN)
    //     WHERE ZN IS IN THE RIGHT HALF PLANE, KIND=(3-MR)/2, MR=+1 OR
    //     -1. HERE ZN=ZR*I OR -ZR*I WHERE ZR=Z if Z IS IN THE RIGHT
    //     HALF PLANE OR ZR=-Z if Z IS IN THE LEFT HALF PLANE. MR INDIC-
    //     ATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION.

    const CR1: Complex64 = Complex64::new(1.0, 1.73205080756887729);
    const CR2: Complex64 = Complex64::new(-0.5, -8.66025403784438647e-01);

    let mut found_one_good_entry = false;
    let mut nz = 0;
    let mut y = c_zeros(n);
    let zr = if z.re < 0.0 { -z } else { z };
    let mut zn = -Complex64::I * zr;
    let mut zb = zr;
    let integer_order = order as usize;
    let order_fract = order.fract();
    let ANG = -FRAC_PI_2 * order_fract;
    let c2 = -Complex64::I * Complex64::from_polar(FRAC_PI_2, ANG);
    let mut cs = CR1 * c2 * CIP[integer_order % 4].conj();
    if zr.im <= 0.0 {
        zn.re = -zn.re;
        zb.im = -zb.im;
    }
    //-----------------------------------------------------------------------
    //     K(FNU,Z) IS COMPUTED FROM H(2,FNU,-I*Z) WHERE Z IS IN THE FIRST
    //     QUADRANT. FOURTH QUADRANT VALUES (YY <= 0.0E0) ARE COMPUTED BY
    //     CONJUGATION SINCE THE K FUNCTION IS REAL ON THE POSITIVE REAL AXIS
    //-----------------------------------------------------------------------

    let mut phi = [c_zero(); 2];
    let mut arg = [c_zero(); 2];
    let mut zeta1 = [c_zero(); 2];
    let mut zeta2 = [c_zero(); 2];
    let mut asum = [None; 2];
    let mut bsum = [None; 2];
    let mut cy = [c_zero(); 2];
    let mut j = 1;
    let mut overflow_state_k = Overflow::None;
    let mut n_elements_set = 0;

    for i in 0..n {
        n_elements_set = i + 1;
        // j flip-flops between 0 and 1 using  = 1-j
        j = 1 - j;
        let modified_order = order + (i as f64);
        (phi[j], arg[j], zeta1[j], zeta2[j], asum[j], bsum[j]) = zunhj(zn, modified_order, false);
        let s1 = -scaling.scale_zetas(zb, modified_order, zeta1[j], zeta2[j]);
        let of = Overflow::find_overflow(s1.re, phi[j], -0.25 * arg[j].abs().ln() - AIC);

        let mut handle_underflow = |of_already: &mut bool, cs_: &mut Complex64| {
            //-----------------------------------------------------------------------
            //     FOR ZR < 0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
            //-----------------------------------------------------------------------
            if z.re < 0.0 {
                return Err(Overflow);
            }
            *of_already = false;
            y[i] = c_zero();
            nz += 1;
            *cs_ *= -Complex64::I;
            if i != 0 && y[i - 1] != c_zero() {
                y[i - 1] = c_zero();
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
                let s1 = s1.exp() * MACHINE_CONSTANTS.scaling_factors[overflow_state_k];
                s2 *= s1;
                if overflow_state_k == Overflow::NearUnder
                    && will_underflow(
                        s2,
                        MACHINE_CONSTANTS.smallness_threshold[0],
                        MACHINE_CONSTANTS.abs_error_tolerance,
                    )
                {
                    handle_underflow(&mut found_one_good_entry, &mut cs)?
                }
                if zr.im <= 0.0 {
                    s2 = s2.conj();
                }
                cy[found_one_good_entry as usize] = s2;
                y[i] = s2 * MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state_k];
                cs = -Complex64::I * cs;
                if found_one_good_entry {
                    break;
                }
                found_one_good_entry = true;
            }
        };
    }

    let rz = 2.0 * zr.conj() / zr.abs().powi(2);
    let mut phid = c_zero();
    let mut argd = c_zero();
    let mut zeta1d = c_zero();
    let mut zeta2d = c_zero();
    let mut asumd = None;
    let mut bsumd = None;
    let do_overflow_check = n_elements_set < n;
    if do_overflow_check {
        //-----------------------------------------------------------------------;
        //     TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW. SET SEQUENCE TO ZERO;
        //     ON UNDERFLOW.;
        //-----------------------------------------------------------------------;
        let modified_order = order + ((n - 1) as f64);
        (phid, argd, zeta1d, zeta2d, asumd, bsumd) =
            zunhj(zn, modified_order, rotation == RotationDirection::None);
        let s1 = -scaling.scale_zetas(zb, modified_order, zeta1d, zeta2d);
        match Overflow::find_overflow(s1.re, phid, 0.0) {
            Overflow::Over(_) => return Err(Overflow),

            Overflow::Under(_) => {
                if z.re < 0.0 {
                    return Err(Overflow);
                }
                return Ok((c_zeros(n), nz));
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
    let sgn = -PI * rotation.signum();
    //-----------------------------------------------------------------------
    //     CSPN AND CSGN ARE COEFF OF K AND I FUNCIONS RESP.
    //-----------------------------------------------------------------------
    let csgn = if zr.im <= 0.0 { -sgn } else { sgn };
    let modified_integer_order = integer_order + n - 1;
    let mut cspn = Complex64::cis(order_fract * sgn);
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
    let cos_sin = Complex64::cis(ANG);
    // let mut cs = Complex64::I * Complex64::from_polar(CSGNI, ANG);
    let mut cs = csgn * Complex64::new(cos_sin.im, cos_sin.re);
    cs *= CIP[modified_integer_order % 4];
    let mut IUF = 0;

    found_one_good_entry = false;
    let mut overflow_state_i = Overflow::None;
    let mut remaining_n = n;
    for (i, yi) in y.iter_mut().enumerate().rev() {
        remaining_n = i;
        let modified_order = order + (i as f64);
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

        let of = Overflow::find_overflow(s1.re, phid, -0.25 * argd.abs().ln() - AIC);
        if !found_one_good_entry {
            overflow_state_i = if matches!(of, Overflow::Under(_)) {
                Overflow::None
            } else {
                of
            };
        }
        let mut s2 = match of {
            Overflow::Over(_) => return Err(Overflow),
            Overflow::Under(_) => c_zero(),
            Overflow::NearOver | Overflow::None | Overflow::NearUnder => {
                let (airy, d_airy) = airy_pair(argd);
                let pt = ((d_airy * bsumd.unwrap()) + (airy * asumd.unwrap())) * phid;
                let mut s2 = pt * cs;
                s1 = s1.exp() * MACHINE_CONSTANTS.scaling_factors[overflow_state_i];
                s2 *= s1;
                if overflow_state_i == Overflow::NearUnder
                    && will_underflow(
                        s2,
                        MACHINE_CONSTANTS.smallness_threshold[0],
                        MACHINE_CONSTANTS.abs_error_tolerance,
                    )
                {
                    s2 = c_zero();
                }
                s2
            }
        };
        if zr.im <= 0.0 {
            s2 = s2.conj();
        }
        cy[found_one_good_entry as usize] = s2;
        let c2 = s2;
        s2 *= MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state_i];
        //-----------------------------------------------------------------------;
        //     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N;
        //-----------------------------------------------------------------------;
        s1 = *yi;
        if scaling == Scaling::Scaled {
            nz += underflow_add_i_k(zr, &mut s1, &mut s2, &mut IUF);
        }
        *yi = s1 * cspn + s2;
        cspn = -cspn;
        cs *= -Complex64::I;
        if c2 == c_zero() {
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

        let mut recip_scale_factor = MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state_i];
        let mut ascle = MACHINE_CONSTANTS.smallness_threshold[overflow_state_i];
        // TODO recurr with assignment fn
        for (i, yi) in y.iter_mut().enumerate().take(remaining_n).rev() {
            let modified_order = order + (i + 1) as f64;
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
                ascle = MACHINE_CONSTANTS.smallness_threshold[overflow_state_i];
                s1 *= recip_scale_factor;
                s2 = old_c2;
                s1 *= MACHINE_CONSTANTS.scaling_factors[overflow_state_i];
                s2 *= MACHINE_CONSTANTS.scaling_factors[overflow_state_i];
                recip_scale_factor = MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state_i];
            }
        }
    }
    Ok((y, nz))
}

fn ZBUNI(
    z: Complex64,
    order: f64,
    KODE: Scaling,
    N: usize,
    NUI: usize,
    y: &mut [Complex64],
) -> Result<(usize, usize), BesselError> {
    // ***BEGIN PROLOGUE  ZBUNI
    // ***REFER TO  ZBESI,ZBESK
    //
    //     ZBUNI COMPUTES THE I BESSEL FUNCTION FOR LARGE CABS(Z) >
    //     FNUL AND FNU+N-1 < FNUL. THE ORDER IS INCREASED FROM
    //     FNU+N-1 GREATER THAN FNUL BY ADDING NUI AND COMPUTING
    //     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR I(FNU,Z)
    //     ON IFORM=1 AND THE EXPANSION FOR J(FNU,Z) ON IFORM=2
    //
    // ***ROUTINES CALLED  ZUNI1,ZUNI2,ZABS,d1mach
    // ***END PROLOGUE  ZBUNI

    let imaginary_dominant = imaginary_dominant(z);
    if NUI != 0 {
        let mut FNUI = NUI as f64;
        let DFNU = order + ((N - 1) as f64);
        let GNU = DFNU + FNUI;
        let mut cy = [c_zero(); 2];
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
            return Ok((0, N));
        }
        //----------------------------------------------------------------------
        //     SCALE BACKWARD RECURRENCE, BRY(3) IS DEFINED BUT NEVER USED
        //----------------------------------------------------------------------
        let (mut overflow_state, mut ASCLE, mut CSCLR) =
            if cy[0].abs() <= MACHINE_CONSTANTS.smallness_threshold[0] {
                (
                    Overflow::NearUnder,
                    MACHINE_CONSTANTS.smallness_threshold[0],
                    1.0 / MACHINE_CONSTANTS.abs_error_tolerance,
                )
            } else if cy[0].abs() >= MACHINE_CONSTANTS.smallness_threshold[1] {
                (
                    Overflow::NearOver,
                    MACHINE_CONSTANTS.smallness_threshold[2],
                    MACHINE_CONSTANTS.abs_error_tolerance,
                )
            } else {
                (
                    Overflow::None,
                    MACHINE_CONSTANTS.smallness_threshold[1],
                    1.0,
                )
            };

        let mut CSCRR = 1.0 / CSCLR;
        let mut s1 = cy[1] * CSCLR;
        let mut s2 = cy[0] * CSCLR;
        // working out rz in multiple steps seems to give different floating point answer.
        let rz = 2.0 * z.conj() / z.abs().pow(2);

        for _ in 0..NUI {
            let st = s2;
            s2 = (DFNU + FNUI) * rz * s2 + s1;
            s1 = st;
            FNUI -= 1.0;
            if overflow_state == Overflow::NearOver {
                continue;
            }
            let st = s2 * CSCRR;
            if max_abs_component(st) <= ASCLE {
                continue;
            }
            overflow_state.increment();
            ASCLE = MACHINE_CONSTANTS.smallness_threshold[overflow_state];
            s1 *= CSCRR;
            s2 = st;
            CSCLR *= MACHINE_CONSTANTS.abs_error_tolerance;
            CSCRR = 1.0 / CSCLR;
            s1 *= CSCLR;
            s2 *= CSCLR;
        }
        y[N - 1] = s2 * CSCRR;
        if N == 1 {
            return Ok((0, NLAST));
        }
        let NL = N - 1;
        FNUI = NL as f64;
        let mut K = NL;
        for _ in 0..NL {
            let st = s2;
            s2 = (order + FNUI) * (rz * s2) + s1;
            s1 = st;
            y[K - 1] = s2 * CSCRR;
            FNUI -= 1.0;
            K -= 1;
            if overflow_state == Overflow::NearOver {
                continue;
            }
            // using K (rather than K-1) below as Amos "saved" the y value before K was decremented
            if max_abs_component(y[K]) <= ASCLE {
                continue;
            }
            overflow_state.increment();
            ASCLE = MACHINE_CONSTANTS.smallness_threshold[overflow_state];
            s1 *= CSCRR;
            s2 = y[K - 1];
            CSCLR *= MACHINE_CONSTANTS.abs_error_tolerance;
            CSCRR = 1.0 / CSCLR;
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
        ZUNI2(z, order, KODE, N, y)?
    } else {
        //-----------------------------------------------------------------------
        //     ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN
        //     -PI/3 <= ARG(Z) <= PI/3
        //-----------------------------------------------------------------------
        ZUNI1(z, order, KODE, N, y)?
    };

    Ok((NW, NLAST))
}

fn ZUNI1(
    z: Complex64,
    order: f64,
    scaling: Scaling,
    n: usize,
    y: &mut [Complex64],
) -> BesselResult<(usize, usize)> {
    // ***BEGIN PROLOGUE  ZUNI1
    // ***REFER TO  ZBESI,ZBESK
    //
    //     ZUNI1 COMPUTES I(FNU,Z)  BY MEANS OF THE UNIFORM ASYMPTOTIC
    //     EXPANSION FOR I(FNU,Z) IN -PI/3 <= ARG Z <= PI/3.
    //
    //     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
    //     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
    //     NLAST != 0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
    //     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1 < FNUL.
    //     Y(I)=CZERO FOR I=NLAST+1,N

    let mut nz = 0;
    let mut n_remaining = n;
    //-----------------------------------------------------------------------
    //     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
    //-----------------------------------------------------------------------
    let mut modified_order = order.max(1.0);
    let (_, zeta1, zeta2, _) = zunik(z, modified_order, IKType::I, true);
    let s1 = scaling.scale_zetas(z, modified_order, zeta1, zeta2);
    // phi is chosen here for refined tests to equal the original tests
    // which don't test refinement
    match Overflow::find_overflow(s1.re, c_one(), 0.0) {
        Overflow::Over(_) => return Err(Overflow),
        Overflow::Under(_) => return Ok((n, 0)),
        _ => (),
    }
    let mut overflow_state = Overflow::None; // this value should never be used
    let mut cy = [c_zero(); 2];
    let mut handle_underflow = |n_remaining: &mut usize,
                                y: &mut [Complex64]|
     -> BesselResult<bool> {
        //-----------------------------------------------------------------------
        //     SET UNDERFLOW AND UPDATE PARAMETERS
        //-----------------------------------------------------------------------
        y[*n_remaining - 1] = c_zero();
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
        let modified_order = order + ((*n_remaining - 1) as f64);
        if modified_order < MACHINE_CONSTANTS.asymptotic_order_limit {
            return Ok(true);
        }
        Ok(false)
    };

    'outer: loop {
        for i in 0..2.min(n_remaining) {
            modified_order = order + ((n_remaining - (i + 1)) as f64);
            let (phi, zeta1, zeta2, sum) = zunik(z, modified_order, IKType::I, false);
            let sum = sum.unwrap();
            let mut s1 = scaling.scale_zetas(z, modified_order, zeta1, zeta2);
            if scaling == Scaling::Scaled {
                s1 += Complex64::new(0.0, z.im);
            }

            let of = Overflow::find_overflow(s1.re, phi, 0.0);
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
            s1 = MACHINE_CONSTANTS.scaling_factors[overflow_state] * s1.exp();
            s2 *= s1;
            if overflow_state == Overflow::NearUnder
                && will_underflow(
                    s2,
                    MACHINE_CONSTANTS.smallness_threshold[0],
                    MACHINE_CONSTANTS.abs_error_tolerance,
                )
            {
                if handle_underflow(&mut n_remaining, y)? {
                    return Ok((nz, n_remaining));
                }
                continue 'outer;
            }
            cy[i] = s2;
            y[n_remaining - i - 1] =
                s2 * MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
        }
        break 'outer;
    }
    if n_remaining > 2 {
        let [s1, s2] = cy;
        recurr(false, order, z, y, n_remaining - 2, s1, s2, overflow_state);
    }
    Ok((nz, 0))
}

fn ZUNI2(
    z: Complex64,
    order: f64,
    scaling: Scaling,
    n: usize,
    y: &mut [Complex64],
) -> BesselResult<(usize, usize)> {
    // ***BEGIN PROLOGUE  ZUNI2
    // ***REFER TO  ZBESI,ZBESK
    //
    //     ZUNI2 COMPUTES I(FNU,Z) IN THE RIGHT HALF PLANE BY MEANS OF
    //     UNIFORM ASYMPTOTIC EXPANSION FOR J(FNU,ZN) WHERE ZN IS Z*I
    //     OR -Z*I AND ZN IS IN THE RIGHT HALF PLANE ALSO.
    //
    //     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
    //     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
    //     NLAST != 0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
    //     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1 < FNUL.
    //     Y(I)=CZERO FOR I=NLAST+1,N
    //
    // ***ROUTINES CALLED  ZAIRY,ZUunderflowCHK,ZUNHJ,ZUOIK,d1mach,ZABS
    // ***END PROLOGUE  ZUNI2

    let mut nz = 0;
    let mut n_remaining = n;

    //-----------------------------------------------------------------------
    //     ZN IS IN THE RIGHT HALF PLANE AFTER ROTATION BY CI OR -CI
    //-----------------------------------------------------------------------
    let mut zn = Complex64::new(z.im, -z.re);
    let mut zb = z;
    let integer_order = order as usize;

    let build_c2 = |effective_n: usize| {
        let index = (integer_order + effective_n - 1) % 4;
        let mut c2 = Complex64::cis(FRAC_PI_2 * order.fract()) * CIP[index];
        if z.im <= 0.0 {
            c2 = c2.conj();
        }
        c2
    };
    let mut c2 = build_c2(n);
    let sign_of_i = if z.im <= 0.0 {
        zn.re = -zn.re;
        zb.im = -zb.im;
        1.0
    } else {
        -1.0
    };
    //-----------------------------------------------------------------------
    //     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
    //-----------------------------------------------------------------------
    let mut modified_order = order.max(1.0);
    let (_, _, zeta1, zeta2, _, _) = zunhj(zn, modified_order, true);

    let s1 = scaling.scale_zetas(zb, modified_order, zeta1, zeta2);

    // phi is chosen here for refined tests to equal the original tests
    // which don't test refinement
    match Overflow::find_overflow(s1.re, c_one(), 0.0) {
        Overflow::Over(_) => return Err(Overflow),
        Overflow::Under(_) => return Ok((n, 0)),
        _ => (),
    }

    debug_assert!(modified_order + (n - 1) as f64 > MACHINE_CONSTANTS.asymptotic_order_limit);

    let mut overflow_state = Overflow::NearUnder;
    let mut cy = [c_zero(); 2];
    let mut handle_underflow = |n_remaining: &mut usize,
                                c2: &mut Complex64,
                                y: &mut [Complex64]|
     -> BesselResult<bool> {
        //-----------------------------------------------------------------------
        //     SET UNDERFLOW AND UPDATE PARAMETERS
        //-----------------------------------------------------------------------
        y[*n_remaining - 1] = c_zero();
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
        let modified_order = order + ((*n_remaining - 1) as f64);
        if modified_order < MACHINE_CONSTANTS.asymptotic_order_limit {
            return Ok(true);
        }
        *c2 = build_c2(*n_remaining);
        Ok(false)
    };
    'outer: loop {
        for i in 0..2.min(n_remaining) {
            modified_order = order + ((n_remaining - (i + 1)) as f64);
            let (phi, arg, zeta1, zeta2, asum, bsum) = zunhj(zn, modified_order, false);
            let asum = asum.unwrap();
            let bsum = bsum.unwrap();
            let mut s1 = scaling.scale_zetas(zb, modified_order, zeta1, zeta2);
            if scaling == Scaling::Scaled {
                s1 += Complex64::I * z.im.abs();
            }

            //-----------------------------------------------------------------------
            //     TEST FOR UNDERFLOW AND OVERFLOW
            //-----------------------------------------------------------------------
            let of = Overflow::find_overflow(s1.re, phi, -0.25 * arg.abs().ln() - AIC);
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
            let s1 = MACHINE_CONSTANTS.scaling_factors[overflow_state] * s1.exp();
            s2 *= s1;
            if overflow_state == Overflow::NearUnder
                && will_underflow(
                    s2,
                    MACHINE_CONSTANTS.smallness_threshold[0],
                    MACHINE_CONSTANTS.abs_error_tolerance,
                )
            {
                if handle_underflow(&mut n_remaining, &mut c2, y)? {
                    return Ok((nz, n_remaining));
                }
                continue 'outer;
            }
            if z.im <= 0.0 {
                s2 = s2.conj();
            }
            s2 *= c2;
            cy[i] = s2;
            y[n_remaining - i - 1] =
                s2 * MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
            c2 *= sign_of_i * Complex64::I;
        }
        break 'outer;
    }
    if n_remaining > 2 {
        let [s1, s2] = cy;
        recurr(false, order, z, y, n_remaining - 2, s1, s2, overflow_state);
    }
    Ok((nz, 0))
}

fn recurr(
    forward: bool,
    order: f64,
    z: Complex64,
    y: &mut [Complex64],
    n_offset: usize,
    mut s1: Complex64,
    mut s2: Complex64,
    mut overflow_state: Overflow,
) {
    let rz = 2.0 * z.conj() / z.abs().pow(2);

    let base_iterator = y.iter_mut().enumerate();
    let iterator = if forward {
        Either::Right(base_iterator.skip(n_offset))
    } else {
        Either::Left(base_iterator.take(n_offset).rev())
    };
    let index_adjustment = if forward { -1.0 } else { 1.0 };

    let mut recip_scale_factor = MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
    let mut boundary = MACHINE_CONSTANTS.smallness_threshold[overflow_state];

    for (i, yi) in iterator {
        let modified_order = order + (i as f64) + index_adjustment;
        (s1, s2) = (s2, s1 + modified_order * rz * s2);
        *yi = s2 * recip_scale_factor;
        if overflow_state != Overflow::NearOver && max_abs_component(*yi) > boundary {
            overflow_state.increment();
            boundary = MACHINE_CONSTANTS.smallness_threshold[overflow_state];
            s1 *= recip_scale_factor;
            s2 = *yi;
            s1 *= MACHINE_CONSTANTS.scaling_factors[overflow_state];
            s2 *= MACHINE_CONSTANTS.scaling_factors[overflow_state];
            recip_scale_factor = MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
        }
    }
}
