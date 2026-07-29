use num::{Complex, Zero, complex::ComplexFloat};

use crate::{
    BesselError, Scaling,
    amos::{
        IKType,
        asymptotics::i_asymp_large_order,
        asymptotics::i_asymptotic,
        gamma_ln, i_power_series,
        limits::{OverflowState, check_underflow_uniform_asymp_params},
        recurrence::{i_miller, scale_controlled_recurrence, scale_k_recurrence},
        utils::{two_over_z_safe, will_underflow},
        wronksian::i_wronksian,
    },
    types::{BesselFloat, BesselResult, BesselValues},
};

/// i_right_half_plane computes the i function in the right half z plane
///
/// i_right_half_plane  acts as a master dispatcher for computing the I_v (z)
/// sequence in the right half of the complex plane (Re (z) ≥ 0).
///
/// Because the Modified Bessel Function
///   I_v (z)
/// has wildly different behavioral regimes (growing exponentially for large z, decaying exponentially for large ν),
/// there is no single algorithm that can compute it accurately everywhere. This function checks the magnitude of
/// z and ν and routes the computation to the most numerically stable algorithm.
///
/// ### 1. The Power Series Regime (Small z)
///
/// If |z| is very small (|z| ≤ 2) or ν is relatively large compared to z, it attempts to use a standard Taylor/Power Series ( [i_power_series] ).
///
/// - If the series converges successfully, it returns.
/// - If the highest orders underflow to  0 , it subtracts them from  remaining_n  and passes the rest of the array down to the next algorithms.
///
/// ### 2. Large Argument Asymptotics (Large z, Small ν)
///
/// If |z| is massive (above the machine's asymptotic limit) and ν is relatively small (2|z| ≥ ν²),
/// it uses Hankel's asymptotic expansion for large arguments ( [i_asymptotic] ).
///
/// ### 3. Underflow Truncation
///
/// Before running any intermediate recurrences, if νₘₐₓ > 1, it runs a quick check ( [check_underflow_uniform_asymp_params] ).
/// Because I_v (z) decays extremely fast as ν → ∞, this calculates exactly how many of the highest-order requested values
/// will mathematically underflow to `0.0`.
/// It fills those with `0.0` , updates `remaining_n`, and saves the system from doing pointless, unstable math.
///
/// ### 4. Large Order Asymptotics (Extremely Large ν or z)
///
/// If the requested order or argument exceeds the machine's strict asymptotic limit, standard recurrences break down.
/// It delegates to [i_asymp_large_order], which uses Debye polynomials (Uniform Asymptotic Expansions) to bootstrap
/// the recurrence from a very high order.
///
/// ### 5. Miller's Algorithm (Intermediate Regime)
///
/// If it reaches this point, we are in the "messy middle" where asymptotics fail. It uses Miller's Algorithm,
/// which runs a backward recurrence from a dynamically calculated starting bound.
/// But because Miller's algorithm only gives relative (unscaled) values, they must be normalized.
/// It chooses between two normalization methods:
///
/// - Series Normalization ( [i_miller] ): If |z| is small enough, it normalizes the sequence by
///   plugging it into the known Neumann series identity ∑Iₙ(z) = eᶻ.
/// - Wronskian Normalization ( [i_wronksian] ): For larger |z|, the Neumann series
///   requires too many terms to converge. Instead, it computes K_v (z)
///   independently, and normalizes the I sequence using the Wronskian cross-product identity:
///   I_v K_v+1    + I_v+1   K_v  = 1/z
///
/// Originally ZBINU
pub(crate) fn i_right_half_plane<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    n: usize,
) -> BesselResult<T, usize> {
    let mut n_zeros = 0;
    let abs_z = z.abs();
    let mut remaining_n: usize = n;
    let mut max_order = order + T::from_usize(n - 1);
    let mut y = T::c_zeros(n);
    if abs_z <= T::two() || abs_z.powi(2) * T::from_f64(0.25) <= max_order + T::one() {
        // Power series for small z
        let n_zeros_inner;
        // i_power_series return *signed* n_zeros. As per the docs
        // n_zeros > 0 means that the last n_zeros components were set to zero
        // due to underflow. (As is the normal convention)
        // n_zeros < 0 means underflow occurred, but the
        // condition z.abs() <= 2*(order+1).sqrt() was violated and the
        // computation must be completed in another routine with n=n-abs(n_zeros).
        (y, n_zeros_inner) = i_power_series(z, order, scaling, remaining_n)?;
        let calculation_finished = n_zeros_inner >= 0;
        let n_to_complete: usize = n_zeros_inner.unsigned_abs();
        n_zeros += n_to_complete;
        remaining_n -= n_to_complete;
        if remaining_n == 0 || calculation_finished {
            return Ok((y, n_zeros));
        }
        max_order = order + (T::from_usize(remaining_n) - T::one());
    }

    if (abs_z >= T::MACHINE_CONSTANTS.asymptotic_z_limit)
        && ((max_order <= T::one()) || (max_order.powi(2) <= abs_z + abs_z))
    {
        // Large Argument Asymptotics (Large z, Small order)
        let (cy, n_zeros_asymptotic) = i_asymptotic(z, order, scaling, remaining_n)?;
        debug_assert!(n_zeros_asymptotic == n_zeros);
        return Ok((cy, n_zeros));
    }

    if max_order > T::one() {
        //-----------------------------------------------------------------------
        //     OVERFLOW AND UNDERFLOW TEST ON I SEQUENCE FOR MILLER ALGORITHM
        //-----------------------------------------------------------------------
        let n_zeros_underflow = check_underflow_uniform_asymp_params(
            z,
            order,
            scaling,
            IKType::I,
            remaining_n,
            &mut y,
        )?;
        n_zeros += n_zeros_underflow;
        remaining_n -= n_zeros_underflow;
        if remaining_n == 0 {
            return Ok((y, n_zeros));
        }
        max_order = order + T::from_usize(remaining_n - 1);
    }

    if (max_order > T::MACHINE_CONSTANTS.asymptotic_order_limit)
        || (abs_z > T::MACHINE_CONSTANTS.asymptotic_order_limit)
    {
        let (n_zeros_asymp_lo, remaining_n) =
            i_asymp_large_order(z, order, scaling, remaining_n, &mut y)?;
        n_zeros += n_zeros_asymp_lo;
        if remaining_n == 0 {
            return Ok((y, n_zeros));
        }
    }

    if max_order <= T::one() && abs_z <= T::MACHINE_CONSTANTS.asymptotic_z_limit {
        // Miller algorithm with series normalisation
        //-----------------------------------------------------------------------
        let y = i_miller(z, order, scaling, remaining_n)?;
        return Ok((y, n_zeros));
    }

    // Miller algorithm nomralised by the Wronksian
    let n_zeros_wr = i_wronksian(z, order, scaling, remaining_n, &mut y)?;
    Ok((y, n_zeros + n_zeros_wr))
}

/// k_right_half_plane computes the k bessel function in the right half z plane.
/// Originally ZBKNU
pub fn k_right_half_plane<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    n: usize,
) -> Result<BesselValues<T, usize>, BesselError<T>> {
    const KMAX: usize = 30;
    let RTFRAC_PI_2: T = T::from_f64(1.253_314_137_315_500_3);
    let SPI: T = T::from_f64(1.909_859_317_102_744);
    let FPI: T = T::from_f64(1.897_699_993_315_177_5);
    let CC: [T; 8] = [
        T::from_f64(5.772_156_649_015_329e-1),
        T::from_f64(-4.200_263_503_409_524e-2),
        T::from_f64(-4.219_773_455_554_433e-2),
        T::from_f64(7.218_943_246_663_1e-3),
        T::from_f64(-2.152_416_741_149_509_8e-4),
        T::from_f64(-2.013_485_478_078_824e-5),
        T::from_f64(1.133_027_231_981_696e-6),
        T::from_f64(6.116_095_104_481_416e-9),
    ];

    let abs_z = z.abs();
    let mut n_zeros = 0;
    let mut underflow_occurred = false;
    let mut overflow_state;
    let two_over_z = two_over_z_safe(z);
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
        let mut smu = two_over_z.ln();
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
            return Ok((vec![y], n_zeros));
        }

        //-----------------------------------------------------------------------;
        //     GENERATE K(DNU,Z) AND K(DNU+1,Z) FOR FORWARD RECURRENCE;
        //-----------------------------------------------------------------------;
        let (mut s1, mut s2) =
            k_right_half_plane_helper(z, frac_order_sqr, signed_fractional_order, f, p, q, ck);

        overflow_state =
            if (order + T::one()) * smu.re.abs() > T::MACHINE_CONSTANTS.approximation_limit {
                OverflowState::NearOver
            } else {
                OverflowState::None
            };
        s2 *= overflow_state.scaling_factor::<T>() * two_over_z;
        s1 *= overflow_state.scaling_factor::<T>();
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
        overflow_state = OverflowState::None;
        if scaling == Scaling::Unscaled {
            if z.re > T::MACHINE_CONSTANTS.approximation_limit {
                underflow_occurred = true;
                overflow_state = OverflowState::NearUnder;
            } else {
                coeff *= overflow_state.scaling_factor::<T>() * (-z).exp();
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
                        return Err(BesselError::DidNotConverge);
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
    let mut ck = (signed_fractional_order + T::one()) * two_over_z;
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

                let mut zd = z;
                let mut IC: isize = -1;
                let mut J = 1;
                for i in 0..integer_order {
                    n_tested = i + 2;
                    // TODO same calculation as other loops - this one is over different range and sets cy
                    // (so is designed to run until cy is set, and record this in INUB)
                    (s1, s2) = (s2, s2 * ck + s1);
                    ck += two_over_z;
                    let abs_ln_s2 = s2.abs().ln();
                    if -zd.re + abs_ln_s2 >= -T::MACHINE_CONSTANTS.exponent_limit {
                        let p1 = (-zd + s2.ln()).exp() / T::MACHINE_CONSTANTS.abs_error_tolerance;
                        if !will_underflow(p1) {
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
                overflow_state = OverflowState::NearUnder;

                s2 = cy[J];
                J = 1 - J;
                s1 = cy[J];
            }

            let mut P1R = overflow_state.reciprocal_scaling_factor::<T>();
            let mut ASCLE = overflow_state.boundary();
            for _ in n_tested..=integer_order {
                // TODO same loop as below?
                // TODO and same recurrence logic used in ZUNKX, ZUNIX?
                (s1, s2) = (s2, ck * s2 + s1);
                ck += two_over_z;
                let p2 = s2 * P1R;
                overflow_state.scale_recurrence(&mut s1, &mut s2, p2, &mut ASCLE, &mut P1R);
            }
        }
        if n == 1 {
            s1 = s2;
        }
    }

    let mut y = T::c_zeros(n);
    let n_completed = if !underflow_occurred {
        // ********* basic setup
        y[0] = s1 * overflow_state.reciprocal_scaling_factor::<T>();
        if n > 1 {
            y[1] = s2 * overflow_state.reciprocal_scaling_factor::<T>();
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
        scale_k_recurrence(z, order, n, &mut y, &mut n_zeros, two_over_z);
        let n_non_zero = (n - n_zeros) as isize;
        if n_non_zero <= 0 {
            return Ok((y, n_zeros));
        }
        let mut working_index = n_zeros;
        s1 = y[working_index];
        y[working_index] *= T::MACHINE_CONSTANTS.abs_error_tolerance;
        if n_non_zero > 1 {
            // if n_non_zero == 1 {
            //     return Ok((y, n_zeros));
            // }
            working_index += 1;
            s2 = y[working_index];
            y[working_index] *= T::MACHINE_CONSTANTS.abs_error_tolerance;
        }
        if n_non_zero > 2 {
            overflow_state = OverflowState::NearUnder;
        }
        working_index + 1
    };
    // End Setup
    if n_completed >= n {
        return Ok((y, n_zeros));
    }
    scale_controlled_recurrence(true, order, z, &mut y, n_completed, s1, s2, overflow_state);
    Ok((y, n_zeros))
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
