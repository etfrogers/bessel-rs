use std::f64::consts::FRAC_PI_2;

use num::{Complex, Integer, complex::ComplexFloat};

use super::consts::{
    AIRY_ASYMP_COEFFS_A, AIRY_ASYMP_COEFFS_B, AIRY_HJ_POLYNOMIAL_COEFFS,
    DEBYE_IK_POLYNOMIAL_COEFFS, IK_NORMALIZATION_FACTORS, TRANSITION_AIRY_A_COEFFS,
    TRANSITION_AIRY_B_COEFFS, TURNING_POINT_ZETA_COEFFS,
};
use crate::{
    BesselError, Scaling,
    amos::{
        IKType, PositiveArg, RotationDirection,
        airy::airy_pair,
        i_pow,
        limits::{OverflowState, check_underflow_uniform_asymp_params, underflow_add_i_k},
        max_abs_component,
        recurrence::scale_controlled_recurrence,
        utils::{AIC, two_over_z_safe, will_underflow},
    },
    types::{BesselFloat, BesselResult},
};

fn geom_underflow<T: BesselFloat>(z: Complex<T>, order: T) -> bool {
    let test = order * T::MACHINE_CONSTANTS.underflow_limit;
    z.re.abs() <= test && z.im.abs() <= test
}

fn underflow_zetas<T: BesselFloat>(order: T) -> (Complex<T>, Complex<T>) {
    let zeta1 = Complex::<T>::new(
        T::two() * T::MACHINE_CONSTANTS.underflow_limit.ln().abs() + order,
        T::zero(),
    );
    let zeta2 = Complex::<T>::new(order, T::zero());
    (zeta1, zeta2)
}

// ***BEGIN PROLOGUE  ZUNIK
// ***REFER TO  ZBESI,ZBESK
//
//        ZUNIK COMPUTES PARAMETERS FOR THE UNIFORM ASYMPTOTIC
//        EXPANSIONS OF THE I AND K FUNCTIONS ON IKFLG= 1 OR 2
//        RESPECTIVELY BY
//
//        W(FNU,ZR) = PHI*EXP(ZETA)*SUM
//
//        WHERE       ZETA=-ZETA1 + ZETA2       OR
//                          ZETA1 - ZETA2
//
//        THE FIRST CALL MUST HAVE INIT=0. SUBSEQUENT CALLS WITH THE
//        SAME ZR AND FNU WILL RETURN THE I OR K FUNCTION ON IKFLG=
//        1 OR 2 WITH NO CHANGE IN INIT. CWRK IS A COMPLEX WORK
//        ARRAY. IPMTR=0 COMPUTES ALL PARAMETERS. IPMTR=1 COMPUTES PHI,
//        ZETA1,ZETA2.
//
// ***ROUTINES CALLED  ZDIV,ZLOG,ZSQRT,d1mach
// ***END PROLOGUE  ZUNIK
//

/// Initaial for the uniform asymptotic Debye expansion of modified Bessel functions $I_\nu$ and $K_\nu$.
#[derive(Clone, Copy, Debug)]
pub(crate) struct DebyeGeometry<T: BesselFloat> {
    pub phi_i: Complex<T>,
    pub phi_k: Complex<T>,
    pub zeta1: Complex<T>,
    pub zeta2: Complex<T>,

    s: Complex<T>,
    sr: Complex<T>,
    reciprocal_order: T,
    is_underflow: bool,
}

impl<T: BesselFloat> DebyeGeometry<T> {
    pub fn compute(z: Complex<T>, order: T) -> Self {
        if geom_underflow(z, order) {
            let (zeta1, zeta2) = underflow_zetas(order);

            return Self {
                phi_i: T::C_ONE,
                phi_k: T::C_ONE,
                zeta1,
                zeta2,
                s: T::C_ONE,
                sr: T::C_ZERO,
                reciprocal_order: T::ZERO,
                is_underflow: true,
            };
        }

        let reciprocal_order = T::one() / order;
        let t = z * reciprocal_order;
        let s = T::C_ONE + t.powi(2);
        let s_root = s.sqrt();
        let zn = (T::C_ONE + s_root) / t;
        let zeta1 = zn.ln() * order;
        let zeta2 = s_root * order;
        let t = T::C_ONE / s_root;
        let sr = t * reciprocal_order;
        let sr_root = sr.sqrt();
        let phi_i = sr_root * T::from_f64(IK_NORMALIZATION_FACTORS[0]);
        let phi_k = sr_root * T::from_f64(IK_NORMALIZATION_FACTORS[1]);

        Self {
            phi_i,
            phi_k,
            zeta1,
            zeta2,
            s,
            sr,
            reciprocal_order,
            is_underflow: false,
        }
    }
}

/// Parameters for the uniform asymptotic Debye expansion of modified Bessel functions $I_\nu$ and $K_\nu$.
#[derive(Clone, Copy, Debug)]
pub(crate) struct DebyeParams<T: BesselFloat> {
    pub phi_i: Complex<T>,
    pub phi_k: Complex<T>,
    pub zeta1: Complex<T>,
    pub zeta2: Complex<T>,
    pub sum_i: Complex<T>, // Olver A(ζ) series
    pub sum_k: Complex<T>, // Olver B(ζ) series
}

impl<T: BesselFloat> DebyeParams<T> {
    pub fn compute(z: Complex<T>, order: T) -> Self {
        let geom = DebyeGeometry::compute(z, order);
        if geom.is_underflow {
            Self {
                phi_i: geom.phi_i,
                phi_k: geom.phi_k,
                zeta1: geom.zeta1,
                zeta2: geom.zeta2,
                sum_i: T::C_ZERO,
                sum_k: T::C_ZERO,
            }
        } else {
            let t2 = T::C_ONE / geom.s;
            let mut sum_i = T::C_ONE;
            let mut sum_k = T::C_ONE;
            let mut term_sign = T::one();
            let mut crfn = T::C_ONE;
            let mut ac = T::one();
            let mut l = 0;
            for k in 1..15 {
                let mut s = T::C_ZERO;
                for _ in 0..=k {
                    l += 1;
                    s = s * t2 + T::from_f64(DEBYE_IK_POLYNOMIAL_COEFFS[l]);
                }
                crfn *= geom.sr;
                let term = crfn * s;
                sum_i += term;
                term_sign = -term_sign;
                sum_k += term * term_sign;

                ac *= geom.reciprocal_order;
                let test = term.re.abs() + term.im.abs();
                if ac < T::MACHINE_CONSTANTS.abs_error_tolerance
                    && test < T::MACHINE_CONSTANTS.abs_error_tolerance
                {
                    break;
                }
            }

            Self {
                phi_i: geom.phi_i,
                phi_k: geom.phi_k,
                zeta1: geom.zeta1,
                zeta2: geom.zeta2,
                sum_i,
                sum_k,
            }
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct AiryGeometry<T: BesselFloat> {
    pub phi: Complex<T>,
    pub arg: Complex<T>, // Argument to the Airy function: ζ · ν^(2/3)
    pub zeta1: Complex<T>,
    pub zeta2: Complex<T>,

    state: AiryState<T>,
}

impl<T: BesselFloat> AiryGeometry<T> {
    const ONE_THIRD: f64 = 3.333_333_333_333_333e-1;
    const TWO_THIRDS: f64 = 6.666_666_666_666_666e-1;
    const THREE_PI_BY_2: f64 = 4.712_388_980_384_69;

    pub fn compute(z: Complex<T>, order: T) -> Self {
        if geom_underflow(z, order) {
            let (zeta1, zeta2) = underflow_zetas(order);
            let phi = T::C_ONE;
            let arg = T::C_ONE;
            return Self {
                phi,
                arg,
                zeta1,
                zeta2,
                state: AiryState::Underflow,
            };
        }
        let reciprocal_order = T::ONE / order;
        let z_over_order = z * reciprocal_order;
        let recip_order_sqr = reciprocal_order * reciprocal_order;
        //-----------------------------------------------------------------------
        //     COMPUTE IN THE FOURTH QUADRANT
        //-----------------------------------------------------------------------
        let fn13 = order.powf(T::from_f64(Self::ONE_THIRD));
        let fn23 = fn13 * fn13;
        let rfn13 = T::one() / fn13;

        let w2 = T::C_ONE - z_over_order.powi(2);
        let aw2 = w2.abs();

        let power_series = aw2 <= T::from_f64(0.25);
        if power_series {
            //-----------------------------------------------------------------------
            //     POWER SERIES FOR CABS(W2) <= 0.25
            //-----------------------------------------------------------------------
            let mut k_max = 1;
            let mut p = [T::C_ZERO; 30];
            let mut abs_p = [T::zero(); 30];
            p[0] = T::C_ONE;
            let mut suma = Complex::<T>::new(T::from_f64(TURNING_POINT_ZETA_COEFFS[0]), T::zero());
            abs_p[0] = T::one();
            if aw2 >= T::MACHINE_CONSTANTS.abs_error_tolerance {
                for k in 1..30 {
                    k_max = k + 1;
                    p[k] = p[k - 1] * w2;
                    suma += p[k] * T::from_f64(TURNING_POINT_ZETA_COEFFS[k]);
                    abs_p[k] = abs_p[k - 1] * aw2;
                    if abs_p[k] < T::MACHINE_CONSTANTS.abs_error_tolerance {
                        break;
                    }
                }
            }
            let zeta = w2 * suma;
            let arg = zeta * fn23;
            let mut za = suma.sqrt();
            let zeta2 = w2.sqrt() * order;
            let zeta1 = (T::C_ONE + T::from_f64(Self::TWO_THIRDS) * zeta * za) * zeta2;
            za *= T::two();
            let phi = za.sqrt() * rfn13;
            Self {
                phi,
                arg,
                zeta1,
                zeta2,
                state: AiryState::Transition {
                    p,
                    abs_p,
                    k_max,
                    recip_order_sqr,
                    reciprocal_order,
                    rfn13,
                },
            }
        } else {
            //-----------------------------------------------------------------------
            //     CABS(W2) > 0.25
            //-----------------------------------------------------------------------
            let mut w = w2.sqrt();
            if w.re < T::zero() {
                w.re = T::zero()
            };
            if w.im < T::zero() {
                w.im = T::zero()
            };

            let za = (T::C_ONE + w) / z_over_order;
            let mut zc = za.ln();
            zc.im = zc.im.clamp(T::zero(), T::from_f64(FRAC_PI_2));
            if zc.re < T::zero() {
                zc.re = T::zero()
            };
            let zth = (zc - w) * T::from_f64(1.5);
            let zeta1 = zc * order;
            let zeta2 = w * order;
            let azth = zth.abs();
            let mut ang = zth.parg();
            ang = ang.clamp(T::zero(), T::from_f64(Self::THREE_PI_BY_2));
            let pp = azth.powf(T::from_f64(Self::TWO_THIRDS));
            ang *= T::from_f64(Self::TWO_THIRDS);
            let mut zeta = Complex::<T>::cis(ang) * pp;
            if zeta.im < T::zero() {
                zeta.im = T::zero()
            };
            let arg = zeta * fn23;
            let rtzt = zth / zeta;
            let za = rtzt / w;
            let tazr = za + za;
            let phi = tazr.sqrt() * rfn13;

            Self {
                phi,
                arg,
                zeta1,
                zeta2,
                state: AiryState::Asymptotic {
                    w2,
                    aw2,
                    w,
                    zth,
                    azth,
                    rtzt,
                    reciprocal_order,
                    recip_order_sqr,
                    rfn13,
                },
            }
        }
    }
}

// ***BEGIN PROLOGUE  ZUNHJ
// ***REFER TO  ZBESI,ZBESK
//
//     REFERENCES
//         HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ AND I.A.
//         STEGUN, AMS55, NATIONAL BUREAU OF STANDARDS, 1965, CHAPTER 9.
//
//         ASYMPTOTICS AND SPECIAL FUNCTIONS BY F.W.J. OLVER, ACADEMIC
//         PRESS, N.Y., 1974, PAGE 420
//
//     ABSTRACT
//         ZUNHJ COMPUTES PARAMETERS FOR BESSEL FUNCTIONS C(FNU,Z) =
//         J(FNU,Z), Y(FNU,Z) OR H(I,FNU,Z) I=1,2 FOR LARGE ORDERS FNU
//         BY MEANS OF THE UNIFORM ASYMPTOTIC EXPANSION
//
//         C(FNU,Z)=C1*PHI*( ASUM*AIRY(ARG) + C2*BSUM*DAIRY(ARG) )
//
//         FOR PROPER CHOICES OF C1, C2, AIRY AND DAIRY WHERE AIRY IS
//         AN AIRY FUNCTION AND DAIRY IS ITS DERIVATIVE.
//
//               (2/3)*FNU*ZETA**1.5 = ZETA1-ZETA2,
//
//         ZETA1=0.5*FNU*CLOG((1+W)/(1-W)), ZETA2=FNU*W FOR SCALING
//         PURPOSES IN AIRY FUNCTIONS FROM CAIRY OR CBIRY.
//
//         MCONJ=SIGN OF AIMAG(Z), BUT IS AMBIGUOUS WHEN Z IS REAL AND
//         MUST BE SPECIFIED. IPMTR=0 RETURNS ALL PARAMETERS. IPMTR=
//         1 COMPUTES ALL EXCEPT ASUM AND BSUM.
//
// ***ROUTINES CALLED  ZABS,ZDIV,ZLOG,ZSQRT,d1mach
// ***END PROLOGUE  ZUNHJ

#[derive(Debug, Clone, Copy)]
pub(crate) struct AiryParams<T: BesselFloat> {
    pub phi: Complex<T>,
    pub arg: Complex<T>,
    pub zeta1: Complex<T>,
    pub zeta2: Complex<T>,
    pub asum: Complex<T>, // Olver A(ζ) series
    pub bsum: Complex<T>, // Olver B(ζ) series
}

#[derive(Debug, Clone, Copy)]
pub(crate) enum AiryState<T: BesselFloat> {
    Underflow,
    Transition {
        p: [Complex<T>; 30],
        abs_p: [T; 30],
        k_max: usize,
        reciprocal_order: T,
        recip_order_sqr: T,
        rfn13: T,
    },
    Asymptotic {
        w2: Complex<T>,
        aw2: T,
        w: Complex<T>,
        zth: Complex<T>,
        azth: T,
        rtzt: Complex<T>,
        reciprocal_order: T,
        recip_order_sqr: T,
        rfn13: T,
    },
}
impl<T: BesselFloat> AiryParams<T> {
    pub fn compute(z: Complex<T>, order: T) -> Self {
        let geom = AiryGeometry::compute(z, order);
        let phi = geom.phi;
        let arg = geom.arg;
        let zeta1 = geom.zeta1;
        let zeta2 = geom.zeta2;

        match geom.state {
            AiryState::Underflow => {
                let (zeta1, zeta2) = underflow_zetas(order);
                let phi = T::C_ONE;
                let arg = T::C_ONE;
                return Self {
                    phi,
                    arg,
                    zeta1,
                    zeta2,
                    asum: T::C_ZERO,
                    bsum: T::C_ZERO,
                };
            }
            AiryState::Transition {
                p,
                abs_p,
                k_max,
                recip_order_sqr,
                reciprocal_order,
                rfn13,
            } => {
                let sumb: Complex<T> = p[..k_max]
                    .iter()
                    .zip(TRANSITION_AIRY_B_COEFFS)
                    .map(|(p, b)| p * T::from_f64(b))
                    .sum();
                let mut asum = T::C_ZERO;
                let mut bsum = sumb;
                let mut l1 = 0;
                let mut l2 = 30;
                let btol =
                    T::MACHINE_CONSTANTS.abs_error_tolerance * (bsum.re.abs() + bsum.im.abs());
                let mut atol = T::MACHINE_CONSTANTS.abs_error_tolerance;
                let mut pp = T::one();
                let mut a_converged = false;
                let mut b_converged = false;
                if recip_order_sqr >= T::MACHINE_CONSTANTS.abs_error_tolerance {
                    for _ in 1..7 {
                        atol /= recip_order_sqr;
                        pp *= recip_order_sqr;
                        if !a_converged {
                            let mut suma = T::C_ZERO;
                            for k in 0..k_max {
                                suma += p[k] * T::from_f64(TRANSITION_AIRY_A_COEFFS[l1 + k]);
                                if abs_p[k] < atol {
                                    break;
                                }
                            }
                            asum += suma * pp;
                            if pp < T::MACHINE_CONSTANTS.abs_error_tolerance {
                                a_converged = true
                            };
                        }
                        if !b_converged {
                            let mut sumb = T::C_ZERO;
                            for k in 0..k_max {
                                sumb += p[k] * T::from_f64(TRANSITION_AIRY_B_COEFFS[l2 + k]);
                                if abs_p[k] < atol {
                                    break;
                                }
                            }
                            bsum += sumb * pp;
                            if pp < btol {
                                b_converged = true;
                            }
                        }
                        if a_converged && b_converged {
                            break;
                        }
                        l1 += 30;
                        l2 += 30;
                    }
                }
                asum += T::one();
                pp = reciprocal_order * rfn13;
                bsum *= pp;
                Self {
                    phi,
                    arg,
                    zeta1,
                    zeta2,
                    asum,
                    bsum,
                }
            }
            AiryState::Asymptotic {
                w2,
                aw2,
                w,
                zth,
                azth,
                rtzt,
                reciprocal_order,
                recip_order_sqr,
                rfn13,
            } => {
                let raw = T::one() / aw2.sqrt();
                let tfn = w.conj() * raw * raw * reciprocal_order;
                let razth = T::one() / azth;
                let rzth = zth.conj() * razth * razth * reciprocal_order;
                let zc = rzth * T::from_f64(AIRY_ASYMP_COEFFS_A[1]);
                let raw2 = T::one() / aw2;
                let t2 = w2.conj() * raw2 * raw2;
                let mut up = T::c_zeros(14);
                up[1] = (t2 * T::from_f64(AIRY_HJ_POLYNOMIAL_COEFFS[1])
                    + T::from_f64(AIRY_HJ_POLYNOMIAL_COEFFS[2]))
                    * tfn;
                let mut bsum = up[1] + zc;
                let mut asum = T::C_ZERO;
                if reciprocal_order >= T::MACHINE_CONSTANTS.abs_error_tolerance {
                    let mut przth = rzth;
                    let mut ptfn = tfn;
                    up[0] = T::C_ONE;
                    let mut pp = T::one();
                    let btol =
                        T::MACHINE_CONSTANTS.abs_error_tolerance * (bsum.re.abs() + bsum.im.abs());
                    let mut ks = 0;
                    let mut kp1 = 2;
                    let mut l = 2; //3;
                    let mut a_converged = false;
                    let mut b_converged = false;
                    let mut cr = T::c_zeros(14);
                    let mut dr = T::c_zeros(14);
                    for lr in (2..=12).step_by(2) {
                        let lrp1 = lr + 1;
                        //-----------------------------------------------------------------------
                        //     COMPUTE TWO ADDITIONAL CR, DR, AND UP FOR TWO MORE TERMS IN
                        //     NEXT SUMA AND SUMB
                        //-----------------------------------------------------------------------
                        for _k in lr..=lrp1 {
                            ks += 1;
                            kp1 += 1;
                            l += 1;
                            let mut za = Complex::<T>::new(
                                T::from_f64(AIRY_HJ_POLYNOMIAL_COEFFS[l]),
                                T::zero(),
                            );
                            for _ in 1..kp1 {
                                l += 1;
                                za = za * t2 + T::from_f64(AIRY_HJ_POLYNOMIAL_COEFFS[l]);
                            }
                            ptfn *= tfn;
                            up[kp1 - 1] = ptfn * za;
                            cr[ks - 1] = przth * T::from_f64(AIRY_ASYMP_COEFFS_B[ks]);
                            przth *= rzth;
                            dr[ks - 1] = przth * T::from_f64(AIRY_ASYMP_COEFFS_A[ks + 1]);
                        }
                        pp *= recip_order_sqr;
                        if !a_converged {
                            let mut suma = up[lrp1 - 1];
                            let mut ju = lrp1;
                            for cr_i in cr.iter().take(lr) {
                                ju -= 1;
                                suma += cr_i * up[ju - 1];
                            }
                            asum += suma;
                            let test = suma.re.abs() + suma.im.abs();
                            if pp < T::MACHINE_CONSTANTS.abs_error_tolerance
                                && test < T::MACHINE_CONSTANTS.abs_error_tolerance
                            {
                                a_converged = true
                            };
                        }
                        if !b_converged {
                            let mut sumb = up[lr + 1] + up[lrp1 - 1] * zc;
                            let mut ju = lrp1;
                            for jr_i in dr.iter().take(lr) {
                                ju -= 1;
                                sumb += jr_i * up[ju - 1];
                            }
                            bsum += sumb;
                            let test = sumb.re.abs() + sumb.im.abs();
                            if pp < btol && test < btol {
                                b_converged = true
                            };
                        }
                        if a_converged && b_converged {
                            break;
                        }
                    }
                }
                asum += T::C_ONE;
                bsum = (-bsum * rfn13) / rtzt;

                Self {
                    phi,
                    arg,
                    zeta1,
                    zeta2,
                    asum,
                    bsum,
                }
            }
        }
    }
}

/// i_uniform_asymp1 computes I(fnu,z)  by means of the uniform asymptotic
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
pub(crate) fn i_uniform_asymp1<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    n: usize,
    y: &mut [Complex<T>],
) -> Result<(usize, usize), BesselError<T>> {
    let mut n_zeros = 0;
    let mut n_remaining = n;
    //-----------------------------------------------------------------------
    //     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
    //-----------------------------------------------------------------------
    let mut modified_order = order.max(T::one());

    let DebyeGeometry { zeta1, zeta2, .. } = DebyeGeometry::compute(z, modified_order);

    // let (_, zeta1, zeta2, _) = ik_uniform_asymp_params(z, modified_order, IKType::I, true);
    let s1 = scaling.scale_zetas(z, modified_order, zeta1, zeta2);
    // phi is chosen here for refined tests to equal the original tests
    // which don't test refinement
    match OverflowState::check(s1.re, T::C_ONE, T::ZERO) {
        OverflowState::Over { .. } => return Err(BesselError::Overflow),
        OverflowState::Under { .. } => return Ok((n, 0)),
        _ => (),
    }
    let mut overflow_state = OverflowState::None; // this value should never be used
    let mut cy = [T::C_ZERO; 2];
    let mut handle_underflow = |n_remaining: &mut usize,
                                y: &mut [Complex<T>]|
     -> Result<bool, BesselError<T>> {
        //-----------------------------------------------------------------------
        //     SET UNDERFLOW AND UPDATE PARAMETERS
        //-----------------------------------------------------------------------
        y[*n_remaining - 1] = T::C_ZERO;
        n_zeros += 1;
        *n_remaining -= 1;
        if *n_remaining == 0 {
            return Ok(true);
        }
        let n_underflow =
            check_underflow_uniform_asymp_params(z, order, scaling, IKType::I, *n_remaining, y)?;
        *n_remaining -= n_underflow;
        n_zeros += n_underflow;
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
            let DebyeParams {
                phi_i: phi,
                zeta1,
                zeta2,
                sum_i: sum,
                ..
            } = DebyeParams::compute(z, modified_order);
            // let (phi, zeta1, zeta2, sum) =
            //     ik_uniform_asymp_params(z, modified_order, IKType::I, false);
            // let sum = sum.unwrap();
            let mut s1 = scaling.scale_zetas(z, modified_order, zeta1, zeta2);
            if scaling == Scaling::Scaled {
                s1 += Complex::<T>::new(T::ZERO, z.im);
            }

            let of = OverflowState::check(s1.re, phi, T::ZERO);
            if i == 0 {
                overflow_state = of;
            }
            match of {
                OverflowState::Over { .. } => return Err(BesselError::Overflow),
                OverflowState::Under { .. } => {
                    if handle_underflow(&mut n_remaining, y)? {
                        return Ok((n_zeros, n_remaining));
                    }
                    continue 'outer;
                }
                _ => (),
            }
            //-----------------------------------------------------------------------
            //     SCALE S1 if CABS(S1) < ASCLE
            //-----------------------------------------------------------------------
            let mut s2 = phi * sum;
            s1 = overflow_state.scaling_factor::<T>() * s1.exp();
            s2 *= s1;
            if overflow_state == OverflowState::NearUnder && will_underflow(s2) {
                if handle_underflow(&mut n_remaining, y)? {
                    return Ok((n_zeros, n_remaining));
                }
                continue 'outer;
            }
            cy[i] = s2;
            y[n_remaining - i - 1] = s2 * overflow_state.reciprocal_scaling_factor::<T>();
        }
        break 'outer;
    }
    if n_remaining > 2 {
        let [s1, s2] = cy;
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
        );
    }
    Ok((n_zeros, 0))
}

/// i_uniform_asymp2 computes I(fnu,z) in the right half plane by means of
/// uniform asymptotic expansion for J(fnu,zn) where zn is z*i
/// or -z*i and zn is in the right half plane also.
///
/// asymptotic_order_limit is the smallest order permitted for the asymptotic
/// expansion. nlast=0 means all of the y values were set.
/// nlast != 0 is the number left to be computed by another
/// formula for orders fnu to fnu+nlast-1 because fnu+nlast-1 < asymptotic_order_limit.
/// y(i)=czero for i=nlast+1,n
///
/// Originally ZUNI2
pub(crate) fn i_uniform_asymp2<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    n: usize,
    y: &mut [Complex<T>],
) -> Result<(usize, usize), BesselError<T>> {
    let mut n_zeros = 0;
    let mut n_remaining = n;

    //-----------------------------------------------------------------------
    //     ZN IS IN THE RIGHT HALF PLANE AFTER ROTATION BY CI OR -CI
    //-----------------------------------------------------------------------
    let mut zn = Complex::<T>::new(z.im, -z.re);
    let mut zb = z;
    let integer_order = order.to_usize().unwrap();

    let build_c2 = |effective_n: usize| {
        let mut c2 = Complex::<T>::cis(T::FRAC_PI_2() * order.fract())
            * i_pow(integer_order + effective_n - 1);
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
    let AiryGeometry { zeta1, zeta2, .. } = AiryGeometry::compute(zn, modified_order);

    let s1 = scaling.scale_zetas(zb, modified_order, zeta1, zeta2);

    // phi is chosen here for refined tests to equal the original tests
    // which don't test refinement
    match OverflowState::check(s1.re, T::C_ONE, T::zero()) {
        OverflowState::Over { .. } => return Err(BesselError::Overflow),
        OverflowState::Under { .. } => return Ok((n, 0)),
        _ => (),
    }

    debug_assert!(
        modified_order + T::from_usize(n - 1) > T::MACHINE_CONSTANTS.asymptotic_order_limit
    );

    let mut overflow_state = OverflowState::NearUnder;
    let mut cy = [T::C_ZERO; 2];
    let mut handle_underflow = |n_remaining: &mut usize,
                                c2: &mut Complex<T>,
                                y: &mut [Complex<T>]|
     -> Result<bool, BesselError<T>> {
        //-----------------------------------------------------------------------
        //     SET UNDERFLOW AND UPDATE PARAMETERS
        //-----------------------------------------------------------------------
        y[*n_remaining - 1] = T::C_ZERO;
        n_zeros += 1;
        *n_remaining -= 1;
        if *n_remaining == 0 {
            return Ok(true);
        }
        let n_underflow =
            check_underflow_uniform_asymp_params(z, order, scaling, IKType::I, *n_remaining, y)?;
        *n_remaining -= n_underflow;
        n_zeros += n_underflow;
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
            let AiryParams {
                phi,
                arg,
                zeta1,
                zeta2,
                asum,
                bsum,
                ..
            } = AiryParams::compute(zn, modified_order);

            // let (phi, arg, zeta1, zeta2, asum, bsum) =
            //     hj_uniform_asymp_params(zn, modified_order, false);
            // let asum = asum.unwrap();
            // let bsum = bsum.unwrap();
            let mut s1 = scaling.scale_zetas(zb, modified_order, zeta1, zeta2);
            if scaling == Scaling::Scaled {
                s1 += T::I * z.im.abs();
            }

            //-----------------------------------------------------------------------
            //     TEST FOR UNDERFLOW AND OVERFLOW
            //-----------------------------------------------------------------------
            let of = OverflowState::check(
                s1.re,
                phi,
                T::from_f64(-0.25) * arg.abs().ln() - T::from_f64(AIC),
            );
            if i == 0 {
                overflow_state = of;
            }
            match of {
                OverflowState::Over { .. } => return Err(BesselError::Overflow),
                OverflowState::Under { .. } => {
                    if handle_underflow(&mut n_remaining, &mut c2, y)? {
                        return Ok((n_zeros, n_remaining));
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
            let s1 = overflow_state.scaling_factor::<T>() * s1.exp();
            s2 *= s1;
            if overflow_state == OverflowState::NearUnder && will_underflow(s2) {
                if handle_underflow(&mut n_remaining, &mut c2, y)? {
                    return Ok((n_zeros, n_remaining));
                }
                continue 'outer;
            }
            if z.im <= T::ZERO {
                s2 = s2.conj();
            }
            s2 *= c2;
            cy[i] = s2;
            y[n_remaining - i - 1] = s2 * overflow_state.reciprocal_scaling_factor::<T>();
            c2 *= sign_of_i * T::I;
        }
        break 'outer;
    }
    if n_remaining > 2 {
        let [s1, s2] = cy;
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
        );
    }
    Ok((n_zeros, 0))
}

/// zunk1 computes K(fnu,z) and its analytic continuation from the
/// right half plane to the left half plane by means of the
/// uniform asymptotic expansion.
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
    let mut found_one_good_entry = false;
    let mut n_zeros = 0;
    //-----------------------------------------------------------------------
    //     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION GREATER THAN
    //     THE UNDERFLOW LIMIT
    //-----------------------------------------------------------------------
    let modified_z = if z.re < T::ZERO { -z } else { z };
    // let mut phi = [T::C_ZERO; 2];
    // let mut zeta1 = [T::C_ZERO; 2];
    // let mut zeta2 = [T::C_ZERO; 2];
    // let mut sum = [T::C_ZERO; 2];
    let mut debye_seeds: [Option<DebyeParams<T>>; 2] = [None, None];
    let mut cy = [T::C_ZERO; 2];
    let mut n_elements_set = 0;
    let mut y = T::c_zeros(n);
    let mut k_overflow_state = OverflowState::NearUnder;

    for i in 0..n {
        n_elements_set = i + 1;
        // j flip-flops between 0 and 1 using j = 1-j
        // j = 1 - j;
        let modified_order = order + T::from_usize(i);

        // Note: use modified_z so Re(z) >= 0
        let params = DebyeParams::compute(modified_z, modified_order);
        if i < 2 {
            debye_seeds[i] = Some(params);
        }

        // Use the K fields:
        let phi = params.phi_k;
        let sum = params.sum_k;
        let mut s1 = -scaling.scale_zetas(modified_z, modified_order, params.zeta1, params.zeta2);
        let of = OverflowState::check(s1.re, phi, T::ZERO);
        if !found_one_good_entry {
            k_overflow_state = of;
        }
        match of {
            OverflowState::Over { .. } => return Err(BesselError::Overflow),
            OverflowState::Under { .. } => {
                if z.re < T::ZERO {
                    return Err(BesselError::Overflow);
                }
                found_one_good_entry = false;
                y[i] = T::C_ZERO;
                n_zeros += 1;
            }
            OverflowState::None | OverflowState::NearOver | OverflowState::NearUnder => {
                //-----------------------------------------------------------------------
                //     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
                //     EXPONENT EXTREMES
                //-----------------------------------------------------------------------
                let mut s2 = phi * sum;
                s1 = k_overflow_state.scaling_factor::<T>() * s1.exp();
                s2 *= s1;
                let will_underflow = will_underflow(s2);
                if k_overflow_state != OverflowState::NearUnder || !will_underflow {
                    cy[found_one_good_entry as usize] = s2;
                    y[i] = s2 * k_overflow_state.reciprocal_scaling_factor::<T>();
                    if found_one_good_entry {
                        break;
                    }
                    found_one_good_entry = true;
                } else if will_underflow {
                    if z.re < T::ZERO {
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

    let two_over_z = two_over_z_safe(modified_z);
    if n_elements_set < n {
        //-----------------------------------------------------------------------
        //     TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW. SET SEQUENCE TO ZERO
        //     ON UNDERFLOW.
        //-----------------------------------------------------------------------
        let max_order = order + T::from_usize(n - 1);
        let DebyeGeometry {
            phi_k: phi,
            zeta1: zet1d,
            zeta2: zet2d,
            ..
        } = DebyeGeometry::compute(modified_z, max_order);
        // let (phi, zet1d, zet2d, _sumd) = ik_uniform_asymp_params(
        //     modified_z,
        //     max_order,
        //     IKType::K,
        //     rotation == RotationDirection::None,
        // );
        let overflow_test = -scaling.scale_zetas(modified_z, max_order, zet1d, zet2d);

        match OverflowState::check(overflow_test.re.abs(), phi, T::ZERO) {
            OverflowState::Over { .. } => return Err(BesselError::Overflow),
            OverflowState::Under { .. } => {
                return if z.re < T::ZERO {
                    Err(BesselError::Overflow)
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
        scale_controlled_recurrence(
            true,
            order,
            modified_z,
            Some(&mut y),
            n_elements_set,
            n,
            s1,
            s2,
            k_overflow_state,
        );
    }
    if rotation == RotationDirection::None {
        return Ok((y, n_zeros));
    }
    //-----------------------------------------------------------------------
    //     ANALYTIC CONTINUATION FOR RE(Z) < 0.0
    //-----------------------------------------------------------------------
    n_zeros = 0;
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
    let mut i_overflow_state = OverflowState::None;
    let mut remaining_n = n;
    for (i, yi) in y.iter_mut().enumerate().rev() {
        remaining_n = i;
        let current_order = order + T::from_usize(i);

        // Reuse from stack if i is 0 or 1, otherwise compute fresh:
        let params = if i < 2 {
            debye_seeds[i].unwrap_or_else(|| DebyeParams::compute(modified_z, current_order))
        } else {
            DebyeParams::compute(modified_z, current_order)
        };

        // Use the I fields:
        let phid = params.phi_i;
        let sumd = params.sum_i;
        //-----------------------------------------------------------------------
        //     LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K
        //     FUNCTION ABOVE
        //-----------------------------------------------------------------------
        // TODO no logic needed! as zunik ccahes the values. Should the other similar functions, too?
        // let (phid, zet1d, zet2d, sumd) =
        //     ik_uniform_asymp_params(modified_z, current_order, IKType::I, false); //, &mut INITD);
        // let sumd = sumd.unwrap();
        let mut s1 = scaling.scale_zetas(modified_z, current_order, params.zeta1, params.zeta2);
        //-----------------------------------------------------------------------
        //     TEST FOR UNDERFLOW AND OVERFLOW
        //-----------------------------------------------------------------------
        let of = OverflowState::check(s1.re, phid, T::ZERO);
        if !found_one_good_entry && !matches!(of, OverflowState::Under { .. }) {
            i_overflow_state = of;
        }
        let mut s2 = match of {
            OverflowState::Over { .. } => {
                return Err(BesselError::Overflow);
            }
            OverflowState::Under { .. } => T::C_ZERO,
            OverflowState::NearOver | OverflowState::NearUnder | OverflowState::None => {
                let st = phid * sumd;
                let mut s2 = T::I * st * rotation_angle;
                s1 = s1.exp() * i_overflow_state.scaling_factor::<T>();
                s2 *= s1;
                if i_overflow_state == OverflowState::NearUnder && will_underflow(s2) {
                    s2 = T::C_ZERO;
                }
                s2
            }
        };
        cy[found_one_good_entry as usize] = s2;
        let c2 = s2;
        s2 *= i_overflow_state.reciprocal_scaling_factor::<T>();
        //-----------------------------------------------------------------------
        //     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N
        //-----------------------------------------------------------------------
        s1 = *yi;
        if scaling == Scaling::Scaled
            && underflow_add_i_k(modified_z, &mut s1, &mut s2, &mut dummy_n_good)
        {
            n_zeros += 1;
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
        let mut reciprocal_scale_factor = i_overflow_state.reciprocal_scaling_factor::<T>();
        let mut absolute_approximation_limit = i_overflow_state.boundary::<T>();
        for (i, yi) in y.iter_mut().enumerate().take(remaining_n).rev() {
            let current_order = order + T::from_usize(i + 1);
            (s1, s2) = (s2, s1 + current_order * (two_over_z * s2));
            let mut unscaled_s2 = s2 * reciprocal_scale_factor;
            let ck = unscaled_s2;

            let mut c1 = *yi;
            if scaling == Scaling::Scaled
                && underflow_add_i_k(modified_z, &mut c1, &mut unscaled_s2, &mut dummy_n_good)
            {
                n_zeros += 1;
            }
            *yi = c1 * cspn + unscaled_s2;
            cspn = -cspn;
            if i_overflow_state == OverflowState::NearOver {
                continue;
            }
            if max_abs_component(unscaled_s2) <= absolute_approximation_limit {
                continue;
            }
            i_overflow_state.increment();
            absolute_approximation_limit = i_overflow_state.boundary::<T>();
            s1 *= reciprocal_scale_factor;
            s2 = ck; // ck is previously calculated s2 * reciprocal_scale_factor
            s1 *= i_overflow_state.scaling_factor::<T>();
            s2 *= i_overflow_state.scaling_factor::<T>();
            reciprocal_scale_factor = i_overflow_state.reciprocal_scaling_factor::<T>();
        }
    }
    Ok((y, n_zeros))
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
pub(crate) fn k_uniform_asymp2<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    rotation: RotationDirection,
    n: usize,
) -> BesselResult<T> {
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
                let s1 = s1.exp() * k_overflow_state.scaling_factor::<T>();
                s2 *= s1;
                if k_overflow_state == OverflowState::NearUnder && will_underflow(s2) {
                    handle_underflow(&mut found_one_good_entry, &mut cs)?
                }
                if zr.im <= T::ZERO {
                    s2 = s2.conj();
                }
                cy[found_one_good_entry as usize] = s2;
                y[i] = s2 * k_overflow_state.reciprocal_scaling_factor::<T>();
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
        match OverflowState::check(s1.re, phid, T::ZERO) {
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
                s1 = s1.exp() * i_overflow_state.scaling_factor::<T>();
                s2 *= s1;
                if i_overflow_state == OverflowState::NearUnder && will_underflow(s2) {
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
        s2 *= i_overflow_state.reciprocal_scaling_factor::<T>();
        //-----------------------------------------------------------------------;
        //     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N;
        //-----------------------------------------------------------------------;
        s1 = *yi;
        if scaling == Scaling::Scaled && underflow_add_i_k(zr, &mut s1, &mut s2, &mut iuf) {
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

        let mut recip_scale_factor = i_overflow_state.reciprocal_scaling_factor::<T>();
        let mut ascle = i_overflow_state.boundary::<T>();
        let mut ck = (order + T::from_usize(remaining_n)) * two_over_z;
        for yi in y.iter_mut().take(remaining_n).rev() {
            (s1, s2) = (s2, s1 + ck * s2);
            ck -= two_over_z;
            let mut c2 = s2 * recip_scale_factor;
            let old_c2 = c2;
            let mut c1 = *yi;
            if scaling == Scaling::Scaled && underflow_add_i_k(zr, &mut c1, &mut c2, &mut iuf) {
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
            );
        }
    }
    Ok((y, n_zeros))
}
