use std::f64::consts::FRAC_PI_2;

use num::{Complex, complex::ComplexFloat};

use crate::{
    amos::{
        PositiveArg,
        asymptotics::consts::{
            AIRY_ASYMP_COEFFS_A, AIRY_ASYMP_COEFFS_B, AIRY_HJ_POLYNOMIAL_COEFFS,
            DEBYE_IK_POLYNOMIAL_COEFFS, IK_NORMALIZATION_FACTORS, TRANSITION_AIRY_A_COEFFS,
            TRANSITION_AIRY_B_COEFFS, TURNING_POINT_ZETA_COEFFS,
        },
    },
    types::BesselFloat,
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
