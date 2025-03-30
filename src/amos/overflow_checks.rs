use std::f64::consts::{FRAC_PI_2, PI};

use num::complex::{Complex64, ComplexFloat};

use crate::amos::{PositiveArg, c_zero, machine::d1mach, utils::will_z_underflow};

use super::{BesselError::*, BesselResult, IKType, MachineConsts, Scaling, c_one, c_zeros};

pub fn zuoik(
    z: Complex64,
    order: f64, //ZR, ZI, FNU,
    kode: Scaling,
    ikflg: IKType,
    n: usize, //YR, YI, NUF,
    mut y: Vec<Complex64>,
    machine_consts: &MachineConsts,
) -> BesselResult {
    // ***BEGIN PROLOGUE  ZUOIK
    // ***REFER TO  ZBESI,ZBESK,ZBESH
    //
    //     ZUOIK COMPUTES THE LEADING TERMS OF THE UNIFORM ASYMPTOTIC
    //     EXPANSIONS FOR THE I AND K FUNCTIONS AND COMPARES THEM
    //     (IN LOGARITHMIC FORM) TO ALIM AND ELIM FOR OVER AND UNDERFLOW
    //     WHERE ALIM < ELIM. if THE MAGNITUDE, BASED ON THE LEADING
    //     EXPONENTIAL, IS LESS THAN ALIM OR GREATER THAN -ALIM, THEN
    //     THE RESULT IS ON SCALE. if NOT, THEN A REFINED TEST USING OTHER
    //     MULTIPLIERS (IN LOGARITHMIC FORM) IS MADE BASED ON ELIM. HERE
    //     EXP(-ELIM)=SMALLEST MACHINE NUMBER*1.0E+3 AND EXP(-ALIM)=
    //     EXP(-ELIM)/TOL
    //
    //     IKFLG=1 MEANS THE I SEQUENCE IS TESTED
    //          =2 MEANS THE K SEQUENCE IS TESTED
    //     NUF = 0 MEANS THE LAST MEMBER OF THE SEQUENCE IS ON SCALE
    //         =-1 MEANS AN OVERFLOW WOULD OCCUR
    //     IKFLG=1 AND NUF > 0 MEANS THE LAST NUF Y VALUES WERE SET TO ZERO
    //             THE FIRST N-NUF VALUES MUST BE SET BY ANOTHER ROUTINE
    //     IKFLG=2 AND NUF == N MEANS ALL Y VALUES WERE SET TO ZERO
    //     IKFLG=2 AND 0 < NUF < N NOT CONSIDERED. Y MUST BE SET BY
    //             ANOTHER ROUTINE
    //
    // ***ROUTINES CALLED  ZUunderflowCHK,ZUNHJ,ZUNIK,d1mach,ZABS,ZLOG
    // ***END PROLOGUE  ZUOIK
    const AIC: f64 = 1.265512123484645396e+00;
    let mut nuf = 0;
    let mut nn = n;
    let zr = if z.re < 0.0 { -z } else { z };
    let zb = zr;
    let ax = z.re.abs() * 1.7321;
    let ay = z.im.abs();
    let iform = if ay > ax { 2 } else { 1 };
    let mut gnu = order.max(1.0);
    if ikflg == IKType::K {
        let fnn = nn as f64;
        let gnn = order + fnn - 1.0;
        gnu = gnn.max(fnn);
    }
    //-----------------------------------------------------------------------;
    //     ONLY THE MAGNITUDE OF ARG AND PHI ARE NEEDED ALONG WITH THE;
    //     REAL PARTS OF ZETA1, ZETA2 AND ZB. NO ATTEMPT IS MADE TO GET;
    //     THE SIGN OF THE IMAGINARY PART CORRECT.;
    //-----------------------------------------------------------------------;
    let mut zn = None;
    let (mut cz, phi, arg, aarg) = if iform != 2 {
        let mut init = 0;
        let (phi, zeta1, zeta2, _) = zunik(zr, gnu, ikflg, true, machine_consts, &mut init);
        (-zeta1 + zeta2, phi, c_zero(), 0.0)
    } else {
        let mut zn_ = Complex64::new(zr.im, -zr.re);
        if z.im <= 0.0 {
            zn_.re = -zn_.re;
        }
        let (phi, arg, zeta1, zeta2, _, _) = zunhj(zn_, gnu, true, machine_consts.tol);
        zn = Some(zn_);
        let cz = -zeta1 + zeta2;
        let aarg = arg.abs();
        (cz, phi, arg, aarg)
    };
    if kode == Scaling::Scaled {
        cz -= zb;
    }
    if ikflg == IKType::K {
        cz = -cz;
    }
    let aphi = phi.abs();
    let mut rcz = cz.re;
    //-----------------------------------------------------------------------;
    //     OVERFLOW TEST;
    //-----------------------------------------------------------------------;
    if rcz > machine_consts.elim {
        return Err(Overflow);
    }
    if rcz >= machine_consts.alim {
        rcz += aphi.ln();
        if iform == 2 {
            rcz = rcz - 0.25 * aarg.ln() - AIC
        };
        if rcz > machine_consts.elim {
            return Err(Overflow);
        }
    } else {
        //-----------------------------------------------------------------------;
        //     UNDERFLOW TEST;
        //-----------------------------------------------------------------------;
        if rcz < (-machine_consts.elim) {
            y[0..nn].iter_mut().for_each(|v| *v = c_zero());
            return Ok((y, nn));
        }
        if rcz <= (-machine_consts.alim) {
            rcz += aphi.ln();
            if iform == 2 {
                rcz = rcz - 0.25 * aarg.ln() - AIC
            };
            if !(rcz > (-machine_consts.elim)) {
                return Ok((y, nn));
            }
            cz += phi.ln();
            if iform != 1 {
                cz -= 0.25 * arg.ln() + AIC
            }
            let ax = rcz.exp() / machine_consts.tol;
            let ay = cz.im;
            cz = ax * Complex64::cis(ay);
            if will_z_underflow(cz, machine_consts.ascle, machine_consts.tol) {
                y[0..nn].iter_mut().for_each(|v| *v = c_zero());
                return Ok((y, nn));
            }
        }
    }
    if ikflg == IKType::K || n == 1 {
        return Ok((y, nuf));
    }
    //-----------------------------------------------------------------------;
    //     SET UNDERFLOWS ON I SEQUENCE;
    //-----------------------------------------------------------------------;
    let mut go_to_180 = false;
    let mut skip_to_190;
    'outer: loop {
        'l140: loop {
            skip_to_190 = false;
            if !go_to_180 {
                gnu = order + ((nn - 1) as f64);
                let (phi, mut cz, aarg) = if iform != 2 {
                    let mut init = 0;
                    let (phi, zeta1, zeta2, _) =
                        zunik(zr, gnu, ikflg, true, machine_consts, &mut init);
                    cz = -zeta1 + zeta2;
                    (phi, cz, 0.0)
                } else {
                    let (phi, arg, zeta1, zeta2, _, _) =
                        zunhj(zn.unwrap(), gnu, true, machine_consts.tol);
                    cz = -zeta1 + zeta2;
                    let aarg = arg.abs();
                    (phi, cz, aarg)
                };
                if kode == Scaling::Scaled {
                    cz -= zb;
                }
                let aphi = phi.abs();
                rcz = cz.re;
                if !(rcz < (-machine_consts.elim)) {
                    if rcz > (-machine_consts.alim) {
                        return Ok((y, nuf));
                    };
                    rcz += aphi.ln();
                    if iform == 2 {
                        rcz = rcz - 0.25 * aarg.ln() - AIC;
                    }
                    if rcz > (-machine_consts.elim) {
                        skip_to_190 = true
                    }
                }
            }
            go_to_180 = false;
            if !skip_to_190 {
                y[nn - 1] = c_zero();
                nn = nn - 1;
                nuf = nuf + 1;
                if nn == 0 {
                    return Ok((y, nuf));
                }
            } else {
                break 'l140;
            }
        }
        let ascle = 1.0e+3 * d1mach(1) / machine_consts.tol;
        cz += phi.ln();
        if iform != 1 {
            cz -= 0.25 * arg.ln() + AIC
        }
        let ax = rcz.exp() / machine_consts.tol;
        let ay = cz.im;
        cz = ax * Complex64::cis(ay);
        if will_z_underflow(cz, ascle, machine_consts.tol) {
            go_to_180 = true;
        } else {
            break 'outer;
        }
    }
    return Ok((y, nuf));
}

fn zunik(
    zr: Complex64,
    order: f64,
    ikflg: IKType,
    only_phi_zeta: bool,
    machine_consts: &MachineConsts,
    init: &mut usize,
) -> (Complex64, Complex64, Complex64, Option<Complex64>) {
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
    let mut working = c_zeros(16);

    if *init != 0 {
        todo!()
    }
    // if (INIT != 0) GO TO 40;
    //-----------------------------------------------------------------------;
    //     INITIALIZE ALL VARIABLES;
    //-----------------------------------------------------------------------;
    let rfn = 1.0 / order;
    //-----------------------------------------------------------------------;
    //     OVERFLOW TEST (ZR/FNU TOO SMALL);
    //-----------------------------------------------------------------------;
    let test = d1mach(1) * 1.0e+3;
    let ac = order * test;
    if !(zr.re.abs() > ac || zr.im.abs() > ac) {
        let zeta1 = Complex64::new(2.0 * test.ln().abs() + order, 0.0);
        let zeta2 = Complex64::new(order, 0.0);
        let phi = c_one();
        return (zeta1, zeta2, phi, None);
    }
    let t = zr * rfn;
    let s = c_one() + t * t;
    let s_root = s.sqrt();
    let zn = (c_one() + s_root) / t;
    let zeta1 = order * zn.ln();
    let zeta2 = order * s_root;
    let t = c_one() / s_root;
    let sr = t * rfn;
    working[15] = sr.sqrt();
    let mut phi = working[15] * CON[ikflg.index() - 1_usize];
    if only_phi_zeta {
        return (phi, zeta1, zeta2, None);
    };
    let t2 = c_one() / s;
    working[0] = c_one();
    let mut crfn = c_one();
    let mut ac = 1.0;
    let mut l = 1;
    let mut k = 0;
    'l20: for k_ in 1..15 {
        k = k_;
        let mut s = c_zero();
        '_l10: for _ in 0..k_ {
            l += 1;
            s = s * t2 + C_ZUNIK[l];
        }
        crfn *= sr;
        working[k_] = crfn * s;
        ac *= rfn;
        let test = working[k_].re.abs() + working[k_].im.abs();
        if ac < machine_consts.tol && test < machine_consts.tol {
            break 'l20;
        } //GO TO 30;
    }
    *init = k;
    //    40 CONTINUE;
    //-----------------------------------------------------------------------;
    // FINISH INIT. Use OnceCell here later?
    //-----------------------------------------------------------------------;
    match ikflg {
        IKType::I => {
            //GO TO 60;
            //-----------------------------------------------------------------------;
            //     COMPUTE SUM FOR THE I FUNCTION;
            //-----------------------------------------------------------------------;
            // s=c_zero();
            // SR = ZEROR;
            // SI = ZEROI;
            // DO 50 I=1,INIT;
            // for i in 0..INIT{
            //       s += working[i]
            //   SR = SR + CWRKR(I);
            //   SI = SI + CWRKI(I);
            // }
            //    50 CONTINUE;
            // SUMR = SR;
            // SUMI = SI;
            let sum = working[..*init].iter().sum();
            phi = working[15] * CON[0];
            // PHIR = CWRKR(16)*CON(1);
            // PHII = CWRKI(16)*CON(1);
            // RETURN;
            (phi, zeta1, zeta2, Some(sum))
        }
        IKType::K => {
            //    60 CONTINUE;
            //-----------------------------------------------------------------------;
            //     COMPUTE SUM FOR THE K FUNCTION;
            //-----------------------------------------------------------------------;
            let mut tr = 1.0;
            let sum = working[..*init]
                .iter()
                .map(|v| {
                    let output = tr * v;
                    tr = -tr;
                    output
                })
                .sum();

            // SR = ZEROR;
            // SI = ZEROI;
            // TR = CONER;
            // DO 70 I=1,INIT;
            //   SR = SR + TR*CWRKR(I);
            //   SI = SI + TR*CWRKI(I);
            //   TR = -TR;
            //    70 CONTINUE;
            //       SUMR = SR;
            //       SUMI = SI;
            // PHIR = CWRKR(16)*CON(2);
            // PHII = CWRKI(16)*CON(2);
            phi = working[15] * CON[1];
            (phi, zeta1, zeta2, Some(sum))

            // RETURN;
            // END;
        }
    }
}

fn zunhj(
    z: Complex64,
    order: f64,
    only_phi_zeta: bool,
    tol: f64,
    //returns: phi, arg, zeta1, zeta2, asum, bsum
) -> (
    Complex64,
    Complex64,
    Complex64,
    Complex64,
    Option<Complex64>,
    Option<Complex64>,
) {
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

    const EX1: f64 = 3.33333333333333333e-01;
    const EX2: f64 = 6.66666666666666667e-01;
    const THREE_PI_BY_2: f64 = 4.71238898038468986e+00;
    let rfnu = 1.0 / order;
    //-----------------------------------------------------------------------
    //     OVERFLOW TEST (Z/FNU TOO SMALL)
    //-----------------------------------------------------------------------
    let test = d1mach(1) * 1.0e+3;
    let ac = order * test;
    if !((z.re).abs() > ac || (z.im).abs() > ac) {
        let zeta1 = Complex64::new(2.0 * test.ln().abs() + order, 0.0);
        let zeta2 = Complex64::new(order, 0.0);
        let phi = c_one();
        let arg = c_one();
        return (phi, arg, zeta1, zeta2, None, None);
    }
    let zb = z * rfnu;
    let rfnu2 = rfnu * rfnu;
    //-----------------------------------------------------------------------;
    //     COMPUTE IN THE FOURTH QUADRANT;
    //-----------------------------------------------------------------------;
    let fn13 = order.powf(EX1);
    let fn23 = fn13 * fn13;
    let rfn13 = 1.0 / fn13;
    let w2 = c_one() - zb * zb;
    let aw2 = w2.abs();

    if aw2 <= 0.25 {
        //-----------------------------------------------------------------------;
        //     POWER SERIES FOR CABS(W2) <= 0.25;
        //-----------------------------------------------------------------------;
        let mut k = 0;
        let mut p = c_zeros(30);
        let mut ap = vec![0.0; 30];
        p[0] = c_one();
        let mut suma = Complex64::new(GAMA[0], 0.0);
        ap[0] = 1.0;
        if aw2 >= tol {
            for k_ in 1..30 {
                k = k_;
                p[k_] = p[k_ - 1] * w2;
                suma += p[k_] * GAMA[k_];
                ap[k_] = ap[k_ - 1] * aw2;
                if ap[k_] < tol {
                    break;
                }
            }
        }
        let kmax = k;
        let zeta = w2 * suma;
        let arg = zeta * fn23;
        let mut za = suma.sqrt();
        let zeta2 = w2.sqrt() * order;
        let zeta1 = (c_one() + EX2 * zeta * za) * zeta2;
        za *= 2.0;
        let phi = za.sqrt() * rfn13;
        if only_phi_zeta {
            return (phi, arg, zeta1, zeta2, None, None);
        }
        //-----------------------------------------------------------------------;
        //     SUM SERIES FOR ASUM AND BSUM;
        //-----------------------------------------------------------------------;
        let sumb: Complex64 = p[..kmax].iter().zip(BETA).map(|(p, b)| p * b).sum();
        let mut asum = c_zero();
        let mut bsum = sumb;
        let mut l1 = 0;
        let mut l2 = 30;
        let btol = tol * (bsum.re.abs() + bsum.im.abs());
        let mut atol = tol;
        let mut pp = 1.0;
        let mut ias = false;
        let mut ibs = false;
        if !(rfnu2 < tol) {
            //GO TO 110;
            // DO 100 IS=2,7;
            for _is in 1..7 {
                atol /= rfnu2;
                pp *= rfnu2;
                if !ias {
                    //GO TO 60;
                    let mut suma = c_zero();
                    //   SUMAR = ZEROR;
                    //   SUMAI = ZEROI;
                    //   DO 40 K=1,KMAX;
                    for k in 0..kmax {
                        //     let M = L1 + K;
                        suma += p[k] * ALFA[l1 + k];
                        //     SUMAR = SUMAR + PR(K)*ALFA(M);
                        //     SUMAI = SUMAI + PI(K)*ALFA(M);
                        if ap[k] < atol {
                            break;
                        } //GO TO 50;
                    }
                    //    40   CONTINUE;
                    //    50   CONTINUE;
                    asum += suma * pp;
                    //   ASUMR = ASUMR + SUMAR*PP;
                    //   ASUMI = ASUMI + SUMAI*PP;
                    if pp < tol {
                        ias = true
                    };
                }
                //    60   CONTINUE;
                if !ibs {
                    //GO TO 90;
                    let mut sumb = c_zero();
                    //   SUMBR = ZEROR;
                    //   SUMBI = ZEROI;
                    //   DO 70 K=1,KMAX;
                    for k in 0..kmax {
                        //     let M = L2 + K;
                        sumb += p[k] * BETA[l2 + k];
                        //     SUMBR = SUMBR + PR(K)*BETA(M);
                        //     SUMBI = SUMBI + PI(K)*BETA(M);
                        if ap[k] < atol {
                            break;
                        } //GO TO 80;
                    }
                    //    70   CONTINUE;
                    //    80   CONTINUE;
                    bsum += sumb * pp;
                    //   BSUMR = BSUMR + SUMBR*PP;
                    //   BSUMI = BSUMI + SUMBI*PP;
                    if pp < btol {
                        ibs = true;
                    }
                }
                //    90   CONTINUE;
                if ias && ibs {
                    break;
                } //GO TO 110;
                l1 = l1 + 30;
                l2 = l2 + 30;
            }
            //   100 CONTINUE;
        }
        //   110 CONTINUE;
        asum += 1.0;
        // ASUMR = ASUMR + CONER;
        pp = rfnu * rfn13;
        bsum *= pp;
        // BSUMR = BSUMR*PP;
        // BSUMI = BSUMI*PP;
        //   120 CONTINUE;
        //       RETURN;
        return (phi, arg, zeta1, zeta2, Some(bsum), Some(bsum));
    } else {
        //-----------------------------------------------------------------------;
        //     CABS(W2) > 0.25;
        //-----------------------------------------------------------------------;
        let mut w = w2.sqrt();
        if w.re < 0.0 {
            w.re = 0.0
        };
        if w.im < 0.0 {
            w.im = 0.0
        };

        let za = (c_one() + w) / zb;
        let mut zc = za.ln();
        if zc.im < 0.0 {
            zc.im = 0.0;
        }

        if zc.im > FRAC_PI_2 {
            zc.im = FRAC_PI_2;
        }
        if zc.re < 0.0 {
            zc.re = 0.0
        };
        let zth = (zc - w) * 1.5;
        let zeta1 = zc * order;
        let zeta2 = w * order;
        let azth = zth.abs();
        let mut ang = zth.parg();
        ang = ang.clamp(0.0, THREE_PI_BY_2);
        let mut pp = azth.powf(EX2);
        ang *= EX2;
        let mut zeta = pp * Complex64::cis(ang);
        if zeta.im < 0.0 {
            zeta.im = 0.0
        };
        let arg = zeta * fn23;
        let rtzt = zth / zeta;
        let za = rtzt / w;
        let tazr = za + za;
        let phi = tazr.sqrt() * rfn13;
        if only_phi_zeta {
            return (phi, arg, zeta1, zeta2, None, None);
        }

        let raw = 1.0 / aw2.sqrt();
        // STR = WR*RAW;
        // STI = -WI*RAW;
        let tfn = w.conj() * raw * raw * rfnu;
        // TFNR = STR*RFNU*RAW;
        // TFNI = STI*RFNU*RAW;
        let razth = 1.0 / azth;
        // STR = ZTHR*RAZTH;
        // STI = -ZTHI*RAZTH;
        let rzth = zth.conj() * razth * razth * rfnu;
        // RZTHR = STR*RAZTH*RFNU;
        // RZTHI = STI*RAZTH*RFNU;
        let zc = razth * AR[1];
        // ZCR = RZTHR*AR(2);
        // ZCI = RZTHI*AR(2);
        let raw2 = 1.0 / aw2;
        // STR = W2R*RAW2;
        // STI = -W2I*RAW2;
        let t2 = w2.conj() * raw2 * raw2;
        // T2R = STR*RAW2;
        // T2I = STI*RAW2;
        // let st = t2 * C_ZUNHJ[1] + C_ZUNHJ[2];
        // STR = T2R*C(2) + C(3);
        // STI = T2I*C(2);
        let mut up = c_zeros(14);
        up[1] = (t2 * C_ZUNHJ[1] + C_ZUNHJ[2]) * tfn;
        // UPR(2) = STR*TFNR - STI*TFNI;
        // UPI(2) = STR*TFNI + STI*TFNR;
        let mut bsum = up[1] + zc;
        // BSUMR = UPR(2) + ZCR;
        // BSUMI = UPI(2) + ZCI;
        let mut asum = c_zero();
        // ASUMR = ZEROR;
        // ASUMI = ZEROI;
        if !(rfnu < tol) {
            //GO TO 220;
            let mut przth = rzth;
            // PRZTHR = RZTHR;
            // PRZTHI = RZTHI;
            let mut ptfn = tfn;
            // PTFNR = TFNR;
            // PTFNI = TFNI;
            up[0] = c_one();
            // UPR(1) = CONER;
            // UPI(1) = CONEI;
            pp = 1.0;
            let btol = tol * (bsum.re.abs() + bsum.im.abs());
            let mut ks = 0;
            let mut kp1 = 2;
            let mut l = 2; //3;
            let mut ias = false;
            let mut ibs = false;
            // DO 210 LR=2,12,2;
            for lr in (2..=12).step_by(2) {
                let lrp1 = lr + 1;
                //-----------------------------------------------------------------------;
                //     COMPUTE TWO ADDITIONAL CR, DR, AND UP FOR TWO MORE TERMS IN;
                //     NEXT SUMA AND SUMB;
                //-----------------------------------------------------------------------;
                //   DO 160 K=LR,LRP1;
                let mut cr = c_zeros(14);
                let mut dr = c_zeros(14);
                for _k in lr..=lrp1 {
                    ks = ks + 1;
                    kp1 = kp1 + 1;
                    l = l + 1;
                    let mut za = Complex64::new(C_ZUNHJ[l - 1], 0.0);
                    //     ZAR = C(L);
                    //     ZAI = ZEROI;
                    //     DO 150 J=2,KP1;
                    for _j in 2..=kp1 {
                        l = l + 1;
                        // STR = ZAR*T2R - T2I*ZAI + C(L);
                        // ZAI = ZAR*T2I + ZAI*T2R;
                        // ZAR = STR;
                        za = za * t2 + C_ZUNHJ[l - 1];
                    }
                    //   150     CONTINUE;
                    ptfn *= tfn;
                    //     STR = PTFNR*TFNR - PTFNI*TFNI;
                    //     PTFNI = PTFNR*TFNI + PTFNI*TFNR;
                    //     PTFNR = STR;
                    up[kp1 - 1] = ptfn * za;
                    //     UPR(KP1) = PTFNR*ZAR - PTFNI*ZAI;
                    //     UPI(KP1) = PTFNI*ZAR + PTFNR*ZAI;
                    cr[ks - 1] = przth * BR[ks];
                    //     CRR(KS) = PRZTHR*BR(KS+1);
                    //     CRI(KS) = PRZTHI*BR(KS+1);
                    przth *= rzth;
                    //     STR = PRZTHR*RZTHR - PRZTHI*RZTHI;
                    //     PRZTHI = PRZTHR*RZTHI + PRZTHI*RZTHR;
                    //     PRZTHR = STR;
                    dr[ks - 1] = przth * AR[ks + 1];
                    //     DRR(KS) = PRZTHR*AR(KS+2);
                    //     DRI(KS) = PRZTHI*AR(KS+2);
                }
                //   160   CONTINUE;
                pp = pp * rfnu2;
                if !ias {
                    //GO TO 180;
                    let mut suma = up[lrp1 - 1];
                    //   SUMAR = UPR(LRP1);
                    //   SUMAI = UPI(LRP1);
                    let mut ju = lrp1;
                    //   DO 170 JR=1,LR;
                    for jr in 0..lr {
                        ju = ju - 1;
                        suma += cr[jr] * up[ju - 1];
                        //     SUMAR = SUMAR + CRR(JR)*UPR(JU) - CRI(JR)*UPI(JU);
                        //     SUMAI = SUMAI + CRR(JR)*UPI(JU) + CRI(JR)*UPR(JU);
                    }
                    //   170   CONTINUE;
                    asum += suma;
                    //   ASUMR = ASUMR + SUMAR;
                    //   ASUMI = ASUMI + SUMAI;
                    let test = suma.re.abs() + suma.im.abs();
                    if pp < tol && test < tol {
                        ias = true
                    };
                }
                //   180   CONTINUE;
                if !ibs {
                    //GO TO 200;
                    let mut sumb = up[lr - 2] + up[lrp1 - 1] * zc;
                    //   SUMBR = UPR(LR+2) + UPR(LRP1)*ZCR - UPI(LRP1)*ZCI;
                    //   SUMBI = UPI(LR+2) + UPR(LRP1)*ZCI + UPI(LRP1)*ZCR;
                    let mut ju = lrp1;
                    for jr in 0..lr {
                        //   DO 190 JR=1,LR;
                        ju = ju - 1;
                        sumb += dr[jr] * up[ju - 1];
                        //     SUMBR = SUMBR + DRR(JR)*UPR(JU) - DRI(JR)*UPI(JU);
                        //     SUMBI = SUMBI + DRR(JR)*UPI(JU) + DRI(JR)*UPR(JU);
                    }
                    //   190   CONTINUE;
                    bsum += sumb;
                    //   BSUMR = BSUMR + SUMBR;
                    //   BSUMI = BSUMI + SUMBI;
                    let test = sumb.re.abs() + sumb.im.abs();
                    if pp < btol && test < btol {
                        ibs = true
                    };
                    //   200   CONTINUE;
                }
                if ias && ibs {
                    break;
                }
                //   210 CONTINUE;
            }
        }
        //   220 CONTINUE;
        asum += c_one();
        // ASUMR = ASUMR + CONER;
        bsum = (-bsum * rfn13) / rtzt;
        // STR = -BSUMR*RFN13;
        // STI = -BSUMI*RFN13;
        // CALL ZDIV(STR, STI, RTZTR, RTZTI, BSUMR, BSUMI);
        // GO TO 120;
        return (phi, arg, zeta1, zeta2, Some(bsum), Some(bsum));

        // END;
    }
}

const CON: [f64; 2] = [3.98942280401432678e-01, 1.25331413731550025e+00];
const C_ZUNIK: [f64; 120] = [
    1.00000000000000000e+00,
    -2.08333333333333333e-01,
    1.25000000000000000e-01,
    3.34201388888888889e-01,
    -4.01041666666666667e-01,
    7.03125000000000000e-02,
    -1.02581259645061728e+00,
    1.84646267361111111e+00,
    -8.91210937500000000e-01,
    7.32421875000000000e-02,
    4.66958442342624743e+00,
    -1.12070026162229938e+01,
    8.78912353515625000e+00,
    -2.36408691406250000e+00,
    1.12152099609375000e-01,
    -2.82120725582002449e+01,
    8.46362176746007346e+01,
    -9.18182415432400174e+01,
    4.25349987453884549e+01,
    -7.36879435947963170e+00,
    2.27108001708984375e-01,
    2.12570130039217123e+02,
    -7.65252468141181642e+02,
    1.05999045252799988e+03,
    -6.99579627376132541e+02,
    2.18190511744211590e+02,
    -2.64914304869515555e+01,
    5.72501420974731445e-01,
    -1.91945766231840700e+03,
    8.06172218173730938e+03,
    -1.35865500064341374e+04,
    1.16553933368645332e+04,
    -5.30564697861340311e+03,
    1.20090291321635246e+03,
    -1.08090919788394656e+02,
    1.72772750258445740e+00,
    2.02042913309661486e+04,
    -9.69805983886375135e+04,
    1.92547001232531532e+05,
    -2.03400177280415534e+05,
    1.22200464983017460e+05,
    -4.11926549688975513e+04,
    7.10951430248936372e+03,
    -4.93915304773088012e+02,
    6.07404200127348304e+00,
    -2.42919187900551333e+05,
    1.31176361466297720e+06,
    -2.99801591853810675e+06,
    3.76327129765640400e+06,
    -2.81356322658653411e+06,
    1.26836527332162478e+06,
    -3.31645172484563578e+05,
    4.52187689813627263e+04,
    -2.49983048181120962e+03,
    2.43805296995560639e+01,
    3.28446985307203782e+06,
    -1.97068191184322269e+07,
    5.09526024926646422e+07,
    -7.41051482115326577e+07,
    6.63445122747290267e+07,
    -3.75671766607633513e+07,
    1.32887671664218183e+07,
    -2.78561812808645469e+06,
    3.08186404612662398e+05,
    -1.38860897537170405e+04,
    1.10017140269246738e+02,
    -4.93292536645099620e+07,
    3.25573074185765749e+08,
    -9.39462359681578403e+08,
    1.55359689957058006e+09,
    -1.62108055210833708e+09,
    1.10684281682301447e+09,
    -4.95889784275030309e+08,
    1.42062907797533095e+08,
    -2.44740627257387285e+07,
    2.24376817792244943e+06,
    -8.40054336030240853e+04,
    5.51335896122020586e+02,
    8.14789096118312115e+08,
    -5.86648149205184723e+09,
    1.86882075092958249e+10,
    -3.46320433881587779e+10,
    4.12801855797539740e+10,
    -3.30265997498007231e+10,
    1.79542137311556001e+10,
    -6.56329379261928433e+09,
    1.55927986487925751e+09,
    -2.25105661889415278e+08,
    1.73951075539781645e+07,
    -5.49842327572288687e+05,
    3.03809051092238427e+03,
    -1.46792612476956167e+10,
    1.14498237732025810e+11,
    -3.99096175224466498e+11,
    8.19218669548577329e+11,
    -1.09837515608122331e+12,
    1.00815810686538209e+12,
    -6.45364869245376503e+11,
    2.87900649906150589e+11,
    -8.78670721780232657e+10,
    1.76347306068349694e+10,
    -2.16716498322379509e+09,
    1.43157876718888981e+08,
    -3.87183344257261262e+06,
    1.82577554742931747e+04,
    2.86464035717679043e+11,
    -2.40629790002850396e+12,
    9.10934118523989896e+12,
    -2.05168994109344374e+13,
    3.05651255199353206e+13,
    -3.16670885847851584e+13,
    2.33483640445818409e+13,
    -1.23204913055982872e+13,
    4.61272578084913197e+12,
    -1.19655288019618160e+12,
    2.05914503232410016e+11,
    -2.18229277575292237e+10,
    1.24700929351271032e+09,
    -2.91883881222208134e+07,
    1.18838426256783253e+05,
];

const AR: [f64; 14] = [
    1.00000000000000000e+00,
    1.04166666666666667e-01,
    8.35503472222222222e-02,
    1.28226574556327160e-01,
    2.91849026464140464e-01,
    8.81627267443757652e-01,
    3.32140828186276754e+00,
    1.49957629868625547e+01,
    7.89230130115865181e+01,
    4.74451538868264323e+02,
    3.20749009089066193e+03,
    2.40865496408740049e+04,
    1.98923119169509794e+05,
    1.79190200777534383e+06,
];

const BR: [f64; 14] = [
    1.00000000000000000e+00,
    -1.45833333333333333e-01,
    -9.87413194444444444e-02,
    -1.43312053915895062e-01,
    -3.17227202678413548e-01,
    -9.42429147957120249e-01,
    -3.51120304082635426e+00,
    -1.57272636203680451e+01,
    -8.22814390971859444e+01,
    -4.92355370523670524e+02,
    -3.31621856854797251e+03,
    -2.48276742452085896e+04,
    -2.04526587315129788e+05,
    -1.83844491706820990e+06,
];
const C_ZUNHJ: [f64; 105] = [
    1.00000000000000000e+00,
    -2.08333333333333333e-01,
    1.25000000000000000e-01,
    3.34201388888888889e-01,
    -4.01041666666666667e-01,
    7.03125000000000000e-02,
    -1.02581259645061728e+00,
    1.84646267361111111e+00,
    -8.91210937500000000e-01,
    7.32421875000000000e-02,
    4.66958442342624743e+00,
    -1.12070026162229938e+01,
    8.78912353515625000e+00,
    -2.36408691406250000e+00,
    1.12152099609375000e-01,
    -2.82120725582002449e+01,
    8.46362176746007346e+01,
    -9.18182415432400174e+01,
    4.25349987453884549e+01,
    -7.36879435947963170e+00,
    2.27108001708984375e-01,
    2.12570130039217123e+02,
    -7.65252468141181642e+02,
    1.05999045252799988e+03,
    -6.99579627376132541e+02,
    2.18190511744211590e+02,
    -2.64914304869515555e+01,
    5.72501420974731445e-01,
    -1.91945766231840700e+03,
    8.06172218173730938e+03,
    -1.35865500064341374e+04,
    1.16553933368645332e+04,
    -5.30564697861340311e+03,
    1.20090291321635246e+03,
    -1.08090919788394656e+02,
    1.72772750258445740e+00,
    2.02042913309661486e+04,
    -9.69805983886375135e+04,
    1.92547001232531532e+05,
    -2.03400177280415534e+05,
    1.22200464983017460e+05,
    -4.11926549688975513e+04,
    7.10951430248936372e+03,
    -4.93915304773088012e+02,
    6.07404200127348304e+00,
    -2.42919187900551333e+05,
    1.31176361466297720e+06,
    -2.99801591853810675e+06,
    3.76327129765640400e+06,
    -2.81356322658653411e+06,
    1.26836527332162478e+06,
    -3.31645172484563578e+05,
    4.52187689813627263e+04,
    -2.49983048181120962e+03,
    2.43805296995560639e+01,
    3.28446985307203782e+06,
    -1.97068191184322269e+07,
    5.09526024926646422e+07,
    -7.41051482115326577e+07,
    6.63445122747290267e+07,
    -3.75671766607633513e+07,
    1.32887671664218183e+07,
    -2.78561812808645469e+06,
    3.08186404612662398e+05,
    -1.38860897537170405e+04,
    1.10017140269246738e+02,
    -4.93292536645099620e+07,
    3.25573074185765749e+08,
    -9.39462359681578403e+08,
    1.55359689957058006e+09,
    -1.62108055210833708e+09,
    1.10684281682301447e+09,
    -4.95889784275030309e+08,
    1.42062907797533095e+08,
    -2.44740627257387285e+07,
    2.24376817792244943e+06,
    -8.40054336030240853e+04,
    5.51335896122020586e+02,
    8.14789096118312115e+08,
    -5.86648149205184723e+09,
    1.86882075092958249e+10,
    -3.46320433881587779e+10,
    4.12801855797539740e+10,
    -3.30265997498007231e+10,
    1.79542137311556001e+10,
    -6.56329379261928433e+09,
    1.55927986487925751e+09,
    -2.25105661889415278e+08,
    1.73951075539781645e+07,
    -5.49842327572288687e+05,
    3.03809051092238427e+03,
    -1.46792612476956167e+10,
    1.14498237732025810e+11,
    -3.99096175224466498e+11,
    8.19218669548577329e+11,
    -1.09837515608122331e+12,
    1.00815810686538209e+12,
    -6.45364869245376503e+11,
    2.87900649906150589e+11,
    -8.78670721780232657e+10,
    1.76347306068349694e+10,
    -2.16716498322379509e+09,
    1.43157876718888981e+08,
    -3.87183344257261262e+06,
    1.82577554742931747e+04,
];
const ALFA: [f64; 180] = [
    -4.44444444444444444e-03,
    -9.22077922077922078e-04,
    -8.84892884892884893e-05,
    1.65927687832449737e-04,
    2.46691372741792910e-04,
    2.65995589346254780e-04,
    2.61824297061500945e-04,
    2.48730437344655609e-04,
    2.32721040083232098e-04,
    2.16362485712365082e-04,
    2.00738858762752355e-04,
    1.86267636637545172e-04,
    1.73060775917876493e-04,
    1.61091705929015752e-04,
    1.50274774160908134e-04,
    1.40503497391269794e-04,
    1.31668816545922806e-04,
    1.23667445598253261e-04,
    1.16405271474737902e-04,
    1.09798298372713369e-04,
    1.03772410422992823e-04,
    9.82626078369363448e-05,
    9.32120517249503256e-05,
    8.85710852478711718e-05,
    8.42963105715700223e-05,
    8.03497548407791151e-05,
    7.66981345359207388e-05,
    7.33122157481777809e-05,
    7.01662625163141333e-05,
    6.72375633790160292e-05,
    6.93735541354588974e-04,
    2.32241745182921654e-04,
    -1.41986273556691197e-05,
    -1.16444931672048640e-04,
    -1.50803558053048762e-04,
    -1.55121924918096223e-04,
    -1.46809756646465549e-04,
    -1.33815503867491367e-04,
    -1.19744975684254051e-04,
    -1.06184319207974020e-04,
    -9.37699549891194492e-05,
    -8.26923045588193274e-05,
    -7.29374348155221211e-05,
    -6.44042357721016283e-05,
    -5.69611566009369048e-05,
    -5.04731044303561628e-05,
    -4.48134868008882786e-05,
    -3.98688727717598864e-05,
    -3.55400532972042498e-05,
    -3.17414256609022480e-05,
    -2.83996793904174811e-05,
    -2.54522720634870566e-05,
    -2.28459297164724555e-05,
    -2.05352753106480604e-05,
    -1.84816217627666085e-05,
    -1.66519330021393806e-05,
    -1.50179412980119482e-05,
    -1.35554031379040526e-05,
    -1.22434746473858131e-05,
    -1.10641884811308169e-05,
    -3.54211971457743841e-04,
    -1.56161263945159416e-04,
    3.04465503594936410e-05,
    1.30198655773242693e-04,
    1.67471106699712269e-04,
    1.70222587683592569e-04,
    1.56501427608594704e-04,
    1.36339170977445120e-04,
    1.14886692029825128e-04,
    9.45869093034688111e-05,
    7.64498419250898258e-05,
    6.07570334965197354e-05,
    4.74394299290508799e-05,
    3.62757512005344297e-05,
    2.69939714979224901e-05,
    1.93210938247939253e-05,
    1.30056674793963203e-05,
    7.82620866744496661e-06,
    3.59257485819351583e-06,
    1.44040049814251817e-07,
    -2.65396769697939116e-06,
    -4.91346867098485910e-06,
    -6.72739296091248287e-06,
    -8.17269379678657923e-06,
    -9.31304715093561232e-06,
    -1.02011418798016441e-05,
    -1.08805962510592880e-05,
    -1.13875481509603555e-05,
    -1.17519675674556414e-05,
    -1.19987364870944141e-05,
    3.78194199201772914e-04,
    2.02471952761816167e-04,
    -6.37938506318862408e-05,
    -2.38598230603005903e-04,
    -3.10916256027361568e-04,
    -3.13680115247576316e-04,
    -2.78950273791323387e-04,
    -2.28564082619141374e-04,
    -1.75245280340846749e-04,
    -1.25544063060690348e-04,
    -8.22982872820208365e-05,
    -4.62860730588116458e-05,
    -1.72334302366962267e-05,
    5.60690482304602267e-06,
    2.31395443148286800e-05,
    3.62642745856793957e-05,
    4.58006124490188752e-05,
    5.24595294959114050e-05,
    5.68396208545815266e-05,
    5.94349820393104052e-05,
    6.06478527578421742e-05,
    6.08023907788436497e-05,
    6.01577894539460388e-05,
    5.89199657344698500e-05,
    5.72515823777593053e-05,
    5.52804375585852577e-05,
    5.31063773802880170e-05,
    5.08069302012325706e-05,
    4.84418647620094842e-05,
    4.60568581607475370e-05,
    -6.91141397288294174e-04,
    -4.29976633058871912e-04,
    1.83067735980039018e-04,
    6.60088147542014144e-04,
    8.75964969951185931e-04,
    8.77335235958235514e-04,
    7.49369585378990637e-04,
    5.63832329756980918e-04,
    3.68059319971443156e-04,
    1.88464535514455599e-04,
    3.70663057664904149e-05,
    -8.28520220232137023e-05,
    -1.72751952869172998e-04,
    -2.36314873605872983e-04,
    -2.77966150694906658e-04,
    -3.02079514155456919e-04,
    -3.12594712643820127e-04,
    -3.12872558758067163e-04,
    -3.05678038466324377e-04,
    -2.93226470614557331e-04,
    -2.77255655582934777e-04,
    -2.59103928467031709e-04,
    -2.39784014396480342e-04,
    -2.20048260045422848e-04,
    -2.00443911094971498e-04,
    -1.81358692210970687e-04,
    -1.63057674478657464e-04,
    -1.45712672175205844e-04,
    -1.29425421983924587e-04,
    -1.14245691942445952e-04,
    1.92821964248775885e-03,
    1.35592576302022234e-03,
    -7.17858090421302995e-04,
    -2.58084802575270346e-03,
    -3.49271130826168475e-03,
    -3.46986299340960628e-03,
    -2.82285233351310182e-03,
    -1.88103076404891354e-03,
    -8.89531718383947600e-04,
    3.87912102631035228e-06,
    7.28688540119691412e-04,
    1.26566373053457758e-03,
    1.62518158372674427e-03,
    1.83203153216373172e-03,
    1.91588388990527909e-03,
    1.90588846755546138e-03,
    1.82798982421825727e-03,
    1.70389506421121530e-03,
    1.55097127171097686e-03,
    1.38261421852276159e-03,
    1.20881424230064774e-03,
    1.03676532638344962e-03,
    8.71437918068619115e-04,
    7.16080155297701002e-04,
    5.72637002558129372e-04,
    4.42089819465802277e-04,
    3.24724948503090564e-04,
    2.20342042730246599e-04,
    1.28412898401353882e-04,
    4.82005924552095464e-05,
];
const BETA: [f64; 210] = [
    1.79988721413553309e-02,
    5.59964911064388073e-03,
    2.88501402231132779e-03,
    1.80096606761053941e-03,
    1.24753110589199202e-03,
    9.22878876572938311e-04,
    7.14430421727287357e-04,
    5.71787281789704872e-04,
    4.69431007606481533e-04,
    3.93232835462916638e-04,
    3.34818889318297664e-04,
    2.88952148495751517e-04,
    2.52211615549573284e-04,
    2.22280580798883327e-04,
    1.97541838033062524e-04,
    1.76836855019718004e-04,
    1.59316899661821081e-04,
    1.44347930197333986e-04,
    1.31448068119965379e-04,
    1.20245444949302884e-04,
    1.10449144504599392e-04,
    1.01828770740567258e-04,
    9.41998224204237509e-05,
    8.74130545753834437e-05,
    8.13466262162801467e-05,
    7.59002269646219339e-05,
    7.09906300634153481e-05,
    6.65482874842468183e-05,
    6.25146958969275078e-05,
    5.88403394426251749e-05,
    -1.49282953213429172e-03,
    -8.78204709546389328e-04,
    -5.02916549572034614e-04,
    -2.94822138512746025e-04,
    -1.75463996970782828e-04,
    -1.04008550460816434e-04,
    -5.96141953046457895e-05,
    -3.12038929076098340e-05,
    -1.26089735980230047e-05,
    -2.42892608575730389e-07,
    8.05996165414273571e-06,
    1.36507009262147391e-05,
    1.73964125472926261e-05,
    1.98672978842133780e-05,
    2.14463263790822639e-05,
    2.23954659232456514e-05,
    2.28967783814712629e-05,
    2.30785389811177817e-05,
    2.30321976080909144e-05,
    2.28236073720348722e-05,
    2.25005881105292418e-05,
    2.20981015361991429e-05,
    2.16418427448103905e-05,
    2.11507649256220843e-05,
    2.06388749782170737e-05,
    2.01165241997081666e-05,
    1.95913450141179244e-05,
    1.90689367910436740e-05,
    1.85533719641636667e-05,
    1.80475722259674218e-05,
    5.52213076721292790e-04,
    4.47932581552384646e-04,
    2.79520653992020589e-04,
    1.52468156198446602e-04,
    6.93271105657043598e-05,
    1.76258683069991397e-05,
    -1.35744996343269136e-05,
    -3.17972413350427135e-05,
    -4.18861861696693365e-05,
    -4.69004889379141029e-05,
    -4.87665447413787352e-05,
    -4.87010031186735069e-05,
    -4.74755620890086638e-05,
    -4.55813058138628452e-05,
    -4.33309644511266036e-05,
    -4.09230193157750364e-05,
    -3.84822638603221274e-05,
    -3.60857167535410501e-05,
    -3.37793306123367417e-05,
    -3.15888560772109621e-05,
    -2.95269561750807315e-05,
    -2.75978914828335759e-05,
    -2.58006174666883713e-05,
    -2.41308356761280200e-05,
    -2.25823509518346033e-05,
    -2.11479656768912971e-05,
    -1.98200638885294927e-05,
    -1.85909870801065077e-05,
    -1.74532699844210224e-05,
    -1.63997823854497997e-05,
    -4.74617796559959808e-04,
    -4.77864567147321487e-04,
    -3.20390228067037603e-04,
    -1.61105016119962282e-04,
    -4.25778101285435204e-05,
    3.44571294294967503e-05,
    7.97092684075674924e-05,
    1.03138236708272200e-04,
    1.12466775262204158e-04,
    1.13103642108481389e-04,
    1.08651634848774268e-04,
    1.01437951597661973e-04,
    9.29298396593363896e-05,
    8.40293133016089978e-05,
    7.52727991349134062e-05,
    6.69632521975730872e-05,
    5.92564547323194704e-05,
    5.22169308826975567e-05,
    4.58539485165360646e-05,
    4.01445513891486808e-05,
    3.50481730031328081e-05,
    3.05157995034346659e-05,
    2.64956119950516039e-05,
    2.29363633690998152e-05,
    1.97893056664021636e-05,
    1.70091984636412623e-05,
    1.45547428261524004e-05,
    1.23886640995878413e-05,
    1.04775876076583236e-05,
    8.79179954978479373e-06,
    7.36465810572578444e-04,
    8.72790805146193976e-04,
    6.22614862573135066e-04,
    2.85998154194304147e-04,
    3.84737672879366102e-06,
    -1.87906003636971558e-04,
    -2.97603646594554535e-04,
    -3.45998126832656348e-04,
    -3.53382470916037712e-04,
    -3.35715635775048757e-04,
    -3.04321124789039809e-04,
    -2.66722723047612821e-04,
    -2.27654214122819527e-04,
    -1.89922611854562356e-04,
    -1.55058918599093870e-04,
    -1.23778240761873630e-04,
    -9.62926147717644187e-05,
    -7.25178327714425337e-05,
    -5.22070028895633801e-05,
    -3.50347750511900522e-05,
    -2.06489761035551757e-05,
    -8.70106096849767054e-06,
    1.13698686675100290e-06,
    9.16426474122778849e-06,
    1.56477785428872620e-05,
    2.08223629482466847e-05,
    2.48923381004595156e-05,
    2.80340509574146325e-05,
    3.03987774629861915e-05,
    3.21156731406700616e-05,
    -1.80182191963885708e-03,
    -2.43402962938042533e-03,
    -1.83422663549856802e-03,
    -7.62204596354009765e-04,
    2.39079475256927218e-04,
    9.49266117176881141e-04,
    1.34467449701540359e-03,
    1.48457495259449178e-03,
    1.44732339830617591e-03,
    1.30268261285657186e-03,
    1.10351597375642682e-03,
    8.86047440419791759e-04,
    6.73073208165665473e-04,
    4.77603872856582378e-04,
    3.05991926358789362e-04,
    1.60315694594721630e-04,
    4.00749555270613286e-05,
    -5.66607461635251611e-05,
    -1.32506186772982638e-04,
    -1.90296187989614057e-04,
    -2.32811450376937408e-04,
    -2.62628811464668841e-04,
    -2.82050469867598672e-04,
    -2.93081563192861167e-04,
    -2.97435962176316616e-04,
    -2.96557334239348078e-04,
    -2.91647363312090861e-04,
    -2.83696203837734166e-04,
    -2.73512317095673346e-04,
    -2.61750155806768580e-04,
    6.38585891212050914e-03,
    9.62374215806377941e-03,
    7.61878061207001043e-03,
    2.83219055545628054e-03,
    -2.09841352012720090e-03,
    -5.73826764216626498e-03,
    -7.70804244495414620e-03,
    -8.21011692264844401e-03,
    -7.65824520346905413e-03,
    -6.47209729391045177e-03,
    -4.99132412004966473e-03,
    -3.45612289713133280e-03,
    -2.01785580014170775e-03,
    -7.59430686781961401e-04,
    2.84173631523859138e-04,
    1.10891667586337403e-03,
    1.72901493872728771e-03,
    2.16812590802684701e-03,
    2.45357710494539735e-03,
    2.61281821058334862e-03,
    2.67141039656276912e-03,
    2.65203073395980430e-03,
    2.57411652877287315e-03,
    2.45389126236094427e-03,
    2.30460058071795494e-03,
    2.13684837686712662e-03,
    1.95896528478870911e-03,
    1.77737008679454412e-03,
    1.59690280765839059e-03,
    1.42111975664438546e-03,
];
const GAMA: [f64; 30] = [
    6.29960524947436582e-01,
    2.51984209978974633e-01,
    1.54790300415655846e-01,
    1.10713062416159013e-01,
    8.57309395527394825e-02,
    6.97161316958684292e-02,
    5.86085671893713576e-02,
    5.04698873536310685e-02,
    4.42600580689154809e-02,
    3.93720661543509966e-02,
    3.54283195924455368e-02,
    3.21818857502098231e-02,
    2.94646240791157679e-02,
    2.71581677112934479e-02,
    2.51768272973861779e-02,
    2.34570755306078891e-02,
    2.19508390134907203e-02,
    2.06210828235646240e-02,
    1.94388240897880846e-02,
    1.83810633800683158e-02,
    1.74293213231963172e-02,
    1.65685837786612353e-02,
    1.57865285987918445e-02,
    1.50729501494095594e-02,
    1.44193250839954639e-02,
    1.38184805735341786e-02,
    1.32643378994276568e-02,
    1.27517121970498651e-02,
    1.22761545318762767e-02,
    1.18338262398482403e-02,
];
