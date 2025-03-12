use std::f64::consts::{FRAC_PI_2, PI};

use num::{
    One, Zero,
    complex::{Complex64, ComplexFloat},
};

use crate::amos::{machine::d1mach, utils::will_z_underflow};

use super::{BesselError, BesselError::*, Scaling};

pub fn ZUOIK(
    z: Complex64,
    order: f64, //ZR, ZI, FNU,
    KODE: Scaling,
    IKFLG: usize,
    N: usize, //YR, YI, NUF,
    mut y: Vec<Complex64>,
    TOL: f64,
    ELIM: f64,
    ALIM: f64,
) -> Result<(Vec<Complex64>, usize), BesselError> {
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
    //     COMPLEX ARG,ASUM,BSUM,CWRK,CZ,CZERO,PHI,SUM,Y,Z,ZB,ZETA1,ZETA2,ZN,
    //    *ZR
    //       EXTERNAL ZABS
    //       DOUBLE PRECISION AARG, AIC, ALIM, APHI, ARGI, ARGR, ASUMI, ASUMR,
    //      * ASCLE, AX, AY, BSUMI, BSUMR, CWRKI, CWRKR, CZI, CZR, ELIM, FNN,
    //      * FNU, GNN, GNU, PHII, PHIR, RCZ, STR, STI, SUMI, SUMR, TOL, YI,
    //      * YR, ZBI, ZBR, ZEROI, ZEROR, ZETA1I, ZETA1R, ZETA2I, ZETA2R, ZI,
    //      * ZNI, ZNR, ZR, ZRI, ZRR, d1mach, ZABS
    //       INTEGER I, IDUM, IFORM, IKFLG, INIT, KODE, N, NN, NUF, NW
    //       DIMENSION YR(N), YI(N), CWRKR(16), CWRKI(16)
    //       DATA ZEROR,ZEROI / 0.0, 0.0 /
    const AIC: f64 = 1.265512123484645396e+00;
    let mut NUF = 0;
    let mut NN = N;
    //     let zr = z;
    // ZRR = ZR;
    // ZRI = ZI;
    let zr = if z.re < 0.0 //GO TO 10;
      {-z}else{z};
    // ZRR = -ZR;
    // ZRI = -ZI;

    //    10 CONTINUE;
    let zb = zr;
    // ZBR = ZRR;
    // ZBI = ZRI;
    let AX = z.re.abs() * 1.7321;
    let AY = z.im.abs();
    let IFORM = if AY > AX { 2 } else { 1 };
    let mut GNU = order.max(1.0);
    if IKFLG != 1
    //GO TO 20;
    {
        let FNN = NN as f64;
        let GNN = order + FNN - 1.0;
        GNU = GNN.max(FNN);
    }
    //    20 CONTINUE;
    //-----------------------------------------------------------------------;
    //     ONLY THE MAGNITUDE OF ARG AND PHI ARE NEEDED ALONG WITH THE;
    //     REAL PARTS OF ZETA1, ZETA2 AND ZB. NO ATTEMPT IS MADE TO GET;
    //     THE SIGN OF THE IMAGINARY PART CORRECT.;
    //-----------------------------------------------------------------------;
    let mut zn = None;
    let (mut cz, phi, arg, AARG) = if IFORM != 2 {
        //GO TO 30;
        let INIT = 0;
        let (phi, zeta1, zeta2, _) = ZUNIK(
            zr, /*ZRR, ZRI*/ GNU, IKFLG, true, TOL,
            INIT, /* , PHIR, PHII,
                 ZETA1R, ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI*/
        );
        (-zeta1 + zeta2, phi, Complex64::zero(), 0.0)
    // CZR = -ZETA1R + ZETA2R;
    // CZI = -ZETA1I + ZETA2I;
    // GO TO 50;
    //    30 CONTINUE;
    } else {
        let mut zn_ = Complex64::new(zr.im, -zr.re);
        // ZNR = ZRI;
        // ZNI = -ZRR;
        if z.im <= 0.0
        //GO TO 40;
        {
            zn_.re = -zn_.re;
        }
        //    40 CONTINUE;
        let (phi, arg, zeta1, zeta2, _, _) = ZUNHJ(
            //ZNR, ZNI,
            zn_, GNU, true,
            TOL,
            //       PHIR, PHII, ARGR, ARGI, ZETA1R,
            // ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI
        );
        zn = Some(zn_);
        // CZR = -ZETA1R + ZETA2R;
        // CZI = -ZETA1I + ZETA2I;
        let cz = -zeta1 + zeta2;
        let AARG = arg.abs(); //ZABS(ARGR,ARGI);
        (cz, phi, arg, AARG)
    };
    //    50 CONTINUE;
    if KODE == Scaling::Scaled {
        //GO TO 60;
        cz -= zb;
        // CZR = CZR - ZBR;
        // CZI = CZI - ZBI;
    }
    //    60 CONTINUE;
    if IKFLG != 1 {
        // GO TO 70;
        cz = -cz;
        // CZR = -CZR;
        // CZI = -CZI;
        //    70 CONTINUE;
    }
    let aphi = phi.abs();
    // APHI = ZABS(PHIR,PHII);
    let mut RCZ = cz.re;
    // RCZ = CZR;
    //-----------------------------------------------------------------------;
    //     OVERFLOW TEST;
    //-----------------------------------------------------------------------;
    if RCZ > ELIM {
        return Err(Overflow);
    } //GO TO 210;
    if RCZ >= ALIM {
        //GO TO 80;
        // RCZ = RCZ + DLOG(APHI);
        RCZ += aphi.ln();
        if IFORM == 2 {
            RCZ = RCZ - 0.25 * AARG.ln() - AIC
        };
        if RCZ > ELIM {
            return Err(Overflow);
        } //GO TO 210;
    // GO TO 130;
    } else {
        //    80 CONTINUE;
        //-----------------------------------------------------------------------;
        //     UNDERFLOW TEST;
        //-----------------------------------------------------------------------;
        if RCZ < (-ELIM) {
            y[0..NN].iter_mut().for_each(|v| *v = Complex64::zero());
            return Ok((y, NN));
        }
        //GO TO 90;
        if RCZ <= (-ALIM) {
            //GO TO 130;
            RCZ += aphi.ln();
            if IFORM == 2 {
                RCZ = RCZ - 0.25 * AARG.ln() - AIC
            };
            if !(RCZ > (-ELIM)) {
                //GO TO 110;

                //    90 CONTINUE;
                //       DO 100 I=1,NN;
                //         YR(I) = ZEROR;
                //         YI(I) = ZEROI;
                //   100 CONTINUE;
                //       NUF = NN;
                //       RETURN;
                return Ok((y, NN));
            }
            //   110 CONTINUE;
            let ASCLE = 1.0e+3 * d1mach(1) / TOL;
            cz += phi.ln();
            // CALL ZLOG(PHIR, PHII, STR, STI, IDUM);
            // CZR = CZR + STR;
            // CZI = CZI + STI;
            if IFORM != 1 {
                //GO TO 120;
                cz -= 0.25 * arg.ln() + AIC
                // CALL ZLOG(ARGR, ARGI, STR, STI, IDUM);
                // CZR = CZR - 0.25*STR - AIC;
                // CZI = CZI - 0.25*STI;
            } //120 CONTINUE;
            let AX = RCZ.exp() / TOL;
            // AX = DEXP(RCZ)/TOL;
            let AY = cz.im;
            cz = AX * Complex64::new(AY.cos(), AY.sin());
            // CZR = AX*DCOS(AY);
            // CZI = AX*DSIN(AY);
            if will_z_underflow(cz, ASCLE, TOL) {
                y[0..NN].iter_mut().for_each(|v| *v = Complex64::zero());
                return Ok((y, NN));
            }
            // if (NW != 0) {}//GO TO 90;
        }
    }
    //   130 CONTINUE;
    if IKFLG == 2 || N == 1 {
        return Ok((y, NUF));
    }
    // if (IKFLG == 2) RETURN;
    // if (N == 1) RETURN;
    //-----------------------------------------------------------------------;
    //     SET UNDERFLOWS ON I SEQUENCE;
    //-----------------------------------------------------------------------;
    //   140 CONTINUE;
    let mut go_to_180 = false;
    let mut skip_to_190;
    'outer: loop {
        'l140: loop {
            skip_to_190 = false;
            if !go_to_180 {
                GNU = order + ((NN - 1) as f64);
                let (phi, mut cz, AARG) = if IFORM != 2 {
                    //GO TO 150;
                    let INIT = 0;
                    let (phi, zeta1, zeta2, _) = ZUNIK(zr, GNU, IKFLG, true, TOL, INIT);
                    //       CALL ZUNIK(ZRR, ZRI, GNU, IKFLG, 1, TOL, INIT, PHIR, PHII,;
                    //      * ZETA1R, ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI);
                    cz = -zeta1 + zeta2;
                    // CZR = -ZETA1R + ZETA2R;
                    // CZI = -ZETA1I + ZETA2I;
                    // GO TO 160;
                    (phi, cz, 0.0)
                } else {
                    //   150 CONTINUE;
                    let (phi, arg, zeta1, zeta2, _, _) = ZUNHJ(zn.unwrap(), GNU, true, TOL);
                    //       CALL ZUNHJ(ZNR, ZNI, GNU, 1, TOL, PHIR, PHII, ARGR, ARGI, ZETA1R,;
                    //      * ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI);
                    cz = -zeta1 + zeta2;
                    // CZR = -ZETA1R + ZETA2R;
                    // CZI = -ZETA1I + ZETA2I;
                    let AARG = arg.abs();
                    // AARG = ZABS(ARGR,ARGI);
                    (phi, cz, AARG)
                };
                //   160 CONTINUE;
                if KODE == Scaling::Scaled {
                    //GO TO 170;
                    cz -= zb;
                    // CZR = CZR - ZBR;
                    // CZI = CZI - ZBI;
                }
                //   170 CONTINUE;
                let aphi = phi.abs();
                // APHI = ZABS(PHIR,PHII);
                RCZ = cz.re;
                if !(RCZ < (-ELIM)) {
                    //GO TO 180;
                    if (RCZ > (-ALIM)) {
                        return Ok((y, NUF));
                    };
                    // RCZ = RCZ + DLOG(APHI);
                    RCZ += aphi.ln();
                    if (IFORM == 2) {
                        RCZ = RCZ - 0.25 * AARG.ln() - AIC;
                    }
                    if (RCZ > (-ELIM)) {
                        skip_to_190 = true
                    } //GO TO 190
                }
                //   180 CONTINUE;
            }
            go_to_180 = false;
            if !skip_to_190 {
                y[NN - 1] = Complex64::zero();
                // YR(NN) = ZEROR;
                // YI(NN) = ZEROI;
                NN = NN - 1;
                NUF = NUF + 1;
                if NN == 0 {
                    return Ok((y, NUF));
                }
                // GO TO 140;
            } else {
                break 'l140;
            }
        }
        //   190 CONTINUE;
        let ASCLE = 1.0e+3 * d1mach(1) / TOL;
        cz += phi.ln();
        // CALL ZLOG(PHIR, PHII, STR, STI, IDUM);
        // CZR = CZR + STR;
        // CZI = CZI + STI;
        if IFORM != 1 {
            //GO TO 200;
            cz -= 0.25 * arg.ln() + AIC
            // CALL ZLOG(ARGR, ARGI, STR, STI, IDUM);
            // CZR = CZR - 0.25*STR - AIC;
            // CZI = CZI - 0.25*STI;
        }
        //   200 CONTINUE;
        let AX = RCZ.exp() / TOL;
        let AY = cz.im;
        cz = AX * Complex64::new(AY.cos(), AY.sin());
        // CZR = AX*DCOS(AY);
        // CZI = AX*DSIN(AY);
        //ZUunderflowCHK(CZR, CZI, NW, ASCLE, TOL);
        if will_z_underflow(cz, ASCLE, TOL) {
            go_to_180 = true;
        } else {
            break 'outer;
        }
    }
    // if (NW != 0) GO TO 180;
    return Ok((y, NUF));
    //   210 CONTINUE;
    //   return Err(Overflow);
    // NUF = -1;
    // RETURN;
    // END;
}

fn ZUNIK(
    /*ZRR, ZRI, FNU,*/
    zr: Complex64,
    order: f64,
    IKFLG: usize,
    only_phi_zeta: bool,
    TOL: f64,
    mut INIT: usize, //PHIR,
                     // PHII, ZETA1R, ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI
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
    //     COMPLEX CFN,CON,CONE,CRFN,CWRK,CZERO,PHI,S,SR,SUM,T,T2,ZETA1,
    //    *ZETA2,ZN,ZR
    //       DOUBLE PRECISION AC, C, CON, CONEI, CONER, CRFNI, CRFNR, CWRKI,
    //      * CWRKR, FNU, PHII, PHIR, RFN, SI, SR, SRI, SRR, STI, STR, SUMI,
    //      * SUMR, TEST, TI, TOL, TR, T2I, T2R, ZEROI, ZEROR, ZETA1I, ZETA1R,
    //      * ZETA2I, ZETA2R, ZNI, ZNR, ZRI, ZRR, d1mach
    // INTEGER I, IDUM, IKFLG, INIT, IPMTR, J, K, L
    // DIMENSION C(120), CWRKR(16), CWRKI(16), CON(2)
    // DATA ZEROR,ZEROI,CONER,CONEI / 0.0, 0.0, 1.0, 0.0 /
    const CON: [f64; 2] = //CON(1), CON(2)  /
        [3.98942280401432678e-01, 1.25331413731550025e+00];
    const C: [f64; 120] = [
        //       DATA C(1), C(2), C(3), C(4), C(5), C(6), C(7), C(8), C(9), C(10),
        //      1     C(11), C(12), C(13), C(14), C(15), C(16), C(17), C(18),
        //      2     C(19), C(20), C(21), C(22), C(23), C(24)/
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
        // DATA C(25), C(26), C(27), C(28), C(29), C(30), C(31), C(32),
        //     C(33), C(34), C(35), C(36), C(37), C(38), C(39), C(40),
        //     C(41), C(42), C(43), C(44), C(45), C(46), C(47), C(48)/
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
        // DATA C(49), C(50), C(51), C(52), C(53), C(54), C(55), C(56),
        //     C(57), C(58), C(59), C(60), C(61), C(62), C(63), C(64),
        //     C(65), C(66), C(67), C(68), C(69), C(70), C(71), C(72)/
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
        // DATA C(73), C(74), C(75), C(76), C(77), C(78), C(79), C(80),
        //     C(81), C(82), C(83), C(84), C(85), C(86), C(87), C(88),
        //     C(89), C(90), C(91), C(92), C(93), C(94), C(95), C(96)/
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
        // DATA C(97), C(98), C(99), C(100), C(101), C(102), C(103), C(104),
        //     C(105), C(106), C(107), C(108), C(109), C(110), C(111),
        //     C(112), C(113), C(114), C(115), C(116), C(117), C(118)/
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
        // DATA C(119), C(120)/
        -2.91883881222208134e+07,
        1.18838426256783253e+05,
    ];
    //
    let mut working = vec![Complex64::zero(); 16];

    // if (INIT != 0) GO TO 40;
    //-----------------------------------------------------------------------;
    //     INITIALIZE ALL VARIABLES;
    //-----------------------------------------------------------------------;
    let RFN = 1.0 / order;
    //-----------------------------------------------------------------------;
    //     OVERFLOW TEST (ZR/FNU TOO SMALL);
    //-----------------------------------------------------------------------;
    let TEST = d1mach(1) * 1.0e+3;
    let AC = order * TEST;
    if !(zr.re.abs() > AC || zr.im.abs() > AC) {
        let zeta1 = Complex64::new(2.0 * TEST.ln().abs() + order, 0.0);
        // ZETA1R = 2.0*DABS(DLOG(TEST))+FNU;
        // ZETA1I = 0.0;
        let zeta2 = Complex64::new(order, 0.0);
        // ZETA2R = FNU;
        // ZETA2I = 0.0;
        let phi = Complex64::one();
        // PHIR = 1.0;
        // PHII = 0.0;
        return (zeta1, zeta2, phi, None);
    }
    //    15 CONTINUE;
    let t = zr * RFN;
    // TR = ZRR*RFN;
    // TI = ZRI*RFN;
    let s = Complex64::one() + t * t;
    // SR = CONER + (TR*TR-TI*TI);
    // SI = CONEI + (TR*TI+TI*TR);
    let sr = s.sqrt();
    // CALL ZSQRT(SR, SI, SRR, SRI);
    let st = Complex64::one() + sr;
    // STR = CONER + SRR;
    // STI = CONEI + SRI;
    let zn = st / t;
    // CALL ZDIV(STR, STI, TR, TI, ZNR, ZNI);
    // CALL ZLOG(ZNR, ZNI, STR, STI, IDUM);
    let st = zn.ln();
    let zeta1 = order * st;
    // ZETA1R = FNU*STR;
    // ZETA1I = FNU*STI;
    let zeta2 = order * sr;
    // ZETA2R = FNU*SRR;
    // ZETA2I = FNU*SRI;
    // CALL ZDIV(CONER, CONEI, SRR, SRI, TR, TI);
    let t = Complex64::one() / sr;
    let sr = t * RFN;
    // SRR = TR*RFN;
    // SRI = TI*RFN;
    working[15] = sr.sqrt();
    // CALL ZSQRT(SRR, SRI, CWRKR(16), CWRKI(16));
    let mut phi = working[15] * CON[IKFLG - 1];
    // PHIR = CWRKR(16)*CON(IKFLG);
    // PHII = CWRKI(16)*CON(IKFLG);
    if only_phi_zeta {
        return (phi, zeta1, zeta2, None);
    };
    let t2 = Complex64::one() / s;
    // CALL ZDIV(CONER, CONEI, SR, SI, T2R, T2I);
    working[0] = Complex64::one();
    // CWRKR(1) = CONER;
    // CWRKI(1) = CONEI;
    let mut crfn = Complex64::one();
    // CRFNR = CONER;
    // CRFNI = CONEI;
    let mut AC = 1.0;
    let mut L = 1;
    // DO 20 K=2,15;
    let mut K = 0;
    'l20: for k in 1..15 {
        K = k;
        let mut s = Complex64::zero();
        //   SR = ZEROR;
        //   SI = ZEROI;
        //   DO 10 J=1,K;
        '_l10: for _ in 0..k {
            L += 1;
            s = s * t2 + C[L];
            //     STR = SR*T2R - SI*T2I + C(L);
            //     s.im = s.re*t2.im + s.im*t2.re;
            //     SR = STR;
        }
        //    10   CONTINUE;
        //   STR = CRFNR*SRR - CRFNI*SRI;
        //   CRFNI = CRFNR*SRI + CRFNI*SRR;
        //   CRFNR = STR;
        crfn *= sr;
        working[k] = crfn * s;
        //   CWRKR(K) = CRFNR*SR - CRFNI*SI;
        //   CWRKI(K) = CRFNR*SI + CRFNI*SR;
        AC *= RFN;
        let TEST = working[k].re.abs() + working[k].im.abs();
        //   TEST = DABS(CWRKR(K)) + DABS(CWRKI(K));
        if AC < TOL && TEST < TOL {
            break 'l20;
        } //GO TO 30;
    }
    //    20 CONTINUE;
    // K = 15;
    //    30 CONTINUE;
    INIT = K;
    //    40 CONTINUE;
    //-----------------------------------------------------------------------;
    // FINISH INIT. Use OnceCell here later?
    //-----------------------------------------------------------------------;
    return if IKFLG == 1 {
        //GO TO 60;
        //-----------------------------------------------------------------------;
        //     COMPUTE SUM FOR THE I FUNCTION;
        //-----------------------------------------------------------------------;
        // s=Complex64::zero();
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
        let sum = working[..INIT].iter().sum();
        phi = working[15] * CON[0];
        // PHIR = CWRKR(16)*CON(1);
        // PHII = CWRKI(16)*CON(1);
        // RETURN;
        (phi, zeta1, zeta2, Some(sum))
    } else if IKFLG == 2 {
        //    60 CONTINUE;
        //-----------------------------------------------------------------------;
        //     COMPUTE SUM FOR THE K FUNCTION;
        //-----------------------------------------------------------------------;
        let mut tr = 1.0;
        let sum = working[..INIT]
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
    } else {
        panic!("Invalid value of IKFLG")
    };
}

fn ZUNHJ(
    z: Complex64,
    order: f64,          //ZR, ZI, FNU,
    only_phi_zeta: bool, //IPMTR,
    TOL: f64,
    //PHIR, PHII, ARGR, ARGI, ZETA1R, ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI
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
    //     COMPLEX ARG,ASUM,BSUM,CFNU,CONE,CR,CZERO,DR,P,PHI,PRZTH,PTFN,
    //    *RFN13,RTZTA,RZTH,SUMA,SUMB,TFN,T2,UP,W,W2,Z,ZA,ZB,ZC,ZETA,ZETA1,
    //    *ZETA2,ZTH
    //       EXTERNAL ZABS
    //       DOUBLE PRECISION ALFA, ANG, AP, AR, ARGI, ARGR, ASUMI, ASUMR,
    //      * ATOL, AW2, AZTH, BETA, BR, BSUMI, BSUMR, BTOL, C, CONEI, CONER,
    //      * CRI, CRR, DRI, DRR, EX1, EX2, FNU, FN13, FN23, GAMA, GPI, FRAC_PI_2,
    //      * PHII, PHIR, PI, PP, PR, PRZTHI, PRZTHR, PTFNI, PTFNR, RAW, RAW2,
    //      * RAZTH, RFNU, RFNU2, RFN13, RTZTI, RTZTR, RZTHI, RZTHR, STI, STR,
    //      * SUMAI, SUMAR, SUMBI, SUMBR, TEST, TFNI, TFNR, TFRAC_PI_2, TOL, TZAI,
    //      * TZAR, T2I, T2R, UPI, UPR, WI, WR, W2I, W2R, ZAI, ZAR, ZBI, ZBR,
    //      * ZCI, ZCR, ZEROI, ZEROR, ZETAI, ZETAR, ZETA1I, ZETA1R, ZETA2I,
    //      * ZETA2R, ZI, ZR, ZTHI, ZTHR, ZABS, AC, d1mach
    //       INTEGER IAS, IBS, IPMTR, IS, J, JR, JU, K, KMAX, KP1, KS, L, LR,
    //      * LRP1, L1, L2, M, IDUM
    //       DIMENSION AR(14), BR(14), C(105), ALFA(180), BETA(210), GAMA(30),
    //      * AP(30), PR(30), PI(30), UPR(14), UPI(14), CRR(14), CRI(14),
    //      * DRR(14), DRI(14)
    //       DATA AR(1), AR(2), AR(3), AR(4), AR(5), AR(6), AR(7), AR(8),
    //      1     AR(9), AR(10), AR(11), AR(12), AR(13), AR(14)/
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
        //       DATA BR(1), BR(2), BR(3), BR(4), BR(5), BR(6), BR(7), BR(8),
        //      1     BR(9), BR(10), BR(11), BR(12), BR(13), BR(14)/
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
    const C: [f64; 105] = [
        //       DATA C(1), C(2), C(3), C(4), C(5), C(6), C(7), C(8), C(9), C(10),
        //      1     C(11), C(12), C(13), C(14), C(15), C(16), C(17), C(18),
        //      2     C(19), C(20), C(21), C(22), C(23), C(24)/
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
        // DATA C(25), C(26), C(27), C(28), C(29), C(30), C(31), C(32),
        //     C(33), C(34), C(35), C(36), C(37), C(38), C(39), C(40),
        //     C(41), C(42), C(43), C(44), C(45), C(46), C(47), C(48),
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
        // DATA C(49), C(50), C(51), C(52), C(53), C(54), C(55), C(56),
        //     C(57), C(58), C(59), C(60), C(61), C(62), C(63), C(64),
        //     C(65), C(66), C(67), C(68), C(69), C(70), C(71), C(72),
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
        // DATA C(73), C(74), C(75), C(76), C(77), C(78), C(79), C(80),
        //     C(81), C(82), C(83), C(84), C(85), C(86), C(87), C(88),
        //     C(89), C(90), C(91), C(92), C(93), C(94), C(95), C(96),
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
        // DATA C(97), C(98), C(99), C(100), C(101), C(102), C(103), C(104),
        //     C(105),
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
        // DATA ALFA(1), ALFA(2), ALFA(3), ALFA(4), ALFA(5), ALFA(6),
        //     ALFA(7), ALFA(8), ALFA(9), ALFA(10), ALFA(11), ALFA(12),
        //     ALFA(13), ALFA(14), ALFA(15), ALFA(16), ALFA(17), ALFA(18),
        //     ALFA(19), ALFA(20), ALFA(21), ALFA(22),
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
        // DATA ALFA(23), ALFA(24), ALFA(25), ALFA(26), ALFA(27), ALFA(28),
        //     ALFA(29), ALFA(30), ALFA(31), ALFA(32), ALFA(33), ALFA(34),
        //     ALFA(35), ALFA(36), ALFA(37), ALFA(38), ALFA(39), ALFA(40),
        //     ALFA(41), ALFA(42), ALFA(43), ALFA(44),
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
        // DATA ALFA(45), ALFA(46), ALFA(47), ALFA(48), ALFA(49), ALFA(50),
        //     ALFA(51), ALFA(52), ALFA(53), ALFA(54), ALFA(55), ALFA(56),
        //     ALFA(57), ALFA(58), ALFA(59), ALFA(60), ALFA(61), ALFA(62),
        //     ALFA(63), ALFA(64), ALFA(65), ALFA(66),
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
        // DATA ALFA(67), ALFA(68), ALFA(69), ALFA(70), ALFA(71), ALFA(72),
        //     ALFA(73), ALFA(74), ALFA(75), ALFA(76), ALFA(77), ALFA(78),
        //     ALFA(79), ALFA(80), ALFA(81), ALFA(82), ALFA(83), ALFA(84),
        //     ALFA(85), ALFA(86), ALFA(87), ALFA(88),
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
        // DATA ALFA(89), ALFA(90), ALFA(91), ALFA(92), ALFA(93), ALFA(94),
        //     ALFA(95), ALFA(96), ALFA(97), ALFA(98), ALFA(99), ALFA(100),
        //     ALFA(101), ALFA(102), ALFA(103), ALFA(104), ALFA(105),
        //     ALFA(106), ALFA(107), ALFA(108), ALFA(109), ALFA(110),
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
        // DATA ALFA(111), ALFA(112), ALFA(113), ALFA(114), ALFA(115),
        //     ALFA(116), ALFA(117), ALFA(118), ALFA(119), ALFA(120),
        //     ALFA(121), ALFA(122), ALFA(123), ALFA(124), ALFA(125),
        //     ALFA(126), ALFA(127), ALFA(128), ALFA(129), ALFA(130),
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
        // DATA ALFA(131), ALFA(132), ALFA(133), ALFA(134), ALFA(135),
        //     ALFA(136), ALFA(137), ALFA(138), ALFA(139), ALFA(140),
        //     ALFA(141), ALFA(142), ALFA(143), ALFA(144), ALFA(145),
        //     ALFA(146), ALFA(147), ALFA(148), ALFA(149), ALFA(150),
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
        // DATA ALFA(151), ALFA(152), ALFA(153), ALFA(154), ALFA(155),
        //     ALFA(156), ALFA(157), ALFA(158), ALFA(159), ALFA(160),
        //     ALFA(161), ALFA(162), ALFA(163), ALFA(164), ALFA(165),
        //     ALFA(166), ALFA(167), ALFA(168), ALFA(169), ALFA(170),
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
        // DATA ALFA(171), ALFA(172), ALFA(173), ALFA(174), ALFA(175),
        //     ALFA(176), ALFA(177), ALFA(178), ALFA(179), ALFA(180),
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
        // DATA BETA(1), BETA(2), BETA(3), BETA(4), BETA(5), BETA(6),
        //     BETA(7), BETA(8), BETA(9), BETA(10), BETA(11), BETA(12),
        //     BETA(13), BETA(14), BETA(15), BETA(16), BETA(17), BETA(18),
        //     BETA(19), BETA(20), BETA(21), BETA(22),
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
        // DATA BETA(23), BETA(24), BETA(25), BETA(26), BETA(27), BETA(28),
        //     BETA(29), BETA(30), BETA(31), BETA(32), BETA(33), BETA(34),
        //     BETA(35), BETA(36), BETA(37), BETA(38), BETA(39), BETA(40),
        //     BETA(41), BETA(42), BETA(43), BETA(44),
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
        // DATA BETA(45), BETA(46), BETA(47), BETA(48), BETA(49), BETA(50),
        //     BETA(51), BETA(52), BETA(53), BETA(54), BETA(55), BETA(56),
        //     BETA(57), BETA(58), BETA(59), BETA(60), BETA(61), BETA(62),
        //     BETA(63), BETA(64), BETA(65), BETA(66),
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
        // DATA BETA(67), BETA(68), BETA(69), BETA(70), BETA(71), BETA(72),
        //     BETA(73), BETA(74), BETA(75), BETA(76), BETA(77), BETA(78),
        //     BETA(79), BETA(80), BETA(81), BETA(82), BETA(83), BETA(84),
        //     BETA(85), BETA(86), BETA(87), BETA(88),
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
        // DATA BETA(89), BETA(90), BETA(91), BETA(92), BETA(93), BETA(94),
        //     BETA(95), BETA(96), BETA(97), BETA(98), BETA(99), BETA(100),
        //     BETA(101), BETA(102), BETA(103), BETA(104), BETA(105),
        //     BETA(106), BETA(107), BETA(108), BETA(109), BETA(110),
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
        // DATA BETA(111), BETA(112), BETA(113), BETA(114), BETA(115),
        //     BETA(116), BETA(117), BETA(118), BETA(119), BETA(120),
        //     BETA(121), BETA(122), BETA(123), BETA(124), BETA(125),
        //     BETA(126), BETA(127), BETA(128), BETA(129), BETA(130),
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
        // DATA BETA(131), BETA(132), BETA(133), BETA(134), BETA(135),
        //     BETA(136), BETA(137), BETA(138), BETA(139), BETA(140),
        //     BETA(141), BETA(142), BETA(143), BETA(144), BETA(145),
        //     BETA(146), BETA(147), BETA(148), BETA(149), BETA(150),
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
        // DATA BETA(151), BETA(152), BETA(153), BETA(154), BETA(155),
        //     BETA(156), BETA(157), BETA(158), BETA(159), BETA(160),
        //     BETA(161), BETA(162), BETA(163), BETA(164), BETA(165),
        //     BETA(166), BETA(167), BETA(168), BETA(169), BETA(170),
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
        // DATA BETA(171), BETA(172), BETA(173), BETA(174), BETA(175),
        //     BETA(176), BETA(177), BETA(178), BETA(179), BETA(180),
        //     BETA(181), BETA(182), BETA(183), BETA(184), BETA(185),
        //     BETA(186), BETA(187), BETA(188), BETA(189), BETA(190),
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
        // DATA BETA(191), BETA(192), BETA(193), BETA(194), BETA(195),
        //     BETA(196), BETA(197), BETA(198), BETA(199), BETA(200),
        //     BETA(201), BETA(202), BETA(203), BETA(204), BETA(205),
        //     BETA(206), BETA(207), BETA(208), BETA(209), BETA(210),
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
        // DATA GAMA(1), GAMA(2), GAMA(3), GAMA(4), GAMA(5), GAMA(6),
        //     GAMA(7), GAMA(8), GAMA(9), GAMA(10), GAMA(11), GAMA(12),
        //     GAMA(13), GAMA(14), GAMA(15), GAMA(16), GAMA(17), GAMA(18),
        //     GAMA(19), GAMA(20), GAMA(21), GAMA(22),
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
        // DATA GAMA(23), GAMA(24), GAMA(25), GAMA(26), GAMA(27), GAMA(28),
        //     GAMA(29), GAMA(30),
        1.57865285987918445e-02,
        1.50729501494095594e-02,
        1.44193250839954639e-02,
        1.38184805735341786e-02,
        1.32643378994276568e-02,
        1.27517121970498651e-02,
        1.22761545318762767e-02,
        1.18338262398482403e-02,
    ];
    // DATA EX1, EX2, FRAC_PI_2, GPI,  ,
    //     3.33333333333333333e-01,     6.66666666666666667e-01,
    //     1.57079632679489662e+00,     3.14159265358979324e+00,
    const EX1: f64 = 3.33333333333333333e-01;
    const EX2: f64 = 6.66666666666666667e-01;
    const THREE_PI_BY_2: f64 = 4.71238898038468986e+00;
    // DATA ZEROR,ZEROI,CONER,CONEI / 0.0, 0.0, 1.0, 0.0 /
    //
    let RFNU = 1.0 / order;
    //-----------------------------------------------------------------------
    //     OVERFLOW TEST (Z/FNU TOO SMALL)
    //-----------------------------------------------------------------------
    let TEST = d1mach(1) * 1.0e+3;
    let AC = order * TEST;
    if !((z.re).abs() > AC || (z.im).abs() > AC) {
        // GO TO 15;
        let zeta1 = Complex64::new(2.0 * TEST.ln().abs() + order, 0.0);
        // ZETA1R = 2.0*DABS(DLOG(TEST))+FNU;
        // ZETA1I = 0.0;
        let zeta2 = Complex64::new(order, 0.0);
        // ZETA2R = FNU;
        // ZETA2I = 0.0;
        let phi = Complex64::one();
        let arg = Complex64::one();
        // PHIR = 1.0;
        // PHII = 0.0;
        // ARGR = 1.0;
        // ARGI = 0.0;
        // RETURN;
        return (phi, arg, zeta1, zeta2, None, None);
    }
    //    15 CONTINUE;
    let zb = z * RFNU;
    // ZBR = ZR*RFNU;
    // ZBI = ZI*RFNU;
    let RFNU2 = RFNU * RFNU;
    //-----------------------------------------------------------------------;
    //     COMPUTE IN THE FOURTH QUADRANT;
    //-----------------------------------------------------------------------;
    let FN13 = order.powf(EX1); //FNU**EX1;
    let FN23 = FN13 * FN13;
    let RFN13 = 1.0 / FN13;
    let w2 = Complex64::one() - zb * zb;
    // W2R = CONER - ZBR*ZBR + ZBI*ZBI;
    // W2I = CONEI - ZBR*ZBI - ZBR*ZBI;
    let AW2 = w2.abs();
    // AW2 = ZABS(W2R,W2I);

    if AW2 <= 0.25 {
        //GO TO 130;
        //-----------------------------------------------------------------------;
        //     POWER SERIES FOR CABS(W2) <= 0.25;
        //-----------------------------------------------------------------------;
        let mut K = 0;
        let mut p = vec![Complex64::zero(); 30];
        let mut ap = vec![0.0; 30];
        p[0] = Complex64::one();
        // let PR(1) = CONER;
        // PI(1) = CONEI;
        let mut suma = Complex64::new(GAMA[0], 0.0);
        // SUMAR = GAMA(1);
        // SUMAI = ZEROI;
        ap[0] = 1.0;
        if AW2 >= TOL {
            //GO TO 20;
            // DO 10 K=2,30;
            // let mut K = 0;
            for k in 1..30 {
                K = k;
                p[k] = p[k - 1] * w2;
                //   PR(K) = PR(K-1)*W2R - PI(K-1)*W2I;
                //   PI(K) = PR(K-1)*W2I + PI(K-1)*W2R;
                suma += p[k] * GAMA[k];
                //   SUMAR = SUMAR + PR(K)*GAMA(K);
                //   SUMAI = SUMAI + PI(K)*GAMA(K);
                ap[k] = ap[k - 1] * AW2;
                //   AP(K) = AP(K-1)*AW2;
                if ap[k] < TOL {
                    break;
                }
                //    10 CONTINUE;
            }
        }
        // K = 30;
        //    20 CONTINUE;
        let KMAX = K;
        let zeta = w2 * suma;
        // ZETAR = W2R*SUMAR - W2I*SUMAI;
        // ZETAI = W2R*SUMAI + W2I*SUMAR;
        let arg = zeta * FN23;
        // ARGR = ZETAR*FN23;
        // ARGI = ZETAI*FN23;
        let mut za = suma.sqrt();
        // CALL ZSQRT(SUMAR, SUMAI, ZAR, ZAI);
        // let st = w2.sqrt();
        // CALL ZSQRT(W2R, W2I, STR, STI);
        let zeta2 = w2.sqrt() * order;
        // ZETA2R = STR*FNU;
        // ZETA2I = STI*FNU;
        //  st = Complex64::one() + EX2 * zeta * za;
        // STR = CONER + EX2*(ZETAR*ZAR-ZETAI*ZAI);
        // STI = CONEI + EX2*(ZETAR*ZAI+ZETAI*ZAR);
        let zeta1 = (Complex64::one() + EX2 * zeta * za) * zeta2;
        // ZETA1R = STR*ZETA2R - STI*ZETA2I;
        // ZETA1I = STR*ZETA2I + STI*ZETA2R;
        za *= 2.0;
        // ZAR = ZAR + ZAR;
        // ZAI = ZAI + ZAI;
        // st = za.sqrt();
        // CALL ZSQRT(ZAR, ZAI, STR, STI);
        let phi = za.sqrt() * RFN13;
        // PHIR = STR*RFN13;
        // PHII = STI*RFN13;
        if only_phi_zeta {
            return (phi, arg, zeta1, zeta2, None, None);
        }
        //-----------------------------------------------------------------------;
        //     SUM SERIES FOR ASUM AND BSUM;
        //-----------------------------------------------------------------------;
        //       SUMBR = ZEROR;
        //       SUMBI = ZEROI;
        //       DO 30 K=1,KMAX;
        //         SUMBR = SUMBR + PR(K)*BETA(K);
        //         SUMBI = SUMBI + PI(K)*BETA(K);
        //    30 CONTINUE;
        let sumb: Complex64 = p[..KMAX].iter().zip(BETA).map(|(p, b)| p * b).sum();
        let mut asum = Complex64::zero();
        // ASUMR = ZEROR;
        // ASUMI = ZEROI;
        let mut bsum = sumb;
        // BSUMR = SUMBR;
        // BSUMI = SUMBI;
        let mut L1 = 0;
        let mut L2 = 30;
        let BTOL = TOL * (bsum.re.abs() + bsum.im.abs());
        let mut ATOL = TOL;
        let mut PP = 1.0;
        let mut IAS = false;
        let mut IBS = false;
        if !(RFNU2 < TOL) {
            //GO TO 110;
            // DO 100 IS=2,7;
            for _IS in 1..7 {
                ATOL /= RFNU2;
                PP *= RFNU2;
                if !IAS {
                    //GO TO 60;
                    let mut suma = Complex64::zero();
                    //   SUMAR = ZEROR;
                    //   SUMAI = ZEROI;
                    //   DO 40 K=1,KMAX;
                    for k in 0..KMAX {
                        //     let M = L1 + K;
                        suma += p[k] * ALFA[L1 + k];
                        //     SUMAR = SUMAR + PR(K)*ALFA(M);
                        //     SUMAI = SUMAI + PI(K)*ALFA(M);
                        if ap[k] < ATOL {
                            break;
                        } //GO TO 50;
                    }
                    //    40   CONTINUE;
                    //    50   CONTINUE;
                    asum += suma * PP;
                    //   ASUMR = ASUMR + SUMAR*PP;
                    //   ASUMI = ASUMI + SUMAI*PP;
                    if PP < TOL {
                        IAS = true
                    };
                }
                //    60   CONTINUE;
                if !IBS {
                    //GO TO 90;
                    let mut sumb = Complex64::zero();
                    //   SUMBR = ZEROR;
                    //   SUMBI = ZEROI;
                    //   DO 70 K=1,KMAX;
                    for k in 0..KMAX {
                        //     let M = L2 + K;
                        sumb += p[k] * BETA[L2 + k];
                        //     SUMBR = SUMBR + PR(K)*BETA(M);
                        //     SUMBI = SUMBI + PI(K)*BETA(M);
                        if ap[k] < ATOL {
                            break;
                        } //GO TO 80;
                    }
                    //    70   CONTINUE;
                    //    80   CONTINUE;
                    bsum += sumb * PP;
                    //   BSUMR = BSUMR + SUMBR*PP;
                    //   BSUMI = BSUMI + SUMBI*PP;
                    if PP < BTOL {
                        IBS = true;
                    }
                }
                //    90   CONTINUE;
                if IAS && IBS {
                    break;
                } //GO TO 110;
                L1 = L1 + 30;
                L2 = L2 + 30;
            }
            //   100 CONTINUE;
        }
        //   110 CONTINUE;
        asum += 1.0;
        // ASUMR = ASUMR + CONER;
        PP = RFNU * RFN13;
        bsum *= PP;
        // BSUMR = BSUMR*PP;
        // BSUMI = BSUMI*PP;
        //   120 CONTINUE;
        //       RETURN;
        return (phi, arg, zeta1, zeta2, Some(bsum), Some(bsum));
    } else {
        //-----------------------------------------------------------------------;
        //     CABS(W2) > 0.25;
        //-----------------------------------------------------------------------;
        //   130 CONTINUE;
        let mut w = w2.sqrt();
        // CALL ZSQRT(W2R, W2I, WR, WI);
        if w.re < 0.0 {
            w.re = 0.0
        };
        if w.im < 0.0 {
            w.im = 0.0
        };
        let st = Complex64::one() + w;
        // STR = CONER + WR;
        // STI = WI;
        let za = st / zb;
        // CALL ZDIV(STR, STI, ZBR, ZBI, ZAR, ZAI);
        // CALL ZLOG(ZAR, ZAI, ZCR, ZCI, IDUM);
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
        // ZTHR = (ZCR-WR)*1.5;
        // ZTHI = (ZCI-WI)*1.5;
        let zeta1 = zc * order;
        // ZETA1R = ZCR*FNU;
        // ZETA1I = ZCI*FNU;
        let zeta2 = w * order;
        // ZETA2R = WR*FNU;
        // ZETA2I = WI*FNU;
        let AZTH = zth.abs(); //ZABS(ZTHR,ZTHI);
        let mut ANG = THREE_PI_BY_2;
        if !(zth.re >= 0.0 && zth.im < 0.0) {
            //GO TO 140;
            ANG = FRAC_PI_2;
            if zth.re != 0.0 {
                //GO TO 140;
                ANG = (zth.im / zth.re).atan(); //DATAN(ZTHI/ZTHR);
                if zth.re < 0.0 {
                    ANG = ANG + PI
                };
            }
        }
        //   140 CONTINUE;

        let mut PP = AZTH.powf(EX2);
        ANG *= EX2;
        let mut zeta = PP * Complex64::new(ANG.cos(), ANG.sin());
        // ZETAR = PP*DCOS(ANG);
        // ZETAI = PP*DSIN(ANG);
        if zeta.im < 0.0 {
            zeta.im = 0.0
        };
        let arg = zeta * FN23;
        // ARGR = ZETAR*FN23;
        // ARGI = ZETAI*FN23;
        let rtzt = zth / zeta;
        // CALL ZDIV(ZTHR, ZTHI, ZETAR, ZETAI, RTZTR, RTZTI);
        let za = rtzt / w;
        // CALL ZDIV(RTZTR, RTZTI, WR, WI, ZAR, ZAI);
        let tazr = za + za;
        // TZAR = ZAR + ZAR;
        // TZAI = ZAI + ZAI;
        // let st = tazr.sqrt();
        // CALL ZSQRT(TZAR, TZAI, STR, STI);
        let phi = tazr.sqrt() * RFN13;
        // PHIR = STR*RFN13;
        // PHII = STI*RFN13;
        // if (IPMTR == 1) GO TO 120;
        if only_phi_zeta {
            return (phi, arg, zeta1, zeta2, None, None);
        }

        let RAW = 1.0 / AW2.sqrt();
        // let st = ;
        // STR = WR*RAW;
        // STI = -WI*RAW;
        let tfn = w.conj() * RAW * RAW * RFNU;
        // TFNR = STR*RFNU*RAW;
        // TFNI = STI*RFNU*RAW;
        let RAZTH = 1.0 / AZTH;
        // STR = ZTHR*RAZTH;
        // STI = -ZTHI*RAZTH;
        let rzth = zth.conj() * RAZTH * RAZTH * RFNU;
        // RZTHR = STR*RAZTH*RFNU;
        // RZTHI = STI*RAZTH*RFNU;
        let zc = RAZTH * AR[1];
        // ZCR = RZTHR*AR(2);
        // ZCI = RZTHI*AR(2);
        let RAW2 = 1.0 / AW2;
        // STR = W2R*RAW2;
        // STI = -W2I*RAW2;
        let t2 = w2.conj() * RAW2 * RAW2;
        // T2R = STR*RAW2;
        // T2I = STI*RAW2;
        let st = t2 * C[1] + C[2];
        // STR = T2R*C(2) + C(3);
        // STI = T2I*C(2);
        let mut up = vec![Complex64::zero(); 14];
        up[1] = st * tfn;
        // UPR(2) = STR*TFNR - STI*TFNI;
        // UPI(2) = STR*TFNI + STI*TFNR;
        let mut bsum = up[1] + zc;
        // BSUMR = UPR(2) + ZCR;
        // BSUMI = UPI(2) + ZCI;
        let mut asum = Complex64::zero();
        // ASUMR = ZEROR;
        // ASUMI = ZEROI;
        if !(RFNU < TOL) {
            //GO TO 220;
            let mut przth = rzth;
            // PRZTHR = RZTHR;
            // PRZTHI = RZTHI;
            let mut ptfn = tfn;
            // PTFNR = TFNR;
            // PTFNI = TFNI;
            up[0] = Complex64::one();
            // UPR(1) = CONER;
            // UPI(1) = CONEI;
            PP = 1.0;
            let BTOL = TOL * (bsum.re.abs() + bsum.im.abs());
            let mut KS = 0;
            let mut KP1 = 2;
            let mut L = 2; //3;
            let mut IAS = false;
            let mut IBS = false;
            // DO 210 LR=2,12,2;
            for LR in (2..=12).step_by(2) {
                let LRP1 = LR + 1;
                //-----------------------------------------------------------------------;
                //     COMPUTE TWO ADDITIONAL CR, DR, AND UP FOR TWO MORE TERMS IN;
                //     NEXT SUMA AND SUMB;
                //-----------------------------------------------------------------------;
                //   DO 160 K=LR,LRP1;
                let mut cr = vec![Complex64::zero(); 14];
                let mut dr = vec![Complex64::zero(); 14];
                for _K in LR..=LRP1 {
                    KS = KS + 1;
                    KP1 = KP1 + 1;
                    L = L + 1;
                    let mut za = Complex64::new(C[L - 1], 0.0);
                    //     ZAR = C(L);
                    //     ZAI = ZEROI;
                    //     DO 150 J=2,KP1;
                    for _J in 2..=KP1 {
                        L = L + 1;
                        // STR = ZAR*T2R - T2I*ZAI + C(L);
                        // ZAI = ZAR*T2I + ZAI*T2R;
                        // ZAR = STR;
                        za = za * t2 + C[L - 1];
                    }
                    //   150     CONTINUE;
                    ptfn *= tfn;
                    //     STR = PTFNR*TFNR - PTFNI*TFNI;
                    //     PTFNI = PTFNR*TFNI + PTFNI*TFNR;
                    //     PTFNR = STR;
                    up[KP1 - 1] = ptfn * za;
                    //     UPR(KP1) = PTFNR*ZAR - PTFNI*ZAI;
                    //     UPI(KP1) = PTFNI*ZAR + PTFNR*ZAI;
                    cr[KS - 1] = przth * BR[KS];
                    //     CRR(KS) = PRZTHR*BR(KS+1);
                    //     CRI(KS) = PRZTHI*BR(KS+1);
                    przth *= rzth;
                    //     STR = PRZTHR*RZTHR - PRZTHI*RZTHI;
                    //     PRZTHI = PRZTHR*RZTHI + PRZTHI*RZTHR;
                    //     PRZTHR = STR;
                    dr[KS - 1] = przth * AR[KS + 1];
                    //     DRR(KS) = PRZTHR*AR(KS+2);
                    //     DRI(KS) = PRZTHI*AR(KS+2);
                }
                //   160   CONTINUE;
                PP = PP * RFNU2;
                if !IAS {
                    //GO TO 180;
                    let mut suma = up[LRP1 - 1];
                    //   SUMAR = UPR(LRP1);
                    //   SUMAI = UPI(LRP1);
                    let mut JU = LRP1;
                    //   DO 170 JR=1,LR;
                    for JR in 0..LR {
                        JU = JU - 1;
                        suma += cr[JR] * up[JU - 1];
                        //     SUMAR = SUMAR + CRR(JR)*UPR(JU) - CRI(JR)*UPI(JU);
                        //     SUMAI = SUMAI + CRR(JR)*UPI(JU) + CRI(JR)*UPR(JU);
                    }
                    //   170   CONTINUE;
                    asum += suma;
                    //   ASUMR = ASUMR + SUMAR;
                    //   ASUMI = ASUMI + SUMAI;
                    let TEST = suma.re.abs() + suma.im.abs();
                    if PP < TOL && TEST < TOL {
                        IAS = true
                    };
                }
                //   180   CONTINUE;
                if !IBS {
                    //GO TO 200;
                    let mut sumb = up[LR - 2] + up[LRP1 - 1] * zc;
                    //   SUMBR = UPR(LR+2) + UPR(LRP1)*ZCR - UPI(LRP1)*ZCI;
                    //   SUMBI = UPI(LR+2) + UPR(LRP1)*ZCI + UPI(LRP1)*ZCR;
                    let mut JU = LRP1;
                    for JR in 0..LR {
                        //   DO 190 JR=1,LR;
                        JU = JU - 1;
                        sumb += dr[JR] * up[JU - 1];
                        //     SUMBR = SUMBR + DRR(JR)*UPR(JU) - DRI(JR)*UPI(JU);
                        //     SUMBI = SUMBI + DRR(JR)*UPI(JU) + DRI(JR)*UPR(JU);
                    }
                    //   190   CONTINUE;
                    bsum += sumb;
                    //   BSUMR = BSUMR + SUMBR;
                    //   BSUMI = BSUMI + SUMBI;
                    let TEST = sumb.re.abs() + sumb.im.abs();
                    if PP < BTOL && TEST < BTOL {
                        IBS = true
                    };
                    //   200   CONTINUE;
                }
                if IAS && IBS {
                    break;
                }
                //   210 CONTINUE;
            }
        }
        //   220 CONTINUE;
        asum += Complex64::one();
        // ASUMR = ASUMR + CONER;
        bsum = (-bsum * RFN13) / rtzt;
        // STR = -BSUMR*RFN13;
        // STI = -BSUMI*RFN13;
        // CALL ZDIV(STR, STI, RTZTR, RTZTI, BSUMR, BSUMI);
        // GO TO 120;
        return (phi, arg, zeta1, zeta2, Some(bsum), Some(bsum));

        // END;
    }
}
