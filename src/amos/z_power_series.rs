use num::{
    One, Zero,
    complex::{Complex64, ComplexFloat},
};

use super::{BesselError, CONEI, CONER, Scaling, gamma_ln, machine::d1mach, utils::ZUunderflowCHK};

pub fn z_power_series(
    z: Complex64, //ZR, ZI,
    order: f64,   //FNU,
    kode: Scaling,
    n: usize,
    //Ouputs: YR, YI, NZ,
    tol: f64,
    elim: f64,
    alim: f64,
) -> Result<(Vec<Complex64>, isize), BesselError> {
    // ***BEGIN PROLOGUE  z_power_series
    // ***REFER TO  ZBESI,ZBESK
    //
    //     z_power_series COMPUTES THE I BESSEL FUNCTION FOR REAL(Z) >= 0.0 BY
    //     MEANS OF THE POWER SERIES FOR LARGE CABS(Z) IN THE
    //     REGION CABS(Z) <= 2*SQRT(FNU+1). NZ=0 IS A NORMAL RETURN.
    //     NZ > 0 MEANS THAT THE LAST NZ COMPONENTS WERE SET TO ZERO
    //     DUE TO UNDERFLOW. NZ < 0 MEANS UNDERFLOW OCCURRED, BUT THE
    //     CONDITION CABS(Z) <= 2*SQRT(FNU+1) WAS VIOLATED AND THE
    //     COMPUTATION MUST BE COMPLETED IN ANOTHER ROUTINE WITH N=N-ABS(NZ).
    //
    // ***ROUTINES CALLED  gamma_ln,d1mach,ZUunderflowCHK,ZABS,ZDIV,ZLOG,ZMLT
    // ***END PROLOGUE  z_power_series
    //     COMPLEX AK1,CK,COEF,CONE,CRSC,CSCL,CZ,CZERO,HZ,RZ,S1,S2,Y,Z
    //       EXTERNAL ZABS
    //       DOUBLE PRECISION AA, ACZ, AK, AK1I, AK1R, ALIM, ARM, ASCLE, ATOL,
    //      * AZ, CKI, CKR, COEFI, COEFR, CONEI, CONER, CRSCR, CZI, CZR, DFNU,
    //      * ELIM, FNU, FNUP, HZI, HZR, RAZ, RS, RTR1, RZI, RZR, S, SS, STI,
    //      * STR, S1I, S1R, S2I, S2R, TOL, YI, YR, WI, WR, ZEROI, ZEROR, ZI,
    //      * ZR, gamma_ln, d1mach, ZABS
    //       INTEGER I, IB, IDUM, IFLAG, IL, K, KODE, L, M, N, NN, NZ, NW
    //       DIMENSION YR(N), YI(N), WR(2), WI(2)
    //       DATA ZEROR,ZEROI,CONER,CONEI / 0.0, 0.0, 1.0, 0.0 /
    //
    let mut NZ = 0;
    let AZ = z.abs(); //ZABS(ZR,ZI)
    let mut y = vec![Complex64::zero(); n];
    let mut w = vec![Complex64::zero(); 2];
    if AZ == 0.0 {
        // Not setting zero below, as was set in initialisation of y above
        // y[1] = Complex64::new(ZEROR, ZEROI);
        // YR[1] = ZEROR;
        // YI[1] = ZEROI;
        if order == 0.0 {
            y[1] = Complex64::new(CONER, CONEI);
            // YR(1) = CONER;
            // YI(1) = CONEI;
        }
        // Not setting zero below, as was set in initialisation of y above
        // if N == 1 {return Ok((y, NZ));}
        // for I in 2..N
        // {
        //       YR(I) = ZEROR;
        //       YI(I) = ZEROI;
        // }
        return Ok((y, NZ));
    }
    let ARM = 1.0e+3 * d1mach(1);
    let RTR1 = ARM.sqrt();
    let mut CRSCR = 1.0;
    let mut IFLAG = false;
    if AZ < ARM {
        //GO TO 150;
        NZ = n.try_into().unwrap();
        if order == 0.0 {
            NZ = NZ - 1;
        }
        // Not setting zero below, as was set in initialisation of y above
        // y[1] = Complex64::new(ZEROR, ZEROI);
        // YR[1] = ZEROR;
        // YI[1] = ZEROI;
        if order == 0.0 {
            y[1] = Complex64::new(CONER, CONEI);
            // YR(1) = CONER;
            // YI(1) = CONEI;
        }
        // Not setting zero below, as was set in initialisation of y above
        // if N == 1 {return Ok((y, NZ));}
        // for I in 2..N
        // {
        //       YR(I) = ZEROR;
        //       YI(I) = ZEROI;
        // }
        return Ok((y, NZ));
    }
    // let HZR = 0.5*z.re;
    // let HZI = 0.5*z.im;
    let hz = 0.5 * z;
    // let CZR = ZEROR;
    // let CZI = ZEROI;
    let mut cz = Complex64::zero();
    if !(AZ <= RTR1) {
        //CALL ZMLT(HZR, HZI, HZR, HZI, CZR, CZI);
        cz = (hz).powu(2);
    }

    let ACZ = cz.abs(); //ZABS(CZR,CZI);
    let mut NN = n;
    // CALL ZLOG(HZR, HZI, CKR, CKI, IDUM);
    let ck = hz.ln();

    let mut ak1;
    let mut FNUP;
    let mut AK;
    let mut st;
    let mut ASCLE = f64::INFINITY;
    'l20: loop {
        // 20 CONTINUE;
        let mut DFNU = order + ((NN - 1) as f64);
        FNUP = DFNU + 1.0;
        //-----------------------------------------------------------------------;
        //     UNDERFLOW TEST
        //-----------------------------------------------------------------------;
        // let AK1R = CKR*DFNU;
        // let AK1I = CKI*DFNU;
        ak1 = ck * DFNU;
        AK = gamma_ln(FNUP).unwrap();
        ak1.re -= AK;
        if kode == Scaling::Scaled {
            ak1.re -= z.re;
        }
        let skip_to_40 = ak1.re > (-elim);
        //    30 CONTINUE;
        'l30: loop {
            if !skip_to_40 {
                NZ = NZ + 1;
                // YR(NN) = ZEROR;
                // YI(NN) = ZEROI;
                y[NN - 1] = Complex64::zero();
                if ACZ > DFNU {
                    return Ok((y, -NZ));
                }
                NN -= 1;
                if NN == 0 {
                    return Ok((y, -NZ));
                }
                continue 'l20;
            }
            // GO TO 20;
            //    40 CONTINUE;
            let mut SS = 0.0;
            if !(ak1.re > (-alim))
            //GO TO 50;
            {
                IFLAG = true;
                SS = 1.0 / tol;
                CRSCR = tol;
                ASCLE = ARM * SS;
            }
            //    50 CONTINUE;
            let mut AA = ak1.re.exp();
            if IFLAG {
                AA *= SS
            };
            let mut coef = AA * Complex64::new(ak1.im.cos(), ak1.im.sin());
            // let COEFR = AA * ak1.im.cos();
            // let COEFI = AA * ak1.im.sin();
            let ATOL = tol * ACZ / FNUP;
            let IL = 2.min(NN);
            for I in 0..IL {
                // DO 90 I=1,IL;
                DFNU = order + ((NN - (I + 1)) as f64);
                FNUP = DFNU + 1.0;
                let mut s1 = Complex64::one();
                // S1R = CONER;
                // S1I = CONEI;
                if !(ACZ < tol * FNUP) {
                    //GO TO 70;
                    //     AK1R = CONER;
                    //     AK1I = CONEI;
                    ak1 = Complex64::one();
                    AK = FNUP + 2.0;
                    let mut S = FNUP;
                    AA = 2.0;
                    //    60   CONTINUE;
                    'l60: loop {
                        let RS = 1.0 / S;
                        st = ak1 * cz;
                        //   STR = AK1R * CZR - AK1I * CZI;
                        //   STI = AK1R * CZI + AK1I * CZR;
                        // AK1R = STR * RS;
                        // AK1I = STI * RS;
                        ak1 = st * RS;
                        s1 += ak1;
                        // S1R = S1R + AK1R;
                        // S1I = S1I + AK1I;
                        S = S + AK;
                        AK = AK + 2.0;
                        AA = AA * ACZ * RS;
                        if !(AA > ATOL) {
                            break;
                        }
                    }
                }
                // 70   CONTINUE;
                let s2 = s1 * coef;
                w[I] = s2;
                //     WR(I) = S2R;
                //     WI(I) = S2I;
                if IFLAG && ZUunderflowCHK(s2, ASCLE, tol)
                //GO TO 80;
                {
                    continue 'l30;
                    //GO TO 30;
                }
                //    80   CONTINUE;
                let M = NN - I - 1;
                //   YR(M) = S2R*CRSCR;
                //   YI(M) = S2I*CRSCR;
                y[M] = s2 * CRSCR;
                if I != (IL - 1)
                //GO TO 90;
                {
                    // CALL ZDIV(COEFR, COEFI, HZR, HZI, STR, STI);
                    st = coef / hz;
                    // COEFR = STR*DFNU;
                    // COEFI = STI*DFNU;
                    coef = st * DFNU;
                }
            }
            break 'l30;
        }
        break 'l20;
    }
    //    90 CONTINUE;
    if NN <= 2 {
        return Ok((y, NZ));
    }
    let mut K = NN - 2;
    AK = K as f64;
    let RAZ = 1.0 / AZ;
    st = z.conj() * RAZ;
    //     STR = ZR * RAZ;
    //     STI = -ZI * RAZ;
    let rz = 2.0 * st * RAZ;
    //     RZR = (STR + STR) * RAZ;
    //     RZI = (STI + STI) * RAZ;
    let mut IB;
    'l100: loop {
        // moved outside if statment, but believe the logic holds
        if !IFLAG {
            //GO TO 120;
            IB = 3;
            // 100 CONTINUE;
            'l110: for _I in IB..=NN {
                //DO 110 I=IB,NN;
                y[K].re = (AK + order) * (rz.re * y[K + 1].re - rz.im * y[K + 1].im) + y[K + 2].re;
                y[K].im = (AK + order) * (rz.re * y[K + 1].im + rz.im * y[K + 1].re) + y[K + 2].im;
                AK = AK - 1.0;
                K = K - 1;
            }
            //   110 CONTINUE;
            return Ok((y, NZ));
        }
        //-----------------------------------------------------------------------;
        //     RECUR BACKWARD WITH SCALED VALUES;
        //-----------------------------------------------------------------------;
        //   120 CONTINUE;
        //-----------------------------------------------------------------------;
        //     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION ABOVE THE;
        //     UNDERFLOW LIMIT = ASCLE = d1mach(1)*SS*1.0D+3;
        //-----------------------------------------------------------------------;
        let mut s1 = w[1];
        let mut s2 = w[2];
        //   S1R = WR(1);
        //   S1I = WI(1);
        //   S2R = WR(2);
        //   S2I = WI(2);
        let mut to_return = true;
        let mut L = 0;
        for L_inner in 3..=NN {
            // DO 130 L=3,NN;
            // CKR = S2R;
            // CKI = S2I;
            let ck = s2;
            s2 = s1 + (AK + order) * (rz * ck);
            // S2R = S1R + (AK + order) * (RZR * CKR - RZI * CKI);
            // S2I = S1I + (AK + order) * (RZR * CKI + RZI * CKR);
            s1 = ck;
            // S1R = CKR;
            // S1I = CKI;
            let ck = s2 * CRSCR;
            // CKR = S2R * CRSCR;
            // CKI = S2I * CRSCR;
            y[K] = ck;
            // YR(K) = CKR;
            // YI(K) = CKI;
            AK = AK - 1.0;
            K = K - 1;
            L = L_inner;
            if ck.abs() > ASCLE {
                // if (ZABS(CKR, CKI) > ASCLE) {
                to_return = false;
                break;
            } //GO TO 140;
        }
        //   130 CONTINUE;
        if to_return {
            return Ok((y, NZ));
        }
        //   140 CONTINUE;
        IB = L + 1;
        if IB > NN {
            return Ok((y, NZ));
        };
    } //GO TO 100;
    //   150 CONTINUE;
    //       NZ = N;
    //       if (FNU == 0.0) NZ = NZ - 1;
    //   160 CONTINUE;
    //       YR(1) = ZEROR;
    //       YI(1) = ZEROI;
    //       if (FNU == 0.0) {;
    //       YR(1) = CONER;
    //       YI(1) = CONEI;
    //       };
    //       if (N == 1) RETURN;
    //       DO 180 I=2,N;
    //         YR(I) = ZEROR;
    //         YI(I) = ZEROI;
    //   180 CONTINUE;
    //       RETURN;
    //-----------------------------------------------------------------------
    //     RETURN WITH NZ < 0 if CABS(Z*Z/4) > FNU+N-NZ-1 COMPLETE
    //     THE CALCULATION IN CBINU WITH N=N-NZ.abs()
    //-----------------------------------------------------------------------
    //   190 CONTINUE
    //       NZ = -NZ
    //       RETURN
    //       END
}
