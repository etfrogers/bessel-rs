use std::f64::consts::PI;

use num::{
    One, Zero,
    complex::{Complex64, ComplexFloat},
};

use super::{BesselError, BesselError::*, Scaling, machine::d1mach, utils::RTPI};

pub fn z_asymptotic_i(
    //ZR, ZI, FNU,
    z: Complex64,
    order: f64,
    kode: Scaling,
    n: usize,
    //YR, YI, NZ,
    rl: f64,
    tol: f64,
    elim: f64,
    alim: f64,
) -> Result<(Vec<Complex64>, usize), BesselError> {
    // ***BEGIN PROLOGUE  ZASYI
    // ***REFER TO  ZBESI,ZBESK
    //
    //     ZASYI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z) >= 0.0 BY
    //     MEANS OF THE ASYMPTOTIC EXPANSION FOR LARGE CABS(Z) IN THE
    //     REGION CABS(Z) > MAX(RL,FNU*FNU/2). NZ=0 IS A NORMAL RETURN.
    //     NZ < 0 INDICATES AN OVERFLOW ON KODE=1.
    //
    // ***ROUTINES CALLED  d1mach,ZABS,ZDIV,ZEXP,ZMLT,ZSQRT
    // ***END PROLOGUE  ZASYI
    //     COMPLEX AK1,CK,CONE,CS1,CS2,CZ,CZERO,DK,EZ,P1,RZ,S2,Y,Z
    //       EXTERNAL ZABS
    //       DOUBLE PRECISION AA, AEZ, AK, AK1I, AK1R, ALIM, ARG, ARM, ATOL,
    //      * AZ, BB, BK, CKI, CKR, CONEI, CONER, CS1I, CS1R, CS2I, CS2R, CZI,
    //      * CZR, DFNU, DKI, DKR, DNU2, ELIM, EZI, EZR, FDN, FNU, PI, P1I,
    //      * P1R, RAZ, RL, RTPI, RTR1, RZI, RZR, S, SGN, SQK, STI, STR, S2I,
    //      * S2R, TOL, TZI, TZR, YI, YR, ZEROI, ZEROR, ZI, ZR, d1mach, ZABS
    //       INTEGER I, IB, IL, INU, J, JL, K, KODE, KODED, M, N, NN, NZ
    //       DIMENSION YR(N), YI(N)
    // DATA PI, RTPI  /3.14159265358979324 , 0.159154943091895336 /
    // DATA ZEROR,ZEROI,CONER,CONEI / 0.0, 0.0, 1.0, 0.0 /
    //
    let nz = 0;
    let mut y = vec![Complex64::zero(); n];
    let az = z.abs();
    let arm = 1.0e+3 * d1mach(1);
    let rtr1 = arm.sqrt();
    let il = 2.min(n);
    let dfnu = order + ((n - il) as f64);
    //-----------------------------------------------------------------------;
    //     OVERFLOW TEST;
    //-----------------------------------------------------------------------;
    let raz = 1.0 / az;
    let mut st = raz * z.conj();
    // STR = ZR*RAZ;
    // STI = -ZI*RAZ;
    let mut ak1 = RTPI * st * raz;
    // AK1R = RTPI*STR*RAZ;
    // AK1I = RTPI*STI*RAZ;
    // CALL ZSQRT(AK1R, AK1I, AK1R, AK1I);
    ak1 = ak1.sqrt();
    // CZR = ZR;
    // CZI = ZI;
    let mut cz = z;
    if kode == Scaling::Scaled
    //GO TO 10;
    {
        cz.re = 0.0;
        // CZR = ZEROR;
        // CZI = ZI;
    }
    //    10 CONTINUE;
    if cz.re.abs() > elim {
        return Err(Overflow);
    }
    let dnu2 = dfnu + dfnu;
    let mut koded = true;
    if !((cz.re.abs() > alim) && (n > 2)) {
        // GO TO 20;
        koded = false;
        st = cz.exp();
        // CALL ZEXP(CZR, CZI, STR, STI);
        ak1 *= st;
        // CALL ZMLT(AK1R, AK1I, STR, STI, AK1R, AK1I);
    }
    //    20 CONTINUE;
    let mut fdn = 0.0;
    if dnu2 > rtr1 {
        fdn = dnu2 * dnu2
    };
    let ez = z * 8.0;
    // EZR = ZR*8.0;
    // EZI = ZI*8.0;
    //-----------------------------------------------------------------------;
    //     WHEN Z IS IMAGINARY, THE ERROR TEST MUST BE MADE RELATIVE TO THE;
    //     FIRST RECIPROCAL POWER SINCE THIS IS THE LEADING TERM OF THE;
    //     EXPANSION FOR THE IMAGINARY PART.;
    //-----------------------------------------------------------------------;
    let aez = 8.0 * az;
    let s = tol / aez;
    let jl = (rl * 2.0) as i32 + 2;
    // JL = INT(SNGL(RL+RL)) + 2;
    // P1R = ZEROR;
    // P1I = ZEROI;
    let mut p1 = Complex64::zero();
    if !(z.im == 0.0) {
        //GO TO 30;
        //-----------------------------------------------------------------------;
        //     CALCULATE EXP(PI*(0.5+FNU+N-IL)*I) TO MINIMIZE LOSSES OF;
        //     SIGNIFICANCE WHEN FNU OR N IS LARGE;
        //-----------------------------------------------------------------------;
        let mut inu = order as usize; //INT(SNGL(FNU));
        let arg = (order - (inu as f64)) * PI;
        inu = inu + n - il;
        let ak = -arg.sin();
        let mut bk = arg.cos();
        if z.im < 0.0 {
            bk = -bk
        };
        p1 = Complex64::new(ak, bk);
        // P1R = AK;
        // P1I = BK;
        if !(inu % 2) == 0 {
            // GO TO 30;
            p1 = -p1;
            // P1R = -P1R;
            // P1I = -P1I;
        }
    }
    //    30 CONTINUE;
    // DO 70 K=1,IL;
    '_l70: for k in 0..il {
        let mut sqk = fdn - 1.0;
        let atol = s * (sqk).abs();
        let mut sign = 1.0;
        let mut cs1 = Complex64::one();
        //   CS1R = CONER;
        //   CS1I = CONEI;
        let mut cs2 = Complex64::one();
        //   CS2R = CONER;
        //   CS2I = CONEI;
        let mut ck = Complex64::one();
        //   CKR = CONER;
        //   CKI = CONEI;
        let mut ak = 0.0;
        let mut aa = 1.0;
        let mut bb = aez;
        let mut dk = ez;
        //   DKR = EZR;
        //   DKI = EZI;
        //   DO 40 J=1,JL;
        let mut converged = false;
        'l40: for _ in 0..jl {
            st = ck / dk;
            //     CALL ZDIV(CKR, CKI, DKR, DKI, STR, STI);
            ck = st * sqk;
            //     CKR = STR*SQK;
            //     CKI = STI*SQK;
            cs2 += ck;
            //     CS2R = CS2R + CKR;
            //     CS2I = CS2I + CKI;
            sign = -sign;
            cs1 += ck * sign;
            //     CS1R = CS1R + CKR*SGN;
            //     CS1I = CS1I + CKI*SGN;
            dk += ez;
            //     DKR = DKR + EZR;
            //     DKI = DKI + EZI;
            aa = aa * sqk.abs() / bb;
            bb = bb + aez;
            ak = ak + 8.0;
            sqk = sqk - ak;
            if aa <= atol {
                //GO TO 50;
                converged = true;
                break 'l40;
            }
        }
        //    40   CONTINUE;
        if !converged {
            return Err(DidNotConverge);
        }
        //   GO TO 110;
        //    50   CONTINUE;
        let mut s2 = cs1;
        //   S2R = CS1R;
        //   S2I = CS1I;
        if !(z.re * 2.0 >= elim) {
            //GO TO 60;
            let tz = z * 2.0;
            //   TZR = ZR + ZR;
            //   TZI = ZI + ZI;
            st = (-tz).exp();
            //   CALL ZEXP(-TZR, -TZI, STR, STI);
            st *= p1;
            //   CALL ZMLT(STR, STI, P1R, P1I, STR, STI);
            st *= cs2;
            //   CALL ZMLT(STR, STI, CS2R, CS2I, STR, STI);
            s2 += st;
            //   S2R = S2R + STR;
            //   S2I = S2I + STI;
        }
        //    60   CONTINUE;
        fdn = fdn + 8.0 * dfnu + 4.0;
        p1 = -p1;
        //   P1R = -P1R;
        //   P1I = -P1I;
        let m = n - il + k;
        y[m] = s2 * ak1;
        //   YR(M) = S2R*AK1R - S2I*AK1I;
        //   YI(M) = S2R*AK1I + S2I*AK1R;
    }
    //    70 CONTINUE;
    if n <= 2 {
        return Ok((y, nz));
    }

    // NOT TESTED by ZBINU. May be used by ZACAI.
    let nn = n;
    let mut k = nn - 2;
    let mut ak = k as f64;
    st = z.conj() * raz;
    // STR = ZR*RAZ;
    // STI = -ZI*RAZ;
    let rz = st * 2.0 * raz;
    // RZR = (STR+STR)*RAZ;
    // RZI = (STI+STI)*RAZ;
    let ib = 3;
    // DO 80 I=IB,NN;
    '_l80: for _ in ib..nn {
        y[k] = (ak + order) * (rz * y[k + 1]) + y[k + 2];
        //   YR(K) = (AK+FNU)*(RZR*YR(K+1)-RZI*YI(K+1)) + YR(K+2);
        //   YI(K) = (AK+FNU)*(RZR*YI(K+1)+RZI*YR(K+1)) + YI(K+2);
        ak = ak - 1.0;
        k = k - 1;
    }
    //    80 CONTINUE;
    if !koded {
        return Ok((y, nz));
    }
    let ck = cz.exp();
    // CALL ZEXP(CZR, CZI, CKR, CKI);
    // DO 90 I=1,NN;
    '_l90: for i in 0..nn {
        //   STR = YR(I)*CKR - YI(I)*CKI;
        //   YI(I) = YR(I)*CKI + YI(I)*CKR;
        //   YR(I) = STR;
        y[i] *= ck;
    }
    //    90 CONTINUE;
    return Ok((y, nz));
    //   100 CONTINUE;
    //       NZ = -1;
    //       RETURN;
    //   110 CONTINUE;
    //       NZ=-2;
    //       RETURN;
    //       END;
}
