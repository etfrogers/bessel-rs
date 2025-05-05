use std::f64::consts::PI;

use num::complex::{Complex64, ComplexFloat};

use super::{
    BesselError::*, BesselResult, MachineConsts, Scaling, c_one, c_zero, c_zeros, utils::RTPI,
};

pub fn z_asymptotic_i(
    z: Complex64,
    order: f64,
    kode: Scaling,
    n: usize,
    machine_consts: &MachineConsts,
) -> BesselResult {
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
    let nz = 0;
    let mut y = c_zeros(n);
    let az = z.abs();
    let rtr1 = machine_consts.underflow_limit.sqrt();
    let il = 2.min(n);
    let dfnu = order + ((n - il) as f64);
    //-----------------------------------------------------------------------;
    //     OVERFLOW TEST;
    //-----------------------------------------------------------------------;
    let raz = 1.0 / az;
    let mut ak1 = (RTPI * z.conj() * raz.powi(2)).sqrt();
    let mut cz = z;
    if kode == Scaling::Scaled {
        cz.re = 0.0;
    }
    if cz.re.abs() > machine_consts.exponent_limit {
        return Err(Overflow);
    }
    let dnu2 = dfnu + dfnu;
    let mut koded = true;
    if !((cz.re.abs() > machine_consts.approximation_limit) && (n > 2)) {
        koded = false;
        ak1 *= cz.exp();
    }
    let mut fdn = 0.0;
    if dnu2 > rtr1 {
        fdn = dnu2 * dnu2
    };
    let ez = z * 8.0;
    //-----------------------------------------------------------------------;
    //     WHEN Z IS IMAGINARY, THE ERROR TEST MUST BE MADE RELATIVE TO THE;
    //     FIRST RECIPROCAL POWER SINCE THIS IS THE LEADING TERM OF THE;
    //     EXPANSION FOR THE IMAGINARY PART.;
    //-----------------------------------------------------------------------;
    let aez = 8.0 * az;
    let s = machine_consts.abs_error_tolerance / aez;
    let jl = (machine_consts.asymptotic_z_limit * 2.0) as i32 + 2;
    let mut p1 = c_zero();
    if z.im != 0.0 {
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
        if inu % 2 != 0 {
            p1 = -p1;
        }
    }
    //    30 CONTINUE;
    // DO 70 K=1,IL;
    '_l70: for k in 0..il {
        let mut sqk = fdn - 1.0;
        let atol = s * (sqk).abs();
        let mut sign = 1.0;
        let mut cs1 = c_one();
        //   CS1R = CONER;
        //   CS1I = CONEI;
        let mut cs2 = c_one();
        //   CS2R = CONER;
        //   CS2I = CONEI;
        let mut ck = c_one();
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
            ck *= sqk / dk;
            cs2 += ck;
            sign = -sign;
            cs1 += ck * sign;
            dk += ez;
            aa *= sqk.abs() / bb;
            bb += aez;
            ak += 8.0;
            sqk -= ak;
            if aa <= atol {
                converged = true;
                break 'l40;
            }
        }
        if !converged {
            return Err(DidNotConverge);
        }
        let mut s2 = cs1;
        if z.re * 2.0 < machine_consts.exponent_limit {
            s2 += (-z * 2.0).exp() * p1 * cs2;
        }
        fdn = fdn + 8.0 * dfnu + 4.0;
        p1 = -p1;
        y[n - il + k] = s2 * ak1;
    }
    if n <= 2 {
        return Ok((y, nz));
    }

    // let mut k = n - 2;
    // let mut ak = k as f64;
    let rz = z.conj() * 2.0 * raz.powi(2);
    // let ib = 3;
    for k in (1..=n - 2).rev() {
        y[k - 1] = ((k as f64) + order) * (rz * y[k]) + y[k + 1];
        // ak -= 1.0;
        // if k > 0 {
        //     k -= 1;
        // }
    }
    if koded {
        let ck = cz.exp();
        for yi in y.iter_mut() {
            *yi *= ck;
        }
    }
    Ok((y, nz))
}
