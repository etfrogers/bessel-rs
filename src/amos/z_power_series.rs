use num::{
    One, Zero,
    complex::{Complex64, ComplexFloat},
};

use super::{BesselError, Scaling, gamma_ln, machine::d1mach, utils::will_z_underflow};

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
    let mut nz = 0;
    let az = z.abs(); //ZABS(ZR,ZI)
    let mut y = vec![Complex64::zero(); n];
    let mut w = vec![Complex64::zero(); 2];
    if az == 0.0 {
        // Not setting zero below, as was set in initialisation of y above
        if order == 0.0 {
            y[1] = Complex64::one();
        }
        return Ok((y, nz));
    }
    let arm = 1.0e+3 * d1mach(1);
    let rtr1 = arm.sqrt();
    let mut crscr = 1.0;
    let mut iflag = false;
    if az < arm {
        nz = n.try_into().unwrap();
        if order == 0.0 {
            nz = nz - 1;
        }
        // Not setting zero below, as was set in initialisation of y above
        if order == 0.0 {
            y[1] = Complex64::one();
        }
        // Not setting zero below, as was set in initialisation of y above
        return Ok((y, nz));
    }
    let hz = 0.5 * z;
    let mut cz = Complex64::zero();
    if !(az <= rtr1) {
        cz = (hz).powu(2);
    }

    let acz = cz.abs();
    let mut nn = n;
    let ck = hz.ln();

    let mut ak1;
    let mut fnup;
    let mut ak;
    let mut st;
    let mut ascle = f64::INFINITY;
    'l20: loop {
        let mut dfnu = order + ((nn - 1) as f64);
        fnup = dfnu + 1.0;
        //-----------------------------------------------------------------------;
        //     UNDERFLOW TEST
        //-----------------------------------------------------------------------;
        ak1 = ck * dfnu;
        ak = gamma_ln(fnup).unwrap();
        ak1.re -= ak;
        if kode == Scaling::Scaled {
            ak1.re -= z.re;
        }
        let mut skip_to_40 = ak1.re > (-elim);
        'l30: loop {
            if !skip_to_40 {
                nz = nz + 1;
                y[nn - 1] = Complex64::zero();
                if acz > dfnu {
                    return Ok((y, -nz));
                }
                nn -= 1;
                if nn == 0 {
                    return Ok((y, -nz));
                }
                continue 'l20;
            } else {
                skip_to_40 = false; // should only skip once until sent back to 'l20
            }
            let mut ss = 0.0;
            if !(ak1.re > (-alim)) {
                iflag = true;
                ss = 1.0 / tol;
                crscr = tol;
                ascle = arm * ss;
            }
            let mut aa = ak1.re.exp();
            if iflag {
                aa *= ss
            };
            let mut coef = aa * Complex64::new(ak1.im.cos(), ak1.im.sin());
            let atol = tol * acz / fnup;
            let il = 2.min(nn);
            for i in 0..il {
                dfnu = order + ((nn - (i + 1)) as f64);
                fnup = dfnu + 1.0;
                let mut s1 = Complex64::one();
                if !(acz < tol * fnup) {
                    ak1 = Complex64::one();
                    ak = fnup + 2.0;
                    let mut s = fnup;
                    aa = 2.0;
                    '_l60: loop {
                        let rs = 1.0 / s;
                        st = ak1 * cz;
                        ak1 = st * rs;
                        s1 += ak1;
                        s = s + ak;
                        ak = ak + 2.0;
                        aa = aa * acz * rs;
                        if !(aa > atol) {
                            break;
                        }
                    }
                }
                let s2 = s1 * coef;
                w[i] = s2;
                if iflag && will_z_underflow(s2, ascle, tol) {
                    continue 'l30;
                }
                let m = nn - i - 1;
                y[m] = s2 * crscr;
                if i != (il - 1) {
                    st = coef / hz;
                    coef = st * dfnu;
                }
            }
            break 'l30;
        }
        break 'l20;
    }
    if nn <= 2 {
        return Ok((y, nz));
    }
    let mut k = nn - 2;
    ak = k as f64;
    let raz = 1.0 / az;
    st = z.conj() * raz;
    let rz = 2.0 * st * raz;
    let mut ib;
    '_l100: loop {
        if !iflag {
            ib = 3;
            '_l110: for _ in ib..=nn {
                y[k - 1] = (ak + order) * (rz * y[k]) + y[k + 1];
                ak = ak - 1.0;
                k = k - 1;
            }
            return Ok((y, nz));
        }
        //-----------------------------------------------------------------------;
        //     RECUR BACKWARD WITH SCALED VALUES;
        //-----------------------------------------------------------------------;
        //   120 CONTINUE;
        //-----------------------------------------------------------------------;
        //     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION ABOVE THE;
        //     UNDERFLOW LIMIT = ASCLE = d1mach(1)*SS*1.0D+3;
        //-----------------------------------------------------------------------;
        let mut s1 = w[0];
        let mut s2 = w[1];
        let mut to_return = true;
        let mut l = 0;
        debug_assert!(nn >= 3);
        for l_inner in 2..nn {
            let ck = s2;
            s2 = s1 + (ak + order) * (rz * ck);
            s1 = ck;
            let ck = s2 * crscr;
            y[k] = ck;
            ak = ak - 1.0;
            if k > 0 {
                k -= 1;
            }
            l = l_inner + 1;
            if ck.abs() > ascle {
                to_return = false;
                break;
            }
        }
        if to_return {
            return Ok((y, nz));
        }
        ib = l + 1;
        if ib > nn {
            return Ok((y, nz));
        };
    }
}
