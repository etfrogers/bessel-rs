use num::complex::{Complex64, ComplexFloat};

use super::{
    BesselResult, MachineConsts, Scaling, c_one, c_zero, c_zeros, gamma_ln, utils::will_z_underflow,
};

pub fn z_power_series(
    z: Complex64,
    order: f64,
    kode: Scaling,
    n: usize,
    machine_consts: &MachineConsts,
) -> BesselResult<isize> {
    // ***BEGIN PROLOGUE  z_power_series - was ZSERI
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
    let mut y = c_zeros(n);
    let mut w = [c_zero(); 2];
    if az == 0.0 {
        // Not setting zero below, as was set in initialisation of y above
        if order == 0.0 {
            y[1] = c_one();
        }
        return Ok((y, nz));
    }
    let rtr1 = machine_consts.arm.sqrt();
    let mut crscr = 1.0;
    let mut underflow_would_occur = false;
    if az < machine_consts.arm {
        nz = n.try_into().unwrap();
        if order == 0.0 {
            nz -= 1;
        }
        // Not setting zero below, as was set in initialisation of y above
        if order == 0.0 {
            y[1] = c_one();
        }
        // Not setting zero below, as was set in initialisation of y above
        return Ok((y, nz));
    }
    let hz = 0.5 * z;
    let mut cz = c_zero();
    if az > rtr1 {
        cz = (hz).powu(2);
    }

    let acz = cz.abs();
    let mut nn = n;
    let ck = hz.ln();

    let mut ak1;
    let mut fnup;
    let mut ak;
    'l20: loop {
        let mut dfnu = order + ((nn - 1) as f64);
        fnup = dfnu + 1.0;
        //-----------------------------------------------------------------------;
        //     UNDERFLOW TEST
        //     recur down (setting y to zero) from N until underflow no longer found, then skip to 40
        //-----------------------------------------------------------------------;
        ak1 = ck * dfnu;
        ak = gamma_ln(fnup).unwrap();
        ak1.re -= ak;
        if kode == Scaling::Scaled {
            ak1.re -= z.re;
        }
        let mut skip_to_40 = ak1.re > (-machine_consts.elim);
        'l30: loop {
            if !skip_to_40 {
                nz += 1;
                y[nn - 1] = c_zero();
                if acz > dfnu {
                    return Ok((y, -nz));
                    // Feels like this should return an error, but the amos code effectively
                    // ignores the negative return value anyway. Treating it like an error
                    // changes the behavior of the code
                    // if nz == 2 {
                    //     return Err(DidNotConverge);
                    // } else {
                    //     return Err(Overflow);
                    // };
                }
                nn -= 1;
                if nn == 0 {
                    return Ok((y, nz));
                }
                continue 'l20;
            } else {
                skip_to_40 = false; // should only skip once until sent back to 'l20
            }
            if ak1.re <= (-machine_consts.alim) {
                underflow_would_occur = true;
                crscr = machine_consts.tol;
            }
            let mut aa = ak1.re.exp();
            if underflow_would_occur {
                aa *= machine_consts.rtol
            };
            let mut coef = aa * Complex64::new(ak1.im.cos(), ak1.im.sin());
            let atol = machine_consts.tol * acz / fnup;
            let il = 2.min(nn);
            for i in 0..il {
                dfnu = order + ((nn - (i + 1)) as f64);
                fnup = dfnu + 1.0;
                let mut s1 = c_one();
                if acz >= machine_consts.tol * fnup {
                    ak1 = c_one();
                    ak = fnup + 2.0;
                    let mut s = fnup;
                    aa = 2.0;
                    '_l60: loop {
                        let rs = 1.0 / s;
                        ak1 = ak1 * cz * rs;
                        s1 += ak1;
                        s += ak;
                        ak += 2.0;
                        aa = aa * acz * rs;
                        if aa <= atol {
                            break;
                        }
                    }
                }
                let s2 = s1 * coef;
                w[i] = s2;
                if underflow_would_occur
                    && will_z_underflow(s2, machine_consts.ascle, machine_consts.tol)
                {
                    continue 'l30;
                }
                let m = nn - i - 1;
                y[m] = s2 * crscr;
                if i != (il - 1) {
                    coef = coef.fdiv(hz) * dfnu; //normal div underflows, fdiv more accurate
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
    let rz = 2.0 * z.conj() / (az.powi(2));
    let ib = if underflow_would_occur {
        //-----------------------------------------------------------------------;
        //     RECUR BACKWARD WITH SCALED VALUES;
        //-----------------------------------------------------------------------;
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
            y[k - 1] = ck;
            ak -= 1.0;
            if k > 0 {
                k -= 1;
            }
            l = l_inner + 1;
            if ck.abs() > machine_consts.ascle {
                to_return = false;
                break;
            }
        }
        if to_return {
            return Ok((y, nz));
        }
        if l + 1 > nn {
            return Ok((y, nz));
        };
        l + 1
    } else {
        3
    };
    for _ in ib..=nn {
        y[k - 1] = (ak + order) * (rz * y[k]) + y[k + 1];
        ak -= 1.0;
        k -= 1;
    }
    return Ok((y, nz));
}
