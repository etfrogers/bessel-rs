use num::complex::Complex64;

pub fn ZUunderflowCHK(
    y: Complex64, //YR, YI,
    // Outputs: NZ
    ascle: f64,
    tol: f64,
) -> bool {
    // ***BEGIN PROLOGUE  ZUunderflowCHK
    // ***REFER TO z_power_series,ZUOIK,ZUNK1,ZUNK2,ZUNI1,ZUNI2,ZKSCL
    //
    //      Y ENTERS AS A SCALED QUANTITY WHOSE MAGNITUDE IS GREATER THAN
    //      EXP(-ALIM)=ASCLE=1.0E+3*d1mach(1)/TOL. THE TEST IS MADE TO SEE
    //      if THE MAGNITUDE OF THE REAL OR IMAGINARY PART WOULD UNDERFLOW
    //      WHEN Y IS SCALED (BY TOL) TO ITS PROPER VALUE. Y IS ACCEPTED
    //      if THE UNDERFLOW IS AT LEAST ONE PRECISION BELOW THE MAGNITUDE
    //      OF THE LARGEST COMPONENT; OTHERWISE THE PHASE ANGLE DOES NOT HAVE
    //      ABSOLUTE ACCURACY AND AN UNDERFLOW IS ASSUMED.
    //
    // ***ROUTINES CALLED  (NONE)
    // ***END PROLOGUE  ZUunderflowCHK
    //
    //     COMPLEX Y
    // DOUBLE PRECISION ASCLE, SS, ST, TOL, WR, WI, YR, YI
    // INTEGER NZ
    let re_abs = y.re.abs();
    let im_abs = y.im.abs();
    let min_abs_component = re_abs.min(im_abs);
    if min_abs_component > ascle {
        false
    } else {
        let max_abs_component = re_abs.max(im_abs);
        max_abs_component < min_abs_component / tol
    }
}
