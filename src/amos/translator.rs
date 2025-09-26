#![allow(non_snake_case, clippy::excessive_precision)]
use super::{
    BesselError, BesselResult, HankelKind, IKType, Scaling, c_one, c_zero, c_zeros, gamma_ln,
    i_power_series,
    overflow_checks::{check_underflow_uniform_asymp_params, underflow_add_i_k, zunik},
    utils::{AIC, TWO_THIRDS, is_sigificance_lost, will_underflow},
};
use crate::amos::{
    BesselError::*,
    MACHINE_CONSTANTS, RotationDirection,
    asymptotic_i::asymptotic_i,
    max_abs_component,
    overflow_checks::{Overflow, zunhj},
    utils::{imaginary_dominant, sanitise_inputs},
};
use num::{
    Integer, Zero,
    complex::{Complex64, ComplexFloat},
    pow::Pow,
};
use std::{
    cmp::min,
    f64::consts::{FRAC_2_PI, FRAC_PI_2, PI},
};

/// complex_bessel_h computes the H-bessel functions (Hankel functions) of a complex argument
///
/// ON KODE=1, ZBESH COMPUTES AN N MEMBER SEQUENCE OF COMPLEX
/// HANKEL (BESSEL) FUNCTIONS CY(J)=H(M,FNU+J-1,Z) FOR KINDS M=1
/// OR 2, REAL, NONNEGATIVE ORDERS FNU+J-1, J=1,...,N, AND COMPLEX
/// Z != CMPLX(0.0,0.0) IN THE CUT PLANE -PI < ARG(Z) <= PI.
/// ON KODE=2, ZBESH RETURNS THE SCALED HANKEL FUNCTIONS
///
/// CY(I)=EXP(-MM*Z*I)*H(M,FNU+J-1,Z)       MM=3-2*M,   I**2=-1.
///
/// WHICH REMOVES THE EXPONENTIAL BEHAVIOR IN BOTH THE UPPER AND
/// LOWER HALF PLANES. DEFINITIONS AND NOTATION ARE FOUND IN THE
/// NBS HANDBOOK OF MATHEMATICAL FUNCTIONS (REF. 1).
///
/// INPUT      ZR,ZI,FNU ARE DOUBLE PRECISION
///   ZR,ZI  - Z=CMPLX(ZR,ZI), Z != CMPLX(0.0,0.0),
///            -PT < ARG(Z) <= PI
///   FNU    - ORDER OF INITIAL H FUNCTION, FNU >= 0.0
///   KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
///            KODE= 1  RETURNS
///                     CY(J)=H(M,FNU+J-1,Z),   J=1,...,N
///                = 2  RETURNS
///                     CY(J)=H(M,FNU+J-1,Z)*EXP(-I*Z*(3-2M))
///                          J=1,...,N  ,  I**2=-1
///   M      - KIND OF HANKEL FUNCTION, M=1 OR 2
///   N      - NUMBER OF MEMBERS IN THE SEQUENCE, N >= 1
///
/// OUTPUT     CYR,CYI ARE DOUBLE PRECISION
///   CYR,CYI- DOUBLE PRECISION VECTORS WHOSE FIRST N COMPONENTS
///            CONTAIN REAL AND IMAGINARY PARTS FOR THE SEQUENCE
///            CY(J)=H(M,FNU+J-1,Z)  OR
///            CY(J)=H(M,FNU+J-1,Z)*EXP(-I*Z*(3-2M))  J=1,...,N
///            DEPENDING ON KODE, I**2=-1.
///   NZ     - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW,
///            NZ= 0   , NORMAL RETURN
///            NZ > 0 , FIRST NZ COMPONENTS OF CY SET TO ZERO DUE
///                      TO UNDERFLOW, CY(J)=CMPLX(0.0,0.0)
///                      J=1,...,NZ WHEN Y > 0.0 AND M=1 OR
///                      Y < 0.0 AND M=2. FOR THE COMPLMENTARY
///                      HALF PLANES, NZ STATES ONLY THE NUMBER
///                      OF UNDERFLOWS.
///
/// ***LONG DESCRIPTION
///
/// THE COMPUTATION IS CARRIED OUT BY THE RELATION
///
/// H(M,FNU,Z)=(1/MP)*EXP(-MP*FNU)*K(FNU,Z*EXP(-MP))
///     MP=MM*FRAC_PI_2*I,  MM=3-2*M,  FRAC_PI_2=PI/2,  I**2=-1
///
/// FOR M=1 OR 2 WHERE THE K BESSEL FUNCTION IS COMPUTED FOR THE
/// RIGHT HALF PLANE RE(Z) >= 0.0. THE K FUNCTION IS CONTINUED
/// TO THE LEFT HALF PLANE BY THE RELATION
///
/// K(FNU,Z*EXP(MP)) = EXP(-MP*FNU)*K(FNU,Z)-MP*I(FNU,Z)
/// MP=MR*PI*I, MR=+1 OR -1, RE(Z) > 0, I**2=-1
///
/// WHERE I(FNU,Z) IS THE I BESSEL FUNCTION.
///
/// EXPONENTIAL DECAY OF H(M,FNU,Z) OCCURS IN THE UPPER HALF Z
/// PLANE FOR M=1 AND THE LOWER HALF Z PLANE FOR M=2.  EXPONENTIAL
/// GROWTH OCCURS IN THE COMPLEMENTARY HALF PLANES.  SCALING
/// BY EXP(-MM*Z*I) REMOVES THE EXPONENTIAL BEHAVIOR IN THE
/// WHOLE Z PLANE FOR Z TO INFINITY.
///
/// FOR NEGATIVE ORDERS,THE FORMULAE
///
///   H(1,-FNU,Z) = H(1,FNU,Z)*CEXP( PI*FNU*I)
///   H(2,-FNU,Z) = H(2,FNU,Z)*CEXP(-PI*FNU*I)
///                 I**2=-1
///
/// CAN BE USED.
///
/// IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
/// MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
/// LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
/// CONSEQUENTLY, if EITHER ONE EXCEEDS U1=SQRT(0.5/UR), THEN
/// LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
/// IERR=3 IS TRIGGERED WHERE UR=DMAX1(d1mach(4),1.0e-18) IS
/// DOUBLE PRECISION UNIT ROUNDOFF LIMITED TO 18 DIGITS PRECISION.
/// if EITHER IS LARGER THAN U2=0.5/UR, THEN ALL SIGNIFICANCE IS
/// LOST AND IERR=4. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS
/// MUST BE FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE
/// INTEGER, U3=i1mach(9). THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS
/// RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, AND U3
/// ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION
/// ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE PRECISION
/// ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN
/// THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT ONE CAN EXPECT
/// TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS
/// IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
/// SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.
//
/// THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
/// BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
/// ROUNDOFF,1.0e-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
/// SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
/// ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
/// ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
/// CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
/// HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
/// ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
/// SEVERAL ORDERS OF MAGNITUDE. if ONE COMPONENT IS 10**K LARGER
/// THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
/// 0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
/// THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
/// COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
/// BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
/// COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
/// MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
/// THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
/// OR -PI/2+P.
///
/// ***REFERENCES
///   HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
///         AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
///         COMMERCE, 1955.
///
///   COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
///         BY D. E. AMOS, SAND83-0083, MAY, 1983.
///
///   COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
///         AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
///
///   A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
///         ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
///         1018, MAY, 1985
///
///   A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
///         ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, ACM
///         TRANS. MATH. SOFTWARE, VOL. 12, NO. 3, SEPTEMBER 1986,
///         PP 265-273.
///
/// Original metadata:
/// - Name:  ZBESH
/// - Date written:   830501   (YYMMDD)
/// - Revision date:  890801, 930101   (YYMMDD)
/// - Keywords:  h-bessel functions,bessel functions of complex argument,
///   bessel functions of third kind, hankel functions
/// - Author:  Amos, Donald E., Sandia National Laboratories
pub fn complex_bessel_h(
    z: Complex64,
    order: f64,
    scaling: Scaling,
    hankel_kind: HankelKind,
    n: usize,
) -> BesselResult {
    sanitise_inputs(z, order, n, true)?;
    let mut nz = 0;

    let modified_order = order + ((n - 1) as f64);

    let rotation = hankel_kind.get_rotation();
    let rotation_f64: f64 = rotation.into();
    let mut zn = -Complex64::I * rotation_f64 * z;
    //-----------------------------------------------------------------------
    //     TEST FOR PROPER RANGE
    //-----------------------------------------------------------------------
    let abs_z = z.abs();
    let partial_loss_of_significance = is_sigificance_lost(abs_z, modified_order, false)?;
    //-----------------------------------------------------------------------
    //     OVERFLOW TEST ON THE LAST MEMBER OF THE SEQUENCE
    //-----------------------------------------------------------------------
    if abs_z < MACHINE_CONSTANTS.underflow_limit {
        return Err(Overflow);
    }
    let (mut cy, NZ) = if order < MACHINE_CONSTANTS.asymptotic_order_limit {
        if modified_order > 1.0 {
            if modified_order > 2.0 {
                let mut cy = c_zeros(n);
                let n_underflow = check_underflow_uniform_asymp_params(
                    zn,
                    order,
                    scaling,
                    IKType::K,
                    n,
                    &mut cy,
                )?;

                nz += n_underflow;

                // Here nn=n or nn=0 since n_undeflow=(0 or nn) on return from
                // check_underflow_uniform_asymp_params (for ik_type = k)
                //
                // if nuf=nn, then cy[i]=c_zero() for all i
                if n == n_underflow {
                    return if zn.re < 0.0 {
                        Err(Overflow)
                    } else if partial_loss_of_significance {
                        Err(BesselError::PartialLossOfSignificance { y: cy, nz })
                    } else {
                        Ok((cy, nz))
                    };
                }
            }
            if abs_z <= MACHINE_CONSTANTS.abs_error_tolerance
                && -modified_order * (0.5 * abs_z).ln() > MACHINE_CONSTANTS.exponent_limit
            {
                return Err(Overflow);
            }
        }
        if !((zn.re < 0.0) || (zn.re == 0.0 && zn.im < 0.0 && hankel_kind == HankelKind::Second)) {
            //-----------------------------------------------------------------------
            //     RIGHT HALF PLANE COMPUTATION, XN >= 0. && (XN != 0. ||
            //     YN >= 0. || M=1)
            //-----------------------------------------------------------------------
            k_right_half_plane(zn, order, scaling, n)?
        } else {
            //-----------------------------------------------------------------------
            //     LEFT HALF PLANE COMPUTATION
            //-----------------------------------------------------------------------
            analytic_continuation(zn, order, scaling, -rotation, n)?
        }
    } else {
        //-----------------------------------------------------------------------
        //     UNIFORM ASYMPTOTIC EXPANSIONS FOR FNU > FNUL
        //-----------------------------------------------------------------------
        let mut asymptotic_rotation = RotationDirection::None;
        if !((zn.re >= 0.0) && (zn.re != 0.0 || zn.im >= 0.0 || hankel_kind != HankelKind::Second))
        {
            asymptotic_rotation = -rotation;
            if !(zn.re != 0.0 || zn.im >= 0.0) {
                zn = -zn;
            }
        }
        let (cy, NW) = ZBUNK(zn, order, scaling, asymptotic_rotation, n)?;
        nz += NW;
        (cy, nz)
    };
    //-----------------------------------------------------------------------
    //     H(M,FNU,Z) = -FMM*(I/FRAC_PI_2)*(ZT**FNU)*K(FNU,-Z*ZT)
    //
    //     ZT=EXP(-FMM*FRAC_PI_2*I) = CMPLX(0.0,-FMM), FMM=3-2*M, M=1,2
    //-----------------------------------------------------------------------
    let sign = -FRAC_PI_2 * rotation.signum();
    //-----------------------------------------------------------------------
    //     CALCULATE EXP(FNU*FRAC_PI_2*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
    //     WHEN FNU IS LARGE
    //-----------------------------------------------------------------------
    let int_order = order as i64;
    let half_int_order = int_order / 2;
    let int_remain = int_order - 2 * half_int_order;
    let arg = (order - ((int_order - int_remain) as f64)) * sign;
    let mut csgn = (1.0 / sign) * Complex64::I * Complex64::cis(arg);
    if half_int_order % 2 != 0 {
        csgn = -csgn;
    }
    for element in cy.iter_mut() {
        let ATOL = if max_abs_component(*element) < MACHINE_CONSTANTS.absolute_approximation_limit {
            *element *= MACHINE_CONSTANTS.rtol;
            MACHINE_CONSTANTS.abs_error_tolerance
        } else {
            1.0
        };
        *element *= csgn * ATOL;
        csgn *= Complex64::I * -rotation_f64;
    }
    if partial_loss_of_significance {
        Err(BesselError::PartialLossOfSignificance { y: cy, nz: NZ })
    } else {
        Ok((cy, NZ))
    }
}

pub fn complex_bessel_i(z: Complex64, order: f64, scaling: Scaling, n: usize) -> BesselResult {
    // ***BEGIN PROLOGUE  ZBESI
    // ***DATE WRITTEN   830501   (YYMMDD)
    // ***REVISION DATE  890801, 930101   (YYMMDD)
    // ***CATEGORY NO.  B5K
    // ***KEYWORDS  I-BESSEL FUNCTION,COMPLEX BESSEL FUNCTION,
    //             MODIFIED BESSEL FUNCTION OF THE FIRST KIND
    // ***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
    // ***PURPOSE  TO COMPUTE I-BESSEL FUNCTIONS OF COMPLEX ARGUMENT
    // ***DESCRIPTION
    //
    //                    ***A DOUBLE PRECISION ROUTINE***
    //         ON KODE=1, ZBESI COMPUTES AN N MEMBER SEQUENCE OF COMPLEX
    //         BESSEL FUNCTIONS CY(J)=I(FNU+J-1,Z) FOR REAL, NONNEGATIVE
    //         ORDERS FNU+J-1, J=1,...,N AND COMPLEX Z IN THE CUT PLANE
    //         -PI < ARG(Z) <= PI. ON KODE=2, ZBESI RETURNS THE SCALED
    //         FUNCTIONS
    //
    //         CY(J)=EXP(-ABS(X))*I(FNU+J-1,Z)   J = 1,...,N , X=REAL(Z)
    //
    //         WITH THE EXPONENTIAL GROWTH REMOVED IN BOTH THE LEFT AND
    //         RIGHT HALF PLANES FOR Z TO INFINITY. DEFINITIONS AND NOTATION
    //         ARE FOUND IN THE NBS HANDBOOK OF MATHEMATICAL FUNCTIONS
    //         (REF. 1).
    //
    //         INPUT      ZR,ZI,FNU ARE DOUBLE PRECISION
    //           ZR,ZI  - Z=CMPLX(ZR,ZI),  -PI < ARG(Z) <= PI
    //           FNU    - ORDER OF INITIAL I FUNCTION, FNU >= 0.0
    //           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
    //                    KODE= 1  RETURNS
    //                             CY(J)=I(FNU+J-1,Z), J=1,...,N
    //                        = 2  RETURNS
    //                             CY(J)=I(FNU+J-1,Z)*EXP(-ABS(X)), J=1,...,N
    //           N      - NUMBER OF MEMBERS OF THE SEQUENCE, N >= 1
    //
    //         OUTPUT     CYR,CYI ARE DOUBLE PRECISION
    //           CYR,CYI- DOUBLE PRECISION VECTORS WHOSE FIRST N COMPONENTS
    //                    CONTAIN REAL AND IMAGINARY PARTS FOR THE SEQUENCE
    //                    CY(J)=I(FNU+J-1,Z)  OR
    //                    CY(J)=I(FNU+J-1,Z)*EXP(-ABS(X))  J=1,...,N
    //                    DEPENDING ON KODE, X=REAL(Z)
    //           NZ     - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW,
    //                    NZ= 0   , NORMAL RETURN
    //                    NZ > 0 , LAST NZ COMPONENTS OF CY SET TO ZERO
    //                              TO UNDERFLOW, CY(J)=CMPLX(0.0,0.0)
    //                              J = N-NZ+1,...,N
    //           IERR   - ERROR FLAG
    //                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
    //                    IERR=1, INPUT ERROR   - NO COMPUTATION
    //                    IERR=2, OVERFLOW      - NO COMPUTATION, REAL(Z) TOO
    //                            LARGE ON KODE=1
    //                    IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
    //                            BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
    //                            REDUCTION PRODUCE LESS THAN HALF OF MACHINE
    //                            ACCURACY
    //                    IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
    //                            TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
    //                            CANCE BY ARGUMENT REDUCTION
    //                    IERR=5, ERROR              - NO COMPUTATION,
    //                            ALGORITHM TERMINATION CONDITION NOT MET
    //
    // ***LONG DESCRIPTION
    //
    //         THE COMPUTATION IS CARRIED OUT BY THE POWER SERIES FOR
    //         SMALL CABS(Z), THE ASYMPTOTIC EXPANSION FOR LARGE CABS(Z),
    //         THE MILLER ALGORITHM NORMALIZED BY THE WRONSKIAN AND A
    //         NEUMANN SERIES FOR IMTERMEDIATE MAGNITUDES, AND THE
    //         UNIFORM ASYMPTOTIC EXPANSIONS FOR I(FNU,Z) AND J(FNU,Z)
    //         FOR LARGE ORDERS. BACKWARD RECURRENCE IS USED TO GENERATE
    //         SEQUENCES OR REDUCE ORDERS WHEN NECESSARY.
    //
    //         THE CALCULATIONS ABOVE ARE DONE IN THE RIGHT HALF PLANE AND
    //         CONTINUED INTO THE LEFT HALF PLANE BY THE FORMULA
    //
    //         I(FNU,Z*EXP(M*PI)) = EXP(M*PI*FNU)*I(FNU,Z)  REAL(Z) > 0.0
    //                       M = +I OR -I,  I**2=-1
    //
    //         FOR NEGATIVE ORDERS,THE FORMULA
    //
    //              I(-FNU,Z) = I(FNU,Z) + (2/PI)*SIN(PI*FNU)*K(FNU,Z)
    //
    //         CAN BE USED. HOWEVER,FOR LARGE ORDERS CLOSE TO INTEGERS, THE
    //         THE FUNCTION CHANGES RADICALLY. WHEN FNU IS A LARGE POSITIVE
    //         INTEGER,THE MAGNITUDE OF I(-FNU,Z)=I(FNU,Z) IS A LARGE
    //         NEGATIVE POWER OF TEN. BUT WHEN FNU IS NOT AN INTEGER,
    //         K(FNU,Z) DOMINATES IN MAGNITUDE WITH A LARGE POSITIVE POWER OF
    //         TEN AND THE MOST THAT THE SECOND TERM CAN BE REDUCED IS BY
    //         UNIT ROUNDOFF FROM THE COEFFICIENT. THUS, WIDE CHANGES CAN
    //         OCCUR WITHIN UNIT ROUNDOFF OF A LARGE INTEGER FOR FNU. HERE,
    //         LARGE MEANS FNU > CABS(Z).
    //
    //         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
    //         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
    //         LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
    //         CONSEQUENTLY, if EITHER ONE EXCEEDS U1=SQRT(0.5/UR), THEN
    //         LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
    //         IERR=3 IS TRIGGERED WHERE UR=DMAX1(d1mach(4),1.0e-18) IS
    //         DOUBLE PRECISION UNIT ROUNDOFF LIMITED TO 18 DIGITS PRECISION.
    //         if EITHER IS LARGER THAN U2=0.5/UR, THEN ALL SIGNIFICANCE IS
    //         LOST AND IERR=4. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS
    //         MUST BE FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE
    //         INTEGER, U3=i1mach(9). THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS
    //         RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, AND U3
    //         ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION
    //         ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE PRECISION
    //         ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN
    //         THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT ONE CAN EXPECT
    //         TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS
    //         IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
    //         SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.
    //
    //         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
    //         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
    //         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
    //         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
    //         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
    //         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
    //         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
    //         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
    //         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
    //         SEVERAL ORDERS OF MAGNITUDE. if ONE COMPONENT IS 10**K LARGER
    //         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
    //         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
    //         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
    //         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
    //         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
    //         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
    //         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
    //         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
    //         OR -PI/2+P.
    //
    // ***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
    //                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
    //                 COMMERCE, 1955.
    //
    //               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
    //                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
    //
    //               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
    //                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
    //
    //               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
    //                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
    //                 1018, MAY, 1985
    //
    //               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
    //                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, ACM
    //                 TRANS. MATH. SOFTWARE, VOL. 12, NO. 3, SEPTEMBER 1986,
    //                 PP 265-273.

    sanitise_inputs(z, order, n, false)?;

    let abs_z = z.abs();
    let modified_order = order + ((n - 1) as f64);
    let partial_significance_loss = is_sigificance_lost(abs_z, modified_order, false)?;

    let (zn, mut csgn) = if z.re >= 0.0 {
        (z, c_one())
    } else {
        //-----------------------------------------------------------------------
        //     CALCULATE CSGN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
        //     WHEN FNU IS LARGE
        //-----------------------------------------------------------------------
        let integer_order = order as usize;
        let arg = order.fract() * PI * if z.im < 0.0 { -1.0 } else { 1.0 };
        let mut csgn = Complex64::cis(arg);
        if !integer_order.is_multiple_of(2) {
            csgn = -csgn;
        }
        (-z, csgn)
    };
    let (mut y, nz) = i_right_half_plane(zn, order, scaling, n)?;
    let remaining_n = n - nz;
    if z.re < 0.0 && remaining_n > 0 {
        //-----------------------------------------------------------------------
        //     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE
        //-----------------------------------------------------------------------
        for yi in y.iter_mut().take(remaining_n) {
            let correction =
                if max_abs_component(*yi) <= MACHINE_CONSTANTS.absolute_approximation_limit {
                    *yi *= MACHINE_CONSTANTS.rtol;
                    MACHINE_CONSTANTS.abs_error_tolerance
                } else {
                    1.0
                };
            *yi *= csgn;
            *yi *= correction;
            csgn = -csgn;
        }
    }

    if partial_significance_loss {
        Err(PartialLossOfSignificance { y, nz })
    } else {
        Ok((y, nz))
    }
}

pub fn complex_bessel_j(z: Complex64, order: f64, scaling: Scaling, n: usize) -> BesselResult {
    // ***BEGIN PROLOGUE  ZBESJ
    // ***DATE WRITTEN   830501   (YYMMDD)
    // ***REVISION DATE  890801, 930101   (YYMMDD)
    // ***CATEGORY NO.  B5K
    // ***KEYWORDS  J-BESSEL FUNCTION,BESSEL FUNCTION OF COMPLEX ARGUMENT,
    //             BESSEL FUNCTION OF FIRST KIND
    // ***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
    // ***PURPOSE  TO COMPUTE THE J-BESSEL FUNCTION OF A COMPLEX ARGUMENT
    // ***DESCRIPTION
    //
    //                      ***A DOUBLE PRECISION ROUTINE***
    //         ON KODE=1, ZBESJ COMPUTES AN N MEMBER  SEQUENCE OF COMPLEX
    //         BESSEL FUNCTIONS CY(I)=J(FNU+I-1,Z) FOR REAL, NONNEGATIVE
    //         ORDERS FNU+I-1, I=1,...,N AND COMPLEX Z IN THE CUT PLANE
    //         -PI < ARG(Z) <= PI. ON KODE=2, ZBESJ RETURNS THE SCALED
    //         FUNCTIONS
    //
    //         CY(I)=EXP(-ABS(Y))*J(FNU+I-1,Z)   I = 1,...,N , Y=AIMAG(Z)
    //
    //         WHICH REMOVE THE EXPONENTIAL GROWTH IN BOTH THE UPPER AND
    //         LOWER HALF PLANES FOR Z TO INFINITY. DEFINITIONS AND NOTATION
    //         ARE FOUND IN THE NBS HANDBOOK OF MATHEMATICAL FUNCTIONS
    //         (REF. 1).
    //
    //         INPUT      ZR,ZI,FNU ARE DOUBLE PRECISION
    //           ZR,ZI  - Z=CMPLX(ZR,ZI),  -PI < ARG(Z) <= PI
    //           FNU    - ORDER OF INITIAL J FUNCTION, FNU >= 0.0
    //           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
    //                    KODE= 1  RETURNS
    //                             CY(I)=J(FNU+I-1,Z), I=1,...,N
    //                        = 2  RETURNS
    //                             CY(I)=J(FNU+I-1,Z)EXP(-ABS(Y)), I=1,...,N
    //           N      - NUMBER OF MEMBERS OF THE SEQUENCE, N >= 1
    //
    //         OUTPUT     CYR,CYI ARE DOUBLE PRECISION
    //           CYR,CYI- DOUBLE PRECISION VECTORS WHOSE FIRST N COMPONENTS
    //                    CONTAIN REAL AND IMAGINARY PARTS FOR THE SEQUENCE
    //                    CY(I)=J(FNU+I-1,Z)  OR
    //                    CY(I)=J(FNU+I-1,Z)EXP(-ABS(Y))  I=1,...,N
    //                    DEPENDING ON KODE, Y=AIMAG(Z).
    //           NZ     - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW,
    //                    NZ= 0   , NORMAL RETURN
    //                    NZ > 0 , LAST NZ COMPONENTS OF CY SET  ZERO DUE
    //                              TO UNDERFLOW, CY(I)=CMPLX(0.0,0.0),
    //                              I = N-NZ+1,...,N
    //           IERR   - ERROR FLAG
    //                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
    //                    IERR=1, INPUT ERROR   - NO COMPUTATION
    //                    IERR=2, OVERFLOW      - NO COMPUTATION, AIMAG(Z)
    //                            TOO LARGE ON KODE=1
    //                    IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
    //                            BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
    //                            REDUCTION PRODUCE LESS THAN HALF OF MACHINE
    //                            ACCURACY
    //                    IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
    //                            TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
    //                            CANCE BY ARGUMENT REDUCTION
    //                    IERR=5, ERROR              - NO COMPUTATION,
    //                            ALGORITHM TERMINATION CONDITION NOT MET
    //
    // ***LONG DESCRIPTION
    //
    //         THE COMPUTATION IS CARRIED OUT BY THE FORMULA
    //
    //         J(FNU,Z)=EXP( FNU*PI*I/2)*I(FNU,-I*Z)    AIMAG(Z) >= 0.0
    //
    //         J(FNU,Z)=EXP(-FNU*PI*I/2)*I(FNU, I*Z)    AIMAG(Z) < 0.0
    //
    //         WHERE I**2 = -1 AND I(FNU,Z) IS THE I BESSEL FUNCTION.
    //
    //         FOR NEGATIVE ORDERS,THE FORMULA
    //
    //              J(-FNU,Z) = J(FNU,Z)*COS(PI*FNU) - Y(FNU,Z)*SIN(PI*FNU)
    //
    //         CAN BE USED. HOWEVER,FOR LARGE ORDERS CLOSE TO INTEGERS, THE
    //         THE FUNCTION CHANGES RADICALLY. WHEN FNU IS A LARGE POSITIVE
    //         INTEGER,THE MAGNITUDE OF J(-FNU,Z)=J(FNU,Z)*COS(PI*FNU) IS A
    //         LARGE NEGATIVE POWER OF TEN. BUT WHEN FNU IS NOT AN INTEGER,
    //         Y(FNU,Z) DOMINATES IN MAGNITUDE WITH A LARGE POSITIVE POWER OF
    //         TEN AND THE MOST THAT THE SECOND TERM CAN BE REDUCED IS BY
    //         UNIT ROUNDOFF FROM THE COEFFICIENT. THUS, WIDE CHANGES CAN
    //         OCCUR WITHIN UNIT ROUNDOFF OF A LARGE INTEGER FOR FNU. HERE,
    //         LARGE MEANS FNU > CABS(Z).
    //
    //         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
    //         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
    //         LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
    //         CONSEQUENTLY, if EITHER ONE EXCEEDS U1=SQRT(0.5/UR), THEN
    //         LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
    //         IERR=3 IS TRIGGERED WHERE UR=DMAX1(d1mach(4),1.0e-18) IS
    //         DOUBLE PRECISION UNIT ROUNDOFF LIMITED TO 18 DIGITS PRECISION.
    //         if EITHER IS LARGER THAN U2=0.5/UR, THEN ALL SIGNIFICANCE IS
    //         LOST AND IERR=4. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS
    //         MUST BE FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE
    //         INTEGER, U3=i1mach(9). THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS
    //         RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, AND U3
    //         ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION
    //         ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE PRECISION
    //         ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN
    //         THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT ONE CAN EXPECT
    //         TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS
    //         IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
    //         SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.
    //
    //         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
    //         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
    //         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
    //         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
    //         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
    //         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
    //         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
    //         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
    //         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
    //         SEVERAL ORDERS OF MAGNITUDE. if ONE COMPONENT IS 10**K LARGER
    //         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
    //         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
    //         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
    //         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
    //         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
    //         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
    //         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
    //         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
    //         OR -PI/2+P.
    //
    // ***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
    //                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
    //                 COMMERCE, 1955.
    //
    //               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
    //                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
    //
    //               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
    //                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
    //
    //               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
    //                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
    //                 1018, MAY, 1985
    //
    //               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
    //                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, ACM
    //                 TRANS. MATH. SOFTWARE, VOL. 12, NO. 3, SEPTEMBER 1986,
    //                 PP 265-273.
    //

    sanitise_inputs(z, order, n, false)?;

    let partial_significance_loss = is_sigificance_lost(z.abs(), order + ((n - 1) as f64), false)?;
    //-----------------------------------------------------------------------
    //     CALCULATE CSGN=EXP(FNU*FRAC_PI_2*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
    //     WHEN FNU IS LARGE
    //-----------------------------------------------------------------------
    let order_int = order as i64;
    let half_order_int = order_int / 2;
    let order_rounded_down_to_even = 2 * half_order_int;
    let arg = (order - (order_rounded_down_to_even as f64)) * FRAC_PI_2;
    let mut csgn = Complex64::cis(arg);
    if (half_order_int % 2) != 0 {
        csgn = -csgn;
    }
    //-----------------------------------------------------------------------
    //     ZN IS IN THE RIGHT HALF PLANE
    //-----------------------------------------------------------------------
    let mut sign_selector = 1.0;
    let mut zn = -Complex64::I * z;
    if z.im < 0.0 {
        zn = -zn;
        csgn.im = -csgn.im;
        sign_selector = -sign_selector;
    }
    let (mut cy, nz) = i_right_half_plane(zn, order, scaling, n)?;
    for cyi in cy.iter_mut().take(n - nz) {
        let mut ATOL = 1.0;
        // TODO is the below a pattern?
        if (max_abs_component(*cyi)) <= MACHINE_CONSTANTS.absolute_approximation_limit {
            *cyi *= MACHINE_CONSTANTS.rtol;
            ATOL = MACHINE_CONSTANTS.abs_error_tolerance;
        }
        *cyi *= csgn * ATOL;
        csgn *= sign_selector * Complex64::I;
    }
    if partial_significance_loss {
        Err(PartialLossOfSignificance { y: cy, nz })
    } else {
        Ok((cy, nz))
    }
}

pub fn complex_bessel_k(z: Complex64, order: f64, scaling: Scaling, n: usize) -> BesselResult {
    // ***BEGIN PROLOGUE  ZBESK
    // ***DATE WRITTEN   830501   (YYMMDD)
    // ***REVISION DATE  890801, 930101   (YYMMDD)
    // ***CATEGORY NO.  B5K
    // ***KEYWORDS  K-BESSEL FUNCTION,COMPLEX BESSEL FUNCTION,
    //             MODIFIED BESSEL FUNCTION OF THE SECOND KIND,
    //             BESSEL FUNCTION OF THE THIRD KIND
    // ***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
    // ***PURPOSE  TO COMPUTE K-BESSEL FUNCTIONS OF COMPLEX ARGUMENT
    // ***DESCRIPTION
    //
    //                      ***A DOUBLE PRECISION ROUTINE***
    //
    //         ON KODE=1, ZBESK COMPUTES AN N MEMBER SEQUENCE OF COMPLEX
    //         BESSEL FUNCTIONS CY(J)=K(FNU+J-1,Z) FOR REAL, NONNEGATIVE
    //         ORDERS FNU+J-1, J=1,...,N AND COMPLEX Z != CMPLX(0.0,0.0)
    //         IN THE CUT PLANE -PI < ARG(Z) <= PI. ON KODE=2, ZBESK
    //         RETURNS THE SCALED K FUNCTIONS,
    //
    //         CY(J)=EXP(Z)*K(FNU+J-1,Z) , J=1,...,N,
    //
    //         WHICH REMOVE THE EXPONENTIAL BEHAVIOR IN BOTH THE LEFT AND
    //         RIGHT HALF PLANES FOR Z TO INFINITY. DEFINITIONS AND
    //         NOTATION ARE FOUND IN THE NBS HANDBOOK OF MATHEMATICAL
    //         FUNCTIONS (REF. 1).
    //
    //         INPUT      ZR,ZI,FNU ARE DOUBLE PRECISION
    //           ZR,ZI  - Z=CMPLX(ZR,ZI), Z != CMPLX(0.0,0.0),
    //                    -PI < ARG(Z) <= PI
    //           FNU    - ORDER OF INITIAL K FUNCTION, FNU >= 0.0
    //           N      - NUMBER OF MEMBERS OF THE SEQUENCE, N >= 1
    //           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
    //                    KODE= 1  RETURNS
    //                             CY(I)=K(FNU+I-1,Z), I=1,...,N
    //                        = 2  RETURNS
    //                             CY(I)=K(FNU+I-1,Z)*EXP(Z), I=1,...,N
    //
    //         OUTPUT     CYR,CYI ARE DOUBLE PRECISION
    //           CYR,CYI- DOUBLE PRECISION VECTORS WHOSE FIRST N COMPONENTS
    //                    CONTAIN REAL AND IMAGINARY PARTS FOR THE SEQUENCE
    //                    CY(I)=K(FNU+I-1,Z), I=1,...,N OR
    //                    CY(I)=K(FNU+I-1,Z)*EXP(Z), I=1,...,N
    //                    DEPENDING ON KODE
    //           NZ     - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW.
    //                    NZ= 0   , NORMAL RETURN
    //                    NZ > 0 , FIRST NZ COMPONENTS OF CY SET TO ZERO DUE
    //                              TO UNDERFLOW, CY(I)=CMPLX(0.0,0.0),
    //                              I=1,...,N WHEN X >= 0.0. WHEN X < 0.0
    //                              NZ STATES ONLY THE NUMBER OF UNDERFLOWS
    //                              IN THE SEQUENCE.
    //
    //           IERR   - ERROR FLAG
    //                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
    //                    IERR=1, INPUT ERROR   - NO COMPUTATION
    //                    IERR=2, OVERFLOW      - NO COMPUTATION, FNU IS
    //                            TOO LARGE OR CABS(Z) IS TOO SMALL OR BOTH
    //                    IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
    //                            BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
    //                            REDUCTION PRODUCE LESS THAN HALF OF MACHINE
    //                            ACCURACY
    //                    IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
    //                            TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
    //                            CANCE BY ARGUMENT REDUCTION
    //                    IERR=5, ERROR              - NO COMPUTATION,
    //                            ALGORITHM TERMINATION CONDITION NOT MET
    //
    // ***LONG DESCRIPTION
    //
    //         EQUATIONS OF THE REFERENCE ARE IMPLEMENTED FOR SMALL ORDERS
    //         DNU AND DNU+1.0 IN THE RIGHT HALF PLANE X >= 0.0. FORWARD
    //         RECURRENCE GENERATES HIGHER ORDERS. K IS CONTINUED TO THE LEFT
    //         HALF PLANE BY THE RELATION
    //
    //         K(FNU,Z*EXP(MP)) = EXP(-MP*FNU)*K(FNU,Z)-MP*I(FNU,Z)
    //         MP=MR*PI*I, MR=+1 OR -1, RE(Z) > 0, I**2=-1
    //
    //         WHERE I(FNU,Z) IS THE I BESSEL FUNCTION.
    //
    //         FOR LARGE ORDERS, FNU > FNUL, THE K FUNCTION IS COMPUTED
    //         BY MEANS OF ITS UNIFORM ASYMPTOTIC EXPANSIONS.
    //
    //         FOR NEGATIVE ORDERS, THE FORMULA
    //
    //                       K(-FNU,Z) = K(FNU,Z)
    //
    //         CAN BE USED.
    //
    //         ZBESK ASSUMES THAT A SIGNIFICANT DIGIT SINH(X) FUNCTION IS
    //         AVAILABLE.
    //
    //         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
    //         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
    //         LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
    //         CONSEQUENTLY, if EITHER ONE EXCEEDS U1=SQRT(0.5/UR), THEN
    //         LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
    //         IERR=3 IS TRIGGERED WHERE UR=DMAX1(d1mach(4),1.0e-18) IS
    //         DOUBLE PRECISION UNIT ROUNDOFF LIMITED TO 18 DIGITS PRECISION.
    //         if EITHER IS LARGER THAN U2=0.5/UR, THEN ALL SIGNIFICANCE IS
    //         LOST AND IERR=4. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS
    //         MUST BE FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE
    //         INTEGER, U3=i1mach(9). THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS
    //         RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, AND U3
    //         ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION
    //         ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE PRECISION
    //         ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN
    //         THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT ONE CAN EXPECT
    //         TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS
    //         IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
    //         SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.
    //
    //         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
    //         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
    //         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
    //         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
    //         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
    //         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
    //         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
    //         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
    //         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
    //         SEVERAL ORDERS OF MAGNITUDE. if ONE COMPONENT IS 10**K LARGER
    //         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
    //         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
    //         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
    //         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
    //         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
    //         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
    //         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
    //         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
    //         OR -PI/2+P.
    //
    // ***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
    //                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
    //                 COMMERCE, 1955.
    //
    //               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
    //                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
    //
    //               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
    //                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983.
    //
    //               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
    //                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
    //                 1018, MAY, 1985
    //
    //               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
    //                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, ACM
    //                 TRANS. MATH. SOFTWARE, VOL. 12, NO. 3, SEPTEMBER 1986,
    //                 PP 265-273.

    sanitise_inputs(z, order, n, true)?;
    //-----------------------------------------------------------------------------;
    //     TEST FOR PROPER RANGE;
    //-----------------------------------------------------------------------;
    let abs_z = z.abs();
    let modified_order = order + ((n - 1) as f64);
    let partial_significance_loss = is_sigificance_lost(abs_z, modified_order, false)?;

    //-----------------------------------------------------------------------;
    //     OVERFLOW TEST ON THE LAST MEMBER OF THE SEQUENCE;
    //-----------------------------------------------------------------------;
    if abs_z < MACHINE_CONSTANTS.underflow_limit {
        return Err(Overflow);
    }

    let mut nz = 0;
    if order > MACHINE_CONSTANTS.asymptotic_order_limit {
        //-----------------------------------------------------------------------
        //     UNIFORM ASYMPTOTIC EXPANSIONS FOR FNU > FNUL
        //-----------------------------------------------------------------------
        let rotation = if z.re >= 0.0 {
            RotationDirection::None
        } else if z.im < 0.0 {
            RotationDirection::Left
        } else {
            RotationDirection::Right
        };

        let (y, nz) = ZBUNK(z, order, scaling, rotation, n)?;
        return if partial_significance_loss {
            Err(PartialLossOfSignificance { y, nz })
        } else {
            Ok((y, nz))
        };
    }

    if modified_order > 2.0 {
        let mut y = c_zeros(n);
        let NUF = check_underflow_uniform_asymp_params(z, order, scaling, IKType::K, n, &mut y)?;
        nz += NUF;

        //-----------------------------------------------------------------------;
        //     HERE NN=N OR NN=0 SINCE NUF=0,NN, OR -1 ON RETURN FROM CUOIK;
        //     if NUF=NN, THEN CY(I)=CZERO FOR ALL I;
        //-----------------------------------------------------------------------;
        if NUF == n {
            return if z.re < 0.0 {
                Err(Overflow)
            } else if partial_significance_loss {
                Err(PartialLossOfSignificance { y, nz })
            } else {
                Ok((y, nz))
            };
        }
    }
    if (modified_order > 1.0) && abs_z <= MACHINE_CONSTANTS.abs_error_tolerance {
        let half_abs_z = 0.5 * abs_z;
        if -modified_order * half_abs_z.ln() > MACHINE_CONSTANTS.exponent_limit {
            return Err(Overflow);
        }
    }
    let (y, nz) = if z.re >= 0.0 {
        //-----------------------------------------------------------------------;
        //     RIGHT HALF PLANE COMPUTATION, REAL(Z) >= 0.;
        //-----------------------------------------------------------------------;
        k_right_half_plane(z, order, scaling, n)?
    } else {
        //-----------------------------------------------------------------------;
        //     LEFT HALF PLANE COMPUTATION;
        //     PI/2 < ARG(Z) <= PI AND -PI < ARG(Z) < -PI/2.;
        //-----------------------------------------------------------------------;
        if nz != 0 {
            return Err(Overflow);
        }
        let rotation = if z.im < 0.0 {
            RotationDirection::Left
        } else {
            RotationDirection::Right
        };
        analytic_continuation(z, order, scaling, rotation, n)?
    };
    if partial_significance_loss {
        Err(PartialLossOfSignificance { y, nz })
    } else {
        Ok((y, nz))
    }
}

pub fn complex_bessel_y(z: Complex64, order: f64, scaling: Scaling, n: usize) -> BesselResult {
    // ***BEGIN PROLOGUE  ZBESY
    // ***DATE WRITTEN   830501   (YYMMDD)
    // ***REVISION DATE  890801, 930101   (YYMMDD)
    // ***CATEGORY NO.  B5K
    // ***KEYWORDS  Y-BESSEL FUNCTION,BESSEL FUNCTION OF COMPLEX ARGUMENT,
    //             BESSEL FUNCTION OF SECOND KIND
    // ***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
    // ***PURPOSE  TO COMPUTE THE Y-BESSEL FUNCTION OF A COMPLEX ARGUMENT
    // ***DESCRIPTION
    //
    //                      ***A DOUBLE PRECISION ROUTINE***
    //
    //         ON KODE=1, ZBESY COMPUTES AN N MEMBER SEQUENCE OF COMPLEX
    //         BESSEL FUNCTIONS CY(I)=Y(FNU+I-1,Z) FOR REAL, NONNEGATIVE
    //         ORDERS FNU+I-1, I=1,...,N AND COMPLEX Z IN THE CUT PLANE
    //         -PI < ARG(Z) <= PI. ON KODE=2, ZBESY RETURNS THE SCALED
    //         FUNCTIONS
    //
    //         CY(I)=EXP(-ABS(Y))*Y(FNU+I-1,Z)   I = 1,...,N , Y=AIMAG(Z)
    //
    //         WHICH REMOVE THE EXPONENTIAL GROWTH IN BOTH THE UPPER AND
    //         LOWER HALF PLANES FOR Z TO INFINITY. DEFINITIONS AND NOTATION
    //         ARE FOUND IN THE NBS HANDBOOK OF MATHEMATICAL FUNCTIONS
    //         (REF. 1).
    //
    //         INPUT      ZR,ZI,FNU ARE DOUBLE PRECISION
    //           ZR,ZI  - Z=CMPLX(ZR,ZI), Z != CMPLX(0.0,0.0),
    //                    -PI < ARG(Z) <= PI
    //           FNU    - ORDER OF INITIAL Y FUNCTION, FNU >= 0.0
    //           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
    //                    KODE= 1  RETURNS
    //                             CY(I)=Y(FNU+I-1,Z), I=1,...,N
    //                        = 2  RETURNS
    //                             CY(I)=Y(FNU+I-1,Z)*EXP(-ABS(Y)), I=1,...,N
    //                             WHERE Y=AIMAG(Z)
    //           N      - NUMBER OF MEMBERS OF THE SEQUENCE, N >= 1
    //           CWRKR, - DOUBLE PRECISION WORK VECTORS OF DIMENSION AT
    //           CWRKI    AT LEAST N
    //
    //         OUTPUT     CYR,CYI ARE DOUBLE PRECISION
    //           CYR,CYI- DOUBLE PRECISION VECTORS WHOSE FIRST N COMPONENTS
    //                    CONTAIN REAL AND IMAGINARY PARTS FOR THE SEQUENCE
    //                    CY(I)=Y(FNU+I-1,Z)  OR
    //                    CY(I)=Y(FNU+I-1,Z)*EXP(-ABS(Y))  I=1,...,N
    //                    DEPENDING ON KODE.
    //           NZ     - NZ=0 , A NORMAL RETURN
    //                    NZ > 0 , NZ COMPONENTS OF CY SET TO ZERO DUE TO
    //                    UNDERFLOW (GENERALLY ON KODE=2)
    //           IERR   - ERROR FLAG
    //                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
    //                    IERR=1, INPUT ERROR   - NO COMPUTATION
    //                    IERR=2, OVERFLOW      - NO COMPUTATION, FNU IS
    //                            TOO LARGE OR CABS(Z) IS TOO SMALL OR BOTH
    //                    IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
    //                            BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
    //                            REDUCTION PRODUCE LESS THAN HALF OF MACHINE
    //                            ACCURACY
    //                    IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
    //                            TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
    //                            CANCE BY ARGUMENT REDUCTION
    //                    IERR=5, ERROR              - NO COMPUTATION,
    //                            ALGORITHM TERMINATION CONDITION NOT MET
    //
    // ***LONG DESCRIPTION
    //
    //         THE COMPUTATION IS CARRIED OUT IN TERMS OF THE I(FNU,Z) AND
    //         K(FNU,Z) BESSEL FUNCTIONS IN THE RIGHT HALF PLANE BY
    //
    //             Y(FNU,Z) = I*CC*I(FNU,ARG) - (2/PI)*CONJG(CC)*K(FNU,ARG)
    //
    //             Y(FNU,Z) = CONJG(Y(FNU,CONJG(Z)))
    //
    //         FOR AIMAG(Z) >= 0 AND AIMAG(Z) < 0 RESPECTIVELY, WHERE
    //         CC=EXP(I*PI*FNU/2), ARG=Z*EXP(-I*PI/2) AND I**2=-1.
    //
    //         FOR NEGATIVE ORDERS,THE FORMULA
    //
    //              Y(-FNU,Z) = Y(FNU,Z)*COS(PI*FNU) + J(FNU,Z)*SIN(PI*FNU)
    //
    //         CAN BE USED. HOWEVER,FOR LARGE ORDERS CLOSE TO HALF ODD
    //         INTEGERS THE FUNCTION CHANGES RADICALLY. WHEN FNU IS A LARGE
    //         POSITIVE HALF ODD INTEGER,THE MAGNITUDE OF Y(-FNU,Z)=J(FNU,Z)*
    //         SIN(PI*FNU) IS A LARGE NEGATIVE POWER OF TEN. BUT WHEN FNU IS
    //         NOT A HALF ODD INTEGER, Y(FNU,Z) DOMINATES IN MAGNITUDE WITH A
    //         LARGE POSITIVE POWER OF TEN AND THE MOST THAT THE SECOND TERM
    //         CAN BE REDUCED IS BY UNIT ROUNDOFF FROM THE COEFFICIENT. THUS,
    //         WIDE CHANGES CAN OCCUR WITHIN UNIT ROUNDOFF OF A LARGE HALF
    //         ODD INTEGER. HERE, LARGE MEANS FNU > CABS(Z).
    //
    //         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
    //         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
    //         LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
    //         CONSEQUENTLY, if EITHER ONE EXCEEDS U1=SQRT(0.5/UR), THEN
    //         LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
    //         IERR=3 IS TRIGGERED WHERE UR=DMAX1(d1mach(4),1.0e-18) IS
    //         DOUBLE PRECISION UNIT ROUNDOFF LIMITED TO 18 DIGITS PRECISION.
    //         if EITHER IS LARGER THAN U2=0.5/UR, THEN ALL SIGNIFICANCE IS
    //         LOST AND IERR=4. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS
    //         MUST BE FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE
    //         INTEGER, U3=i1mach(9). THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS
    //         RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, AND U3
    //         ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION
    //         ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE PRECISION
    //         ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN
    //         THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT ONE CAN EXPECT
    //         TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS
    //         IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
    //         SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.
    //
    //         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
    //         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
    //         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
    //         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
    //         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
    //         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
    //         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
    //         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
    //         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
    //         SEVERAL ORDERS OF MAGNITUDE. if ONE COMPONENT IS 10**K LARGER
    //         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
    //         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
    //         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
    //         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
    //         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
    //         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
    //         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
    //         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
    //         OR -PI/2+P.
    //
    // ***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
    //                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
    //                 COMMERCE, 1955.
    //
    //               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
    //                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
    //
    //               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
    //                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
    //
    //               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
    //                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
    //                 1018, MAY, 1985
    //
    //               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
    //                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, ACM
    //                 TRANS. MATH. SOFTWARE, VOL. 12, NO. 3, SEPTEMBER 1986,
    //                 PP 265-273.

    sanitise_inputs(z, order, n, true)?;
    let zz = if z.im < 0.0 { z.conj() } else { z };
    let zn = -Complex64::I * zz;
    let mut partial_loss_of_significance = false;

    let mut unwrap_psl = |result: BesselResult| match result {
        Ok((y_, nz_)) => Ok((y_, nz_)),
        Err(PartialLossOfSignificance { y: y_, nz: nz_ }) => {
            partial_loss_of_significance = true;
            Ok((y_, nz_))
        }
        err => err,
    };

    let (bess_i, nz_i) = unwrap_psl(complex_bessel_i(zn, order, scaling, n))?;
    let (bess_k, nz_k) = unwrap_psl(complex_bessel_k(zn, order, scaling, n))?;

    let mut nz = nz_i.min(nz_k);
    let frac_order = order.fract();
    let integer_order = order as usize;
    let ARG = FRAC_PI_2 * frac_order;
    let mut csgn = Complex64::cis(ARG);
    let index = integer_order % 4;
    csgn *= CIP[index];
    let mut cspn = csgn.conj() * FRAC_2_PI;
    csgn *= Complex64::I;

    let mut ey = 1.0;
    if scaling == Scaling::Scaled {
        let ex = Complex64::cis(z.re);
        let two_abs_z = 2.0 * z.im.abs();
        ey = if two_abs_z < MACHINE_CONSTANTS.exponent_limit {
            (-two_abs_z).exp()
        } else {
            0.0
        };
        cspn *= ex * ey;
        nz = 0;
    }
    let mut y: Vec<Complex64> = bess_i
        .iter()
        .zip(bess_k)
        .map(|(&z_i, z_k)| {
            //----------------------------------------------------------------------;
            //       CY(I) = CSGN*CY(I)-CSPN*CWRK(I): PRODUCTS ARE COMPUTED IN;
            //       SCALED MODE if CY(I) OR CWRK(I) ARE CLOSE TO UNDERFLOW TO;
            //       PREVENT UNDERFLOW IN AN INTERMEDIATE COMPUTATION.;
            //----------------------------------------------------------------------;
            let z_k = scaled_multiply(z_k, cspn, scaling);
            let z_i = scaled_multiply(z_i, csgn, scaling);
            let val = z_i - z_k;
            if scaling == Scaling::Scaled && val == c_zero() && ey == 0.0 {
                nz += 1;
            }
            csgn *= Complex64::I;
            cspn *= -Complex64::I;
            val
        })
        .collect();

    if z.im < 0.0 {
        y.iter_mut().for_each(|v| *v = v.conj());
    }
    if partial_loss_of_significance {
        Err(PartialLossOfSignificance { y, nz })
    } else {
        Ok((y, nz))
    }
}

fn scaled_multiply(mut z: Complex64, coeff: Complex64, scaling: Scaling) -> Complex64 {
    match scaling {
        Scaling::Unscaled => z * coeff,
        Scaling::Scaled => {
            let atol = if max_abs_component(z) <= MACHINE_CONSTANTS.absolute_approximation_limit {
                z *= MACHINE_CONSTANTS.rtol;
                MACHINE_CONSTANTS.abs_error_tolerance
            } else {
                1.0
            };
            (z * coeff) * atol
        }
    }
}

pub fn complex_airy(
    z: Complex64,
    return_derivative: bool,
    scaling: Scaling,
) -> BesselResult<(Complex64, usize)> {
    // ***BEGIN PROLOGUE  ZAIRY
    // ***DATE WRITTEN   830501   (YYMMDD)
    // ***REVISION DATE  890801, 930101   (YYMMDD)
    // ***CATEGORY NO.  B5K
    // ***KEYWORDS  AIRY FUNCTION,BESSEL FUNCTIONS OF ORDER ONE THIRD
    // ***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
    // ***PURPOSE  TO COMPUTE AIRY FUNCTIONS AI(Z) AND DAI(Z) FOR COMPLEX Z
    // ***DESCRIPTION
    //
    //                      ***A DOUBLE PRECISION ROUTINE***
    //         ON KODE=1, ZAIRY COMPUTES THE COMPLEX AIRY FUNCTION AI(Z) OR
    //         ITS DERIVATIVE DAI(Z)/DZ ON ID=0 OR ID=1 RESPECTIVELY. ON
    //         KODE=2, A SCALING OPTION CEXP(ZTA)*AI(Z) OR CEXP(ZTA)*
    //         DAI(Z)/DZ IS PROVIDED TO REMOVE THE EXPONENTIAL DECAY IN
    //         -PI/3 < ARG(Z) < PI/3 AND THE EXPONENTIAL GROWTH IN
    //         PI/3 < ABS(ARG(Z)) < PI WHERE ZTA=(2/3)*Z*CSQRT(Z).
    //
    //         WHILE THE AIRY FUNCTIONS AI(Z) AND DAI(Z)/DZ ARE ANALYTIC IN
    //         THE WHOLE Z PLANE, THE CORRESPONDING SCALED FUNCTIONS DEFINED
    //         FOR KODE=2 HAVE A CUT ALONG THE NEGATIVE REAL AXIS.
    //         DEFINTIONS AND NOTATION ARE FOUND IN THE NBS HANDBOOK OF
    //         MATHEMATICAL FUNCTIONS (REF. 1).
    //
    //         INPUT      ZR,ZI ARE DOUBLE PRECISION
    //           ZR,ZI  - Z=CMPLX(ZR,ZI)
    //           ID     - ORDER OF DERIVATIVE, ID=0 OR ID=1
    //           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
    //                    KODE= 1  RETURNS
    //                             AI=AI(Z)                ON ID=0 OR
    //                             AI=DAI(Z)/DZ            ON ID=1
    //                        = 2  RETURNS
    //                             AI=CEXP(ZTA)*AI(Z)       ON ID=0 OR
    //                             AI=CEXP(ZTA)*DAI(Z)/DZ   ON ID=1 WHERE
    //                             ZTA=(2/3)*Z*CSQRT(Z)
    //
    //         OUTPUT     AIR,AII ARE DOUBLE PRECISION
    //           AIR,AII- COMPLEX ANSWER DEPENDING ON THE CHOICES FOR ID AND
    //                    KODE
    //           NZ     - UNDERFLOW INDICATOR
    //                    NZ= 0   , NORMAL RETURN
    //                    NZ= 1   , AI=CMPLX(0.0,0.0) DUE TO UNDERFLOW IN
    //                              -PI/3 < ARG(Z) < PI/3 ON KODE=1
    //           IERR   - ERROR FLAG
    //                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
    //                    IERR=1, INPUT ERROR   - NO COMPUTATION
    //                    IERR=2, OVERFLOW      - NO COMPUTATION, REAL(ZTA)
    //                            TOO LARGE ON KODE=1
    //                    IERR=3, CABS(Z) LARGE      - COMPUTATION COMPLETED
    //                            LOSSES OF SIGNIFCANCE BY ARGUMENT REDUCTION
    //                            PRODUCE LESS THAN HALF OF MACHINE ACCURACY
    //                    IERR=4, CABS(Z) TOO LARGE  - NO COMPUTATION
    //                            COMPLETE LOSS OF ACCURACY BY ARGUMENT
    //                            REDUCTION
    //                    IERR=5, ERROR              - NO COMPUTATION,
    //                            ALGORITHM TERMINATION CONDITION NOT MET
    //
    // ***LONG DESCRIPTION
    //
    //         AI AND DAI ARE COMPUTED FOR CABS(Z) > 1.0 FROM THE K BESSEL
    //         FUNCTIONS BY
    //
    //            AI(Z)=C*SQRT(Z)*K(1/3,ZTA) , DAI(Z)=-C*Z*K(2/3,ZTA)
    //                           C=1.0/(PI*SQRT(3.0))
    //                            ZTA=(2/3)*Z**(3/2)
    //
    //         WITH THE POWER SERIES FOR CABS(Z) <= 1.0.
    //
    //         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
    //         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z IS LARGE, LOSSES
    //         OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR. CONSEQUENTLY, if
    //         THE MAGNITUDE OF ZETA=(2/3)*Z**1.5 EXCEEDS U1=SQRT(0.5/UR),
    //         THEN LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR
    //         FLAG IERR=3 IS TRIGGERED WHERE UR=DMAX1(d1mach(4),1.0e-18) IS
    //         DOUBLE PRECISION UNIT ROUNDOFF LIMITED TO 18 DIGITS PRECISION.
    //         ALSO, if THE MAGNITUDE OF ZETA IS LARGER THAN U2=0.5/UR, THEN
    //         ALL SIGNIFICANCE IS LOST AND IERR=4. IN ORDER TO USE THE INT
    //         FUNCTION, ZETA MUST BE FURTHER RESTRICTED NOT TO EXCEED THE
    //         LARGEST INTEGER, U3=i1mach(9). THUS, THE MAGNITUDE OF ZETA
    //         MUST BE RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2,
    //         AND U3 ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE
    //         PRECISION ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE
    //         PRECISION ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMIT-
    //         ING IN THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT THE MAG-
    //         NITUDE OF Z CANNOT EXCEED 3.1E+4 IN SINGLE AND 2.1E+6 IN
    //         DOUBLE PRECISION ARITHMETIC. THIS ALSO MEANS THAT ONE CAN
    //         EXPECT TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES,
    //         NO DIGITS IN SINGLE PRECISION AND ONLY 7 DIGITS IN DOUBLE
    //         PRECISION ARITHMETIC. SIMILAR CONSIDERATIONS HOLD FOR OTHER
    //         MACHINES.
    //
    //         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
    //         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
    //         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
    //         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
    //         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
    //         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
    //         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
    //         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
    //         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
    //         SEVERAL ORDERS OF MAGNITUDE. if ONE COMPONENT IS 10**K LARGER
    //         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
    //         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
    //         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
    //         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
    //         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
    //         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
    //         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
    //         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
    //         OR -PI/2+P.
    //
    // ***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
    //                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
    //                 COMMERCE, 1955.
    //
    //               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
    //                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
    //
    //               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
    //                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
    //                 1018, MAY, 1985
    //
    //               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
    //                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, ACM
    //                 TRANS. MATH. SOFTWARE, VOL. 12, NO. 3, SEPTEMBER 1986,
    //                 PP 265-273.

    const POWER_SERIES_COEFFS: (f64, f64) = (3.55028053887817240e-01, 2.58819403792806799e-01);
    const COEFF: f64 = 1.83776298473930683e-01;

    let abs_z = z.abs();
    let float_is_derivative = if return_derivative { 1.0 } else { 0.0 };
    //--------------------------------------------------------------------------
    //     TEST FOR PROPER RANGE
    //-----------------------------------------------------------------------
    // significance loss only tested against z, not order, so 0.0 is used to never cause significance loss
    let partial_loss_of_significance = is_sigificance_lost(abs_z, 0.0, true)?;

    let retval = if abs_z <= 1.0 {
        //-----------------------------------------------------------------------
        //     POWER SERIES FOR CABS(Z) <= 1.
        //-----------------------------------------------------------------------
        let ai = airy_power_series(z, return_derivative, POWER_SERIES_COEFFS);
        (
            match scaling {
                Scaling::Scaled => ai * (TWO_THIRDS * z * z.sqrt()).exp(),
                Scaling::Unscaled => ai,
            },
            0,
        )
    } else {
        //-----------------------------------------------------------------------
        //     CASE FOR CABS(Z) > 1.0
        //-----------------------------------------------------------------------
        let order = (1.0 + float_is_derivative) / 3.0;
        let ln_abs_z = abs_z.ln();

        let sqrt_z = z.sqrt();
        let mut zta = TWO_THIRDS * z * sqrt_z;
        //-----------------------------------------------------------------------
        //     RE(ZTA) <= 0 WHEN RE(Z) < 0, ESPECIALLY WHEN IM(Z) IS SMALL
        //-----------------------------------------------------------------------
        let mut scale_factor = 1.0;
        if z.re < 0.0 {
            zta.re = -zta.re.abs();
        }
        if z.im == 0.0 && z.re <= 0.0 {
            zta.re = 0.0;
        }
        let re_zta = zta.re;
        let (cy, NZ) = if re_zta < 0.0 || z.re <= 0.0 {
            //-----------------------------------------------------------------------
            //     OVERFLOW TEST
            //-----------------------------------------------------------------------
            if scaling == Scaling::Unscaled && re_zta <= -MACHINE_CONSTANTS.approximation_limit {
                scale_factor = MACHINE_CONSTANTS.abs_error_tolerance;
                if (-re_zta + 0.25 * ln_abs_z) > MACHINE_CONSTANTS.exponent_limit {
                    return Err(Overflow);
                }
            }
            //-----------------------------------------------------------------------
            //     CBKNU AND CACON RETURN EXP(ZTA)*K(FNU,ZTA) ON KODE=2
            //-----------------------------------------------------------------------
            let rotation = if z.im < 0.0 {
                RotationDirection::Left
            } else {
                RotationDirection::Right
            };
            ZACAI(zta, order, scaling, rotation, 1)?
        } else {
            //-----------------------------------------------------------------------
            //     UNDERFLOW TEST
            //-----------------------------------------------------------------------
            let mut retval = None;
            if scaling == Scaling::Unscaled && re_zta > MACHINE_CONSTANTS.approximation_limit {
                scale_factor = 1.0 / MACHINE_CONSTANTS.abs_error_tolerance;
                if (-re_zta - 0.25 * ln_abs_z) < -MACHINE_CONSTANTS.exponent_limit {
                    retval = Some(Ok((c_zeros(1), 1)));
                }
            }
            retval.unwrap_or_else(|| k_right_half_plane(zta, order, scaling, 1))?
        };

        let mut s1 = cy[0] * COEFF * scale_factor;
        s1 *= if return_derivative { -z } else { sqrt_z };
        (s1 / scale_factor, NZ)
    };
    if partial_loss_of_significance {
        Err(PartialLossOfSignificance {
            y: vec![retval.0],
            nz: retval.1,
        })
    } else {
        Ok(retval)
    }
}

pub fn complex_airy_b(
    z: Complex64,
    return_derivative: bool,
    scaling: Scaling,
) -> BesselResult<Complex64> {
    // ***BEGIN PROLOGUE  ZBIRY
    // ***DATE WRITTEN   830501   (YYMMDD)
    // ***REVISION DATE  890801, 930101   (YYMMDD)
    // ***CATEGORY NO.  B5K
    // ***KEYWORDS  AIRY FUNCTION,BESSEL FUNCTIONS OF ORDER ONE THIRD
    // ***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
    // ***PURPOSE  TO COMPUTE AIRY FUNCTIONS BI(Z) AND DBI(Z) FOR COMPLEX Z
    // ***DESCRIPTION
    //
    //                      ***A DOUBLE PRECISION ROUTINE***
    //         ON KODE=1, CBIRY COMPUTES THE COMPLEX AIRY FUNCTION BI(Z) OR
    //         ITS DERIVATIVE DBI(Z)/DZ ON ID=0 OR ID=1 RESPECTIVELY. ON
    //         KODE=2, A SCALING OPTION CEXP(-AXZTA)*BI(Z) OR CEXP(-AXZTA)*
    //         DBI(Z)/DZ IS PROVIDED TO REMOVE THE EXPONENTIAL BEHAVIOR IN
    //         BOTH THE LEFT AND RIGHT HALF PLANES WHERE
    //         ZTA=(2/3)*Z*CSQRT(Z)=CMPLX(XZTA,YZTA) AND AXZTA=ABS(XZTA).
    //         DEFINTIONS AND NOTATION ARE FOUND IN THE NBS HANDBOOK OF
    //         MATHEMATICAL FUNCTIONS (REF. 1).
    //
    //         INPUT      ZR,ZI ARE DOUBLE PRECISION
    //           ZR,ZI  - Z=CMPLX(ZR,ZI)
    //           ID     - ORDER OF DERIVATIVE, ID=0 OR ID=1
    //           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
    //                    KODE= 1  RETURNS
    //                             BI=BI(Z)                 ON ID=0 OR
    //                             BI=DBI(Z)/DZ             ON ID=1
    //                        = 2  RETURNS
    //                             BI=CEXP(-AXZTA)*BI(Z)     ON ID=0 OR
    //                             BI=CEXP(-AXZTA)*DBI(Z)/DZ ON ID=1 WHERE
    //                             ZTA=(2/3)*Z*CSQRT(Z)=CMPLX(XZTA,YZTA)
    //                             AND AXZTA=ABS(XZTA)
    //
    //         OUTPUT     BIR,BII ARE DOUBLE PRECISION
    //           BIR,BII- COMPLEX ANSWER DEPENDING ON THE CHOICES FOR ID AND
    //                    KODE
    //           IERR   - ERROR FLAG
    //                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
    //                    IERR=1, INPUT ERROR   - NO COMPUTATION
    //                    IERR=2, OVERFLOW      - NO COMPUTATION, REAL(Z)
    //                            TOO LARGE ON KODE=1
    //                    IERR=3, CABS(Z) LARGE      - COMPUTATION COMPLETED
    //                            LOSSES OF SIGNIFCANCE BY ARGUMENT REDUCTION
    //                            PRODUCE LESS THAN HALF OF MACHINE ACCURACY
    //                    IERR=4, CABS(Z) TOO LARGE  - NO COMPUTATION
    //                            COMPLETE LOSS OF ACCURACY BY ARGUMENT
    //                            REDUCTION
    //                    IERR=5, ERROR              - NO COMPUTATION,
    //                            ALGORITHM TERMINATION CONDITION NOT MET
    //
    // ***LONG DESCRIPTION
    //
    //         BI AND DBI ARE COMPUTED FOR CABS(Z) > 1.0 FROM THE I BESSEL
    //         FUNCTIONS BY
    //
    //                BI(Z)=C*SQRT(Z)*( I(-1/3,ZTA) + I(1/3,ZTA) )
    //               DBI(Z)=C *  Z  * ( I(-2/3,ZTA) + I(2/3,ZTA) )
    //                               C=1.0/SQRT(3.0)
    //                             ZTA=(2/3)*Z**(3/2)
    //
    //         WITH THE POWER SERIES FOR CABS(Z) <= 1.0.
    //
    //         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
    //         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z IS LARGE, LOSSES
    //         OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR. CONSEQUENTLY, if
    //         THE MAGNITUDE OF ZETA=(2/3)*Z**1.5 EXCEEDS U1=SQRT(0.5/UR),
    //         THEN LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR
    //         FLAG IERR=3 IS TRIGGERED WHERE UR=DMAX1(d1mach(4),1.0e-18) IS
    //         DOUBLE PRECISION UNIT ROUNDOFF LIMITED TO 18 DIGITS PRECISION.
    //         ALSO, if THE MAGNITUDE OF ZETA IS LARGER THAN U2=0.5/UR, THEN
    //         ALL SIGNIFICANCE IS LOST AND IERR=4. IN ORDER TO USE THE INT
    //         FUNCTION, ZETA MUST BE FURTHER RESTRICTED NOT TO EXCEED THE
    //         LARGEST INTEGER, U3=i1mach(9). THUS, THE MAGNITUDE OF ZETA
    //         MUST BE RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2,
    //         AND U3 ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE
    //         PRECISION ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE
    //         PRECISION ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMIT-
    //         ING IN THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT THE MAG-
    //         NITUDE OF Z CANNOT EXCEED 3.1E+4 IN SINGLE AND 2.1E+6 IN
    //         DOUBLE PRECISION ARITHMETIC. THIS ALSO MEANS THAT ONE CAN
    //         EXPECT TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES,
    //         NO DIGITS IN SINGLE PRECISION AND ONLY 7 DIGITS IN DOUBLE
    //         PRECISION ARITHMETIC. SIMILAR CONSIDERATIONS HOLD FOR OTHER
    //         MACHINES.
    //
    //         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
    //         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
    //         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
    //         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
    //         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
    //         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
    //         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
    //         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
    //         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
    //         SEVERAL ORDERS OF MAGNITUDE. if ONE COMPONENT IS 10**K LARGER
    //         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
    //         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
    //         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
    //         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
    //         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
    //         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
    //         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
    //         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
    //         OR -PI/2+P.
    //
    // ***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
    //                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
    //                 COMMERCE, 1955.
    //
    //               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
    //                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
    //
    //               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
    //                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
    //                 1018, MAY, 1985
    //
    //               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
    //                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, ACM
    //                 TRANS. MATH. SOFTWARE, VOL. 12, NO. 3, SEPTEMBER 1986,
    //                 PP 265-273.
    const POWER_SERIES_COEFFS: (f64, f64) = (6.14926627446000736e-01, -4.48288357353826359e-01);
    const COEF: f64 = 5.77350269189625765e-01;

    let abs_z = z.abs();
    let float_is_derivative = if return_derivative { 1.0 } else { 0.0 };
    let mut partial_loss_of_significance = false;

    let bi = if abs_z <= 1.0 {
        //-----------------------------------------------------------------------
        //     POWER SERIES FOR CABS(Z) <= 1.
        //-----------------------------------------------------------------------
        let bi = airy_power_series(z, return_derivative, POWER_SERIES_COEFFS);
        match scaling {
            Scaling::Scaled => {
                //TODO ZTA used many places with similar definition
                let zta = TWO_THIRDS * (z * z.sqrt());
                bi * (-(zta.re.abs())).exp()
            }
            Scaling::Unscaled => bi,
        }
    } else {
        //-----------------------------------------------------------------------;
        //     CASE FOR CABS(Z) > 1.0;
        //-----------------------------------------------------------------------;
        let order = (1.0 + float_is_derivative) / 3.0;
        //-----------------------------------------------------------------------;
        //     TEST FOR RANGE;
        //-----------------------------------------------------------------------;
        // significance loss only tested against z, not order, so 0.0 is used to never cause significance loss
        partial_loss_of_significance = is_sigificance_lost(abs_z, 0.0, true)?;
        let mut scale_factor = 1.0;
        let mut zta = TWO_THIRDS * (z * z.sqrt());

        //-----------------------------------------------------------------------;
        //     RE(ZTA) <= 0 WHEN RE(Z) < 0, ESPECIALLY WHEN IM(Z) IS SMALL;
        //-----------------------------------------------------------------------;
        if z.re < 0.0 {
            zta.re = -zta.re.abs();
        }
        if z.im == 0.0 && z.re < 0.0 {
            zta.re = 0.0;
        }
        if scaling == Scaling::Unscaled {
            //-----------------------------------------------------------------------;
            //     OVERFLOW TEST;
            //-----------------------------------------------------------------------;
            let re_zta = zta.re.abs();
            if re_zta > MACHINE_CONSTANTS.approximation_limit {
                scale_factor = MACHINE_CONSTANTS.abs_error_tolerance;
                if re_zta + 0.25 * abs_z.ln() > MACHINE_CONSTANTS.exponent_limit {
                    return Err(Overflow);
                }
            }
        }
        let mut rotation_angle = 0.0;
        if zta.re < 0.0 || z.re <= 0.0 {
            rotation_angle = PI;
            if z.im < 0.0 {
                rotation_angle = -PI;
            }
            zta *= -1.0;
        }
        //-----------------------------------------------------------------------;
        //     AA=FACTOR FOR ANALYTIC CONTINUATION OF I(FNU,ZTA);
        //     KODE=2 RETURNS EXP(-ABS(XZTA))*I(FNU,ZTA) FROM ZBESI;
        //-----------------------------------------------------------------------;
        let (cy, _) = i_right_half_plane(zta, order, scaling, 1)?;
        let mut s1 = Complex64::cis(rotation_angle * order) * cy[0] * scale_factor;
        let order = (2.0 - float_is_derivative) / 3.0;
        let (mut cy, _) = i_right_half_plane(zta, order, scaling, 2)?;
        cy[0] *= scale_factor;
        cy[1] *= scale_factor;

        //-----------------------------------------------------------------------;
        //     BACKWARD RECUR ONE STEP FOR ORDERS -1/3 OR -2/3;
        //-----------------------------------------------------------------------;
        let s2 = (2.0 * order) * (cy[0] / zta) + cy[1];
        s1 = COEF * (s1 + s2 * Complex64::cis(rotation_angle * (order - 1.0)));
        let z_factor = if return_derivative { z } else { z.sqrt() };
        s1 * z_factor / scale_factor
    };
    if partial_loss_of_significance {
        Err(PartialLossOfSignificance { y: vec![bi], nz: 0 })
    } else {
        Ok(bi)
    }
}

fn airy_power_series(z: Complex64, return_derivative: bool, coeffs: (f64, f64)) -> Complex64 {
    let float_is_derivative = if return_derivative { 1.0 } else { 0.0 };

    let abs_z = z.abs();
    let z_floor = if abs_z < MACHINE_CONSTANTS.underflow_limit {
        c_zero()
    } else {
        z
    };
    let (s1, s2) = if abs_z < MACHINE_CONSTANTS.abs_error_tolerance {
        (c_one(), c_one())
    } else {
        let abs_z_sq = abs_z * abs_z;
        let mut s1 = c_one();
        let mut s2 = c_one();

        if abs_z_sq >= MACHINE_CONSTANTS.abs_error_tolerance / abs_z {
            let mut term1 = c_one();
            let mut term2 = c_one();
            let mut a_term = 1.0;
            let z3 = z.pow(3.0);
            let AZ3 = abs_z * abs_z_sq;
            let (AK, BK, CK, DK) = (
                2.0 + float_is_derivative,
                3.0 - 2.0 * float_is_derivative,
                4.0 - float_is_derivative,
                3.0 + 2.0 * float_is_derivative,
            );
            let mut D1: f64 = AK * DK;
            let mut D2 = BK * CK;
            let mut AD = D1.min(D2);
            let mut AK = 24.0 + 9.0 * float_is_derivative;
            let mut BK = 30.0 - 9.0 * float_is_derivative;
            for _ in 0..25 {
                term1 = term1 * z3 / D1;
                s1 += term1;
                term2 = term2 * z3 / D2;
                s2 += term2;
                a_term = a_term * AZ3 / AD;
                D1 += AK;
                D2 += BK;
                AD = D1.min(D2);
                if a_term < MACHINE_CONSTANTS.abs_error_tolerance * AD {
                    break;
                }
                AK += 18.0;
                BK += 18.0;
            }
        }
        (s1, s2)
    };
    let (c1, c2) = coeffs;
    if return_derivative {
        (c1 / 2.0) * z_floor.pow(2.0) * s1 - s2 * c2
    } else {
        s1 * c1 - c2 * z_floor * s2
    }
}

/// zbknu computes the k bessel function in the right half z plane.
/// Originally ZBKNU
fn k_right_half_plane(z: Complex64, order: f64, scaling: Scaling, n: usize) -> BesselResult {
    const KMAX: usize = 30;
    const RTFRAC_PI_2: f64 = 1.25331413731550025;
    const SPI: f64 = 1.90985931710274403;
    const FPI: f64 = 1.89769999331517738;
    const CC: [f64; 8] = [
        5.77215664901532861e-01,
        -4.20026350340952355e-02,
        -4.21977345555443367e-02,
        7.21894324666309954e-03,
        -2.15241674114950973e-04,
        -2.01348547807882387e-05,
        1.13302723198169588e-06,
        6.11609510448141582e-09,
    ];

    let abs_z = z.abs();
    let mut nz = 0;
    let mut underflow_occurred = false;
    let mut overflow_state;
    let rz = 2.0 * z.conj() / abs_z.powi(2);
    let mut integer_order = (order + 0.5) as isize; // round to nearest int
    let simple_case = integer_order == 0 && n == 1;

    let signed_fractional_order = order - (integer_order as f64); // signed fractional part (-0.5 < DNU < 0.5 )
    let frac_order_sqr = if signed_fractional_order.abs() > MACHINE_CONSTANTS.abs_error_tolerance {
        signed_fractional_order * signed_fractional_order
    } else {
        0.0
    };

    let (mut s1, mut s2) = if (signed_fractional_order.abs() != 0.5) && (abs_z <= 2.0) {
        // series for (z.abs() <= 2.0) and not half integer order
        let mut fc = 1.0;
        let mut smu = rz.ln();
        let fmu = smu * signed_fractional_order;
        let csh = fmu.sinh();
        let cch = fmu.cosh();
        if signed_fractional_order != 0.0 {
            fc = signed_fractional_order * PI;
            fc /= fc.sin();
            smu = csh / signed_fractional_order;
        }
        //-----------------------------------------------------------------------;
        //     GAM(1-Z)*GAM(1+Z)=PI*Z/SIN(PI*Z), T1=1/GAM(1-DNU), T2=1/GAM(1+DNU);
        //-----------------------------------------------------------------------;
        let t2 = (-gamma_ln(1.0 + signed_fractional_order).unwrap()).exp();
        let t1 = 1.0 / (t2 * fc);

        let g1 = if signed_fractional_order.abs() <= 0.1 {
            //-----------------------------------------------------------------------;
            //     SERIES FOR F0 TO RESOLVE INDETERMINACY FOR SMALL ABS(DNU);
            //-----------------------------------------------------------------------;
            let mut ak = 1.0;
            let mut sum = CC[0];
            for cc in CC[1..].iter() {
                ak *= frac_order_sqr;
                let tm = cc * ak;
                sum += tm;
                if tm.abs() < MACHINE_CONSTANTS.abs_error_tolerance {
                    break;
                }
            }
            -sum
        } else {
            (t1 - t2) / (2.0 * signed_fractional_order)
        };
        let g2 = (t1 + t2) * 0.5;
        let f = fc * (g1 * cch + g2 * smu);
        let p = 0.5 * fmu.exp() / t2;
        let q = (0.5 / fmu.exp()) / t1;

        let ck = c_one();
        if simple_case {
            //-----------------------------------------------------------------------;
            //     SPECIAL CASE
            //     GENERATE K(FNU,Z), 0.0  <=  FNU  <  0.5 AND N=1;
            //-----------------------------------------------------------------------;
            let (s1, _) =
                k_right_half_plane_helper(z, frac_order_sqr, signed_fractional_order, f, p, q, ck);

            let mut y = s1;
            if scaling == Scaling::Scaled {
                y *= z.exp();
            }
            return Ok((vec![y], nz));
        }

        //-----------------------------------------------------------------------;
        //     GENERATE K(DNU,Z) AND K(DNU+1,Z) FOR FORWARD RECURRENCE;
        //-----------------------------------------------------------------------;
        let (mut s1, mut s2) =
            k_right_half_plane_helper(z, frac_order_sqr, signed_fractional_order, f, p, q, ck);

        overflow_state = if (order + 1.0) * smu.re.abs() > MACHINE_CONSTANTS.approximation_limit {
            Overflow::NearOver
        } else {
            Overflow::None
        };
        s2 *= MACHINE_CONSTANTS.scaling_factors[overflow_state] * rz;
        s1 *= MACHINE_CONSTANTS.scaling_factors[overflow_state];
        if scaling == Scaling::Scaled {
            let z_exp = z.exp();
            s1 *= z_exp;
            s2 *= z_exp;
        }
        (s1, s2)
    } else {
        // alternative series for z.abs() > 2.0 or half integer order
        //-----------------------------------------------------------------------;
        //     underflow_occured=0 MEANS NO UNDERFLOW OCCURRED;
        //     underflow_occured=1 MEANS AN UNDERFLOW OCCURRED- COMPUTATION PROCEEDS WITH;
        //     KODED=2 AND A TEST FOR ON SCALE VALUES IS MADE DURING FORWARD;
        //     RECURSION;
        //-----------------------------------------------------------------------;
        let mut coeff = Complex64::new(RTFRAC_PI_2, 0.0) / z.sqrt();
        overflow_state = Overflow::None;
        if scaling == Scaling::Unscaled {
            if z.re > MACHINE_CONSTANTS.approximation_limit {
                underflow_occurred = true;
                overflow_state = Overflow::NearUnder;
            } else {
                coeff *= MACHINE_CONSTANTS.scaling_factors[overflow_state] * (-z).exp();
            }
        }
        let mut AK = (signed_fractional_order * PI).cos().abs();
        let mut FHS = (0.25 - frac_order_sqr).abs();

        if signed_fractional_order.abs() == 0.5 || AK == 0.0 || FHS == 0.0 {
            (coeff, coeff)
        } else {
            //-----------------------------------------------------------------------
            //     MILLER ALGORITHM FOR CABS(Z) > R1;
            //-----------------------------------------------------------------------
            //-----------------------------------------------------------------------
            //     COMPUTE R2=F(E). if CABS(Z) >= R2, USE FORWARD RECURRENCE TO
            //     DETERMINE THE BACKWARD INDEX K. R2=F(E) IS A STRAIGHT LINE ON
            //     12 <= E <= 60. E IS COMPUTED FROM 2**(-E)=B**(1-i1mach(14))=
            //     TOL WHERE B IS THE BASE OF THE ARITHMETIC.
            //-----------------------------------------------------------------------
            let f64_significant_digits =
                (f64::MANTISSA_DIGITS - 1) as f64 * (f64::RADIX as f64).log10();
            let determiner = (f64_significant_digits * std::f64::consts::LOG2_10).clamp(12.0, 60.0);
            let recurrence_threshold = TWO_THIRDS * determiner - 6.0;
            let arg_z = z.arg();

            let (FK, FHS) = if abs_z > recurrence_threshold {
                //-----------------------------------------------------------------------;
                //     FORWARD RECURRENCE LOOP WHEN CABS(Z) >= R2;
                //-----------------------------------------------------------------------;
                let convergence_test = AK / (PI * abs_z * MACHINE_CONSTANTS.abs_error_tolerance);
                let mut FK = 1.0;
                if convergence_test >= 1.0 {
                    let mut FKS = 2.0;
                    let mut CKR = abs_z + abs_z + 2.0;
                    let mut p1 = 0.0;
                    let mut p2 = 1.0;
                    let mut converged = false;
                    for _ in 0..KMAX {
                        let AK = FHS / FKS;
                        let CBR = CKR / (FK + 1.0);
                        let pt = p2;
                        p2 = CBR * p2 - AK * p1;
                        p1 = pt;
                        CKR += 2.0;
                        FKS += FK + FK + 2.0;
                        FHS += FK + FK;
                        FK += 1.0;
                        if convergence_test < p2.abs() * FK {
                            converged = true;
                            break;
                        }
                    }
                    if !converged {
                        return Err(DidNotConverge);
                    }
                    FK += SPI * arg_z * (recurrence_threshold / abs_z).sqrt();
                    FHS = (0.25 - frac_order_sqr).abs();
                }
                (FK, FHS)
            } else {
                //-----------------------------------------------------------------------;
                //     COMPUTE BACKWARD INDEX K FOR CABS(Z) < R2;
                //-----------------------------------------------------------------------;
                AK *= FPI / (MACHINE_CONSTANTS.abs_error_tolerance * abs_z.sqrt().sqrt());
                let AA = 3.0 * arg_z / (1.0 + abs_z);
                let BB = 14.7 * arg_z / (28.0 + abs_z);
                AK = (AK.ln() + abs_z * AA.cos() / (1.0 + 0.008 * abs_z)) / BB.cos();
                let FK = 0.12125 * AK * AK / abs_z + 1.5;
                (FK, FHS)
            };
            //-----------------------------------------------------------------------;
            //     BACKWARD RECURRENCE LOOP FOR MILLER ALGORITHM;
            //-----------------------------------------------------------------------;
            let K = FK as usize;
            let mut k_squared = FK.floor().pow(2);
            let mut p1 = Complex64::zero();
            let mut p2 = Complex64::new(MACHINE_CONSTANTS.abs_error_tolerance, 0.0);
            let mut cs = p2;
            for i in (0..K).rev() {
                let k_f64 = (i + 1) as f64;
                let cb = (z + k_f64) * 2.0 / (k_f64 + 1.0);
                (p1, p2) = (
                    p2,
                    (p2 * cb - p1) * (k_squared + k_f64) / (k_squared - k_f64 + FHS),
                );
                cs += p2;
                k_squared -= (2.0 * k_f64) - 1.0;
            }
            //-----------------------------------------------------------------------;
            //     COMPUTE (P2/CS)=(P2/CABS(CS))*(CONJG(CS)/CABS(CS)) FOR BETTER;
            //     SCALING;
            //-----------------------------------------------------------------------;
            let mut s1 = p2 / cs.abs();
            let mut s2 = Complex64::zero();
            cs = cs.conj() / cs.abs();
            s1 *= coeff * cs;
            if !simple_case {
                //-----------------------------------------------------------------------;
                //     COMPUTE P1/P2=(P1/CABS(P2)*CONJG(P2)/CABS(P2) FOR SCALING;
                //-----------------------------------------------------------------------;
                p1 /= p2.abs();
                p2 = p2.conj() / p2.abs();
                s2 = (((signed_fractional_order + 0.5 - (p1 * p2)) / z) + 1.0) * s1;
            }
            (s1, s2)
        }
    };

    // Now s1, s2 set up, we can go to recurrence

    //-----------------------------------------------------------------------
    //     FORWARD RECURSION ON THE THREE TERM RECURSION WITH RELATION WITH
    //     SCALING NEAR EXPONENT EXTREMES ON KFLAG=1 OR KFLAG=3
    //-----------------------------------------------------------------------
    let mut ck = (signed_fractional_order + 1.0) * rz;
    if n == 1 {
        integer_order -= 1
    };

    if !simple_case {
        if integer_order > 0 {
            let mut n_tested = 1;
            if underflow_occurred {
                underflow_occurred = false;
                //-----------------------------------------------------------------------;
                //     underflow_occured=1 CASES, FORWARD RECURRENCE ON SCALED VALUES ON UNDERFLOW;
                //-----------------------------------------------------------------------;
                let mut cy = [c_zero(); 2];
                let half_exponent_limit = 0.5 * MACHINE_CONSTANTS.exponent_limit;

                let abs_limit = (-MACHINE_CONSTANTS.exponent_limit).exp();
                let ASCLE = MACHINE_CONSTANTS.smallness_threshold[0];
                let mut zd = z;
                let mut IC: isize = -1;
                let mut J = 1;
                for i in 0..integer_order {
                    n_tested = i + 2;
                    // TODO same calculation as other loops - this one is over different range and sets cy
                    // (so is designed to run until cy is set, and record this in INUB)
                    (s1, s2) = (s2, s2 * ck + s1);
                    ck += rz;
                    let abs_ln_s2 = s2.abs().ln();
                    if -zd.re + abs_ln_s2 >= -MACHINE_CONSTANTS.exponent_limit {
                        let p1 = (-zd + s2.ln()).exp() / MACHINE_CONSTANTS.abs_error_tolerance;
                        if !will_underflow(p1, ASCLE, MACHINE_CONSTANTS.abs_error_tolerance) {
                            J = 1 - J;
                            cy[J] = p1;
                            // below implies we got here twice in a row
                            if IC == i - 1 {
                                // underflow_occurred = true; //implies 270
                                break;
                            } else {
                                IC = i;
                                continue;
                            }
                        }
                        if abs_ln_s2 < half_exponent_limit {
                            continue;
                        }
                        zd.re -= MACHINE_CONSTANTS.exponent_limit;
                        s1 *= abs_limit;
                        s2 *= abs_limit;
                    }
                }
                overflow_state = Overflow::NearUnder;

                s2 = cy[J];
                J = 1 - J;
                s1 = cy[J];
            }

            let mut P1R = MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
            let mut ASCLE = MACHINE_CONSTANTS.smallness_threshold[overflow_state];
            for _ in n_tested..=integer_order {
                // TODO same loop as below?
                (s1, s2) = (s2, ck * s2 + s1);
                ck += rz;
                if overflow_state == Overflow::NearOver {
                    continue;
                }
                let p2 = s2 * P1R;
                if max_abs_component(p2) <= ASCLE {
                    continue;
                }
                overflow_state.increment();
                ASCLE = MACHINE_CONSTANTS.smallness_threshold[overflow_state];
                s1 *= P1R;
                s2 = p2;
                s1 *= MACHINE_CONSTANTS.scaling_factors[overflow_state];
                s2 *= MACHINE_CONSTANTS.scaling_factors[overflow_state];
                P1R = MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
            }
        }
        if n == 1 {
            s1 = s2;
        }
    }

    let mut y = c_zeros(n);
    let n_completed = if !underflow_occurred {
        // ********* basic setup
        y[0] = s1 * MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
        if n > 1 {
            y[1] = s2 * MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
            2
        } else {
            1
        }

        // ********* End Basic Setup
    } else {
        // ******** Alternative setup if underflow_occured

        y[0] = s1;
        if n > 1 {
            y[1] = s2;
        }
        ZKSCL(
            z,
            order,
            n,
            &mut y,
            &mut nz,
            rz,
            MACHINE_CONSTANTS.absolute_approximation_limit,
        );
        let n_non_zero = (n - nz) as isize;
        if n_non_zero <= 0 {
            return Ok((y, nz));
        }
        let mut working_index = nz;
        s1 = y[working_index];
        y[working_index] *= MACHINE_CONSTANTS.reciprocal_scaling_factors[0];
        if n_non_zero > 1 {
            // if n_non_zero == 1 {
            //     return Ok((y, nz));
            // }
            working_index += 1;
            s2 = y[working_index];
            y[working_index] *= MACHINE_CONSTANTS.reciprocal_scaling_factors[0];
        }
        if n_non_zero > 2 {
            ck = (order + (working_index as f64)) * rz;
            overflow_state = Overflow::NearUnder;
        }
        working_index + 1
    };
    // End Setup
    if n_completed >= n {
        return Ok((y, nz));
    }
    let mut P1R = MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
    let mut ASCLE = MACHINE_CONSTANTS.smallness_threshold[overflow_state];
    for y_elem in y.iter_mut().skip(n_completed) {
        // TODO same loops as above
        (s1, s2) = (s2, ck * s2 + s1);
        ck += rz;
        *y_elem = s2 * P1R;
        if overflow_state == Overflow::NearOver {
            continue;
        };
        if max_abs_component(*y_elem) <= ASCLE {
            continue;
        }
        overflow_state.increment();
        ASCLE = MACHINE_CONSTANTS.smallness_threshold[overflow_state];
        s1 *= P1R;
        s2 = *y_elem;
        s1 *= MACHINE_CONSTANTS.scaling_factors[overflow_state];
        s2 *= MACHINE_CONSTANTS.scaling_factors[overflow_state];
        P1R = MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
    }
    Ok((y, nz))
}

fn k_right_half_plane_helper(
    z: Complex64,
    frac_order_sqr: f64,
    signed_fractional_order: f64,
    mut f: Complex64,
    mut p: Complex64,
    mut q: Complex64,
    mut ck: Complex64,
) -> (Complex64, Complex64) {
    let mut a1 = 1.0;
    let cz_sqr_over_4 = 0.25 * z.powu(2);
    let abs_z = z.abs();
    let abs_z_sqr_over_4 = 0.25 * abs_z * abs_z;
    let mut ak = 1.0;
    let mut bk = 1.0 - frac_order_sqr;

    let mut s1 = f;
    let mut s2 = p;
    if abs_z >= MACHINE_CONSTANTS.abs_error_tolerance {
        while a1 > MACHINE_CONSTANTS.abs_error_tolerance {
            f = (f * ak + p + q) / bk;
            p /= ak - signed_fractional_order;
            q /= ak + signed_fractional_order;
            ck *= cz_sqr_over_4 / ak;
            s1 += ck * f;
            s2 += ck * (p - ak * f);
            a1 *= abs_z_sqr_over_4 / ak;
            bk += (2.0 * ak) + 1.0;
            ak += 1.0;
        }
    }
    (s1, s2)
}

/// Set k functions to zero on underflow, continue recurrence
/// on scaled functions until two members come on scale, then
/// return with min(nz+2,n) values scaled by 1/tol.
fn ZKSCL(
    zr: Complex64,
    order: f64,
    n: usize,
    y: &mut [Complex64],
    nz: &mut usize,
    rz: Complex64,
    ASCLE: f64,
) {
    *nz = 0;
    // let NN = min(2, n);
    let mut cy = [c_zero(); 2];
    let mut i_completed = 0;
    // repeats twice, unless n < 2
    for i in 0..min(2, n) {
        let s1 = y[i];
        cy[i] = s1;
        *nz += 1;
        y[i] = c_zero();
        if -zr.re + s1.abs().ln() < -MACHINE_CONSTANTS.exponent_limit {
            continue;
        }

        let cs = (s1.ln() - zr).exp() / MACHINE_CONSTANTS.abs_error_tolerance;
        if will_underflow(cs, ASCLE, MACHINE_CONSTANTS.abs_error_tolerance) {
            continue;
        }
        y[i] = cs;
        i_completed = i;
        *nz -= 1;
    }
    if n <= 2 || *nz == 0 {
        return;
    }
    // if i_completed < 1 {
    //     y[0] = c_zero();
    //     *nz = 2;
    // }
    // if n == 2 {
    //     return;
    // }
    // if *nz == 0 {
    //     return;
    // }
    let FN = order + 1.0;
    let mut ck = FN * rz;
    let mut s1 = cy[0];
    let mut s2 = cy[1];
    let half_elim = 0.5 * MACHINE_CONSTANTS.exponent_limit;
    let ELM = (-MACHINE_CONSTANTS.exponent_limit).exp();
    let CELMR = ELM;
    let mut zd = zr;
    //     FIND TWO CONSECUTIVE Y VALUES ON SCALE. SCALE RECURRENCE if
    //     S2 GETS LARGER THAN EXP(ELIM/2)
    let mut skip_to_40 = false;
    let mut I = 0;
    for i in 2..n {
        I = i;
        let mut cs = s2;
        s2 = cs * ck + s1;
        s1 = cs;
        ck += rz;
        let ALAS = s2.abs().ln();
        *nz += 1;
        y[i] = Complex64::zero();
        if -zd.re + s2.abs().ln() >= -MACHINE_CONSTANTS.exponent_limit {
            cs = s2.ln() - zd;
            cs = cs.exp() / MACHINE_CONSTANTS.abs_error_tolerance;
            if !will_underflow(cs, ASCLE, MACHINE_CONSTANTS.abs_error_tolerance) {
                y[i] = cs;
                *nz -= 1;
                if i_completed == i - 1 {
                    skip_to_40 = true;
                    break;
                }
                i_completed = i;
                continue;
            }
        }

        if ALAS < half_elim {
            continue;
        }
        zd -= MACHINE_CONSTANTS.exponent_limit;
        s1 *= CELMR;
        s2 *= CELMR;
    }
    if !skip_to_40 {
        *nz = n;
        if i_completed == n {
            *nz = n - 1
        };
    } else {
        *nz = I - 2;
    }
    for element in y.iter_mut().take(*nz) {
        *element = c_zero();
    }
}

/// ratios_i computes ratios of I bessel functions by backward
/// recurrence. The starting index is determined by forward
/// recurrence as described in J. Res. of Nat. Bur. of Standards-B,
/// Mathematical Sciences, vol 77b, p111-114, September, 1973,
/// Bessel functions I and J of complex argument and integer order,
/// by D. J. Sookne.
///
/// Originally ZRATI
fn ratios_i(z: Complex64, order: f64, n: usize) -> Vec<Complex64> {
    let abs_z = z.abs();
    let integer_order = order as usize;
    let modified_int_order = integer_order + n - 1;
    let int_abs_z = abs_z as isize;
    let FNUP = (int_abs_z + 1).max(modified_int_order as isize) as f64;
    let ID_ = modified_int_order as isize - int_abs_z - 1;
    let ID = if ID_ > 0 { 0 } else { ID_ };

    let rz = 2.0 * z.conj() / abs_z.powi(2);
    let mut K = 1;
    let mut abs_p2;
    {
        let mut t1 = rz * FNUP;
        let mut p2 = -t1;
        let mut p1 = c_one();
        t1 += rz;

        abs_p2 = p2.abs();
        let mut abs_p1 = p1.abs();
        //-----------------------------------------------------------------------
        //     THE OVERFLOW TEST ON K(FNU+I-1,Z) BEFORE THE CALL TO CBKNU
        //     GUARANTEES THAT P2 IS ON SCALE. SCALE TEST1 AND ALL SUBSEQUENT
        //     P2 VALUES BY AP1 TO ENSURE THAT AN OVERFLOW DOES NOT OCCUR
        //     PREMATURELY.
        //-----------------------------------------------------------------------
        let ARG = (abs_p2 + abs_p2) / (abs_p1 * MACHINE_CONSTANTS.abs_error_tolerance);
        let TEST1 = ARG.sqrt();
        let mut TEST = TEST1;
        p1 /= abs_p1;
        p2 /= abs_p1;
        abs_p2 /= abs_p1;
        let mut first_pass = true;
        'l10: loop {
            K += 1;
            abs_p1 = abs_p2;
            (p1, p2) = (p2, p1 - (t1 * p2));
            t1 += rz;
            abs_p2 = p2.abs();
            if abs_p1 <= TEST {
                continue;
            }
            if !first_pass {
                break 'l10;
            }
            {
                let ak = t1.abs() / 2.0;
                let flam = ak + (ak.powi(2) - 1.0).sqrt();
                let rho = abs_p2 / abs_p1.min(flam);
                TEST = TEST1 * (rho / (rho.powi(2) - 1.0)).sqrt();
            }
            first_pass = false;
        }
    }

    let mut p1 = Complex64::new(1.0 / abs_p2, 0.0);
    let mut p2 = c_zero();

    {
        let kk: usize = (K as isize + 1 - ID).try_into().unwrap();
        let mut t1 = Complex64::new(kk as f64, 0.0);
        let modified_order = order + ((n - 1) as f64);
        for _ in 0..kk {
            (p1, p2) = (p1 * (rz * (modified_order + t1.re)) + p2, p1);
            t1.re -= 1.0;
        }
        if p1.re == 0.0 && p1.im == 0.0 {
            p1 = Complex64::new(
                MACHINE_CONSTANTS.abs_error_tolerance,
                MACHINE_CONSTANTS.abs_error_tolerance,
            );
        }
    }
    let mut cy = c_zeros(n);
    cy[n - 1] = p2 / p1;
    if n > 1 {
        let mut t1 = Complex64::new((n - 1) as f64, 0.0);
        let cdfnu = order * rz;
        for k in (1..n).rev() {
            let mut pt = cdfnu + t1 * rz + cy[k];
            let mut abs_pt = pt.abs();
            if abs_pt == 0.0 {
                pt = Complex64::new(
                    MACHINE_CONSTANTS.abs_error_tolerance,
                    MACHINE_CONSTANTS.abs_error_tolerance,
                );
                abs_pt = pt.abs();
            }
            cy[k - 1] = pt.conj() / abs_pt.powi(2);
            t1 -= 1.0;
        }
    }
    cy
}

fn ZBUNK(
    z: Complex64,
    order: f64,
    scaling: Scaling,
    rotation: RotationDirection,
    N: usize,
) -> BesselResult {
    // ***BEGIN PROLOGUE  ZBUNK
    // ***REFER TO  ZBESK,ZBESH
    //
    //     ZBUNK COMPUTES THE K BESSEL FUNCTION FOR FNU > FNUL.
    //     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR K(FNU,Z)
    //     IN ZUNK1 AND THE EXPANSION FOR H(2,FNU,Z) IN ZUNK2
    //
    if imaginary_dominant(z) {
        //-----------------------------------------------------------------------
        //     ASYMPTOTIC EXPANSION FOR H(2,FNU,Z*EXP(M*FRAC_PI_2)) FOR LARGE FNU
        //     APPLIED IN PI/3 < ABS(ARG(Z)) <= PI/2 WHERE M=+I OR -I
        //     AND FRAC_PI_2=PI/2
        //-----------------------------------------------------------------------
        ZUNK2(z, order, scaling, rotation, N)
    } else {
        //-----------------------------------------------------------------------
        //     ASYMPTOTIC EXPANSION FOR K(FNU,Z) FOR LARGE FNU APPLIED IN
        //     -PI/3 <= ARG(Z) <= PI/3
        //-----------------------------------------------------------------------
        ZUNK1(z, order, scaling, rotation, N)
    }
}

fn i_miller(z: Complex64, order: f64, KODE: Scaling, N: usize) -> BesselResult {
    // ***BEGIN PROLOGUE  ZMLRI
    // ***REFER TO  ZBESI,ZBESK
    //
    //     ZMLRI COMPUTES THE I BESSEL FUNCTION FOR RE(Z) >= 0.0 BY THE
    //     MILLER ALGORITHM NORMALIZED BY A NEUMANN SERIES.
    //

    let SCLE: f64 = 2.0 * f64::MIN_POSITIVE / MACHINE_CONSTANTS.abs_error_tolerance;
    let NZ = 0;
    let AZ = z.abs();
    let IAZ = AZ as usize;
    let IFNU = order as usize;
    let INU = IFNU + N - 1;
    let AT = (IAZ as f64) + 1.0;
    let RAZ = 1.0 / AZ;
    let mut ck = z.conj() * RAZ * RAZ * AT;
    let rz = z.conj() * 2.0 * RAZ * RAZ;
    let mut p1 = c_zero();
    let mut p2 = c_one();
    let mut ACK = (AT + 1.0) * RAZ;
    let RHO = ACK + (ACK * ACK - 1.0).sqrt();
    let RHO2 = RHO * RHO;
    let mut TST = (RHO2 + RHO2) / ((RHO2 - 1.0) * (RHO - 1.0));
    TST /= MACHINE_CONSTANTS.abs_error_tolerance;
    //-----------------------------------------------------------------------
    //     COMPUTE RELATIVE TRUNCATION ERROR INDEX FOR SERIES
    //-----------------------------------------------------------------------
    let mut AK = AT;
    let mut converged = false;
    let mut I = 0;
    for i in 0..80 {
        I = i + 1;
        let pt = p2;
        p2 = p1 - ck * p2;
        p1 = pt;
        ck += rz;
        if p2.abs() > TST * AK * AK {
            converged = true;
            break;
        }
        AK += 1.0;
    }
    if !converged {
        return Err(DidNotConverge);
    }
    I += 1;
    let mut K = 0;
    if INU >= IAZ {
        //-----------------------------------------------------------------------
        //     COMPUTE RELATIVE TRUNCATION ERROR FOR RATIOS
        //-----------------------------------------------------------------------
        p1 = c_zero();
        p2 = c_one();
        let AT = (INU as f64) + 1.0;
        ck = z.conj() * RAZ * RAZ * AT;
        ACK = AT * RAZ;
        TST = (ACK / MACHINE_CONSTANTS.abs_error_tolerance).sqrt();
        let mut hit_loop_end = false;
        converged = false;
        for k in 0..80 {
            K = k;
            let pt = p2;
            p2 = p1 - ck * pt;
            p1 = pt;
            ck += rz;
            let AP = p2.abs();
            if AP < TST {
                continue;
            }
            if hit_loop_end {
                converged = true;
                break;
            }
            ACK = ck.abs();
            let FLAM = ACK + (ACK * ACK - 1.0).sqrt();
            let FKAP = AP / p1.abs();
            let RHO = FLAM.min(FKAP);
            TST *= (RHO / (RHO * RHO - 1.0)).sqrt();
            hit_loop_end = true;
        }
        if !converged {
            return Err(DidNotConverge);
        }
    }
    //-----------------------------------------------------------------------
    //     BACKWARD RECURRENCE AND SUM NORMALIZING RELATION
    //-----------------------------------------------------------------------
    K += 1;
    let KK = (I + IAZ).max(K + INU);
    let mut FKK = KK as f64;
    let mut p1 = c_zero();
    //-----------------------------------------------------------------------
    //     SCALE P2 AND SUM BY SCLE
    //-----------------------------------------------------------------------
    let mut p2 = Complex64::new(SCLE, 0.0);
    let FNF = order - (IFNU as f64);
    let TFNF = FNF + FNF;
    let mut BK = (gamma_ln(FKK + TFNF + 1.0).unwrap()
        - gamma_ln(FKK + 1.0).unwrap()
        - gamma_ln(TFNF + 1.0).unwrap())
    .exp();
    let mut sumr = c_zero();
    for _ in 0..(KK - INU) {
        let pt = p2;
        p2 = p1 + (FKK + FNF) * (rz * p2);
        p1 = pt;
        AK = 1.0 - TFNF / (FKK + TFNF);
        ACK = BK * AK;
        sumr += (ACK + BK) * p1;
        BK = ACK;
        FKK -= 1.0;
    }
    let mut y = c_zeros(N);
    y[N - 1] = p2;
    if N != 1 {
        for i in 1..N {
            let pt = p2;
            p2 = p1 + (FKK + FNF) * (rz * pt);
            p1 = pt;
            AK = 1.0 - TFNF / (FKK + TFNF);
            ACK = BK * AK;
            sumr += (ACK + BK) * p1;
            BK = ACK;
            FKK -= 1.0;
            y[N - (i + 1)] = p2;
        }
    }
    if IFNU > 0 {
        for _i in 0..IFNU {
            let pt = p2;
            p2 = p1 + (FKK + FNF) * (rz * pt);
            p1 = pt;
            AK = 1.0 - TFNF / (FKK + TFNF);
            ACK = BK * AK;
            sumr += (ACK + BK) * p1;
            BK = ACK;
            FKK -= 1.0;
        }
    }

    let mut pt = z;
    if KODE == Scaling::Scaled {
        pt.re = 0.0;
    }
    p1 = -FNF * rz.ln() + pt;
    let AP = gamma_ln(1.0 + FNF).unwrap();
    p1 -= AP;
    //-----------------------------------------------------------------------
    //     THE DIVISION CEXP(PT)/(SUM+P2) IS ALTERED TO AVOID OVERFLOW
    //     IN THE DENOMINATOR BY SQUARING LARGE QUANTITIES
    //-----------------------------------------------------------------------
    p2 += sumr;
    let AP = p2.abs();
    ck = p1.exp() / AP;
    let cnorm = ck * p2.conj() / AP;
    for element in y.iter_mut() {
        *element *= cnorm;
    }
    Ok((y, NZ))
}

fn i_wronksian(
    zr: Complex64,
    order: f64,
    KODE: Scaling,
    N: usize,
    y: &mut [Complex64],
) -> BesselResult<usize> {
    // ***BEGIN PROLOGUE  ZWRSK
    // ***REFER TO  ZBESI,ZBESK
    //
    //     ZWRSK COMPUTES THE I BESSEL FUNCTION FOR RE(Z) >= 0.0 BY
    //     NORMALIZING THE I FUNCTION RATIOS FROM ZRATI BY THE WRONSKIAN
    //

    //-----------------------------------------------------------------------
    //     I(FNU+I-1,Z) BY BACKWARD RECURRENCE FOR RATIOS
    //     Y(I)=I(FNU+I,Z)/I(FNU+I-1,Z) FROM CRATI NORMALIZED BY THE
    //     WRONSKIAN WITH K(FNU,Z) AND K(FNU+1,Z) FROM CBKNU.
    //-----------------------------------------------------------------------
    let NZ = 0;
    let (cw, _) = k_right_half_plane(zr, order, KODE, 2)?;
    let y_ratios = ratios_i(zr, order, N);
    //-----------------------------------------------------------------------
    //     RECUR FORWARD ON I(FNU+1,Z) = R(FNU,Z)*I(FNU,Z),
    //     R(FNU+J-1,Z)=Y(J),  J=1,...,N
    //-----------------------------------------------------------------------
    let mut cinu = c_one();
    if KODE == Scaling::Scaled {
        cinu = Complex64::cis(zr.im);
    }
    //-----------------------------------------------------------------------
    //     ON LOW EXPONENT MACHINES THE K FUNCTIONS CAN BE CLOSE TO BOTH
    //     THE UNDER AND OVERFLOW LIMITS AND THE NORMALIZATION MUST BE
    //     SCALED TO PREVENT OVER OR UNDERFLOW. CUOIK HAS DETERMINED THAT
    //     THE RESULT IS ON SCALE.
    //-----------------------------------------------------------------------
    let acw = cw[1].abs();
    let CSCLR = if acw <= MACHINE_CONSTANTS.absolute_approximation_limit {
        1.0 / MACHINE_CONSTANTS.abs_error_tolerance
    } else if acw >= 1.0 / MACHINE_CONSTANTS.absolute_approximation_limit {
        MACHINE_CONSTANTS.abs_error_tolerance
    } else {
        1.0
    };
    let c1 = cw[0] * CSCLR;
    let c2 = cw[1] * CSCLR;
    //-----------------------------------------------------------------------
    //     CINU=CINU*(CONJG(CT)/CABS(CT))*(1.0/CABS(CT) PREVENTS
    //     UNDER- OR OVERFLOW PREMATURELY BY SQUARING CABS(CT)
    //-----------------------------------------------------------------------
    let mut ct = zr * (y_ratios[0] * c1 + c2);
    let ct_abs = ct.abs();
    ct = ct.conj() / ct_abs;
    cinu = (cinu / ct_abs) * ct;
    y[0] = cinu * CSCLR;
    for i in 1..N {
        cinu *= y_ratios[i - 1];
        y[i] = cinu * CSCLR;
    }
    Ok(NZ)
}

fn analytic_continuation(
    z: Complex64,
    order: f64,
    scaling: Scaling,
    rotation: RotationDirection,
    N: usize,
) -> BesselResult {
    // ***BEGIN PROLOGUE  ZACON
    // ***REFER TO  ZBESK,ZBESH
    //
    //     ZACON APPLIES THE ANALYTIC CONTINUATION FORMULA
    //
    //         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
    //                 MP=PI*MR*CMPLX(0.0,1.0)
    //
    //     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT
    //     HALF Z PLANE
    //

    let mut NZ = 0;
    let zn = -z;
    let (mut y, _) = i_right_half_plane(zn, order, scaling, N)?;
    //-----------------------------------------------------------------------
    //     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
    //-----------------------------------------------------------------------
    let NN = 2.min(N);
    let (cy, NW) = k_right_half_plane(zn, order, scaling, NN)?;
    if NW > 0 {
        return Err(Overflow);
        // the NW = -1 or -2 is handled by ZBNKU returning an error,
        // but the amos code defaults to an overflow, if NW != 0
    }
    let mut s1 = cy[0];
    let SGN = -PI * rotation.signum();
    let mut csgn = Complex64::new(0.0, SGN);
    if scaling == Scaling::Scaled {
        csgn *= Complex64::cis(-zn.im);
    }
    //-----------------------------------------------------------------------
    //     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
    //     WHEN FNU IS LARGE
    //-----------------------------------------------------------------------
    let mut cspn = Complex64::cis(order.fract() * SGN);
    if order as i64 % 2 != 0 {
        cspn = -cspn;
    }
    let mut IUF = 0;
    let mut c1 = s1;
    let mut c2 = y[0];
    let mut sc1;
    if scaling == Scaling::Scaled {
        let NW = underflow_add_i_k(zn, &mut c1, &mut c2, &mut IUF);
        NZ += NW;
    }
    let st = cspn * c1;
    let pt = csgn * c2;
    y[0] = st + pt;
    if N == 1 {
        return Ok((y, NZ));
    }
    cspn = -cspn;
    let mut s2 = cy[1];
    c1 = s2;
    c2 = y[1];
    // this value never used, as initialised and used if scaling is needed
    let mut sc2 = c_zero() * f64::NAN;
    if scaling == Scaling::Scaled {
        let NW = underflow_add_i_k(zn, &mut c1, &mut c2, &mut IUF);
        NZ += NW;
        sc2 = c1;
    }
    let st = cspn * c1;
    let pt = csgn * c2;
    y[1] = st + pt;

    if N == 2 {
        return Ok((y, NZ));
    }
    cspn = -cspn;
    let rz = 2.0 * (zn.conj()) / zn.abs().powi(2);
    let FN = order + 1.0;
    let mut ck = FN * rz;
    //-----------------------------------------------------------------------
    //     SCALE NEAR EXPONENT EXTREMES DURING RECURRENCE ON K FUNCTIONS
    //-----------------------------------------------------------------------
    let abs_s2 = s2.abs();
    let mut overflow_state = if abs_s2 <= MACHINE_CONSTANTS.smallness_threshold[0] {
        Overflow::NearUnder
    } else if abs_s2 > MACHINE_CONSTANTS.smallness_threshold[1] {
        Overflow::NearOver
    } else {
        Overflow::None
    };
    let mut b_scale = MACHINE_CONSTANTS.smallness_threshold[overflow_state];
    s1 *= MACHINE_CONSTANTS.scaling_factors[overflow_state];
    s2 *= MACHINE_CONSTANTS.scaling_factors[overflow_state];
    let mut CSR = MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
    for yi in y.iter_mut().skip(2) {
        //TODO common pattern below
        (s1, s2) = (s2, ck * s2 + s1);
        c1 = s2 * CSR;
        let mut st = c1;
        c2 = *yi;
        if scaling == Scaling::Scaled && IUF >= 0 {
            let NW = underflow_add_i_k(zn, &mut c1, &mut c2, &mut IUF);
            NZ += NW;
            sc1 = sc2;
            sc2 = c1;
            if IUF == 3 {
                IUF = -4;
                s1 = sc1 * MACHINE_CONSTANTS.scaling_factors[overflow_state];
                s2 = sc2 * MACHINE_CONSTANTS.scaling_factors[overflow_state];
                st = sc2;
            }
        }
        *yi = cspn * c1 + csgn * c2;
        ck += rz;
        cspn = -cspn;
        if overflow_state == Overflow::NearOver {
            continue;
        }
        if max_abs_component(c1) > b_scale {
            overflow_state.increment();
            b_scale = MACHINE_CONSTANTS.smallness_threshold[overflow_state];
            s1 *= CSR;
            s2 = st;
            s1 *= MACHINE_CONSTANTS.scaling_factors[overflow_state];
            s2 *= MACHINE_CONSTANTS.scaling_factors[overflow_state];
            CSR = MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state]; //CSRR(KFLAG);
        }
    }
    Ok((y, NZ))
}

/// i_right_half_plane computes the i function in the right half z plane
/// Originally ZBINU
fn i_right_half_plane(z: Complex64, order: f64, KODE: Scaling, N: usize) -> BesselResult {
    let mut NZ = 0;
    let AZ = z.abs();
    let mut NN: usize = N;
    let mut DFNU = order + ((N - 1) as f64);
    let mut cy = c_zeros(N);
    if AZ <= 2.0 || AZ * AZ * 0.25 <= DFNU + 1.0 {
        //-----------------------------------------------------------------------
        //     POWER SERIES
        //-----------------------------------------------------------------------
        let NW;
        (cy, NW) = i_power_series(z, order, KODE, NN)?;
        let INW: usize = NW.unsigned_abs();
        NZ += INW;
        NN -= INW;
        if NN == 0 || NW >= 0 {
            return Ok((cy, NZ));
        }

        DFNU = order + ((NN as f64) - 1.0);
    }

    if (AZ >= MACHINE_CONSTANTS.asymptotic_z_limit) && ((DFNU <= 1.0) || (AZ + AZ >= DFNU * DFNU)) {
        //-----------------------------------------------------------------------
        //     ASYMPTOTIC EXPANSION FOR LARGE Z
        //-----------------------------------------------------------------------
        let (cy, nw) = asymptotic_i(z, order, KODE, NN)?;
        debug_assert!(nw == NZ);
        return Ok((cy, NZ));
    }
    let mut skip_az_rl_check = true;
    if DFNU > 1.0 {
        skip_az_rl_check = false;
        //-----------------------------------------------------------------------
        //     OVERFLOW AND UNDERFLOW TEST ON I SEQUENCE FOR MILLER ALGORITHM
        //-----------------------------------------------------------------------
        let nw = check_underflow_uniform_asymp_params(z, order, KODE, IKType::I, NN, &mut cy)?;
        NZ += nw;
        NN -= nw;
        if NN == 0 {
            return Ok((cy, NZ));
        }
        DFNU = order + ((NN - 1) as f64);
    }
    if (DFNU > MACHINE_CONSTANTS.asymptotic_order_limit)
        || (AZ > MACHINE_CONSTANTS.asymptotic_order_limit)
    {
        //-----------------------------------------------------------------------
        //     INCREMENT FNU+NN-1 UP TO FNUL, COMPUTE AND RECUR BACKWARD
        //-----------------------------------------------------------------------
        let NUI_isize = (MACHINE_CONSTANTS.asymptotic_order_limit - DFNU) as isize + 1;
        let NUI = NUI_isize.max(0) as usize;
        let (NW, NLAST) = ZBUNI(z, order, KODE, NN, NUI, &mut cy)?;
        NZ += NW;
        if NLAST == 0 {
            return Ok((cy, NZ));
        }
        NN = NLAST;
    }
    if !skip_az_rl_check && AZ <= MACHINE_CONSTANTS.asymptotic_z_limit {
        //-----------------------------------------------------------------------
        //     MILLER ALGORITHM NORMALIZED BY THE SERIES
        //-----------------------------------------------------------------------
        let (cy, _) = i_miller(z, order, KODE, NN)?;
        return Ok((cy, NZ)); //}
    }
    //-----------------------------------------------------------------------
    //     MILLER ALGORITHM NORMALIZED BY THE WRONSKIAN
    //-----------------------------------------------------------------------
    //-----------------------------------------------------------------------
    //     OVERFLOW TEST ON K FUNCTIONS USED IN WRONSKIAN
    //-----------------------------------------------------------------------
    if let Ok(NW) =
        check_underflow_uniform_asymp_params(z, order, KODE, IKType::K, 2, &mut [c_one(); 2])
    {
        if NW > 0 {
            Err(Overflow)
        } else {
            let nz = i_wronksian(z, order, KODE, NN, &mut cy)?;
            Ok((cy, nz))
        }
    } else {
        Ok((vec![c_one(); NN], NN))
    }
}

fn ZACAI(
    z: Complex64,
    order: f64,
    KODE: Scaling,
    rotation: RotationDirection,
    N: usize,
) -> BesselResult {
    // ***BEGIN PROLOGUE  ZACAI
    // ***REFER TO  ZAIRY
    //
    //     ZACAI APPLIES THE ANALYTIC CONTINUATION FORMULA
    //
    //         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
    //                 MP=PI*MR*CMPLX(0.0,1.0)
    //
    //     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT
    //     HALF Z PLANE FOR USE WITH ZAIRY WHERE FNU=1/3 OR 2/3 AND N=1.
    //     ZACAI IS THE SAME AS ZACON WITH THE PARTS FOR LARGER ORDERS AND
    //     RECURRENCE REMOVED. A RECURSIVE CALL TO ZACON CAN RESULT if ZACON
    //     IS CALLED FROM ZAIRY.
    //

    let mut NZ = 0;
    let zn = -z;
    let AZ = z.abs();
    let NN = N;
    let DFNU = order + ((N - 1) as f64);
    let (mut y, _) = if (AZ * AZ * 0.25 <= DFNU + 1.0) || (AZ <= 2.0) {
        //-----------------------------------------------------------------------
        //     POWER SERIES FOR THE I FUNCTION
        //-----------------------------------------------------------------------
        let (y, NW_signed) = i_power_series(zn, order, KODE, NN)?;
        debug_assert!(NW_signed >= 0);
        (y, NW_signed.unsigned_abs())
    } else if AZ >= MACHINE_CONSTANTS.asymptotic_z_limit {
        //-----------------------------------------------------------------------
        //     ASYMPTOTIC EXPANSION FOR LARGE Z FOR THE I FUNCTION
        //-----------------------------------------------------------------------
        asymptotic_i(zn, order, KODE, NN)?
    //-----------------------------------------------------------------------
    //     MILLER ALGORITHM NORMALIZED BY THE SERIES FOR THE I FUNCTION
    //-----------------------------------------------------------------------
    } else {
        i_miller(zn, order, KODE, NN)?
    };
    //-----------------------------------------------------------------------
    //     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
    //-----------------------------------------------------------------------s
    let (cy, nz) = k_right_half_plane(zn, order, KODE, 1)?;
    if nz != 0 {
        return Err(Overflow);
    }
    let SGN = -PI * rotation.signum();
    let mut csgn = Complex64::new(0.0, SGN);
    if KODE == Scaling::Scaled {
        csgn = Complex64::I * csgn.im * Complex64::cis(-zn.im);
    }
    //-----------------------------------------------------------------------
    //     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
    //     WHEN FNU IS LARGE
    //-----------------------------------------------------------------------
    let INU = order as usize;
    let mut cspn = Complex64::cis(order.fract() * SGN);
    if !INU.is_multiple_of(2) {
        cspn = -cspn;
    }
    let mut c1 = cy[0];
    let mut c2 = y[0];
    if KODE == Scaling::Scaled {
        let mut IUF = 0;
        let NW = underflow_add_i_k(zn, &mut c1, &mut c2, &mut IUF);
        NZ += NW;
    }
    y[0] = cspn * c1 + csgn * c2;
    Ok((y, NZ))
}

fn ZUNK1(
    z: Complex64,
    order: f64,
    scaling: Scaling,
    rotation: RotationDirection,
    N: usize,
) -> BesselResult {
    //     ZUNK1 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE
    //     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE
    //     UNIFORM ASYMPTOTIC EXPANSION.
    //     MR INDICATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION.

    // TODO better name for KDFLAG
    let mut KDFLG = false; // = 1;
    let mut NZ = 0;
    //-----------------------------------------------------------------------
    //     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION GREATER THAN
    //     THE UNDERFLOW LIMIT
    //-----------------------------------------------------------------------
    let zr = if z.re < 0.0 { -z } else { z };
    let mut J = 1;
    let mut phi = [c_zero(); 2];
    let mut zeta1 = [c_zero(); 2];
    let mut zeta2 = [c_zero(); 2];
    let mut sum = [c_zero(); 2];
    let mut cy = [c_zero(); 2];
    let mut I = 0;
    let mut modified_order = 0.0;
    let mut y = c_zeros(N);
    let mut k_overflow_state = Overflow::NearUnder;

    for i in 0..N {
        I = i;
        //-----------------------------------------------------------------------
        //     J FLIP FLOPS BETWEEN 1 AND 2 IN J = 3 - J
        //     In Rust this is 0 and 1
        //-----------------------------------------------------------------------
        J = 1 - J;
        modified_order = order + (i as f64);

        let sum_opt;
        (phi[J], zeta1[J], zeta2[J], sum_opt) = zunik(zr, modified_order, IKType::K, false);
        sum[J] = sum_opt.unwrap();
        let mut s1 = -scaling.scale_zetas(zr, modified_order, zeta1[J], zeta2[J]);
        let of = Overflow::find_overflow(s1.re, phi[J], 0.0);
        if !KDFLG {
            k_overflow_state = of;
        }
        match of {
            Overflow::Over(_) => return Err(Overflow),
            Overflow::Under(_) => {
                if z.re < 0.0 {
                    return Err(Overflow);
                }
                KDFLG = false;
                y[i] = c_zero();
                NZ += 1;
            }
            Overflow::None | Overflow::NearOver | Overflow::NearUnder => {
                //-----------------------------------------------------------------------
                //     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
                //     EXPONENT EXTREMES
                //-----------------------------------------------------------------------
                let mut s2 = phi[J] * sum[J];
                s1 = MACHINE_CONSTANTS.scaling_factors[k_overflow_state] * s1.exp();
                s2 *= s1;
                let will_underflow = will_underflow(
                    s2,
                    MACHINE_CONSTANTS.smallness_threshold[0],
                    MACHINE_CONSTANTS.abs_error_tolerance,
                );
                if k_overflow_state != Overflow::NearUnder || !will_underflow {
                    cy[KDFLG as usize] = s2;
                    y[i] = s2 * MACHINE_CONSTANTS.reciprocal_scaling_factors[k_overflow_state];
                    if KDFLG {
                        break;
                    }
                    KDFLG = true;
                } else if will_underflow {
                    if z.re < 0.0 {
                        return Err(Overflow);
                    }
                    y[i] = c_zero();
                    NZ += 1;
                    if i > 0 && y[i - 1] != c_zero() {
                        y[i - 1] = c_zero();
                        NZ += 1
                    }
                }
            }
        };
    }
    if !KDFLG {
        I = N;
    }
    let rz = 2.0 * zr.conj() / zr.abs().powi(2);
    let mut ck = modified_order * rz;
    let IB = I + 1;
    if N > IB {
        //-----------------------------------------------------------------------
        //     TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW. SET SEQUENCE TO ZERO
        //     ON UNDERFLOW.
        //-----------------------------------------------------------------------
        modified_order = order + ((N - 1) as f64);
        let (phi, zet1d, zet2d, _sumd) = zunik(
            zr,
            modified_order,
            IKType::K,
            rotation == RotationDirection::None,
        );
        let overflow_test = -scaling.scale_zetas(zr, modified_order, zet1d, zet2d);

        match Overflow::find_overflow(overflow_test.re.abs(), phi, 0.0) {
            Overflow::Over(_) => return Err(Overflow),
            Overflow::Under(_) => {
                return if z.re < 0.0 {
                    Err(Overflow)
                } else {
                    Ok((vec![c_zero(); N], N))
                };
            }
            _ => (),
        }
        //---------------------------------------------------------------------------
        //     FORWARD RECUR FOR REMAINDER OF THE SEQUENCE
        //----------------------------------------------------------------------------
        let [mut s1, mut s2] = cy;
        let mut C1R = MACHINE_CONSTANTS.reciprocal_scaling_factors[k_overflow_state];
        let mut ASCLE = MACHINE_CONSTANTS.smallness_threshold[k_overflow_state];
        for item in y[IB..N].iter_mut() {
            (s1, s2) = (s2, ck * s2 + s1);
            ck += rz;
            *item = s2 * C1R;
            if k_overflow_state == Overflow::NearOver {
                continue;
            }
            if max_abs_component(*item) <= ASCLE {
                continue;
            }
            k_overflow_state.increment();
            ASCLE = MACHINE_CONSTANTS.smallness_threshold[k_overflow_state];
            s1 *= C1R;
            s2 = *item;
            s1 *= MACHINE_CONSTANTS.scaling_factors[k_overflow_state];
            s2 *= MACHINE_CONSTANTS.scaling_factors[k_overflow_state];
            C1R = MACHINE_CONSTANTS.reciprocal_scaling_factors[k_overflow_state];
        }
        if rotation == RotationDirection::None {
            return Ok((y, NZ));
        }
        //-----------------------------------------------------------------------
        //     ANALYTIC CONTINUATION FOR RE(Z) < 0.0
        //-----------------------------------------------------------------------
        NZ = 0;
        let SGN = -PI * rotation.signum();
        //-----------------------------------------------------------------------
        //     CSPN AND CSGN ARE COEFF OF K AND I FUNCTIONS RESP.
        //-----------------------------------------------------------------------
        let CSGNI = SGN;

        let INU = order as i64;
        let order_frac = order.fract();
        let IFN = INU + (N as i64) - 1;
        let ANG = order_frac * SGN;
        let mut cspn = Complex64::cis(ANG);
        if (IFN % 2) != 0 {
            cspn = -cspn;
        }
        let mut IUF = 0;
        KDFLG = false;
        let mut i_overflow_state = Overflow::None;
        let mut left_early = false;
        let mut completed_k = 0;
        for k in (0..N).rev() {
            completed_k = N - k;
            modified_order = order + (k as f64);
            //-----------------------------------------------------------------------
            //     LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K
            //     FUNCTION ABOVE
            //-----------------------------------------------------------------------
            let (phid, zet1d, zet2d, sumd) = zunik(zr, modified_order, IKType::I, false); //, &mut INITD);
            let sumd = sumd.unwrap();
            let mut s1 = scaling.scale_zetas(zr, modified_order, zet1d, zet2d);
            //-----------------------------------------------------------------------
            //     TEST FOR UNDERFLOW AND OVERFLOW
            //-----------------------------------------------------------------------
            let of = Overflow::find_overflow(s1.re, phid, 0.0);
            if !KDFLG && !matches!(of, Overflow::Under(_)) {
                i_overflow_state = of;
            }
            let mut s2 = match of {
                Overflow::Over(_) => {
                    return Err(Overflow);
                }
                Overflow::Under(_) => c_zero(),
                Overflow::NearOver | Overflow::NearUnder | Overflow::None => {
                    let st = phid * sumd;
                    let mut s2 = Complex64::I * st * CSGNI;
                    s1 = s1.exp() * MACHINE_CONSTANTS.scaling_factors[i_overflow_state];
                    s2 *= s1;
                    if i_overflow_state == Overflow::NearUnder
                        && will_underflow(
                            s2,
                            MACHINE_CONSTANTS.smallness_threshold[0],
                            MACHINE_CONSTANTS.abs_error_tolerance,
                        )
                    {
                        s2 = c_zero();
                    }
                    s2
                }
            };
            cy[KDFLG as usize] = s2;
            let c2 = s2;
            s2 *= MACHINE_CONSTANTS.reciprocal_scaling_factors[i_overflow_state];
            //-----------------------------------------------------------------------
            //     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N
            //-----------------------------------------------------------------------
            s1 = y[k];
            if scaling == Scaling::Scaled {
                let NW = underflow_add_i_k(zr, &mut s1, &mut s2, &mut IUF);
                nz += NW;
            }
            y[k] = s1 * cspn + s2;
            cspn = -cspn;
            if c2 == c_zero() {
                KDFLG = false;
                continue;
            }
            if KDFLG {
                left_early = true;
            }
            KDFLG = true;
        }
        if !left_early {
            completed_k = N;
        }
        let IL = N - completed_k;
        if IL == 0 {
            return Ok((y, NZ));
        }
        //-----------------------------------------------------------------------
        //     RECUR BACKWARD FOR REMAINDER OF I SEQUENCE AND ADD IN THE
        //     K FUNCTIONS, SCALING THE I SEQUENCE DURING RECURRENCE TO KEEP
        //     INTERMEDIATE ARITHMETIC ON SCALE NEAR EXPONENT EXTREMES.
        //-----------------------------------------------------------------------
        let mut s1 = cy[0];
        let mut s2 = cy[1];
        let mut csr = MACHINE_CONSTANTS.scaling_factors[i_overflow_state];
        let mut ASCLE = MACHINE_CONSTANTS.smallness_threshold[i_overflow_state];
        modified_order = INU as f64 + IL as f64;
        let mut KK = N - completed_k;
        for _ in 0..IL {
            let mut c2 = s2;
            s2 = s1 + (modified_order * order_frac) * (rz * c2);
            s1 = c2;
            modified_order -= 1.0;
            c2 = s2 * csr;
            let ck = c2;

            let mut c1 = y[KK];
            if scaling == Scaling::Scaled {
                let NW = ZS1S2(zr, &mut c1, &mut c2, &mut IUF);
                NZ += NW;
            }
            y[KK] = c1 * cspn + c2;
            KK -= 1;
            cspn = -cspn;
            if i_overflow_state == Overflow::NearOver {
                continue;
            }
            if max_abs_component(c2) <= ASCLE {
                continue;
            }
            i_overflow_state.increment();
            ASCLE = MACHINE_CONSTANTS.smallness_threshold[i_overflow_state];
            s1 *= csr;
            s2 = ck;
            s1 *= MACHINE_CONSTANTS.scaling_factors[i_overflow_state];
            s2 *= MACHINE_CONSTANTS.scaling_factors[i_overflow_state];
            csr = MACHINE_CONSTANTS.scaling_factors[i_overflow_state];
        }
    }
    Ok((y, NZ))
}

fn airy_pair(z: Complex64) -> (Complex64, Complex64) {
    //note that ZAIRY calls in fortran code ignore IERR (using IDUM)
    let airy = match complex_airy(z, false, Scaling::Scaled) {
        Ok((y, _)) => y,
        Err(PartialLossOfSignificance { y, nz: _ }) => y[0],
        // If loss of significance, Fortran code would continue with un-initialised y,
        // which is usually ~=0. Also long as it is << d_airy, the logic below means
        // it will not matter what the precise value is
        Err(LossOfSignificance) => c_zero(),
        Err(err) => {
            panic!(
                "An error {:?} was generated, which is not handled by the Amos code",
                err
            )
        }
    };
    let d_airy = match complex_airy(z, true, Scaling::Scaled) {
        Ok((y, _)) => y,
        Err(PartialLossOfSignificance { y, nz: _ }) => y[0],
        // If loss of significance, Fortran code would continue with un-initialised y,
        // which is usually ~=0. Also long as it is << a_airy, the logic below means
        // it will not matter what the precise value is
        Err(LossOfSignificance) => c_zero(),
        Err(err) => {
            panic!(
                "An error {:?} was generated, which is not handled by the Amos code",
                err
            )
        }
    };
    (airy, d_airy)
}

const CIP: [Complex64; 4] = [
    Complex64::new(1.0, 0.0),
    Complex64::new(0.0, 1.0),
    Complex64::new(-1.0, 0.0),
    Complex64::new(0.0, -1.0),
];

fn ZUNK2(
    z: Complex64,
    order: f64,
    scaling: Scaling,
    rotation: RotationDirection,
    n: usize,
) -> BesselResult {
    //     ZUNK2 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE
    //     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE
    //     UNIFORM ASYMPTOTIC EXPANSIONS FOR H(KIND,FNU,ZN) AND J(FNU,ZN)
    //     WHERE ZN IS IN THE RIGHT HALF PLANE, KIND=(3-MR)/2, MR=+1 OR
    //     -1. HERE ZN=ZR*I OR -ZR*I WHERE ZR=Z if Z IS IN THE RIGHT
    //     HALF PLANE OR ZR=-Z if Z IS IN THE LEFT HALF PLANE. MR INDIC-
    //     ATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION.

    const CR1: Complex64 = Complex64::new(1.0, 1.73205080756887729);
    const CR2: Complex64 = Complex64::new(-0.5, -8.66025403784438647e-01);

    let mut underflowed_already = false;
    let mut nz = 0;
    let mut y = c_zeros(n);
    let zr = if z.re < 0.0 { -z } else { z };
    let mut zn = -Complex64::I * zr;
    let mut zb = zr;
    let integer_order = order as usize;
    let order_fract = order.fract();
    let ANG = -FRAC_PI_2 * order_fract;
    let mut c2 = -Complex64::I * Complex64::from_polar(FRAC_PI_2, ANG);
    let mut cs = CR1 * c2 * CIP[integer_order % 4].conj();
    if zr.im <= 0.0 {
        zn.re = -zn.re;
        zb.im = -zb.im;
    }
    //-----------------------------------------------------------------------
    //     K(FNU,Z) IS COMPUTED FROM H(2,FNU,-I*Z) WHERE Z IS IN THE FIRST
    //     QUADRANT. FOURTH QUADRANT VALUES (YY <= 0.0E0) ARE COMPUTED BY
    //     CONJUGATION SINCE THE K FUNCTION IS REAL ON THE POSITIVE REAL AXIS
    //-----------------------------------------------------------------------

    let mut phi = [c_zero(); 2];
    let mut arg = [c_zero(); 2];
    let mut zeta1 = [c_zero(); 2];
    let mut zeta2 = [c_zero(); 2];
    let mut asum = [None; 2];
    let mut bsum = [None; 2];
    let mut cy = [c_zero(); 2];
    let mut j = 1;
    let mut overflow_state_k = Overflow::None;
    let mut modified_order = 0.0;
    let mut n_elements_set = 0;

    for i in 0..n {
        n_elements_set = i + 1;
        // J flip-flops between 0 and 1 using J = 1-J
        j = 1 - j;
        modified_order = order + (i as f64);
        (phi[j], arg[j], zeta1[j], zeta2[j], asum[j], bsum[j]) = zunhj(zn, modified_order, false);
        let s1 = -scaling.scale_zetas(zb, modified_order, zeta1[j], zeta2[j]);
        let of = Overflow::find_overflow(s1.re, phi[j], -0.25 * arg[j].abs().ln() - AIC);

        let mut handle_underflow = |of_already: &mut bool, cs_: &mut Complex64| {
            //-----------------------------------------------------------------------
            //     FOR ZR < 0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
            //-----------------------------------------------------------------------
            if z.re < 0.0 {
                return Err(Overflow);
            }
            *of_already = false;
            y[i] = c_zero();
            nz += 1;
            *cs_ *= -Complex64::I;
            if i != 0 && y[i - 1] != c_zero() {
                y[i - 1] = c_zero();
                nz += 1;
            }
            Ok(())
        };

        if !underflowed_already {
            overflow_state_k = of;
        }

        match of {
            Overflow::Over(_) => return Err(Overflow),

            Overflow::Under(_) => handle_underflow(&mut underflowed_already, &mut cs)?,
            Overflow::NearOver | Overflow::NearUnder | Overflow::None => {
                //-----------------------------------------------------------------------;
                //     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR;
                //     EXPONENT EXTREMES;
                //-----------------------------------------------------------------------;
                let c2 = CR2 * arg[j];

                let (airy, d_airy) = airy_pair(c2);
                let pt = ((d_airy * bsum[j].unwrap()) * CR2 + (airy * asum[j].unwrap())) * phi[j];
                let mut s2 = pt * cs;
                let s1 = s1.exp() * MACHINE_CONSTANTS.scaling_factors[overflow_state_k];
                s2 *= s1;
                if overflow_state_k == Overflow::NearUnder
                    && will_underflow(
                        s2,
                        MACHINE_CONSTANTS.smallness_threshold[0],
                        MACHINE_CONSTANTS.abs_error_tolerance,
                    )
                {
                    handle_underflow(&mut underflowed_already, &mut cs)?
                }
                if zr.im <= 0.0 {
                    s2 = s2.conj();
                }
                cy[underflowed_already as usize] = s2;
                y[i] = s2 * MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state_k];
                cs = -Complex64::I * cs;
                if underflowed_already {
                    break;
                }
                underflowed_already = true;
            }
        };
    }

    let rz = 2.0 * zr.conj() / zr.abs().powi(2);
    let mut ck = modified_order * rz;
    let mut phid = c_zero();
    let mut argd = c_zero();
    let mut zeta1d = c_zero();
    let mut zeta2d = c_zero();
    let mut asumd = None;
    let mut bsumd = None;
    let do_overflow_check = n_elements_set < n;
    if do_overflow_check {
        //-----------------------------------------------------------------------;
        //     TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW. SET SEQUENCE TO ZERO;
        //     ON UNDERFLOW.;
        //-----------------------------------------------------------------------;
        modified_order = order + ((n - 1) as f64);
        (phid, argd, zeta1d, zeta2d, asumd, bsumd) =
            zunhj(zn, modified_order, rotation == RotationDirection::None);
        let s1 = -scaling.scale_zetas(zb, modified_order, zeta1d, zeta2d);
        match Overflow::find_overflow(s1.re, phid, 0.0) {
            Overflow::Over(_) => return Err(Overflow),

            Overflow::Under(_) => {
                if z.re < 0.0 {
                    return Err(Overflow);
                }
                return Ok((c_zeros(n), nz));
            }
            Overflow::NearOver | Overflow::None | Overflow::NearUnder => (),
        }
        let [mut s1, mut s2] = cy;
        let mut recip_scaling = MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state_k];
        let mut ascle = MACHINE_CONSTANTS.smallness_threshold[overflow_state_k];

        for yi in y.iter_mut().skip(n_elements_set) {
            (s1, s2) = (s2, s2 * ck + s1);
            ck += rz;
            c2 = s2 * recip_scaling;
            *yi = c2;
            if overflow_state_k == Overflow::NearOver {
                continue;
            }
            if max_abs_component(c2) <= ascle {
                continue;
            }
            overflow_state_k.increment();
            ascle = MACHINE_CONSTANTS.smallness_threshold[overflow_state_k];
            s1 *= recip_scaling;
            s2 = c2;
            s1 *= MACHINE_CONSTANTS.scaling_factors[overflow_state_k];
            s2 *= MACHINE_CONSTANTS.scaling_factors[overflow_state_k];
            recip_scaling = MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state_k];
        }
    }
    if rotation == RotationDirection::None {
        return Ok((y, nz));
    }
    //-----------------------------------------------------------------------
    //     ANALYTIC CONTINUATION FOR RE(Z) < 0.0
    //-----------------------------------------------------------------------
    nz = 0;
    let sgn = -PI * rotation.signum();
    //-----------------------------------------------------------------------
    //     CSPN AND CSGN ARE COEFF OF K AND I FUNCIONS RESP.
    //-----------------------------------------------------------------------
    let csgn = if zr.im <= 0.0 { -sgn } else { sgn };
    let modified_integer_order = integer_order + n - 1;
    let mut cspn = Complex64::cis(order_fract * sgn);
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
    let cos_sin = Complex64::cis(ANG);
    // let mut cs = Complex64::I * Complex64::from_polar(CSGNI, ANG);
    let mut cs = csgn * Complex64::new(cos_sin.im, cos_sin.re);
    cs *= CIP[modified_integer_order % 4];
    let mut IUF = 0;

    underflowed_already = false;
    let mut overflow_state_i = Overflow::None;
    let mut remaining_n = n;
    for (kk, yi) in y.iter_mut().enumerate().rev() {
        remaining_n = kk;
        modified_order = order + (kk as f64);
        //-----------------------------------------------------------------------
        //     LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K
        //     FUNCTION ABOVE
        //-----------------------------------------------------------------------
        // Note that, is the overflow check was done, the ___d are already set, and
        // valid for kk == n-1. Also that kk == n-1 on the first pas through this loop.
        let use_preset_overflow = (kk == n - 1) && do_overflow_check;
        // these where the last two kk values where phi etc where recorded in the previous run.
        // Would it be better to store all of them?!
        let in_last_two_set = (kk == n_elements_set - 1) || (kk == n_elements_set - 2);
        if n <= 2 || (!use_preset_overflow) && in_last_two_set {
            phid = phi[j];
            argd = arg[j];
            zeta1d = zeta1[j];
            zeta2d = zeta2[j];
            asumd = asum[j];
            bsumd = bsum[j];
            j = 1 - j;
        } else if !(use_preset_overflow || in_last_two_set) {
            (phid, argd, zeta1d, zeta2d, asumd, bsumd) = zunhj(zn, modified_order, false);
        } else {
            // Case were overflow check has already set the ___d variables ?
        }
        let mut s1 = scaling.scale_zetas(zb, modified_order, zeta1d, zeta2d);

        let of = Overflow::find_overflow(s1.re, phid, -0.25 * argd.abs().ln() - AIC);
        if !underflowed_already {
            overflow_state_i = if matches!(of, Overflow::Under(_)) {
                Overflow::None
            } else {
                of
            };
        }
        let mut s2 = match of {
            Overflow::Over(_) => return Err(Overflow),
            Overflow::Under(_) => c_zero(),
            Overflow::NearOver | Overflow::None | Overflow::NearUnder => {
                let (airy, d_airy) = airy_pair(argd);
                let pt = ((d_airy * bsumd.unwrap()) + (airy * asumd.unwrap())) * phid;
                let mut s2 = pt * cs;
                s1 = s1.exp() * MACHINE_CONSTANTS.scaling_factors[overflow_state_i];
                s2 *= s1;
                if overflow_state_i == Overflow::NearUnder
                    && will_underflow(
                        s2,
                        MACHINE_CONSTANTS.smallness_threshold[0],
                        MACHINE_CONSTANTS.abs_error_tolerance,
                    )
                {
                    s2 = c_zero();
                }
                s2
            }
        };
        if zr.im <= 0.0 {
            s2 = s2.conj();
        }
        cy[underflowed_already as usize] = s2;
        let c2 = s2;
        s2 *= MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state_i];
        //-----------------------------------------------------------------------;
        //     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N;
        //-----------------------------------------------------------------------;
        s1 = *yi;
        if scaling == Scaling::Scaled {
            nz += underflow_add_i_k(zr, &mut s1, &mut s2, &mut IUF);
        }
        *yi = s1 * cspn + s2;
        cspn = -cspn;
        cs *= -Complex64::I;
        if c2 == c_zero() {
            underflowed_already = false;
        } else {
            if underflowed_already {
                break;
            }
            underflowed_already = true;
        }
    }

    if remaining_n == 0 {
        return Ok((y, nz));
    }
    //-----------------------------------------------------------------------
    //     RECUR BACKWARD FOR REMAINDER OF I SEQUENCE AND ADD IN THE
    //     K FUNCTIONS, SCALING THE I SEQUENCE DURING RECURRENCE TO KEEP
    //     INTERMEDIATE ARITHMETIC ON SCALE NEAR EXPONENT EXTREMES.
    //-----------------------------------------------------------------------
    let [mut s1, mut s2] = cy;

    let mut recip_scale_factor = MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state_i];
    let mut ascle = MACHINE_CONSTANTS.smallness_threshold[overflow_state_i];
    modified_order = (integer_order + remaining_n) as f64;
    for yi in y.iter_mut().take(remaining_n).rev() {
        (s1, s2) = (s2, s1 + (modified_order + order_fract) * (rz * s2));
        modified_order -= 1.0;
        let mut c2 = s2 * recip_scale_factor;
        let old_c2 = c2;
        let mut c1 = *yi;
        if scaling == Scaling::Scaled {
            nz += underflow_add_i_k(zr, &mut c1, &mut c2, &mut IUF);
        }
        *yi = c1 * cspn + c2;
        cspn = -cspn;
        if overflow_state_i != Overflow::NearOver && max_abs_component(c2) > ascle {
            overflow_state_i.increment();
            ascle = MACHINE_CONSTANTS.smallness_threshold[overflow_state_i];
            s1 *= recip_scale_factor;
            s2 = old_c2;
            s1 *= MACHINE_CONSTANTS.scaling_factors[overflow_state_i];
            s2 *= MACHINE_CONSTANTS.scaling_factors[overflow_state_i];
            recip_scale_factor = MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state_i];
        }
    }
    Ok((y, nz))
}

fn ZBUNI(
    z: Complex64,
    order: f64,
    KODE: Scaling,
    N: usize,
    NUI: usize,
    y: &mut [Complex64],
) -> Result<(usize, usize), BesselError> {
    // ***BEGIN PROLOGUE  ZBUNI
    // ***REFER TO  ZBESI,ZBESK
    //
    //     ZBUNI COMPUTES THE I BESSEL FUNCTION FOR LARGE CABS(Z) >
    //     FNUL AND FNU+N-1 < FNUL. THE ORDER IS INCREASED FROM
    //     FNU+N-1 GREATER THAN FNUL BY ADDING NUI AND COMPUTING
    //     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR I(FNU,Z)
    //     ON IFORM=1 AND THE EXPANSION FOR J(FNU,Z) ON IFORM=2
    //
    // ***ROUTINES CALLED  ZUNI1,ZUNI2,ZABS,d1mach
    // ***END PROLOGUE  ZBUNI

    let imaginary_dominant = imaginary_dominant(z);
    if NUI != 0 {
        let mut FNUI = NUI as f64;
        let DFNU = order + ((N - 1) as f64);
        let GNU = DFNU + FNUI;
        let mut cy = [c_zero(); 2];
        let (NW, NLAST) = if imaginary_dominant {
            //-----------------------------------------------------------------------
            //     ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*FRAC_PI_2)) FOR LARGE FNU
            //     APPLIED IN PI/3 < ABS(ARG(Z)) <= PI/2 WHERE M=+I OR -I
            //     AND FRAC_PI_2=PI/2
            //-----------------------------------------------------------------------
            ZUNI2(z, GNU, KODE, 2, &mut cy)?
        } else {
            //-----------------------------------------------------------------------
            //     ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN
            //     -PI/3 <= ARG(Z) <= PI/3
            //-----------------------------------------------------------------------
            ZUNI1(z, GNU, KODE, 2, &mut cy)?
        };
        if NW != 0 {
            return Ok((0, N));
        }
        //----------------------------------------------------------------------
        //     SCALE BACKWARD RECURRENCE, BRY(3) IS DEFINED BUT NEVER USED
        //----------------------------------------------------------------------
        let (mut overflow_state, mut ASCLE, mut CSCLR) =
            if cy[0].abs() <= MACHINE_CONSTANTS.smallness_threshold[0] {
                (
                    Overflow::NearUnder,
                    MACHINE_CONSTANTS.smallness_threshold[0],
                    1.0 / MACHINE_CONSTANTS.abs_error_tolerance,
                )
            } else if cy[0].abs() >= MACHINE_CONSTANTS.smallness_threshold[1] {
                (
                    Overflow::NearOver,
                    MACHINE_CONSTANTS.smallness_threshold[2],
                    MACHINE_CONSTANTS.abs_error_tolerance,
                )
            } else {
                (
                    Overflow::None,
                    MACHINE_CONSTANTS.smallness_threshold[1],
                    1.0,
                )
            };

        let mut CSCRR = 1.0 / CSCLR;
        let mut s1 = cy[1] * CSCLR;
        let mut s2 = cy[0] * CSCLR;
        // working out rz in multiple steps seems to give different floating point answer.
        let rz = 2.0 * z.conj() / z.abs().pow(2);

        for _ in 0..NUI {
            let st = s2;
            s2 = (DFNU + FNUI) * rz * s2 + s1;
            s1 = st;
            FNUI -= 1.0;
            if overflow_state == Overflow::NearOver {
                continue;
            }
            let st = s2 * CSCRR;
            if max_abs_component(st) <= ASCLE {
                continue;
            }
            overflow_state.increment();
            ASCLE = MACHINE_CONSTANTS.smallness_threshold[overflow_state];
            s1 *= CSCRR;
            s2 = st;
            CSCLR *= MACHINE_CONSTANTS.abs_error_tolerance;
            CSCRR = 1.0 / CSCLR;
            s1 *= CSCLR;
            s2 *= CSCLR;
        }
        y[N - 1] = s2 * CSCRR;
        if N == 1 {
            return Ok((0, NLAST));
        }
        let NL = N - 1;
        FNUI = NL as f64;
        let mut K = NL;
        for _ in 0..NL {
            let st = s2;
            s2 = (order + FNUI) * (rz * s2) + s1;
            s1 = st;
            y[K - 1] = s2 * CSCRR;
            FNUI -= 1.0;
            K -= 1;
            if overflow_state == Overflow::NearOver {
                continue;
            }
            // using K (rather than K-1) below as Amos "saved" the y value before K was decremented
            if max_abs_component(y[K]) <= ASCLE {
                continue;
            }
            overflow_state.increment();
            ASCLE = MACHINE_CONSTANTS.smallness_threshold[overflow_state];
            s1 *= CSCRR;
            s2 = y[K - 1];
            CSCLR *= MACHINE_CONSTANTS.abs_error_tolerance;
            CSCRR = 1.0 / CSCLR;
            s1 *= CSCLR;
            s2 *= CSCLR;
        }
        return Ok((0, NLAST));
    }
    let (NW, NLAST) = if imaginary_dominant {
        //-----------------------------------------------------------------------
        //     ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*FRAC_PI_2)) FOR LARGE FNU
        //     APPLIED IN PI/3 < ABS(ARG(Z)) <= PI/2 WHERE M=+I OR -I
        //     AND FRAC_PI_2=PI/2
        //-----------------------------------------------------------------------
        ZUNI2(z, order, KODE, N, y)?
    } else {
        //-----------------------------------------------------------------------
        //     ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN
        //     -PI/3 <= ARG(Z) <= PI/3
        //-----------------------------------------------------------------------
        ZUNI1(z, order, KODE, N, y)?
    };

    Ok((NW, NLAST))
}

fn ZUNI1(
    z: Complex64,
    order: f64,
    scaling: Scaling,
    N: usize,
    y: &mut [Complex64],
) -> BesselResult<(usize, usize)> {
    // ***BEGIN PROLOGUE  ZUNI1
    // ***REFER TO  ZBESI,ZBESK
    //
    //     ZUNI1 COMPUTES I(FNU,Z)  BY MEANS OF THE UNIFORM ASYMPTOTIC
    //     EXPANSION FOR I(FNU,Z) IN -PI/3 <= ARG Z <= PI/3.
    //
    //     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
    //     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
    //     NLAST != 0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
    //     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1 < FNUL.
    //     Y(I)=CZERO FOR I=NLAST+1,N

    let mut NZ = 0;
    let mut ND = N;
    let NLAST = 0;
    //-----------------------------------------------------------------------
    //     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
    //-----------------------------------------------------------------------
    let mut modified_order = order.max(1.0);
    let (_, zeta1, zeta2, _) = zunik(z, modified_order, IKType::I, true);
    let s1 = scaling.scale_zetas(z, modified_order, zeta1, zeta2);
    // phi is chosen here for refined tests to equal the original tests
    // which don't test refinement
    match Overflow::find_overflow(s1.re, c_one(), 0.0) {
        Overflow::Over(_) => return Err(Overflow),
        Overflow::Under(_) => return Ok((N, NLAST)),
        _ => (),
    }
    let mut overflow_state = Overflow::None; // this value should never be used
    let mut cy = [c_zero(); 2];
    let mut set_underflow_and_update = false;
    // Unroll loop here and ZUNI2. Check ZUNK1, 2 as well.
    'l30: loop {
        if set_underflow_and_update {
            y[ND - 1] = c_zero();
            NZ += 1;
            ND -= 1;
            if ND == 0 {
                return Ok((NZ, NLAST));
            }
            let NUF = check_underflow_uniform_asymp_params(z, order, scaling, IKType::I, ND, y)?;
            ND -= NUF;
            NZ += NUF;
            if ND == 0 {
                return Ok((NZ, NLAST));
            }
            modified_order = order + ((ND - 1) as f64);
            if modified_order < MACHINE_CONSTANTS.asymptotic_order_limit {
                return Ok((NZ, ND));
            }
        }

        for i in 0..2.min(ND) {
            modified_order = order + ((ND - (i + 1)) as f64);
            let (phi, zeta1, zeta2, sum) = zunik(z, modified_order, IKType::I, false);
            let sum = sum.unwrap();
            let mut s1 = scaling.scale_zetas(z, modified_order, zeta1, zeta2);
            if scaling == Scaling::Scaled {
                s1 += Complex64::new(0.0, z.im);
            }

            let of = Overflow::find_overflow(s1.re, phi, 0.0);
            if i == 0 {
                overflow_state = of;
            }
            match of {
                Overflow::Over(_) => return Err(Overflow),
                Overflow::Under(_) => {
                    set_underflow_and_update = true;
                    continue 'l30;
                }
                _ => (),
            }
            //-----------------------------------------------------------------------
            //     SCALE S1 if CABS(S1) < ASCLE
            //-----------------------------------------------------------------------
            let mut s2 = phi * sum;
            s1 = MACHINE_CONSTANTS.scaling_factors[overflow_state] * s1.exp();
            s2 *= s1;
            if overflow_state == Overflow::NearUnder
                && will_underflow(
                    s2,
                    MACHINE_CONSTANTS.smallness_threshold[0],
                    MACHINE_CONSTANTS.abs_error_tolerance,
                )
            {
                set_underflow_and_update = true;
                continue 'l30;
            }
            cy[i] = s2;
            y[ND - i - 1] = s2 * MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
        }
        break 'l30;
    }
    if ND <= 2 {
        return Ok((NZ, NLAST));
    }
    let rz = 2.0 * z.conj() / z.abs().pow(2);
    let [mut s1, mut s2] = cy;
    let mut C1R = MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
    let mut ASCLE = MACHINE_CONSTANTS.smallness_threshold[overflow_state];
    let mut K = ND - 2;
    modified_order = K as f64;
    for _ in 2..ND {
        let mut c2 = s2;
        s2 = s1 + (order + modified_order) * (rz * c2);
        s1 = c2;
        c2 = s2 * C1R;
        y[K - 1] = c2;
        K -= 1;
        modified_order -= 1.0;
        if overflow_state == Overflow::NearOver {
            continue;
        }
        if max_abs_component(c2) <= ASCLE {
            continue;
        }
        overflow_state.increment();
        ASCLE = MACHINE_CONSTANTS.smallness_threshold[overflow_state];
        s1 *= C1R;
        s2 = c2;
        s1 *= MACHINE_CONSTANTS.scaling_factors[overflow_state];
        s2 *= MACHINE_CONSTANTS.scaling_factors[overflow_state];
        C1R = MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
    }
    Ok((NZ, NLAST))
}

fn ZUNI2(
    z: Complex64,
    order: f64,
    scaling: Scaling,
    N: usize,
    y: &mut [Complex64],
) -> BesselResult<(usize, usize)> {
    // ***BEGIN PROLOGUE  ZUNI2
    // ***REFER TO  ZBESI,ZBESK
    //
    //     ZUNI2 COMPUTES I(FNU,Z) IN THE RIGHT HALF PLANE BY MEANS OF
    //     UNIFORM ASYMPTOTIC EXPANSION FOR J(FNU,ZN) WHERE ZN IS Z*I
    //     OR -Z*I AND ZN IS IN THE RIGHT HALF PLANE ALSO.
    //
    //     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
    //     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
    //     NLAST != 0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
    //     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1 < FNUL.
    //     Y(I)=CZERO FOR I=NLAST+1,N
    //
    // ***ROUTINES CALLED  ZAIRY,ZUunderflowCHK,ZUNHJ,ZUOIK,d1mach,ZABS
    // ***END PROLOGUE  ZUNI2

    let mut NZ = 0;
    let mut ND = N;
    let NLAST = 0;
    //-----------------------------------------------------------------------
    //     ZN IS IN THE RIGHT HALF PLANE AFTER ROTATION BY CI OR -CI
    //-----------------------------------------------------------------------
    let mut zn = Complex64::new(z.im, -z.re);
    let mut zb = z;
    let mut CIDI = -1.0;
    let INU = order as usize;
    let ANG = FRAC_PI_2 * order.fract();
    let mut c2 = Complex64::cis(ANG);
    let CAR = c2.re;
    let SAR = c2.im;
    let index = (INU + N - 1) % 4;
    c2 *= CIP[index];
    if z.im <= 0.0 {
        zn.re = -zn.re;
        zb.im = -zb.im;
        CIDI = -CIDI;
        c2.im = -c2.im;
    }
    //-----------------------------------------------------------------------
    //     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
    //-----------------------------------------------------------------------
    let mut modified_order = order.max(1.0);
    let (_, _, zeta1, zeta2, _, _) = zunhj(zn, modified_order, true);

    let s1 = scaling.scale_zetas(zb, modified_order, zeta1, zeta2);

    // phi is chosen here for refined tests to equal the original tests
    // which don't test refinement
    match Overflow::find_overflow(s1.re, c_one(), 0.0) {
        Overflow::Over(_) => return Err(Overflow),
        Overflow::Under(_) => return Ok((N, NLAST)),
        _ => (),
    }
    let mut set_underflow_and_update = false;
    'l40: loop {
        if set_underflow_and_update {
            //-----------------------------------------------------------------------
            //     SET UNDERFLOW AND UPDATE PARAMETERS
            //-----------------------------------------------------------------------
            y[ND - 1] = c_zero();
            NZ += 1;
            ND -= 1;
            if ND == 0 {
                return Ok((NZ, NLAST));
            }
            let NUF = check_underflow_uniform_asymp_params(z, order, scaling, IKType::I, ND, y)?;
            ND -= NUF;
            NZ += NUF;
            if ND == 0 {
                return Ok((NZ, NLAST));
            }
            modified_order = order + ((ND - 1) as f64);
            if modified_order < MACHINE_CONSTANTS.asymptotic_order_limit {
                return Ok((NZ, ND));
            }
            let index = (INU + ND - 1) % 4;
            c2 = Complex64::new(CAR, SAR) * CIP[index];
            if z.im <= 0.0 {
                c2 = c2.conj();
            }
        }
        let mut overflow_state = Overflow::NearUnder;
        let mut cy = [c_zero(); 2];
        for i in 0..ND.min(2) {
            modified_order = order + ((ND - (i + 1)) as f64);
            let (phi, arg, zeta1, zeta2, asum, bsum) = zunhj(zn, modified_order, false);
            let asum = asum.unwrap();
            let bsum = bsum.unwrap();
            let mut s1 = scaling.scale_zetas(zb, modified_order, zeta1, zeta2);
            if scaling == Scaling::Scaled {
                s1 += Complex64::I * z.im.abs();
            }

            //-----------------------------------------------------------------------
            //     TEST FOR UNDERFLOW AND OVERFLOW
            //-----------------------------------------------------------------------
            let of = Overflow::find_overflow(s1.re, phi, -0.25 * arg.abs().ln() - AIC);
            if i == 0 {
                overflow_state = of;
            }
            match of {
                Overflow::Over(_) => return Err(Overflow),
                Overflow::Under(_) => {
                    set_underflow_and_update = true;
                    continue 'l40;
                }
                _ => (),
            }
            //-----------------------------------------------------------------------
            //     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
            //     EXPONENT EXTREMES
            //-----------------------------------------------------------------------
            let (a_airy, d_airy) = airy_pair(arg);

            let mut s2 = phi * (d_airy * bsum + a_airy * asum);
            let s1 = MACHINE_CONSTANTS.scaling_factors[overflow_state] * s1.exp();
            s2 *= s1;
            if overflow_state == Overflow::NearUnder
                && will_underflow(
                    s2,
                    MACHINE_CONSTANTS.smallness_threshold[0],
                    MACHINE_CONSTANTS.abs_error_tolerance,
                )
            {
                set_underflow_and_update = true;
                continue 'l40;
            }
            if z.im <= 0.0 {
                s2 = s2.conj();
            }
            s2 *= c2;
            cy[i] = s2;
            y[ND - i - 1] = s2 * MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
            c2 *= CIDI * Complex64::I;
        }
        if ND <= 2 {
            break 'l40;
        }
        let rz = 2.0 * z.conj() / z.abs().pow(2);
        let [mut s1, mut s2] = cy;
        let mut C1R = MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
        let mut ASCLE = MACHINE_CONSTANTS.smallness_threshold[overflow_state];
        for K in (0..(ND - 2)).rev() {
            let st = s2;
            s2 = s1 + (order + ((K + 1) as f64)) * rz * s2;
            s1 = st;
            y[K] = s2 * C1R;
            if overflow_state == Overflow::NearOver {
                continue;
            }
            if max_abs_component(y[K]) <= ASCLE {
                continue;
            }
            overflow_state.increment();
            ASCLE = MACHINE_CONSTANTS.smallness_threshold[overflow_state];
            s1 *= C1R;
            s2 = y[K];
            s1 *= MACHINE_CONSTANTS.scaling_factors[overflow_state];
            s2 *= MACHINE_CONSTANTS.scaling_factors[overflow_state];
            C1R = MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
        }
        break 'l40;
    }
    Ok((NZ, NLAST))
}
