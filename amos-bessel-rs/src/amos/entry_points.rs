use std::f64::consts::{FRAC_2_PI, FRAC_PI_2, PI};

use num::complex::{Complex64, ComplexFloat};

use crate::{
    Scaling,
    amos::{
        BesselError::*,
        BesselResult, CIP, HankelKind, IKType, MACHINE_CONSTANTS, RotationDirection, c_one, c_zero,
        c_zeros, max_abs_component,
        overflow_checks::check_underflow_uniform_asymp_params,
        translator::{
            ZACAI, ZBUNK, airy_power_series, analytic_continuation, i_right_half_plane,
            k_right_half_plane,
        },
        utils::{TWO_THIRDS, is_sigificance_lost, sanitise_inputs},
    },
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
    let (mut cy, nz) = if order < MACHINE_CONSTANTS.asymptotic_order_limit {
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
                        Err(PartialLossOfSignificance { y: cy, nz })
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
        let (cy, nw) = ZBUNK(zn, order, scaling, asymptotic_rotation, n)?;
        nz += nw;
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
        let scaling =
            if max_abs_component(*element) < MACHINE_CONSTANTS.absolute_approximation_limit {
                *element *= MACHINE_CONSTANTS.rtol;
                MACHINE_CONSTANTS.abs_error_tolerance
            } else {
                1.0
            };
        *element *= csgn * scaling;
        csgn *= Complex64::I * -rotation_f64;
    }
    if partial_loss_of_significance {
        Err(PartialLossOfSignificance { y: cy, nz })
    } else {
        Ok((cy, nz))
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
        let mut scaling = 1.0;
        // TODO is the below a pattern?
        if (max_abs_component(*cyi)) <= MACHINE_CONSTANTS.absolute_approximation_limit {
            *cyi *= MACHINE_CONSTANTS.rtol;
            scaling = MACHINE_CONSTANTS.abs_error_tolerance;
        }
        *cyi *= csgn * scaling;
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
        let n_underflow =
            check_underflow_uniform_asymp_params(z, order, scaling, IKType::K, n, &mut y)?;
        nz += n_underflow;

        //-----------------------------------------------------------------------;
        //     HERE NN=N OR NN=0 SINCE NUF=0,NN, OR -1 ON RETURN FROM CUOIK;
        //     if NUF=NN, THEN CY(I)=CZERO FOR ALL I;
        //-----------------------------------------------------------------------;
        if n_underflow == n {
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
    let mut csgn = Complex64::cis(FRAC_PI_2 * frac_order);
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
        let (cy, nz) = if re_zta < 0.0 || z.re <= 0.0 {
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
        (s1 / scale_factor, nz)
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
