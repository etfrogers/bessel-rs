#![allow(non_snake_case, clippy::excessive_precision)]
use super::{
    BesselError, BesselResult, HankelKind, IKType, Scaling, c_one, c_zero, c_zeros, gamma_ln,
    i_power_series,
    overflow_checks::{zunik, zuoik},
    utils::{TWO_THIRDS, is_sigificance_lost, will_z_underflow},
};
use crate::amos::{
    BesselError::*,
    MACHINE_CONSTANTS, max_abs_component,
    overflow_checks::{Overflow, zunhj},
    utils::imaginary_dominant,
    z_asymptotic_i::z_asymptotic_i,
};
use num::{
    Zero,
    complex::{Complex64, ComplexFloat},
    pow::Pow,
};
use std::{
    cmp::min,
    f64::consts::{FRAC_PI_2, PI},
};

pub fn zbesh(z: Complex64, order: f64, KODE: Scaling, M: HankelKind, N: usize) -> BesselResult {
    // ***BEGIN PROLOGUE  ZBESH
    // ***DATE WRITTEN   830501   (YYMMDD)
    // ***REVISION DATE  890801, 930101   (YYMMDD)
    // ***CATEGORY NO.  B5K
    // ***KEYWORDS  H-BESSEL FUNCTIONS,BESSEL FUNCTIONS OF COMPLEX ARGUMENT,
    //             BESSEL FUNCTIONS OF THIRD KIND,HANKEL FUNCTIONS
    // ***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
    // ***PURPOSE  TO COMPUTE THE H-BESSEL FUNCTIONS OF A COMPLEX ARGUMENT
    // ***DESCRIPTION
    //
    //                      ***A DOUBLE PRECISION ROUTINE***
    //         ON KODE=1, ZBESH COMPUTES AN N MEMBER SEQUENCE OF COMPLEX
    //         HANKEL (BESSEL) FUNCTIONS CY(J)=H(M,FNU+J-1,Z) FOR KINDS M=1
    //         OR 2, REAL, NONNEGATIVE ORDERS FNU+J-1, J=1,...,N, AND COMPLEX
    //         Z != CMPLX(0.0,0.0) IN THE CUT PLANE -PI < ARG(Z) <= PI.
    //         ON KODE=2, ZBESH RETURNS THE SCALED HANKEL FUNCTIONS
    //
    //         CY(I)=EXP(-MM*Z*I)*H(M,FNU+J-1,Z)       MM=3-2*M,   I**2=-1.
    //
    //         WHICH REMOVES THE EXPONENTIAL BEHAVIOR IN BOTH THE UPPER AND
    //         LOWER HALF PLANES. DEFINITIONS AND NOTATION ARE FOUND IN THE
    //         NBS HANDBOOK OF MATHEMATICAL FUNCTIONS (REF. 1).
    //
    //         INPUT      ZR,ZI,FNU ARE DOUBLE PRECISION
    //           ZR,ZI  - Z=CMPLX(ZR,ZI), Z != CMPLX(0.0,0.0),
    //                    -PT < ARG(Z) <= PI
    //           FNU    - ORDER OF INITIAL H FUNCTION, FNU >= 0.0
    //           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
    //                    KODE= 1  RETURNS
    //                             CY(J)=H(M,FNU+J-1,Z),   J=1,...,N
    //                        = 2  RETURNS
    //                             CY(J)=H(M,FNU+J-1,Z)*EXP(-I*Z*(3-2M))
    //                                  J=1,...,N  ,  I**2=-1
    //           M      - KIND OF HANKEL FUNCTION, M=1 OR 2
    //           N      - NUMBER OF MEMBERS IN THE SEQUENCE, N >= 1
    //
    //         OUTPUT     CYR,CYI ARE DOUBLE PRECISION
    //           CYR,CYI- DOUBLE PRECISION VECTORS WHOSE FIRST N COMPONENTS
    //                    CONTAIN REAL AND IMAGINARY PARTS FOR THE SEQUENCE
    //                    CY(J)=H(M,FNU+J-1,Z)  OR
    //                    CY(J)=H(M,FNU+J-1,Z)*EXP(-I*Z*(3-2M))  J=1,...,N
    //                    DEPENDING ON KODE, I**2=-1.
    //           NZ     - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW,
    //                    NZ= 0   , NORMAL RETURN
    //                    NZ > 0 , FIRST NZ COMPONENTS OF CY SET TO ZERO DUE
    //                              TO UNDERFLOW, CY(J)=CMPLX(0.0,0.0)
    //                              J=1,...,NZ WHEN Y > 0.0 AND M=1 OR
    //                              Y < 0.0 AND M=2. FOR THE COMPLMENTARY
    //                              HALF PLANES, NZ STATES ONLY THE NUMBER
    //                              OF UNDERFLOWS.
    //           IERR   - ERROR FLAG
    //                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
    //                    IERR=1, INPUT ERROR   - NO COMPUTATION
    //                    IERR=2, OVERFLOW      - NO COMPUTATION, FNU TOO
    //                            LARGE OR CABS(Z) TOO SMALL OR BOTH
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
    //         THE COMPUTATION IS CARRIED OUT BY THE RELATION
    //
    //         H(M,FNU,Z)=(1/MP)*EXP(-MP*FNU)*K(FNU,Z*EXP(-MP))
    //             MP=MM*FRAC_PI_2*I,  MM=3-2*M,  FRAC_PI_2=PI/2,  I**2=-1
    //
    //         FOR M=1 OR 2 WHERE THE K BESSEL FUNCTION IS COMPUTED FOR THE
    //         RIGHT HALF PLANE RE(Z) >= 0.0. THE K FUNCTION IS CONTINUED
    //         TO THE LEFT HALF PLANE BY THE RELATION
    //
    //         K(FNU,Z*EXP(MP)) = EXP(-MP*FNU)*K(FNU,Z)-MP*I(FNU,Z)
    //         MP=MR*PI*I, MR=+1 OR -1, RE(Z) > 0, I**2=-1
    //
    //         WHERE I(FNU,Z) IS THE I BESSEL FUNCTION.
    //
    //         EXPONENTIAL DECAY OF H(M,FNU,Z) OCCURS IN THE UPPER HALF Z
    //         PLANE FOR M=1 AND THE LOWER HALF Z PLANE FOR M=2.  EXPONENTIAL
    //         GROWTH OCCURS IN THE COMPLEMENTARY HALF PLANES.  SCALING
    //         BY EXP(-MM*Z*I) REMOVES THE EXPONENTIAL BEHAVIOR IN THE
    //         WHOLE Z PLANE FOR Z TO INFINITY.
    //
    //         FOR NEGATIVE ORDERS,THE FORMULAE
    //
    //               H(1,-FNU,Z) = H(1,FNU,Z)*CEXP( PI*FNU*I)
    //               H(2,-FNU,Z) = H(2,FNU,Z)*CEXP(-PI*FNU*I)
    //                         I**2=-1
    //
    //         CAN BE USED.
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
    //         ROUNDOFF,1.0e-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
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
    // ***ROUTINES CALLED  ZACON,ZBKNU,ZBUNK,ZUOIK,ZABS,i1mach,d1mach
    // ***END PROLOGUE  ZBESH
    let mut err = None;
    if z.re == 0.0 && z.im == 0.0 {
        err = Some("z must not be zero");
    }
    if order < 0.0_f64 {
        err = Some("order must be positive");
    };
    if N < 1 {
        err = Some("N must be >= 1");
    };
    if let Some(details) = err {
        return Err(BesselError::InvalidInput {
            details: details.to_owned(),
        });
    }
    let mut NZ = 0;

    let mut NN = N;
    let FN = order + ((NN - 1) as f64);
    let MM: i64 = 3_i64 - 2 * <HankelKind as Into<i64>>::into(M);
    let FMM = MM as f64;
    let mut zn = -Complex64::I * FMM * z;
    //-----------------------------------------------------------------------
    //     TEST FOR PROPER RANGE
    //-----------------------------------------------------------------------
    let abs_z = z.abs();
    let partial_loss_of_significance = is_sigificance_lost(abs_z, FN, false)?;
    //-----------------------------------------------------------------------
    //     OVERFLOW TEST ON THE LAST MEMBER OF THE SEQUENCE
    //-----------------------------------------------------------------------
    if abs_z < MACHINE_CONSTANTS.underflow_limit {
        return Err(Overflow);
    }
    let (mut cy, NZ) = if order < MACHINE_CONSTANTS.asymptotic_order_limit {
        if FN > 1.0 {
            if FN > 2.0 {
                let mut cy = c_zeros(N);
                let NUF = zuoik(zn, order, KODE, IKType::K, NN, &mut cy)?;

                NZ += NUF;
                NN -= NUF;
                //-----------------------------------------------------------------------
                //     HERE NN=N OR NN=0 SINCE NUF=0,NN, OR -1 ON RETURN FROM CUOIK
                //     if NUF=NN, THEN CY(I)=CZERO FOR ALL I
                //-----------------------------------------------------------------------
                if NN == 0 {
                    return if zn.re < 0.0 {
                        Err(Overflow)
                    } else if partial_loss_of_significance {
                        Err(BesselError::PartialLossOfSignificance { y: cy, nz: NZ })
                    } else {
                        Ok((cy, NZ))
                    };
                }
            }
            if abs_z <= MACHINE_CONSTANTS.abs_error_tolerance
                && -FN * (0.5 * abs_z).ln() > MACHINE_CONSTANTS.exponent_limit
            {
                return Err(Overflow);
            }
        }
        if !((zn.re < 0.0) || (zn.re == 0.0 && zn.im < 0.0 && M == HankelKind::Second)) {
            //-----------------------------------------------------------------------
            //     RIGHT HALF PLANE COMPUTATION, XN >= 0. && (XN != 0. ||
            //     YN >= 0. || M=1)
            //-----------------------------------------------------------------------
            ZBKNU(zn, order, KODE, NN)?
        } else {
            //-----------------------------------------------------------------------
            //     LEFT HALF PLANE COMPUTATION
            //-----------------------------------------------------------------------
            analytic_continuation(zn, order, KODE, -MM, NN)?
        }
    } else {
        //-----------------------------------------------------------------------
        //     UNIFORM ASYMPTOTIC EXPANSIONS FOR FNU > FNUL
        //-----------------------------------------------------------------------
        let mut MR = 0;
        if !((zn.re >= 0.0) && (zn.re != 0.0 || zn.im >= 0.0 || M != HankelKind::Second)) {
            MR = -MM;
            if !(zn.re != 0.0 || zn.im >= 0.0) {
                zn = -zn;
            }
        }
        let (cy, NW) = ZBUNK(zn, order, KODE, MR, NN)?;
        NZ += NW;
        (cy, NZ)
    };
    //-----------------------------------------------------------------------
    //     H(M,FNU,Z) = -FMM*(I/FRAC_PI_2)*(ZT**FNU)*K(FNU,-Z*ZT)
    //
    //     ZT=EXP(-FMM*FRAC_PI_2*I) = CMPLX(0.0,-FMM), FMM=3-2*M, M=1,2
    //-----------------------------------------------------------------------
    let sign = -FRAC_PI_2 * FMM.signum();
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
    for element in cy.iter_mut().take(NN) {
        let ATOL = if max_abs_component(*element) < MACHINE_CONSTANTS.absolute_approximation_limit {
            *element *= MACHINE_CONSTANTS.rtol;
            MACHINE_CONSTANTS.abs_error_tolerance
        } else {
            1.0
        };
        *element *= csgn * ATOL;
        csgn *= Complex64::I * -FMM;
    }
    if partial_loss_of_significance {
        Err(BesselError::PartialLossOfSignificance { y: cy, nz: NZ })
    } else {
        Ok((cy, NZ))
    }
}

pub fn zbesi(z: Complex64, order: f64, KODE: Scaling, N: usize) -> BesselResult {
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

    let mut err = None;
    if order < 0.0_f64 {
        err = Some("order must be positive")
    };
    if N < 1 {
        err = Some("N must be >= 1")
    };
    if let Some(details) = err {
        return Err(BesselError::InvalidInput {
            details: details.to_owned(),
        });
    }

    //-----------------------------------------------------------------------------
    //     TEST FOR PROPER RANGE
    //-----------------------------------------------------------------------
    let AZ = z.abs();
    let FN = order + ((N - 1) as f64);
    let partial_significance_loss = is_sigificance_lost(AZ, FN, false)?;

    let (zn, mut csgn) = if z.re >= 0.0 {
        (z, c_one())
    } else {
        //-----------------------------------------------------------------------
        //     CALCULATE CSGN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
        //     WHEN FNU IS LARGE
        //-----------------------------------------------------------------------
        let INU = order as usize;
        let ARG = order.fract() * PI * if z.im < 0.0 { -1.0 } else { 1.0 };
        let mut csgn = Complex64::cis(ARG);
        if !INU.is_multiple_of(2) {
            csgn = -csgn;
        }
        (-z, csgn)
    };
    let (mut cy, nz) = ZBINU(zn, order, KODE, N)?;
    let remaining_n = N - nz;
    if z.re < 0.0 && remaining_n > 0 {
        //-----------------------------------------------------------------------
        //     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE
        //-----------------------------------------------------------------------
        for element in cy.iter_mut().take(remaining_n) {
            let correction =
                if max_abs_component(*element) <= MACHINE_CONSTANTS.absolute_approximation_limit {
                    *element *= MACHINE_CONSTANTS.rtol;
                    MACHINE_CONSTANTS.abs_error_tolerance
                } else {
                    1.0
                };
            *element *= csgn;
            *element *= correction;
            csgn = -csgn;
        }
    }

    if partial_significance_loss {
        Err(PartialLossOfSignificance { y: cy, nz })
    } else {
        Ok((cy, nz))
    }
}

pub fn zbesj(
    z: Complex64, //ZR, ZI,
    order: f64,   //FNU,
    KODE: Scaling,
    n: usize,
) -> BesselResult {
    // ouputs: CYR, CYI, NZ, IERR)
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

    let mut err = None;
    if order < 0.0_f64 {
        err = Some("order must be positive")
    };
    if n < 1 {
        err = Some("N must be >= 1")
    };
    if let Some(details) = err {
        return Err(BesselError::InvalidInput {
            details: details.to_owned(),
        });
    }

    //-----------------------------------------------------------------------
    //     TEST FOR PROPER RANGE
    //-----------------------------------------------------------------------
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
    let (mut cy, nz) = ZBINU(zn, order, KODE, n)?;
    for i in 0..n - nz {
        let mut cyi = cy[i];
        let mut ATOL = 1.0;
        // TODO is the below a pattern?
        if (max_abs_component(cyi)) <= MACHINE_CONSTANTS.absolute_approximation_limit {
            cyi *= MACHINE_CONSTANTS.rtol;
            ATOL = MACHINE_CONSTANTS.abs_error_tolerance;
        }
        let st = cyi * csgn;
        cy[i] = st * ATOL;
        csgn *= sign_selector * Complex64::I;
    }
    if partial_significance_loss {
        Err(PartialLossOfSignificance { y: cy, nz })
    } else {
        Ok((cy, nz))
    }
}

/*
fn ZBESK(ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, IERR)
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
//
// ***ROUTINES CALLED  ZACON,ZBKNU,ZBUNK,ZUOIK,ZABS,i1mach,d1mach
// ***END PROLOGUE  ZBESK
//
//     COMPLEX CY,Z
      EXTERNAL ZABS
      DOUBLE PRECISION AA, ALIM, ALN, ARG, AZ, CYI, CYR, DIG, ELIM, FN,
     * FNU, FNUL, RL, R1M5, TOL, UFL, ZI, ZR, d1mach, ZABS, BB
      INTEGER IERR, K, KODE, K1, K2, MR, N, NN, NUF, NW, NZ, i1mach
      DIMENSION CYR(N), CYI(N)
// ***FIRST EXECUTABLE STATEMENT  ZBESK
      IERR = 0
      NZ=0
      if (ZI == 0.0E0 && ZR == 0.0E0) IERR=1
      if (FNU < 0.0) IERR=1
      if (KODE < 1 || KODE > 2) IERR=1
      if (N < 1) IERR=1
      if (IERR != 0) RETURN
      NN = N
//-----------------------------------------------------------------------
//     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
//     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
//     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
//     EXP(-ELIM) < EXP(-ALIM)=EXP(-ELIM)/TOL    AND
//     EXP(ELIM) > EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
//     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
//     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
//     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
//     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU
//-----------------------------------------------------------------------
      TOL = DMAX1(d1mach(4),1.0e-18)
      K1 = i1mach(15)
      K2 = i1mach(16)
      R1M5 = d1mach(5)
      K = MIN0(K1.abs(),K2.abs())
      ELIM = 2.303*(K as f64)*R1M5-3.0)
      K1 = i1mach(14) - 1
      AA = R1M5*(K1 as f64)
      DIG = DMIN1(AA,18.0)
      AA = AA*2.303
      ALIM = ELIM + DMAX1(-AA,-41.45)
      FNUL = 10.0 + 6.0*(DIG-3.0)
      RL = 1.2*DIG + 3.0
//-----------------------------------------------------------------------------
//     TEST FOR PROPER RANGE
//-----------------------------------------------------------------------
      AZ = ZABS(ZR,ZI)
      FN = FNU + ((NN-1) as f64)
      AA = 0.5/TOL
      BB=DBLE(FLOAT(i1mach(9)))*0.5
      AA = DMIN1(AA,BB)
      if (AZ > AA) {return Err(LossOfSignificance);}
      if (FN > AA) {return Err(LossOfSignificance);}
      AA = DSQRT(AA)
      if (AZ > AA) IERR=3
      if (FN > AA) IERR=3
//-----------------------------------------------------------------------
//     OVERFLOW TEST ON THE LAST MEMBER OF THE SEQUENCE
//-----------------------------------------------------------------------
//     UFL = DEXP(-ELIM)
      UFL = d1mach(1)*1.0e+3
      if (AZ < UFL) GO TO 180
      if (FNU > FNUL) GO TO 80
      if (FN <= 1.0) GO TO 60
      if (FN > 2.0) GO TO 50
      if (AZ > TOL) GO TO 60
      ARG = 0.5*AZ
      ALN = -FN*DLOG(ARG)
      if (ALN > ELIM) GO TO 180
      GO TO 60
   50 CONTINUE
      CALL ZUOIK(ZR, ZI, FNU, KODE, 2, NN, CYR, CYI, NUF, TOL, ELIM,
     * ALIM)
      if (NUF < 0) GO TO 180
      NZ = NZ + NUF
      NN = NN - NUF
//-----------------------------------------------------------------------
//     HERE NN=N OR NN=0 SINCE NUF=0,NN, OR -1 ON RETURN FROM CUOIK
//     if NUF=NN, THEN CY(I)=CZERO FOR ALL I
//-----------------------------------------------------------------------
      if (NN == 0) GO TO 100
   60 CONTINUE
      if (ZR < 0.0) GO TO 70
//-----------------------------------------------------------------------
//     RIGHT HALF PLANE COMPUTATION, REAL(Z) >= 0.
//-----------------------------------------------------------------------
      CALL ZBKNU(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, TOL, ELIM, ALIM)
      if (NW < 0) GO TO 200
      NZ=NW
      RETURN
//-----------------------------------------------------------------------
//     LEFT HALF PLANE COMPUTATION
//     PI/2 < ARG(Z) <= PI AND -PI < ARG(Z) < -PI/2.
//-----------------------------------------------------------------------
   70 CONTINUE
      if (NZ != 0) GO TO 180
      MR = 1
      if (ZI < 0.0) MR = -1
      CALL ZACON(ZR, ZI, FNU, KODE, MR, NN, CYR, CYI, NW, RL, FNUL,
     * TOL, ELIM, ALIM)
      if (NW < 0) GO TO 200
      NZ=NW
      RETURN
//-----------------------------------------------------------------------
//     UNIFORM ASYMPTOTIC EXPANSIONS FOR FNU > FNUL
//-----------------------------------------------------------------------
   80 CONTINUE
      MR = 0
      if (ZR >= 0.0) GO TO 90
      MR = 1
      if (ZI < 0.0) MR = -1
   90 CONTINUE
      CALL ZBUNK(ZR, ZI, FNU, KODE, MR, NN, CYR, CYI, NW, TOL, ELIM,
     * ALIM)
      if (NW < 0) GO TO 200
      NZ = NZ + NW
      RETURN
  100 CONTINUE
      if (ZR < 0.0) GO TO 180
      RETURN
  180 CONTINUE
      NZ = 0
      IERR=2
      RETURN
  200 CONTINUE
      if(NW == (-1)) GO TO 180
      NZ=0
      IERR=5
      RETURN
  260 CONTINUE
      NZ=0
      IERR=4
      RETURN
      END
fn ZBESY(ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, CWRKR,
     *           CWRKI, IERR)
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
//
// ***ROUTINES CALLED  ZBESI,ZBESK,i1mach,d1mach
// ***END PROLOGUE  ZBESY
//
//     COMPLEX CWRK,CY,C1,C2,EX,HCI,Z,ZU,ZV
      DOUBLE PRECISION ARG, ASCLE, CIPI, CIPR, CSGNI, CSGNR, CSPNI,
     * CSPNR, CWRKI, CWRKR, CYI, CYR, D1M5, d1mach, ELIM, EXI, EXR, EY,
     * FNU, FFNU, FRAC_PI_2, RFRAC_PI_2, STR, STI, TAY, TOL, ATOL, RTOL, ZI, ZR,
     * ZNI, ZNR, ZUI, ZUR, ZVI, ZVR, ZZI, ZZR
      INTEGER I, IERR, IFNU, I4, K, KODE, K1, K2, N, NZ, NZ1, NZ2,
     * i1mach
      DIMENSION CYR(N), CYI(N), CWRKR(N), CWRKI(N), CIPR(4), CIPI(4)
      DATA CIPR(1),CIPR(2),CIPR(3),CIPR(4)/1.0, 0.0, -1.0, 0.0/
      DATA CIPI(1),CIPI(2),CIPI(3),CIPI(4)/0.0, 1.0, 0.0, -1.0/
      DATA FRAC_PI_2 / 1.57079632679489662 /
// ***FIRST EXECUTABLE STATEMENT  ZBESY
      IERR = 0
      NZ=0
      if (ZR == 0.0 && ZI == 0.0) IERR=1
      if (FNU < 0.0) IERR=1
      if (KODE < 1 || KODE > 2) IERR=1
      if (N < 1) IERR=1
      if (IERR != 0) RETURN
      ZZR = ZR
      ZZI = ZI
      if (ZI < 0.0) ZZI = -ZZI
      ZNR = ZZI
      ZNI = -ZZR
      CALL ZBESI(ZNR, ZNI, FNU, KODE, N, CYR, CYI, NZ1, IERR)
      if (IERR != 0&&IERR != 3) GO TO 90
      CALL ZBESK(ZNR, ZNI, FNU, KODE, N, CWRKR, CWRKI, NZ2, IERR)
      if (IERR != 0&&IERR != 3) GO TO 90
      NZ = MIN(NZ1,NZ2)
      IFNU = INT(SNGL(FNU))
      FFNU = FNU - (IFNU as f64)
      ARG = FRAC_PI_2*FFNU
      CSGNR = COS(ARG)
      CSGNI = SIN(ARG)
      I4 = MOD(IFNU,4) + 1
      STR = CSGNR*CIPR(I4) - CSGNI*CIPI(I4)
      CSGNI = CSGNR*CIPI(I4) + CSGNI*CIPR(I4)
      CSGNR = STR
      RFRAC_PI_2 = 1.0/FRAC_PI_2
      CSPNR = CSGNR*RFRAC_PI_2
      CSPNI = -CSGNI*RFRAC_PI_2
      STR = -CSGNI
      CSGNI = CSGNR
      CSGNR = STR
      if (KODE == 2) GO TO 60
      DO 50 I=1,N
//       CY(I) = CSGN*CY(I)-CSPN*CWRK(I)
        STR = CSGNR*CYR(I) - CSGNI*CYI(I)
        STR = STR - (CSPNR*CWRKR(I) - CSPNI*CWRKI(I))
        STI = CSGNR*CYI(I) + CSGNI*CYR(I)
        STI = STI - (CSPNR*CWRKI(I) + CSPNI*CWRKR(I))
        CYR(I) = STR
        CYI(I) = STI
        STR = - CSGNI
        CSGNI = CSGNR
        CSGNR = STR
        STR = CSPNI
        CSPNI = -CSPNR
        CSPNR = STR
   50 CONTINUE
      if (ZI < 0.0) THEN
        DO 55 I=1,N
          CYI(I) = -CYI(I)
   55   CONTINUE
      ENDIF
      RETURN
   60 CONTINUE
      EXR = COS(ZR)
      EXI = SIN(ZR)
      TOL = MAX(d1mach(4),1.0e-18)
      K1 = i1mach(15)
      K2 = i1mach(16)
      K = MIN(K1.abs(),K2.abs())
      D1M5 = d1mach(5)
//-----------------------------------------------------------------------
//     ELIM IS THE APPROXIMATE EXPONENTIAL UNDER- AND OVERFLOW LIMIT
//-----------------------------------------------------------------------
      ELIM = 2.303*(K as f64)*D1M5-3.0)
      EY = 0.0
      TAY = ABS(ZI+ZI)
      if (TAY < ELIM) EY = EXP(-TAY)
      STR = (EXR*CSPNR - EXI*CSPNI)*EY
      CSPNI = (EXR*CSPNI + EXI*CSPNR)*EY
      CSPNR = STR
      NZ = 0
      RTOL = 1.0/TOL
      ASCLE = d1mach(1)*RTOL*1.0e+3
      DO 80 I=1,N
//----------------------------------------------------------------------
//       CY(I) = CSGN*CY(I)-CSPN*CWRK(I): PRODUCTS ARE COMPUTED IN
//       SCALED MODE if CY(I) OR CWRK(I) ARE CLOSE TO UNDERFLOW TO
//       PREVENT UNDERFLOW IN AN INTERMEDIATE COMPUTATION.
//----------------------------------------------------------------------
        ZVR = CWRKR(I)
        ZVI = CWRKI(I)
        ATOL=1.0
        if (MAX(ABS(ZVR),ABS(ZVI)) > ASCLE) GO TO 75
          ZVR = ZVR*RTOL
          ZVI = ZVI*RTOL
          ATOL = TOL
   75   CONTINUE
        STR = (ZVR*CSPNR - ZVI*CSPNI)*ATOL
        ZVI = (ZVR*CSPNI + ZVI*CSPNR)*ATOL
        ZVR = STR
        ZUR = CYR(I)
        ZUI = CYI(I)
        ATOL=1.0
        if (MAX(ABS(ZUR),ABS(ZUI)) > ASCLE) GO TO 85
          ZUR = ZUR*RTOL
          ZUI = ZUI*RTOL
          ATOL = TOL
   85   CONTINUE
        STR = (ZUR*CSGNR - ZUI*CSGNI)*ATOL
        ZUI = (ZUR*CSGNI + ZUI*CSGNR)*ATOL
        ZUR = STR
        CYR(I) = ZUR - ZVR
        CYI(I) = ZUI - ZVI
        if (ZI < 0.0) CYI(I) = -CYI(I)
        if (CYR(I) == 0.0 && CYI(I) == 0.0 && EY == 0.0)
     &      NZ = NZ + 1
        STR = -CSGNI
        CSGNI = CSGNR
        CSGNR = STR
        STR = CSPNI
        CSPNI = -CSPNR
        CSPNR = STR
   80 CONTINUE
      RETURN
   90 CONTINUE
      NZ = 0
      RETURN
      END
      */

fn ZAIRY(
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
    //
    // ***ROUTINES CALLED  ZACAI,ZBKNU,ZEXP,ZSQRT,ZABS,i1mach,d1mach
    // ***END PROLOGUE  ZAIRY
    const C1: f64 = 3.55028053887817240e-01;
    const C2: f64 = 2.58819403792806799e-01;
    const COEFF: f64 = 1.83776298473930683e-01;

    let abs_z = z.abs();
    let float_is_derivative = if return_derivative { 1.0 } else { 0.0 };
    if abs_z <= 1.0 {
        //-----------------------------------------------------------------------
        //     POWER SERIES FOR CABS(Z) <= 1.
        //-----------------------------------------------------------------------
        let mut s1 = c_one();
        let mut s2 = c_one();
        if abs_z < MACHINE_CONSTANTS.abs_error_tolerance {
            s1 = c_zero();
            return if return_derivative {
                let mut ai = Complex64::new(-C2, 0.0);

                if abs_z > MACHINE_CONSTANTS.underflow_limit.sqrt() {
                    s1 = z.pow(2.0) / 2.0;
                }
                ai += C1 * s1;
                Ok((ai, 0))
            } else {
                if abs_z > MACHINE_CONSTANTS.underflow_limit {
                    s1 = C2 * z;
                }
                let ai = C1 - s1;
                Ok((ai, 0))
            };
        }
        let abs_z_sq = abs_z * abs_z;
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
        let mut ai = if return_derivative {
            let mut ai_inner = -s2 * C2;
            if abs_z > MACHINE_CONSTANTS.abs_error_tolerance {
                let CC = C1 / (1.0 + float_is_derivative);
                ai_inner += CC * z.pow(2.0) * s1;
            }
            ai_inner
        } else {
            s1 * C1 - C2 * z * s2
        };
        if scaling == Scaling::Scaled {
            ai *= (TWO_THIRDS * z * z.sqrt()).exp();
        }
        Ok((ai, 0))
    } else {
        //-----------------------------------------------------------------------
        //     CASE FOR CABS(Z) > 1.0
        //-----------------------------------------------------------------------
        let FNU = (1.0 + float_is_derivative) / 3.0;
        let ln_abs_z = abs_z.ln();
        //--------------------------------------------------------------------------
        //     TEST FOR PROPER RANGE
        //-----------------------------------------------------------------------
        // significance loss only tested against z, not order, so 0.0 is used to never cause significance loss
        let partial_loss_of_significance = is_sigificance_lost(abs_z, 0.0, true)?;

        let csq = z.sqrt();
        let mut zta = TWO_THIRDS * z * csq;
        //-----------------------------------------------------------------------
        //     RE(ZTA) <= 0 WHEN RE(Z) < 0, ESPECIALLY WHEN IM(Z) IS SMALL
        //-----------------------------------------------------------------------
        let mut overflow_state = Overflow::None;
        let mut scale_factor = 1.0;
        if z.re < 0.0 {
            zta.re = -zta.re.abs();
        }
        if z.im == 0.0 && z.re <= 0.0 {
            zta.re = 0.0;
        }
        let mut AA = zta.re;
        let (cy, NZ) = if !(AA >= 0.0 && z.re > 0.0) {
            //-----------------------------------------------------------------------
            //     OVERFLOW TEST
            //-----------------------------------------------------------------------
            if scaling == Scaling::Unscaled && AA <= -MACHINE_CONSTANTS.approximation_limit {
                AA = -AA + 0.25 * ln_abs_z;
                overflow_state = Overflow::NearOver;
                scale_factor = MACHINE_CONSTANTS.abs_error_tolerance;
                if AA > MACHINE_CONSTANTS.exponent_limit {
                    return Err(Overflow);
                }
            }
            //-----------------------------------------------------------------------
            //     CBKNU AND CACON RETURN EXP(ZTA)*K(FNU,ZTA) ON KODE=2
            //-----------------------------------------------------------------------
            let MR = if z.im < 0.0 { -1 } else { 1 };
            ZACAI(zta, FNU, scaling, MR, 1)?
        } else {
            //-----------------------------------------------------------------------
            //     UNDERFLOW TEST
            //-----------------------------------------------------------------------
            if scaling == Scaling::Unscaled && AA > MACHINE_CONSTANTS.approximation_limit {
                AA = -AA - 0.25 * ln_abs_z;
                overflow_state = Overflow::NearUnder;
                scale_factor = 1.0 / MACHINE_CONSTANTS.abs_error_tolerance;
                if AA < -MACHINE_CONSTANTS.exponent_limit {
                    return Ok((c_zero(), 1));
                }
            }
            ZBKNU(zta, FNU, scaling, 1)?
        };
        let mut s1 = cy[0] * COEFF;
        let retval = if overflow_state == Overflow::None {
            let ai = if return_derivative { -z * s1 } else { csq * s1 };
            (ai, NZ)
        } else {
            s1 *= scale_factor;
            s1 *= if return_derivative { -z } else { csq };
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
}
/*
fn ZBIRY(ZR, ZI, ID, KODE, BIR, BII, IERR)
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
//
// ***ROUTINES CALLED  ZBINU,ZABS,ZDIV,ZSQRT,d1mach,i1mach
// ***END PROLOGUE  ZBIRY
//     COMPLEX BI,CONE,CSQ,CY,S1,S2,TRM1,TRM2,Z,ZTA,Z3
      EXTERNAL ZABS
      DOUBLE PRECISION AA, AD, AK, ALIM, ATRM, AZ, AZ3, BB, BII, BIR,
     * BK, CC, CK, COEF, CONEI, CONER, CSQI, CSQR, CYI, CYR, C1, C2,
     * DIG, DK, D1, D2, EAA, ELIM, FID, FMR, FNU, FNUL, PI, RL, R1M5,
     * SFAC, STI, STR, S1I, S1R, S2I, S2R, TOL, TRM1I, TRM1R, TRM2I,
     * TRM2R, TWO_THIRDS, ZI, ZR, ZTAI, ZTAR, Z3I, Z3R, d1mach, ZABS
      INTEGER ID, IERR, K, KODE, K1, K2, NZ, i1mach
      DIMENSION CYR(2), CYI(2)
      DATA TWO_THIRDS, C1, C2, COEF, PI /6.66666666666666667e-01,
     * 6.14926627446000736e-01,4.48288357353826359e-01,
     * 5.77350269189625765e-01,3.14159265358979324e+00/
      DATA CONER, CONEI /1.0,0.0/
// ***FIRST EXECUTABLE STATEMENT  ZBIRY
      IERR = 0
      NZ=0
      if (ID < 0 || ID > 1) IERR=1
      if (KODE < 1 || KODE > 2) IERR=1
      if (IERR != 0) RETURN
      AZ = ZABS(ZR,ZI)
      TOL = DMAX1(d1mach(4),1.0e-18)
      FID = (ID as f64)
      if (AZ > 1.0E0) GO TO 70
//-----------------------------------------------------------------------
//     POWER SERIES FOR CABS(Z) <= 1.
//-----------------------------------------------------------------------
      S1R = CONER
      S1I = CONEI
      S2R = CONER
      S2I = CONEI
      if (AZ < TOL) GO TO 130
      AA = AZ*AZ
      if (AA < TOL/AZ) GO TO 40
      TRM1R = CONER
      TRM1I = CONEI
      TRM2R = CONER
      TRM2I = CONEI
      ATRM = 1.0
      STR = ZR*ZR - ZI*ZI
      STI = ZR*ZI + ZI*ZR
      Z3R = STR*ZR - STI*ZI
      Z3I = STR*ZI + STI*ZR
      AZ3 = AZ*AA
      AK = 2.0 + FID
      BK = 3.0 - FID - FID
      CK = 4.0 - FID
      DK = 3.0 + FID + FID
      D1 = AK*DK
      D2 = BK*CK
      AD = DMIN1(D1,D2)
      AK = 24.0 + 9.0*FID
      BK = 30.0 - 9.0*FID
      DO 30 K=1,25
        STR = (TRM1R*Z3R-TRM1I*Z3I)/D1
        TRM1I = (TRM1R*Z3I+TRM1I*Z3R)/D1
        TRM1R = STR
        S1R = S1R + TRM1R
        S1I = S1I + TRM1I
        STR = (TRM2R*Z3R-TRM2I*Z3I)/D2
        TRM2I = (TRM2R*Z3I+TRM2I*Z3R)/D2
        TRM2R = STR
        S2R = S2R + TRM2R
        S2I = S2I + TRM2I
        ATRM = ATRM*AZ3/AD
        D1 = D1 + AK
        D2 = D2 + BK
        AD = DMIN1(D1,D2)
        if (ATRM < TOL*AD) GO TO 40
        AK = AK + 18.0
        BK = BK + 18.0
   30 CONTINUE
   40 CONTINUE
      if (ID == 1) GO TO 50
      BIR = C1*S1R + C2*(ZR*S2R-ZI*S2I)
      BII = C1*S1I + C2*(ZR*S2I+ZI*S2R)
      if (KODE == 1) RETURN
      CALL ZSQRT(ZR, ZI, STR, STI)
      ZTAR = TWO_THIRDS*(ZR*STR-ZI*STI)
      ZTAI = TWO_THIRDS*(ZR*STI+ZI*STR)
      AA = ZTAR
      AA = -(AA).abs()
      EAA = DEXP(AA)
      BIR = BIR*EAA
      BII = BII*EAA
      RETURN
   50 CONTINUE
      BIR = S2R*C2
      BII = S2I*C2
      if (AZ <= TOL) GO TO 60
      CC = C1/(1.0+FID)
      STR = S1R*ZR - S1I*ZI
      STI = S1R*ZI + S1I*ZR
      BIR = BIR + CC*(STR*ZR-STI*ZI)
      BII = BII + CC*(STR*ZI+STI*ZR)
   60 CONTINUE
      if (KODE == 1) RETURN
      CALL ZSQRT(ZR, ZI, STR, STI)
      ZTAR = TWO_THIRDS*(ZR*STR-ZI*STI)
      ZTAI = TWO_THIRDS*(ZR*STI+ZI*STR)
      AA = ZTAR
      AA = -(AA).abs()
      EAA = DEXP(AA)
      BIR = BIR*EAA
      BII = BII*EAA
      RETURN
//-----------------------------------------------------------------------
//     CASE FOR CABS(Z) > 1.0
//-----------------------------------------------------------------------
   70 CONTINUE
      FNU = (1.0+FID)/3.0
//-----------------------------------------------------------------------
//     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
//     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
//     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
//     EXP(-ELIM) < EXP(-ALIM)=EXP(-ELIM)/TOL    AND
//     EXP(ELIM) > EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
//     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
//     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
//     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
//     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU.
//-----------------------------------------------------------------------
      K1 = i1mach(15)
      K2 = i1mach(16)
      R1M5 = d1mach(5)
      K = MIN0(K1.abs(),K2.abs())
      ELIM = 2.303*(K as f64)*R1M5-3.0)
      K1 = i1mach(14) - 1
      AA = R1M5*(K1 as f64)
      DIG = DMIN1(AA,18.0)
      AA = AA*2.303
      ALIM = ELIM + DMAX1(-AA,-41.45)
      RL = 1.2*DIG + 3.0
      FNUL = 10.0 + 6.0*(DIG-3.0)
//-----------------------------------------------------------------------
//     TEST FOR RANGE
//-----------------------------------------------------------------------
      AA=0.5/TOL
      BB=DBLE(FLOAT(i1mach(9)))*0.5
      AA=DMIN1(AA,BB)
      AA=AA**TWO_THIRDS
      if (AZ > AA) {return Err(LossOfSignificance);}
      AA=DSQRT(AA)
      if (AZ > AA) IERR=3
      CALL ZSQRT(ZR, ZI, CSQR, CSQI)
      ZTAR = TWO_THIRDS*(ZR*CSQR-ZI*CSQI)
      ZTAI = TWO_THIRDS*(ZR*CSQI+ZI*CSQR)
//-----------------------------------------------------------------------
//     RE(ZTA) <= 0 WHEN RE(Z) < 0, ESPECIALLY WHEN IM(Z) IS SMALL
//-----------------------------------------------------------------------
      SFAC = 1.0
      AK = ZTAI
      if (ZR >= 0.0) GO TO 80
      BK = ZTAR
      CK = -(BK).abs()
      ZTAR = CK
      ZTAI = AK
   80 CONTINUE
      if (ZI != 0.0 || ZR > 0.0) GO TO 90
      ZTAR = 0.0
      ZTAI = AK
   90 CONTINUE
      AA = ZTAR
      if (KODE == 2) GO TO 100
//-----------------------------------------------------------------------
//     OVERFLOW TEST
//-----------------------------------------------------------------------
      BB = (AA).abs()
      if (BB < ALIM) GO TO 100
      BB = BB + 0.25*DLOG(AZ)
      SFAC = TOL
      if (BB > ELIM) GO TO 190
  100 CONTINUE
      FMR = 0.0
      if (AA >= 0.0 && ZR > 0.0) GO TO 110
      FMR = PI
      if (ZI < 0.0) FMR = -PI
      ZTAR = -ZTAR
      ZTAI = -ZTAI
  110 CONTINUE
//-----------------------------------------------------------------------
//     AA=FACTOR FOR ANALYTIC CONTINUATION OF I(FNU,ZTA)
//     KODE=2 RETURNS EXP(-ABS(XZTA))*I(FNU,ZTA) FROM ZBESI
//-----------------------------------------------------------------------
      CALL ZBINU(ZTAR, ZTAI, FNU, KODE, 1, CYR, CYI, NZ, RL, FNUL, TOL,
     * ELIM, ALIM)
      if (NZ < 0) GO TO 200
      AA = FMR*FNU
      Z3R = SFAC
      STR = DCOS(AA)
      STI = DSIN(AA)
      S1R = (STR*CYR(1)-STI*CYI(1))*Z3R
      S1I = (STR*CYI(1)+STI*CYR(1))*Z3R
      FNU = (2.0-FID)/3.0
      CALL ZBINU(ZTAR, ZTAI, FNU, KODE, 2, CYR, CYI, NZ, RL, FNUL, TOL,
     * ELIM, ALIM)
      CYR(1) = CYR(1)*Z3R
      CYI(1) = CYI(1)*Z3R
      CYR(2) = CYR(2)*Z3R
      CYI(2) = CYI(2)*Z3R
//-----------------------------------------------------------------------
//     BACKWARD RECUR ONE STEP FOR ORDERS -1/3 OR -2/3
//-----------------------------------------------------------------------
      CALL ZDIV(CYR(1), CYI(1), ZTAR, ZTAI, STR, STI)
      S2R = (FNU+FNU)*STR + CYR(2)
      S2I = (FNU+FNU)*STI + CYI(2)
      AA = FMR*(FNU-1.0)
      STR = DCOS(AA)
      STI = DSIN(AA)
      S1R = COEF*(S1R+S2R*STR-S2I*STI)
      S1I = COEF*(S1I+S2R*STI+S2I*STR)
      if (ID == 1) GO TO 120
      STR = CSQR*S1R - CSQI*S1I
      S1I = CSQR*S1I + CSQI*S1R
      S1R = STR
      BIR = S1R/SFAC
      BII = S1I/SFAC
      RETURN
  120 CONTINUE
      STR = ZR*S1R - ZI*S1I
      S1I = ZR*S1I + ZI*S1R
      S1R = STR
      BIR = S1R/SFAC
      BII = S1I/SFAC
      RETURN
  130 CONTINUE
      AA = C1*(1.0-FID) + FID*C2
      BIR = AA
      BII = 0.0
      RETURN
  190 CONTINUE
      IERR=2
      NZ=0
      RETURN
  200 CONTINUE
      if(NZ == (-1)) GO TO 190
      NZ=0
      IERR=5
      RETURN
  260 CONTINUE
      IERR=4
      NZ=0
      RETURN
      END
*/

fn ZBKNU(z: Complex64, order: f64, KODE: Scaling, N: usize) -> BesselResult {
    // ***BEGIN PROLOGUE  ZBKNU
    // ***REFER TO  ZBESI,ZBESK,ZAIRY,ZBESH
    //
    //     ZBKNU COMPUTES THE K BESSEL FUNCTION IN THE RIGHT HALF Z PLANE.
    //
    // ***ROUTINES CALLED  gamma_ln,i1mach,d1mach,ZKSCL,ZSHCH,ZUunderflowCHK,ZABS,ZDIV,
    //                    ZEXP,ZLOG,ZMLT,ZSQRT
    // ***END PROLOGUE  ZBKNU
    const CTWOR: f64 = 2.0;
    const R1: f64 = 2.0;
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

    let CAZ = z.abs();
    let mut NZ = 0;
    let mut underflow_occurred = false;
    let mut overflow_state;
    let rz = 2.0 * z.conj() / CAZ.powi(2);
    let mut INU = (order + 0.5) as isize; // round to nearest int
    let DNU = order - (INU as f64); // signed fractional part (-0.5 < DNU < 0.5 )
    let DNU2 = if DNU.abs() > MACHINE_CONSTANTS.abs_error_tolerance {
        DNU * DNU
    } else {
        0.0
    };
    let mut skip_to_240 = false;
    let (mut s1, mut s2) = if (DNU.abs() != 0.5) && (CAZ <= R1) {
        //-----------------------------------------------------------------------;
        //     SERIES FOR CABS(Z) <= R1; and not half integer order
        //-----------------------------------------------------------------------;
        let mut FC = 1.0;
        let mut smu = rz.ln();
        let fmu = smu * DNU;
        let csh = fmu.sinh();
        let cch = fmu.cosh();
        if DNU != 0.0 {
            FC = DNU * PI;
            FC /= FC.sin();
            smu = csh / DNU;
        }
        //-----------------------------------------------------------------------;
        //     GAM(1-Z)*GAM(1+Z)=PI*Z/SIN(PI*Z), T1=1/GAM(1-DNU), T2=1/GAM(1+DNU);
        //-----------------------------------------------------------------------;
        let T2 = (-gamma_ln(1.0 + DNU).unwrap()).exp();
        let T1 = 1.0 / (T2 * FC);

        let G1 = if (DNU).abs() <= 0.1 {
            //-----------------------------------------------------------------------;
            //     SERIES FOR F0 TO RESOLVE INDETERMINACY FOR SMALL ABS(DNU);
            //-----------------------------------------------------------------------;
            let mut ak = 1.0;
            let mut sum = CC[0];
            for cc in CC[1..].iter() {
                ak *= DNU2;
                let TM = cc * ak;
                sum += TM;
                if TM.abs() < MACHINE_CONSTANTS.abs_error_tolerance {
                    break;
                }
            }
            -sum
        } else {
            (T1 - T2) / (DNU + DNU)
        };
        let G2 = (T1 + T2) * 0.5;
        let mut f = FC * (G1 * cch + G2 * smu);
        let mut p = 0.5 * fmu.exp() / T2;
        let mut q = (0.5 / fmu.exp()) / T1;
        let mut s1 = f;
        let mut s2 = p;
        let mut AK = 1.0;
        let mut A1 = 1.0;
        let mut ck = c_one();
        let mut BK = 1.0 - DNU2;
        if INU == 0 && N == 1 {
            //-----------------------------------------------------------------------;
            //     SPECIAL CASE
            //     GENERATE K(FNU,Z), 0.0  <=  FNU  <  0.5 AND N=1;
            //-----------------------------------------------------------------------;
            if CAZ >= MACHINE_CONSTANTS.abs_error_tolerance {
                let cz = 0.25 * z.powu(2);
                let T1 = 0.25 * CAZ * CAZ;
                '_l60: loop {
                    f = (f * AK + p + q) / BK;
                    p /= AK - DNU;
                    q /= AK + DNU;
                    ck *= cz / AK;
                    s1 += ck * f;
                    A1 *= T1 / AK;
                    BK += (2.0 * AK) + 1.0;
                    AK += 1.0;
                    if A1 <= MACHINE_CONSTANTS.abs_error_tolerance {
                        break;
                    }
                }
            }
            let mut y = s1;
            if KODE == Scaling::Scaled {
                y *= z.exp();
            }
            return Ok((vec![y], NZ));
            //-----------------------------------------------------------------------;
            //     GENERATE K(DNU,Z) AND K(DNU+1,Z) FOR FORWARD RECURRENCE;
            //-----------------------------------------------------------------------;
        }

        if CAZ >= MACHINE_CONSTANTS.abs_error_tolerance {
            let cz = 0.25 * z.powu(2);
            let T1 = 0.25 * CAZ * CAZ;
            '_l90: loop {
                f = (f * AK + p + q) / BK;
                p /= AK - DNU;
                q /= AK + DNU;
                ck *= cz / AK;
                s1 += ck * f;

                // ... TODO this is the only bit that differs from the loop above...
                s2 += ck * (p - AK * f);
                A1 *= T1 / AK;
                BK += AK + AK + 1.0;
                AK += 1.0;

                if A1 <= MACHINE_CONSTANTS.abs_error_tolerance {
                    break;
                }
            }
        }
        AK = (order + 1.0) * smu.re.abs();
        overflow_state = if AK > MACHINE_CONSTANTS.approximation_limit {
            Overflow::NearOver
        } else {
            Overflow::None
        };
        s2 *= MACHINE_CONSTANTS.scaling_factors[overflow_state] * rz;
        s1 *= MACHINE_CONSTANTS.scaling_factors[overflow_state];
        if KODE == Scaling::Scaled {
            let z_exp = z.exp();
            s1 *= z_exp;
            s2 *= z_exp;
        }
        (s1, s2)
    } else {
        // alternative to SERIES FOR CABS(Z) <= R1; Or half integer order
        //-----------------------------------------------------------------------;
        //     underflow_occured=0 MEANS NO UNDERFLOW OCCURRED;
        //     underflow_occured=1 MEANS AN UNDERFLOW OCCURRED- COMPUTATION PROCEEDS WITH;
        //     KODED=2 AND A TEST FOR ON SCALE VALUES IS MADE DURING FORWARD;
        //     RECURSION;
        //-----------------------------------------------------------------------;
        let mut coef = Complex64::new(RTFRAC_PI_2, 0.0) / z.sqrt();
        overflow_state = Overflow::None;
        if KODE == Scaling::Unscaled {
            if z.re > MACHINE_CONSTANTS.approximation_limit {
                underflow_occurred = true;
                overflow_state = Overflow::NearUnder;
            } else {
                coef *= MACHINE_CONSTANTS.scaling_factors[overflow_state] * (-z).exp();
            }
        }
        let mut AK = (DNU * PI).cos().abs();
        let mut FHS = (0.25 - DNU2).abs();

        if DNU.abs() == 0.5 || AK == 0.0 || FHS == 0.0 {
            (coef, coef)
        } else {
            //-----------------------------------------------------------------------;
            //     MILLER ALGORITHM FOR CABS(Z) > R1;
            //-----------------------------------------------------------------------;
            //-----------------------------------------------------------------------;
            //     COMPUTE R2=F(E). if CABS(Z) >= R2, USE FORWARD RECURRENCE TO;
            //     DETERMINE THE BACKWARD INDEX K. R2=F(E) IS A STRAIGHT LINE ON;
            //     12 <= E <= 60. E IS COMPUTED FROM 2**(-E)=B**(1-i1mach(14))=;
            //     TOL WHERE B IS THE BASE OF THE ARITHMETIC.;
            //-----------------------------------------------------------------------;
            // TODO is this a pattern?
            let mut T1 = ((f64::MANTISSA_DIGITS - 1) as f64
                * (f64::RADIX as f64).log10()
                * std::f64::consts::LOG2_10)
                .clamp(12.0, 60.0);
            let T2 = TWO_THIRDS * T1 - 6.0;
            T1 = if z.re == 0.0 {
                FRAC_PI_2
            } else {
                (z.im / z.re).atan().abs()
            };

            let (FK, FHS) = if T2 <= CAZ {
                //-----------------------------------------------------------------------;
                //     FORWARD RECURRENCE LOOP WHEN CABS(Z) >= R2;
                //-----------------------------------------------------------------------;
                let ETEST = AK / (PI * CAZ * MACHINE_CONSTANTS.abs_error_tolerance);
                let mut FK = 1.0;
                if ETEST >= 1.0 {
                    let mut FKS = CTWOR;
                    let mut CKR = CAZ + CAZ + CTWOR;
                    let mut p1 = 0.0;
                    let mut p2 = 1.0;
                    let mut converged = false;
                    for _ in 0..KMAX {
                        let AK = FHS / FKS;
                        let CBR = CKR / (FK + 1.0);
                        let pt = p2;
                        p2 = CBR * p2 - AK * p1;
                        p1 = pt;
                        CKR += CTWOR;
                        FKS += FK + FK + CTWOR;
                        FHS += FK + FK;
                        FK += 1.0;
                        if ETEST < p2.abs() * FK {
                            converged = true;
                            break;
                        }
                    }
                    if !converged {
                        return Err(DidNotConverge);
                    }
                    FK += SPI * T1 * (T2 / CAZ).sqrt();
                    FHS = (0.25 - DNU2).abs();
                }
                (FK, FHS)
            } else {
                //-----------------------------------------------------------------------;
                //     COMPUTE BACKWARD INDEX K FOR CABS(Z) < R2;
                //-----------------------------------------------------------------------;
                AK *= FPI / (MACHINE_CONSTANTS.abs_error_tolerance * CAZ.sqrt().sqrt());
                let AA = 3.0 * T1 / (1.0 + CAZ);
                let BB = 14.7 * T1 / (28.0 + CAZ);
                AK = (AK.ln() + CAZ * AA.cos() / (1.0 + 0.008 * CAZ)) / BB.cos();
                let FK = 0.12125 * AK * AK / CAZ + 1.5;
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
            for i in (1..=K).rev() {
                let k_f64 = i as f64;
                let cb = (z + k_f64) * 2.0 / (k_f64 + 1.0);
                let pt = p2;
                p2 = (p2 * cb - p1) * (k_squared + k_f64) / (k_squared - k_f64 + FHS);
                p1 = pt;
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
            s1 *= coef * cs;
            if INU <= 0 && N <= 1 {
                skip_to_240 = true;
            } else {
                //-----------------------------------------------------------------------;
                //     COMPUTE P1/P2=(P1/CABS(P2)*CONJG(P2)/CABS(P2) FOR SCALING;
                //-----------------------------------------------------------------------;
                p1 /= p2.abs();
                p2 = p2.conj() / p2.abs();
                s2 = (((DNU + 0.5 - (p1 * p2)) / z) + 1.0) * s1;
            }
            (s1, s2)
        }
    };

    // Now s1, s2 set up, we can go to recurrence

    let zd = z;
    let mut INUB = 1;

    //-----------------------------------------------------------------------
    //     FORWARD RECURSION ON THE THREE TERM RECURSION WITH RELATION WITH
    //     SCALING NEAR EXPONENT EXTREMES ON KFLAG=1 OR KFLAG=3
    //-----------------------------------------------------------------------
    let mut ck = (DNU + 1.0) * rz;
    if N == 1 {
        INU -= 1
    };

    'l225: loop {
        if !skip_to_240 {
            if INU > 0 {
                if underflow_occurred {
                    underflow_occurred = false;
                    //-----------------------------------------------------------------------;
                    //     underflow_occured=1 CASES, FORWARD RECURRENCE ON SCALED VALUES ON UNDERFLOW;
                    //-----------------------------------------------------------------------;
                    let mut cy = c_zeros(2);
                    let HELIM = 0.5 * MACHINE_CONSTANTS.exponent_limit;
                    let ELM = (-MACHINE_CONSTANTS.exponent_limit).exp();
                    let CELMR = ELM;
                    let ASCLE = MACHINE_CONSTANTS.bry[0];
                    let mut zd = z;
                    let mut IC: isize = -1;
                    let mut J = 1;
                    let mut I = 0;
                    for i in 0..INU {
                        I = i + 1;
                        let st = s2;
                        s2 = s2 * ck + s1;
                        s1 = st;
                        ck += rz;
                        let ALAS = s2.abs().ln();
                        if -zd.re + ALAS >= -MACHINE_CONSTANTS.exponent_limit {
                            let p2 = -zd + s2.ln();
                            let p1 = (p2.re.exp() / MACHINE_CONSTANTS.abs_error_tolerance)
                                * Complex64::cis(p2.im);
                            if !will_z_underflow(p1, ASCLE, MACHINE_CONSTANTS.abs_error_tolerance) {
                                J = 1 - J;
                                cy[J] = p1;
                                // IF(IC.EQ.(I-1)) GO TO 264
                                // below implies we got here twice in a row
                                if IC == I - 1 {
                                    // underflow_occurred = true; //implies 270
                                    break;
                                } else {
                                    IC = I;
                                    continue;
                                }
                            }
                            if ALAS < HELIM {
                                continue;
                            }
                            zd.re -= MACHINE_CONSTANTS.exponent_limit;
                            s1 *= CELMR;
                            s2 *= CELMR;
                        }
                    }
                    overflow_state = Overflow::NearUnder;
                    INUB = I + 1;
                    s2 = cy[J];
                    J = 1 - J;
                    s1 = cy[J];
                    if INUB <= INU {
                        continue 'l225;
                    }
                } else {
                    let mut P1R = MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
                    let mut ASCLE = MACHINE_CONSTANTS.bry[overflow_state];
                    for _ in INUB..=INU {
                        let st = s2;
                        s2 = ck * s2 + s1;
                        s1 = st;
                        ck += rz;
                        if overflow_state == Overflow::NearOver {
                            continue;
                        }
                        let p2 = s2 * P1R;
                        if max_abs_component(p2) <= ASCLE {
                            continue;
                        }
                        overflow_state.increment();
                        ASCLE = MACHINE_CONSTANTS.bry[overflow_state];
                        s1 *= P1R;
                        s2 = p2;
                        s1 *= MACHINE_CONSTANTS.scaling_factors[overflow_state];
                        s2 *= MACHINE_CONSTANTS.scaling_factors[overflow_state];
                        P1R = MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
                    }
                }
            }
            if N == 1 {
                s1 = s2;
            }
        }
        // ********* basic setup
        let (mut KK, mut y) = if !underflow_occurred {
            let mut y = c_zeros(N);
            y[0] = s1 * MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
            if N == 1 {
                return Ok((y, NZ));
            }
            y[1] = s2 * MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
            if N == 2 {
                return Ok((y, NZ));
            }
            let KK = 1;
            (KK, y)
        // ********* End Basic Setup
        } else {
            //Complex setup from 270 onwards
            // ******** Alternative setup if underflow_occured
            let mut y = c_zeros(N);
            y[0] = s1;
            if N > 1 {
                y[1] = s2;
            }
            ZKSCL(zd, order, N, &mut y, &mut NZ, rz, MACHINE_CONSTANTS.bry[0]);
            INU = (N - NZ) as isize;
            if INU <= 0 {
                return Ok((y, NZ));
            }
            let mut KK = NZ; // + 1;
            s1 = y[KK];
            y[KK] *= MACHINE_CONSTANTS.reciprocal_scaling_factors[0];
            if INU == 1 {
                return Ok((y, NZ));
            }
            KK = NZ + 1;
            s2 = y[KK];
            y[KK] *= MACHINE_CONSTANTS.reciprocal_scaling_factors[0];
            if INU == 2 {
                return Ok((y, NZ));
            }
            ck = (order + (KK as f64)) * rz;
            overflow_state = Overflow::NearUnder;
            (KK, y)
        };
        KK += 1;
        if KK >= N {
            return Ok((y, NZ));
        }
        let mut P1R = MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
        let mut ASCLE = MACHINE_CONSTANTS.bry[overflow_state];
        let mut I;
        for i in KK..N {
            I = i;
            let mut p2 = s2;
            s2 = ck * s2 + s1;
            s1 = p2;
            ck += rz;
            p2 = s2 * P1R;
            y[I] = p2;
            if overflow_state == Overflow::NearOver {
                continue;
            };
            if max_abs_component(p2) <= ASCLE {
                continue;
            }
            overflow_state.increment();
            ASCLE = MACHINE_CONSTANTS.bry[overflow_state];
            s1 *= P1R;
            s2 = p2;
            s1 *= MACHINE_CONSTANTS.scaling_factors[overflow_state];
            s2 *= MACHINE_CONSTANTS.scaling_factors[overflow_state];
            P1R = MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
        }
        return Ok((y, NZ));
    }
}

fn ZKSCL(
    //ZRR,ZRI,FNU,
    zr: Complex64,
    order: f64,
    N: usize, //YR,YI,NZ,
    y: &mut [Complex64],
    NZ: &mut usize,
    rz: Complex64, //RZR,RZI,
    ASCLE: f64,
) //-> BesselResult//ASCLE,TOL,ELIM)
{
    // ***BEGIN PROLOGUE  ZKSCL
    // ***REFER TO  ZBESK
    //
    //     SET K FUNCTIONS TO ZERO ON UNDERFLOW, CONTINUE RECURRENCE
    //     ON SCALED FUNCTIONS UNTIL TWO MEMBERS COME ON SCALE, THEN
    //     RETURN WITH MIN(NZ+2,N) VALUES SCALED BY 1/TOL.

    *NZ = 0;
    let NN = min(2, N);
    let mut cy = c_zeros(2);
    let mut IC = 0;
    for i in 0..NN {
        let s1 = y[i];
        cy[i] = s1;
        *NZ += 1;
        y[i] = c_zero();
        if -zr.re + s1.abs().ln() < (-MACHINE_CONSTANTS.exponent_limit) {
            continue;
        }
        let mut cs = s1.ln() - zr;
        cs = cs.exp() / MACHINE_CONSTANTS.abs_error_tolerance;

        if will_z_underflow(cs, ASCLE, MACHINE_CONSTANTS.abs_error_tolerance) {
            continue;
        }
        y[i] = cs;
        IC = i;
        *NZ -= 1;
    }
    if N == 1 {
        return;
    }
    if IC < 1 {
        y[0] = c_zero();
        *NZ = 2;
    }
    if N == 2 {
        return;
    }
    if *NZ == 0 {
        return;
    }
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
    for i in 2..N {
        I = i;
        let mut cs = s2;
        s2 = cs * ck + s1;
        s1 = cs;
        ck += rz;
        let ALAS = s2.abs().ln();
        *NZ += 1;
        y[i] = Complex64::zero();
        if -zd.re + s2.abs().ln() >= -MACHINE_CONSTANTS.exponent_limit {
            cs = s2.ln() - zd;
            cs = cs.exp() / MACHINE_CONSTANTS.abs_error_tolerance;
            if !will_z_underflow(cs, ASCLE, MACHINE_CONSTANTS.abs_error_tolerance) {
                y[i] = cs;
                *NZ -= 1;
                if IC == i - 1 {
                    skip_to_40 = true;
                    break;
                }
                IC = i;
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
        *NZ = N;
        if IC == N {
            *NZ = N - 1
        };
    } else {
        *NZ = I - 2;
    }
    for element in y.iter_mut().take(*NZ) {
        *element = c_zero();
    }
}

fn i_ratios(z: Complex64, order: f64, N: usize) -> Vec<Complex64> {
    // ***BEGIN PROLOGUE  ZRATI
    // ***REFER TO  ZBESI,ZBESK,ZBESH
    //
    //     ZRATI COMPUTES RATIOS OF I BESSEL FUNCTIONS BY BACKWARD
    //     RECURRENCE.  THE STARTING INDEX IS DETERMINED BY FORWARD
    //     RECURRENCE AS DESCRIBED IN J. RES. OF NAT. BUR. OF STANDARDS-B,
    //     MATHEMATICAL SCIENCES, VOL 77B, P111-114, SEPTEMBER, 1973,
    //     BESSEL FUNCTIONS I AND J OF COMPLEX ARGUMENT AND INTEGER ORDER,
    //     BY D. J. SOOKNE.
    //
    // ***ROUTINES CALLED  ZABS,ZDIV
    // ***END PROLOGUE  ZRATI
    let AZ = z.abs();
    let INU = order as usize;
    let IDNU = INU + N - 1;
    let MAGZ = AZ as isize;
    let AMAGZ = (MAGZ + 1) as f64;
    let FDNU = IDNU as f64;
    let FNUP = AMAGZ.max(FDNU);
    let mut ID = IDNU as isize - MAGZ - 1;
    let mut K = 1;
    let rz = 2.0 * z.conj() / AZ.powi(2);
    let mut t1 = rz * FNUP;
    let mut p2 = -t1;
    let mut p1 = c_one();
    t1 += rz;
    if ID > 0 {
        ID = 0;
    }
    let mut AP2 = p2.abs();
    let mut AP1 = p1.abs();
    //-----------------------------------------------------------------------
    //     THE OVERFLOW TEST ON K(FNU+I-1,Z) BEFORE THE CALL TO CBKNU
    //     GUARANTEES THAT P2 IS ON SCALE. SCALE TEST1 AND ALL SUBSEQUENT
    //     P2 VALUES BY AP1 TO ENSURE THAT AN OVERFLOW DOES NOT OCCUR
    //     PREMATURELY.
    //-----------------------------------------------------------------------
    let ARG = (AP2 + AP2) / (AP1 * MACHINE_CONSTANTS.abs_error_tolerance);
    let TEST1 = ARG.sqrt();
    let mut TEST = TEST1;
    p1 /= AP1;
    p2 /= AP1;
    AP2 /= AP1;
    let mut first_pass = true;
    'l10: loop {
        K += 1;
        AP1 = AP2;
        let pt = p2;
        p2 = p1 - (t1 * p2);
        p1 = pt;
        t1 += rz;
        AP2 = p2.abs();
        if AP1 <= TEST {
            continue;
        }
        if !first_pass {
            break 'l10;
        }
        let ak = t1.abs() / 2.0;
        let FLAM = ak + (ak.powi(2) - 1.0).sqrt();
        let rho = AP2 / AP1.min(FLAM);
        TEST = TEST1 * (rho / (rho.powi(2) - 1.0)).sqrt();
        first_pass = false;
    }
    let KK: usize = (K as isize + 1 - ID).try_into().unwrap();
    let mut t1 = Complex64::new(KK as f64, 0.0);
    let DFNU = order + ((N - 1) as f64);
    p1 = Complex64::new(1.0 / AP2, 0.0);
    p2 = c_zero();
    for _ in 0..KK {
        let pt = p1;
        let tt = rz * (DFNU + t1.re);
        p1 = p1 * tt + p2;
        p2 = pt;
        t1.re -= 1.0;
    }
    if p1.re == 0.0 && p1.im == 0.0 {
        p1 = Complex64::new(
            MACHINE_CONSTANTS.abs_error_tolerance,
            MACHINE_CONSTANTS.abs_error_tolerance,
        );
    }
    let mut cy = c_zeros(N);
    cy[N - 1] = p2 / p1;
    if N > 1 {
        t1 = Complex64::new((N - 1) as f64, 0.0);
        let cdfnu = order * rz;
        for k in (1..N).rev() {
            let mut pt = cdfnu + t1 * rz + cy[k];
            let mut AK = pt.abs();
            if AK == 0.0 {
                pt = Complex64::new(
                    MACHINE_CONSTANTS.abs_error_tolerance,
                    MACHINE_CONSTANTS.abs_error_tolerance,
                );
                AK = pt.abs();
            }
            cy[k - 1] = pt.conj() / AK.powi(2);
            t1 -= 1.0;
        }
    }
    cy
}

fn ZS1S2(zr: Complex64, s1: &mut Complex64, s2: &mut Complex64, IUF: &mut isize) -> usize {
    // ***BEGIN PROLOGUE  ZS1S2
    // ***REFER TO  ZBESK,ZAIRY
    //
    //     ZS1S2 TESTS FOR A POSSIBLE UNDERFLOW RESULTING FROM THE
    //     ADDITION OF THE I AND K FUNCTIONS IN THE ANALYTIC CON-
    //     TINUATION FORMULA WHERE S1=K FUNCTION AND S2=I FUNCTION.
    //     ON KODE=1 THE I AND K FUNCTIONS ARE DIFFERENT ORDERS OF
    //     MAGNITUDE, BUT FOR KODE=2 THEY CAN BE OF THE SAME ORDER
    //     OF MAGNITUDE AND THE MAXIMUM MUST BE AT LEAST ONE
    //     PRECISION ABOVE THE UNDERFLOW LIMIT.

    let NZ = 0;
    let mut abs_s1 = s1.abs();
    let abs_s2 = s2.abs();
    if (s1.re != 0.0 || s1.im != 0.0) && (abs_s1 != 0.0) {
        let ALN = (-2.0 * zr.re) + abs_s1.ln();
        let s1d = *s1;
        *s1 = c_zero();
        abs_s1 = 0.0;
        if ALN >= (-MACHINE_CONSTANTS.approximation_limit) {
            *s1 = (s1d.ln() - 2.0 * zr).exp();
            abs_s1 = s1.abs();
            *IUF += 1;
        }
    }
    if abs_s1.max(abs_s2) > MACHINE_CONSTANTS.absolute_approximation_limit {
        NZ
    } else {
        *s1 = c_zero();
        *s2 = c_zero();
        *IUF = 0;
        1
    }
}

fn ZBUNK(z: Complex64, order: f64, KODE: Scaling, MR: i64, N: usize) -> BesselResult {
    //ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, TOL, ELIM,
    // * ALIM)
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
        // ZUNK2(z, order, KODE, MR, N)
        Err(BesselError::NotYetImplemented)
        // ZUNK2(ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, TOL, ELIM, ALIM)
    } else {
        //-----------------------------------------------------------------------
        //     ASYMPTOTIC EXPANSION FOR K(FNU,Z) FOR LARGE FNU APPLIED IN
        //     -PI/3 <= ARG(Z) <= PI/3
        //-----------------------------------------------------------------------
        ZUNK1(z, order, KODE, MR, N)
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
    let (cw, _) = ZBKNU(zr, order, KODE, 2)?;
    let y_ratios = i_ratios(zr, order, N);
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
    KODE: Scaling,
    MR: i64,
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
    let (mut y, _) = ZBINU(zn, order, KODE, N)?;
    //-----------------------------------------------------------------------
    //     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
    //-----------------------------------------------------------------------
    let NN = 2.min(N);
    let (cy, NW) = ZBKNU(zn, order, KODE, NN)?;
    if NW > 0 {
        return Err(Overflow);
        // the NW = -1 or -2 is handled by ZBNKU returning an error,
        // but the amos code defaults to an overflow, if NW != 0
    }
    let mut s1 = cy[0];
    let FMR = MR as f64;
    let SGN = -PI * FMR.signum();
    let mut csgn = Complex64::new(0.0, SGN);
    if KODE == Scaling::Scaled {
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
    if KODE == Scaling::Scaled {
        let NW = ZS1S2(zn, &mut c1, &mut c2, &mut IUF);
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
    if KODE == Scaling::Scaled {
        let NW = ZS1S2(zn, &mut c1, &mut c2, &mut IUF);
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
    let mut overflow_state = if abs_s2 <= MACHINE_CONSTANTS.bry[0] {
        Overflow::NearUnder
    } else if abs_s2 > MACHINE_CONSTANTS.bry[1] {
        Overflow::NearOver
    } else {
        Overflow::None
    };
    let mut b_scale = MACHINE_CONSTANTS.bry[overflow_state];
    s1 *= MACHINE_CONSTANTS.scaling_factors[overflow_state];
    s2 *= MACHINE_CONSTANTS.scaling_factors[overflow_state];
    let mut CSR = MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state];
    for i in 2..N {
        let st = s2;
        s2 = ck * s2 + s1;
        s1 = st;
        c1 = s2 * CSR;
        let mut st = c1;
        c2 = y[i];
        if KODE == Scaling::Scaled && IUF >= 0 {
            let NW = ZS1S2(zn, &mut c1, &mut c2, &mut IUF);
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
        y[i] = cspn * c1 + csgn * c2;
        ck += rz;
        cspn = -cspn;
        if overflow_state == Overflow::NearOver {
            continue;
        }
        if max_abs_component(c1) > b_scale {
            overflow_state.increment();
            b_scale = MACHINE_CONSTANTS.bry[overflow_state];
            s1 *= CSR;
            s2 = st;
            s1 *= MACHINE_CONSTANTS.scaling_factors[overflow_state];
            s2 *= MACHINE_CONSTANTS.scaling_factors[overflow_state];
            CSR = MACHINE_CONSTANTS.reciprocal_scaling_factors[overflow_state]; //CSRR(KFLAG);
        }
    }
    Ok((y, NZ))
}

fn ZBINU(z: Complex64, order: f64, KODE: Scaling, N: usize) -> BesselResult {
    // ***BEGIN PROLOGUE  ZBINU
    // ***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZAIRY,ZBIRY
    //
    //     ZBINU COMPUTES THE I FUNCTION IN THE RIGHT HALF Z PLANE
    //

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
        let INW: usize = NW.abs().try_into().unwrap();
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
        let (cy, nw) = z_asymptotic_i(z, order, KODE, NN)?;
        debug_assert!(nw == NZ);
        return Ok((cy, NZ));
    }
    let mut skip_az_rl_check = true;
    if DFNU > 1.0 {
        skip_az_rl_check = false;
        //-----------------------------------------------------------------------
        //     OVERFLOW AND UNDERFLOW TEST ON I SEQUENCE FOR MILLER ALGORITHM
        //-----------------------------------------------------------------------
        let nw = zuoik(z, order, KODE, IKType::I, NN, &mut cy)?;
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
    if let Ok(NW) = zuoik(z, order, KODE, IKType::K, 2, &mut [c_one(); 2]) {
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

fn ZACAI(z: Complex64, order: f64, KODE: Scaling, MR: i64, N: usize) -> BesselResult {
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
        z_asymptotic_i(zn, order, KODE, NN)?
    //-----------------------------------------------------------------------
    //     MILLER ALGORITHM NORMALIZED BY THE SERIES FOR THE I FUNCTION
    //-----------------------------------------------------------------------
    } else {
        i_miller(zn, order, KODE, NN)?
    };
    //-----------------------------------------------------------------------
    //     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
    //-----------------------------------------------------------------------s
    let (cy, _) = ZBKNU(zn, order, KODE, 1)?;
    let SGN = -PI * (MR as f64).signum();
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
        let NW = ZS1S2(zn, &mut c1, &mut c2, &mut IUF);
        NZ += NW;
    }
    y[0] = cspn * c1 + csgn * c2;
    Ok((y, NZ))
}

fn ZUNK1(z: Complex64, order: f64, scaling: Scaling, MR: i64, N: usize) -> BesselResult {
    // ***BEGIN PROLOGUE  ZUNK1
    // ***REFER TO  ZBESK
    //
    //     ZUNK1 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE
    //     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE
    //     UNIFORM ASYMPTOTIC EXPANSION.
    //     MR INDICATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION.
    //     NZ=-1 MEANS AN OVERFLOW WILL OCCUR
    //

    // TODO better name for KDFLAG
    let mut KDFLG = false; // = 1;
    let mut NZ = 0;
    //-----------------------------------------------------------------------
    //     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION GREATER THAN
    //     THE UNDERFLOW LIMIT
    //-----------------------------------------------------------------------
    let zr = if z.re < 0.0 { -z } else { z };
    let mut J = 1;
    let mut phi = c_zeros(2);
    let mut zeta1 = c_zeros(2);
    let mut zeta2 = c_zeros(2);
    let mut sum = c_zeros(2);
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
                let will_underflow = will_z_underflow(
                    s2,
                    MACHINE_CONSTANTS.bry[0],
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
        let (phi, zet1d, zet2d, _sumd) = zunik(zr, modified_order, IKType::K, MR == 0);
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
        let mut ASCLE = MACHINE_CONSTANTS.bry[k_overflow_state];
        for i in IB..N {
            (s1, s2) = (s2, ck * s2 + s1);
            ck += rz;
            y[i] = s2 * C1R;
            if k_overflow_state == Overflow::NearOver {
                continue;
            }
            if max_abs_component(y[i]) <= ASCLE {
                continue;
            }
            k_overflow_state.increment();
            ASCLE = MACHINE_CONSTANTS.bry[k_overflow_state];
            s1 *= C1R;
            s2 = y[i];
            s1 *= MACHINE_CONSTANTS.scaling_factors[k_overflow_state];
            s2 *= MACHINE_CONSTANTS.scaling_factors[k_overflow_state];
            C1R = MACHINE_CONSTANTS.reciprocal_scaling_factors[k_overflow_state];
        }
        if MR == 0 {
            return Ok((y, NZ));
        }
        //-----------------------------------------------------------------------
        //     ANALYTIC CONTINUATION FOR RE(Z) < 0.0
        //-----------------------------------------------------------------------
        NZ = 0;
        let FMR = MR as f64;
        let SGN = -PI * FMR.signum();
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
                        && will_z_underflow(
                            s2,
                            MACHINE_CONSTANTS.bry[0],
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
                let NW = ZS1S2(zr, &mut s1, &mut s2, &mut IUF);
                NZ += NW;
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
        let mut ASCLE = MACHINE_CONSTANTS.bry[i_overflow_state];
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
            ASCLE = MACHINE_CONSTANTS.bry[i_overflow_state];
            s1 *= csr;
            s2 = ck;
            s1 *= MACHINE_CONSTANTS.scaling_factors[i_overflow_state];
            s2 *= MACHINE_CONSTANTS.scaling_factors[i_overflow_state];
            csr = MACHINE_CONSTANTS.scaling_factors[i_overflow_state];
        }
    }
    Ok((y, NZ))
}

/*

fn ZUNK2(z: Complex64,order: f64,KODE: Scaling,MR: i64,N: usize,)->BesselResult{
    //ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, TOL, ELIM,
     // * ALIM)
// ***BEGIN PROLOGUE  ZUNK2
// ***REFER TO  ZBESK
//
//     ZUNK2 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE
//     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE
//     UNIFORM ASYMPTOTIC EXPANSIONS FOR H(KIND,FNU,ZN) AND J(FNU,ZN)
//     WHERE ZN IS IN THE RIGHT HALF PLANE, KIND=(3-MR)/2, MR=+1 OR
//     -1. HERE ZN=ZR*I OR -ZR*I WHERE ZR=Z if Z IS IN THE RIGHT
//     HALF PLANE OR ZR=-Z if Z IS IN THE LEFT HALF PLANE. MR INDIC-
//     ATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION.
//     NZ=-1 MEANS AN OVERFLOW WILL OCCUR
//
// ***ROUTINES CALLED  ZAIRY,ZKSCL,ZS1S2,ZUunderflowCHK,ZUNHJ,d1mach,ZABS
// ***END PROLOGUE  ZUNK2
//     COMPLEX AI,ARG,ARGD,ASUM,ASUMD,BSUM,BSUMD,CFN,CI,CIP,CK,CONE,CRSC,
//    *CR1,CR2,CS,CSCL,CSGN,CSPN,CSR,CSS,CY,CZERO,C1,C2,DAI,PHI,PHID,RZ,
//    *S1,S2,Y,Z,ZB,ZETA1,ZETA1D,ZETA2,ZETA2D,ZN,ZR
    //   EXTERNAL ZABS
    //   DOUBLE PRECISION AARG, AIC, AII, AIR, ALIM, ANG, APHI, ARGDI,
    //  * ARGDR, ARGI, ARGR, ASC, ASCLE, ASUMDI, ASUMDR, ASUMI, ASUMR,
    //  * BRY, BSUMDI, BSUMDR, BSUMI, BSUMR, CAR, CIPI, CIPR, CKI, CKR,
    //  * CONER, CRSC, CR1I, CR1R, CR2I, CR2R, CSCL, CSGNI, CSI,
    //  * CSPNI, CSPNR, CSR, CSRR, CSSR, CYI, CYR, C1I, C1R, C2I, C2M,
    //  * C2R, DAII, DAIR, ELIM, FMR, FN, FNF, FNU, FRAC_PI_2, PHIDI, PHIDR,
    //  * PHII, PHIR, PI, PTI, PTR, RAST, RAZR, RS1, RZI, RZR, SAR, SGN,
    //  * STI, STR, S1I, S1R, S2I, S2R, TOL, YI, YR, YY, ZBI, ZBR, ZEROI,
    //  * ZEROR, ZETA1I, ZETA1R, ZETA2I, ZETA2R, ZET1DI, ZET1DR, ZET2DI,
    //  * ZET2DR, ZI, ZNI, ZNR, ZR, ZRI, ZRR, d1mach, ZABS
    //   INTEGER I, IB, IFLAG, IFN, IL, IN, INU, IUF, K, KDFLG, KFLAG, KK,
    //  * KODE, MR, N, NAI, NDAI, NW, NZ, IDUM, J, IPARD, IC
    //   DIMENSION BRY(3), YR(N), YI(N), ASUMR(2), ASUMI(2), BSUMR(2),
    //  * BSUMI(2), PHIR(2), PHII(2), ARGR(2), ARGI(2), ZETA1R(2),
    //  * ZETA1I(2), ZETA2R(2), ZETA2I(2), CYR(2), CYI(2), CIPR(4),
    //  * CIPI(4), CSSR(3), CSRR(3)
    //   DATA ZEROR,ZEROI,CONER,CR1R,CR1I,CR2R,CR2I /
    //  1         0.0, 0.0, 1.0,
    // const CR1R = 1.73205080756887729
    // const CR1I = -0.5
     const CR1: Complex64 = Complex64::new(1.0, 1.73205080756887729);
     const CR2: Complex64 = Complex64::new( -0.5, -8.66025403784438647e-01 );
    // const CR1I = -0.5

    // -8.66025403784438647e-01 /
    //   DATA FRAC_PI_2, PI, AIC /
    //  1     1.57079632679489662e+00,     3.14159265358979324e+00,
     const AIC:f64 =      1.26551212348464539e+00;
    //   DATA CIPR(1),CIPI(1),CIPR(2),CIPI(2),CIPR(3),CIPI(3),CIPR(4),
    //  * CIPI(4) /
    //  1  1.0,0.0 ,  0.0,-1.0 ,  -1.0,0.0 ,  0.0,1.0 /
    const CIP:[Complex64;4] = [
        Complex64::new(1.0, 0.0),
        Complex64::new(0.0, -1.0),
        Complex64::new(-1.0, 0.0),
        Complex64::new(0.0, 1.0),
    ];
//
      KDFLG = 1;
      NZ = 0;
//-----------------------------------------------------------------------;
//     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION GREATER THAN;
//     THE UNDERFLOW LIMIT;
//-----------------------------------------------------------------------;
      CSCL = 1.0/TOL;
      CRSC = TOL;
      CSSR(1) = CSCL;
      CSSR(2) = CONER;
      CSSR(3) = CRSC;
      CSRR(1) = CRSC;
      CSRR(2) = CONER;
      CSRR(3) = CSCL;
      BRY(1) = 1.0e+3*d1mach(1)/TOL;
      BRY(2) = 1.0/BRY(1);
      BRY(3) = d1mach(2);
      ZRR = ZR;
      ZRI = ZI;
      if (ZR >= 0.0) GO TO 10;
      ZRR = -ZR;
      ZRI = -ZI;
   10 CONTINUE;
      YY = ZRI;
      ZNR = ZRI;
      ZNI = -ZRR;
      ZBR = ZRR;
      ZBI = ZRI;
      INU = INT(SNGL(FNU));
      FNF = FNU - (INU as f64);
      ANG = -FRAC_PI_2*FNF;
      CAR = DCOS(ANG);
      SAR = DSIN(ANG);
      C2R = FRAC_PI_2*SAR;
      C2I = -FRAC_PI_2*CAR;
      KK = MOD(INU,4) + 1;
      STR = C2R*CIPR(KK) - C2I*CIPI(KK);
      STI = C2R*CIPI(KK) + C2I*CIPR(KK);
      CSR = CR1R*STR - CR1I*STI;
      CSI = CR1R*STI + CR1I*STR;
      if (YY > 0.0) GO TO 20;
      ZNR = -ZNR;
      ZBI = -ZBI;
   20 CONTINUE;
//-----------------------------------------------------------------------;
//     K(FNU,Z) IS COMPUTED FROM H(2,FNU,-I*Z) WHERE Z IS IN THE FIRST;
//     QUADRANT. FOURTH QUADRANT VALUES (YY <= 0.0E0) ARE COMPUTED BY;
//     CONJUGATION SINCE THE K FUNCTION IS REAL ON THE POSITIVE REAL AXIS;
//-----------------------------------------------------------------------;
      J = 2;
      DO 80 I=1,N;
//-----------------------------------------------------------------------;
//     J FLIP FLOPS BETWEEN 1 AND 2 IN J = 3 - J;
//-----------------------------------------------------------------------;
        J = 3 - J;
        FN = FNU + ((I-1) as f64);
        CALL ZUNHJ(ZNR, ZNI, FN, 0, TOL, PHIR(J), PHII(J), ARGR(J),;
     *   ARGI(J), ZETA1R(J), ZETA1I(J), ZETA2R(J), ZETA2I(J), ASUMR(J),;
     *   ASUMI(J), BSUMR(J), BSUMI(J));
        if (KODE == 1) GO TO 30;
        STR = ZBR + ZETA2R(J);
        STI = ZBI + ZETA2I(J);
        RAST = FN/ZABS(STR,STI);
        STR = STR*RAST*RAST;
        STI = -STI*RAST*RAST;
        S1R = ZETA1R(J) - STR;
        S1I = ZETA1I(J) - STI;
        GO TO 40;
   30   CONTINUE;
        S1R = ZETA1R(J) - ZETA2R(J);
        S1I = ZETA1I(J) - ZETA2I(J);
   40   CONTINUE;
//-----------------------------------------------------------------------;
//     TEST FOR UNDERFLOW AND OVERFLOW;
//-----------------------------------------------------------------------;
        RS1 = S1R;
        if ((RS1).abs() > ELIM) GO TO 70;
        if (KDFLG == 1) KFLAG = 2;
        if ((RS1).abs() < ALIM) GO TO 50;
//-----------------------------------------------------------------------;
//     REFINE  TEST AND SCALE;
//-----------------------------------------------------------------------;
        APHI = ZABS(PHIR(J),PHII(J));
        AARG = ZABS(ARGR(J),ARGI(J));
        RS1 = RS1 + DLOG(APHI) - 0.25*DLOG(AARG) - AIC;
        if ((RS1).abs() > ELIM) GO TO 70;
        if (KDFLG == 1) KFLAG = 1;
        if (RS1 < 0.0) GO TO 50;
        if (KDFLG == 1) KFLAG = 3;
   50   CONTINUE;
//-----------------------------------------------------------------------;
//     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR;
//     EXPONENT EXTREMES;
//-----------------------------------------------------------------------;
        C2R = ARGR(J)*CR2R - ARGI(J)*CR2I;
        C2I = ARGR(J)*CR2I + ARGI(J)*CR2R;
        CALL ZAIRY(C2R, C2I, 0, 2, AIR, AII, NAI, IDUM);
        CALL ZAIRY(C2R, C2I, 1, 2, DAIR, DAII, NDAI, IDUM);
        STR = DAIR*BSUMR(J) - DAII*BSUMI(J);
        STI = DAIR*BSUMI(J) + DAII*BSUMR(J);
        PTR = STR*CR2R - STI*CR2I;
        PTI = STR*CR2I + STI*CR2R;
        STR = PTR + (AIR*ASUMR(J)-AII*ASUMI(J));
        STI = PTI + (AIR*ASUMI(J)+AII*ASUMR(J));
        PTR = STR*PHIR(J) - STI*PHII(J);
        PTI = STR*PHII(J) + STI*PHIR(J);
        S2R = PTR*CSR - PTI*CSI;
        S2I = PTR*CSI + PTI*CSR;
        STR = DEXP(S1R)*CSSR(KFLAG);
        S1R = STR*DCOS(S1I);
        S1I = STR*DSIN(S1I);
        STR = S2R*S1R - S2I*S1I;
        S2I = S1R*S2I + S2R*S1I;
        S2R = STR;
        if (KFLAG != 1) GO TO 60;
        CALL ZUunderflowCHK(S2R, S2I, NW, BRY(1), TOL);
        if (NW != 0) GO TO 70;
   60   CONTINUE;
        if (YY <= 0.0) S2I = -S2I;
        CYR(KDFLG) = S2R;
        CYI(KDFLG) = S2I;
        YR(I) = S2R*CSRR(KFLAG);
        YI(I) = S2I*CSRR(KFLAG);
        STR = CSI;
        CSI = -CSR;
        CSR = STR;
        if (KDFLG == 2) GO TO 85;
        KDFLG = 2;
        GO TO 80;
   70   CONTINUE;
        if (RS1 > 0.0) GO TO 320;
//-----------------------------------------------------------------------;
//     FOR ZR < 0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW;
//-----------------------------------------------------------------------;
        if (ZR < 0.0) GO TO 320;
        KDFLG = 1;
        YR(I)=ZEROR;
        YI(I)=ZEROI;
        NZ=NZ+1;
        STR = CSI;
        CSI =-CSR;
        CSR = STR;
        if (I == 1) GO TO 80;
        if ((YR(I-1) == ZEROR)&&(YI(I-1) == ZEROI)) GO TO 80;
        YR(I-1)=ZEROR;
        YI(I-1)=ZEROI;
        NZ=NZ+1;
   80 CONTINUE;
      I = N;
   85 CONTINUE;
      RAZR = 1.0/ZABS(ZRR,ZRI);
      STR = ZRR*RAZR;
      STI = -ZRI*RAZR;
      RZR = (STR+STR)*RAZR;
      RZI = (STI+STI)*RAZR;
      CKR = FN*RZR;
      CKI = FN*RZI;
      IB = I + 1;
      if (N < IB) GO TO 180;
//-----------------------------------------------------------------------;
//     TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW. SET SEQUENCE TO ZERO;
//     ON UNDERFLOW.;
//-----------------------------------------------------------------------;
      FN = FNU + ((N-1) as f64);
      IPARD = 1;
      if (MR != 0) IPARD = 0;
      CALL ZUNHJ(ZNR, ZNI, FN, IPARD, TOL, PHIDR, PHIDI, ARGDR, ARGDI,;
     * ZET1DR, ZET1DI, ZET2DR, ZET2DI, ASUMDR, ASUMDI, BSUMDR, BSUMDI);
      if (KODE == 1) GO TO 90;
      STR = ZBR + ZET2DR;
      STI = ZBI + ZET2DI;
      RAST = FN/ZABS(STR,STI);
      STR = STR*RAST*RAST;
      STI = -STI*RAST*RAST;
      S1R = ZET1DR - STR;
      S1I = ZET1DI - STI;
      GO TO 100;
   90 CONTINUE;
      S1R = ZET1DR - ZET2DR;
      S1I = ZET1DI - ZET2DI;
  100 CONTINUE;
      RS1 = S1R;
      if ((RS1).abs() > ELIM) GO TO 105;
      if ((RS1).abs() < ALIM) GO TO 120;
//----------------------------------------------------------------------------;
//     REFINE ESTIMATE AND TEST;
//-------------------------------------------------------------------------;
      APHI = ZABS(PHIDR,PHIDI);
      RS1 = RS1+DLOG(APHI);
      if ((RS1).abs() < ELIM) GO TO 120;
  105 CONTINUE;
      if (RS1 > 0.0) GO TO 320;
//-----------------------------------------------------------------------;
//     FOR ZR < 0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW;
//-----------------------------------------------------------------------;
      if (ZR < 0.0) GO TO 320;
      NZ = N;
      DO 106 I=1,N;
        YR(I) = ZEROR;
        YI(I) = ZEROI;
  106 CONTINUE;
      RETURN;
  120 CONTINUE;
      S1R = CYR(1);
      S1I = CYI(1);
      S2R = CYR(2);
      S2I = CYI(2);
      C1R = CSRR(KFLAG);
      ASCLE = BRY(KFLAG);
      DO 130 I=IB,N;
        C2R = S2R;
        C2I = S2I;
        S2R = CKR*C2R - CKI*C2I + S1R;
        S2I = CKR*C2I + CKI*C2R + S1I;
        S1R = C2R;
        S1I = C2I;
        CKR = CKR + RZR;
        CKI = CKI + RZI;
        C2R = S2R*C1R;
        C2I = S2I*C1R;
        YR(I) = C2R;
        YI(I) = C2I;
        if (KFLAG >= 3) GO TO 130;
        STR = (C2R).abs();
        STI = (C2I).abs();
        C2M = DMAX1(STR,STI);
        if (C2M <= ASCLE) GO TO 130;
        KFLAG = KFLAG + 1;
        ASCLE = BRY(KFLAG);
        S1R = S1R*C1R;
        S1I = S1I*C1R;
        S2R = C2R;
        S2I = C2I;
        S1R = S1R*CSSR(KFLAG);
        S1I = S1I*CSSR(KFLAG);
        S2R = S2R*CSSR(KFLAG);
        S2I = S2I*CSSR(KFLAG);
        C1R = CSRR(KFLAG);
  130 CONTINUE;
  180 CONTINUE;
      if (MR == 0) RETURN;
//-----------------------------------------------------------------------;
//     ANALYTIC CONTINUATION FOR RE(Z) < 0.0;
//-----------------------------------------------------------------------;
      NZ = 0;
      FMR = (MR as f64);
      SGN = -DSIGN(PI,FMR);
//-----------------------------------------------------------------------;
//     CSPN AND CSGN ARE COEFF OF K AND I FUNCIONS RESP.;
//-----------------------------------------------------------------------;
      CSGNI = SGN;
      if (YY <= 0.0) CSGNI = -CSGNI;
      IFN = INU + N - 1;
      ANG = FNF*SGN;
      CSPNR = DCOS(ANG);
      CSPNI = DSIN(ANG);
      if (MOD(IFN,2) == 0) GO TO 190;
      CSPNR = -CSPNR;
      CSPNI = -CSPNI;
  190 CONTINUE;
//-----------------------------------------------------------------------;
//     CS=COEFF OF THE J FUNCTION TO GET THE I FUNCTION. I(FNU,Z) IS;
//     COMPUTED FROM EXP(I*FNU*FRAC_PI_2)*J(FNU,-I*Z) WHERE Z IS IN THE FIRST;
//     QUADRANT. FOURTH QUADRANT VALUES (YY <= 0.0E0) ARE COMPUTED BY;
//     CONJUGATION SINCE THE I FUNCTION IS REAL ON THE POSITIVE REAL AXIS;
//-----------------------------------------------------------------------;
      CSR = SAR*CSGNI;
      CSI = CAR*CSGNI;
      IN = MOD(IFN,4) + 1;
      C2R = CIPR(IN);
      C2I = CIPI(IN);
      STR = CSR*C2R + CSI*C2I;
      CSI = -CSR*C2I + CSI*C2R;
      CSR = STR;
      ASC = BRY(1);
      IUF = 0;
      KK = N;
      KDFLG = 1;
      IB = IB - 1;
      IC = IB - 1;
      DO 290 K=1,N;
        FN = FNU + ((KK-1) as f64);
//-----------------------------------------------------------------------;
//     LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K;
//     FUNCTION ABOVE;
//-----------------------------------------------------------------------;
        if (N > 2) GO TO 175;
  172   CONTINUE;
        PHIDR = PHIR(J);
        PHIDI = PHII(J);
        ARGDR = ARGR(J);
        ARGDI = ARGI(J);
        ZET1DR = ZETA1R(J);
        ZET1DI = ZETA1I(J);
        ZET2DR = ZETA2R(J);
        ZET2DI = ZETA2I(J);
        ASUMDR = ASUMR(J);
        ASUMDI = ASUMI(J);
        BSUMDR = BSUMR(J);
        BSUMDI = BSUMI(J);
        J = 3 - J;
        GO TO 210;
  175   CONTINUE;
        if ((KK == N)&&(IB < N)) GO TO 210;
        if ((KK == IB)||(KK == IC)) GO TO 172;
        CALL ZUNHJ(ZNR, ZNI, FN, 0, TOL, PHIDR, PHIDI, ARGDR,;
     *   ARGDI, ZET1DR, ZET1DI, ZET2DR, ZET2DI, ASUMDR,;
     *   ASUMDI, BSUMDR, BSUMDI);
  210   CONTINUE;
        if (KODE == 1) GO TO 220;
        STR = ZBR + ZET2DR;
        STI = ZBI + ZET2DI;
        RAST = FN/ZABS(STR,STI);
        STR = STR*RAST*RAST;
        STI = -STI*RAST*RAST;
        S1R = -ZET1DR + STR;
        S1I = -ZET1DI + STI;
        GO TO 230;
  220   CONTINUE;
        S1R = -ZET1DR + ZET2DR;
        S1I = -ZET1DI + ZET2DI;
  230   CONTINUE;
//-----------------------------------------------------------------------;
//     TEST FOR UNDERFLOW AND OVERFLOW;
//-----------------------------------------------------------------------;
        RS1 = S1R;
        if ((RS1).abs() > ELIM) GO TO 280;
        if (KDFLG == 1) IFLAG = 2;
        if ((RS1).abs() < ALIM) GO TO 240;
//-----------------------------------------------------------------------;
//     REFINE  TEST AND SCALE;
//-----------------------------------------------------------------------;
        APHI = ZABS(PHIDR,PHIDI);
        AARG = ZABS(ARGDR,ARGDI);
        RS1 = RS1 + DLOG(APHI) - 0.25*DLOG(AARG) - AIC;
        if ((RS1).abs() > ELIM) GO TO 280;
        if (KDFLG == 1) IFLAG = 1;
        if (RS1 < 0.0) GO TO 240;
        if (KDFLG == 1) IFLAG = 3;
  240   CONTINUE;
        CALL ZAIRY(ARGDR, ARGDI, 0, 2, AIR, AII, NAI, IDUM);
        CALL ZAIRY(ARGDR, ARGDI, 1, 2, DAIR, DAII, NDAI, IDUM);
        STR = DAIR*BSUMDR - DAII*BSUMDI;
        STI = DAIR*BSUMDI + DAII*BSUMDR;
        STR = STR + (AIR*ASUMDR-AII*ASUMDI);
        STI = STI + (AIR*ASUMDI+AII*ASUMDR);
        PTR = STR*PHIDR - STI*PHIDI;
        PTI = STR*PHIDI + STI*PHIDR;
        S2R = PTR*CSR - PTI*CSI;
        S2I = PTR*CSI + PTI*CSR;
        STR = DEXP(S1R)*CSSR(IFLAG);
        S1R = STR*DCOS(S1I);
        S1I = STR*DSIN(S1I);
        STR = S2R*S1R - S2I*S1I;
        S2I = S2R*S1I + S2I*S1R;
        S2R = STR;
        if (IFLAG != 1) GO TO 250;
        CALL ZUunderflowCHK(S2R, S2I, NW, BRY(1), TOL);
        if (NW == 0) GO TO 250;
        S2R = ZEROR;
        S2I = ZEROI;
  250   CONTINUE;
        if (YY <= 0.0) S2I = -S2I;
        CYR(KDFLG) = S2R;
        CYI(KDFLG) = S2I;
        C2R = S2R;
        C2I = S2I;
        S2R = S2R*CSRR(IFLAG);
        S2I = S2I*CSRR(IFLAG);
//-----------------------------------------------------------------------;
//     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N;
//-----------------------------------------------------------------------;
        S1R = YR(KK);
        S1I = YI(KK);
        if (KODE == 1) GO TO 270;
        CALL ZS1S2(ZRR, ZRI, S1R, S1I, S2R, S2I, NW, ASC, ALIM, IUF);
        NZ = NZ + NW;
  270   CONTINUE;
        YR(KK) = S1R*CSPNR - S1I*CSPNI + S2R;
        YI(KK) = S1R*CSPNI + S1I*CSPNR + S2I;
        KK = KK - 1;
        CSPNR = -CSPNR;
        CSPNI = -CSPNI;
        STR = CSI;
        CSI = -CSR;
        CSR = STR;
        if (C2R != 0.0 || C2I != 0.0) GO TO 255;
        KDFLG = 1;
        GO TO 290;
  255   CONTINUE;
        if (KDFLG == 2) GO TO 295;
        KDFLG = 2;
        GO TO 290;
  280   CONTINUE;
        if (RS1 > 0.0) GO TO 320;
        S2R = ZEROR;
        S2I = ZEROI;
        GO TO 250;
  290 CONTINUE;
      K = N;
  295 CONTINUE;
      IL = N - K;
      if (IL == 0) RETURN;
//-----------------------------------------------------------------------;
//     RECUR BACKWARD FOR REMAINDER OF I SEQUENCE AND ADD IN THE;
//     K FUNCTIONS, SCALING THE I SEQUENCE DURING RECURRENCE TO KEEP;
//     INTERMEDIATE ARITHMETIC ON SCALE NEAR EXPONENT EXTREMES.;
//-----------------------------------------------------------------------;
      S1R = CYR(1);
      S1I = CYI(1);
      S2R = CYR(2);
      S2I = CYI(2);
      CSR = CSRR(IFLAG);
      ASCLE = BRY(IFLAG);
      FN = DBLE(FLOAT(INU+IL));
      DO 310 I=1,IL;
        C2R = S2R;
        C2I = S2I;
        S2R = S1R + (FN+FNF)*(RZR*C2R-RZI*C2I);
        S2I = S1I + (FN+FNF)*(RZR*C2I+RZI*C2R);
        S1R = C2R;
        S1I = C2I;
        FN = FN - 1.0;
        C2R = S2R*CSR;
        C2I = S2I*CSR;
        CKR = C2R;
        CKI = C2I;
        C1R = YR(KK);
        C1I = YI(KK);
        if (KODE == 1) GO TO 300;
        CALL ZS1S2(ZRR, ZRI, C1R, C1I, C2R, C2I, NW, ASC, ALIM, IUF);
        NZ = NZ + NW;
  300   CONTINUE;
        YR(KK) = C1R*CSPNR - C1I*CSPNI + C2R;
        YI(KK) = C1R*CSPNI + C1I*CSPNR + C2I;
        KK = KK - 1;
        CSPNR = -CSPNR;
        CSPNI = -CSPNI;
        if (IFLAG >= 3) GO TO 310;
        C2R = (CKR).abs();
        C2I = (CKI).abs();
        C2M = DMAX1(C2R,C2I);
        if (C2M <= ASCLE) GO TO 310;
        IFLAG = IFLAG + 1;
        ASCLE = BRY(IFLAG);
        S1R = S1R*CSR;
        S1I = S1I*CSR;
        S2R = CKR;
        S2I = CKI;
        S1R = S1R*CSSR(IFLAG);
        S1I = S1I*CSSR(IFLAG);
        S2R = S2R*CSSR(IFLAG);
        S2I = S2I*CSSR(IFLAG);
        CSR = CSRR(IFLAG);
  310 CONTINUE;
      RETURN;
  320 CONTINUE;
      NZ = -1;
      RETURN;
      END;
        }

*/
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
        let mut cy = c_zeros(2);
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
        let (mut overflow_state, mut ASCLE, mut CSCLR) = if cy[0].abs() <= MACHINE_CONSTANTS.bry[0]
        {
            (
                Overflow::NearUnder,
                MACHINE_CONSTANTS.bry[0],
                1.0 / MACHINE_CONSTANTS.abs_error_tolerance,
            )
        } else if cy[0].abs() >= MACHINE_CONSTANTS.bry[1] {
            (
                Overflow::NearOver,
                MACHINE_CONSTANTS.bry[2],
                MACHINE_CONSTANTS.abs_error_tolerance,
            )
        } else {
            (Overflow::None, MACHINE_CONSTANTS.bry[1], 1.0)
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
            ASCLE = MACHINE_CONSTANTS.bry[overflow_state];
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
            ASCLE = MACHINE_CONSTANTS.bry[overflow_state];
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
            let NUF = zuoik(z, order, scaling, IKType::I, ND, y)?;
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
                && will_z_underflow(
                    s2,
                    MACHINE_CONSTANTS.bry[0],
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
    let mut ASCLE = MACHINE_CONSTANTS.bry[overflow_state];
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
        ASCLE = MACHINE_CONSTANTS.bry[overflow_state];
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
    let CIP = [
        Complex64::new(1.0, 0.0),
        Complex64::new(0.0, 1.0),
        Complex64::new(-1.0, 0.0),
        Complex64::new(0.0, -1.0),
    ];
    const AIC: f64 = 1.265512123484645396;
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
    let (_, _, zeta1, zeta2, _, _) = zunhj(
        zn,
        modified_order,
        true,
        MACHINE_CONSTANTS.abs_error_tolerance,
    );

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
            let NUF = zuoik(z, order, scaling, IKType::I, ND, y)?;
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
            let (phi, arg, zeta1, zeta2, asum, bsum) = zunhj(
                zn,
                modified_order,
                false,
                MACHINE_CONSTANTS.abs_error_tolerance,
            );
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
            //note that ZAIRY calls in fortran code ignore IERR (using IDUM)
            let a_airy = match ZAIRY(arg, false, Scaling::Scaled) {
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
            let d_airy = match ZAIRY(arg, true, Scaling::Scaled) {
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

            let mut s2 = phi * (d_airy * bsum + a_airy * asum);
            let s1 = MACHINE_CONSTANTS.scaling_factors[overflow_state] * s1.exp();
            s2 *= s1;
            if overflow_state == Overflow::NearUnder
                && will_z_underflow(
                    s2,
                    MACHINE_CONSTANTS.bry[0],
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
        let mut ASCLE = MACHINE_CONSTANTS.bry[overflow_state];
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
            ASCLE = MACHINE_CONSTANTS.bry[overflow_state];
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
