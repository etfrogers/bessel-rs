#![allow(non_snake_case)]
use super::{
    BesselError, BesselResult, IKType, MachineConsts, Scaling, c_one, c_zero, c_zeros, gamma_ln,
    i_power_series,
    machine::i1mach,
    overflow_checks::{zunik, zuoik},
    utils::will_z_underflow,
};
use crate::amos::{
    BesselError::*, max_abs_component, overflow_checks::zunhj, z_asymptotic_i::z_asymptotic_i,
};
use num::{
    Zero,
    complex::{Complex64, ComplexFloat},
    pow::Pow,
};
use std::{
    cmp::min,
    f64::{
        self,
        consts::{FRAC_PI_2, PI},
    },
};

const TWO_THIRDS: f64 = 6.66666666666666666e-01;

/*
fn ZBESH(ZR, ZI, FNU, KODE, M, N, CYR, CYI, NZ, IERR)
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
//
//     COMPLEX CY,Z,ZN,ZT,CSGN
//       EXTERNAL ZABS
//       DOUBLE PRECISION AA, ALIM, ALN, ARG, AZ, CYI, CYR, DIG, ELIM,
//      * FMM, FN, FNU, FNUL, FRAC_PI_2, RFRAC_PI_2, RL, R1M5, SGN, STR, TOL, UFL, ZI,
//      * ZNI, ZNR, ZR, ZTI, d1mach, ZABS, BB, ASCLE, RTOL, ATOL, STI,
//      * CSGNR, CSGNI
//       INTEGER I, IERR, INU, INUH, IR, K, KODE, K1, K2, M,
//      * MM, MR, N, NN, NUF, NW, NZ, i1mach
//       DIMENSION CYR(N), CYI(N)
//
//      DATA FRAC_PI_2 /1.57079632679489662/
//
// ***FIRST EXECUTABLE STATEMENT  ZBESH
      IERR = 0
      NZ=0
      if (ZR == 0.0 && ZI == 0.0) IERR=1
      if (FNU < 0.0) IERR=1
      if (M < 1 || M > 2) IERR=1
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
      FN = FNU + ((NN-1) as f64)
      MM = 3 - M - M
      FMM = (MM as f64)
      ZNR = FMM*ZI
      ZNI = -FMM*ZR
//-----------------------------------------------------------------------
//     TEST FOR PROPER RANGE
//-----------------------------------------------------------------------
      AZ = ZABS(ZR,ZI)
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
      UFL = d1mach(1)*1.0e+3
      if (AZ < UFL) GO TO 230
      if (FNU > FNUL) GO TO 90
      if (FN <= 1.0) GO TO 70
      if (FN > 2.0) GO TO 60
      if (AZ > TOL) GO TO 70
      ARG = 0.5*AZ
      ALN = -FN*DLOG(ARG)
      if (ALN > ELIM) GO TO 230
      GO TO 70
   60 CONTINUE
      CALL ZUOIK(ZNR, ZNI, FNU, KODE, 2, NN, CYR, CYI, NUF, TOL, ELIM,
     * ALIM)
      if (NUF < 0) GO TO 230
      NZ = NZ + NUF
      NN = NN - NUF
//-----------------------------------------------------------------------
//     HERE NN=N OR NN=0 SINCE NUF=0,NN, OR -1 ON RETURN FROM CUOIK
//     if NUF=NN, THEN CY(I)=CZERO FOR ALL I
//-----------------------------------------------------------------------
      if (NN == 0) GO TO 140
   70 CONTINUE
      if ((ZNR < 0.0) || (ZNR == 0.0 && ZNI < 0.0 &&
     * M == 2)) GO TO 80
//-----------------------------------------------------------------------
//     RIGHT HALF PLANE COMPUTATION, XN >= 0. && (XN != 0. ||
//     YN >= 0. || M=1)
//-----------------------------------------------------------------------
      CALL ZBKNU(ZNR, ZNI, FNU, KODE, NN, CYR, CYI, NZ, TOL, ELIM, ALIM)
      GO TO 110
//-----------------------------------------------------------------------
//     LEFT HALF PLANE COMPUTATION
//-----------------------------------------------------------------------
   80 CONTINUE
      MR = -MM
      CALL ZACON(ZNR, ZNI, FNU, KODE, MR, NN, CYR, CYI, NW, RL, FNUL,
     * TOL, ELIM, ALIM)
      if (NW < 0) GO TO 240
      NZ=NW
      GO TO 110
   90 CONTINUE
//-----------------------------------------------------------------------
//     UNIFORM ASYMPTOTIC EXPANSIONS FOR FNU > FNUL
//-----------------------------------------------------------------------
      MR = 0
      if ((ZNR >= 0.0) && (ZNR != 0.0 || ZNI >= 0.0 ||
     * M != 2)) GO TO 100
      MR = -MM
      if (ZNR != 0.0 || ZNI >= 0.0) GO TO 100
      ZNR = -ZNR
      ZNI = -ZNI
  100 CONTINUE
      CALL ZBUNK(ZNR, ZNI, FNU, KODE, MR, NN, CYR, CYI, NW, TOL, ELIM,
     * ALIM)
      if (NW < 0) GO TO 240
      NZ = NZ + NW
  110 CONTINUE
//-----------------------------------------------------------------------
//     H(M,FNU,Z) = -FMM*(I/FRAC_PI_2)*(ZT**FNU)*K(FNU,-Z*ZT)
//
//     ZT=EXP(-FMM*FRAC_PI_2*I) = CMPLX(0.0,-FMM), FMM=3-2*M, M=1,2
//-----------------------------------------------------------------------
      SGN = DSIGN(FRAC_PI_2,-FMM)
//-----------------------------------------------------------------------
//     CALCULATE EXP(FNU*FRAC_PI_2*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
//     WHEN FNU IS LARGE
//-----------------------------------------------------------------------
      INU = INT(SNGL(FNU))
      INUH = INU/2
      IR = INU - 2*INUH
      ARG = (FNU-((INU-IR) as f64))*SGN
      RFRAC_PI_2 = 1.0/SGN
//     ZNI = RFRAC_PI_2*DCOS(ARG)
//     ZNR = -RFRAC_PI_2*DSIN(ARG)
      CSGNI = RFRAC_PI_2*DCOS(ARG)
      CSGNR = -RFRAC_PI_2*DSIN(ARG)
      if (MOD(INUH,2) == 0) GO TO 120
//     ZNR = -ZNR
//     ZNI = -ZNI
      CSGNR = -CSGNR
      CSGNI = -CSGNI
  120 CONTINUE
      ZTI = -FMM
      RTOL = 1.0/TOL
      ASCLE = UFL*RTOL
      DO 130 I=1,NN
//       STR = CYR(I)*ZNR - CYI(I)*ZNI
//       CYI(I) = CYR(I)*ZNI + CYI(I)*ZNR
//       CYR(I) = STR
//       STR = -ZNI*ZTI
//       ZNI = ZNR*ZTI
//       ZNR = STR
        AA = CYR(I)
        BB = CYI(I)
        ATOL = 1.0
        if (DMAX1((AA).abs(),(BB).abs()) > ASCLE) GO TO 135
          AA = AA*RTOL
          BB = BB*RTOL
          ATOL = TOL
  135 CONTINUE
      STR = AA*CSGNR - BB*CSGNI
      STI = AA*CSGNI + BB*CSGNR
      CYR(I) = STR*ATOL
      CYI(I) = STI*ATOL
      STR = -CSGNI*ZTI
      CSGNI = CSGNR*ZTI
      CSGNR = STR
  130 CONTINUE
      RETURN
  140 CONTINUE
      if (ZNR < 0.0) GO TO 230
      RETURN
  230 CONTINUE
      NZ=0
      IERR=2
      RETURN
  240 CONTINUE
      if(NW == (-1)) GO TO 230
      NZ=0
      IERR=5
      RETURN
  260 CONTINUE
      NZ=0
      IERR=4
      RETURN
      END
fn ZBESI(ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, IERR)
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
//
// ***ROUTINES CALLED  ZBINU,ZABS,i1mach,d1mach
// ***END PROLOGUE  ZBESI
//     COMPLEX CONE,CSGN,CW,CY,CZERO,Z,ZN
      EXTERNAL ZABS
      DOUBLE PRECISION AA, ALIM, ARG, CONEI, CONER, CSGNI, CSGNR, CYI,
     * CYR, DIG, ELIM, FNU, FNUL, PI, RL, R1M5, STR, TOL, ZI, ZNI, ZNR,
     * ZR, d1mach, AZ, BB, FN, ZABS, ASCLE, RTOL, ATOL, STI
      INTEGER I, IERR, INU, K, KODE, K1,K2,N,NZ,NN, i1mach
      DIMENSION CYR(N), CYI(N)
      DATA PI /3.14159265358979324/
      DATA CONER, CONEI /1.0,0.0/
//
// ***FIRST EXECUTABLE STATEMENT  ZBESI
      IERR = 0
      NZ=0
      if (FNU < 0.0) IERR=1
      if (KODE < 1 || KODE > 2) IERR=1
      if (N < 1) IERR=1
      if (IERR != 0) RETURN
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
      RL = 1.2*DIG + 3.0
      FNUL = 10.0 + 6.0*(DIG-3.0)
//-----------------------------------------------------------------------------
//     TEST FOR PROPER RANGE
//-----------------------------------------------------------------------
      AZ = ZABS(ZR,ZI)
      FN = FNU+((N-1) as f64)
      AA = 0.5/TOL
      BB=DBLE(FLOAT(i1mach(9)))*0.5
      AA = DMIN1(AA,BB)
      if (AZ > AA) {return Err(LossOfSignificance);}
      if (FN > AA) {return Err(LossOfSignificance);}
      AA = DSQRT(AA)
      if (AZ > AA) IERR=3
      if (FN > AA) IERR=3
      ZNR = ZR
      ZNI = ZI
      CSGNR = CONER
      CSGNI = CONEI
      if (ZR >= 0.0) GO TO 40
      ZNR = -ZR
      ZNI = -ZI
//-----------------------------------------------------------------------
//     CALCULATE CSGN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
//     WHEN FNU IS LARGE
//-----------------------------------------------------------------------
      INU = INT(SNGL(FNU))
      ARG = (FNU-(INU as f64))*PI
      if (ZI < 0.0) ARG = -ARG
      CSGNR = DCOS(ARG)
      CSGNI = DSIN(ARG)
      if (MOD(INU,2) == 0) GO TO 40
      CSGNR = -CSGNR
      CSGNI = -CSGNI
   40 CONTINUE
      CALL ZBINU(ZNR, ZNI, FNU, KODE, N, CYR, CYI, NZ, RL, FNUL, TOL,
     * ELIM, ALIM)
      if (NZ < 0) GO TO 120
      if (ZR >= 0.0) RETURN
//-----------------------------------------------------------------------
//     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE
//-----------------------------------------------------------------------
      NN = N - NZ
      if (NN == 0) RETURN
      RTOL = 1.0/TOL
      ASCLE = d1mach(1)*RTOL*1.0e+3
      DO 50 I=1,NN
//       STR = CYR(I)*CSGNR - CYI(I)*CSGNI
//       CYI(I) = CYR(I)*CSGNI + CYI(I)*CSGNR
//       CYR(I) = STR
        AA = CYR(I)
        BB = CYI(I)
        ATOL = 1.0
        if (DMAX1((AA).abs(),(BB).abs()) > ASCLE) GO TO 55
          AA = AA*RTOL
          BB = BB*RTOL
          ATOL = TOL
   55   CONTINUE
        STR = AA*CSGNR - BB*CSGNI
        STI = AA*CSGNI + BB*CSGNR
        CYR(I) = STR*ATOL
        CYI(I) = STI*ATOL
        CSGNR = -CSGNR
        CSGNI = -CSGNI
   50 CONTINUE
      RETURN
  120 CONTINUE
      if(NZ == (-2)) GO TO 130
      NZ = 0
      IERR=2
      RETURN
  130 CONTINUE
      NZ=0
      IERR=5
      RETURN
  260 CONTINUE
      NZ=0
      IERR=4
      RETURN
      END
*/
pub fn zbesj(
    z: Complex64, //ZR, ZI,
    order: f64,   //FNU,
    KODE: Scaling,
    N: usize,
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
    // ***ROUTINES CALLED  ZBINU,ZABS,i1mach,d1mach
    // ***END PROLOGUE  ZBESJ
    //
    //     COMPLEX CI,CSGN,CY,Z,ZN
    //       EXTERNAL ZABS
    //       DOUBLE PRECISION AA, ALIM, ARG, CII, CSGNI, CSGNR, CYI, CYR, DIG,
    //      * ELIM, FNU, FNUL, FRAC_PI_2, RL, R1M5, STR, TOL, ZI, ZNI, ZNR, ZR,
    //      * d1mach, BB, FN, AZ, ZABS, ASCLE, RTOL, ATOL, STI
    //       INTEGER I, IERR, INU, INUH, IR, K, KODE, K1, K2, N, NL, NZ, i1mach
    //       DIMENSION CYR(N), CYI(N)
    //      DATA FRAC_PI_2 /1.57079632679489662/
    //
    // ***FIRST EXECUTABLE STATEMENT  ZBESJ
    // IERR = 0
    let mut significance_loss = false;
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
    let machine_consts = MachineConsts::new();
    //-----------------------------------------------------------------------
    //     TEST FOR PROPER RANGE
    //-----------------------------------------------------------------------
    let az = z.abs(); //ZABS(ZR, ZI);
    let FN = order + ((N - 1) as f64);
    let mut AA = 0.5 / machine_consts.tol;
    let mut bb = (i1mach(9) as f64) * 0.5;
    AA = AA.min(bb);
    if az > AA {
        return Err(LossOfSignificance);
    }
    if FN > AA {
        return Err(LossOfSignificance);
    }
    AA = AA.sqrt();
    if (az > AA) || (FN > AA) {
        significance_loss = true;
    }
    //-----------------------------------------------------------------------
    //     CALCULATE CSGN=EXP(FNU*FRAC_PI_2*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
    //     WHEN FNU IS LARGE
    //-----------------------------------------------------------------------
    let mut CII = 1.0;
    let INU = order as i64;
    let INUH = INU / 2;
    let IR = INU - 2 * INUH;
    let ARG = (order - ((INU - IR) as f64)) * FRAC_PI_2;
    let mut CSGNR = ARG.cos();
    let mut CSGNI = ARG.sin();
    if (INUH % 2) != 0 {
        CSGNR = -CSGNR;
        CSGNI = -CSGNI;
    }
    //-----------------------------------------------------------------------
    //     ZN IS IN THE RIGHT HALF PLANE
    //-----------------------------------------------------------------------
    let mut ZNR = z.im;
    let mut ZNI = -z.re;
    if !(z.im >= 0.0) {
        ZNR = -ZNR;
        ZNI = -ZNI;
        CSGNI = -CSGNI;
        CII = -CII;
    }
    let (mut cy, NZ) = ZBINU(Complex64::new(ZNR, ZNI), order, KODE, N, &machine_consts)?;
    let NL = N - NZ;
    for i in 0..NL {
        AA = cy[i].re;
        bb = cy[i].im;
        let mut ATOL = 1.0;
        if !((AA.abs().max(bb.abs())) > machine_consts.ascle) {
            AA = AA * machine_consts.rtol;
            bb = bb * machine_consts.rtol;
            ATOL = machine_consts.tol;
        }
        let mut STR = AA * CSGNR - bb * CSGNI;
        let STI = AA * CSGNI + bb * CSGNR;
        cy[i] = Complex64::new(STR * ATOL, STI * ATOL);
        STR = -CSGNI * CII;
        CSGNI = CSGNR * CII;
        CSGNR = STR;
    }
    if significance_loss {
        Err(PartialLossOfSignificance { y: cy, nz: NZ })
    } else {
        Ok((cy, NZ))
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

fn ZAIRY(//ZR, ZI, ID, KODE, AIR, AII, NZ, IERR){
    z:Complex64, return_derivative: bool, KODE: Scaling) -> BesselResult<(Complex64, usize)> {
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
//     COMPLEX AI,CONE,CSQ,CY,S1,S2,TRM1,TRM2,Z,ZTA,Z3
    //   EXTERNAL ZABS
    //   DOUBLE PRECISION AA, AD, AII, AIR, AK, ALIM, ATRM, AZ, AZ3, BK,
    //  * CC, CK, COEF, CONEI, CONER, CSQI, CSQR, CYI, CYR, C1, C2, DIG,
    //  * DK, D1, D2, ELIM, FID, FNU, PTR, RL, R1M5, SFAC, STI, STR,
    //  * S1I, S1R, S2I, S2R, TOL, TRM1I, TRM1R, TRM2I, TRM2R, TWO_THIRDS, ZEROI,
    //  * ZEROR, ZI, ZR, ZTAI, ZTAR, Z3I, Z3R, d1mach, ZABS, ALAZ, BB
    //   INTEGER ID, IERR, IFLAG, K, KODE, K1, K2, MR, NN, NZ, i1mach
    //   DIMENSION CYR(1), CYI(1)
    //   DATA TWO_THIRDS, C1, C2, COEF /6.66666666666666667e-01,
     const C1:f64= 3.55028053887817240e-01;
     const C2:f64 = 2.58819403792806799e-01;
     const COEFF:f64 = 1.83776298473930683e-01;
    //   DATA ZEROR, ZEROI, CONER, CONEI /0.0,0.0,1.0,0.0/
// ***FIRST EXECUTABLE STATEMENT  ZAIRY
    //   IERR = 0;
    //   let NZ=0;
    //   if (ID < 0 || ID > 1) IERR=1;
    //   if (KODE < 1 || KODE > 2) IERR=1;
    //   if (IERR != 0) RETURN;
      let AZ = z.abs();//ZABS(ZR,ZI);
    //   TOL = DMAX1(d1mach(4),1.0e-18);
    let machine_consts = MachineConsts::new();
      let FID = if return_derivative{1.0} else{0.0};//(ID as f64);
      let mut significance_loss = false;
      if AZ <= 1.0 {//GO TO 70;
//-----------------------------------------------------------------------;
//     POWER SERIES FOR CABS(Z) <= 1.;
//-----------------------------------------------------------------------;
let mut s1 = c_one();
let mut s2 = c_one();
    //   S1R = CONER;
    //   S1I = CONEI;
    //   S2R = CONER;
    //   S2I = CONEI;
      if AZ < machine_consts.tol {//GO TO 170;
        // 170 CONTINUE;
        // AA = 1.0e+3*d1mach(1);
        // S1R = ZEROR;
        // S1I = ZEROI;
        s1 = c_zero();
        return if return_derivative {//GO TO 190;
            // 190 CONTINUE;
            let mut ai = Complex64::new(-C2, 0.0);
        // AIR = -C2;
        // AII = 0.0;
        // AA = DSQRT(AA);

        if AZ > machine_consts.arm.sqrt() {//GO TO 200;
            s1 = z.pow(2.0)/2.0;
        // S1R = 0.5*(ZR*ZR-ZI*ZI);
        // S1I = ZR*ZI;
        }
    // 200 CONTINUE;
    ai+= C1 * s1;
        // AIR = AIR + C1*S1R;
        // AII = AII + C1*S1I;
        // RETURN;
        Ok((ai, 0))
    }else{
        if AZ > machine_consts.arm {//GO TO 180;
            s1 = C2*z;
        // S1R = C2*ZR;
        // S1I = C2*ZI;
        }
    // 180 CONTINUE;
    let ai = C1-s1;
        // AIR = C1 - S1R;
        // AII = -S1I;
        // RETURN;
        Ok((ai, 0))

    }
      }
      let AA = AZ*AZ;
      if AA >= machine_consts.tol/AZ {//GO TO 40;
        let mut term1 = c_one();
        let mut term2 = c_one();
    //   TRM1R = CONER;
    //   TRM1I = CONEI;
    //   TRM2R = CONER;
    //   TRM2I = CONEI;
      let mut a_term = 1.0;
      let z3 = z.pow(3.0);
    //   STR = ZR*ZR - ZI*ZI;
    //   STI = ZR*ZI + ZI*ZR;
    //   Z3R = STR*ZR - STI*ZI;
    //   Z3I = STR*ZI + STI*ZR;
      let AZ3 = AZ*AA;
      let (mut AK, mut BK,CK, DK) = (2.0+FID, 3.0-2.0*FID, 4.0-FID,3.0-2.0*FID);
    //   AK = 2.0 + FID;
    //   BK = 3.0 - FID - FID;
    //   CK = 4.0 - FID;
    //   DK = 3.0 + FID + FID;
      let mut D1:f64 = AK*DK;
      let mut D2 = BK*CK;
      let mut AD = D1.min(D2);//DMIN1(D1,D2);
      AK = 24.0 + 9.0*FID;
      BK = 30.0 - 9.0*FID;
    //   DO 30 K=1,25;
      for _ in 0..25{
        // STR = (TRM1R*Z3R-TRM1I*Z3I)/D1;
        // TRM1I = (TRM1R*Z3I+TRM1I*Z3R)/D1;
        // TRM1R = STR;
        term1 = term1*z3/D1;
        s1 += term1;
        // S1R = S1R + TRM1R;
        // S1I = S1I + TRM1I;
        term2 = term2*z3 /D2;
        // STR = (TRM2R*Z3R-TRM2I*Z3I)/D2;
        // TRM2I = (TRM2R*Z3I+TRM2I*Z3R)/D2;
        // TRM2R = STR;
        s2 += term2;
        // S2R = S2R + TRM2R;
        // S2I = S2I + TRM2I;
        a_term = a_term*AZ3/AD;
        D1 += AK;
        D2 += BK;
        AD = D1.min(D2);
        // AD = DMIN1(D1,D2);
        if a_term < machine_consts.tol*AD {break;}//GO TO 40;
        AK = AK + 18.0;
        BK = BK + 18.0;
//    30 CONTINUE;
      }
    }
      return if return_derivative {//GO TO 50;
        // 50 CONTINUE;
        let mut ai = -s2*C2;
        // AIR = -S2R*C2;
        // AII = -S2I*C2;
        if AZ > machine_consts.tol {//GO TO 60;
        // STR = ZR*S1R - ZI*S1I;
        // STI = ZR*S1I + ZI*S1R;
        let CC = C1/(1.0+FID);
        ai += CC * z.pow(2.0)*s1;
        // AIR = AIR + CC*(STR*ZR-STI*ZI);
        // AII = AII + CC*(STR*ZI+STI*ZR);
        Ok((ai, 0))
        }else{
    //  60 CONTINUE;
        if KODE == Scaling::Scaled
        {
            ai *= (TWO_THIRDS*z*z.sqrt()).exp();
        // CALL ZSQRT(ZR, ZI, STR, STI);
        // ZTAR = TWO_THIRDS*(ZR*STR-ZI*STI);
        // ZTAI = TWO_THIRDS*(ZR*STI+ZI*STR);
        // CALL ZEXP(ZTAR, ZTAI, STR, STI);
        // PTR = STR*AIR - STI*AII;
        // AII = STR*AII + STI*AIR;
        // AIR = PTR;
        }
        Ok((ai, 0))

    }
        // RETURN;

      }else{
//    40 CONTINUE;
        let mut ai = s1*C1 - C2*z*s2;
//      AIR = S1R*C1 - C2*(ZR*S2R-ZI*S2I);
//       AII = S1I*C1 - C2*(ZR*S2I+ZI*S2R);
      if KODE == Scaling::Scaled {//RETURN;

        ai *= (TWO_THIRDS*z*z.sqrt()).exp();
    //   CALL ZSQRT(ZR, ZI, STR, STI);
    //   ZTAR = TWO_THIRDS*(ZR*STR-ZI*STI);
    //   ZTAI = TWO_THIRDS*(ZR*STI+ZI*STR);
    //   CALL ZEXP(ZTAR, ZTAI, STR, STI);
    //   PTR = AIR*STR - AII*STI;
    //   AII = AIR*STI + AII*STR;
    //   AIR = PTR;
    //   RETURN;
      }
      Ok((ai, 0))
      }
    }else{
//-----------------------------------------------------------------------;
//     CASE FOR CABS(Z) > 1.0;
//-----------------------------------------------------------------------;
//    70 CONTINUE;
      let FNU = (1.0+FID)/3.0;
//-----------------------------------------------------------------------;
//     SET PARAMETERS RELATED TO MACHINE CONSTANTS.;
//     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0e-18.;
//     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.;
//     EXP(-ELIM) < EXP(-ALIM)=EXP(-ELIM)/TOL    AND;
//     EXP(ELIM) > EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR;
//     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.;
//     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.;
//     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).;
//-----------------------------------------------------------------------;
    //   K1 = i1mach(15);
    //   K2 = i1mach(16);
    //   R1M5 = d1mach(5);
    //   K = MIN0(K1.abs(),K2.abs());
    //   ELIM = 2.303*(K as f64)*R1M5-3.0);
    //   K1 = i1mach(14) - 1;
    //   AA = R1M5*(K1 as f64);
    //   DIG = DMIN1(AA,18.0);
    //   AA = AA*2.303;
    //   ALIM = ELIM + DMAX1(-AA,-41.45);
    //   RL = 1.2*DIG + 3.0;
    //   ALAZ = DLOG(AZ);
    let ALAZ = AZ.ln();
//--------------------------------------------------------------------------;
//     TEST FOR PROPER RANGE;
//-----------------------------------------------------------------------;
      let mut AA= 0.5/machine_consts.tol;
    //   BB=DBLE(FLOAT(i1mach(9)))*0.5;
     AA=AA.min(i32::MAX as f64/2.0);
      AA=AA.pow(TWO_THIRDS);
      if AZ > AA {
        return Err(LossOfSignificance);
    };
    //   AA=DSQRT(AA);
        AA = AA.sqrt();
        significance_loss = AZ>AA;
    //   if (AZ > AA) IERR=3;
    let csq = z.sqrt();
    //   CALL ZSQRT(ZR, ZI, CSQR, CSQI);
    let mut zta = TWO_THIRDS*z*csq;
    //   ZTAR = TWO_THIRDS*(ZR*CSQR-ZI*CSQI);
    //   ZTAI = TWO_THIRDS*(ZR*CSQI+ZI*CSQR);
//-----------------------------------------------------------------------;
//     RE(ZTA) <= 0 WHEN RE(Z) < 0, ESPECIALLY WHEN IM(Z) IS SMALL;
//-----------------------------------------------------------------------;
      let mut IFLAG = 0;
      let mut SFAC = 1.0;
    //   let AK = zta.im;//ZTAI;
      if z.re < 0.0 {//GO TO 80;
    //   BK = ZTAR;
    //   CK = -(BK).abs();
    //   ZTAR = CK;
    //   ZTAI = AK;
        zta.re = -zta.re.abs();
      }
//    80 CONTINUE;
    if z.im ==0.0 && z.re<=0.0{
    //   if (ZI != 0.0) GO TO 90;
    //   if (ZR > 0.0) GO TO 90;
    //   ZTAR = 0.0;
    //   ZTAI = AK;
    zta.re = 0.0;
    }
//    90 CONTINUE;
      let mut AA = zta.re;//ZTAR;
      let (cy, NZ) = if !(AA >= 0.0 && z.re > 0.0) {//GO TO 110;
    //   if (KODE == 2) GO TO 100;
//-----------------------------------------------------------------------;
//     OVERFLOW TEST;
//-----------------------------------------------------------------------;
        if KODE == Scaling::Unscaled && AA<= -machine_consts.alim{
    //   if (AA > (-ALIM)) GO TO 100;
      AA = -AA + 0.25*ALAZ;
      IFLAG = 1;
      SFAC = machine_consts.tol;
      if AA > machine_consts.elim {return Err(Overflow);}//GO TO 270;
        }
//   100 CONTINUE;
//-----------------------------------------------------------------------;
//     CBKNU AND CACON RETURN EXP(ZTA)*K(FNU,ZTA) ON KODE=2;
//-----------------------------------------------------------------------;
      let MR = if z.im < 0.0 { -1} else{1};
      ZACAI(zta, FNU, KODE, MR, 1, &machine_consts)?
    //    ZACAI(ZTAR, ZTAI, FNU, KODE, MR, 1, CYR, CYI, NN, RL, TOL,
    //   ELIM, ALIM)?;
    //   if (NN < 0) GO TO 280;
    //   NZ = NZ + NN;
    //   GO TO 130;
    }else{
//   110 CONTINUE;
    //   if (KODE == 2) GO TO 120;
//-----------------------------------------------------------------------;
//     UNDERFLOW TEST;
//-----------------------------------------------------------------------;
    //   if (AA < ALIM) GO TO 120;
    if KODE == Scaling::Unscaled && AA>machine_consts.alim{
      AA = -AA - 0.25*ALAZ;
      IFLAG = 2;
      SFAC = 1.0/machine_consts.tol;
      if AA < -machine_consts.elim {return Ok((c_zero(), 1));} //GO TO 210;
    }
//   120 CONTINUE;
    ZBKNU(z, FNU, KODE, 1, &machine_consts)?
    //   CALL ZBKNU(ZTAR, ZTAI, FNU, KODE, 1, CYR, CYI, NZ, TOL, ELIM,
    //  * ALIM);
      };
//   130 CONTINUE;
      let mut s1 = cy[0]*COEFF;
    //   S1R = CYR(1)*COEF;
    //   S1I = CYI(1)*COEF;
    let retval = if IFLAG == 0 {//GO TO 150;
      let ai = if return_derivative //GO TO 140;
      {
        // 140 CONTINUE;
    //   AIR = -(ZR*S1R-ZI*S1I);
    //   AII = -(ZR*S1I+ZI*S1R);
        -z*s1
      }else{
        csq*s1
    //   AIR = CSQR*S1R - CSQI*S1I;
    //   AII = CSQR*S1I + CSQI*S1R;
      };
 (ai, NZ)
    //   RETURN;

    //   RETURN;
    }else{

//   150 CONTINUE;
        s1 *=SFAC;
    //   S1R = S1R*SFAC;
    //   S1I = S1I*SFAC;
    s1 *= if return_derivative{
        // 160 CONTINUE;
        -z
        // STR = -(S1R*ZR-S1I*ZI);
        // S1I = -(S1R*ZI+S1I*ZR);
        // S1R = STR;
        // AIR = S1R/SFAC;
        // AII = S1I/SFAC;
        // RETURN;

    }else{
    //   if (ID == 1) GO TO 160;
    //   STR = S1R*CSQR - S1I*CSQI;
    //   S1I = S1R*CSQI + S1I*CSQR;
    //   S1R = STR;
    //   AIR = S1R/SFAC;
    //   AII = S1I/SFAC;
    //   RETURN;
        csq
    };
    // ai = s1/SFAC;
    (s1/SFAC, NZ)

};
if significance_loss{
    Err(PartialLossOfSignificance { y: vec![retval.0], nz: retval.1 })

}else {
    Ok(retval)
}
//   170 CONTINUE;
//       AA = 1.0e+3*d1mach(1);
//       S1R = ZEROR;
//       S1I = ZEROI;
//       if (ID == 1) GO TO 190;
//       if (AZ <= AA) GO TO 180;
//       S1R = C2*ZR;
//       S1I = C2*ZI;
//   180 CONTINUE;
//       AIR = C1 - S1R;
//       AII = -S1I;
//       RETURN;
//   190 CONTINUE;
//       AIR = -C2;
//       AII = 0.0;
//       AA = DSQRT(AA);
//       if (AZ <= AA) GO TO 200;
//       S1R = 0.5*(ZR*ZR-ZI*ZI);
//       S1I = ZR*ZI;
//   200 CONTINUE;
//       AIR = AIR + C1*S1R;
//       AII = AII + C1*S1I;
//       RETURN;
//   210 CONTINUE;
//       NZ = 1;
//       AIR = ZEROR;
//       AII = ZEROI;
//       RETURN;
//   270 CONTINUE;
//       NZ = 0;
//       IERR=2;
//       RETURN;
//   280 CONTINUE;
//       if(NN == (-1)) GO TO 270;
//       NZ=0;
//       IERR=5;
//       RETURN;
//   260 CONTINUE;
//       IERR=4;
//       NZ=0;
//       RETURN;
//       END;
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
fn ZBKNU(
    z: Complex64,
    order: f64,
    KODE: Scaling,
    N: usize,
    machine_consts: &MachineConsts,
) -> BesselResult {
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
    let CSCLR = machine_consts.rtol;
    let CRSCR = machine_consts.tol;
    let CSSR = [CSCLR, 1.0, CRSCR];
    let CSRR = [CRSCR, 1.0, CSCLR];
    let BRY = [
        machine_consts.ascle,
        1.0 / machine_consts.ascle,
        f64::MAX / 2.0,
    ];
    let mut NZ = 0;
    let mut underflow_occurred = false;
    let mut KFLAG;
    let mut KODED = KODE;
    let rz = 2.0 * z.conj() / CAZ.powi(2);
    let mut INU = (order + 0.5) as isize; // round to nearest int
    let DNU = order - (INU as f64); // signed fractional part (-0.5 < DNU < 0.5 )
    let DNU2 = if DNU.abs() > machine_consts.tol {
        DNU * DNU
    } else {
        0.0
    };
    let mut skip_to_240= false;
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
            for k in 1..8 {
                ak *= DNU2;
                let TM = CC[k] * ak;
                sum += TM;
                if TM.abs() < machine_consts.tol {
                    break;
                }
            }
            -sum
        } else {
            (T1 - T2) / (DNU + DNU)
        };
        //    50 CONTINUE;
        let G2 = (T1 + T2) * 0.5;
        let f = FC * (G1 * cch + G2 * smu);
        // FR = FC*(CCHR*G1+SMUR*G2);
        // FI = FC*(CCHI*G1+SMUI*G2);
        let mut p = 0.5 * fmu.exp() / T2;
        // CALL ZEXP(FMUR, FMUI, STR, STI);
        // PR = 0.5*STR/T2;
        // PI = 0.5*STI/T2;
        // let pt = 0.5 / fmu.exp();
        // CALL ZDIV(0.5, 0.0, STR, STI, PTR, PTI);
        let mut q = (0.5 / fmu.exp()) / T1;
        // QR = PTR/T1;
        // QI = PTI/T1;
        let mut s1 = f;
        let mut s2 = p;
        // S1R = FR;
        // S1I = FI;
        // S2R = PR;
        // S2I = PI;
        let mut AK = 1.0;
        let mut A1 = 1.0;
        let mut ck = c_one();
        // CKR = CONER;
        // CKI = CONEI;
        let mut BK = 1.0 - DNU2;
        if INU == 0 && N == 1 {
            //GO TO 80;
            //-----------------------------------------------------------------------;
            //     SPECIAL CASE
            //     GENERATE K(FNU,Z), 0.0  <=  FNU  <  0.5 AND N=1;
            //-----------------------------------------------------------------------;
            if CAZ >= machine_consts.tol {
                //GO TO 70;
                let cz = 0.25 * z.powu(2);
                // CALL ZMLT(ZR, ZI, ZR, ZI, CZR, CZI);
                // CZR = 0.25*CZR;
                // CZI = 0.25*CZI;
                let T1 = 0.25 * CAZ * CAZ;
                //    60 CONTINUE;
                '_l60: loop {
                    let f = (f * AK + p + q) / BK;
                    // FR = (FR*AK+PR+QR)/BK;
                    // FI = (FI*AK+PI+QI)/BK;
                    p /= AK - DNU;
                    // STR = 1.0/(AK-DNU);
                    // PR = PR*STR;
                    // PI = PI*STR;
                    q /= AK + DNU;
                    // STR = 1.0/(AK+DNU);
                    // QR = QR*STR;
                    // QI = QI*STR;
                    ck *= cz / AK;
                    // STR = CKR*CZR - CKI*CZI;
                    // RAK = 1.0/AK;
                    // CKI = (CKR*CZI+CKI*CZR)*RAK;
                    // CKR = STR*RAK;
                    s1 += ck * f;
                    // S1R = CKR*FR - CKI*FI + S1R;
                    // S1I = CKR*FI + CKI*FR + S1I;
                    // A1 = A1*T1*RAK;
                    A1 *= T1 / AK;
                    BK += (2.0 * AK) + 1.0;
                    AK += 1.0;
                    if A1 <= machine_consts.tol {
                        break;
                    } //GO TO 60;
                }
            }
            //    70 CONTINUE;
            // y[0] = s1;
            // YR(1) = S1R;
            // YI(1) = S1I;
            let mut y = s1;
            // let y = if (KODED == Scaling::Unscaled) {s1}else{
            //       s1 * z.exp();
            if KODED == Scaling::Scaled {
                y *= z.exp();
            }
            // CALL ZEXP(ZR, ZI, STR, STI);
            // CALL ZMLT(S1R, S1I, STR, STI, YR(1), YI(1));
            // RETURN;}
            return Ok((vec![y], NZ));
            //-----------------------------------------------------------------------;
            //     GENERATE K(DNU,Z) AND K(DNU+1,Z) FOR FORWARD RECURRENCE;
            //-----------------------------------------------------------------------;
        }

        //    80 CONTINUE;
        if CAZ >= machine_consts.tol {
            //GO TO 100;
            // CALL ZMLT(ZR, ZI, ZR, ZI, CZR, CZI);
            let cz = 0.25 * z.powu(2);
            // CZR = 0.25*CZR;
            // CZI = 0.25*CZI;
            let T1 = 0.25 * CAZ * CAZ;
            //    90 CONTINUE;
            '_l90: loop {
                let f = (f * AK + p + q) / BK;
                // FR = (FR*AK+PR+QR)/BK;
                // FI = (FI*AK+PI+QI)/BK;
                p /= AK - DNU;
                // STR = 1.0/(AK-DNU);
                // PR = PR*STR;
                // PI = PI*STR;
                q /= AK + DNU;
                // STR = 1.0/(AK+DNU);
                // QR = QR*STR;
                // QI = QI*STR;
                ck *= cz / AK;
                // STR = CKR*CZR - CKI*CZI;
                // RAK = 1.0/AK;
                // CKI = (CKR*CZI+CKI*CZR)*RAK;
                // CKR = STR*RAK;
                s1 += ck * f;

                // S1R = CKR*FR - CKI*FI + S1R;
                // S1I = CKR*FI + CKI*FR + S1I;
                // ... TODO this is the only bit that differs from the loop above...
                s2 += ck * (p - AK * f);
                // STR = PR - FR*AK;
                // STI = PI - FI*AK;
                // S2R = CKR*STR - CKI*STI + S2R;
                // S2I = CKR*STI + CKI*STR + S2I;
                A1 = A1 * T1 / AK;
                BK = BK + AK + AK + 1.0;
                AK = AK + 1.0;

                if A1 <= machine_consts.tol {
                    break;
                } //GO TO 90;
            }
        }
        //   100 CONTINUE;
        // let KFLAG = 2;
        // A1 = order + 1.0;
        AK = (order + 1.0) * smu.re.abs();
        KFLAG = if AK > machine_consts.alim { 3 } else { 2 };
        // let p2 = s2 * CSSR[KFLAG-1];
        // STR = CSSR(KFLAG);
        // P2R = S2R*STR;
        // P2I = S2I*STR;
        // CALL ZMLT(P2R, P2I, RZR, RZI, S2R, S2I);
        s2 *= CSSR[KFLAG - 1] * rz;
        s1 *= CSSR[KFLAG - 1];
        // S1R = S1R*STR;
        // S1I = S1I*STR;
        if KODED == Scaling::Scaled {
            //GO TO 210;
            let z_exp = z.exp();
            s1 *= z_exp;
            s2 *= z_exp;
        }
        // CALL ZEXP(ZR, ZI, FR, FI);
        // CALL ZMLT(S1R, S1I, FR, FI, S1R, S1I);
        // CALL ZMLT(S2R, S2I, FR, FI, S2R, S2I);
        // GO TO 210;
        (s1, s2)
    } else {
        // alterntive to SERIES FOR CABS(Z) <= R1; Or half integer order
        //-----------------------------------------------------------------------;
        //     underflow_occured=0 MEANS NO UNDERFLOW OCCURRED;
        //     underflow_occured=1 MEANS AN UNDERFLOW OCCURRED- COMPUTATION PROCEEDS WITH;
        //     KODED=2 AND A TEST FOR ON SCALE VALUES IS MADE DURING FORWARD;
        //     RECURSION;
        //-----------------------------------------------------------------------;
        let mut coef = Complex64::new(RTFRAC_PI_2, 0.0) / z.sqrt();
        KFLAG = 2;
        if KODED == Scaling::Unscaled {
            if z.re > machine_consts.alim {
                KODED = Scaling::Scaled;
                underflow_occurred = true;
                KFLAG = 2;
            } else {
                coef *= (-z.re).exp() * CSSR[KFLAG - 1] * Complex64::cis(z.im).conj();
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
            let mut T1 = ((f64::MANTISSA_DIGITS - 1) as f64 * (f64::RADIX as f64).log10() * 3.321928094)
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
                let ETEST = AK / (PI * CAZ * machine_consts.tol);
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
                    FK +=  SPI * T1 * (T2 / CAZ).sqrt();
                    FHS = (0.25 - DNU2).abs();
                }
                (FK, FHS)
            } else {
                //-----------------------------------------------------------------------;
                //     COMPUTE BACKWARD INDEX K FOR CABS(Z) < R2;
                //-----------------------------------------------------------------------;
                AK *= FPI / (machine_consts.tol * CAZ.sqrt().sqrt());
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
            let mut p2 = Complex64::new(machine_consts.tol, 0.0);
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
                if underflow_occurred {
                    // skip_to_270 = true;
                } //GO TO 270;
            // GO TO 240;
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
    // if INU > 0 {
    //     if underflow_occurred {} //GO TO 261;
    // }
    // if underflow_occurred {} //{break 'l270;}

    'l225: loop {
        if !skip_to_240{
        if !(INU <= 0 && N <= 1) {
            //   220 CONTINUE;
            // let INUB = 1;
            // if !underflow_occurred GO TO 261;

            //   225 CONTINUE;
            let mut P1R = CSRR[KFLAG - 1];
            let mut ASCLE = BRY[KFLAG - 1];
            for _ in INUB..=INU {
                let st = s2;
                s2 = ck * s2 + s1;
                s1 = st;
                ck += rz;
                if KFLAG >= 3 {
                    continue;
                }
                let p2 = s2 * P1R;
                if max_abs_component(p2) <= ASCLE {
                    continue;
                }
                KFLAG += 1;
                ASCLE = BRY[KFLAG - 1];
                s1 *= P1R;
                s2 = p2;
                s1 *= CSSR[KFLAG - 1];
                s2 *= CSSR[KFLAG - 1];
                P1R = CSRR[KFLAG - 1];
            }
        } else {
            //-----------------------------------------------------------------------;
            //     underflow_occured=1 CASES, FORWARD RECURRENCE ON SCALED VALUES ON UNDERFLOW;
            //-----------------------------------------------------------------------;
            //   261 CONTINUE;
            let mut cy = c_zeros(2);
            let HELIM = 0.5 * machine_consts.elim;
            let ELM = (-machine_consts.elim).exp();
            let CELMR = ELM;
            let ASCLE = BRY[0]; //BRY(1);
            let mut zd = z;
            // ZDR = ZR;
            // ZDI = ZI;
            let mut IC: isize = -2;
            let mut J = 2;
            // DO 262 I=1,INU;
            // let skip_to_264=false;
            let mut I = 0;
            for i in 0..INU {
                I = i + 1;
                let st = s2;
                s2 = s2 * ck + s1;
                //   STR = S2R;
                //   STI = S2I;
                //   S2R = STR*CKR-STI*CKI+S1R;
                //   S2I = STI*CKR+STR*CKI+S1I;
                //   S1R = STR;
                //   S1I = STI;
                s1 = st;
                ck += rz;
                //   CKR = CKR+RZR;
                //   CKI = CKI+RZI;
                //   AS = ZABS(S2R,S2I);
                let ALAS = s2.abs().ln(); //DLOG(AS);
                //   P2R = -ZDR+ALAS;
                if !((-zd.re + ALAS) < (-machine_consts.elim)) {
                    //GO TO 263;

                    //   CALL ZLOG(S2R,S2I,STR,STI,IDUM);
                    let p2 = -zd + s2.ln();
                    //   P2R = -ZDR+STR;
                    //   P2I = -ZDI+STI;
                    //   P2M = DEXP(P2R)/TOL;
                    //   P1R = P2M*DCOS(P2I);
                    //   P1I = P2M*DSIN(P2I);
                    let p1 = (p2.re.exp() / machine_consts.tol) * Complex64::cis(p2.im);

                    //   CALL ZUunderflowCHK(P1R,P1I,NW,ASCLE,TOL);
                    if will_z_underflow(p1, ASCLE, machine_consts.tol) {
                        //GO TO 263;
                        J = 3 - J;
                        cy[J - 1] = p1;
                        //   CYR(J) = P1R;
                        //   CYI(J) = P1I;
                        // IF(IC.EQ.(I-1)) GO TO 264
                        // below implies we got here twice in a row
                        if IC == (i as isize) - 1 {
                            underflow_occurred = true; //implies 270}//{skip_to_264 = true;break;}//GO TO 264;
                            IC = i as isize;
                            continue;
                        }
                    }
                    //   263   CONTINUE;
                    if ALAS < HELIM {
                        continue;
                    } //GO TO 262;
                    zd.re -= machine_consts.elim;
                    //   ZDR = ZDR-ELIM;
                    s1 *= CELMR;
                    s2 *= CELMR;
                    //   S1R = S1R*CELMR;
                    //   S1I = S1I*CELMR;
                    //   S2R = S2R*CELMR;
                    //   S2I = S2I*CELMR;
                }
            }
            // if !skip_to_264{
            //   262 CONTINUE;
            // if(N == 1) {//GO TO 270;
            //       s1 = s2;
            // }
            // S1R = S2R;
            // S1I = S2I;
            // GO TO 270;
            // }
            //   264 CONTINUE;
            KFLAG = 1;
            INUB = I + 1;
            s2 = cy[J - 1];
            // S2R = CYR(J);
            // S2I = CYI(J);
            J = 3 - J;
            s1 = cy[J - 1];
            // S1R = CYR(J);
            // S1I = CYI(J);
            if INUB <= INU {
                continue 'l225;
            } //GO TO 225;
            // }
            // if(N == 1) {//GO TO 240;
            // s1 = s2;
            // }
            // S1R = S2R;
            // S1I = S2I;
            // GO TO 240;
            // break 'l270;
            //}//loop 270
            //}
        }
        // }
        if N == 1 {
            s1 = s2;
        }
    }
        // ********* basic setup
        let (mut KK, mut y) = if !underflow_occurred {
            let mut y = c_zeros(N);
            y[0] = s1 * CSRR[KFLAG - 1];
            if N == 1 {
                return Ok((y, NZ));
            }
            y[1] = s2 * CSRR[KFLAG - 1];
            if N == 2 {
                return Ok((y, NZ));
            }
            let KK = 1;
            (KK, y)
        // ********* End Basic Setup
        } else {
            //Complex setup from 270 onwards
            // ******** Alternative setup if underflow_occured
            // 270 CONTINUE;
            let mut y = c_zeros(N);
            y[0] = s1;
            if N > 1 {
                y[1] = s2;
            }
            // YR(1) = S1R;
            // YI(1) = S1I;
            // if(N == 1) GO TO 280;
            // YR(2) = S2R;
            // YI(2) = S2I;
            //   280 CONTINUE;
            // ASCLE = BRY[0];
            ZKSCL(zd, order, N, &mut y, &mut NZ, rz, BRY[0], machine_consts);
            // CALL ZKSCL(ZDR,ZDI,FNU,N,YR,YI,NZ,RZR,RZI,ASCLE,TOL,ELIM);
            INU = (N - NZ) as isize;
            if INU <= 0 {
                return Ok((y, NZ));
            } //;RETURN;}
            let mut KK = NZ; // + 1;
            y[KK - 1] *= CSRR[0];
            // S1R = YR(KK);
            // S1I = YI(KK);
            // YR(KK) = S1R*CSRR(1);
            // YI(KK) = S1I*CSRR(1);
            if INU == 1 {
                return Ok((y, NZ));
            }
            KK = NZ + 1; //+ 2;
            y[KK] *= CSRR[0];
            // S2R = YR(KK);
            // S2I = YI(KK);
            // YR(KK) = S2R*CSRR(1);
            // YI(KK) = S2I*CSRR(1);
            if INU == 2 {
                return Ok((y, NZ));
            }
            ck = (order + ((KK - 1) as f64)) * rz;
            // T2 = order + ((KK-1) as f64);
            // CKR = T2*RZR;
            // CKI = T2*RZI;
            KFLAG = 1;
            // GO TO 250;
            (KK, y)
        };
        //   250 CONTINUE;
        KK += 1;
        if KK >= N {
            return Ok((y, NZ));
        }
        let mut P1R = CSRR[KFLAG - 1];
        let mut ASCLE = BRY[KFLAG - 1];
        // DO 260 I=KK,N;
        let mut I;
        for i in KK..N {
            I = i;
            let mut p2 = s2;
            //   P2R = S2R;
            //   P2I = S2I;
            s2 = ck * s2 + s1;
            //   S2R = CKR*P2R - CKI*P2I + S1R;
            //   S2I = CKI*P2R + CKR*P2I + S1I;
            let mut s1 = p2;
            //   S1R = P2R;
            //   S1I = P2I;
            ck *= rz;
            //   CKR = CKR + RZR;
            //   CKI = CKI + RZI;
            p2 = s2 * P1R;
            //   P2R = S2R*P1R;
            //   P2I = S2I*P1R;
            y[I - 1] = p2;
            //   YR(I) = P2R;
            //   YI(I) = P2I;
            if KFLAG >= 3 {
                return Err(LossOfSignificance);
            };
            //   STR = (P2R).abs();
            //   STI = (P2I).abs();
            //   P2M = DMAX1(STR,STI);
            if max_abs_component(p2) <= ASCLE {
                continue;
            } //GO TO 260;
            KFLAG += 1;
            ASCLE = BRY[KFLAG - 1];
            s1 *= P1R;
            //   S1R = S1R*P1R;
            //   S1I = S1I*P1R;
            s2 = p2;
            //   S2R = P2R;
            //   S2I = P2I;
            s1 *= CSSR[KFLAG - 1];
            s2 *= CSSR[KFLAG - 1];
            //   STR = CSSR(KFLAG);
            //   S1R = S1R*STR;
            //   S1I = S1I*STR;
            //   S2R = S2R*STR;
            //   S2I = S2I*STR;
            P1R = CSRR[KFLAG - 1];
            //   P1R = CSRR(KFLAG);
        }
        //   260 CONTINUE;
        // RETURN;
        return Ok((y, NZ));
    }
}

fn ZKSCL(
    //ZRR,ZRI,FNU,
    zr: Complex64,
    order: f64,
    N: usize, //YR,YI,NZ,
    y: &mut Vec<Complex64>,
    NZ: &mut usize,
    rz: Complex64, //RZR,RZI,
    ASCLE: f64,
    machine_consts: &MachineConsts,
) //-> BesselResult//ASCLE,TOL,ELIM)
{
    // ***BEGIN PROLOGUE  ZKSCL
    // ***REFER TO  ZBESK
    //
    //     SET K FUNCTIONS TO ZERO ON UNDERFLOW, CONTINUE RECURRENCE
    //     ON SCALED FUNCTIONS UNTIL TWO MEMBERS COME ON SCALE, THEN
    //     RETURN WITH MIN(NZ+2,N) VALUES SCALED BY 1/TOL.
    //
    // ***ROUTINES CALLED  ZUunderflowCHK,ZABS,ZLOG
    // ***END PROLOGUE  ZKSCL
    //     COMPLEX CK,CS,CY,CZERO,RZ,S1,S2,Y,ZR,ZD,CELM
    //       EXTERNAL ZABS
    //       DOUBLE PRECISION ACS, AS, ASCLE, CKI, CKR, CSI, CSR, CYI,
    //      * CYR, ELIM, FN, FNU, RZI, RZR, STR, S1I, S1R, S2I,
    //      * S2R, TOL, YI, YR, ZEROI, ZEROR, ZRI, ZRR, ZABS,
    //      * ZDR, ZDI, CELMR, ELM, HELIM, ALAS
    //       INTEGER I, IC, IDUM, KK, N, NN, NW, NZ
    //       DIMENSION YR(N), YI(N), CYR(2), CYI(2)
    //       DATA ZEROR,ZEROI / 0.0 , 0.0 /
    //
    *NZ = 0;
    // let mut IC = 0;
    //NN = MIN0(2,N);
    let NN = min(2, N);
    // DO 10 I=1,NN;
    let mut cy = c_zeros(2);
    let mut IC = 0;
    for i in 0..NN {
        let s1 = y[i];
        //   S1R = YR(I);
        //   S1I = YI(I);
        cy[i] = s1;
        //   CYR(I) = S1R;
        //   CYI(I) = S1I;
        //   let AS = s1.abs();//ZABS(S1R,S1I);
        //   let acs = -zr.re +  s1.abs().ln();
        //   ACS = -ZRR + DLOG(AS);
        *NZ += 1;
        //   NZ = NZ + 1;
        y[i] = c_zero();
        //   YR(I) = ZEROR;
        //   YI(I) = ZEROI;
        if -zr.re + s1.abs().ln() < (-machine_consts.elim) {
            continue;
        } //GO TO 10;
        let mut cs = s1.ln() - zr;
        //   CALL ZLOG(S1R, S1I, CSR, CSI, IDUM);
        //   CSR = CSR - ZRR;
        //   CSI = CSI - ZRI;
        cs = (cs.re.exp() / machine_consts.tol) * Complex64::cis(cs.im);
        //   STR = DEXP(CSR)/TOL;
        //   CSR = STR*DCOS(CSI);
        //   CSI = STR*DSIN(CSI);
        //   CALL ZUunderflowCHK(CSR, CSI, NW, ASCLE, TOL);

        if will_z_underflow(cs, ASCLE, machine_consts.tol) {
            continue;
        } //GO TO 10;
        y[i] = cs;
        //   YR(I) = CSR;
        //   YI(I) = CSI;
        IC = i;
        *NZ -= 1;
    }
    //    10 CONTINUE;
    if N == 1 {
        return;
    }
    if IC <= 1 {
        //GO TO 20;
        y[0] = c_zero();
        // YR(1) = ZEROR;
        // YI(1) = ZEROI;
        *NZ = 2;
    }
    //    20 CONTINUE;
    if N == 2 {
        return;
    }
    if *NZ == 0 {
        return;
    }
    let FN = order + 1.0;
    let mut ck = FN * rz;
    // CKR = FN*RZR;
    // CKI = FN*RZI;
    let mut s1 = cy[0];
    // S1R = CYR(1);
    // S1I = CYI(1);
    let mut s2 = cy[1];
    // S2R = CYR(2);
    // S2I = CYI(2);
    let half_elim = 0.5 * machine_consts.elim;
    let ELM = (-machine_consts.elim).exp();
    let CELMR = ELM;
    let mut zd = zr;
    // ZDR = ZRR;
    // ZDI = ZRI;
    //;
    //     FIND TWO CONSECUTIVE Y VALUES ON SCALE. SCALE RECURRENCE if;
    //     S2 GETS LARGER THAN EXP(ELIM/2);
    //;
    let mut skip_to_40 = false;
    // DO 30 I=3,N;
    let mut KK = 0;
    for i in 2..N {
        KK = i + 1;
        let mut cs = s2;
        //   CSR = S2R;
        //   CSI = S2I;
        s2 = cs * ck + s1;
        //   S2R = CKR*CSR - CKI*CSI + S1R;
        //   S2I = CKI*CSR + CKR*CSI + S1I;
        s1 = cs;
        //   S1R = CSR;
        //   S1I = CSI;
        ck += rz;
        //   CKR = CKR + RZR;
        //   CKI = CKI + RZI;
        //   AS = ZABS(S2R,S2I);
        //   ALAS = DLOG(AS);
        let ALAS = s2.abs().ln();
        //   ACS = -ZDR + ALAS;
        *NZ += 1;
        y[i] = Complex64::zero();
        //   YR(I) = ZEROR;
        //   YI(I) = ZEROI;
        if !(-zd.re + s2.abs().ln() < (-machine_consts.elim)) {
            //GO TO 25;
            //   CALL ZLOG(S2R, S2I, CSR, CSI, IDUM);
            cs = s2.ln() - zd;
            //   CSR = CSR - ZDR;
            //   CSI = CSI - ZDI;
            cs = (cs.exp() / machine_consts.tol) * Complex64::cis(cs.im);
            //   STR = DEXP(CSR)/TOL;
            //   CSR = STR*DCOS(CSI);
            //   CSI = STR*DSIN(CSI);
            //   CALL ZUunderflowCHK(CSR, CSI, NW, ASCLE, TOL);
            if !will_z_underflow(cs, ASCLE, machine_consts.tol) {
                //GO TO 25;
                y[i] = cs;
                //   YR(I) = CSR;
                //   YI(I) = CSI;
                *NZ -= 1;
                //   NZ = NZ - 1;
                if IC == KK - 1 {
                    skip_to_40 = true;
                    break;
                } //GO TO 40;
                IC = KK;
                //   GO TO 30;
                continue;
            }
        }

        //    25   CONTINUE;
        if ALAS < half_elim {
            continue;
        } //GO TO 30;}
        zd -= machine_consts.elim;
        //   ZDR = ZDR - ELIM;
        s1 *= CELMR;
        s2 *= CELMR;
        //   S1R = S1R*CELMR;
        //   S1I = S1I*CELMR;
        //   S2R = S2R*CELMR;
        //   S2I = S2I*CELMR;
    }
    //    30 CONTINUE;
    if !skip_to_40 {
        *NZ = N;
        if IC == N {
            *NZ = N - 1
        };

        // GO TO 45;
    } else {
        //    40 CONTINUE;
        *NZ = KK - 2;
    }
    //    45 CONTINUE;
    for i in 0..*NZ {
        // DO 50 I=1,NZ;
        y[i] = c_zero();
        //   YR(I) = ZEROR;
        //   YI(I) = ZEROI;
    }
    //    50 CONTINUE;
    // RETURN;
    // END;
}
/*
fn ZSHCH(ZR, ZI, CSHR, CSHI, CCHR, CCHI)
// ***BEGIN PROLOGUE  ZSHCH
// ***REFER TO  ZBESK,ZBESH
//
//     ZSHCH COMPUTES THE COMPLEX HYPERBOLIC FUNCTIONS CSH=SINH(X+I*Y)
//     AND CCH=COSH(X+I*Y), WHERE I**2=-1.
//
// ***ROUTINES CALLED  (NONE)
// ***END PROLOGUE  ZSHCH
//
      DOUBLE PRECISION CCHI, CCHR, CH, CN, CSHI, CSHR, SH, SN, ZI, ZR,
     * DCOSH, DSINH
      SH = DSINH(ZR)
      CH = DCOSH(ZR)
      SN = DSIN(ZI)
      CN = DCOS(ZI)
      CSHR = SH*CN
      CSHI = CH*SN
      CCHR = CH*CN
      CCHI = SH*SN
      RETURN
      END
      */
fn i_ratios(
    //ZR, ZI, FNU,
    z: Complex64,
    order: f64,
    N: usize, //CYR, CYI, TOL
    machine_consts: &MachineConsts,
) -> Vec<Complex64> {
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
    //-----------------------------------------------------------------------;
    //     THE OVERFLOW TEST ON K(FNU+I-1,Z) BEFORE THE CALL TO CBKNU
    //     GUARANTEES THAT P2 IS ON SCALE. SCALE TEST1 AND ALL SUBSEQUENT
    //     P2 VALUES BY AP1 TO ENSURE THAT AN OVERFLOW DOES NOT OCCUR
    //     PREMATURELY.;
    //-----------------------------------------------------------------------;
    let ARG = (AP2 + AP2) / (AP1 * machine_consts.tol);
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
        p1 = Complex64::new(machine_consts.tol, machine_consts.tol);
    }
    let mut cy = c_zeros(N);
    cy[N - 1] = p2 / p1;
    if N == 1 {
        return cy;
    }
    t1 = Complex64::new((N - 1) as f64, 0.0);
    let cdfnu = order * rz;
    for k in (1..N).rev() {
        let mut pt = cdfnu + t1 * rz + cy[k];
        let mut AK = pt.abs();
        if AK == 0.0 {
            pt = Complex64::new(machine_consts.tol, machine_consts.tol);
            AK = pt.abs();
        }
        cy[k - 1] = pt.conj() / AK.powi(2);
        t1 -= 1.0;
    }
    return cy;
}

fn ZS1S2(//ZRR, ZRI, S1R, S1I, S2R, S2I, NZ, ASCLE, ALIM,
     //* IUF)
     zr: Complex64, s1:&mut Complex64, s2:&mut Complex64,  IUF:&mut usize, machine_consts: &MachineConsts
)-> usize{
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
//
// ***ROUTINES CALLED  ZABS,ZEXP,ZLOG
// ***END PROLOGUE  ZS1S2
//     COMPLEX CZERO,C1,S1,S1D,S2,ZR
    //   EXTERNAL ZABS
    //   DOUBLE PRECISION AA, ALIM, ALN, ASCLE, AS1, AS2, C1I, C1R, S1DI,
    //  * S1DR, S1I, S1R, S2I, S2R, ZEROI, ZEROR, ZRI, ZRR, ZABS
    //   INTEGER IUF, IDUM, NZ
    //   DATA ZEROR,ZEROI  / 0.0 , 0.0 /
      let NZ = 0;
      let mut abs_s1 = s1.abs();
      let abs_s2  =s2.abs();
    //   AS1 = ZABS(S1R,S1I);
    //   AS2 = ZABS(S2R,S2I);
    if (s1.re != 0.0 || s1.im != 0.0) && (abs_s1 !=0.0){


    //   if (S1R == 0.0 && S1I == 0.0) GO TO 10;
    //   if (AS1 == 0.0) GO TO 10;
      let ALN = (-2.0*zr.re) + abs_s1.ln();//-ZRR - ZRR + DLOG(AS1);
        let s1d = *s1;
        *s1 = c_zero();
        abs_s1 = 0.0;
    //   S1DR = S1R;
    //   S1DI = S1I;
    //   S1R = ZEROR;
    //   S1I = ZEROI;
    //   AS1 = ZEROR;
      if ALN >= (-machine_consts.alim) {//} GO TO 10;
        *s1 = (s1d.ln() - 2.0*zr).exp();
        abs_s1 = s1.abs();
    //   CALL ZLOG(S1DR, S1DI, C1R, C1I, IDUM);
    //   C1R = C1R - ZRR - ZRR;
    //   C1I = C1I - ZRI - ZRI;
    //   CALL ZEXP(C1R, C1I, S1R, S1I);
    //   AS1 = ZABS(S1R,S1I);
      *IUF += 1;
      }
    }
//    10 CONTINUE;

    //   AA = DMAX1(AS1,AS2);
      if abs_s1.max(abs_s2) > machine_consts.ascle { NZ} //RETURN;
      else{
        *s1 = c_zero();
        *s2 = c_zero();
    //   S1R = ZEROR;
    //   S1I = ZEROI;
    //   S2R = ZEROR;
    //   S2I = ZEROI;
    //   NZ = 1;
      *IUF = 0;
    //   RETURN;
    //   END;
      1
    }
}
      /*
fn ZBUNK(ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, TOL, ELIM,
     * ALIM)
// ***BEGIN PROLOGUE  ZBUNK
// ***REFER TO  ZBESK,ZBESH
//
//     ZBUNK COMPUTES THE K BESSEL FUNCTION FOR FNU > FNUL.
//     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR K(FNU,Z)
//     IN ZUNK1 AND THE EXPANSION FOR H(2,FNU,Z) IN ZUNK2
//
// ***ROUTINES CALLED  ZUNK1,ZUNK2
// ***END PROLOGUE  ZBUNK
//     COMPLEX Y,Z
      DOUBLE PRECISION ALIM, AX, AY, ELIM, FNU, TOL, YI, YR, ZI, ZR
      INTEGER KODE, MR, N, NZ
      DIMENSION YR(N), YI(N)
      NZ = 0
      AX = (ZR).abs()*1.7321
      AY = (ZI).abs()
      if (AY > AX) GO TO 10
//-----------------------------------------------------------------------
//     ASYMPTOTIC EXPANSION FOR K(FNU,Z) FOR LARGE FNU APPLIED IN
//     -PI/3 <= ARG(Z) <= PI/3
//-----------------------------------------------------------------------
      CALL ZUNK1(ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, TOL, ELIM, ALIM)
      GO TO 20
   10 CONTINUE
//-----------------------------------------------------------------------
//     ASYMPTOTIC EXPANSION FOR H(2,FNU,Z*EXP(M*FRAC_PI_2)) FOR LARGE FNU
//     APPLIED IN PI/3 < ABS(ARG(Z)) <= PI/2 WHERE M=+I OR -I
//     AND FRAC_PI_2=PI/2
//-----------------------------------------------------------------------
      CALL ZUNK2(ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, TOL, ELIM, ALIM)
   20 CONTINUE
      RETURN
      END
      */
fn i_miller(
    z: Complex64,
    order: f64, //ZR, ZI, FNU,
    KODE: Scaling,
    N: usize, //YR, YI, NZ,
    machine_consts: &MachineConsts,
) -> BesselResult {
    // ***BEGIN PROLOGUE  ZMLRI
    // ***REFER TO  ZBESI,ZBESK
    //
    //     ZMLRI COMPUTES THE I BESSEL FUNCTION FOR RE(Z) >= 0.0 BY THE
    //     MILLER ALGORITHM NORMALIZED BY A NEUMANN SERIES.
    //
    // ***ROUTINES CALLED  gamma_ln,d1mach,ZABS,ZEXP,ZLOG,ZMLT
    // ***END PROLOGUE  ZMLRI

    let SCLE: f64 = 2.0 * f64::MIN_POSITIVE / machine_consts.tol;
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
    TST /= machine_consts.tol;
    //-----------------------------------------------------------------------;
    //     COMPUTE RELATIVE TRUNCATION ERROR INDEX FOR SERIES;
    //-----------------------------------------------------------------------;
    let mut AK = AT;
    let mut converged = false;
    let mut I = 0;
    for i in 0..80 {
        I = i;
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
        //-----------------------------------------------------------------------;
        //     COMPUTE RELATIVE TRUNCATION ERROR FOR RATIOS;
        //-----------------------------------------------------------------------;
        p1 = c_zero();
        p2 = c_one();
        let AT = (INU as f64) + 1.0;
        ck = z.conj() * RAZ * RAZ * AT;
        ACK = AT * RAZ;
        TST = (ACK / machine_consts.tol).sqrt();
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
    //-----------------------------------------------------------------------;
    //     BACKWARD RECURRENCE AND SUM NORMALIZING RELATION;
    //-----------------------------------------------------------------------;
    K += 1;
    let KK = (I + IAZ).max(K + INU);
    let mut FKK = KK as f64;
    let mut p1 = c_zero();
    //-----------------------------------------------------------------------;
    //     SCALE P2 AND SUM BY SCLE;
    //-----------------------------------------------------------------------;
    let mut p2 = Complex64::new(SCLE, 0.0);
    let FNF = order - (IFNU as f64);
    let TFNF = FNF + FNF;
    let mut BK = (gamma_ln(FKK + TFNF + 1.0).unwrap()
        - gamma_ln(FKK + 1.0).unwrap()
        - gamma_ln(TFNF + 1.0).unwrap())
    .exp();
    let mut sumr = c_zero();
    for _i in 0..(KK - INU) {
        let pt = p2;
        p2 = p1 + (FKK + FNF) * (rz * pt);
        p1 = pt;
        AK = 1.0 - TFNF / (FKK + TFNF);
        ACK = BK * AK;
        sumr += (ACK + BK) * p1;
        BK = ACK;
        FKK = FKK - 1.0;
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
            FKK = FKK - 1.0;
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
            FKK = FKK - 1.0;
        }
    }

    let mut pt = z;
    if KODE == Scaling::Scaled {
        pt.re = 0.0;
    }
    p1 = -FNF * rz.ln() + pt;
    let AP = gamma_ln(1.0 + FNF).unwrap();
    p1 -= AP;
    //-----------------------------------------------------------------------;
    //     THE DIVISION CEXP(PT)/(SUM+P2) IS ALTERED TO AVOID OVERFLOW;
    //     IN THE DENOMINATOR BY SQUARING LARGE QUANTITIES;
    //-----------------------------------------------------------------------;
    p2 += sumr;
    let AP = p2.abs();
    ck = p1.exp() / AP;
    let cnorm = ck * p2.conj() / AP;
    for i in 0..N {
        y[i] *= cnorm;
    }
    return Ok((y, NZ));
}

fn i_wronksian(
    zr: Complex64,
    order: f64,
    KODE: Scaling,
    N: usize,
    y: &mut Vec<Complex64>,
    machine_consts: &MachineConsts,
) -> BesselResult<usize> {
    // ***BEGIN PROLOGUE  ZWRSK
    // ***REFER TO  ZBESI,ZBESK
    //
    //     ZWRSK COMPUTES THE I BESSEL FUNCTION FOR RE(Z) >= 0.0 BY
    //     NORMALIZING THE I FUNCTION RATIOS FROM ZRATI BY THE WRONSKIAN
    //
    // ***ROUTINES CALLED  d1mach,ZBKNU,ZRATI,ZABS
    // ***END PROLOGUE  ZWRSK
    //-----------------------------------------------------------------------
    //     I(FNU+I-1,Z) BY BACKWARD RECURRENCE FOR RATIOS
    //     Y(I)=I(FNU+I,Z)/I(FNU+I-1,Z) FROM CRATI NORMALIZED BY THE
    //     WRONSKIAN WITH K(FNU,Z) AND K(FNU+1,Z) FROM CBKNU.
    //-----------------------------------------------------------------------
    let NZ = 0;
    let (cw, _) = ZBKNU(zr, order, KODE, 2, machine_consts)?;
    let y_ratios = i_ratios(zr, order, N, machine_consts);
    //-----------------------------------------------------------------------;
    //     RECUR FORWARD ON I(FNU+1,Z) = R(FNU,Z)*I(FNU,Z),;
    //     R(FNU+J-1,Z)=Y(J),  J=1,...,N;
    //-----------------------------------------------------------------------;
    let mut cinu = c_one();
    if KODE == Scaling::Scaled {
        cinu = Complex64::cis(zr.im);
    }
    //-----------------------------------------------------------------------;
    //     ON LOW EXPONENT MACHINES THE K FUNCTIONS CAN BE CLOSE TO BOTH;
    //     THE UNDER AND OVERFLOW LIMITS AND THE NORMALIZATION MUST BE;
    //     SCALED TO PREVENT OVER OR UNDERFLOW. CUOIK HAS DETERMINED THAT;
    //     THE RESULT IS ON SCALE.;
    //-----------------------------------------------------------------------;
    let acw = cw[1].abs();
    let CSCLR = if acw <= machine_consts.ascle {
        1.0 / machine_consts.tol
    } else if acw >= 1.0 / machine_consts.ascle {
        machine_consts.tol
    } else {
        1.0
    };
    let c1 = cw[0] * CSCLR;
    let c2 = cw[1] * CSCLR;
    //-----------------------------------------------------------------------;
    //     CINU=CINU*(CONJG(CT)/CABS(CT))*(1.0/CABS(CT) PREVENTS;
    //     UNDER- OR OVERFLOW PREMATURELY BY SQUARING CABS(CT);
    //-----------------------------------------------------------------------;
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

/*
fn ZACON(ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, RL, FNUL,
     * TOL, ELIM, ALIM)
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
// ***ROUTINES CALLED  ZBINU,ZBKNU,ZS1S2,d1mach,ZABS,ZMLT
// ***END PROLOGUE  ZACON
//     COMPLEX CK,CONE,CSCL,CSCR,CSGN,CSPN,CY,CZERO,C1,C2,RZ,SC1,SC2,ST,
//    *S1,S2,Y,Z,ZN
      EXTERNAL ZABS
      DOUBLE PRECISION ALIM, ARG, ASCLE, AS2, AZN, BRY, BSCLE, CKI,
     * CKR, CONER, CPN, CSCL, CSCR, CSGNI, CSGNR, CSPNI, CSPNR,
     * CSR, CSRR, CSSR, CYI, CYR, C1I, C1M, C1R, C2I, C2R, ELIM, FMR,
     * FN, FNU, FNUL, PI, PTI, PTR, RAZN, RL, RZI, RZR, SC1I, SC1R,
     * SC2I, SC2R, SGN, SPN, STI, STR, S1I, S1R, S2I, S2R, TOL, YI, YR,
     * YY, ZEROR, ZI, ZNI, ZNR, ZR, d1mach, ZABS
      INTEGER I, INU, IUF, KFLAG, KODE, MR, N, NN, NW, NZ
      DIMENSION YR(N), YI(N), CYR(2), CYI(2), CSSR(3), CSRR(3), BRY(3)
      DATA PI / 3.14159265358979324 /
      DATA ZEROR,CONER / 0.0,1.0 /
      NZ = 0
      ZNR = -ZR
      ZNI = -ZI
      NN = N
      CALL ZBINU(ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, RL, FNUL, TOL,
     * ELIM, ALIM)
      if (NW < 0) GO TO 90
//-----------------------------------------------------------------------
//     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
//-----------------------------------------------------------------------
      NN = MIN0(2,N)
      CALL ZBKNU(ZNR, ZNI, FNU, KODE, NN, CYR, CYI, NW, TOL, ELIM, ALIM)
      if (NW != 0) GO TO 90
      S1R = CYR(1)
      S1I = CYI(1)
      FMR = (MR as f64)
      SGN = -DSIGN(PI,FMR)
      CSGNR = ZEROR
      CSGNI = SGN
      if (KODE == 1) GO TO 10
      YY = -ZNI
      CPN = DCOS(YY)
      SPN = DSIN(YY)
      CALL ZMLT(CSGNR, CSGNI, CPN, SPN, CSGNR, CSGNI)
   10 CONTINUE
//-----------------------------------------------------------------------
//     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
//     WHEN FNU IS LARGE
//-----------------------------------------------------------------------
      INU = INT(SNGL(FNU))
      ARG = (FNU-(INU as f64))*SGN
      CPN = DCOS(ARG)
      SPN = DSIN(ARG)
      CSPNR = CPN
      CSPNI = SPN
      if (MOD(INU,2) == 0) GO TO 20
      CSPNR = -CSPNR
      CSPNI = -CSPNI
   20 CONTINUE
      IUF = 0
      C1R = S1R
      C1I = S1I
      C2R = YR(1)
      C2I = YI(1)
      ASCLE = 1.0e+3*d1mach(1)/TOL
      if (KODE == 1) GO TO 30
      CALL ZS1S2(ZNR, ZNI, C1R, C1I, C2R, C2I, NW, ASCLE, ALIM, IUF)
      NZ = NZ + NW
      SC1R = C1R
      SC1I = C1I
   30 CONTINUE
      CALL ZMLT(CSPNR, CSPNI, C1R, C1I, STR, STI)
      CALL ZMLT(CSGNR, CSGNI, C2R, C2I, PTR, PTI)
      YR(1) = STR + PTR
      YI(1) = STI + PTI
      if (N == 1) RETURN
      CSPNR = -CSPNR
      CSPNI = -CSPNI
      S2R = CYR(2)
      S2I = CYI(2)
      C1R = S2R
      C1I = S2I
      C2R = YR(2)
      C2I = YI(2)
      if (KODE == 1) GO TO 40
      CALL ZS1S2(ZNR, ZNI, C1R, C1I, C2R, C2I, NW, ASCLE, ALIM, IUF)
      NZ = NZ + NW
      SC2R = C1R
      SC2I = C1I
   40 CONTINUE
      CALL ZMLT(CSPNR, CSPNI, C1R, C1I, STR, STI)
      CALL ZMLT(CSGNR, CSGNI, C2R, C2I, PTR, PTI)
      YR(2) = STR + PTR
      YI(2) = STI + PTI
      if (N == 2) RETURN
      CSPNR = -CSPNR
      CSPNI = -CSPNI
      AZN = ZABS(ZNR,ZNI)
      RAZN = 1.0/AZN
      STR = ZNR*RAZN
      STI = -ZNI*RAZN
      RZR = (STR+STR)*RAZN
      RZI = (STI+STI)*RAZN
      FN = FNU + 1.0
      CKR = FN*RZR
      CKI = FN*RZI
//-----------------------------------------------------------------------
//     SCALE NEAR EXPONENT EXTREMES DURING RECURRENCE ON K FUNCTIONS
//-----------------------------------------------------------------------
      CSCL = 1.0/TOL
      CSCR = TOL
      CSSR(1) = CSCL
      CSSR(2) = CONER
      CSSR(3) = CSCR
      CSRR(1) = CSCR
      CSRR(2) = CONER
      CSRR(3) = CSCL
      BRY(1) = ASCLE
      BRY(2) = 1.0/ASCLE
      BRY(3) = d1mach(2)
      AS2 = ZABS(S2R,S2I)
      KFLAG = 2
      if (AS2 > BRY(1)) GO TO 50
      KFLAG = 1
      GO TO 60
   50 CONTINUE
      if (AS2 < BRY(2)) GO TO 60
      KFLAG = 3
   60 CONTINUE
      BSCLE = BRY(KFLAG)
      S1R = S1R*CSSR(KFLAG)
      S1I = S1I*CSSR(KFLAG)
      S2R = S2R*CSSR(KFLAG)
      S2I = S2I*CSSR(KFLAG)
      CSR = CSRR(KFLAG)
      DO 80 I=3,N
        STR = S2R
        STI = S2I
        S2R = CKR*STR - CKI*STI + S1R
        S2I = CKR*STI + CKI*STR + S1I
        S1R = STR
        S1I = STI
        C1R = S2R*CSR
        C1I = S2I*CSR
        STR = C1R
        STI = C1I
        C2R = YR(I)
        C2I = YI(I)
        if (KODE == 1) GO TO 70
        if (IUF < 0) GO TO 70
        CALL ZS1S2(ZNR, ZNI, C1R, C1I, C2R, C2I, NW, ASCLE, ALIM, IUF)
        NZ = NZ + NW
        SC1R = SC2R
        SC1I = SC2I
        SC2R = C1R
        SC2I = C1I
        if (IUF != 3) GO TO 70
        IUF = -4
        S1R = SC1R*CSSR(KFLAG)
        S1I = SC1I*CSSR(KFLAG)
        S2R = SC2R*CSSR(KFLAG)
        S2I = SC2I*CSSR(KFLAG)
        STR = SC2R
        STI = SC2I
   70   CONTINUE
        PTR = CSPNR*C1R - CSPNI*C1I
        PTI = CSPNR*C1I + CSPNI*C1R
        YR(I) = PTR + CSGNR*C2R - CSGNI*C2I
        YI(I) = PTI + CSGNR*C2I + CSGNI*C2R
        CKR = CKR + RZR
        CKI = CKI + RZI
        CSPNR = -CSPNR
        CSPNI = -CSPNI
        if (KFLAG >= 3) GO TO 80
        PTR = (C1R).abs()
        PTI = (C1I).abs()
        C1M = DMAX1(PTR,PTI)
        if (C1M <= BSCLE) GO TO 80
        KFLAG = KFLAG + 1
        BSCLE = BRY(KFLAG)
        S1R = S1R*CSR
        S1I = S1I*CSR
        S2R = STR
        S2I = STI
        S1R = S1R*CSSR(KFLAG)
        S1I = S1I*CSSR(KFLAG)
        S2R = S2R*CSSR(KFLAG)
        S2I = S2I*CSSR(KFLAG)
        CSR = CSRR(KFLAG)
   80 CONTINUE
      RETURN
   90 CONTINUE
      NZ = -1
      if(NW == (-2)) NZ=-2
      RETURN
      END
*/
fn ZBINU(
    z: Complex64, //ZR, ZI,
    order: f64,
    KODE: Scaling,
    N: usize,
    machine_consts: &MachineConsts,
) -> BesselResult {
    //output?: CYR, CYI, NZ,)
    // ***BEGIN PROLOGUE  ZBINU
    // ***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZAIRY,ZBIRY
    //
    //     ZBINU COMPUTES THE I FUNCTION IN THE RIGHT HALF Z PLANE
    //
    // ***ROUTINES CALLED  ZABS,ZASYI,ZBUNI,ZMLRI,z_power_series,ZUOIK,ZWRSK
    // ***END PROLOGUE  ZBINU
    //       EXTERNAL ZABS
    //       DOUBLE PRECISION ALIM, AZ, CWI, CWR, CYI, CYR, DFNU, ELIM, FNU,
    //      * FNUL, RL, TOL, ZEROI, ZEROR, ZI, ZR, ZABS
    //       INTEGER I, INW, KODE, N, NLAST, NN, NUI, NW, NZ
    //       DIMENSION CYR(N), CYI(N), CWR(2), CWI(2)
    //       DATA ZEROR,ZEROI / 0.0, 0.0 /
    //
    let mut NZ = 0;
    let AZ = z.abs(); //ZABS(ZR,ZI)
    let mut NN: usize = N;
    let mut DFNU = order + ((N - 1) as f64);
    let mut cy = c_zeros(N);
    if AZ <= 2.0 || AZ * AZ * 0.25 <= DFNU + 1.0 {
        //-----------------------------------------------------------------------
        //     POWER SERIES
        //-----------------------------------------------------------------------
        let NW;
        (cy, NW) = i_power_series(z, order, KODE, NN, machine_consts)?;
        let INW: usize = NW.abs().try_into().unwrap();
        NZ = NZ + INW;
        NN = NN - INW;
        if NN == 0 || NW >= 0 {
            return Ok((cy, NZ.try_into().unwrap()));
        }

        DFNU = order + ((NN as f64) - 1.0);
    }

    if (AZ >=  machine_consts.rl)
          && ((DFNU <= 1.0) //GO TO 30
          || (AZ+AZ >= DFNU*DFNU))
    //GO TO 50 //equiv to go to 40 as (DFNU <= 1.0) is not true to get here
    {
        //     if (AZ < RL) //GO TO 40
        //     if (DFNU <= 1.0) GO TO 30
        //     if (AZ+AZ < DFNU*DFNU) GO TO 50 //equiv to go to 40 as (DFNU <= 1.0) is not true to get here
        //-----------------------------------------------------------------------
        //     ASYMPTOTIC EXPANSION FOR LARGE Z
        //-----------------------------------------------------------------------
        //  30 CONTINUE
        let (cy, nw) = z_asymptotic_i(z, order, KODE, NN, machine_consts)?;
        //     if (NW < 0) GO TO 130
        debug_assert!(nw == NZ);
        return Ok((cy, NZ));
    }
    //  40 CONTINUE
    let mut skip_az_rl_check = true;
    if !(DFNU <= 1.0) {
        //GO TO 70
        skip_az_rl_check = false;
        //  50 CONTINUE

        //-----------------------------------------------------------------------
        //     OVERFLOW AND UNDERFLOW TEST ON I SEQUENCE FOR MILLER ALGORITHM
        //-----------------------------------------------------------------------
        let nw = zuoik(z, order, KODE, IKType::I, NN, &mut cy, machine_consts)?;

        NZ = NZ + nw;
        NN = NN - nw;
        if NN == 0 {
            return Ok((cy, NZ));
        }
        DFNU = order + ((NN - 1) as f64);
        //     if (DFNU > FNUL) GO TO 110
        //     if (AZ > FNUL) GO TO 110
    }
    if (DFNU > machine_consts.fnul) || (AZ > machine_consts.fnul) {
        //-----------------------------------------------------------------------
        //     INCREMENT FNU+NN-1 UP TO FNUL, COMPUTE AND RECUR BACKWARD
        //-----------------------------------------------------------------------
        let NUI_isize = (machine_consts.fnul - DFNU) as isize + 1;
        let NUI = NUI_isize.max(0) as usize;
        let ( NW, NLAST) = ZBUNI(
            //ZR, ZI, FNU,
            z,
            order,
            KODE,
            NN,
            NUI, //CYR, CYI, NW, NUI, NLAST, FNUL,
            &mut cy,
            machine_consts,
        )?;
        //    * TOL, ELIM, ALIM)?
        // if (NW < 0) GO TO 130
        NZ = NZ + NW;
        if NLAST == 0 {
            return Ok((cy, NZ));
        }
        NN = NLAST;
        // GO TO 60
    }
    //  60 CONTINUE
    // 'l60: loop{
    if !skip_az_rl_check & !(AZ > machine_consts.rl) {
        // GO TO 80
        //  70 CONTINUE
        //-----------------------------------------------------------------------
        //     MILLER ALGORITHM NORMALIZED BY THE SERIES
        //-----------------------------------------------------------------------
        let (cy, _) = i_miller(
            //ZR, ZI, FNU,
            z,
            order,
            KODE,
            NN,
            /*CYR, CYI, NW,*/ machine_consts,
        )?;
        //      return{if(NW < 0) { if(NW == (-2)) { Err(DidNotConverge);}// NZ=-2
        //     else{Err(Overflow)}}}else{
        return Ok((cy, NZ)); //}
    }
    //    80 CONTINUE
    //-----------------------------------------------------------------------
    //     MILLER ALGORITHM NORMALIZED BY THE WRONSKIAN
    //-----------------------------------------------------------------------
    //-----------------------------------------------------------------------
    //     OVERFLOW TEST ON K FUNCTIONS USED IN WRONSKIAN
    //-----------------------------------------------------------------------
    //       CALL ZUOIK(ZR, ZI, FNU, KODE, 2, 2, CWR, CWI, NW, TOL, ELIM,
    //      * ALIM)

    if let Ok( NW) = zuoik(
        z,
        order,
        KODE,
        IKType::K,
        2,
        &mut vec![c_one(); 2],
        machine_consts,
    ) {
        if NW > 0 {
            return Err(Overflow);
        } else {
            let nz= i_wronksian(z, order, KODE, NN, &mut cy, machine_consts)?;
            return Ok((cy, nz));
        }
    } else {
        return Ok((vec![c_one(); NN], NN));
    }

    /*
    // 110 CONTINUE

    //   120 CONTINUE
    //       RETURN
      130 CONTINUE
      //     NZ = -1
          return if(NW == (-2)) { Err(DidNotConverge);}// NZ=-2
          else{Err(Overflow)}

          RETURN
          // END
          */
}


fn ZACAI(//ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, RL, TOL,
     //* ELIM, ALIM)
     z: Complex64, order: f64, KODE: Scaling, MR: i32, N:usize, machine_consts: &MachineConsts)-> BesselResult{
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
// ***ROUTINES CALLED  ZASYI,ZBKNU,ZMLRI,z_power_series,ZS1S2,d1mach,ZABS
// ***END PROLOGUE  ZACAI
//     COMPLEX CSGN,CSPN,C1,C2,Y,Z,ZN,CY
    //   EXTERNAL ZABS
    //   DOUBLE PRECISION ALIM, ARG, ASCLE, AZ, CSGNR, CSGNI, CSPNR,
    //  * CSPNI, C1R, C1I, C2R, C2I, CYR, CYI, DFNU, ELIM, FMR, FNU, PI,
    //  * RL, SGN, TOL, YY, YR, YI, ZR, ZI, ZNR, ZNI, d1mach, ZABS
    //   INTEGER INU, IUF, KODE, MR, N, NN, NW, NZ
    //   DIMENSION YR(N), YI(N), CYR(2), CYI(2)
    //   DATA PI / 3.14159265358979324 /
      let mut NZ = 0;
      let zn = -z;
    //   ZNR = -ZR;
    //   ZNI = -ZI;
    let AZ = z.abs();
    //   AZ = ZABS(ZR,ZI);
      let NN = N;
      let DFNU = order + ((N-1) as f64);
    //   if (AZ <= 2.0) GO TO 10;
    let (mut y, _) = if (AZ*AZ*0.25 <= DFNU+1.0) || (AZ <= 2.0) {//GO TO 20;
//    10 CONTINUE;
//-----------------------------------------------------------------------;
//     POWER SERIES FOR THE I FUNCTION;
//-----------------------------------------------------------------------;
let (y, NW_signed) = i_power_series(zn, order, KODE, NN, machine_consts)?;
debug_assert!(NW_signed>=0);
(y, NW_signed.abs() as usize)
    //   CALL z_power_series(ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, TOL, ELIM, ALIM);
      }else if AZ >=machine_consts.rl{
    //   GO TO 40;
//    20 CONTINUE;
    //   if (AZ < RL) GO TO 30;
//-----------------------------------------------------------------------;
//     ASYMPTOTIC EXPANSION FOR LARGE Z FOR THE I FUNCTION;
//-----------------------------------------------------------------------;
    //   CALL ZASYI(ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, RL, TOL, ELIM,;
    //  * ALIM);
     z_asymptotic_i(zn, order, KODE,NN, machine_consts)?
    //   if (NW < 0) GO TO 80;
    //   GO TO 40;
//    30 CONTINUE;
//-----------------------------------------------------------------------;
//     MILLER ALGORITHM NORMALIZED BY THE SERIES FOR THE I FUNCTION;
//-----------------------------------------------------------------------;
      }else{
        i_miller(z, order, KODE, NN, machine_consts)?
    //   CALL ZMLRI(ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, TOL);
    //   if(NW < 0) GO TO 80;
     };
//    40 CONTINUE;
//-----------------------------------------------------------------------;
//     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION;
//-----------------------------------------------------------------------;
     let (cy, _) = ZBKNU(zn, order, KODE, 1, machine_consts)?;
    //   CALL ZBKNU(ZNR, ZNI, FNU, KODE, 1, CYR, CYI, NW, TOL, ELIM, ALIM);
    //   if (NW != 0) GO TO 80;
    //   let FMR = (MR as f64);s
      let SGN = -PI * (MR as f64).signum();//-DSIGN(PI,FMR);
      let mut csgn = Complex64::new(0.0, SGN);
    //   CSGNR = 0.0;
    //   CSGNI = SGN;
      if KODE == Scaling::Scaled {// GO TO 50;
        csgn = Complex64::I * csgn.im * Complex64::cis(-zn.im);
        // csgn.re = -csgn.re;
    //   YY = -ZNI;
    //   CSGNR = -CSGNI*DSIN(YY);
    //   CSGNI = CSGNI*DCOS(YY);
      }
//    50 CONTINUE;
//-----------------------------------------------------------------------;
//     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE;
//     WHEN FNU IS LARGE;
//-----------------------------------------------------------------------;
    //   INU = INT(SNGL(FNU));
      let INU = order as usize;
    //   ARG = (FNU-(INU as f64))*SGN;
    // let ARG = order.fract() * SGN;
      let mut cspn = Complex64::cis(order.fract()*SGN);
    //   CSPNR = DCOS(ARG);
    //   CSPNI = DSIN(ARG);
      if INU%2 != 0 {//GO TO 60;
        cspn = - cspn;
    //   CSPNR = -CSPNR;
    //   CSPNI = -CSPNI;
      }
//    60 CONTINUE;
      let mut c1 = cy[0];
      let mut c2 = y[0];
    //   C1R = CYR(1);
    //   C1I = CYI(1);
    //   C2R = YR(1);
    //   C2I = YI(1);
      if KODE == Scaling::Scaled {//GO TO 70;
      let mut IUF = 0;
    //   ASCLE = 1.0e+3*d1mach(1)/TOL;
        let NW = ZS1S2(zn, &mut c1, &mut c2, &mut IUF, machine_consts);
    //   CALL ZS1S2(ZNR, ZNI, C1R, C1I, C2R, C2I, NW, ASCLE, ALIM, IUF);
      NZ +=  NW;
      }
//    70 CONTINUE;
      y[0] = cspn*c1 + csgn*c2;
    //   YR(1) = CSPNR*C1R - CSPNI*C1I + CSGNR*C2R - CSGNI*C2I;
    //   YI(1) = CSPNR*C1I + CSPNI*C1R + CSGNR*C2I + CSGNI*C2R;
      Ok((y, NZ))
    //   RETURN;
//    80 CONTINUE;
//       NZ = -1;
//       if(NW == (-2)) NZ=-2;
//       RETURN;
//       END;
     }

/*
fn ZUNK1(ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, TOL, ELIM,
     * ALIM)
// ***BEGIN PROLOGUE  ZUNK1
// ***REFER TO  ZBESK
//
//     ZUNK1 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE
//     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE
//     UNIFORM ASYMPTOTIC EXPANSION.
//     MR INDICATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION.
//     NZ=-1 MEANS AN OVERFLOW WILL OCCUR
//
// ***ROUTINES CALLED  ZKSCL,ZS1S2,ZUunderflowCHK,ZUNIK,d1mach,ZABS
// ***END PROLOGUE  ZUNK1
//     COMPLEX CFN,CK,CONE,CRSC,CS,CSCL,CSGN,CSPN,CSR,CSS,CWRK,CY,CZERO,
//    *C1,C2,PHI,PHID,RZ,SUM,SUMD,S1,S2,Y,Z,ZETA1,ZETA1D,ZETA2,ZETA2D,ZR
      EXTERNAL ZABS
      DOUBLE PRECISION ALIM, ANG, APHI, ASC, ASCLE, BRY, CKI, CKR,
     * CONER, CRSC, CSCL, CSGNI, CSPNI, CSPNR, CSR, CSRR, CSSR,
     * CWRKI, CWRKR, CYI, CYR, C1I, C1R, C2I, C2M, C2R, ELIM, FMR, FN,
     * FNF, FNU, PHIDI, PHIDR, PHII, PHIR, PI, RAST, RAZR, RS1, RZI,
     * RZR, SGN, STI, STR, SUMDI, SUMDR, SUMI, SUMR, S1I, S1R, S2I,
     * S2R, TOL, YI, YR, ZEROI, ZEROR, ZETA1I, ZETA1R, ZETA2I, ZETA2R,
     * ZET1DI, ZET1DR, ZET2DI, ZET2DR, ZI, ZR, ZRI, ZRR, d1mach, ZABS
      INTEGER I, IB, IFLAG, IFN, IL, INIT, INU, IUF, K, KDFLG, KFLAG,
     * KK, KODE, MR, N, NW, NZ, INITD, IC, IPARD, J, M
      DIMENSION BRY(3), INIT(2), YR(N), YI(N), SUMR(2), SUMI(2),
     * ZETA1R(2), ZETA1I(2), ZETA2R(2), ZETA2I(2), CYR(2), CYI(2),
     * CWRKR(16,3), CWRKI(16,3), CSSR(3), CSRR(3), PHIR(2), PHII(2)
      DATA ZEROR,ZEROI,CONER / 0.0, 0.0, 1.0 /
      DATA PI / 3.14159265358979324 /
//
      KDFLG = 1
      NZ = 0
//-----------------------------------------------------------------------
//     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION GREATER THAN
//     THE UNDERFLOW LIMIT
//-----------------------------------------------------------------------
      CSCL = 1.0/TOL
      CRSC = TOL
      CSSR(1) = CSCL
      CSSR(2) = CONER
      CSSR(3) = CRSC
      CSRR(1) = CRSC
      CSRR(2) = CONER
      CSRR(3) = CSCL
      BRY(1) = 1.0e+3*d1mach(1)/TOL
      BRY(2) = 1.0/BRY(1)
      BRY(3) = d1mach(2)
      ZRR = ZR
      ZRI = ZI
      if (ZR >= 0.0) GO TO 10
      ZRR = -ZR
      ZRI = -ZI
   10 CONTINUE
      J = 2
      DO 70 I=1,N
//-----------------------------------------------------------------------
//     J FLIP FLOPS BETWEEN 1 AND 2 IN J = 3 - J
//-----------------------------------------------------------------------
        J = 3 - J
        FN = FNU + ((I-1) as f64)
        INIT(J) = 0
        CALL ZUNIK(ZRR, ZRI, FN, 2, 0, TOL, INIT(J), PHIR(J), PHII(J),
     *   ZETA1R(J), ZETA1I(J), ZETA2R(J), ZETA2I(J), SUMR(J), SUMI(J),
     *   CWRKR(1,J), CWRKI(1,J))
        if (KODE == 1) GO TO 20
        STR = ZRR + ZETA2R(J)
        STI = ZRI + ZETA2I(J)
        RAST = FN/ZABS(STR,STI)
        STR = STR*RAST*RAST
        STI = -STI*RAST*RAST
        S1R = ZETA1R(J) - STR
        S1I = ZETA1I(J) - STI
        GO TO 30
   20   CONTINUE
        S1R = ZETA1R(J) - ZETA2R(J)
        S1I = ZETA1I(J) - ZETA2I(J)
   30   CONTINUE
        RS1 = S1R
//-----------------------------------------------------------------------
//     TEST FOR UNDERFLOW AND OVERFLOW
//-----------------------------------------------------------------------
        if ((RS1).abs() > ELIM) GO TO 60
        if (KDFLG == 1) KFLAG = 2
        if ((RS1).abs() < ALIM) GO TO 40
//-----------------------------------------------------------------------
//     REFINE  TEST AND SCALE
//-----------------------------------------------------------------------
        APHI = ZABS(PHIR(J),PHII(J))
        RS1 = RS1 + DLOG(APHI)
        if ((RS1).abs() > ELIM) GO TO 60
        if (KDFLG == 1) KFLAG = 1
        if (RS1 < 0.0) GO TO 40
        if (KDFLG == 1) KFLAG = 3
   40   CONTINUE
//-----------------------------------------------------------------------
//     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
//     EXPONENT EXTREMES
//-----------------------------------------------------------------------
        S2R = PHIR(J)*SUMR(J) - PHII(J)*SUMI(J)
        S2I = PHIR(J)*SUMI(J) + PHII(J)*SUMR(J)
        STR = DEXP(S1R)*CSSR(KFLAG)
        S1R = STR*DCOS(S1I)
        S1I = STR*DSIN(S1I)
        STR = S2R*S1R - S2I*S1I
        S2I = S1R*S2I + S2R*S1I
        S2R = STR
        if (KFLAG != 1) GO TO 50
        CALL ZUunderflowCHK(S2R, S2I, NW, BRY(1), TOL)
        if (NW != 0) GO TO 60
   50   CONTINUE
        CYR(KDFLG) = S2R
        CYI(KDFLG) = S2I
        YR(I) = S2R*CSRR(KFLAG)
        YI(I) = S2I*CSRR(KFLAG)
        if (KDFLG == 2) GO TO 75
        KDFLG = 2
        GO TO 70
   60   CONTINUE
        if (RS1 > 0.0) GO TO 300
//-----------------------------------------------------------------------
//     FOR ZR < 0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
//-----------------------------------------------------------------------
        if (ZR < 0.0) GO TO 300
        KDFLG = 1
        YR(I)=ZEROR
        YI(I)=ZEROI
        NZ=NZ+1
        if (I == 1) GO TO 70
        if ((YR(I-1) == ZEROR)&&(YI(I-1) == ZEROI)) GO TO 70
        YR(I-1)=ZEROR
        YI(I-1)=ZEROI
        NZ=NZ+1
   70 CONTINUE
      I = N
   75 CONTINUE
      RAZR = 1.0/ZABS(ZRR,ZRI)
      STR = ZRR*RAZR
      STI = -ZRI*RAZR
      RZR = (STR+STR)*RAZR
      RZI = (STI+STI)*RAZR
      CKR = FN*RZR
      CKI = FN*RZI
      IB = I + 1
      if (N < IB) GO TO 160
//-----------------------------------------------------------------------
//     TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW. SET SEQUENCE TO ZERO
//     ON UNDERFLOW.
//-----------------------------------------------------------------------
      FN = FNU + ((N-1) as f64)
      IPARD = 1
      if (MR != 0) IPARD = 0
      INITD = 0
      CALL ZUNIK(ZRR, ZRI, FN, 2, IPARD, TOL, INITD, PHIDR, PHIDI,
     * ZET1DR, ZET1DI, ZET2DR, ZET2DI, SUMDR, SUMDI, CWRKR(1,3),
     * CWRKI(1,3))
      if (KODE == 1) GO TO 80
      STR = ZRR + ZET2DR
      STI = ZRI + ZET2DI
      RAST = FN/ZABS(STR,STI)
      STR = STR*RAST*RAST
      STI = -STI*RAST*RAST
      S1R = ZET1DR - STR
      S1I = ZET1DI - STI
      GO TO 90
   80 CONTINUE
      S1R = ZET1DR - ZET2DR
      S1I = ZET1DI - ZET2DI
   90 CONTINUE
      RS1 = S1R
      if ((RS1).abs() > ELIM) GO TO 95
      if ((RS1).abs() < ALIM) GO TO 100
//----------------------------------------------------------------------------
//     REFINE ESTIMATE AND TEST
//-------------------------------------------------------------------------
      APHI = ZABS(PHIDR,PHIDI)
      RS1 = RS1+DLOG(APHI)
      if ((RS1).abs() < ELIM) GO TO 100
   95 CONTINUE
      if ((RS1).abs() > 0.0) GO TO 300
//-----------------------------------------------------------------------
//     FOR ZR < 0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
//-----------------------------------------------------------------------
      if (ZR < 0.0) GO TO 300
      NZ = N
      DO 96 I=1,N
        YR(I) = ZEROR
        YI(I) = ZEROI
   96 CONTINUE
      RETURN
//---------------------------------------------------------------------------
//     FORWARD RECUR FOR REMAINDER OF THE SEQUENCE
//----------------------------------------------------------------------------
  100 CONTINUE
      S1R = CYR(1)
      S1I = CYI(1)
      S2R = CYR(2)
      S2I = CYI(2)
      C1R = CSRR(KFLAG)
      ASCLE = BRY(KFLAG)
      DO 120 I=IB,N
        C2R = S2R
        C2I = S2I
        S2R = CKR*C2R - CKI*C2I + S1R
        S2I = CKR*C2I + CKI*C2R + S1I
        S1R = C2R
        S1I = C2I
        CKR = CKR + RZR
        CKI = CKI + RZI
        C2R = S2R*C1R
        C2I = S2I*C1R
        YR(I) = C2R
        YI(I) = C2I
        if (KFLAG >= 3) GO TO 120
        STR = (C2R).abs()
        STI = (C2I).abs()
        C2M = DMAX1(STR,STI)
        if (C2M <= ASCLE) GO TO 120
        KFLAG = KFLAG + 1
        ASCLE = BRY(KFLAG)
        S1R = S1R*C1R
        S1I = S1I*C1R
        S2R = C2R
        S2I = C2I
        S1R = S1R*CSSR(KFLAG)
        S1I = S1I*CSSR(KFLAG)
        S2R = S2R*CSSR(KFLAG)
        S2I = S2I*CSSR(KFLAG)
        C1R = CSRR(KFLAG)
  120 CONTINUE
  160 CONTINUE
      if (MR == 0) RETURN
//-----------------------------------------------------------------------
//     ANALYTIC CONTINUATION FOR RE(Z) < 0.0
//-----------------------------------------------------------------------
      NZ = 0
      FMR = (MR as f64)
      SGN = -DSIGN(PI,FMR)
//-----------------------------------------------------------------------
//     CSPN AND CSGN ARE COEFF OF K AND I FUNCTIONS RESP.
//-----------------------------------------------------------------------
      CSGNI = SGN
      INU = INT(SNGL(FNU))
      FNF = FNU - (INU as f64)
      IFN = INU + N - 1
      ANG = FNF*SGN
      CSPNR = DCOS(ANG)
      CSPNI = DSIN(ANG)
      if (MOD(IFN,2) == 0) GO TO 170
      CSPNR = -CSPNR
      CSPNI = -CSPNI
  170 CONTINUE
      ASC = BRY(1)
      IUF = 0
      KK = N
      KDFLG = 1
      IB = IB - 1
      IC = IB - 1
      DO 270 K=1,N
        FN = FNU + ((KK-1) as f64)
//-----------------------------------------------------------------------
//     LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K
//     FUNCTION ABOVE
//-----------------------------------------------------------------------
        M=3
        if (N > 2) GO TO 175
  172   CONTINUE
        INITD = INIT(J)
        PHIDR = PHIR(J)
        PHIDI = PHII(J)
        ZET1DR = ZETA1R(J)
        ZET1DI = ZETA1I(J)
        ZET2DR = ZETA2R(J)
        ZET2DI = ZETA2I(J)
        SUMDR = SUMR(J)
        SUMDI = SUMI(J)
        M = J
        J = 3 - J
        GO TO 180
  175   CONTINUE
        if ((KK == N)&&(IB < N)) GO TO 180
        if ((KK == IB)||(KK == IC)) GO TO 172
        INITD = 0
  180   CONTINUE
        CALL ZUNIK(ZRR, ZRI, FN, 1, 0, TOL, INITD, PHIDR, PHIDI,
     *   ZET1DR, ZET1DI, ZET2DR, ZET2DI, SUMDR, SUMDI,
     *   CWRKR(1,M), CWRKI(1,M))
        if (KODE == 1) GO TO 200
        STR = ZRR + ZET2DR
        STI = ZRI + ZET2DI
        RAST = FN/ZABS(STR,STI)
        STR = STR*RAST*RAST
        STI = -STI*RAST*RAST
        S1R = -ZET1DR + STR
        S1I = -ZET1DI + STI
        GO TO 210
  200   CONTINUE
        S1R = -ZET1DR + ZET2DR
        S1I = -ZET1DI + ZET2DI
  210   CONTINUE
//-----------------------------------------------------------------------
//     TEST FOR UNDERFLOW AND OVERFLOW
//-----------------------------------------------------------------------
        RS1 = S1R
        if ((RS1).abs() > ELIM) GO TO 260
        if (KDFLG == 1) IFLAG = 2
        if ((RS1).abs() < ALIM) GO TO 220
//-----------------------------------------------------------------------
//     REFINE  TEST AND SCALE
//-----------------------------------------------------------------------
        APHI = ZABS(PHIDR,PHIDI)
        RS1 = RS1 + DLOG(APHI)
        if ((RS1).abs() > ELIM) GO TO 260
        if (KDFLG == 1) IFLAG = 1
        if (RS1 < 0.0) GO TO 220
        if (KDFLG == 1) IFLAG = 3
  220   CONTINUE
        STR = PHIDR*SUMDR - PHIDI*SUMDI
        STI = PHIDR*SUMDI + PHIDI*SUMDR
        S2R = -CSGNI*STI
        S2I = CSGNI*STR
        STR = DEXP(S1R)*CSSR(IFLAG)
        S1R = STR*DCOS(S1I)
        S1I = STR*DSIN(S1I)
        STR = S2R*S1R - S2I*S1I
        S2I = S2R*S1I + S2I*S1R
        S2R = STR
        if (IFLAG != 1) GO TO 230
        CALL ZUunderflowCHK(S2R, S2I, NW, BRY(1), TOL)
        if (NW == 0) GO TO 230
        S2R = ZEROR
        S2I = ZEROI
  230   CONTINUE
        CYR(KDFLG) = S2R
        CYI(KDFLG) = S2I
        C2R = S2R
        C2I = S2I
        S2R = S2R*CSRR(IFLAG)
        S2I = S2I*CSRR(IFLAG)
//-----------------------------------------------------------------------
//     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N
//-----------------------------------------------------------------------
        S1R = YR(KK)
        S1I = YI(KK)
        if (KODE == 1) GO TO 250
        CALL ZS1S2(ZRR, ZRI, S1R, S1I, S2R, S2I, NW, ASC, ALIM, IUF)
        NZ = NZ + NW
  250   CONTINUE
        YR(KK) = S1R*CSPNR - S1I*CSPNI + S2R
        YI(KK) = CSPNR*S1I + CSPNI*S1R + S2I
        KK = KK - 1
        CSPNR = -CSPNR
        CSPNI = -CSPNI
        if (C2R != 0.0 || C2I != 0.0) GO TO 255
        KDFLG = 1
        GO TO 270
  255   CONTINUE
        if (KDFLG == 2) GO TO 275
        KDFLG = 2
        GO TO 270
  260   CONTINUE
        if (RS1 > 0.0) GO TO 300
        S2R = ZEROR
        S2I = ZEROI
        GO TO 230
  270 CONTINUE
      K = N
  275 CONTINUE
      IL = N - K
      if (IL == 0) RETURN
//-----------------------------------------------------------------------
//     RECUR BACKWARD FOR REMAINDER OF I SEQUENCE AND ADD IN THE
//     K FUNCTIONS, SCALING THE I SEQUENCE DURING RECURRENCE TO KEEP
//     INTERMEDIATE ARITHMETIC ON SCALE NEAR EXPONENT EXTREMES.
//-----------------------------------------------------------------------
      S1R = CYR(1)
      S1I = CYI(1)
      S2R = CYR(2)
      S2I = CYI(2)
      CSR = CSRR(IFLAG)
      ASCLE = BRY(IFLAG)
      FN = DBLE(FLOAT(INU+IL))
      DO 290 I=1,IL
        C2R = S2R
        C2I = S2I
        S2R = S1R + (FN+FNF)*(RZR*C2R-RZI*C2I)
        S2I = S1I + (FN+FNF)*(RZR*C2I+RZI*C2R)
        S1R = C2R
        S1I = C2I
        FN = FN - 1.0
        C2R = S2R*CSR
        C2I = S2I*CSR
        CKR = C2R
        CKI = C2I
        C1R = YR(KK)
        C1I = YI(KK)
        if (KODE == 1) GO TO 280
        CALL ZS1S2(ZRR, ZRI, C1R, C1I, C2R, C2I, NW, ASC, ALIM, IUF)
        NZ = NZ + NW
  280   CONTINUE
        YR(KK) = C1R*CSPNR - C1I*CSPNI + C2R
        YI(KK) = C1R*CSPNI + C1I*CSPNR + C2I
        KK = KK - 1
        CSPNR = -CSPNR
        CSPNI = -CSPNI
        if (IFLAG >= 3) GO TO 290
        C2R = (CKR).abs()
        C2I = (CKI).abs()
        C2M = DMAX1(C2R,C2I)
        if (C2M <= ASCLE) GO TO 290
        IFLAG = IFLAG + 1
        ASCLE = BRY(IFLAG)
        S1R = S1R*CSR
        S1I = S1I*CSR
        S2R = CKR
        S2I = CKI
        S1R = S1R*CSSR(IFLAG)
        S1I = S1I*CSSR(IFLAG)
        S2R = S2R*CSSR(IFLAG)
        S2I = S2I*CSSR(IFLAG)
        CSR = CSRR(IFLAG)
  290 CONTINUE
      RETURN
  300 CONTINUE
      NZ = -1
      RETURN
      END
fn ZUNK2(ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, TOL, ELIM,
     * ALIM)
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
      EXTERNAL ZABS
      DOUBLE PRECISION AARG, AIC, AII, AIR, ALIM, ANG, APHI, ARGDI,
     * ARGDR, ARGI, ARGR, ASC, ASCLE, ASUMDI, ASUMDR, ASUMI, ASUMR,
     * BRY, BSUMDI, BSUMDR, BSUMI, BSUMR, CAR, CIPI, CIPR, CKI, CKR,
     * CONER, CRSC, CR1I, CR1R, CR2I, CR2R, CSCL, CSGNI, CSI,
     * CSPNI, CSPNR, CSR, CSRR, CSSR, CYI, CYR, C1I, C1R, C2I, C2M,
     * C2R, DAII, DAIR, ELIM, FMR, FN, FNF, FNU, FRAC_PI_2, PHIDI, PHIDR,
     * PHII, PHIR, PI, PTI, PTR, RAST, RAZR, RS1, RZI, RZR, SAR, SGN,
     * STI, STR, S1I, S1R, S2I, S2R, TOL, YI, YR, YY, ZBI, ZBR, ZEROI,
     * ZEROR, ZETA1I, ZETA1R, ZETA2I, ZETA2R, ZET1DI, ZET1DR, ZET2DI,
     * ZET2DR, ZI, ZNI, ZNR, ZR, ZRI, ZRR, d1mach, ZABS
      INTEGER I, IB, IFLAG, IFN, IL, IN, INU, IUF, K, KDFLG, KFLAG, KK,
     * KODE, MR, N, NAI, NDAI, NW, NZ, IDUM, J, IPARD, IC
      DIMENSION BRY(3), YR(N), YI(N), ASUMR(2), ASUMI(2), BSUMR(2),
     * BSUMI(2), PHIR(2), PHII(2), ARGR(2), ARGI(2), ZETA1R(2),
     * ZETA1I(2), ZETA2R(2), ZETA2I(2), CYR(2), CYI(2), CIPR(4),
     * CIPI(4), CSSR(3), CSRR(3)
      DATA ZEROR,ZEROI,CONER,CR1R,CR1I,CR2R,CR2I /
     1         0.0, 0.0, 1.0,
     1 1.0,1.73205080756887729 , -0.5,-8.66025403784438647e-01 /
      DATA FRAC_PI_2, PI, AIC /
     1     1.57079632679489662e+00,     3.14159265358979324e+00,
     1     1.26551212348464539e+00/
      DATA CIPR(1),CIPI(1),CIPR(2),CIPI(2),CIPR(3),CIPI(3),CIPR(4),
     * CIPI(4) /
     1  1.0,0.0 ,  0.0,-1.0 ,  -1.0,0.0 ,  0.0,1.0 /
//
      KDFLG = 1
      NZ = 0
//-----------------------------------------------------------------------
//     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION GREATER THAN
//     THE UNDERFLOW LIMIT
//-----------------------------------------------------------------------
      CSCL = 1.0/TOL
      CRSC = TOL
      CSSR(1) = CSCL
      CSSR(2) = CONER
      CSSR(3) = CRSC
      CSRR(1) = CRSC
      CSRR(2) = CONER
      CSRR(3) = CSCL
      BRY(1) = 1.0e+3*d1mach(1)/TOL
      BRY(2) = 1.0/BRY(1)
      BRY(3) = d1mach(2)
      ZRR = ZR
      ZRI = ZI
      if (ZR >= 0.0) GO TO 10
      ZRR = -ZR
      ZRI = -ZI
   10 CONTINUE
      YY = ZRI
      ZNR = ZRI
      ZNI = -ZRR
      ZBR = ZRR
      ZBI = ZRI
      INU = INT(SNGL(FNU))
      FNF = FNU - (INU as f64)
      ANG = -FRAC_PI_2*FNF
      CAR = DCOS(ANG)
      SAR = DSIN(ANG)
      C2R = FRAC_PI_2*SAR
      C2I = -FRAC_PI_2*CAR
      KK = MOD(INU,4) + 1
      STR = C2R*CIPR(KK) - C2I*CIPI(KK)
      STI = C2R*CIPI(KK) + C2I*CIPR(KK)
      CSR = CR1R*STR - CR1I*STI
      CSI = CR1R*STI + CR1I*STR
      if (YY > 0.0) GO TO 20
      ZNR = -ZNR
      ZBI = -ZBI
   20 CONTINUE
//-----------------------------------------------------------------------
//     K(FNU,Z) IS COMPUTED FROM H(2,FNU,-I*Z) WHERE Z IS IN THE FIRST
//     QUADRANT. FOURTH QUADRANT VALUES (YY <= 0.0E0) ARE COMPUTED BY
//     CONJUGATION SINCE THE K FUNCTION IS REAL ON THE POSITIVE REAL AXIS
//-----------------------------------------------------------------------
      J = 2
      DO 80 I=1,N
//-----------------------------------------------------------------------
//     J FLIP FLOPS BETWEEN 1 AND 2 IN J = 3 - J
//-----------------------------------------------------------------------
        J = 3 - J
        FN = FNU + ((I-1) as f64)
        CALL ZUNHJ(ZNR, ZNI, FN, 0, TOL, PHIR(J), PHII(J), ARGR(J),
     *   ARGI(J), ZETA1R(J), ZETA1I(J), ZETA2R(J), ZETA2I(J), ASUMR(J),
     *   ASUMI(J), BSUMR(J), BSUMI(J))
        if (KODE == 1) GO TO 30
        STR = ZBR + ZETA2R(J)
        STI = ZBI + ZETA2I(J)
        RAST = FN/ZABS(STR,STI)
        STR = STR*RAST*RAST
        STI = -STI*RAST*RAST
        S1R = ZETA1R(J) - STR
        S1I = ZETA1I(J) - STI
        GO TO 40
   30   CONTINUE
        S1R = ZETA1R(J) - ZETA2R(J)
        S1I = ZETA1I(J) - ZETA2I(J)
   40   CONTINUE
//-----------------------------------------------------------------------
//     TEST FOR UNDERFLOW AND OVERFLOW
//-----------------------------------------------------------------------
        RS1 = S1R
        if ((RS1).abs() > ELIM) GO TO 70
        if (KDFLG == 1) KFLAG = 2
        if ((RS1).abs() < ALIM) GO TO 50
//-----------------------------------------------------------------------
//     REFINE  TEST AND SCALE
//-----------------------------------------------------------------------
        APHI = ZABS(PHIR(J),PHII(J))
        AARG = ZABS(ARGR(J),ARGI(J))
        RS1 = RS1 + DLOG(APHI) - 0.25*DLOG(AARG) - AIC
        if ((RS1).abs() > ELIM) GO TO 70
        if (KDFLG == 1) KFLAG = 1
        if (RS1 < 0.0) GO TO 50
        if (KDFLG == 1) KFLAG = 3
   50   CONTINUE
//-----------------------------------------------------------------------
//     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
//     EXPONENT EXTREMES
//-----------------------------------------------------------------------
        C2R = ARGR(J)*CR2R - ARGI(J)*CR2I
        C2I = ARGR(J)*CR2I + ARGI(J)*CR2R
        CALL ZAIRY(C2R, C2I, 0, 2, AIR, AII, NAI, IDUM)
        CALL ZAIRY(C2R, C2I, 1, 2, DAIR, DAII, NDAI, IDUM)
        STR = DAIR*BSUMR(J) - DAII*BSUMI(J)
        STI = DAIR*BSUMI(J) + DAII*BSUMR(J)
        PTR = STR*CR2R - STI*CR2I
        PTI = STR*CR2I + STI*CR2R
        STR = PTR + (AIR*ASUMR(J)-AII*ASUMI(J))
        STI = PTI + (AIR*ASUMI(J)+AII*ASUMR(J))
        PTR = STR*PHIR(J) - STI*PHII(J)
        PTI = STR*PHII(J) + STI*PHIR(J)
        S2R = PTR*CSR - PTI*CSI
        S2I = PTR*CSI + PTI*CSR
        STR = DEXP(S1R)*CSSR(KFLAG)
        S1R = STR*DCOS(S1I)
        S1I = STR*DSIN(S1I)
        STR = S2R*S1R - S2I*S1I
        S2I = S1R*S2I + S2R*S1I
        S2R = STR
        if (KFLAG != 1) GO TO 60
        CALL ZUunderflowCHK(S2R, S2I, NW, BRY(1), TOL)
        if (NW != 0) GO TO 70
   60   CONTINUE
        if (YY <= 0.0) S2I = -S2I
        CYR(KDFLG) = S2R
        CYI(KDFLG) = S2I
        YR(I) = S2R*CSRR(KFLAG)
        YI(I) = S2I*CSRR(KFLAG)
        STR = CSI
        CSI = -CSR
        CSR = STR
        if (KDFLG == 2) GO TO 85
        KDFLG = 2
        GO TO 80
   70   CONTINUE
        if (RS1 > 0.0) GO TO 320
//-----------------------------------------------------------------------
//     FOR ZR < 0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
//-----------------------------------------------------------------------
        if (ZR < 0.0) GO TO 320
        KDFLG = 1
        YR(I)=ZEROR
        YI(I)=ZEROI
        NZ=NZ+1
        STR = CSI
        CSI =-CSR
        CSR = STR
        if (I == 1) GO TO 80
        if ((YR(I-1) == ZEROR)&&(YI(I-1) == ZEROI)) GO TO 80
        YR(I-1)=ZEROR
        YI(I-1)=ZEROI
        NZ=NZ+1
   80 CONTINUE
      I = N
   85 CONTINUE
      RAZR = 1.0/ZABS(ZRR,ZRI)
      STR = ZRR*RAZR
      STI = -ZRI*RAZR
      RZR = (STR+STR)*RAZR
      RZI = (STI+STI)*RAZR
      CKR = FN*RZR
      CKI = FN*RZI
      IB = I + 1
      if (N < IB) GO TO 180
//-----------------------------------------------------------------------
//     TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW. SET SEQUENCE TO ZERO
//     ON UNDERFLOW.
//-----------------------------------------------------------------------
      FN = FNU + ((N-1) as f64)
      IPARD = 1
      if (MR != 0) IPARD = 0
      CALL ZUNHJ(ZNR, ZNI, FN, IPARD, TOL, PHIDR, PHIDI, ARGDR, ARGDI,
     * ZET1DR, ZET1DI, ZET2DR, ZET2DI, ASUMDR, ASUMDI, BSUMDR, BSUMDI)
      if (KODE == 1) GO TO 90
      STR = ZBR + ZET2DR
      STI = ZBI + ZET2DI
      RAST = FN/ZABS(STR,STI)
      STR = STR*RAST*RAST
      STI = -STI*RAST*RAST
      S1R = ZET1DR - STR
      S1I = ZET1DI - STI
      GO TO 100
   90 CONTINUE
      S1R = ZET1DR - ZET2DR
      S1I = ZET1DI - ZET2DI
  100 CONTINUE
      RS1 = S1R
      if ((RS1).abs() > ELIM) GO TO 105
      if ((RS1).abs() < ALIM) GO TO 120
//----------------------------------------------------------------------------
//     REFINE ESTIMATE AND TEST
//-------------------------------------------------------------------------
      APHI = ZABS(PHIDR,PHIDI)
      RS1 = RS1+DLOG(APHI)
      if ((RS1).abs() < ELIM) GO TO 120
  105 CONTINUE
      if (RS1 > 0.0) GO TO 320
//-----------------------------------------------------------------------
//     FOR ZR < 0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
//-----------------------------------------------------------------------
      if (ZR < 0.0) GO TO 320
      NZ = N
      DO 106 I=1,N
        YR(I) = ZEROR
        YI(I) = ZEROI
  106 CONTINUE
      RETURN
  120 CONTINUE
      S1R = CYR(1)
      S1I = CYI(1)
      S2R = CYR(2)
      S2I = CYI(2)
      C1R = CSRR(KFLAG)
      ASCLE = BRY(KFLAG)
      DO 130 I=IB,N
        C2R = S2R
        C2I = S2I
        S2R = CKR*C2R - CKI*C2I + S1R
        S2I = CKR*C2I + CKI*C2R + S1I
        S1R = C2R
        S1I = C2I
        CKR = CKR + RZR
        CKI = CKI + RZI
        C2R = S2R*C1R
        C2I = S2I*C1R
        YR(I) = C2R
        YI(I) = C2I
        if (KFLAG >= 3) GO TO 130
        STR = (C2R).abs()
        STI = (C2I).abs()
        C2M = DMAX1(STR,STI)
        if (C2M <= ASCLE) GO TO 130
        KFLAG = KFLAG + 1
        ASCLE = BRY(KFLAG)
        S1R = S1R*C1R
        S1I = S1I*C1R
        S2R = C2R
        S2I = C2I
        S1R = S1R*CSSR(KFLAG)
        S1I = S1I*CSSR(KFLAG)
        S2R = S2R*CSSR(KFLAG)
        S2I = S2I*CSSR(KFLAG)
        C1R = CSRR(KFLAG)
  130 CONTINUE
  180 CONTINUE
      if (MR == 0) RETURN
//-----------------------------------------------------------------------
//     ANALYTIC CONTINUATION FOR RE(Z) < 0.0
//-----------------------------------------------------------------------
      NZ = 0
      FMR = (MR as f64)
      SGN = -DSIGN(PI,FMR)
//-----------------------------------------------------------------------
//     CSPN AND CSGN ARE COEFF OF K AND I FUNCIONS RESP.
//-----------------------------------------------------------------------
      CSGNI = SGN
      if (YY <= 0.0) CSGNI = -CSGNI
      IFN = INU + N - 1
      ANG = FNF*SGN
      CSPNR = DCOS(ANG)
      CSPNI = DSIN(ANG)
      if (MOD(IFN,2) == 0) GO TO 190
      CSPNR = -CSPNR
      CSPNI = -CSPNI
  190 CONTINUE
//-----------------------------------------------------------------------
//     CS=COEFF OF THE J FUNCTION TO GET THE I FUNCTION. I(FNU,Z) IS
//     COMPUTED FROM EXP(I*FNU*FRAC_PI_2)*J(FNU,-I*Z) WHERE Z IS IN THE FIRST
//     QUADRANT. FOURTH QUADRANT VALUES (YY <= 0.0E0) ARE COMPUTED BY
//     CONJUGATION SINCE THE I FUNCTION IS REAL ON THE POSITIVE REAL AXIS
//-----------------------------------------------------------------------
      CSR = SAR*CSGNI
      CSI = CAR*CSGNI
      IN = MOD(IFN,4) + 1
      C2R = CIPR(IN)
      C2I = CIPI(IN)
      STR = CSR*C2R + CSI*C2I
      CSI = -CSR*C2I + CSI*C2R
      CSR = STR
      ASC = BRY(1)
      IUF = 0
      KK = N
      KDFLG = 1
      IB = IB - 1
      IC = IB - 1
      DO 290 K=1,N
        FN = FNU + ((KK-1) as f64)
//-----------------------------------------------------------------------
//     LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K
//     FUNCTION ABOVE
//-----------------------------------------------------------------------
        if (N > 2) GO TO 175
  172   CONTINUE
        PHIDR = PHIR(J)
        PHIDI = PHII(J)
        ARGDR = ARGR(J)
        ARGDI = ARGI(J)
        ZET1DR = ZETA1R(J)
        ZET1DI = ZETA1I(J)
        ZET2DR = ZETA2R(J)
        ZET2DI = ZETA2I(J)
        ASUMDR = ASUMR(J)
        ASUMDI = ASUMI(J)
        BSUMDR = BSUMR(J)
        BSUMDI = BSUMI(J)
        J = 3 - J
        GO TO 210
  175   CONTINUE
        if ((KK == N)&&(IB < N)) GO TO 210
        if ((KK == IB)||(KK == IC)) GO TO 172
        CALL ZUNHJ(ZNR, ZNI, FN, 0, TOL, PHIDR, PHIDI, ARGDR,
     *   ARGDI, ZET1DR, ZET1DI, ZET2DR, ZET2DI, ASUMDR,
     *   ASUMDI, BSUMDR, BSUMDI)
  210   CONTINUE
        if (KODE == 1) GO TO 220
        STR = ZBR + ZET2DR
        STI = ZBI + ZET2DI
        RAST = FN/ZABS(STR,STI)
        STR = STR*RAST*RAST
        STI = -STI*RAST*RAST
        S1R = -ZET1DR + STR
        S1I = -ZET1DI + STI
        GO TO 230
  220   CONTINUE
        S1R = -ZET1DR + ZET2DR
        S1I = -ZET1DI + ZET2DI
  230   CONTINUE
//-----------------------------------------------------------------------
//     TEST FOR UNDERFLOW AND OVERFLOW
//-----------------------------------------------------------------------
        RS1 = S1R
        if ((RS1).abs() > ELIM) GO TO 280
        if (KDFLG == 1) IFLAG = 2
        if ((RS1).abs() < ALIM) GO TO 240
//-----------------------------------------------------------------------
//     REFINE  TEST AND SCALE
//-----------------------------------------------------------------------
        APHI = ZABS(PHIDR,PHIDI)
        AARG = ZABS(ARGDR,ARGDI)
        RS1 = RS1 + DLOG(APHI) - 0.25*DLOG(AARG) - AIC
        if ((RS1).abs() > ELIM) GO TO 280
        if (KDFLG == 1) IFLAG = 1
        if (RS1 < 0.0) GO TO 240
        if (KDFLG == 1) IFLAG = 3
  240   CONTINUE
        CALL ZAIRY(ARGDR, ARGDI, 0, 2, AIR, AII, NAI, IDUM)
        CALL ZAIRY(ARGDR, ARGDI, 1, 2, DAIR, DAII, NDAI, IDUM)
        STR = DAIR*BSUMDR - DAII*BSUMDI
        STI = DAIR*BSUMDI + DAII*BSUMDR
        STR = STR + (AIR*ASUMDR-AII*ASUMDI)
        STI = STI + (AIR*ASUMDI+AII*ASUMDR)
        PTR = STR*PHIDR - STI*PHIDI
        PTI = STR*PHIDI + STI*PHIDR
        S2R = PTR*CSR - PTI*CSI
        S2I = PTR*CSI + PTI*CSR
        STR = DEXP(S1R)*CSSR(IFLAG)
        S1R = STR*DCOS(S1I)
        S1I = STR*DSIN(S1I)
        STR = S2R*S1R - S2I*S1I
        S2I = S2R*S1I + S2I*S1R
        S2R = STR
        if (IFLAG != 1) GO TO 250
        CALL ZUunderflowCHK(S2R, S2I, NW, BRY(1), TOL)
        if (NW == 0) GO TO 250
        S2R = ZEROR
        S2I = ZEROI
  250   CONTINUE
        if (YY <= 0.0) S2I = -S2I
        CYR(KDFLG) = S2R
        CYI(KDFLG) = S2I
        C2R = S2R
        C2I = S2I
        S2R = S2R*CSRR(IFLAG)
        S2I = S2I*CSRR(IFLAG)
//-----------------------------------------------------------------------
//     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N
//-----------------------------------------------------------------------
        S1R = YR(KK)
        S1I = YI(KK)
        if (KODE == 1) GO TO 270
        CALL ZS1S2(ZRR, ZRI, S1R, S1I, S2R, S2I, NW, ASC, ALIM, IUF)
        NZ = NZ + NW
  270   CONTINUE
        YR(KK) = S1R*CSPNR - S1I*CSPNI + S2R
        YI(KK) = S1R*CSPNI + S1I*CSPNR + S2I
        KK = KK - 1
        CSPNR = -CSPNR
        CSPNI = -CSPNI
        STR = CSI
        CSI = -CSR
        CSR = STR
        if (C2R != 0.0 || C2I != 0.0) GO TO 255
        KDFLG = 1
        GO TO 290
  255   CONTINUE
        if (KDFLG == 2) GO TO 295
        KDFLG = 2
        GO TO 290
  280   CONTINUE
        if (RS1 > 0.0) GO TO 320
        S2R = ZEROR
        S2I = ZEROI
        GO TO 250
  290 CONTINUE
      K = N
  295 CONTINUE
      IL = N - K
      if (IL == 0) RETURN
//-----------------------------------------------------------------------
//     RECUR BACKWARD FOR REMAINDER OF I SEQUENCE AND ADD IN THE
//     K FUNCTIONS, SCALING THE I SEQUENCE DURING RECURRENCE TO KEEP
//     INTERMEDIATE ARITHMETIC ON SCALE NEAR EXPONENT EXTREMES.
//-----------------------------------------------------------------------
      S1R = CYR(1)
      S1I = CYI(1)
      S2R = CYR(2)
      S2I = CYI(2)
      CSR = CSRR(IFLAG)
      ASCLE = BRY(IFLAG)
      FN = DBLE(FLOAT(INU+IL))
      DO 310 I=1,IL
        C2R = S2R
        C2I = S2I
        S2R = S1R + (FN+FNF)*(RZR*C2R-RZI*C2I)
        S2I = S1I + (FN+FNF)*(RZR*C2I+RZI*C2R)
        S1R = C2R
        S1I = C2I
        FN = FN - 1.0
        C2R = S2R*CSR
        C2I = S2I*CSR
        CKR = C2R
        CKI = C2I
        C1R = YR(KK)
        C1I = YI(KK)
        if (KODE == 1) GO TO 300
        CALL ZS1S2(ZRR, ZRI, C1R, C1I, C2R, C2I, NW, ASC, ALIM, IUF)
        NZ = NZ + NW
  300   CONTINUE
        YR(KK) = C1R*CSPNR - C1I*CSPNI + C2R
        YI(KK) = C1R*CSPNI + C1I*CSPNR + C2I
        KK = KK - 1
        CSPNR = -CSPNR
        CSPNI = -CSPNI
        if (IFLAG >= 3) GO TO 310
        C2R = (CKR).abs()
        C2I = (CKI).abs()
        C2M = DMAX1(C2R,C2I)
        if (C2M <= ASCLE) GO TO 310
        IFLAG = IFLAG + 1
        ASCLE = BRY(IFLAG)
        S1R = S1R*CSR
        S1I = S1I*CSR
        S2R = CKR
        S2I = CKI
        S1R = S1R*CSSR(IFLAG)
        S1I = S1I*CSSR(IFLAG)
        S2R = S2R*CSSR(IFLAG)
        S2I = S2I*CSSR(IFLAG)
        CSR = CSRR(IFLAG)
  310 CONTINUE
      RETURN
  320 CONTINUE
      NZ = -1
      RETURN
      END
      */
fn ZBUNI(
    //ZR, ZI, FNU, KODE, N,
    z: Complex64,
    order: f64,
    KODE: Scaling,
    N: usize,
    NUI: usize,
    y: &mut Vec<Complex64>,
    // YR, YI, NZ, NUI, NLAST,
    machine_consts: &MachineConsts,
) -> Result<(usize, usize), BesselError> {
    // * FNUL, TOL, ELIM, ALIM)
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
    //     COMPLEX CSCL,CSCR,CY,RZ,ST,S1,S2,Y,Z
    //   EXTERNAL ZABS
    //   DOUBLE PRECISION ALIM, AX, AY, CSCLR, CSCRR, CYI, CYR, DFNU,
    //  * ELIM, FNU, FNUI, FNUL, GNU, RAZ, RZI, RZR, STI, STR, S1I, S1R,
    //  * S2I, S2R, TOL, YI, YR, ZI, ZR, ZABS, ASCLE, BRY, C1R, C1I, C1M,
    //  * d1mach
    //   INTEGER I, IFLAG, IFORM, K, KODE, N, NL, NLAST, NUI, NW, NZ
    //   DIMENSION YR(N), YI(N), CYR(2), CYI(2), BRY(3)
    let mut NZ = 0;
    let AX = z.re.abs() * 1.7321;
    let AY = z.im.abs();
    //   let mut IFORM = 1;
    let IFORM = if AY > AX { 2 } else { 1 };
    if NUI != 0 {
        //GO TO 60;
        let mut FNUI = NUI as f64;
        let DFNU = order + ((N - 1) as f64);
        let GNU = DFNU + FNUI;
        let mut cy = c_zeros(2);
        let ( NW, NLAST) = if IFORM != 2 {
            //GO TO 10;
            //-----------------------------------------------------------------------;
            //     ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN;
            //     -PI/3 <= ARG(Z) <= PI/3;
            //-----------------------------------------------------------------------;
            ZUNI1(
                //ZR, ZI, GNU, KODE, 2,
                z,
                GNU,
                KODE,
                2,
                &mut cy,
                machine_consts,
            )?
            //CYR, CYI, NW, NLAST, machine_consts)
            //FNUL, TOL, ELIM, ALIM);
        } else {
            //   GO TO 20;
            //    10 CONTINUE;
            //-----------------------------------------------------------------------;
            //     ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*FRAC_PI_2)) FOR LARGE FNU;
            //     APPLIED IN PI/3 < ABS(ARG(Z)) <= PI/2 WHERE M=+I OR -I;
            //     AND FRAC_PI_2=PI/2;
            //-----------------------------------------------------------------------;
            ZUNI2(
                //ZR, ZI, GNU, KODE, 2,
                z,
                GNU,
                KODE,
                2,
                &mut cy,
                machine_consts,
            )?
            //CYR, CYI, NW, NLAST, machine_consts)
            //FNUL, TOL, ELIM, ALIM);
        };
        //    20 CONTINUE;
        //   if (NW < 0) GO TO 50;
        if NW != 0 {
            return Ok((NZ, N /*NLAST=N*/));
        } //GO TO 90;
        //   STR = ZABS(CYR(1),CYI(1));
        //----------------------------------------------------------------------;
        //     SCALE BACKWARD RECURRENCE, BRY(3) IS DEFINED BUT NEVER USED;
        //----------------------------------------------------------------------;
        //   BRY(1)= machine_consts.ascle//1.0e+3*d1mach(1)/TOL;
        //   BRY(2) = 1.0/BRY(1);
        //   BRY(3) = BRY(2);
        let BRY = [
            machine_consts.ascle,
            1.0 / machine_consts.ascle,
            1.0 / machine_consts.ascle,
        ];
        let (mut IFLAG, mut ASCLE, mut CSCLR) = if cy[0].abs() <= BRY[0] {
            // GO TO 21;
            (1, BRY[0], 1.0 / machine_consts.tol)
        //   IFLAG = 1;
        //   ASCLE = BRY(1);
        //   CSCLR = 1.0/TOL;
        } else if cy[0].abs() >= BRY[1] {
            //   GO TO 25;
            //    21 CONTINUE{// GO TO 25;
            (3, BRY[2], machine_consts.tol)
        //   IFLAG = 3;
        //   ASCLE=BRY(3);
        //   CSCLR = TOL;
        } else {
            (2, BRY[1], 1.0)
        };

        //    25 CONTINUE;
        let mut CSCRR = 1.0 / CSCLR;
        let mut s1 = cy[1] * CSCLR;
        let mut s2 = cy[0] * CSCLR;
        //   S1R = CYR(2)*CSCLR;
        //   S1I = CYI(2)*CSCLR;
        //   S2R = CYR(1)*CSCLR;
        //   S2I = CYI(1)*CSCLR;
        let RAZ = 1.0 / z.abs(); //ZABS(ZR,ZI);
        //   STR = ZR*RAZ;
        //   STI = -ZI*RAZ;
        //   RZR = (STR+STR)*RAZ;
        //   RZI = (STI+STI)*RAZ;
        let rz = 2.0 * z.conj() * RAZ.pow(2);
        //   DO 30 I=1,NUI;
        for _ in 0..NUI {
            let st = s2;
            // STR = S2R;
            // STI = S2I;
            s2 = (DFNU + FNUI) * rz * s2 + s1;
            // S2R = (DFNU+FNUI)*(RZR*STR-RZI*STI) + S1R;
            // S2I = (DFNU+FNUI)*(RZR*STI+RZI*STR) + S1I;
            s1 = st;
            // S1R = STR;
            // S1I = STI;
            FNUI -= 1.0;
            if IFLAG >= 3 {
                continue;
            } //GO TO 30;
            // let c1 = (s2*CSCRR).abs();
            // STR = S2R*CSCRR;
            // STI = S2I*CSCRR;
            // C1R = (STR).abs();
            // C1I = (STI).abs();
            // let c1m = max_abs_component(s2*CSCRR);
            // C1M = DMAX1(C1R,C1I);
            let st = s2 * CSCRR;
            if max_abs_component(st) <= ASCLE {
                continue;
            } //GO TO 30;
            IFLAG += 1;
            ASCLE = BRY[IFLAG - 1];
            s1 *= CSCRR;
            // S1R = S1R*CSCRR;
            // S1I = S1I*CSCRR;
            // S2R = STR;
            // S2I = STI;
            s2 = st;
            CSCLR *= machine_consts.tol;
            CSCRR = 1.0 / CSCLR;
            s1 *= CSCLR;
            s2 *= CSCLR;
            // S1R = S1R*CSCLR;
            // S1I = S1I*CSCLR;
            // S2R = S2R*CSCLR;
            // S2I = S2I*CSCLR;
        }
        //    30 CONTINUE;
        y[N - 1] = s2 * CSCRR;
        //   YR(N) = S2R*CSCRR;
        //   YI(N) = S2I*CSCRR;
        if N == 1 {
            return Ok(( NZ, NLAST));
        } //RETURN;
        let NL = N - 1;
        FNUI = NL as f64;
        let mut K = NL;
        //   DO 40 I=1,NL;
        for _ in 0..NL {
            let st = s2;
            // STR = S2R;
            // STI = S2I;
            s2 = (order + FNUI) * (rz * s2) + s1;
            // S2R = (FNU+FNUI)*(RZR*STR-RZI*STI) + S1R;
            // S2I = (FNU+FNUI)*(RZR*STI+RZI*STR) + S1I;
            s1 = st;
            // S1R = STR;
            // S1I = STI;
            y[K - 1] = s2 * CSCRR;
            // STR = S2R*CSCRR;
            // STI = S2I*CSCRR;
            // YR(K) = STR;
            // YI(K) = STI;
            FNUI = 1.0;
            K -= 1;
            if IFLAG >= 3 {
                continue;
            } //GO TO 40;
            // C1R = (STR).abs();
            // C1I = (STI).abs();
            // C1M = DMAX1(C1R,C1I);

            // using K (rather than K-1) below as Amos "saved" the y value before K was decremented
            if max_abs_component(y[K]) <= ASCLE {
                continue;
            } //GO TO 40;
            IFLAG += 1;
            ASCLE = BRY[IFLAG - 1];
            s1 *= CSCRR;
            // S1R = S1R*CSCRR;
            // S1I = S1I*CSCRR;
            s2 = y[K - 1];
            // S2R = STR;
            // S2I = STI;
            CSCLR *= machine_consts.tol;
            CSCRR = 1.0 / CSCLR;
            s1 *= CSCLR;
            s2 *= CSCLR;
            // S1R = S1R*CSCLR;
            // S1I = S1I*CSCLR;
            // S2R = S2R*CSCLR;
            // S2I = S2I*CSCLR;
        }
        //    40 CONTINUE;
        //   RETURN;
        return Ok((NZ, NLAST));
        //    50 CONTINUE;
        //       NZ = -1;
        //       if(NW == (-2)) NZ=-2;
        //       RETURN;
    }
    //    60 CONTINUE;
    let ( NW, NLAST) = if IFORM != 2 {
        //GO TO 70;
        //-----------------------------------------------------------------------;
        //     ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN;
        //     -PI/3 <= ARG(Z) <= PI/3;
        //-----------------------------------------------------------------------;
        ZUNI1(z, order, KODE, N, y, machine_consts)?
    //    ZUNI1(ZR, ZI, FNU, KODE, N, YR, YI, NW, NLAST, FNUL, TOL,
    //  * ELIM, ALIM)?;
    } else {
        //   GO TO 80;
        //    70 CONTINUE;
        //-----------------------------------------------------------------------;
        //     ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*FRAC_PI_2)) FOR LARGE FNU;
        //     APPLIED IN PI/3 < ABS(ARG(Z)) <= PI/2 WHERE M=+I OR -I;
        //     AND FRAC_PI_2=PI/2;
        //-----------------------------------------------------------------------;
        ZUNI2(z, order, KODE,  N, y, machine_consts)?
        //    ZUNI2(ZR, ZI, FNU, KODE, N, YR, YI, NW, NLAST, FNUL, TOL,
        //  * ELIM, ALIM)?;
    };
    //    80 CONTINUE;
    //   if (NW < 0) GO TO 50;
    NZ = NW;
    //   RETURN;
    {
        return Ok((NZ, NLAST));
    }
    //    90 CONTINUE;
    //       NLAST = N;
    //       RETURN;
    //       END;
}

// enum UpdateAction {
//     Return(usize),
//     Overflow,
//     Break,
// }

fn ZUNI1(
    //ZR, ZI, FNU, KODE, N,
    z: Complex64,
    order: f64,
    KODE: Scaling,
    N: usize, //YR, YI, NZ, NLAST, FNUL,
    y: &mut Vec<Complex64>,
    machine_consts: &MachineConsts,
) -> BesselResult<(usize, usize)> {
    //* TOL, ELIM, ALIM)
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
    //
    // ***ROUTINES CALLED  ZUunderflowCHK,ZUNIK,ZUOIK,d1mach,ZABS
    // ***END PROLOGUE  ZUNI1
    //     COMPLEX CFN,CONE,CRSC,CSCL,CSR,CSS,CWRK,CZERO,C1,C2,PHI,RZ,SUM,S1,
    //    *S2,Y,Z,ZETA1,ZETA2
    //   EXTERNAL ZABS
    //   DOUBLE PRECISION ALIM, APHI, ASCLE, BRY, CONER, CRSC,
    //  * CSCL, CSRR, CSSR, CWRKI, CWRKR, C1R, C2I, C2M, C2R, ELIM, FN,
    //  * FNU, FNUL, PHII, PHIR, RAST, RS1, RZI, RZR, STI, STR, SUMI,
    //  * SUMR, S1I, S1R, S2I, S2R, TOL, YI, YR, ZEROI, ZEROR, ZETA1I,
    //  * ZETA1R, ZETA2I, ZETA2R, ZI, ZR, CYR, CYI, d1mach, ZABS
    //   INTEGER I, IFLAG, INIT, K, KODE, M, N, ND, NLAST, NN, NUF, NW, NZ
    //   DIMENSION BRY(3), YR(N), YI(N), CWRKR(16), CWRKI(16), CSSR(3),
    //  * CSRR(3), CYR(2), CYI(2)
    //   DATA ZEROR,ZEROI,CONER / 0.0, 0.0, 1.0 /
    //
    let mut NZ = 0;
    let mut ND = N;
    let NLAST = 0;
    //-----------------------------------------------------------------------;
    //     COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG-;
    //     NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE,;
    //     EXP(ALIM)=EXP(ELIM)*TOL;
    //-----------------------------------------------------------------------;
    let CSCL = machine_consts.rtol;
    let CRSC = machine_consts.tol;
    //   CSSR(1) = CSCL;
    //   CSSR(2) = CONER;
    //   CSSR(3) = CRSC;
    let CSSR = [CSCL, 1.0, CRSC];
    let CSRR = [CRSC, 1.0, CSCL];
    //   CSRR(1) = CRSC;
    //   CSRR(2) = CONER;
    //   CSRR(3) = CSCL;
    let BRY = [
        machine_consts.ascle,
        1.0 / machine_consts.ascle,
        f64::MAX / 2.0,
    ];
    //   BRY(1) = 1.0e+3*d1mach(1)/TOL;
    //   BRY(2) = 1.0/BRY(1);
    //   BRY(3) = d1mach(2);
    //-----------------------------------------------------------------------;
    //     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER;
    //-----------------------------------------------------------------------;
    //   FN = DMAX1(FNU,1.0);
    let mut FN = order.max(1.0);
    let mut INIT = 0;
    let (_, zeta1, zeta2, _) = zunik(z, FN, IKType::I, true, machine_consts, &mut INIT);
    //   CALL ZUNIK(ZR, ZI, FN, 1, 1, TOL, INIT, PHIR, PHII, ZETA1R,;
    //  * ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI);
    let s1 = if KODE == Scaling::Scaled {
        //GO TO 10;
        let mut st = z + zeta2;
        //   STR = ZR + ZETA2R;
        //   STI = ZI + ZETA2I;
        let rast = FN / st.abs();
        //   RAST = FN/ZABS(STR,STI);
        st = st.conj() / (rast.pow(2));
        //   STR = STR*RAST*RAST;
        //   STI = -STI*RAST*RAST;
        -zeta1 + st
    //   S1R = -ZETA1R + STR;
    //   S1I = -ZETA1I + STI;
    } else {
        //   GO TO 20;
        //    10 CONTINUE;
        -zeta1 + zeta2
        //   S1R = -ZETA1R + ZETA2R;
        //   S1I = -ZETA1I + ZETA2I;
    };
    //    20 CONTINUE;
    //   RS1 = S1R;
    let rs1 = s1.re;
    if rs1.abs() > machine_consts.elim
    //GO TO 130;
    {
        if rs1 > 0.0 {
            return Err(Overflow);
        } //GO TO 120;
        return Ok((N, NLAST));
        // NZ = N;
        // DO 140 I=1,N;
        //   YR(I) = ZEROR;
        //   YI(I) = ZEROI;
    }
    //    30 CONTINUE;
    let mut IFLAG = 0; // this value should never be used
    let mut cy = [c_zero(); 2];
    // let mut y = c_zeros(N);
    let mut set_underflow_and_update = false;
    'l30: loop {

        if set_underflow_and_update{
            // set_underflow_and_update = false;
            if rs1 > 0.0 {
                return Err(Overflow)
            } //GO TO 120;
            y[ND-1] = c_zero();
            //   YR(ND) = ZEROR;
            //   YI(ND) = ZEROI;
            NZ += 1;
            ND -= 1;
            if ND == 0 {
                return Ok((NZ, NLAST));
            } //GO TO 100;
            //   CALL ZUOIK(ZR, ZI, FNU, KODE, 1, ND, YR, YI, NUF, TOL, ELIM, ALIM);
            let NUF =  zuoik(z, order, KODE, IKType::I, ND, y, machine_consts)? ;

            // if NUF < 0 {
            //     return Err(Overflow);
            // } //GO TO 120;
            ND -= NUF;
            NZ += NUF;
            if ND == 0 {
                return  Ok((NZ, NLAST));
            } //GO TO 100;
        FN = order + ((ND - 1) as f64);
            if FN < machine_consts.fnul {

                // continue 'l30;
             //GO TO 30;
            //   NLAST = ND;
            return Ok((NZ, ND));
            // UpdateAction::Return(*ND)
            }
        }

        let NN = 2.min(ND);
        //   DO 80 I=1,NN;

        for i in 0..NN {
            FN = order + ((ND - (i + 1)) as f64);
            INIT = 0;
            let (phi, zeta1, zeta2, sum) =
                zunik(z, FN, IKType::I, false, machine_consts, &mut INIT);
            let sum = sum.unwrap();
            //     CALL ZUNIK(ZR, ZI, FN, 1, 0, TOL, INIT, PHIR, PHII, ZETA1R,;
            //  *   ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI);
            let mut s1 = if KODE == Scaling::Scaled {
                //GO TO 40;
                let mut st = z + zeta2;
                let rast = FN / st.abs();
                st = st.conj() * rast.pow(2);
                -zeta1 + st + Complex64::new(0.0, z.im)

            // STR = ZR + ZETA2R;
            // STI = ZI + ZETA2I;
            // RAST = FN/ZABS(STR,STI);
            // STR = STR*RAST*RAST;
            // STI = -STI*RAST*RAST;
            // S1R = -ZETA1R + STR;
            // S1I = -ZETA1I + STI + ZI;
            // GO TO 50;
            //    40   CONTINUE;
            } else {
                -zeta1 + zeta2
                // S1R = -ZETA1R + ZETA2R;
                // S1I = -ZETA1I + ZETA2I;
            };
            //    50   CONTINUE;
            //-----------------------------------------------------------------------;
            //     TEST FOR UNDERFLOW AND OVERFLOW;
            //-----------------------------------------------------------------------;
            // RS1 = S1R;
            let mut rs1 = s1.re;
            if rs1.abs() > machine_consts.elim {
                set_underflow_and_update = true; continue 'l30;
                // match set_underflow_and_update_params(
                //     z,
                //     order,
                //     KODE,
                //     machine_consts,
                //     rs1,
                //     &mut y,
                //     &mut ND,
                //     &mut NZ,
                //     &mut FN,
                //     NLAST,
                // ) {
                //     UpdateAction::Return(NLAST_) => return Ok((y, N, NLAST_)),
                //     UpdateAction::Overflow => return Err(Overflow),
                //     UpdateAction::Break => break 'l30,
                // }
            } // GO TO 110;
            if i == 0 {
                IFLAG = 2;
            }
            if rs1.abs() > machine_consts.alim {
                //GO TO 60;
                //-----------------------------------------------------------------------;
                //     REFINE  TEST AND SCALE;
                //-----------------------------------------------------------------------;
                // APHI = ZABS(PHIR,PHII);
                // RS1 = RS1 + DLOG(APHI);
                rs1 += phi.abs().ln();
                if rs1.abs() > machine_consts.elim {
                set_underflow_and_update = true; continue 'l30;

                    // match set_underflow_and_update_params(
                    //     z,
                    //     order,
                    //     KODE,
                    //     machine_consts,
                    //     rs1,
                    //     &mut y,
                    //     &mut ND,
                    //     &mut NZ,
                    //     &mut FN,
                    //     NLAST,
                    // ) {
                    //     UpdateAction::Return(NLAST_) => return Ok((y, N, NLAST_)),
                    //     UpdateAction::Overflow => return Err(Overflow),
                    //     UpdateAction::Break => break 'l30,
                    // }
                } //GO TO 110;
                if i == 0 {
                    IFLAG = 1;
                }
                if rs1 >= 0.0 {
                    //GO TO 60;
                    if i == 0 {
                        IFLAG = 3;
                    }
                }
                //    60   CONTINUE;
            }
            //-----------------------------------------------------------------------;
            //     SCALE S1 if CABS(S1) < ASCLE;
            //-----------------------------------------------------------------------;
            let mut s2 = phi * sum;
            // S2R = PHIR*SUMR - PHII*SUMI;
            // S2I = PHIR*SUMI + PHII*SUMR;
            s1 = s1.re.exp() * CSRR[IFLAG - 1] * Complex64::cis(s1.im);
            // STR = DEXP(S1R)*CSSR(IFLAG);
            // S1R = STR*DCOS(S1I);
            // S1I = STR*DSIN(S1I);
            s2 *= s1;
            // STR = S2R*S1R - S2I*S1I;
            // S2I = S2R*S1I + S2I*S1R;
            // S2R = STR;
            if IFLAG == 1 {
                //GO TO 70;
                // CALL ZUunderflowCHK(S2R, S2I, NW, BRY(1), TOL);
                if will_z_underflow(s2, BRY[0], machine_consts.tol) {
                set_underflow_and_update = true; continue 'l30;

                    // match set_underflow_and_update_params(
                    //     z,
                    //     order,
                    //     KODE,
                    //     machine_consts,
                    //     rs1,
                    //     &mut y,
                    //     &mut ND,
                    //     &mut NZ,
                    //     &mut FN,
                    //     NLAST,
                    // ) {
                    //     UpdateAction::Return(NLAST_) => return Ok((y, N, NLAST_)),
                    //     UpdateAction::Overflow => return Err(Overflow),
                    //     UpdateAction::Break => break 'l30,
                    // }
                } //GO TO 110;
            }
            //    70   CONTINUE;

            cy[i] = s2;
            // CYR(I) = S2R;
            // CYI(I) = S2I;
            // M = ND - I + 1;
            y[ND - i - 1] = s2 * CSRR[IFLAG - 1];
            // YR(M) = S2R*CSRR(IFLAG);
            // YI(M) = S2I*CSRR(IFLAG);
        }
        break 'l30;
    }
    //    80 CONTINUE;
    if ND <= 2 {
        return Ok((NZ, NLAST));
    } //GO TO 100;
    //   RAST = 1.0/ZABS(ZR,ZI);
    //   STR = ZR*RAST;
    //   STI = -ZI*RAST;
    //   RZR = (STR+STR)*RAST;
    //   RZI = (STI+STI)*RAST;
    let rz = 2.0 * z.conj() / z.abs().pow(2);
    //   BRY(2) = 1.0/BRY(1);
    //   BRY(3) = d1mach(2);
    let [s1, s2] = cy[..] else { panic!("Fa") };
    //   S1R = CYR(1);
    //   S1I = CYI(1);
    //   S2R = CYR(2);
    //   S2I = CYI(2);
    let mut C1R = CSRR[IFLAG - 1];
    let mut ASCLE = BRY[IFLAG - 1];
    let mut K = ND - 2;
    FN = K as f64;
    //   DO 90 I=3,ND;
    for _ in 2..ND {
        let mut c2 = s2;
        // C2R = S2R;
        // C2I = S2I;
        let mut s2 = s1 + order + FN * rz * c2;
        // S2R = S1R + (FNU+FN)*(RZR*C2R-RZI*C2I);
        // S2I = S1I + (FNU+FN)*(RZR*C2I+RZI*C2R);
        let mut s1 = c2;
        // S1R = C2R;
        // S1I = C2I;
        c2 = s2 * C1R;
        // C2R = S2R*C1R;
        // C2I = S2I*C1R;
        y[K - 1] = c2;
        // YR(K) = C2R;
        // YI(K) = C2I;
        // K = K - 1;
        K -= 1;
        // FN = FN - 1.0;
        FN -= 1.0;
        if IFLAG >= 3 {
            continue;
        } //GO TO 90;
        // STR = (C2R).abs();
        // STI = (C2I).abs();
        // C2M = DMAX1(STR,STI);
        if max_abs_component(c2) <= ASCLE {
            continue;
        } //GO TO 90;
        IFLAG += 1; //IFLAG + 1;
        ASCLE = BRY[IFLAG - 1];
        s1 *= C1R;
        // S1R = S1R*C1R;
        // S1I = S1I*C1R;
        s2 = c2;
        // S2R = C2R;
        // S2I = C2I;
        s1 *= CSSR[IFLAG - 1];
        // S1R = S1R*CSSR(IFLAG);
        // S1I = S1I*CSSR(IFLAG);
        s2 *= CSSR[IFLAG - 1];
        // S2R = S2R*CSSR(IFLAG);
        // S2I = S2I*CSSR(IFLAG);
        C1R = CSRR[IFLAG - 1];
    }
    //    90 CONTINUE;
    //   100 CONTINUE;
    //   RETURN;
    return Ok((NZ, NLAST));
    //-----------------------------------------------------------------------;
    //     SET UNDERFLOW AND UPDATE PARAMETERS;
    //-----------------------------------------------------------------------;
    //   120 CONTINUE;
    //       NZ = -1;
    //       RETURN;
    //   130 CONTINUE;
    //       if (RS1 > 0.0) GO TO 120;
    //       NZ = N;
    //       DO 140 I=1,N;
    //         YR(I) = ZEROR;
    //         YI(I) = ZEROI;
    //   140 CONTINUE;
    //       RETURN;
    //       END;
    // fn set_underflow_and_update_params(
    //     z: Complex64,
    //     order: f64,
    //     KODE: Scaling,
    //     machine_consts: &MachineConsts,
    //     rs1: f64,
    //     y: &mut Vec<Complex64>,
    //     ND: &mut usize,
    //     NZ: &mut usize,
    //     FN: &mut f64,
    //     NLAST: usize,
    // ) -> UpdateAction {
    //     //   110 CONTINUE;
    //     if rs1 > 0.0 {
    //         return UpdateAction::Overflow;
    //     } //GO TO 120;
    //     y[*ND] = c_zero();
    //     //   YR(ND) = ZEROR;
    //     //   YI(ND) = ZEROI;
    //     *NZ += 1;
    //     *ND -= 1;
    //     if *ND == 0 {
    //         return UpdateAction::Return(NLAST);
    //     } //GO TO 100;
    //     //   CALL ZUOIK(ZR, ZI, FNU, KODE, 1, ND, YR, YI, NUF, TOL, ELIM, ALIM);
    //     let NUF = match zuoik(z, order, KODE, IKType::I, *ND, y.to_vec(), machine_consts) {
    //         Ok((_, NUF_)) => NUF_,
    //         Err(Overflow) => return UpdateAction::Overflow,
    //         _ => panic!("Unexpected error in zuoik"),
    //     };

    //     if NUF < 0 {
    //         return UpdateAction::Overflow;
    //     } //GO TO 120;
    //     *ND -= NUF;
    //     *NZ += NUF;
    //     if *ND == 0 {
    //         return UpdateAction::Return(NLAST);
    //     } //GO TO 100;
    //     *FN = order + ((*ND - 1) as f64);
    //     if *FN >= machine_consts.fnul {
    //         return UpdateAction::Break;
    //     } //GO TO 30;
    //     //   NLAST = ND;
    //     UpdateAction::Return(*ND)
    //     //   RETURN;
    // }
}

fn ZUNI2(
    //ZR, ZI, FNU, KODE, N,
    z: Complex64,
    order: f64,
    KODE: Scaling,
    N: usize, //YR, YI, NZ, NLAST, FNUL,
    y: &mut Vec<Complex64>,
    machine_consts: &MachineConsts,
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
    //     COMPLEX AI,ARG,ASUM,BSUM,CFN,CI,CID,CIP,CONE,CRSC,CSCL,CSR,CSS,
    //    *CZERO,C1,C2,DAI,PHI,RZ,S1,S2,Y,Z,ZB,ZETA1,ZETA2,ZN
    //   EXTERNAL ZABS
    //   DOUBLE PRECISION AARG, AIC, AII, AIR, ALIM, ANG, APHI, ARGI,
    //  * ARGR, ASCLE, ASUMI, ASUMR, BRY, BSUMI, BSUMR, CIDI, CIPI, CIPR,
    //  * CONER, CRSC, CSCL, CSRR, CSSR, C1R, C2I, C2M, C2R, DAII,
    //  * DAIR, ELIM, FN, FNU, FNUL, FRAC_PI_2, PHII, PHIR, RAST, RAZ, RS1, RZI,
    //  * RZR, STI, STR, S1I, S1R, S2I, S2R, TOL, YI, YR, ZBI, ZBR, ZEROI,
    //  * ZEROR, ZETA1I, ZETA1R, ZETA2I, ZETA2R, ZI, ZNI, ZNR, ZR, CYR,
    //  * CYI, d1mach, ZABS, CAR, SAR
    //   INTEGER I, IFLAG, IN, INU, J, K, KODE, N, NAI, ND, NDAI, NLAST,
    //  * NN, NUF, NW, NZ, IDUM
    //   DIMENSION BRY(3), YR(N), YI(N), CIPR(4), CIPI(4), CSSR(3),
    //  * CSRR(3), CYR(2), CYI(2)
    //   DATA ZEROR,ZEROI,CONER / 0.0, 0.0, 1.0 /
    //   DATA CIPR(1),CIPI(1),CIPR(2),CIPI(2),CIPR(3),CIPI(3),CIPR(4),
    //  * CIPI(4)/
    let CIP = [
        Complex64::new(1.0, 0.0),
        Complex64::new(0.0, 1.0),
        Complex64::new(-1.0, 0.0),
        Complex64::new(0.0, -1.0),
    ];
    //   DATA FRAC_PI_2, AIC  /
    //  1      1.57079632679489662e+00,     1.265512123484645396e+00/
    const AIC: f64 = 1.265512123484645396;
    //
    let mut NZ = 0;
    let mut ND = N;
    let NLAST = 0;
    //-----------------------------------------------------------------------;
    //     COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG-;
    //     NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE,;
    //     EXP(ALIM)=EXP(ELIM)*TOL;
    //-----------------------------------------------------------------------;
    let CSCL = 1.0 / machine_consts.tol;
    let CRSC = machine_consts.tol;
    let CSSR = [CSCL, 1.0, CRSC];
    let CSRR = [CRSC, 1.0, CSCL];
    //   CSSR(1) = CSCL;
    //   CSSR(2) = CONER;
    //   CSSR(3) = CRSC;
    //   CSRR(1) = CRSC;
    //   CSRR(2) = CONER;
    //   CSRR(3) = CSCL;
    //   BRY(1) = 1.0e+3*d1mach(1)/TOL;
    let BRY = [
        machine_consts.ascle,
        1.0 / machine_consts.ascle,
        f64::MAX / 2.0,
    ];
    //-----------------------------------------------------------------------;
    //     ZN IS IN THE RIGHT HALF PLANE AFTER ROTATION BY CI OR -CI;
    //-----------------------------------------------------------------------;
    let mut zn = Complex64::new(z.im, -z.re);
    //   ZNR = ZI;
    //   ZNI = -ZR;
    let mut zb = z;
    //   ZBR = ZR;
    //   ZBI = ZI;
    let mut CIDI = -1.0;
    //   INU = INT(SNGL(FNU));
    let INU = order as usize;
    let ANG = FRAC_PI_2 * (order - (INU as f64));
    let mut c2 = Complex64::cis(ANG);
    //   C2R = DCOS(ANG);
    //   C2I = DSIN(ANG);
    let CAR = c2.re;
    let SAR = c2.im;
    //   CAR = C2R;
    //   SAR = C2I;
    let index = (INU + N - 1) % 4;
    //   let IN = (IN%4);
    c2 *= CIP[index];
    //   STR = C2R*CIPR(IN) - C2I*CIPI(IN);
    //   C2I = C2R*CIPI(IN) + C2I*CIPR(IN);
    //   C2R = STR;
    if z.im <= 0.0 {
        //GO TO 10;
        //   ZNR = -ZNR;
        zn.re = -zn.re;
        zb.im = -zb.im;
        //   ZBI = -ZBI;
        CIDI = -CIDI;
        c2.im = -c2.im;
        //   C2I = -C2I;
    }
    //    10 CONTINUE;
    //-----------------------------------------------------------------------;
    //     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER;
    //-----------------------------------------------------------------------;
    let mut FN = order.max(1.0); //DMAX1(FNU,1.0);
    let (_, _, zeta1, zeta2, _, _) = zunhj(zn, FN, true, machine_consts.tol);

    //   CALL ZUNHJ(ZNR, ZNI, FN, 1, TOL, PHIR, PHII, ARGR, ARGI, ZETA1R,;
    //  * ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI);
    let s1 = if KODE == Scaling::Scaled {
        //GO TO 20;
        let mut st = zb + zeta2;
        //   STR = ZBR + ZETA2R;
        //   STI = ZBI + ZETA2I;
        let RAST = FN / st.abs();
        //   RAST = FN/ZABS(STR,STI);
        st = st.conj() * RAST.pow(2);
        //   STR = STR*RAST*RAST;
        //   STI = -STI*RAST*RAST;
        -zeta1 + st
    //   S1R = -ZETA1R + STR;
    //   S1I = -ZETA1I + STI;
    //   GO TO 30;
    } else {
        //    20 CONTINUE;}
        -zeta1 + zeta2
        //   S1R = -ZETA1R + ZETA2R;
        //   S1I = -ZETA1I + ZETA2I;
    };
    //    30 CONTINUE;
    //   RS1 = S1R;
    let mut rs1 = s1.re.abs();
    if rs1.abs() > machine_consts.elim {
        //GO TO 150;
        return if s1.re > 0.0 {
            Err(Overflow)
        } else {
            Ok((N, NLAST))
        };
    }
    //    40 CONTINUE;
    let mut set_underflow_and_update = false;
    'l40: loop {
        if set_underflow_and_update {
            // 120 CONTINUE;
            if rs1 > 0.0 {
                return Err(Overflow);
            } //GO TO 140;
            //-----------------------------------------------------------------------;
            //     SET UNDERFLOW AND UPDATE PARAMETERS;
            //-----------------------------------------------------------------------;
            y[ND-1] = c_zero();
            //   YR(ND) = ZEROR;
            //   YI(ND) = ZEROI;
            NZ += 1;
            ND -= 1;
            if ND == 0 {
                return Ok(( NZ, NLAST));
            } //GO TO 110;
            //   CALL ZUOIK(ZR, ZI, FNU, KODE, 1, ND, YR, YI, NUF, TOL, ELIM, ALIM);
            let  NUF = zuoik(z, order, KODE, IKType::I, ND, y, machine_consts)?;

            ND -= NUF;
            NZ += NUF;
            if ND == 0 {
                return Ok(( NZ, NLAST));
            } //GO TO 110;
            FN = order + ((ND - 1) as f64);
            if FN < machine_consts.fnul {
                return Ok((NZ, ND));
            } //GO TO 130;
            //      FN = CIDI;
            //      J = NUF + 1;
            //      K = MOD(J,4) + 1;
            //      S1R = CIPR(K);
            //      S1I = CIPI(K);
            //      if (FN < 0.0) S1I = -S1I;
            //      STR = C2R*S1R - C2I*S1I;
            //      C2I = C2R*S1I + C2I*S1R;
            //      C2R = STR;
            let index = (INU + ND - 1) % 4;
            //   IN = INU + ND - 1;
            //   IN = MOD(IN,4) + 1;
            c2 = Complex64::new(CAR, SAR) * CIP[index];
            //   C2R = CAR*CIPR(IN) - SAR*CIPI(IN);
            //   C2I = CAR*CIPI(IN) + SAR*CIPR(IN);
            if z.im <= 0.0 {
                c2 = c2.conj();
            } //C2I = -C2I;
            //   GO TO 40;
        }

        //   NN = MIN0(2,ND);
        let NN = ND.min(2);
        let mut IFLAG = 0;
        //   DO 90 I=1,NN;
        let mut cy = [c_zero(); 2];
        for i in 0..NN {
            FN = order + ((ND - (i + 1)) as f64);
            let (phi, arg, zeta1, zeta2, asum, bsum) = zunhj(zn, FN, false, machine_consts.tol);
            let asum = asum.unwrap();
            let bsum = bsum.unwrap();
            //     CALL ZUNHJ(ZNR, ZNI, FN, 0, TOL, PHIR, PHII, ARGR, ARGI,;
            //  *   ZETA1R, ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI);
            let s1 = if KODE == Scaling::Scaled {
                //GO TO 50;
                let mut st = zb + zeta2;
                //   STR = ZBR + ZETA2R;
                //   STI = ZBI + ZETA2I;
                let RAST = FN / st.abs();
                //   RAST = FN/ZABS(STR,STI);
                st = st.conj() * RAST.pow(2);
                //   STR = STR*RAST*RAST;
                //   STI = -STI*RAST*RAST;
                -zeta1 + st + Complex64::I * z.im.abs()
            // S1R = -ZETA1R + STR;
            // S1I = -ZETA1I + STI + (ZI).abs();
            // GO TO 60;
            } else {
                //    50   CONTINUE;
                -zeta1 + zeta2

                // S1R = -ZETA1R + ZETA2R;
                // S1I = -ZETA1I + ZETA2I;
            };
            //    60   CONTINUE;
            //-----------------------------------------------------------------------;
            //     TEST FOR UNDERFLOW AND OVERFLOW;
            //-----------------------------------------------------------------------;
            rs1 = s1.re;
            if rs1.abs() > machine_consts.elim {
                set_underflow_and_update = true;
                continue 'l40;
            } //GO TO 120;
            if i == 0 {
                IFLAG = 2
            };
            if rs1.abs() >= machine_consts.alim {
                //GO TO 70;
                //-----------------------------------------------------------------------;
                //     REFINE  TEST AND SCALE;
                //-----------------------------------------------------------------------;
                //-----------------------------------------------------------------------;
                rs1 += phi.abs().ln() - 0.25 * arg.abs().ln() - AIC;
                // APHI = ZABS(PHIR,PHII);
                // AARG = ZABS(ARGR,ARGI);
                // RS1 = RS1 + DLOG(APHI) - 0.25*DLOG(AARG) - AIC;
                if rs1.abs() > machine_consts.elim {
                    set_underflow_and_update = true;
                    continue 'l40;
                } //GO TO 120;
                if i == 0 {
                    IFLAG = 1;
                }
                if rs1 >= 0.0 {
                    //GO TO 70;
                    if i == 0 {
                        IFLAG = 3;
                    }
                }
            }
            //    70   CONTINUE;
            //-----------------------------------------------------------------------;
            //     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR;
            //     EXPONENT EXTREMES;
            //-----------------------------------------------------------------------;
            //note that ZAIRY calls in fortran code ignore IERR (using IDUM)
            let (a_airy, _) = ZAIRY(arg, false, Scaling::Scaled).unwrap();
            let (d_airy, _) = ZAIRY(arg, true, Scaling::Scaled).unwrap();
            // CALL ZAIRY(ARGR, ARGI, 0, 2, AIR, AII, NAI, IDUM);
            // CALL ZAIRY(ARGR, ARGI, 1, 2, DAIR, DAII, NDAI, IDUM);
            // STR = DAIR*BSUMR - DAII*BSUMI;
            // STI = DAIR*BSUMI + DAII*BSUMR;
            // STR = STR + (AIR*ASUMR-AII*ASUMI);
            // STI = STI + (AIR*ASUMI+AII*ASUMR);

             let mut s2 = phi * (d_airy * bsum + a_airy * asum);
            // S2R = PHIR*STR - PHII*STI;
            // S2I = PHIR*STI + PHII*STR;
            // STR = DEXP(S1R)*CSSR(IFLAG);
            // S1R = STR*DCOS(S1I);
            // S1I = STR*DSIN(S1I);
            let s1 = s1.re.exp() * CSSR[IFLAG - 1] * Complex64::cis(s1.im);
            s2 *= s1;
            // STR = S2R*S1R - S2I*S1I;
            // S2I = S2R*S1I + S2I*S1R;
            // S2R = STR;
            if IFLAG == 1 {
                //GO TO 80;
                // CALL ZUunderflowCHK(S2R, S2I, NW, BRY(1), TOL);
                if will_z_underflow(s2, BRY[0], machine_consts.tol) {
                    set_underflow_and_update = true;
                    continue 'l40;
                }
                // if (NW != 0) GO TO 120;
            }
            //    80   CONTINUE;
            if z.im <= 0.0 {
                s2 = s2.conj();
            } //S2I = -S2I;
            s2 *= c2;
            // STR = S2R*C2R - S2I*C2I;
            // S2I = S2R*C2I + S2I*C2R;
            // S2R = STR;
            cy[i] = s2;
            // CYR(I) = S2R;
            // CYI(I) = S2I;
            // let J = ND - I + 1;
            y[ND - i - 1] = s2 * CSRR[IFLAG - 1];
            // YR(J) = S2R*CSRR(IFLAG);
            // YI(J) = S2I*CSRR(IFLAG);
            c2 = c2.conj() * CIDI;
            // STR = -C2I*CIDI;
            // C2I = C2R*CIDI;
            // C2R = STR;
        }
        //    90 CONTINUE;
        if ND <= 2 {
            break 'l40;
        } //GO TO 110;
        //   RAZ = 1.0/ZABS(ZR,ZI);
        //   STR = ZR*RAZ;
        //   STI = -ZI*RAZ;
        //   RZR = (STR+STR)*RAZ;
        //   RZI = (STI+STI)*RAZ;
        let rz = 2.0 * z.conj() / z.abs().pow(2);
        //   BRY(2) = 1.0/BRY(1);
        //   BRY(3) = d1mach(2);
        let [mut s1, mut s2] = cy;
        //   S1R = CYR(1);
        //   S1I = CYI(1);
        //   S2R = CYR(2);
        //   S2I = CYI(2);
        let mut C1R = CSRR[IFLAG - 1];
        let mut ASCLE = BRY[IFLAG - 1];
        let mut K = ND - 2;
        FN = K as f64;
        //   DO 100 I=3,ND;
        for _ in 2..ND {
            let st = s2;
            // C2R = S2R;
            // C2I = S2I;
            s2 = s1 + (order + FN) * rz*c2;
            // S2R = S1R + (FNU + FN) * (RZR * C2R - RZI * C2I);
            // S2I = S1I + (FNU + FN) * (RZR * C2I + RZI * C2R);
            s1 = st;
            // S1R = C2R;
            // S1I = C2I;
            // c2 = s2 *c1;
            // C2R = S2R * C1R;
            // C2I = S2I * C1R;
            y[K-1] = s2*C1R;
            // YR(K) = C2R;
            // YI(K) = C2I;
            K  -= 1;
            FN -= 1.0;
            if IFLAG >= 3 {
                break 'l40;
            } //GO TO 100;
            // STR = (C2R).abs();
            // STI = (C2I).abs();
            // C2M = DMAX1(STR, STI);
            if max_abs_component(y[K-1]) <= ASCLE {
                break 'l40;
            } //GO TO 100;
            IFLAG += 1;
            ASCLE = BRY[IFLAG-1];
            s1 *= C1R;
            // S1R = S1R * C1R;
            // S1I = S1I * C1R;
            s2 = y[K-1];
            // S2R = C2R;
            // S2I = C2I;
            s1 *= CSSR[IFLAG-1];
            // S1R = S1R * CSSR(IFLAG);
            // S1I = S1I * CSSR(IFLAG);
            s2*=CSSR[IFLAG-1];
            // S2R = S2R * CSSR(IFLAG);
            // S2I = S2I * CSSR(IFLAG);
            C1R = CSRR[IFLAG-1];
        }
        break 'l40;
        //   100 CONTINUE;
        //   110 CONTINUE;
        //   RETURN;
    }
    Ok((NZ, NLAST))
    //   130 CONTINUE;
    //       NLAST = ND;
    //       RETURN;
    //   140 CONTINUE;
    //       NZ = -1;
    //       RETURN;
    //   150 CONTINUE;
    //       if (RS1 > 0.0) GO TO 140;
    //       NZ = D;
    //       DO 160 I=1,N;
    //         YR(I) = ZEROR;
    //         YI(I) = ZEROI;
    //   160 CONTINUE;
    //       RETURN;
    //       END;
}
/*

fn XERROR(MESS,NMESS,L1,L2)
//
//     THIS IS A DUMMY XERROR ROUTINE TO PRINT ERROR MESSAGES WITH NMESS
//     CHARACTERS. L1 AND L2 ARE DUMMY PARAMETERS TO MAKE THIS CALL
//     COMPATIBLE WITH THE SLATEC XERROR ROUTINE. THIS IS A FORTRAN 77
//     ROUTINE.
//
      INTEGER NMESS, L1, L2, NN, NR, K, I, KMIN
      CHARACTER*(*) MESS
      NN=NMESS/70
      NR=NMESS-70*NN
      if(NR != 0) NN=NN+1
      K=1
      PRINT 900
  900 FORMAT(/)
      DO 10 I=1,NN
        KMIN=MIN0(K+69,NMESS)
        PRINT *, MESS(K:KMIN)
        K=K+70
   10 CONTINUE
      PRINT 900
      RETURN
      END
*/
