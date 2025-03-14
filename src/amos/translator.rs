#![allow(non_snake_case)]
use super::{
    BesselError, IKType, MachineConsts, Scaling, c_one, c_zero, gamma_ln,
    machine::{d1mach, i1mach},
    overflow_checks::zuoik,
    z_power_series,
};
use crate::amos::{BesselError::*, z_asymptotic_i::z_asymptotic_i};
use num::complex::{Complex64, ComplexFloat};
use std::f64::{self, consts::FRAC_PI_2};

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
) -> Result<(Vec<Complex64>, usize), BesselError> {
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
    if order == 0.0_f64 {
        err = Some("order cannot be zero")
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
    if NL == 0 {
        return Ok((cy, NZ));
    }
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
    Ok((cy, NZ))
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
fn ZAIRY(ZR, ZI, ID, KODE, AIR, AII, NZ, IERR)
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
      EXTERNAL ZABS
      DOUBLE PRECISION AA, AD, AII, AIR, AK, ALIM, ATRM, AZ, AZ3, BK,
     * CC, CK, COEF, CONEI, CONER, CSQI, CSQR, CYI, CYR, C1, C2, DIG,
     * DK, D1, D2, ELIM, FID, FNU, PTR, RL, R1M5, SFAC, STI, STR,
     * S1I, S1R, S2I, S2R, TOL, TRM1I, TRM1R, TRM2I, TRM2R, TTH, ZEROI,
     * ZEROR, ZI, ZR, ZTAI, ZTAR, Z3I, Z3R, d1mach, ZABS, ALAZ, BB
      INTEGER ID, IERR, IFLAG, K, KODE, K1, K2, MR, NN, NZ, i1mach
      DIMENSION CYR(1), CYI(1)
      DATA TTH, C1, C2, COEF /6.66666666666666667e-01,
     * 3.55028053887817240e-01,2.58819403792806799e-01,
     * 1.83776298473930683e-01/
      DATA ZEROR, ZEROI, CONER, CONEI /0.0,0.0,1.0,0.0/
// ***FIRST EXECUTABLE STATEMENT  ZAIRY
      IERR = 0
      NZ=0
      if (ID < 0 || ID > 1) IERR=1
      if (KODE < 1 || KODE > 2) IERR=1
      if (IERR != 0) RETURN
      AZ = ZABS(ZR,ZI)
      TOL = DMAX1(d1mach(4),1.0e-18)
      FID = (ID as f64)
      if (AZ > 1.0) GO TO 70
//-----------------------------------------------------------------------
//     POWER SERIES FOR CABS(Z) <= 1.
//-----------------------------------------------------------------------
      S1R = CONER
      S1I = CONEI
      S2R = CONER
      S2I = CONEI
      if (AZ < TOL) GO TO 170
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
      AIR = S1R*C1 - C2*(ZR*S2R-ZI*S2I)
      AII = S1I*C1 - C2*(ZR*S2I+ZI*S2R)
      if (KODE == 1) RETURN
      CALL ZSQRT(ZR, ZI, STR, STI)
      ZTAR = TTH*(ZR*STR-ZI*STI)
      ZTAI = TTH*(ZR*STI+ZI*STR)
      CALL ZEXP(ZTAR, ZTAI, STR, STI)
      PTR = AIR*STR - AII*STI
      AII = AIR*STI + AII*STR
      AIR = PTR
      RETURN
   50 CONTINUE
      AIR = -S2R*C2
      AII = -S2I*C2
      if (AZ <= TOL) GO TO 60
      STR = ZR*S1R - ZI*S1I
      STI = ZR*S1I + ZI*S1R
      CC = C1/(1.0+FID)
      AIR = AIR + CC*(STR*ZR-STI*ZI)
      AII = AII + CC*(STR*ZI+STI*ZR)
   60 CONTINUE
      if (KODE == 1) RETURN
      CALL ZSQRT(ZR, ZI, STR, STI)
      ZTAR = TTH*(ZR*STR-ZI*STI)
      ZTAI = TTH*(ZR*STI+ZI*STR)
      CALL ZEXP(ZTAR, ZTAI, STR, STI)
      PTR = STR*AIR - STI*AII
      AII = STR*AII + STI*AIR
      AIR = PTR
      RETURN
//-----------------------------------------------------------------------
//     CASE FOR CABS(Z) > 1.0
//-----------------------------------------------------------------------
   70 CONTINUE
      FNU = (1.0+FID)/3.0
//-----------------------------------------------------------------------
//     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
//     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0e-18.
//     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
//     EXP(-ELIM) < EXP(-ALIM)=EXP(-ELIM)/TOL    AND
//     EXP(ELIM) > EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
//     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
//     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
//     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
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
      ALAZ = DLOG(AZ)
//--------------------------------------------------------------------------
//     TEST FOR PROPER RANGE
//-----------------------------------------------------------------------
      AA=0.5/TOL
      BB=DBLE(FLOAT(i1mach(9)))*0.5
      AA=DMIN1(AA,BB)
      AA=AA**TTH
      if (AZ > AA) {return Err(LossOfSignificance);}
      AA=DSQRT(AA)
      if (AZ > AA) IERR=3
      CALL ZSQRT(ZR, ZI, CSQR, CSQI)
      ZTAR = TTH*(ZR*CSQR-ZI*CSQI)
      ZTAI = TTH*(ZR*CSQI+ZI*CSQR)
//-----------------------------------------------------------------------
//     RE(ZTA) <= 0 WHEN RE(Z) < 0, ESPECIALLY WHEN IM(Z) IS SMALL
//-----------------------------------------------------------------------
      IFLAG = 0
      SFAC = 1.0
      AK = ZTAI
      if (ZR >= 0.0) GO TO 80
      BK = ZTAR
      CK = -(BK).abs()
      ZTAR = CK
      ZTAI = AK
   80 CONTINUE
      if (ZI != 0.0) GO TO 90
      if (ZR > 0.0) GO TO 90
      ZTAR = 0.0
      ZTAI = AK
   90 CONTINUE
      AA = ZTAR
      if (AA >= 0.0 && ZR > 0.0) GO TO 110
      if (KODE == 2) GO TO 100
//-----------------------------------------------------------------------
//     OVERFLOW TEST
//-----------------------------------------------------------------------
      if (AA > (-ALIM)) GO TO 100
      AA = -AA + 0.25*ALAZ
      IFLAG = 1
      SFAC = TOL
      if (AA > ELIM) GO TO 270
  100 CONTINUE
//-----------------------------------------------------------------------
//     CBKNU AND CACON RETURN EXP(ZTA)*K(FNU,ZTA) ON KODE=2
//-----------------------------------------------------------------------
      MR = 1
      if (ZI < 0.0) MR = -1
      CALL ZACAI(ZTAR, ZTAI, FNU, KODE, MR, 1, CYR, CYI, NN, RL, TOL,
     * ELIM, ALIM)
      if (NN < 0) GO TO 280
      NZ = NZ + NN
      GO TO 130
  110 CONTINUE
      if (KODE == 2) GO TO 120
//-----------------------------------------------------------------------
//     UNDERFLOW TEST
//-----------------------------------------------------------------------
      if (AA < ALIM) GO TO 120
      AA = -AA - 0.25*ALAZ
      IFLAG = 2
      SFAC = 1.0/TOL
      if (AA < (-ELIM)) GO TO 210
  120 CONTINUE
      CALL ZBKNU(ZTAR, ZTAI, FNU, KODE, 1, CYR, CYI, NZ, TOL, ELIM,
     * ALIM)
  130 CONTINUE
      S1R = CYR(1)*COEF
      S1I = CYI(1)*COEF
      if (IFLAG != 0) GO TO 150
      if (ID == 1) GO TO 140
      AIR = CSQR*S1R - CSQI*S1I
      AII = CSQR*S1I + CSQI*S1R
      RETURN
  140 CONTINUE
      AIR = -(ZR*S1R-ZI*S1I)
      AII = -(ZR*S1I+ZI*S1R)
      RETURN
  150 CONTINUE
      S1R = S1R*SFAC
      S1I = S1I*SFAC
      if (ID == 1) GO TO 160
      STR = S1R*CSQR - S1I*CSQI
      S1I = S1R*CSQI + S1I*CSQR
      S1R = STR
      AIR = S1R/SFAC
      AII = S1I/SFAC
      RETURN
  160 CONTINUE
      STR = -(S1R*ZR-S1I*ZI)
      S1I = -(S1R*ZI+S1I*ZR)
      S1R = STR
      AIR = S1R/SFAC
      AII = S1I/SFAC
      RETURN
  170 CONTINUE
      AA = 1.0e+3*d1mach(1)
      S1R = ZEROR
      S1I = ZEROI
      if (ID == 1) GO TO 190
      if (AZ <= AA) GO TO 180
      S1R = C2*ZR
      S1I = C2*ZI
  180 CONTINUE
      AIR = C1 - S1R
      AII = -S1I
      RETURN
  190 CONTINUE
      AIR = -C2
      AII = 0.0
      AA = DSQRT(AA)
      if (AZ <= AA) GO TO 200
      S1R = 0.5*(ZR*ZR-ZI*ZI)
      S1I = ZR*ZI
  200 CONTINUE
      AIR = AIR + C1*S1R
      AII = AII + C1*S1I
      RETURN
  210 CONTINUE
      NZ = 1
      AIR = ZEROR
      AII = ZEROI
      RETURN
  270 CONTINUE
      NZ = 0
      IERR=2
      RETURN
  280 CONTINUE
      if(NN == (-1)) GO TO 270
      NZ=0
      IERR=5
      RETURN
  260 CONTINUE
      IERR=4
      NZ=0
      RETURN
      END
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
     * TRM2R, TTH, ZI, ZR, ZTAI, ZTAR, Z3I, Z3R, d1mach, ZABS
      INTEGER ID, IERR, K, KODE, K1, K2, NZ, i1mach
      DIMENSION CYR(2), CYI(2)
      DATA TTH, C1, C2, COEF, PI /6.66666666666666667e-01,
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
      ZTAR = TTH*(ZR*STR-ZI*STI)
      ZTAI = TTH*(ZR*STI+ZI*STR)
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
      ZTAR = TTH*(ZR*STR-ZI*STI)
      ZTAI = TTH*(ZR*STI+ZI*STR)
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
      AA=AA**TTH
      if (AZ > AA) {return Err(LossOfSignificance);}
      AA=DSQRT(AA)
      if (AZ > AA) IERR=3
      CALL ZSQRT(ZR, ZI, CSQR, CSQI)
      ZTAR = TTH*(ZR*CSQR-ZI*CSQI)
      ZTAI = TTH*(ZR*CSQI+ZI*CSQR)
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
fn ZMLT(AR, AI, BR, BI, CR, CI)
// ***BEGIN PROLOGUE  ZMLT
// ***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY
//
//     DOUBLE PRECISION COMPLEX MULTIPLY, C=A*B.
//
// ***ROUTINES CALLED  (NONE)
// ***END PROLOGUE  ZMLT
      DOUBLE PRECISION AR, AI, BR, BI, CR, CI, CA, CB
      CA = AR*BR - AI*BI
      CB = AR*BI + AI*BR
      CR = CA
      CI = CB
      RETURN
      END
fn ZDIV(AR, AI, BR, BI, CR, CI)
// ***BEGIN PROLOGUE  ZDIV
// ***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY
//
//     DOUBLE PRECISION COMPLEX DIVIDE C=A/B.
//
// ***ROUTINES CALLED  ZABS
// ***END PROLOGUE  ZDIV
      EXTERNAL ZABS
      DOUBLE PRECISION AR, AI, BR, BI, CR, CI, BM, CA, CB, CC, CD
      DOUBLE PRECISION ZABS
      BM = 1.0/ZABS(BR,BI)
      CC = BR*BM
      CD = BI*BM
      CA = (AR*CC+AI*CD)*BM
      CB = (AI*CC-AR*CD)*BM
      CR = CA
      CI = CB
      RETURN
      END
fn ZSQRT(AR, AI, BR, BI)
// ***BEGIN PROLOGUE  ZSQRT
// ***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY
//
//     DOUBLE PRECISION COMPLEX SQUARE ROOT, B=CSQRT(A)
//
// ***ROUTINES CALLED  ZABS
// ***END PROLOGUE  ZSQRT
      EXTERNAL ZABS
      DOUBLE PRECISION AR, AI, BR, BI, ZM, DTHETA, DPI, DRT
      DOUBLE PRECISION ZABS
      DATA DRT , DPI / 7.071067811865475244008443621e-1,
     1                 3.141592653589793238462643383e+0/
      ZM = ZABS(AR,AI)
      ZM = DSQRT(ZM)
      if (AR == 0.0e+0) GO TO 10
      if (AI == 0.0e+0) GO TO 20
      DTHETA = DATAN(AI/AR)
      if (DTHETA <= 0.0e+0) GO TO 40
      if (AR < 0.0e+0) DTHETA = DTHETA - DPI
      GO TO 50
   10 if (AI > 0.0e+0) GO TO 60
      if (AI < 0.0e+0) GO TO 70
      BR = 0.0e+0
      BI = 0.0e+0
      RETURN
   20 if (AR > 0.0e+0) GO TO 30
      BR = 0.0e+0
      BI = DSQRT((AR).abs())
      RETURN
   30 BR = DSQRT(AR)
      BI = 0.0e+0
      RETURN
   40 if (AR < 0.0e+0) DTHETA = DTHETA + DPI
   50 DTHETA = DTHETA*0.5e+0
      BR = ZM*DCOS(DTHETA)
      BI = ZM*DSIN(DTHETA)
      RETURN
   60 BR = ZM*DRT
      BI = ZM*DRT
      RETURN
   70 BR = ZM*DRT
      BI = -ZM*DRT
      RETURN
      END
fn ZEXP(AR, AI, BR, BI)
// ***BEGIN PROLOGUE  ZEXP
// ***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY
//
//     DOUBLE PRECISION COMPLEX EXPONENTIAL FUNCTION B=EXP(A)
//
// ***ROUTINES CALLED  (NONE)
// ***END PROLOGUE  ZEXP
      DOUBLE PRECISION AR, AI, BR, BI, ZM, CA, CB
      ZM = DEXP(AR)
      CA = ZM*DCOS(AI)
      CB = ZM*DSIN(AI)
      BR = CA
      BI = CB
      RETURN
      END
fn ZLOG(AR, AI, BR, BI, IERR)
// ***BEGIN PROLOGUE  ZLOG
// ***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY
//
//     DOUBLE PRECISION COMPLEX LOGARITHM B=CLOG(A)
//     IERR=0,NORMAL RETURN      IERR=1, Z=CMPLX(0.0,0.0)
// ***ROUTINES CALLED  ZABS
// ***END PROLOGUE  ZLOG
      EXTERNAL ZABS
      DOUBLE PRECISION AR, AI, BR, BI, ZM, DTHETA, DPI, DFRAC_PI_2
      DOUBLE PRECISION ZABS
      INTEGER IERR
      DATA DPI , DFRAC_PI_2  / 3.141592653589793238462643383e+0,
     1                   1.570796326794896619231321696e+0/
//
      IERR=0
      if (AR == 0.0e+0) GO TO 10
      if (AI == 0.0e+0) GO TO 20
      DTHETA = DATAN(AI/AR)
      if (DTHETA <= 0.0e+0) GO TO 40
      if (AR < 0.0e+0) DTHETA = DTHETA - DPI
      GO TO 50
   10 if (AI == 0.0e+0) GO TO 60
      BI = DFRAC_PI_2
      BR = DLOG((AI).abs())
      if (AI < 0.0e+0) BI = -BI
      RETURN
   20 if (AR > 0.0e+0) GO TO 30
      BR = DLOG((AR).abs())
      BI = DPI
      RETURN
   30 BR = DLOG(AR)
      BI = 0.0e+0
      RETURN
   40 if (AR < 0.0e+0) DTHETA = DTHETA + DPI
   50 ZM = ZABS(AR,AI)
      BR = DLOG(ZM)
      BI = DTHETA
      RETURN
   60 CONTINUE
      IERR=1
      RETURN
      END
      DOUBLE PRECISION FUNCTION ZABS(ZR, ZI)
// ***BEGIN PROLOGUE  ZABS
// ***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY
//
//     ZABS COMPUTES THE ABSOLUTE VALUE OR MAGNITUDE OF A DOUBLE
//     PRECISION COMPLEX VARIABLE CMPLX(ZR,ZI)
//
// ***ROUTINES CALLED  (NONE)
// ***END PROLOGUE  ZABS
      DOUBLE PRECISION ZR, ZI, U, V, Q, S
      U = (ZR).abs()
      V = (ZI).abs()
      S = U + V
//-----------------------------------------------------------------------
//     S*1.0 MAKES AN UNNORMALIZED UNDERFLOW ON CDC MACHINES INTO A
//     TRUE FLOATING ZERO
//-----------------------------------------------------------------------
      S = S*1.0e+0
      if (S == 0.0e+0) GO TO 20
      if (U > V) GO TO 10
      Q = U/V
      ZABS = V*DSQRT(1.D+0+Q*Q)
      RETURN
   10 Q = V/U
      ZABS = U*DSQRT(1.D+0+Q*Q)
      RETURN
   20 ZABS = 0.0e+0
      RETURN
      END
fn ZBKNU(ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL, ELIM,
     * ALIM)
// ***BEGIN PROLOGUE  ZBKNU
// ***REFER TO  ZBESI,ZBESK,ZAIRY,ZBESH
//
//     ZBKNU COMPUTES THE K BESSEL FUNCTION IN THE RIGHT HALF Z PLANE.
//
// ***ROUTINES CALLED  gamma_ln,i1mach,d1mach,ZKSCL,ZSHCH,ZUunderflowCHK,ZABS,ZDIV,
//                    ZEXP,ZLOG,ZMLT,ZSQRT
// ***END PROLOGUE  ZBKNU
//
      EXTERNAL ZABS
      DOUBLE PRECISION AA, AK, ALIM, ASCLE, A1, A2, BB, BK, BRY, CAZ,
     * CBI, CBR, CC, CCHI, CCHR, CKI, CKR, COEFI, COEFR, CONEI, CONER,
     * CRSCR, CSCLR, CSHI, CSHR, CSI, CSR, CSRR, CSSR, CTWOR,
     * CZEROI, CZEROR, CZI, CZR, DNU, DNU2, DPI, ELIM, ETEST, FC, FHS,
     * FI, FK, FKS, FMUI, FMUR, FNU, FPI, FR, G1, G2, FRAC_PI_2, PI, PR, PTI,
     * PTR, P1I, P1R, P2I, P2M, P2R, QI, QR, RAK, RCAZ, RTFRAC_PI_2, RZI,
     * RZR, R1, S, SMUI, SMUR, SPI, STI, STR, S1I, S1R, S2I, S2R, TM,
     * TOL, TTH, T1, T2, YI, YR, ZI, ZR, gamma_ln, d1mach, ZABS, ELM,
     * CELMR, ZDR, ZDI, AS, ALAS, HELIM, CYR, CYI
      INTEGER I, IFLAG, INU, K, KFLAG, KK, KMAX, KODE, KODED, N, NZ,
     * IDUM, i1mach, J, IC, INUB, NW
      DIMENSION YR(N), YI(N), CC(8), CSSR(3), CSRR(3), BRY(3), CYR(2),
     * CYI(2)
//     COMPLEX Z,Y,A,B,RZ,SMU,FU,FMU,F,FLRZ,CZ,S1,S2,CSH,CCH
//     COMPLEX CK,P,Q,COEF,P1,P2,CBK,PT,CZERO,CONE,CTWO,ST,EZ,CS,DK
//
      DATA KMAX / 30 /
      DATA CZEROR,CZEROI,CONER,CONEI,CTWOR,R1/
     1  0.0 , 0.0 , 1.0 , 0.0 , 2.0 , 2.0 /
      DATA DPI, RTFRAC_PI_2, SPI ,FRAC_PI_2, FPI, TTH /
     1     3.14159265358979324,       1.25331413731550025,
     2     1.90985931710274403,       1.57079632679489662,
     3     1.89769999331517738,       6.66666666666666666e-01/
      DATA CC(1), CC(2), CC(3), CC(4), CC(5), CC(6), CC(7), CC(8)/
     1     5.77215664901532861e-01,    -4.20026350340952355e-02,
     2    -4.21977345555443367e-02,     7.21894324666309954e-03,
     3    -2.15241674114950973e-04,    -2.01348547807882387e-05,
     4     1.13302723198169588e-06,     6.11609510448141582e-09/
//
      CAZ = ZABS(ZR,ZI)
      CSCLR = 1.0/TOL
      CRSCR = TOL
      CSSR(1) = CSCLR
      CSSR(2) = 1.0
      CSSR(3) = CRSCR
      CSRR(1) = CRSCR
      CSRR(2) = 1.0
      CSRR(3) = CSCLR
      BRY(1) = 1.0e+3*d1mach(1)/TOL
      BRY(2) = 1.0/BRY(1)
      BRY(3) = d1mach(2)
      NZ = 0
      IFLAG = 0
      KODED = KODE
      RCAZ = 1.0/CAZ
      STR = ZR*RCAZ
      STI = -ZI*RCAZ
      RZR = (STR+STR)*RCAZ
      RZI = (STI+STI)*RCAZ
      INU = INT(FNU+0.5)
      DNU = FNU - (INU as f64)
      if ((DNU).abs() == 0.5) GO TO 110
      DNU2 = 0.0
      if ((DNU).abs() > TOL) DNU2 = DNU*DNU
      if (CAZ > R1) GO TO 110
//-----------------------------------------------------------------------
//     SERIES FOR CABS(Z) <= R1
//-----------------------------------------------------------------------
      FC = 1.0
      CALL ZLOG(RZR, RZI, SMUR, SMUI, IDUM)
      FMUR = SMUR*DNU
      FMUI = SMUI*DNU
      CALL ZSHCH(FMUR, FMUI, CSHR, CSHI, CCHR, CCHI)
      if (DNU == 0.0) GO TO 10
      FC = DNU*DPI
      FC = FC/DSIN(FC)
      SMUR = CSHR/DNU
      SMUI = CSHI/DNU
   10 CONTINUE
      A2 = 1.0 + DNU
//-----------------------------------------------------------------------
//     GAM(1-Z)*GAM(1+Z)=PI*Z/SIN(PI*Z), T1=1/GAM(1-DNU), T2=1/GAM(1+DNU)
//-----------------------------------------------------------------------
      T2 = DEXP(-gamma_ln(A2,IDUM))
      T1 = 1.0/(T2*FC)
      if ((DNU).abs() > 0.1) GO TO 40
//-----------------------------------------------------------------------
//     SERIES FOR F0 TO RESOLVE INDETERMINACY FOR SMALL ABS(DNU)
//-----------------------------------------------------------------------
      AK = 1.0
      S = CC(1)
      DO 20 K=2,8
        AK = AK*DNU2
        TM = CC(K)*AK
        S = S + TM
        if ((TM).abs() < TOL) GO TO 30
   20 CONTINUE
   30 G1 = -S
      GO TO 50
   40 CONTINUE
      G1 = (T1-T2)/(DNU+DNU)
   50 CONTINUE
      G2 = (T1+T2)*0.5
      FR = FC*(CCHR*G1+SMUR*G2)
      FI = FC*(CCHI*G1+SMUI*G2)
      CALL ZEXP(FMUR, FMUI, STR, STI)
      PR = 0.5*STR/T2
      PI = 0.5*STI/T2
      CALL ZDIV(0.5, 0.0, STR, STI, PTR, PTI)
      QR = PTR/T1
      QI = PTI/T1
      S1R = FR
      S1I = FI
      S2R = PR
      S2I = PI
      AK = 1.0
      A1 = 1.0
      CKR = CONER
      CKI = CONEI
      BK = 1.0 - DNU2
      if (INU > 0 || N > 1) GO TO 80
//-----------------------------------------------------------------------
//     GENERATE K(FNU,Z), 0.0  <=  FNU  <  0.5 AND N=1
//-----------------------------------------------------------------------
      if (CAZ < TOL) GO TO 70
      CALL ZMLT(ZR, ZI, ZR, ZI, CZR, CZI)
      CZR = 0.25*CZR
      CZI = 0.25*CZI
      T1 = 0.25*CAZ*CAZ
   60 CONTINUE
      FR = (FR*AK+PR+QR)/BK
      FI = (FI*AK+PI+QI)/BK
      STR = 1.0/(AK-DNU)
      PR = PR*STR
      PI = PI*STR
      STR = 1.0/(AK+DNU)
      QR = QR*STR
      QI = QI*STR
      STR = CKR*CZR - CKI*CZI
      RAK = 1.0/AK
      CKI = (CKR*CZI+CKI*CZR)*RAK
      CKR = STR*RAK
      S1R = CKR*FR - CKI*FI + S1R
      S1I = CKR*FI + CKI*FR + S1I
      A1 = A1*T1*RAK
      BK = BK + AK + AK + 1.0
      AK = AK + 1.0
      if (A1 > TOL) GO TO 60
   70 CONTINUE
      YR(1) = S1R
      YI(1) = S1I
      if (KODED == 1) RETURN
      CALL ZEXP(ZR, ZI, STR, STI)
      CALL ZMLT(S1R, S1I, STR, STI, YR(1), YI(1))
      RETURN
//-----------------------------------------------------------------------
//     GENERATE K(DNU,Z) AND K(DNU+1,Z) FOR FORWARD RECURRENCE
//-----------------------------------------------------------------------
   80 CONTINUE
      if (CAZ < TOL) GO TO 100
      CALL ZMLT(ZR, ZI, ZR, ZI, CZR, CZI)
      CZR = 0.25*CZR
      CZI = 0.25*CZI
      T1 = 0.25*CAZ*CAZ
   90 CONTINUE
      FR = (FR*AK+PR+QR)/BK
      FI = (FI*AK+PI+QI)/BK
      STR = 1.0/(AK-DNU)
      PR = PR*STR
      PI = PI*STR
      STR = 1.0/(AK+DNU)
      QR = QR*STR
      QI = QI*STR
      STR = CKR*CZR - CKI*CZI
      RAK = 1.0/AK
      CKI = (CKR*CZI+CKI*CZR)*RAK
      CKR = STR*RAK
      S1R = CKR*FR - CKI*FI + S1R
      S1I = CKR*FI + CKI*FR + S1I
      STR = PR - FR*AK
      STI = PI - FI*AK
      S2R = CKR*STR - CKI*STI + S2R
      S2I = CKR*STI + CKI*STR + S2I
      A1 = A1*T1*RAK
      BK = BK + AK + AK + 1.0
      AK = AK + 1.0
      if (A1 > TOL) GO TO 90
  100 CONTINUE
      KFLAG = 2
      A1 = FNU + 1.0
      AK = A1*(SMUR).abs()
      if (AK > ALIM) KFLAG = 3
      STR = CSSR(KFLAG)
      P2R = S2R*STR
      P2I = S2I*STR
      CALL ZMLT(P2R, P2I, RZR, RZI, S2R, S2I)
      S1R = S1R*STR
      S1I = S1I*STR
      if (KODED == 1) GO TO 210
      CALL ZEXP(ZR, ZI, FR, FI)
      CALL ZMLT(S1R, S1I, FR, FI, S1R, S1I)
      CALL ZMLT(S2R, S2I, FR, FI, S2R, S2I)
      GO TO 210
//-----------------------------------------------------------------------
//     IFLAG=0 MEANS NO UNDERFLOW OCCURRED
//     IFLAG=1 MEANS AN UNDERFLOW OCCURRED- COMPUTATION PROCEEDS WITH
//     KODED=2 AND A TEST FOR ON SCALE VALUES IS MADE DURING FORWARD
//     RECURSION
//-----------------------------------------------------------------------
  110 CONTINUE
      CALL ZSQRT(ZR, ZI, STR, STI)
      CALL ZDIV(RTFRAC_PI_2, CZEROI, STR, STI, COEFR, COEFI)
      KFLAG = 2
      if (KODED == 2) GO TO 120
      if (ZR > ALIM) GO TO 290
//     BLANK LINE
      STR = DEXP(-ZR)*CSSR(KFLAG)
      STI = -STR*DSIN(ZI)
      STR = STR*DCOS(ZI)
      CALL ZMLT(COEFR, COEFI, STR, STI, COEFR, COEFI)
  120 CONTINUE
      if ((DNU).abs() == 0.5) GO TO 300
//-----------------------------------------------------------------------
//     MILLER ALGORITHM FOR CABS(Z) > R1
//-----------------------------------------------------------------------
      AK = DCOS(DPI*DNU)
      AK = (AK).abs()
      if (AK == CZEROR) GO TO 300
      FHS = (0.25-DNU2).abs()
      if (FHS == CZEROR) GO TO 300
//-----------------------------------------------------------------------
//     COMPUTE R2=F(E). if CABS(Z) >= R2, USE FORWARD RECURRENCE TO
//     DETERMINE THE BACKWARD INDEX K. R2=F(E) IS A STRAIGHT LINE ON
//     12 <= E <= 60. E IS COMPUTED FROM 2**(-E)=B**(1-i1mach(14))=
//     TOL WHERE B IS THE BASE OF THE ARITHMETIC.
//-----------------------------------------------------------------------
      T1 = DBLE(FLOAT(i1mach(14)-1))
      T1 = T1*d1mach(5)*3.321928094
      T1 = DMAX1(T1,12.0)
      T1 = DMIN1(T1,60.0)
      T2 = TTH*T1 - 6.0
      if (ZR != 0.0) GO TO 130
      T1 = FRAC_PI_2
      GO TO 140
  130 CONTINUE
      T1 = DATAN(ZI/ZR)
      T1 = (T1).abs()
  140 CONTINUE
      if (T2 > CAZ) GO TO 170
//-----------------------------------------------------------------------
//     FORWARD RECURRENCE LOOP WHEN CABS(Z) >= R2
//-----------------------------------------------------------------------
      ETEST = AK/(DPI*CAZ*TOL)
      FK = CONER
      if (ETEST < CONER) GO TO 180
      FKS = CTWOR
      CKR = CAZ + CAZ + CTWOR
      P1R = CZEROR
      P2R = CONER
      DO 150 I=1,KMAX
        AK = FHS/FKS
        CBR = CKR/(FK+CONER)
        PTR = P2R
        P2R = CBR*P2R - P1R*AK
        P1R = PTR
        CKR = CKR + CTWOR
        FKS = FKS + FK + FK + CTWOR
        FHS = FHS + FK + FK
        FK = FK + CONER
        STR = (P2R).abs()*FK
        if (ETEST < STR) GO TO 160
  150 CONTINUE
      GO TO 310
  160 CONTINUE
      FK = FK + SPI*T1*DSQRT(T2/CAZ)
      FHS = (0.25-DNU2).abs()
      GO TO 180
  170 CONTINUE
//-----------------------------------------------------------------------
//     COMPUTE BACKWARD INDEX K FOR CABS(Z) < R2
//-----------------------------------------------------------------------
      A2 = DSQRT(CAZ)
      AK = FPI*AK/(TOL*DSQRT(A2))
      AA = 3.0*T1/(1.0+CAZ)
      BB = 14.7*T1/(28.0+CAZ)
      AK = (DLOG(AK)+CAZ*DCOS(AA)/(1.0+0.008*CAZ))/DCOS(BB)
      FK = 0.12125*AK*AK/CAZ + 1.5
  180 CONTINUE
//-----------------------------------------------------------------------
//     BACKWARD RECURRENCE LOOP FOR MILLER ALGORITHM
//-----------------------------------------------------------------------
      K = INT(SNGL(FK))
      FK = (K as f64)
      FKS = FK*FK
      P1R = CZEROR
      P1I = CZEROI
      P2R = TOL
      P2I = CZEROI
      CSR = P2R
      CSI = P2I
      DO 190 I=1,K
        A1 = FKS - FK
        AK = (FKS+FK)/(A1+FHS)
        RAK = 2.0/(FK+CONER)
        CBR = (FK+ZR)*RAK
        CBI = ZI*RAK
        PTR = P2R
        PTI = P2I
        P2R = (PTR*CBR-PTI*CBI-P1R)*AK
        P2I = (PTI*CBR+PTR*CBI-P1I)*AK
        P1R = PTR
        P1I = PTI
        CSR = CSR + P2R
        CSI = CSI + P2I
        FKS = A1 - FK + CONER
        FK = FK - CONER
  190 CONTINUE
//-----------------------------------------------------------------------
//     COMPUTE (P2/CS)=(P2/CABS(CS))*(CONJG(CS)/CABS(CS)) FOR BETTER
//     SCALING
//-----------------------------------------------------------------------
      TM = ZABS(CSR,CSI)
      PTR = 1.0/TM
      S1R = P2R*PTR
      S1I = P2I*PTR
      CSR = CSR*PTR
      CSI = -CSI*PTR
      CALL ZMLT(COEFR, COEFI, S1R, S1I, STR, STI)
      CALL ZMLT(STR, STI, CSR, CSI, S1R, S1I)
      if (INU > 0 || N > 1) GO TO 200
      ZDR = ZR
      ZDI = ZI
      if(IFLAG == 1) GO TO 270
      GO TO 240
  200 CONTINUE
//-----------------------------------------------------------------------
//     COMPUTE P1/P2=(P1/CABS(P2)*CONJG(P2)/CABS(P2) FOR SCALING
//-----------------------------------------------------------------------
      TM = ZABS(P2R,P2I)
      PTR = 1.0/TM
      P1R = P1R*PTR
      P1I = P1I*PTR
      P2R = P2R*PTR
      P2I = -P2I*PTR
      CALL ZMLT(P1R, P1I, P2R, P2I, PTR, PTI)
      STR = DNU + 0.5 - PTR
      STI = -PTI
      CALL ZDIV(STR, STI, ZR, ZI, STR, STI)
      STR = STR + 1.0
      CALL ZMLT(STR, STI, S1R, S1I, S2R, S2I)
//-----------------------------------------------------------------------
//     FORWARD RECURSION ON THE THREE TERM RECURSION WITH RELATION WITH
//     SCALING NEAR EXPONENT EXTREMES ON KFLAG=1 OR KFLAG=3
//-----------------------------------------------------------------------
  210 CONTINUE
      STR = DNU + 1.0
      CKR = STR*RZR
      CKI = STR*RZI
      if (N == 1) INU = INU - 1
      if (INU > 0) GO TO 220
      if (N > 1) GO TO 215
      S1R = S2R
      S1I = S2I
  215 CONTINUE
      ZDR = ZR
      ZDI = ZI
      if(IFLAG == 1) GO TO 270
      GO TO 240
  220 CONTINUE
      INUB = 1
      if(IFLAG == 1) GO TO 261
  225 CONTINUE
      P1R = CSRR(KFLAG)
      ASCLE = BRY(KFLAG)
      DO 230 I=INUB,INU
        STR = S2R
        STI = S2I
        S2R = CKR*STR - CKI*STI + S1R
        S2I = CKR*STI + CKI*STR + S1I
        S1R = STR
        S1I = STI
        CKR = CKR + RZR
        CKI = CKI + RZI
        if (KFLAG >= 3) GO TO 230
        P2R = S2R*P1R
        P2I = S2I*P1R
        STR = (P2R).abs()
        STI = (P2I).abs()
        P2M = DMAX1(STR,STI)
        if (P2M <= ASCLE) GO TO 230
        KFLAG = KFLAG + 1
        ASCLE = BRY(KFLAG)
        S1R = S1R*P1R
        S1I = S1I*P1R
        S2R = P2R
        S2I = P2I
        STR = CSSR(KFLAG)
        S1R = S1R*STR
        S1I = S1I*STR
        S2R = S2R*STR
        S2I = S2I*STR
        P1R = CSRR(KFLAG)
  230 CONTINUE
      if (N != 1) GO TO 240
      S1R = S2R
      S1I = S2I
  240 CONTINUE
      STR = CSRR(KFLAG)
      YR(1) = S1R*STR
      YI(1) = S1I*STR
      if (N == 1) RETURN
      YR(2) = S2R*STR
      YI(2) = S2I*STR
      if (N == 2) RETURN
      KK = 2
  250 CONTINUE
      KK = KK + 1
      if (KK > N) RETURN
      P1R = CSRR(KFLAG)
      ASCLE = BRY(KFLAG)
      DO 260 I=KK,N
        P2R = S2R
        P2I = S2I
        S2R = CKR*P2R - CKI*P2I + S1R
        S2I = CKI*P2R + CKR*P2I + S1I
        S1R = P2R
        S1I = P2I
        CKR = CKR + RZR
        CKI = CKI + RZI
        P2R = S2R*P1R
        P2I = S2I*P1R
        YR(I) = P2R
        YI(I) = P2I
        if (KFLAG >= 3) {return Err(LossOfSignificance);}
        STR = (P2R).abs()
        STI = (P2I).abs()
        P2M = DMAX1(STR,STI)
        if (P2M <= ASCLE) GO TO 260
        KFLAG = KFLAG + 1
        ASCLE = BRY(KFLAG)
        S1R = S1R*P1R
        S1I = S1I*P1R
        S2R = P2R
        S2I = P2I
        STR = CSSR(KFLAG)
        S1R = S1R*STR
        S1I = S1I*STR
        S2R = S2R*STR
        S2I = S2I*STR
        P1R = CSRR(KFLAG)
  260 CONTINUE
      RETURN
//-----------------------------------------------------------------------
//     IFLAG=1 CASES, FORWARD RECURRENCE ON SCALED VALUES ON UNDERFLOW
//-----------------------------------------------------------------------
  261 CONTINUE
      HELIM = 0.5*ELIM
      ELM = DEXP(-ELIM)
      CELMR = ELM
      ASCLE = BRY(1)
      ZDR = ZR
      ZDI = ZI
      IC = -1
      J = 2
      DO 262 I=1,INU
        STR = S2R
        STI = S2I
        S2R = STR*CKR-STI*CKI+S1R
        S2I = STI*CKR+STR*CKI+S1I
        S1R = STR
        S1I = STI
        CKR = CKR+RZR
        CKI = CKI+RZI
        AS = ZABS(S2R,S2I)
        ALAS = DLOG(AS)
        P2R = -ZDR+ALAS
        if(P2R < (-ELIM)) GO TO 263
        CALL ZLOG(S2R,S2I,STR,STI,IDUM)
        P2R = -ZDR+STR
        P2I = -ZDI+STI
        P2M = DEXP(P2R)/TOL
        P1R = P2M*DCOS(P2I)
        P1I = P2M*DSIN(P2I)
        CALL ZUunderflowCHK(P1R,P1I,NW,ASCLE,TOL)
        if(NW != 0) GO TO 263
        J = 3 - J
        CYR(J) = P1R
        CYI(J) = P1I
        if(IC == (I-1)) GO TO 264
        IC = I
        GO TO 262
  263   CONTINUE
        if(ALAS < HELIM) GO TO 262
        ZDR = ZDR-ELIM
        S1R = S1R*CELMR
        S1I = S1I*CELMR
        S2R = S2R*CELMR
        S2I = S2I*CELMR
  262 CONTINUE
      if(N != 1) GO TO 270
      S1R = S2R
      S1I = S2I
      GO TO 270
  264 CONTINUE
      KFLAG = 1
      INUB = I+1
      S2R = CYR(J)
      S2I = CYI(J)
      J = 3 - J
      S1R = CYR(J)
      S1I = CYI(J)
      if(INUB <= INU) GO TO 225
      if(N != 1) GO TO 240
      S1R = S2R
      S1I = S2I
      GO TO 240
  270 CONTINUE
      YR(1) = S1R
      YI(1) = S1I
      if(N == 1) GO TO 280
      YR(2) = S2R
      YI(2) = S2I
  280 CONTINUE
      ASCLE = BRY(1)
      CALL ZKSCL(ZDR,ZDI,FNU,N,YR,YI,NZ,RZR,RZI,ASCLE,TOL,ELIM)
      INU = N - NZ
      if (INU <= 0) RETURN
      KK = NZ + 1
      S1R = YR(KK)
      S1I = YI(KK)
      YR(KK) = S1R*CSRR(1)
      YI(KK) = S1I*CSRR(1)
      if (INU == 1) RETURN
      KK = NZ + 2
      S2R = YR(KK)
      S2I = YI(KK)
      YR(KK) = S2R*CSRR(1)
      YI(KK) = S2I*CSRR(1)
      if (INU == 2) RETURN
      T2 = FNU + ((KK-1) as f64)
      CKR = T2*RZR
      CKI = T2*RZI
      KFLAG = 1
      GO TO 250
  290 CONTINUE
//-----------------------------------------------------------------------
//     SCALE BY DEXP(Z), IFLAG = 1 CASES
//-----------------------------------------------------------------------
      KODED = 2
      IFLAG = 1
      KFLAG = 2
      GO TO 120
//-----------------------------------------------------------------------
//     FNU=HALF ODD INTEGER CASE, DNU=-0.5
//-----------------------------------------------------------------------
  300 CONTINUE
      S1R = COEFR
      S1I = COEFI
      S2R = COEFR
      S2I = COEFI
      GO TO 210
//
//
  310 CONTINUE
      NZ=-2
      RETURN
      END
fn ZKSCL(ZRR,ZRI,FNU,N,YR,YI,NZ,RZR,RZI,ASCLE,TOL,ELIM)
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
      EXTERNAL ZABS
      DOUBLE PRECISION ACS, AS, ASCLE, CKI, CKR, CSI, CSR, CYI,
     * CYR, ELIM, FN, FNU, RZI, RZR, STR, S1I, S1R, S2I,
     * S2R, TOL, YI, YR, ZEROI, ZEROR, ZRI, ZRR, ZABS,
     * ZDR, ZDI, CELMR, ELM, HELIM, ALAS
      INTEGER I, IC, IDUM, KK, N, NN, NW, NZ
      DIMENSION YR(N), YI(N), CYR(2), CYI(2)
      DATA ZEROR,ZEROI / 0.0 , 0.0 /
//
      NZ = 0
      IC = 0
      NN = MIN0(2,N)
      DO 10 I=1,NN
        S1R = YR(I)
        S1I = YI(I)
        CYR(I) = S1R
        CYI(I) = S1I
        AS = ZABS(S1R,S1I)
        ACS = -ZRR + DLOG(AS)
        NZ = NZ + 1
        YR(I) = ZEROR
        YI(I) = ZEROI
        if (ACS < (-ELIM)) GO TO 10
        CALL ZLOG(S1R, S1I, CSR, CSI, IDUM)
        CSR = CSR - ZRR
        CSI = CSI - ZRI
        STR = DEXP(CSR)/TOL
        CSR = STR*DCOS(CSI)
        CSI = STR*DSIN(CSI)
        CALL ZUunderflowCHK(CSR, CSI, NW, ASCLE, TOL)
        if (NW != 0) GO TO 10
        YR(I) = CSR
        YI(I) = CSI
        IC = I
        NZ = NZ - 1
   10 CONTINUE
      if (N == 1) RETURN
      if (IC > 1) GO TO 20
      YR(1) = ZEROR
      YI(1) = ZEROI
      NZ = 2
   20 CONTINUE
      if (N == 2) RETURN
      if (NZ == 0) RETURN
      FN = FNU + 1.0
      CKR = FN*RZR
      CKI = FN*RZI
      S1R = CYR(1)
      S1I = CYI(1)
      S2R = CYR(2)
      S2I = CYI(2)
      HELIM = 0.5*ELIM
      ELM = DEXP(-ELIM)
      CELMR = ELM
      ZDR = ZRR
      ZDI = ZRI
//
//     FIND TWO CONSECUTIVE Y VALUES ON SCALE. SCALE RECURRENCE if
//     S2 GETS LARGER THAN EXP(ELIM/2)
//
      DO 30 I=3,N
        KK = I
        CSR = S2R
        CSI = S2I
        S2R = CKR*CSR - CKI*CSI + S1R
        S2I = CKI*CSR + CKR*CSI + S1I
        S1R = CSR
        S1I = CSI
        CKR = CKR + RZR
        CKI = CKI + RZI
        AS = ZABS(S2R,S2I)
        ALAS = DLOG(AS)
        ACS = -ZDR + ALAS
        NZ = NZ + 1
        YR(I) = ZEROR
        YI(I) = ZEROI
        if (ACS < (-ELIM)) GO TO 25
        CALL ZLOG(S2R, S2I, CSR, CSI, IDUM)
        CSR = CSR - ZDR
        CSI = CSI - ZDI
        STR = DEXP(CSR)/TOL
        CSR = STR*DCOS(CSI)
        CSI = STR*DSIN(CSI)
        CALL ZUunderflowCHK(CSR, CSI, NW, ASCLE, TOL)
        if (NW != 0) GO TO 25
        YR(I) = CSR
        YI(I) = CSI
        NZ = NZ - 1
        if (IC == KK-1) GO TO 40
        IC = KK
        GO TO 30
   25   CONTINUE
        if(ALAS < HELIM) GO TO 30
        ZDR = ZDR - ELIM
        S1R = S1R*CELMR
        S1I = S1I*CELMR
        S2R = S2R*CELMR
        S2I = S2I*CELMR
   30 CONTINUE
      NZ = N
      if(IC == N) NZ=N-1
      GO TO 45
   40 CONTINUE
      NZ = KK - 2
   45 CONTINUE
      DO 50 I=1,NZ
        YR(I) = ZEROR
        YI(I) = ZEROI
   50 CONTINUE
      RETURN
      END
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
fn ZRATI(ZR, ZI, FNU, N, CYR, CYI, TOL)
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
//     COMPLEX Z,CY(1),CONE,CZERO,P1,P2,T1,RZ,PT,CDFNU
      EXTERNAL ZABS
      DOUBLE PRECISION AK, AMAGZ, AP1, AP2, ARG, AZ, CDFNUI, CDFNUR,
     * CONEI, CONER, CYI, CYR, CZEROI, CZEROR, DFNU, FDNU, FLAM, FNU,
     * FNUP, PTI, PTR, P1I, P1R, P2I, P2R, RAK, RAP1, RHO, RT2, RZI,
     * RZR, TEST, TEST1, TOL, TTI, TTR, T1I, T1R, ZI, ZR, ZABS
      INTEGER I, ID, IDNU, INU, ITIME, K, KK, MAGZ, N
      DIMENSION CYR(N), CYI(N)
      DATA CZEROR,CZEROI,CONER,CONEI,RT2/
     1 0.0, 0.0, 1.0, 0.0, 1.41421356237309505 /
      AZ = ZABS(ZR,ZI)
      INU = INT(SNGL(FNU))
      IDNU = INU + N - 1
      MAGZ = INT(SNGL(AZ))
      AMAGZ = DBLE(FLOAT(MAGZ+1))
      FDNU = (IDNU as f64)
      FNUP = DMAX1(AMAGZ,FDNU)
      ID = IDNU - MAGZ - 1
      ITIME = 1
      K = 1
      PTR = 1.0/AZ
      RZR = PTR*(ZR+ZR)*PTR
      RZI = -PTR*(ZI+ZI)*PTR
      T1R = RZR*FNUP
      T1I = RZI*FNUP
      P2R = -T1R
      P2I = -T1I
      P1R = CONER
      P1I = CONEI
      T1R = T1R + RZR
      T1I = T1I + RZI
      if (ID > 0) ID = 0
      AP2 = ZABS(P2R,P2I)
      AP1 = ZABS(P1R,P1I)
//-----------------------------------------------------------------------
//     THE OVERFLOW TEST ON K(FNU+I-1,Z) BEFORE THE CALL TO CBKNU
//     GUARANTEES THAT P2 IS ON SCALE. SCALE TEST1 AND ALL SUBSEQUENT
//     P2 VALUES BY AP1 TO ENSURE THAT AN OVERFLOW DOES NOT OCCUR
//     PREMATURELY.
//-----------------------------------------------------------------------
      ARG = (AP2+AP2)/(AP1*TOL)
      TEST1 = DSQRT(ARG)
      TEST = TEST1
      RAP1 = 1.0/AP1
      P1R = P1R*RAP1
      P1I = P1I*RAP1
      P2R = P2R*RAP1
      P2I = P2I*RAP1
      AP2 = AP2*RAP1
   10 CONTINUE
      K = K + 1
      AP1 = AP2
      PTR = P2R
      PTI = P2I
      P2R = P1R - (T1R*PTR-T1I*PTI)
      P2I = P1I - (T1R*PTI+T1I*PTR)
      P1R = PTR
      P1I = PTI
      T1R = T1R + RZR
      T1I = T1I + RZI
      AP2 = ZABS(P2R,P2I)
      if (AP1 <= TEST) GO TO 10
      if (ITIME == 2) GO TO 20
      AK = ZABS(T1R,T1I)*0.5
      FLAM = AK + DSQRT(AK*AK-1.0)
      RHO = DMIN1(AP2/AP1,FLAM)
      TEST = TEST1*DSQRT(RHO/(RHO*RHO-1.0))
      ITIME = 2
      GO TO 10
   20 CONTINUE
      KK = K + 1 - ID
      AK = (KK as f64)
      T1R = AK
      T1I = CZEROI
      DFNU = FNU + ((N-1) as f64)
      P1R = 1.0/AP2
      P1I = CZEROI
      P2R = CZEROR
      P2I = CZEROI
      DO 30 I=1,KK
        PTR = P1R
        PTI = P1I
        RAP1 = DFNU + T1R
        TTR = RZR*RAP1
        TTI = RZI*RAP1
        P1R = (PTR*TTR-PTI*TTI) + P2R
        P1I = (PTR*TTI+PTI*TTR) + P2I
        P2R = PTR
        P2I = PTI
        T1R = T1R - CONER
   30 CONTINUE
      if (P1R != CZEROR || P1I != CZEROI) GO TO 40
      P1R = TOL
      P1I = TOL
   40 CONTINUE
      CALL ZDIV(P2R, P2I, P1R, P1I, CYR(N), CYI(N))
      if (N == 1) RETURN
      K = N - 1
      AK = (K as f64)
      T1R = AK
      T1I = CZEROI
      CDFNUR = FNU*RZR
      CDFNUI = FNU*RZI
      DO 60 I=2,N
        PTR = CDFNUR + (T1R*RZR-T1I*RZI) + CYR(K+1)
        PTI = CDFNUI + (T1R*RZI+T1I*RZR) + CYI(K+1)
        AK = ZABS(PTR,PTI)
        if (AK != CZEROR) GO TO 50
        PTR = TOL
        PTI = TOL
        AK = TOL*RT2
   50   CONTINUE
        RAK = CONER/AK
        CYR(K) = RAK*PTR*RAK
        CYI(K) = -RAK*PTI*RAK
        T1R = T1R - CONER
        K = K - 1
   60 CONTINUE
      RETURN
      END
fn ZS1S2(ZRR, ZRI, S1R, S1I, S2R, S2I, NZ, ASCLE, ALIM,
     * IUF)
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
      EXTERNAL ZABS
      DOUBLE PRECISION AA, ALIM, ALN, ASCLE, AS1, AS2, C1I, C1R, S1DI,
     * S1DR, S1I, S1R, S2I, S2R, ZEROI, ZEROR, ZRI, ZRR, ZABS
      INTEGER IUF, IDUM, NZ
      DATA ZEROR,ZEROI  / 0.0 , 0.0 /
      NZ = 0
      AS1 = ZABS(S1R,S1I)
      AS2 = ZABS(S2R,S2I)
      if (S1R == 0.0 && S1I == 0.0) GO TO 10
      if (AS1 == 0.0) GO TO 10
      ALN = -ZRR - ZRR + DLOG(AS1)
      S1DR = S1R
      S1DI = S1I
      S1R = ZEROR
      S1I = ZEROI
      AS1 = ZEROR
      if (ALN < (-ALIM)) GO TO 10
      CALL ZLOG(S1DR, S1DI, C1R, C1I, IDUM)
      C1R = C1R - ZRR - ZRR
      C1I = C1I - ZRI - ZRI
      CALL ZEXP(C1R, C1I, S1R, S1I)
      AS1 = ZABS(S1R,S1I)
      IUF = IUF + 1
   10 CONTINUE
      AA = DMAX1(AS1,AS2)
      if (AA > ASCLE) RETURN
      S1R = ZEROR
      S1I = ZEROI
      S2R = ZEROR
      S2I = ZEROI
      NZ = 1
      IUF = 0
      RETURN
      END
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
fn ZMLRI(
    z: Complex64,
    order: f64, //ZR, ZI, FNU,
    KODE: Scaling,
    N: usize, //YR, YI, NZ,
    machine_consts: &MachineConsts,
) -> Result<(Vec<Complex64>, usize), BesselError> {
    // ***BEGIN PROLOGUE  ZMLRI
    // ***REFER TO  ZBESI,ZBESK
    //
    //     ZMLRI COMPUTES THE I BESSEL FUNCTION FOR RE(Z) >= 0.0 BY THE
    //     MILLER ALGORITHM NORMALIZED BY A NEUMANN SERIES.
    //
    // ***ROUTINES CALLED  gamma_ln,d1mach,ZABS,ZEXP,ZLOG,ZMLT
    // ***END PROLOGUE  ZMLRI
    // //     COMPLEX CK,CNORM,CONE,CTWO,CZERO,PT,P1,P2,RZ,SUM,Y,Z
    //       EXTERNAL ZABS
    //       DOUBLE PRECISION ACK, AK, AP, AT, AZ, BK, CKI, CKR, CNORMI,
    //      * CNORMR, CONEI, CONER, FKAP, FKK, FLAM, FNF, FNU, PTI, PTR, P1I,
    //      * P1R, P2I, P2R, RAZ, RHO, RHO2, RZI, RZR, SCLE, STI, STR, SUMI,
    //      * SUMR, TFNF, TOL, TST, YI, YR, ZEROI, ZEROR, ZI, ZR, gamma_ln,
    //      * d1mach, ZABS
    //       INTEGER I, IAZ, IDUM, IFNU, INU, ITIME, K, KK, KM, KODE, M, N, NZ
    //       DIMENSION YR(N), YI(N)
    //       DATA ZEROR,ZEROI,CONER,CONEI / 0.0, 0.0, 1.0, 0.0 /
    let SCLE: f64 = d1mach(1) / machine_consts.tol;
    let NZ = 0;
    let AZ = z.abs(); //ZABS(ZR,ZI);
    let IAZ = AZ as usize;
    let IFNU = order as usize; //INT(SNGL(FNU));
    let INU = IFNU + N - 1;
    let AT = (IAZ as f64) + 1.0;
    let RAZ = 1.0 / AZ;
    // STR = ZR*RAZ;
    // STI = -ZI*RAZ;
    // CKR = STR*AT*RAZ;
    // CKI = STI*AT*RAZ;
    let mut ck = z.conj() * RAZ * RAZ * AT;
    let rz = z.conj() * 2.0 * RAZ;
    // RZR = (STR+STR)*RAZ;
    // RZI = (STI+STI)*RAZ;
    let mut p1 = c_zero();
    let mut p2 = c_one();
    // P1R = ZEROR;
    // P1I = ZEROI;
    // P2R = CONER;
    // P2I = CONEI;
    let mut ACK = (AT + 1.0) * RAZ;
    let RHO = ACK + (ACK * ACK - 1.0).sqrt();
    let RHO2 = RHO * RHO;
    let mut TST = (RHO2 + RHO2) / ((RHO2 - 1.0) * (RHO - 1.0));
    TST /= machine_consts.tol;
    //-----------------------------------------------------------------------;
    //     COMPUTE RELATIVE TRUNCATION ERROR INDEX FOR SERIES;
    //-----------------------------------------------------------------------;
    let mut AK = AT;
    // DO 10 I=1,80;
    let mut converged = false;
    let mut I = 0;
    for i in 0..80 {
        I = i;
        let pt = p2;
        //   PTR = P2R;
        //   PTI = P2I;
        p2 = p1 - ck * p2;
        //   P2R = P1R - (CKR*PTR-CKI*PTI);
        //   P2I = P1I - (CKI*PTR+CKR*PTI);
        p1 = pt;
        //   P1R = PTR;
        //   P1I = PTI;
        ck += rz;
        //   CKR = CKR + RZR;
        //   CKI = CKI + RZI;
        //   AP = ZABS(P2R,P2I);
        if p2.abs() > TST * AK * AK {
            converged = true;
            break;
        }
        AK += 1.0;
    }
    //    10 CONTINUE;
    if !converged {
        return Err(DidNotConverge);
    }
    // GO TO 110;
    //    20 CONTINUE;
    I += 1;
    let mut K = 0;
    if INU >= IAZ {
        //GO TO 40;
        //-----------------------------------------------------------------------;
        //     COMPUTE RELATIVE TRUNCATION ERROR FOR RATIOS;
        //-----------------------------------------------------------------------;
        p1 = c_zero();
        p2 = c_one();
        // P1R = ZEROR;
        // P1I = ZEROI;
        // P2R = CONER;
        // P2I = CONEI;
        let AT = (INU as f64) + 1.0;
        // STR = ZR*RAZ;
        // STI = -ZI*RAZ;
        // CKR = STR*AT*RAZ;
        // CKI = STI*AT*RAZ;
        ck = z.conj() * RAZ * RAZ * AT;
        ACK = AT * RAZ;
        TST = (ACK / machine_consts.tol).sqrt();
        let mut hit_loop_end = false;
        // DO 30 K=1,80;
        converged = false;
        for k in 0..80 {
            K = k;
            let pt = p2;
            //   PTR = P2R;
            //   PTI = P2I;
            p2 = p1 - ck * pt;
            //   P2R = P1R - (CKR*PTR-CKI*PTI);
            //   P2I = P1I - (CKR*PTI+CKI*PTR);
            p1 = pt;
            //   P1R = PTR;
            //   P1I = PTI;
            ck += rz;
            //   CKR = CKR + RZR;
            //   CKI = CKI + RZI;
            //   AP = ZABS(P2R,P2I);
            let AP = p2.abs();
            if AP < TST {
                continue;
            }
            if hit_loop_end {
                converged = true;
                break;
            }
            ACK = ck.abs(); //ZABS(CKR,CKI);
            let FLAM = ACK + (ACK * ACK - 1.0).sqrt();
            let FKAP = AP / p1.abs(); //ZABS(P1R,P1I);
            let RHO = FLAM.min(FKAP); //DMIN1(FLAM,FKAP);
            TST *= (RHO / (RHO * RHO - 1.0)).sqrt();
            hit_loop_end = true;
        }
        //    30 CONTINUE;

        if !converged {
            return Err(DidNotConverge);
        }

        // GO TO 110;
    }
    //    40 CONTINUE;
    //-----------------------------------------------------------------------;
    //     BACKWARD RECURRENCE AND SUM NORMALIZING RELATION;
    //-----------------------------------------------------------------------;
    K += 1;
    // KK = MAX0(I+IAZ,K+INU);
    let KK = (I + IAZ).max(K + INU);
    let mut FKK = KK as f64;
    let mut p1 = c_zero();
    // P1R = ZEROR;
    // P1I = ZEROI;
    //-----------------------------------------------------------------------;
    //     SCALE P2 AND SUM BY SCLE;
    //-----------------------------------------------------------------------;
    let mut p2 = Complex64::new(SCLE, 0.0);
    // P2R = SCLE;
    // P2I = ZEROI;
    let FNF = order - (IFNU as f64);
    let TFNF = FNF + FNF;
    let mut BK = (gamma_ln(FKK + TFNF + 1.0).unwrap()
        - gamma_ln(FKK + 1.0).unwrap()
        - gamma_ln(TFNF + 1.0).unwrap())
    .exp();
    // BK = DEXP(BK);
    let mut sumr = c_zero();
    // SUMR = ZEROR;
    // SUMI = ZEROI;
    // let KM = (KK - INU);
    // DO 50 I=1,KM;
    for _i in 0..(KK - INU) {
        let pt = p2;
        //   PTR = P2R;
        //   PTI = P2I;
        p2 = p1 + (FKK + FNF) * (rz * pt);
        //   P2R = P1R + (FKK+FNF)*(RZR*PTR-RZI*PTI);
        //   P2I = P1I + (FKK+FNF)*(RZI*PTR+RZR*PTI);
        p1 = pt;
        //   P1R = PTR;
        //   P1I = PTI;
        AK = 1.0 - TFNF / (FKK + TFNF);
        ACK = BK * AK;
        sumr += (ACK + BK) * p1;
        //   SUMR = SUMR + (ACK+BK)*P1R;
        //   SUMI = SUMI + (ACK+BK)*P1I;
        BK = ACK;
        FKK = FKK - 1.0;
    }
    //    50 CONTINUE;
    let mut y = vec![c_zero(); N];
    y[N] = p2;
    // YR(N) = P2R;
    // YI(N) = P2I;
    if N != 1 {
        //GO TO 70;
        // DO 60 I=2,N;
        for i in 1..N {
            let pt = p2;
            p2 = p1 + (FKK + FNF) * (rz * pt);
            //   PTR = P2R;
            //   PTI = P2I;
            //   P2R = P1R + (FKK+FNF)*(RZR*PTR-RZI*PTI);
            //   P2I = P1I + (FKK+FNF)*(RZI*PTR+RZR*PTI);
            p1 = pt;
            //   P1R = PTR;
            //   P1I = PTI;
            AK = 1.0 - TFNF / (FKK + TFNF);
            ACK = BK * AK;
            sumr += (ACK + BK) * p1;
            //   SUMR = SUMR + (ACK+BK)*P1R;
            //   SUMI = SUMI + (ACK+BK)*P1I;
            BK = ACK;
            FKK = FKK - 1.0;
            //   M = N - I + 1;
            y[N - i] = p2;
            //   YR(M) = P2R;
            //   YI(M) = P2I;
        }
        //    60 CONTINUE;
        //    70 CONTINUE;
    }
    if IFNU > 0 {
        //GO TO 90;
        // DO 80 I=1,IFNU;
        for _i in 0..IFNU {
            let pt = p2;
            //   PTR = P2R;
            //   PTI = P2I;
            p2 = p1 + (FKK + FNF) * (rz * pt);

            //   P2R = P1R + (FKK+FNF)*(RZR*PTR-RZI*PTI);
            //   P2I = P1I + (FKK+FNF)*(RZR*PTI+RZI*PTR);
            p1 = pt;
            //   P1R = PTR;
            //   P1I = PTI;
            AK = 1.0 - TFNF / (FKK + TFNF);
            ACK = BK * AK;
            sumr += (ACK + BK) * p1;
            //   SUMR = SUMR + (ACK+BK)*P1R;
            //   SUMI = SUMI + (ACK+BK)*P1I;
            BK = ACK;
            FKK = FKK - 1.0;
        }
        //    80 CONTINUE;
        //    90 CONTINUE;
    }

    // PTR = ZR;
    // PTI = ZI;
    let pt = if KODE == Scaling::Scaled { c_zero() } else { z };
    // if (KODE == 2) PTR = ZEROR;
    // CALL ZLOG(RZR, RZI, STR, STI, IDUM);
    p1 = -FNF * rz.ln() + pt;
    // P1R = -FNF*STR + PTR;
    // P1I = -FNF*STI + PTI;
    let AP = gamma_ln(1.0 + FNF).unwrap();
    p1 -= AP;
    // let pt = p1 - AP;
    // PTR = P1R - AP;
    // PTI = P1I;
    //-----------------------------------------------------------------------;
    //     THE DIVISION CEXP(PT)/(SUM+P2) IS ALTERED TO AVOID OVERFLOW;
    //     IN THE DENOMINATOR BY SQUARING LARGE QUANTITIES;
    //-----------------------------------------------------------------------;
    p2 += sumr;
    // P2R = P2R + SUMR;
    // P2I = P2I + SUMI;
    // AP = ZABS(P2R,P2I);
    let AP = p2.abs();
    // P1R = 1.0/AP;
    // CALL ZEXP(PTR, PTI, STR, STI);
    ck = p1.exp() / AP;
    // CKR = STR*P1R;
    // CKI = STI*P1R;
    // PTR = P2R*P1R;
    // PTI = -P2I*P1R;
    // CALL ZMLT(CKR, CKI, PTR, PTI, CNORMR, CNORMI);
    let cnorm = ck * p2.conj() / AP;
    // DO 100 I=1,N;
    for i in 0..N {
        y[i] *= cnorm;
        //   STR = YR(I)*CNORMR - YI(I)*CNORMI;
        //   YI(I) = YR(I)*CNORMI + YI(I)*CNORMR;
        //   YR(I) = STR;
    }
    //   100 CONTINUE;
    // RETURN;
    return Ok((y, NZ));
    //   110 CONTINUE;
    //       NZ=-2; //falied to converge
    //       RETURN;
    //       END;
}
/*
fn ZWRSK(ZRR, ZRI, FNU, KODE, N, YR, YI, NZ, CWR, CWI,
     * TOL, ELIM, ALIM)
// ***BEGIN PROLOGUE  ZWRSK
// ***REFER TO  ZBESI,ZBESK
//
//     ZWRSK COMPUTES THE I BESSEL FUNCTION FOR RE(Z) >= 0.0 BY
//     NORMALIZING THE I FUNCTION RATIOS FROM ZRATI BY THE WRONSKIAN
//
// ***ROUTINES CALLED  d1mach,ZBKNU,ZRATI,ZABS
// ***END PROLOGUE  ZWRSK
//     COMPLEX CINU,CSCL,CT,CW,C1,C2,RCT,ST,Y,ZR
      EXTERNAL ZABS
      DOUBLE PRECISION ACT, ACW, ALIM, ASCLE, CINUI, CINUR, CSCLR, CTI,
     * CTR, CWI, CWR, C1I, C1R, C2I, C2R, ELIM, FNU, PTI, PTR, RACT,
     * STI, STR, TOL, YI, YR, ZRI, ZRR, ZABS, d1mach
      INTEGER I, KODE, N, NW, NZ
      DIMENSION YR(N), YI(N), CWR(2), CWI(2)
//-----------------------------------------------------------------------
//     I(FNU+I-1,Z) BY BACKWARD RECURRENCE FOR RATIOS
//     Y(I)=I(FNU+I,Z)/I(FNU+I-1,Z) FROM CRATI NORMALIZED BY THE
//     WRONSKIAN WITH K(FNU,Z) AND K(FNU+1,Z) FROM CBKNU.
//-----------------------------------------------------------------------
      NZ = 0
      CALL ZBKNU(ZRR, ZRI, FNU, KODE, 2, CWR, CWI, NW, TOL, ELIM, ALIM)
      if (NW != 0) GO TO 50
      CALL ZRATI(ZRR, ZRI, FNU, N, YR, YI, TOL)
//-----------------------------------------------------------------------
//     RECUR FORWARD ON I(FNU+1,Z) = R(FNU,Z)*I(FNU,Z),
//     R(FNU+J-1,Z)=Y(J),  J=1,...,N
//-----------------------------------------------------------------------
      CINUR = 1.0
      CINUI = 0.0
      if (KODE == 1) GO TO 10
      CINUR = DCOS(ZRI)
      CINUI = DSIN(ZRI)
   10 CONTINUE
//-----------------------------------------------------------------------
//     ON LOW EXPONENT MACHINES THE K FUNCTIONS CAN BE CLOSE TO BOTH
//     THE UNDER AND OVERFLOW LIMITS AND THE NORMALIZATION MUST BE
//     SCALED TO PREVENT OVER OR UNDERFLOW. CUOIK HAS DETERMINED THAT
//     THE RESULT IS ON SCALE.
//-----------------------------------------------------------------------
      ACW = ZABS(CWR(2),CWI(2))
      ASCLE = 1.0e+3*d1mach(1)/TOL
      CSCLR = 1.0
      if (ACW > ASCLE) GO TO 20
      CSCLR = 1.0/TOL
      GO TO 30
   20 CONTINUE
      ASCLE = 1.0/ASCLE
      if (ACW < ASCLE) GO TO 30
      CSCLR = TOL
   30 CONTINUE
      C1R = CWR(1)*CSCLR
      C1I = CWI(1)*CSCLR
      C2R = CWR(2)*CSCLR
      C2I = CWI(2)*CSCLR
      STR = YR(1)
      STI = YI(1)
//-----------------------------------------------------------------------
//     CINU=CINU*(CONJG(CT)/CABS(CT))*(1.0/CABS(CT) PREVENTS
//     UNDER- OR OVERFLOW PREMATURELY BY SQUARING CABS(CT)
//-----------------------------------------------------------------------
      PTR = STR*C1R - STI*C1I
      PTI = STR*C1I + STI*C1R
      PTR = PTR + C2R
      PTI = PTI + C2I
      CTR = ZRR*PTR - ZRI*PTI
      CTI = ZRR*PTI + ZRI*PTR
      ACT = ZABS(CTR,CTI)
      RACT = 1.0/ACT
      CTR = CTR*RACT
      CTI = -CTI*RACT
      PTR = CINUR*RACT
      PTI = CINUI*RACT
      CINUR = PTR*CTR - PTI*CTI
      CINUI = PTR*CTI + PTI*CTR
      YR(1) = CINUR*CSCLR
      YI(1) = CINUI*CSCLR
      if (N == 1) RETURN
      DO 40 I=2,N
        PTR = STR*CINUR - STI*CINUI
        CINUI = STR*CINUI + STI*CINUR
        CINUR = PTR
        STR = YR(I)
        STI = YI(I)
        YR(I) = CINUR*CSCLR
        YI(I) = CINUI*CSCLR
   40 CONTINUE
      RETURN
   50 CONTINUE
      NZ = -1
      if(NW == (-2)) NZ=-2
      RETURN
      END
*/

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
) -> Result<(Vec<Complex64>, usize), BesselError> {
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
    let mut cy = vec![c_zero(); N];
    if !(AZ <= 2.0) {
        if !(AZ * AZ * 0.25 > DFNU + 1.0) {
            //-----------------------------------------------------------------------
            //     POWER SERIES
            //-----------------------------------------------------------------------
            let NW;
            (cy, NW) = z_power_series(z, order, KODE, NN, machine_consts)?;
            let INW: usize = NW.abs().try_into().unwrap();
            NZ = NZ + INW;
            NN = NN - INW;
            if NN == 0 || NW >= 0 {
                return Ok((cy, NZ.try_into().unwrap()));
            }

            DFNU = order + ((NN as f64) - 1.0);
        }
    }

    if !(AZ <  machine_consts.rl)
          && (DFNU <= 1.0) //GO TO 30
          && !(AZ+AZ < DFNU*DFNU)
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
    let mut skip_az_rl_check = false;
    if !(DFNU <= 1.0) {
        //GO TO 70
        skip_az_rl_check = true;
        //  50 CONTINUE

        //-----------------------------------------------------------------------
        //     OVERFLOW AND UNDERFLOW TEST ON I SEQUENCE FOR MILLER ALGORITHM
        //-----------------------------------------------------------------------
        let (cy, nw) = zuoik(z, order, KODE, IKType::I, NN, cy, machine_consts)?;

        NZ = NZ + nw;
        NN = NN - nw;
        if NN == 0 {
            return Ok((cy, NZ));
        }
        DFNU = order + ((NN - 1) as f64);
        //     if (DFNU > FNUL) GO TO 110
        //     if (AZ > FNUL) GO TO 110
        if !((DFNU > machine_consts.fnul) || (AZ > machine_consts.fnul))
        //GO TO 110
        {
            //  60 CONTINUE
            // 'l60: loop{
            if !skip_az_rl_check && !(AZ > machine_consts.rl) {
                // GO TO 80
                //  70 CONTINUE
                //-----------------------------------------------------------------------
                //     MILLER ALGORITHM NORMALIZED BY THE SERIES
                //-----------------------------------------------------------------------
                let (cy, _) = ZMLRI(
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
            return Err(NotYetImplemented);
            /*
               80 CONTINUE
            //-----------------------------------------------------------------------
            //     MILLER ALGORITHM NORMALIZED BY THE WRONSKIAN
            //-----------------------------------------------------------------------
            //-----------------------------------------------------------------------
            //     OVERFLOW TEST ON K FUNCTIONS USED IN WRONSKIAN
            //-----------------------------------------------------------------------
            //       CALL ZUOIK(ZR, ZI, FNU, KODE, 2, 2, CWR, CWI, NW, TOL, ELIM,
            //      * ALIM)
            let (cy, NW) = zuoik(z, order, kode, IKType::K, 2, vec![c_one();2], machine_consts)
                  if (NW >= 0) GO TO 100
                  NZ = NN
                  DO 90 I=1,NN
                    CYR(I) = ZEROR
                    CYI(I) = ZEROI
               90 CONTINUE
                  RETURN
              100 CONTINUE
                  if (NW > 0) GO TO 130
                  CALL ZWRSK(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, CWR, CWI, TOL,
                 * ELIM, ALIM)
                  if (NW < 0) GO TO 130
                  return Ok(cy, NZ);
                  */
        }
        // 110 CONTINUE
        /*
        //-----------------------------------------------------------------------
        //     INCREMENT FNU+NN-1 UP TO FNUL, COMPUTE AND RECUR BACKWARD
        //-----------------------------------------------------------------------
              NUI = INT(SNGL(FNUL-DFNU)) + 1
              NUI = MAX0(NUI,0)
              CALL ZBUNI(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, NUI, NLAST, FNUL,
             * TOL, ELIM, ALIM)
              if (NW < 0) GO TO 130
              NZ = NZ + NW
              if (NLAST == 0) {return Ok(cy, NZ);}
              NN = NLAST
              GO TO 60
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
    return Err(NotYetImplemented);
}

/*
fn ZACAI(ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, RL, TOL,
     * ELIM, ALIM)
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
      EXTERNAL ZABS
      DOUBLE PRECISION ALIM, ARG, ASCLE, AZ, CSGNR, CSGNI, CSPNR,
     * CSPNI, C1R, C1I, C2R, C2I, CYR, CYI, DFNU, ELIM, FMR, FNU, PI,
     * RL, SGN, TOL, YY, YR, YI, ZR, ZI, ZNR, ZNI, d1mach, ZABS
      INTEGER INU, IUF, KODE, MR, N, NN, NW, NZ
      DIMENSION YR(N), YI(N), CYR(2), CYI(2)
      DATA PI / 3.14159265358979324 /
      NZ = 0
      ZNR = -ZR
      ZNI = -ZI
      AZ = ZABS(ZR,ZI)
      NN = N
      DFNU = FNU + ((N-1) as f64)
      if (AZ <= 2.0) GO TO 10
      if (AZ*AZ*0.25 > DFNU+1.0) GO TO 20
   10 CONTINUE
//-----------------------------------------------------------------------
//     POWER SERIES FOR THE I FUNCTION
//-----------------------------------------------------------------------
      CALL z_power_series(ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, TOL, ELIM, ALIM)
      GO TO 40
   20 CONTINUE
      if (AZ < RL) GO TO 30
//-----------------------------------------------------------------------
//     ASYMPTOTIC EXPANSION FOR LARGE Z FOR THE I FUNCTION
//-----------------------------------------------------------------------
      CALL ZASYI(ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, RL, TOL, ELIM,
     * ALIM)
      if (NW < 0) GO TO 80
      GO TO 40
   30 CONTINUE
//-----------------------------------------------------------------------
//     MILLER ALGORITHM NORMALIZED BY THE SERIES FOR THE I FUNCTION
//-----------------------------------------------------------------------
      CALL ZMLRI(ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, TOL)
      if(NW < 0) GO TO 80
   40 CONTINUE
//-----------------------------------------------------------------------
//     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
//-----------------------------------------------------------------------
      CALL ZBKNU(ZNR, ZNI, FNU, KODE, 1, CYR, CYI, NW, TOL, ELIM, ALIM)
      if (NW != 0) GO TO 80
      FMR = (MR as f64)
      SGN = -DSIGN(PI,FMR)
      CSGNR = 0.0
      CSGNI = SGN
      if (KODE == 1) GO TO 50
      YY = -ZNI
      CSGNR = -CSGNI*DSIN(YY)
      CSGNI = CSGNI*DCOS(YY)
   50 CONTINUE
//-----------------------------------------------------------------------
//     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
//     WHEN FNU IS LARGE
//-----------------------------------------------------------------------
      INU = INT(SNGL(FNU))
      ARG = (FNU-(INU as f64))*SGN
      CSPNR = DCOS(ARG)
      CSPNI = DSIN(ARG)
      if (MOD(INU,2) == 0) GO TO 60
      CSPNR = -CSPNR
      CSPNI = -CSPNI
   60 CONTINUE
      C1R = CYR(1)
      C1I = CYI(1)
      C2R = YR(1)
      C2I = YI(1)
      if (KODE == 1) GO TO 70
      IUF = 0
      ASCLE = 1.0e+3*d1mach(1)/TOL
      CALL ZS1S2(ZNR, ZNI, C1R, C1I, C2R, C2I, NW, ASCLE, ALIM, IUF)
      NZ = NZ + NW
   70 CONTINUE
      YR(1) = CSPNR*C1R - CSPNI*C1I + CSGNR*C2R - CSGNI*C2I
      YI(1) = CSPNR*C1I + CSPNI*C1R + CSGNR*C2I + CSGNI*C2R
      RETURN
   80 CONTINUE
      NZ = -1
      if(NW == (-2)) NZ=-2
      RETURN
      END
      */

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
fn ZBUNI(ZR, ZI, FNU, KODE, N, YR, YI, NZ, NUI, NLAST,
     * FNUL, TOL, ELIM, ALIM)
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
      EXTERNAL ZABS
      DOUBLE PRECISION ALIM, AX, AY, CSCLR, CSCRR, CYI, CYR, DFNU,
     * ELIM, FNU, FNUI, FNUL, GNU, RAZ, RZI, RZR, STI, STR, S1I, S1R,
     * S2I, S2R, TOL, YI, YR, ZI, ZR, ZABS, ASCLE, BRY, C1R, C1I, C1M,
     * d1mach
      INTEGER I, IFLAG, IFORM, K, KODE, N, NL, NLAST, NUI, NW, NZ
      DIMENSION YR(N), YI(N), CYR(2), CYI(2), BRY(3)
      NZ = 0
      AX = (ZR).abs()*1.7321
      AY = (ZI).abs()
      IFORM = 1
      if (AY > AX) IFORM = 2
      if (NUI == 0) GO TO 60
      FNUI = (NUI as f64)
      DFNU = FNU + ((N-1) as f64)
      GNU = DFNU + FNUI
      if (IFORM == 2) GO TO 10
//-----------------------------------------------------------------------
//     ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN
//     -PI/3 <= ARG(Z) <= PI/3
//-----------------------------------------------------------------------
      CALL ZUNI1(ZR, ZI, GNU, KODE, 2, CYR, CYI, NW, NLAST, FNUL, TOL,
     * ELIM, ALIM)
      GO TO 20
   10 CONTINUE
//-----------------------------------------------------------------------
//     ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*FRAC_PI_2)) FOR LARGE FNU
//     APPLIED IN PI/3 < ABS(ARG(Z)) <= PI/2 WHERE M=+I OR -I
//     AND FRAC_PI_2=PI/2
//-----------------------------------------------------------------------
      CALL ZUNI2(ZR, ZI, GNU, KODE, 2, CYR, CYI, NW, NLAST, FNUL, TOL,
     * ELIM, ALIM)
   20 CONTINUE
      if (NW < 0) GO TO 50
      if (NW != 0) GO TO 90
      STR = ZABS(CYR(1),CYI(1))
//----------------------------------------------------------------------
//     SCALE BACKWARD RECURRENCE, BRY(3) IS DEFINED BUT NEVER USED
//----------------------------------------------------------------------
      BRY(1)=1.0e+3*d1mach(1)/TOL
      BRY(2) = 1.0/BRY(1)
      BRY(3) = BRY(2)
      IFLAG = 2
      ASCLE = BRY(2)
      CSCLR = 1.0
      if (STR > BRY(1)) GO TO 21
      IFLAG = 1
      ASCLE = BRY(1)
      CSCLR = 1.0/TOL
      GO TO 25
   21 CONTINUE
      if (STR < BRY(2)) GO TO 25
      IFLAG = 3
      ASCLE=BRY(3)
      CSCLR = TOL
   25 CONTINUE
      CSCRR = 1.0/CSCLR
      S1R = CYR(2)*CSCLR
      S1I = CYI(2)*CSCLR
      S2R = CYR(1)*CSCLR
      S2I = CYI(1)*CSCLR
      RAZ = 1.0/ZABS(ZR,ZI)
      STR = ZR*RAZ
      STI = -ZI*RAZ
      RZR = (STR+STR)*RAZ
      RZI = (STI+STI)*RAZ
      DO 30 I=1,NUI
        STR = S2R
        STI = S2I
        S2R = (DFNU+FNUI)*(RZR*STR-RZI*STI) + S1R
        S2I = (DFNU+FNUI)*(RZR*STI+RZI*STR) + S1I
        S1R = STR
        S1I = STI
        FNUI = FNUI - 1.0
        if (IFLAG >= 3) GO TO 30
        STR = S2R*CSCRR
        STI = S2I*CSCRR
        C1R = (STR).abs()
        C1I = (STI).abs()
        C1M = DMAX1(C1R,C1I)
        if (C1M <= ASCLE) GO TO 30
        IFLAG = IFLAG+1
        ASCLE = BRY(IFLAG)
        S1R = S1R*CSCRR
        S1I = S1I*CSCRR
        S2R = STR
        S2I = STI
        CSCLR = CSCLR*TOL
        CSCRR = 1.0/CSCLR
        S1R = S1R*CSCLR
        S1I = S1I*CSCLR
        S2R = S2R*CSCLR
        S2I = S2I*CSCLR
   30 CONTINUE
      YR(N) = S2R*CSCRR
      YI(N) = S2I*CSCRR
      if (N == 1) RETURN
      NL = N - 1
      FNUI = (NL as f64)
      K = NL
      DO 40 I=1,NL
        STR = S2R
        STI = S2I
        S2R = (FNU+FNUI)*(RZR*STR-RZI*STI) + S1R
        S2I = (FNU+FNUI)*(RZR*STI+RZI*STR) + S1I
        S1R = STR
        S1I = STI
        STR = S2R*CSCRR
        STI = S2I*CSCRR
        YR(K) = STR
        YI(K) = STI
        FNUI = FNUI - 1.0
        K = K - 1
        if (IFLAG >= 3) GO TO 40
        C1R = (STR).abs()
        C1I = (STI).abs()
        C1M = DMAX1(C1R,C1I)
        if (C1M <= ASCLE) GO TO 40
        IFLAG = IFLAG+1
        ASCLE = BRY(IFLAG)
        S1R = S1R*CSCRR
        S1I = S1I*CSCRR
        S2R = STR
        S2I = STI
        CSCLR = CSCLR*TOL
        CSCRR = 1.0/CSCLR
        S1R = S1R*CSCLR
        S1I = S1I*CSCLR
        S2R = S2R*CSCLR
        S2I = S2I*CSCLR
   40 CONTINUE
      RETURN
   50 CONTINUE
      NZ = -1
      if(NW == (-2)) NZ=-2
      RETURN
   60 CONTINUE
      if (IFORM == 2) GO TO 70
//-----------------------------------------------------------------------
//     ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN
//     -PI/3 <= ARG(Z) <= PI/3
//-----------------------------------------------------------------------
      CALL ZUNI1(ZR, ZI, FNU, KODE, N, YR, YI, NW, NLAST, FNUL, TOL,
     * ELIM, ALIM)
      GO TO 80
   70 CONTINUE
//-----------------------------------------------------------------------
//     ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*FRAC_PI_2)) FOR LARGE FNU
//     APPLIED IN PI/3 < ABS(ARG(Z)) <= PI/2 WHERE M=+I OR -I
//     AND FRAC_PI_2=PI/2
//-----------------------------------------------------------------------
      CALL ZUNI2(ZR, ZI, FNU, KODE, N, YR, YI, NW, NLAST, FNUL, TOL,
     * ELIM, ALIM)
   80 CONTINUE
      if (NW < 0) GO TO 50
      NZ = NW
      RETURN
   90 CONTINUE
      NLAST = N
      RETURN
      END
fn ZUNI1(ZR, ZI, FNU, KODE, N, YR, YI, NZ, NLAST, FNUL,
     * TOL, ELIM, ALIM)
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
      EXTERNAL ZABS
      DOUBLE PRECISION ALIM, APHI, ASCLE, BRY, CONER, CRSC,
     * CSCL, CSRR, CSSR, CWRKI, CWRKR, C1R, C2I, C2M, C2R, ELIM, FN,
     * FNU, FNUL, PHII, PHIR, RAST, RS1, RZI, RZR, STI, STR, SUMI,
     * SUMR, S1I, S1R, S2I, S2R, TOL, YI, YR, ZEROI, ZEROR, ZETA1I,
     * ZETA1R, ZETA2I, ZETA2R, ZI, ZR, CYR, CYI, d1mach, ZABS
      INTEGER I, IFLAG, INIT, K, KODE, M, N, ND, NLAST, NN, NUF, NW, NZ
      DIMENSION BRY(3), YR(N), YI(N), CWRKR(16), CWRKI(16), CSSR(3),
     * CSRR(3), CYR(2), CYI(2)
      DATA ZEROR,ZEROI,CONER / 0.0, 0.0, 1.0 /
//
      NZ = 0
      ND = N
      NLAST = 0
//-----------------------------------------------------------------------
//     COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG-
//     NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE,
//     EXP(ALIM)=EXP(ELIM)*TOL
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
//-----------------------------------------------------------------------
//     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
//-----------------------------------------------------------------------
      FN = DMAX1(FNU,1.0)
      INIT = 0
      CALL ZUNIK(ZR, ZI, FN, 1, 1, TOL, INIT, PHIR, PHII, ZETA1R,
     * ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI)
      if (KODE == 1) GO TO 10
      STR = ZR + ZETA2R
      STI = ZI + ZETA2I
      RAST = FN/ZABS(STR,STI)
      STR = STR*RAST*RAST
      STI = -STI*RAST*RAST
      S1R = -ZETA1R + STR
      S1I = -ZETA1I + STI
      GO TO 20
   10 CONTINUE
      S1R = -ZETA1R + ZETA2R
      S1I = -ZETA1I + ZETA2I
   20 CONTINUE
      RS1 = S1R
      if ((RS1).abs() > ELIM) GO TO 130
   30 CONTINUE
      NN = MIN0(2,ND)
      DO 80 I=1,NN
        FN = FNU + ((ND-I) as f64)
        INIT = 0
        CALL ZUNIK(ZR, ZI, FN, 1, 0, TOL, INIT, PHIR, PHII, ZETA1R,
     *   ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI)
        if (KODE == 1) GO TO 40
        STR = ZR + ZETA2R
        STI = ZI + ZETA2I
        RAST = FN/ZABS(STR,STI)
        STR = STR*RAST*RAST
        STI = -STI*RAST*RAST
        S1R = -ZETA1R + STR
        S1I = -ZETA1I + STI + ZI
        GO TO 50
   40   CONTINUE
        S1R = -ZETA1R + ZETA2R
        S1I = -ZETA1I + ZETA2I
   50   CONTINUE
//-----------------------------------------------------------------------
//     TEST FOR UNDERFLOW AND OVERFLOW
//-----------------------------------------------------------------------
        RS1 = S1R
        if ((RS1).abs() > ELIM) GO TO 110
        if (I == 1) IFLAG = 2
        if ((RS1).abs() < ALIM) GO TO 60
//-----------------------------------------------------------------------
//     REFINE  TEST AND SCALE
//-----------------------------------------------------------------------
        APHI = ZABS(PHIR,PHII)
        RS1 = RS1 + DLOG(APHI)
        if ((RS1).abs() > ELIM) GO TO 110
        if (I == 1) IFLAG = 1
        if (RS1 < 0.0) GO TO 60
        if (I == 1) IFLAG = 3
   60   CONTINUE
//-----------------------------------------------------------------------
//     SCALE S1 if CABS(S1) < ASCLE
//-----------------------------------------------------------------------
        S2R = PHIR*SUMR - PHII*SUMI
        S2I = PHIR*SUMI + PHII*SUMR
        STR = DEXP(S1R)*CSSR(IFLAG)
        S1R = STR*DCOS(S1I)
        S1I = STR*DSIN(S1I)
        STR = S2R*S1R - S2I*S1I
        S2I = S2R*S1I + S2I*S1R
        S2R = STR
        if (IFLAG != 1) GO TO 70
        CALL ZUunderflowCHK(S2R, S2I, NW, BRY(1), TOL)
        if (NW != 0) GO TO 110
   70   CONTINUE
        CYR(I) = S2R
        CYI(I) = S2I
        M = ND - I + 1
        YR(M) = S2R*CSRR(IFLAG)
        YI(M) = S2I*CSRR(IFLAG)
   80 CONTINUE
      if (ND <= 2) GO TO 100
      RAST = 1.0/ZABS(ZR,ZI)
      STR = ZR*RAST
      STI = -ZI*RAST
      RZR = (STR+STR)*RAST
      RZI = (STI+STI)*RAST
      BRY(2) = 1.0/BRY(1)
      BRY(3) = d1mach(2)
      S1R = CYR(1)
      S1I = CYI(1)
      S2R = CYR(2)
      S2I = CYI(2)
      C1R = CSRR(IFLAG)
      ASCLE = BRY(IFLAG)
      K = ND - 2
      FN = (K as f64)
      DO 90 I=3,ND
        C2R = S2R
        C2I = S2I
        S2R = S1R + (FNU+FN)*(RZR*C2R-RZI*C2I)
        S2I = S1I + (FNU+FN)*(RZR*C2I+RZI*C2R)
        S1R = C2R
        S1I = C2I
        C2R = S2R*C1R
        C2I = S2I*C1R
        YR(K) = C2R
        YI(K) = C2I
        K = K - 1
        FN = FN - 1.0
        if (IFLAG >= 3) GO TO 90
        STR = (C2R).abs()
        STI = (C2I).abs()
        C2M = DMAX1(STR,STI)
        if (C2M <= ASCLE) GO TO 90
        IFLAG = IFLAG + 1
        ASCLE = BRY(IFLAG)
        S1R = S1R*C1R
        S1I = S1I*C1R
        S2R = C2R
        S2I = C2I
        S1R = S1R*CSSR(IFLAG)
        S1I = S1I*CSSR(IFLAG)
        S2R = S2R*CSSR(IFLAG)
        S2I = S2I*CSSR(IFLAG)
        C1R = CSRR(IFLAG)
   90 CONTINUE
  100 CONTINUE
      RETURN
//-----------------------------------------------------------------------
//     SET UNDERFLOW AND UPDATE PARAMETERS
//-----------------------------------------------------------------------
  110 CONTINUE
      if (RS1 > 0.0) GO TO 120
      YR(ND) = ZEROR
      YI(ND) = ZEROI
      NZ = NZ + 1
      ND = ND - 1
      if (ND == 0) GO TO 100
      CALL ZUOIK(ZR, ZI, FNU, KODE, 1, ND, YR, YI, NUF, TOL, ELIM, ALIM)
      if (NUF < 0) GO TO 120
      ND = ND - NUF
      NZ = NZ + NUF
      if (ND == 0) GO TO 100
      FN = FNU + ((ND-1) as f64)
      if (FN >= FNUL) GO TO 30
      NLAST = ND
      RETURN
  120 CONTINUE
      NZ = -1
      RETURN
  130 CONTINUE
      if (RS1 > 0.0) GO TO 120
      NZ = N
      DO 140 I=1,N
        YR(I) = ZEROR
        YI(I) = ZEROI
  140 CONTINUE
      RETURN
      END
fn ZUNI2(ZR, ZI, FNU, KODE, N, YR, YI, NZ, NLAST, FNUL,
     * TOL, ELIM, ALIM)
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
      EXTERNAL ZABS
      DOUBLE PRECISION AARG, AIC, AII, AIR, ALIM, ANG, APHI, ARGI,
     * ARGR, ASCLE, ASUMI, ASUMR, BRY, BSUMI, BSUMR, CIDI, CIPI, CIPR,
     * CONER, CRSC, CSCL, CSRR, CSSR, C1R, C2I, C2M, C2R, DAII,
     * DAIR, ELIM, FN, FNU, FNUL, FRAC_PI_2, PHII, PHIR, RAST, RAZ, RS1, RZI,
     * RZR, STI, STR, S1I, S1R, S2I, S2R, TOL, YI, YR, ZBI, ZBR, ZEROI,
     * ZEROR, ZETA1I, ZETA1R, ZETA2I, ZETA2R, ZI, ZNI, ZNR, ZR, CYR,
     * CYI, d1mach, ZABS, CAR, SAR
      INTEGER I, IFLAG, IN, INU, J, K, KODE, N, NAI, ND, NDAI, NLAST,
     * NN, NUF, NW, NZ, IDUM
      DIMENSION BRY(3), YR(N), YI(N), CIPR(4), CIPI(4), CSSR(3),
     * CSRR(3), CYR(2), CYI(2)
      DATA ZEROR,ZEROI,CONER / 0.0, 0.0, 1.0 /
      DATA CIPR(1),CIPI(1),CIPR(2),CIPI(2),CIPR(3),CIPI(3),CIPR(4),
     * CIPI(4)/ 1.0,0.0, 0.0,1.0, -1.0,0.0, 0.0,-1.0/
      DATA FRAC_PI_2, AIC  /
     1      1.57079632679489662e+00,     1.265512123484645396e+00/
//
      NZ = 0
      ND = N
      NLAST = 0
//-----------------------------------------------------------------------
//     COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG-
//     NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE,
//     EXP(ALIM)=EXP(ELIM)*TOL
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
//-----------------------------------------------------------------------
//     ZN IS IN THE RIGHT HALF PLANE AFTER ROTATION BY CI OR -CI
//-----------------------------------------------------------------------
      ZNR = ZI
      ZNI = -ZR
      ZBR = ZR
      ZBI = ZI
      CIDI = -CONER
      INU = INT(SNGL(FNU))
      ANG = FRAC_PI_2*(FNU-(INU as f64))
      C2R = DCOS(ANG)
      C2I = DSIN(ANG)
      CAR = C2R
      SAR = C2I
      IN = INU + N - 1
      IN = MOD(IN,4) + 1
      STR = C2R*CIPR(IN) - C2I*CIPI(IN)
      C2I = C2R*CIPI(IN) + C2I*CIPR(IN)
      C2R = STR
      if (ZI > 0.0) GO TO 10
      ZNR = -ZNR
      ZBI = -ZBI
      CIDI = -CIDI
      C2I = -C2I
   10 CONTINUE
//-----------------------------------------------------------------------
//     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
//-----------------------------------------------------------------------
      FN = DMAX1(FNU,1.0)
      CALL ZUNHJ(ZNR, ZNI, FN, 1, TOL, PHIR, PHII, ARGR, ARGI, ZETA1R,
     * ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI)
      if (KODE == 1) GO TO 20
      STR = ZBR + ZETA2R
      STI = ZBI + ZETA2I
      RAST = FN/ZABS(STR,STI)
      STR = STR*RAST*RAST
      STI = -STI*RAST*RAST
      S1R = -ZETA1R + STR
      S1I = -ZETA1I + STI
      GO TO 30
   20 CONTINUE
      S1R = -ZETA1R + ZETA2R
      S1I = -ZETA1I + ZETA2I
   30 CONTINUE
      RS1 = S1R
      if ((RS1).abs() > ELIM) GO TO 150
   40 CONTINUE
      NN = MIN0(2,ND)
      DO 90 I=1,NN
        FN = FNU + ((ND-I) as f64)
        CALL ZUNHJ(ZNR, ZNI, FN, 0, TOL, PHIR, PHII, ARGR, ARGI,
     *   ZETA1R, ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI)
        if (KODE == 1) GO TO 50
        STR = ZBR + ZETA2R
        STI = ZBI + ZETA2I
        RAST = FN/ZABS(STR,STI)
        STR = STR*RAST*RAST
        STI = -STI*RAST*RAST
        S1R = -ZETA1R + STR
        S1I = -ZETA1I + STI + (ZI).abs()
        GO TO 60
   50   CONTINUE
        S1R = -ZETA1R + ZETA2R
        S1I = -ZETA1I + ZETA2I
   60   CONTINUE
//-----------------------------------------------------------------------
//     TEST FOR UNDERFLOW AND OVERFLOW
//-----------------------------------------------------------------------
        RS1 = S1R
        if ((RS1).abs() > ELIM) GO TO 120
        if (I == 1) IFLAG = 2
        if ((RS1).abs() < ALIM) GO TO 70
//-----------------------------------------------------------------------
//     REFINE  TEST AND SCALE
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
        APHI = ZABS(PHIR,PHII)
        AARG = ZABS(ARGR,ARGI)
        RS1 = RS1 + DLOG(APHI) - 0.25*DLOG(AARG) - AIC
        if ((RS1).abs() > ELIM) GO TO 120
        if (I == 1) IFLAG = 1
        if (RS1 < 0.0) GO TO 70
        if (I == 1) IFLAG = 3
   70   CONTINUE
//-----------------------------------------------------------------------
//     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
//     EXPONENT EXTREMES
//-----------------------------------------------------------------------
        CALL ZAIRY(ARGR, ARGI, 0, 2, AIR, AII, NAI, IDUM)
        CALL ZAIRY(ARGR, ARGI, 1, 2, DAIR, DAII, NDAI, IDUM)
        STR = DAIR*BSUMR - DAII*BSUMI
        STI = DAIR*BSUMI + DAII*BSUMR
        STR = STR + (AIR*ASUMR-AII*ASUMI)
        STI = STI + (AIR*ASUMI+AII*ASUMR)
        S2R = PHIR*STR - PHII*STI
        S2I = PHIR*STI + PHII*STR
        STR = DEXP(S1R)*CSSR(IFLAG)
        S1R = STR*DCOS(S1I)
        S1I = STR*DSIN(S1I)
        STR = S2R*S1R - S2I*S1I
        S2I = S2R*S1I + S2I*S1R
        S2R = STR
        if (IFLAG != 1) GO TO 80
        CALL ZUunderflowCHK(S2R, S2I, NW, BRY(1), TOL)
        if (NW != 0) GO TO 120
   80   CONTINUE
        if (ZI <= 0.0) S2I = -S2I
        STR = S2R*C2R - S2I*C2I
        S2I = S2R*C2I + S2I*C2R
        S2R = STR
        CYR(I) = S2R
        CYI(I) = S2I
        J = ND - I + 1
        YR(J) = S2R*CSRR(IFLAG)
        YI(J) = S2I*CSRR(IFLAG)
        STR = -C2I*CIDI
        C2I = C2R*CIDI
        C2R = STR
   90 CONTINUE
      if (ND <= 2) GO TO 110
      RAZ = 1.0/ZABS(ZR,ZI)
      STR = ZR*RAZ
      STI = -ZI*RAZ
      RZR = (STR+STR)*RAZ
      RZI = (STI+STI)*RAZ
      BRY(2) = 1.0/BRY(1)
      BRY(3) = d1mach(2)
      S1R = CYR(1)
      S1I = CYI(1)
      S2R = CYR(2)
      S2I = CYI(2)
      C1R = CSRR(IFLAG)
      ASCLE = BRY(IFLAG)
      K = ND - 2
      FN = (K as f64)
      DO 100 I=3,ND
        C2R = S2R
        C2I = S2I
        S2R = S1R + (FNU+FN)*(RZR*C2R-RZI*C2I)
        S2I = S1I + (FNU+FN)*(RZR*C2I+RZI*C2R)
        S1R = C2R
        S1I = C2I
        C2R = S2R*C1R
        C2I = S2I*C1R
        YR(K) = C2R
        YI(K) = C2I
        K = K - 1
        FN = FN - 1.0
        if (IFLAG >= 3) GO TO 100
        STR = (C2R).abs()
        STI = (C2I).abs()
        C2M = DMAX1(STR,STI)
        if (C2M <= ASCLE) GO TO 100
        IFLAG = IFLAG + 1
        ASCLE = BRY(IFLAG)
        S1R = S1R*C1R
        S1I = S1I*C1R
        S2R = C2R
        S2I = C2I
        S1R = S1R*CSSR(IFLAG)
        S1I = S1I*CSSR(IFLAG)
        S2R = S2R*CSSR(IFLAG)
        S2I = S2I*CSSR(IFLAG)
        C1R = CSRR(IFLAG)
  100 CONTINUE
  110 CONTINUE
      RETURN
  120 CONTINUE
      if (RS1 > 0.0) GO TO 140
//-----------------------------------------------------------------------
//     SET UNDERFLOW AND UPDATE PARAMETERS
//-----------------------------------------------------------------------
      YR(ND) = ZEROR
      YI(ND) = ZEROI
      NZ = NZ + 1
      ND = ND - 1
      if (ND == 0) GO TO 110
      CALL ZUOIK(ZR, ZI, FNU, KODE, 1, ND, YR, YI, NUF, TOL, ELIM, ALIM)
      if (NUF < 0) GO TO 140
      ND = ND - NUF
      NZ = NZ + NUF
      if (ND == 0) GO TO 110
      FN = FNU + ((ND-1) as f64)
      if (FN < FNUL) GO TO 130
//      FN = CIDI
//      J = NUF + 1
//      K = MOD(J,4) + 1
//      S1R = CIPR(K)
//      S1I = CIPI(K)
//      if (FN < 0.0) S1I = -S1I
//      STR = C2R*S1R - C2I*S1I
//      C2I = C2R*S1I + C2I*S1R
//      C2R = STR
      IN = INU + ND - 1
      IN = MOD(IN,4) + 1
      C2R = CAR*CIPR(IN) - SAR*CIPI(IN)
      C2I = CAR*CIPI(IN) + SAR*CIPR(IN)
      if (ZI <= 0.0) C2I = -C2I
      GO TO 40
  130 CONTINUE
      NLAST = ND
      RETURN
  140 CONTINUE
      NZ = -1
      RETURN
  150 CONTINUE
      if (RS1 > 0.0) GO TO 140
      NZ = N
      DO 160 I=1,N
        YR(I) = ZEROR
        YI(I) = ZEROI
  160 CONTINUE
      RETURN
      END
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
