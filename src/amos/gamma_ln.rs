use thiserror::Error;

use super::machine::{d1mach, i1mach};

#[derive(Error, Debug, PartialEq, Eq)]
pub enum GammaError {
    #[error("Gamma ln can only be calculated for input z > 0.0")]
    ZLessThanZero,
}

// DOUBLE PRECISION FUNCTION DGAMLN(Z,IERR)
pub fn gamma_ln(z: f64) -> Result<f64, GammaError> {
    // ***BEGIN PROLOGUE  DGAMLN
    // ***DATE WRITTEN   830501   (YYMMDD)
    // ***REVISION DATE  830501   (YYMMDD)
    // ***CATEGORY NO.  B5F
    // ***KEYWORDS  GAMMA FUNCTION,LOGARITHM OF GAMMA FUNCTION
    // ***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
    // ***PURPOSE  TO COMPUTE THE LOGARITHM OF THE GAMMA FUNCTION
    // ***DESCRIPTION
    //
    //               **** A DOUBLE PRECISION ROUTINE ****
    //         DGAMLN COMPUTES THE NATURAL LOG OF THE GAMMA FUNCTION FOR
    //         Z > 0.  THE ASYMPTOTIC EXPANSION IS USED TO GENERATE VALUES
    //         GREATER THAN ZMIN WHICH ARE ADJUSTED BY THE RECURSION
    //         G(Z+1)=Z*G(Z) FOR Z <= ZMIN.  THE FUNCTION WAS MADE AS
    //         PORTABLE AS POSSIBLE BY COMPUTIMG ZMIN FROM THE NUMBER OF BASE
    //         10 DIGITS IN A WORD, RLN=AMAX1(-ALOG10(R1MACH(4)),0.5E-18)
    //         LIMITED TO 18 DIGITS OF (RELATIVE) ACCURACY.
    //
    //         SINCE INTEGER ARGUMENTS ARE COMMON, A TABLE LOOK UP ON 100
    //         VALUES IS USED FOR SPEED OF EXECUTION.
    //
    //     DESCRIPTION OF ARGUMENTS
    //
    //         INPUT      Z IS UBLE PRECISION
    //           Z      - ARGUMENT, Z > 0.0
    //
    //         OUTPUT      DGAMLN IS DOUBLE PRECISION
    //           DGAMLN  - NATURAL LOG OF THE GAMMA FUNCTION AT Z != 0.0
    //           IERR    - ERROR FLAG
    //                     IERR=0, NORMAL RETURN, COMPUTATION COMPLETED
    //                     IERR=1, Z <= 0.0,    NO COMPUTATION
    //
    //
    // ***REFERENCES  COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
    //                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
    // ***ROUTINES CALLED  i1mach,d1mach
    // ***END PROLOGUE  DGAMLN
    //       DOUBLE PRECISION CF, CON, FLN, FZ, GLN, RLN, S, TLG, TRM, TST,
    //      * T1, WDTOL, Z, ZDMY, ZINC, ZM, ZMIN, ZP, ZSQ, d1mach
    //       INTEGER I, IERR, I1M, K, MZ, NZ, i1mach
    //       DIMENSION CF(22), GLN(100)
    //           LNGAMMA(N), N=1,100
    //       DATA GLN(1), GLN(2), GLN(3), GLN(4), GLN(5), GLN(6), GLN(7),
    //      1     GLN(8), GLN(9), GLN(10), GLN(11), GLN(12), GLN(13), GLN(14),
    //      2     GLN(15), GLN(16), GLN(17), GLN(18), GLN(19), GLN(20),
    //      3     GLN(21), GLN(22)/
    const INT_GAMMA_LOG_DATA: [f64; 101] = [
        0.0, // value at GLN[0] - can never be reached.
        0.00000000000000000e+00,
        0.00000000000000000e+00,
        6.93147180559945309e-01,
        1.79175946922805500e+00,
        3.17805383034794562e+00,
        4.78749174278204599e+00,
        6.57925121201010100e+00,
        8.52516136106541430e+00,
        1.06046029027452502e+01,
        1.28018274800814696e+01,
        1.51044125730755153e+01,
        1.75023078458738858e+01,
        1.99872144956618861e+01,
        2.25521638531234229e+01,
        2.51912211827386815e+01,
        2.78992713838408916e+01,
        3.06718601060806728e+01,
        3.35050734501368889e+01,
        3.63954452080330536e+01,
        3.93398841871994940e+01,
        4.23356164607534850e+01,
        4.53801388984769080e+01,
        //       DATA GLN(23), GLN(24), GLN(25), GLN(26), GLN(27), GLN(28),
        //      1     GLN(29), GLN(30), GLN(31), GLN(32), GLN(33), GLN(34),
        //      2     GLN(35), GLN(36), GLN(37), GLN(38), GLN(39), GLN(40),
        //      3     GLN(41), GLN(42), GLN(43), GLN(44)/
        4.84711813518352239e+01,
        5.16066755677643736e+01,
        5.47847293981123192e+01,
        5.80036052229805199e+01,
        6.12617017610020020e+01,
        6.45575386270063311e+01,
        6.78897431371815350e+01,
        7.12570389671680090e+01,
        7.46582363488301644e+01,
        7.80922235533153106e+01,
        8.15579594561150372e+01,
        8.50544670175815174e+01,
        8.85808275421976788e+01,
        9.21361756036870925e+01,
        9.57196945421432025e+01,
        9.93306124547874269e+01,
        1.02968198614513813e+02,
        1.06631760260643459e+02,
        1.10320639714757395e+02,
        1.14034211781461703e+02,
        1.17771881399745072e+02,
        1.21533081515438634e+02,
        //       DATA GLN(45), GLN(46), GLN(47), GLN(48), GLN(49), GLN(50),
        //      1     GLN(51), GLN(52), GLN(53), GLN(54), GLN(55), GLN(56),
        //      2     GLN(57), GLN(58), GLN(59), GLN(60), GLN(61), GLN(62),
        //      3     GLN(63), GLN(64), GLN(65), GLN(66)/
        1.25317271149356895e+02,
        1.29123933639127215e+02,
        1.32952575035616310e+02,
        1.36802722637326368e+02,
        1.40673923648234259e+02,
        1.44565743946344886e+02,
        1.48477766951773032e+02,
        1.52409592584497358e+02,
        1.56360836303078785e+02,
        1.60331128216630907e+02,
        1.64320112263195181e+02,
        1.68327445448427652e+02,
        1.72352797139162802e+02,
        1.76395848406997352e+02,
        1.80456291417543771e+02,
        1.84533828861449491e+02,
        1.88628173423671591e+02,
        1.92739047287844902e+02,
        1.96866181672889994e+02,
        2.01009316399281527e+02,
        2.05168199482641199e+02,
        2.09342586752536836e+02,
        //       DATA GLN(67), GLN(68), GLN(69), GLN(70), GLN(71), GLN(72),
        //      1     GLN(73), GLN(74), GLN(75), GLN(76), GLN(77), GLN(78),
        //      2     GLN(79), GLN(80), GLN(81), GLN(82), GLN(83), GLN(84),
        //      3     GLN(85), GLN(86), GLN(87), GLN(88)/
        2.13532241494563261e+02,
        2.17736934113954227e+02,
        2.21956441819130334e+02,
        2.26190548323727593e+02,
        2.30439043565776952e+02,
        2.34701723442818268e+02,
        2.38978389561834323e+02,
        2.43268849002982714e+02,
        2.47572914096186884e+02,
        2.51890402209723194e+02,
        2.56221135550009525e+02,
        2.60564940971863209e+02,
        2.64921649798552801e+02,
        2.69291097651019823e+02,
        2.73673124285693704e+02,
        2.78067573440366143e+02,
        2.82474292687630396e+02,
        2.86893133295426994e+02,
        2.91323950094270308e+02,
        2.95766601350760624e+02,
        3.00220948647014132e+02,
        3.04686856765668715e+02,
        //       DATA GLN(89), GLN(90), GLN(91), GLN(92), GLN(93), GLN(94),
        //      1     GLN(95), GLN(96), GLN(97), GLN(98), GLN(99), GLN(100)/
        3.09164193580146922e+02,
        3.13652829949879062e+02,
        3.18152639620209327e+02,
        3.22663499126726177e+02,
        3.27185287703775217e+02,
        3.31717887196928473e+02,
        3.36261181979198477e+02,
        3.40815058870799018e+02,
        3.45379407062266854e+02,
        3.49954118040770237e+02,
        3.54539085519440809e+02,
        3.59134205369575399e+02,
    ];
    //             COEFFICIENTS OF ASYMPTOTIC EXPANSION
    //       DATA CF(1), CF(2), CF(3), CF(4), CF(5), CF(6), CF(7), CF(8),
    //      CF(9), CF(10), CF(11), CF(12), CF(13), CF(14), CF(15),
    //      CF(16), CF(17), CF(18), CF(19), CF(20), CF(21), CF(22)/
    const CF: [f64; 22] = [
        8.33333333333333333e-02,
        -2.77777777777777778e-03,
        7.93650793650793651e-04,
        -5.95238095238095238e-04,
        8.41750841750841751e-04,
        -1.91752691752691753e-03,
        6.41025641025641026e-03,
        -2.95506535947712418e-02,
        1.79644372368830573e-01,
        -1.39243221690590112e+00,
        1.34028640441683920e+01,
        -1.56848284626002017e+02,
        2.19310333333333333e+03,
        -3.61087712537249894e+04,
        6.91472268851313067e+05,
        -1.52382215394074162e+07,
        3.82900751391414141e+08,
        -1.08822660357843911e+10,
        3.47320283765002252e+11,
        -1.23696021422692745e+13,
        4.88788064793079335e+14,
        -2.13203339609193739e+16,
    ];
    //
    //             LN(2*PI)
    const CON: f64 = 1.83787706640934548e+00;
    //
    // ***FIRST EXECUTABLE STATEMENT  DGAMLN
    // IERR=0
    if z <= 0.0 {
        return Err(GammaError::ZLessThanZero);
    }
    let z_int = z as usize;

    if z <= 100.1 {
        //GO TO 10
        // let NZ = z as u64;
        // let FZ = z - (NZ as f64);

        if z.fract() == 0.0 {
            return Ok(INT_GAMMA_LOG_DATA[z_int]);
            // if (NZ > 100) GO TO 10
            // DGAMLN = GLN(NZ)
            // RETURN
        }
    }
    //    10 CONTINUE
    let wdtol = d1mach(4).max(0.5e-18);
    let i1m = i1mach(14);
    let rln = d1mach(5) * i1m as f64;
    let mut fln = rln.min(20.0);
    fln = fln.max(3.0);
    fln = fln - 3.0;
    let zm = 1.8000 + 0.3875 * fln;
    let mz = zm as u64 + 1;
    let zmin = mz as f64;
    let mut zdmy = z;
    let mut zinc = 0.0;
    if !(z >= zmin) {
        //GO TO 20;
        zinc = zmin - (z_int as f64);
        zdmy = z + zinc;
    }
    //    20 CONTINUE;
    let mut zp = 1.0 / zdmy;
    let t1 = CF[0] * zp;
    let mut s = t1;
    if !(zp < wdtol) {
        //GO TO 40;
        let zsq = zp * zp;
        let tst = t1 * wdtol;
        // DO 30 K=2,22;
        for k in 1..22 {
            zp = zp * zsq;
            let trm = CF[k] * zp;
            if trm.abs() < tst {
                break;
            } //GO TO 40;
            s = s + trm;
        } // 30 CONTINUE;
    } // 40 CONTINUE;
    if zinc == 0.0 {
        //GO TO 50;
        let tlg = z.ln(); //DLOG(Z);
        return Ok(z * (tlg - 1.0) + 0.5 * (CON - tlg) + s);
    }
    //    50 CONTINUE;
    let mut zp = 1.0;
    let nz = zinc as usize;
    // DO 60 I=1,NZ;
    for i in 0..nz {
        zp = zp * (z + (i as f64));
    } // 60 CONTINUE;
    let tlg = zdmy.ln(); //DLOG(ZDMY);
    return Ok(zdmy * (tlg - 1.0) - zp.ln() + 0.5 * (CON - tlg) + s);
    //     RETURN;
    //
    //
    //    70 CONTINUE
    //       IERR=1
    //       RETURN
    //       END
}
