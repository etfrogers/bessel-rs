use std::f64;

use thiserror::Error;

use crate::types::BesselFloat;

#[derive(Error, Debug, PartialEq, Eq)]
pub(crate) enum GammaError {
    #[error("Gamma ln can only be calculated for input z > 0.0")]
    ZLessThanZero,
}

/// Computes the natural log of the Gamma function for `z > 0`.
///
/// # Description
///
/// The asymptotic expansion is used to generate values greater than `z_min`,
/// which are adjusted by the recursion `G(z+1) = z * G(z)` for `z <= z_min`.
///
/// Since integer arguments are common, a table lookup on 100 values is used
/// for speed of execution.
///
/// # Arguments
///
/// * `z` - Argument, must be `> 0.0`.
///
/// # Returns
///
/// * `Ok(T)` - The natural log of the Gamma function at `z`.
/// * `Err(GammaError::ZLessThanZero)` - If `z <= 0.0`, computation fails.
///
/// # Original Amos Fortran Source
/// * **Author:** Donald E. Amos, Sandia National Laboratories
/// * **Date Written:** May 1, 1983
/// * **Reference:** Computation of Bessel Functions of Complex Argument by D.E. Amos, SAND83-0083, May, 1983.
pub(crate) fn gamma_ln<T: BesselFloat>(z: T) -> Result<T, GammaError> {
    if z <= T::zero() {
        return Err(GammaError::ZLessThanZero);
    }

    if z <= T::from_f64(100.1) && z.fract() == T::zero() {
        return Ok(T::from_f64(INT_GAMMA_LOG_DATA[z.to_usize().unwrap()]));
    }

    let working_tolerance = T::EPSILON.max(T::from_f64(0.5e-18));
    let digits_per_bit = (T::from_f64(T::RADIX as f64)).log10();
    let mantissa_base10_digits = T::from_f64(T::MANTISSA_DIGITS as f64) * digits_per_bit;
    let fln = mantissa_base10_digits.clamp(T::from_f64(3.0), T::from_f64(20.0)) - T::from_f64(3.0);
    let z_min = (T::from_f64(1.8000) + T::from_f64(0.3875) * fln).round() + T::one();

    let (z_increment, z_modified) = if z < z_min {
        let zinc = z_min - z.round();
        (zinc.to_usize().unwrap(), z + zinc)
    } else {
        (0, z)
    };

    let mut zp = T::one() / z_modified;
    let term1 = T::from_f64(COEFFS[0]) * zp;
    let mut s = term1;
    if zp >= working_tolerance {
        let zsq = zp * zp;
        let test = term1 * working_tolerance;
        for coeff in COEFFS.iter().skip(1) {
            zp *= zsq;
            let term = T::from_f64(*coeff) * zp;
            if term.abs() < test {
                break;
            }
            s += term;
        }
    }

    let tlg = z_modified.ln();
    let mut return_value =
        z_modified * (tlg - T::one()) + T::half() * (T::from_f64(LN_2_PI) - tlg) + s;

    if z_increment != 0 {
        let mut zp = T::one();
        for i in 0..z_increment {
            zp *= z + T::from_usize(i);
        }
        return_value -= zp.ln();
    }

    Ok(return_value)
}

// ln(2*PI)carg
const LN_2_PI: f64 = 1.837_877_066_409_345_6;

// COEFFICIENTS OF ASYMPTOTIC EXPANSION
const COEFFS: [f64; 22] = [
    8.333_333_333_333_333e-2,
    -2.777_777_777_777_778e-3,
    7.936_507_936_507_937e-4,
    -5.952_380_952_380_953e-4,
    8.417_508_417_508_417e-4,
    -1.917_526_917_526_917_6e-3,
    6.410_256_410_256_41e-3,
    -2.955_065_359_477_124_2e-2,
    1.796_443_723_688_305_7e-1,
    -1.392_432_216_905_901_1,
    1.340_286_404_416_839_3e1,
    -1.568_482_846_260_020_3e2,
    2.193_103_333_333_333_5e3,
    -3.610_877_125_372_499e4,
    6.914_722_688_513_13e5,
    -1.523_822_153_940_741_5e7,
    3.829_007_513_914_141_7e8,
    -1.088_226_603_578_439_1e10,
    3.473_202_837_650_022_6e11,
    -1.236_960_214_226_927_5e13,
    4.887_880_647_930_793e14,
    -2.132_033_396_091_937_2e16,
];

const INT_GAMMA_LOG_DATA: [f64; 101] = [
    0.0, // value at GLN[0] - can never be reached.
    0.00000000000000000e+00,
    0.00000000000000000e+00,
    f64::consts::LN_2,
    1.791_759_469_228_055,
    3.178_053_830_347_945_8,
    4.787_491_742_782_046,
    6.579_251_212_010_101,
    8.525_161_361_065_415,
    1.060_460_290_274_525e1,
    1.280_182_748_008_146_9e1,
    1.510_441_257_307_551_6e1,
    1.750_230_784_587_388_7e1,
    1.998_721_449_566_188_5e1,
    2.255_216_385_312_342_5e1,
    2.519_122_118_273_868e1,
    2.789_927_138_384_089e1,
    3.067_186_010_608_067_2e1,
    3.350_507_345_013_689e1,
    3.639_544_520_803_305e1,
    3.933_988_418_719_949_5e1,
    4.233_561_646_075_348_5e1,
    4.538_013_889_847_691e1,
    4.847_118_135_183_523e1,
    5.160_667_556_776_438e1,
    5.478_472_939_811_232e1,
    5.800_360_522_298_052e1,
    6.126_170_176_100_2e1,
    6.455_753_862_700_634e1,
    6.788_974_313_718_154e1,
    7.125_703_896_716_801e1,
    7.465_823_634_883_016e1,
    7.809_222_355_331_53e1,
    8.155_795_945_611_504e1,
    8.505_446_701_758_152e1,
    8.858_082_754_219_768e1,
    9.213_617_560_368_71e1,
    9.571_969_454_214_32e1,
    9.933_061_245_478_743e1,
    1.029_681_986_145_138_1e2,
    1.066_317_602_606_434_6e2,
    1.103_206_397_147_573_9e2,
    1.140_342_117_814_617e2,
    1.177_718_813_997_450_7e2,
    1.215_330_815_154_386_4e2,
    1.253_172_711_493_569e2,
    1.291_239_336_391_272_2e2,
    1.329_525_750_356_163_2e2,
    1.368_027_226_373_263_5e2,
    1.406_739_236_482_342_5e2,
    1.445_657_439_463_449e2,
    1.484_777_669_517_730_2e2,
    1.524_095_925_844_973_5e2,
    1.563_608_363_030_788e2,
    1.603_311_282_166_309e2,
    1.643_201_122_631_951_7e2,
    1.683_274_454_484_276_5e2,
    1.723_527_971_391_628e2,
    1.763_958_484_069_973_5e2,
    1.804_562_914_175_437_8e2,
    1.845_338_288_614_494_8e2,
    1.886_281_734_236_716e2,
    1.927_390_472_878_449e2,
    1.968_661_816_728_9e2,
    2.010_093_163_992_815_2e2,
    2.051_681_994_826_412e2,
    2.093_425_867_525_368_5e2,
    2.135_322_414_945_632_7e2,
    2.177_369_341_139_542_2e2,
    2.219_564_418_191_303_3e2,
    2.261_905_483_237_276e2,
    2.304_390_435_657_769_6e2,
    2.347_017_234_428_182_6e2,
    2.389_783_895_618_343_2e2,
    2.432_688_490_029_827e2,
    2.475_729_140_961_868_8e2,
    2.518_904_022_097_232e2,
    2.562_211_355_500_095_4e2,
    2.605_649_409_718_632e2,
    2.649_216_497_985_528e2,
    2.692_910_976_510_198e2,
    2.736_731_242_856_937e2,
    2.780_675_734_403_661e2,
    2.824_742_926_876_304e2,
    2.868_931_332_954_27e2,
    2.913_239_500_942_703e2,
    2.957_666_013_507_606_5e2,
    3.002_209_486_470_141_5e2,
    3.046_868_567_656_687e2,
    3.091_641_935_801_469e2,
    3.136_528_299_498_790_5e2,
    3.181_526_396_202_093e2,
    3.226_634_991_267_261_5e2,
    3.271_852_877_037_752e2,
    3.317_178_871_969_285e2,
    3.362_611_819_791_984_5e2,
    3.408_150_588_707_99e2,
    3.453_794_070_622_668_6e2,
    3.499_541_180_407_702_5e2,
    3.545_390_855_194_408e2,
    3.591_342_053_695_754e2,
];
