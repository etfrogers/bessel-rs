use std::{
    collections::HashMap,
    fmt::Debug,
    ops::{AddAssign, Div, DivAssign, Mul, MulAssign, RemAssign, SubAssign},
    sync::{LazyLock, Mutex},
};

use crate::amos::{MACHINE_CONSTANTS_32, MACHINE_CONSTANTS_64, MachineConsts};
use num::{
    Complex, Float,
    complex::{Complex64, ComplexFloat},
    traits::{ConstOne, ConstZero, FloatConst},
};
use thiserror::Error;

pub(crate) trait BesselFloat:
    Float
    + CachableUAP<Self>
    + Debug
    + FloatConst
    + ConstZero
    + ConstOne
    + MulAssign
    + AddAssign
    + SubAssign
    + DivAssign
    + RemAssign
    + Mul<Complex<Self>, Output = Complex<Self>>
    + Div<Complex<Self>, Output = Complex<Self>>
    + 'static
// where Complex<Self>: MulAssign<Complex<Self>>
// where Complex<Self>: Pow<Self, Output = Complex<Self>>
{
    const RADIX: u32;
    const MANTISSA_DIGITS: u32;
    const MIN_EXP: i32;
    const MAX_EXP: i32;
    const EPSILON: Self;
    const C_ONE: Complex<Self> = Complex::<Self>::ONE;
    const C_ZERO: Complex<Self> = Complex::<Self>::ZERO;
    const I: Complex<Self> = Complex::<Self>::I;

    fn from_f64(value: f64) -> Self;
    fn half() -> Self;
    fn two() -> Self;
    fn to_bits(self) -> u64;

    #[inline]
    fn c_zeros(n: usize) -> Vec<Complex<Self>> {
        vec![Complex::<Self>::ZERO; n]
    }

    const MACHINE_CONSTANTS: &'static LazyLock<MachineConsts<Self>>;
}

impl BesselFloat for f64 {
    const RADIX: u32 = f64::RADIX;
    const MANTISSA_DIGITS: u32 = f64::MANTISSA_DIGITS;
    const MIN_EXP: i32 = f64::MIN_EXP;
    const MAX_EXP: i32 = f64::MAX_EXP;
    const EPSILON: Self = f64::EPSILON;

    const MACHINE_CONSTANTS: &'static LazyLock<MachineConsts<Self>> = &MACHINE_CONSTANTS_64;

    #[inline]
    fn from_f64(value: f64) -> Self {
        value
    }

    #[inline]
    fn half() -> Self {
        0.5
    }

    #[inline]
    fn two() -> Self {
        2.0
    }

    #[inline]
    fn to_bits(self) -> u64 {
        f64::to_bits(self)
    }
}

impl BesselFloat for f32 {
    const RADIX: u32 = f32::RADIX;
    const MANTISSA_DIGITS: u32 = f32::MANTISSA_DIGITS;
    const MIN_EXP: i32 = f32::MIN_EXP;
    const MAX_EXP: i32 = f32::MAX_EXP;
    const EPSILON: Self = f32::EPSILON;

    const MACHINE_CONSTANTS: &'static LazyLock<MachineConsts<Self>> = &MACHINE_CONSTANTS_32;

    #[inline]
    fn from_f64(value: f64) -> Self {
        value as f32
    }

    #[inline]
    fn half() -> Self {
        0.5
    }

    #[inline]
    fn two() -> Self {
        2.0
    }

    #[inline]
    fn to_bits(self) -> u64 {
        f32::to_bits(self) as u64
    }
}

pub(crate) trait CachableUAP<T: BesselFloat>: 'static {
    const UNIFORM_ASSYMPTOTIC_PARAMETERS_CACHE: &'static LazyLock<
        Mutex<HashMap<CacheKey, UniformAssymptoticParameters<T>>>,
    >;
}

impl CachableUAP<f64> for f64 {
    const UNIFORM_ASSYMPTOTIC_PARAMETERS_CACHE: &'static LazyLock<
        Mutex<HashMap<CacheKey, UniformAssymptoticParameters<Self>>>,
    > = &UNIFORM_ASSYMPTOTIC_PARAMETERS_CACHE_64;
}

impl CachableUAP<f32> for f32 {
    const UNIFORM_ASSYMPTOTIC_PARAMETERS_CACHE: &'static LazyLock<
        Mutex<HashMap<CacheKey, UniformAssymptoticParameters<Self>>>,
    > = &UNIFORM_ASSYMPTOTIC_PARAMETERS_CACHE_32;
}

pub(crate) struct UniformAssymptoticParameters<T: BesselFloat> {
    pub(crate) phi_i: Complex<T>,
    pub(crate) phi_k: Complex<T>,
    pub(crate) zeta1: Complex<T>,
    pub(crate) zeta2: Complex<T>,
    pub(crate) sum_i: Option<Complex<T>>,
    pub(crate) sum_k: Option<Complex<T>>,
    pub(crate) working: Option<Vec<Complex<T>>>,
}

pub(crate) type CacheKey = (u64, u64, u64);

static UNIFORM_ASSYMPTOTIC_PARAMETERS_CACHE_64: LazyLock<
    Mutex<HashMap<CacheKey, UniformAssymptoticParameters<f64>>>,
> = LazyLock::new(|| Mutex::new(HashMap::new()));

static UNIFORM_ASSYMPTOTIC_PARAMETERS_CACHE_32: LazyLock<
    Mutex<HashMap<CacheKey, UniformAssymptoticParameters<f32>>>,
> = LazyLock::new(|| Mutex::new(HashMap::new()));

pub(crate) fn cache_key<T: BesselFloat>(z: Complex<T>, order: T) -> CacheKey {
    (z.re.to_bits(), z.im.to_bits(), order.to_bits())
}

#[allow(type_alias_bounds)]
pub(crate) type BesselValues<FT: BesselFloat = f64, NT = usize> = (Vec<Complex<FT>>, NT);
pub(crate) type BesselResult<T = BesselValues> = Result<T, BesselError>;

/// A trait for converting back from a type `T` into a `BesselResult<Self>`.
/// Used for allowing both real and complex inputs to the Bessel functions,
/// with real output for real input (provided the answer is real)
/// Implemented for `f64` and `Complex<f64>`
pub trait BackFrom<T>: Sized {
    /// Converts a value (number or `Result<number, BesselError>` wrapping that number)
    /// into the number itself (wrapped in a `Result<number, BesselError>`, due to the
    /// possibility of complex output for real input)
    fn back_from(val: &T) -> BesselResult<Self>;
}

impl BackFrom<Complex64> for Complex64 {
    #[inline]
    fn back_from(val: &Complex64) -> BesselResult<Self> {
        Ok(*val)
    }
}

impl BackFrom<f64> for f64 {
    #[inline]
    fn back_from(val: &f64) -> BesselResult<Self> {
        Ok(*val)
    }
}

impl BackFrom<Complex64> for f64 {
    #[inline]
    fn back_from(val: &Complex64) -> BesselResult<Self> {
        const MARGIN: f64 = 1000.0;
        let tol = MARGIN * f64::MACHINE_CONSTANTS.abs_error_tolerance;
        // if the imainary part is small, pass the value on
        // if the imaginary part is small compared to the real part, pass the value on
        // if the real part is small, the imaginary part is likely inaccurate, so pass the value on
        if val.im().abs() < tol || val.im().abs() < val.re().abs() * tol || val.re() < tol {
            Ok(val.re())
        } else {
            Err(BesselError::ComplexOutputForRealInput { output: *val })
        }
    }
}

impl BackFrom<BesselResult<Complex64>> for f64 {
    #[inline]
    fn back_from(val: &BesselResult<Complex<f64>>) -> BesselResult<Self> {
        match val {
            Ok(cpx) => f64::back_from(cpx),
            // below we can assume that y has one element, as the input type is BesselResult<Complex<f64>> not
            // BesselResult<Vec<Complex<f64>>>
            Err(BesselError::PartialLossOfSignificance { y, nz: _ }) => f64::back_from(&y[0]),
            Err(err) => Err((*err).clone()),
        }
    }
}

impl BackFrom<BesselResult<Complex64>> for Complex<f64> {
    #[inline]
    fn back_from(val: &BesselResult<Complex<f64>>) -> BesselResult<Self> {
        match val {
            Ok(cpx) => Ok(*cpx),
            // below we can assume that y has one element, as the input type is BesselResult<Complex<f64>> not
            // BesselResult<Vec<Complex<f64>>>
            Err(BesselError::PartialLossOfSignificance { y, nz: _ }) => Ok(y[0]),
            Err(err) => Err((*err).clone()),
        }
    }
}

// Original Docs:
// IERR   - ERROR FLAG
//         IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
//         IERR=1, INPUT ERROR   - NO COMPUTATION
//         IERR=2, OVERFLOW      - NO COMPUTATION, order TOO
//                 LARGE OR CABS(Z) TOO SMALL OR BOTH
//         IERR=3, CABS(Z) OR order+N-1 LARGE - COMPUTATION DONE
//                BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
//                 REDUCTION PRODUCE LESS THAN HALF OF MACHINE
//                 ACCURACY
//         IERR=4, CABS(Z) OR order+N-1 TOO LARGE - NO COMPUTA-
//                 TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
//                 CANCE BY ARGUMENT REDUCTION
//         IERR=5, ERROR              - NO COMPUTATION,
//                 ALGORITHM TERMINATION CONDITION NOT MET
/// Error struct returned by Bessel function calculations indicating the
/// nature of the error
#[derive(Error, Debug, PartialEq, Clone)]
#[repr(i32)]
pub enum BesselError {
    /// Indicates that the input is invalid (usually out of bounds) in some way.
    /// Documentation for each function lists valid and invalid inputs
    #[error("Invalid input: {details}")]
    InvalidInput {
        /// Explanation of why the input was invalid.
        details: String,
    } = 1,
    /// Overflow (or underflow) error in calculation: a valid answer cannot be calculated
    /// Usually caused by a (very) large `order`, or small `z.abs()`.
    #[error("Overflow: order too large or z.abs() too small or both")]
    Overflow = 2, //{ too_large: bool },
    /// Calculation is done, and a valkue returned wrapped in this error,
    /// however the value is lower in accuracy than normally expected from these algorithms.
    /// As `z.abs()` or `order` are large, losses of significnace produce
    /// less than half of machine accuracy. This error is conservative, in
    /// that it assume argument reduction causes problems that may not occur
    /// in some architectures.
    /// Not returned by the reduced API `bessel_...` functions, as they unwrap this and
    /// return the value. To detect partial loss of significance, the `complex_bessel_..`
    /// function must be used.
    #[error("Partial loss of significance in output. Losssy values returned.")]
    PartialLossOfSignificance {
        /// Value(s) of Bessel function (reduced accuracy)
        y: Vec<Complex64>,
        /// Number of entries in `y` explicitly set to zero (as per the `complex_bessel_...` docs`)
        nz: usize,
    } = 3,
    #[error("Loss of too much significance in output")]
    /// Complete loss of significance in output. No value could be calculated
    LossOfSignificance = 4,
    /// Algorithm failed to converge to a an answer
    #[error("Algorithm failed to terminate")]
    DidNotConverge = 5,
    /// Returned only when the input `z` to the `bessel_...` functions is real.
    /// As these function return a real output for a real input, the output is
    /// only valid if the imaginary part is small. If the imaginary part of the
    /// answer is siginificant this error is returned. The complex answer is returned
    /// in the output field, if that is wanted.
    #[error("Real input returned complex output. Output value {output}")]
    ComplexOutputForRealInput {
        /// Complex result of the Bessel function calculation.
        output: Complex<f64>,
    } = 6,
}

impl BesselError {
    /// A numeric form of the error equivalent to the error codes returned by the Amos
    /// Fortran code (where equivalence exists)
    pub fn error_code(&self) -> i32 {
        match self {
            BesselError::InvalidInput { .. } => 1,
            BesselError::Overflow => 2,
            BesselError::PartialLossOfSignificance { .. } => 3,
            BesselError::LossOfSignificance => 4,
            BesselError::DidNotConverge => 5,
            BesselError::ComplexOutputForRealInput { .. } => 6,
        }
    }
}

macro_rules! simple_bessel_wrapper {
    (
        $(#[$meta:meta])*
        $base_func:ident // We only pass the base function name now!
    ) => {
        // The paste! macro allows us to create new identifiers
        paste! {
            $(#[$meta])*
            // [<simple_ $base_func>] concatenates into simple_bessel_j
            #[inline]
            fn [<$base_func _single>](order: f64, z: Complex64) -> Result<Complex64, BesselError> {
                let (result_vec, _nz) = [<complex_$base_func>](z, order, Scaling::Unscaled, 1)?;
                Ok(result_vec[0])
            }
        }
    };
}

pub(crate) use simple_bessel_wrapper;
