use std::{
    fmt::Debug,
    ops::{AddAssign, Div, DivAssign, Mul, MulAssign, RemAssign, SubAssign},
    sync::LazyLock,
};

use crate::amos::{MACHINE_CONSTANTS_32, MACHINE_CONSTANTS_64, MachineConsts};
use num::{
    Complex, Float,
    complex::ComplexFloat,
    traits::{ConstOne, ConstZero, FloatConst},
};
use thiserror::Error;

/// A trait defining the mathematical and floating-point constraints required to compute 
/// Bessel and Airy functions.
/// 
/// This trait is implemented for `f32` and `f64`. Downstream users can use this trait 
/// as a generic bound to write their own generic functions over real or complex Bessel evaluations.
pub trait BesselFloat:
    Float
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
    + PartialOrd
    + 'static
{
    /// The radix or base of the internal representation of this type.
    const RADIX: u32;
    /// The number of significant digits in base-`RADIX` for this type.
    const MANTISSA_DIGITS: u32;
    /// The lowest possible power of 2 exponent for a valid normalized float.
    const MIN_EXP: i32;
    /// The highest possible power of 2 exponent for a valid normalized float.
    const MAX_EXP: i32;
    /// The difference between `1.0` and the next larger representable number.
    const EPSILON: Self;
    /// The smallest positive normal floating point number.
    const MIN_POSITIVE: Self;
    /// Not a Number (NaN).
    const NAN: Self;
    /// Pre-computed value of `2.0 / 3.0` in this precision.
    const TWO_THIRDS: Self;
    /// The complex number `1.0 + 0.0i`.
    const C_ONE: Complex<Self> = Complex::<Self>::ONE;
    /// The complex number `0.0 + 0.0i`.
    const C_ZERO: Complex<Self> = Complex::<Self>::ZERO;
    /// The complex imaginary unit `0.0 + 1.0i`.
    const I: Complex<Self> = Complex::<Self>::I;

    /// Casts an `f64` value to this type.
    fn from_f64(value: f64) -> Self;
    /// Casts a `usize` value to this type.
    fn from_usize(value: usize) -> Self;
    /// Casts an `isize` value to this type.
    fn from_isize(value: isize) -> Self;
    /// Casts a `Complex<f64>` value to a `Complex` of this type.
    fn from_cpx64(value: Complex<f64>) -> Complex<Self>;
    
    /// Pre-computed value of `0.5` in this precision.
    fn half() -> Self;
    /// Pre-computed value of `2.0` in this precision.
    fn two() -> Self;
    
    /// Returns the raw bitwise representation of this float.
    fn to_bits(self) -> u64;

    /// Creates a vector of length `n` containing complex zeros.
    #[inline]
    fn c_zeros(n: usize) -> Vec<Complex<Self>> {
        vec![Complex::<Self>::ZERO; n]
    }

    /// Environmental machine constants used for scaling, underflow detection, and iteration bounds 
    /// specific to the AMOS algorithms for this precision.
    const MACHINE_CONSTANTS: &'static LazyLock<MachineConsts<Self>>;
}

impl BesselFloat for f64 {
    const RADIX: u32 = f64::RADIX;
    const MANTISSA_DIGITS: u32 = f64::MANTISSA_DIGITS;
    const MIN_EXP: i32 = f64::MIN_EXP;
    const MAX_EXP: i32 = f64::MAX_EXP;
    const EPSILON: Self = f64::EPSILON;
    const MIN_POSITIVE: Self = f64::MIN_POSITIVE;
    const NAN: Self = f64::NAN;
    const TWO_THIRDS: Self = 2.0 / 3.0;

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

    #[inline]
    fn from_cpx64(value: Complex<f64>) -> Complex<Self> {
        value
    }

    #[inline]
    fn from_usize(value: usize) -> Self {
        value as f64
    }

    #[inline]
    fn from_isize(value: isize) -> Self {
        value as f64
    }
}

impl BesselFloat for f32 {
    const RADIX: u32 = f32::RADIX;
    const MANTISSA_DIGITS: u32 = f32::MANTISSA_DIGITS;
    const MIN_EXP: i32 = f32::MIN_EXP;
    const MAX_EXP: i32 = f32::MAX_EXP;
    const EPSILON: Self = f32::EPSILON;
    const MIN_POSITIVE: Self = f32::MIN_POSITIVE;
    const NAN: Self = f32::NAN;
    const TWO_THIRDS: Self = 2.0 / 3.0;

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

    #[inline]
    fn from_cpx64(value: Complex<f64>) -> Complex<Self> {
        Complex::new(value.re as f32, value.im as f32)
    }

    #[inline]
    fn from_usize(value: usize) -> Self {
        value as f32
    }

    #[inline]
    fn from_isize(value: isize) -> Self {
        value as f32
    }
}

#[allow(type_alias_bounds)]
pub(crate) type BesselValues<FT: BesselFloat = f64, NT = usize> = (Vec<Complex<FT>>, NT);
#[allow(type_alias_bounds)]
pub(crate) type BesselResult<FT: BesselFloat = f64, NT = usize> =
    Result<BesselValues<FT, NT>, BesselError<FT>>;

/// A trait for converting back from a type `T` into a `BesselResult<Self>`.
/// Used for allowing both real and complex inputs to the Bessel functions,
/// with real output for real input (provided the answer is real)
/// Implemented for `f64` and `Complex<f64>`
pub trait BackFrom<T, FT: BesselFloat>: Sized {
    /// Converts a value (number or `Result<number, BesselError>` wrapping that number)
    /// into the number itself (wrapped in a `Result<number, BesselError>`, due to the
    /// possibility of complex output for real input)
    fn back_from(val: &T) -> Result<Self, BesselError<FT>>;
}

impl<T: BesselFloat> BackFrom<Complex<T>, T> for Complex<T> {
    #[inline]
    fn back_from(val: &Complex<T>) -> Result<Complex<T>, BesselError<T>> {
        Ok(*val)
    }
}

impl<T: BesselFloat> BackFrom<T, T> for T {
    #[inline]
    fn back_from(val: &T) -> Result<T, BesselError<T>> {
        Ok(*val)
    }
}

impl<T: BesselFloat> BackFrom<Complex<T>, T> for T {
    #[inline]
    fn back_from(val: &Complex<T>) -> Result<T, BesselError<T>> {
        let margin = T::from_f64(1000.0);
        let tol = margin * T::MACHINE_CONSTANTS.abs_error_tolerance;
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

impl<T: BesselFloat> BackFrom<Result<Complex<T>, BesselError<T>>, T> for T {
    #[inline]
    fn back_from(
        val: &Result<Complex<Self>, BesselError<Self>>,
    ) -> Result<Self, BesselError<Self>> {
        match val {
            Ok(cpx) => T::back_from(cpx),
            // below we can assume that y has one element, as the input type is BesselResult<Complex<T>> not
            // BesselResult<Vec<Complex<T>>>
            Err(BesselError::PartialLossOfSignificance { y, n_zeros: _ }) => T::back_from(&y[0]),
            Err(err) => Err((*err).clone()),
        }
    }
}

impl<T: BesselFloat> BackFrom<Result<Self, BesselError<T>>, T> for Complex<T> {
    #[inline]
    fn back_from(val: &Result<Self, BesselError<T>>) -> Result<Self, BesselError<T>> {
        match val {
            Ok(cpx) => Ok(*cpx),
            // below we can assume that y has one element, as the input type is BesselResult<Complex<T>> not
            // BesselResult<Vec<Complex<T>>>
            Err(BesselError::PartialLossOfSignificance { y, n_zeros: _ }) => Ok(y[0]),
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
pub enum BesselError<T: BesselFloat = f64> {
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
        y: Vec<Complex<T>>,
        /// Number of entries in `y` explicitly set to zero (as per the `complex_bessel_...` docs`)
        n_zeros: usize,
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
        output: Complex<T>,
    } = 6,
}

impl<T: BesselFloat> BesselError<T> {
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

    #[doc(hidden)]
    pub fn from_i32(code: i32) -> Option<BesselError<T>> {
        match code {
            1 => Some(BesselError::InvalidInput {
                details: "from i32".to_string(),
            }),
            2 => Some(BesselError::Overflow),
            3 => Some(BesselError::PartialLossOfSignificance {
                y: vec![],
                n_zeros: 0,
            }),
            4 => Some(BesselError::LossOfSignificance),
            5 => Some(BesselError::DidNotConverge),
            6 => Some(BesselError::ComplexOutputForRealInput {
                output: Complex::new(T::NAN, T::NAN),
            }),
            _ => None,
        }
    }

    #[doc(hidden)]
    pub fn to_f32(&self) -> BesselError<f32> {
        match self {
            BesselError::InvalidInput { details } => BesselError::InvalidInput {
                details: details.clone(),
            },
            BesselError::Overflow => BesselError::Overflow,
            BesselError::PartialLossOfSignificance { y, n_zeros } => {
                BesselError::PartialLossOfSignificance {
                    y: y.iter()
                        .map(|c| Complex::new(c.re.to_f32().unwrap(), c.im.to_f32().unwrap()))
                        .collect(),
                    n_zeros: *n_zeros,
                }
            }
            BesselError::LossOfSignificance => BesselError::LossOfSignificance,
            BesselError::DidNotConverge => BesselError::DidNotConverge,
            BesselError::ComplexOutputForRealInput { output } => {
                BesselError::ComplexOutputForRealInput {
                    output: Complex::new(output.re.to_f32().unwrap(), output.im.to_f32().unwrap()),
                }
            }
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
            fn [<$base_func _single>]<T:BesselFloat>(order: T, z: Complex<T>) -> Result<Complex<T>, BesselError<T>> {
                let (result_vec, _n_zeros) = [<complex_$base_func>](z, order, Scaling::Unscaled, 1)?;
                Ok(result_vec[0])
            }
        }
    };
}

pub(crate) use simple_bessel_wrapper;
