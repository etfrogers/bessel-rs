use alloc::vec::Vec;
use core::{
    fmt::Debug,
    ops::{AddAssign, Deref, DerefMut, Div, DivAssign, Mul, MulAssign, RemAssign, SubAssign},
};
use num::{
    Complex, Float,
    traits::{ConstOne, ConstZero, FloatConst},
};
use thiserror::Error;

#[cfg(feature = "std")]
use std::sync::LazyLock;

use crate::amos::{MACHINE_CONSTANTS_32, MACHINE_CONSTANTS_64, MachineConsts};

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

    /// Pre-computed value of `1.0 / 3.0` in this precision.
    const ONE_THIRD: Self;
    /// Pre-computed value of `2.0 / 3.0` in this precision.
    const TWO_THIRDS: Self;
    /// Pre-computed value of `0.5` in this precision.
    const HALF: Self;
    /// Pre-computed value of `2.0` in this precision.
    const TWO: Self;

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

    /// Returns the raw bitwise representation of this float.
    fn to_bits(self) -> u64;

    /// Creates a vector of length `n` containing complex zeros.
    #[inline]
    fn c_zeros(n: usize) -> Vec<Complex<Self>> {
        vec![Complex::<Self>::ZERO; n]
    }

    #[cfg(feature = "std")]
    /// Environmental machine constants used for scaling, underflow detection, and iteration bounds
    /// specific to the AMOS algorithms for this precision.
    const MACHINE_CONSTANTS: &'static LazyLock<MachineConsts<Self>>;
    #[cfg(not(feature = "std"))]
    /// Environmental machine constants used for scaling, underflow detection, and iteration bounds
    /// specific to the AMOS algorithms for this precision.
    const MACHINE_CONSTANTS: &MachineConsts<Self>;
}

impl BesselFloat for f64 {
    const RADIX: u32 = f64::RADIX;
    const MANTISSA_DIGITS: u32 = f64::MANTISSA_DIGITS;
    const MIN_EXP: i32 = f64::MIN_EXP;
    const MAX_EXP: i32 = f64::MAX_EXP;
    const EPSILON: Self = f64::EPSILON;
    const MIN_POSITIVE: Self = f64::MIN_POSITIVE;
    const NAN: Self = f64::NAN;

    const ONE_THIRD: Self = 1.0 / 3.0;
    const TWO_THIRDS: Self = 2.0 / 3.0;
    const HALF: Self = 0.5;
    const TWO: Self = 2.0;

    #[cfg(feature = "std")]
    const MACHINE_CONSTANTS: &'static LazyLock<MachineConsts<Self>> = &MACHINE_CONSTANTS_64;
    #[cfg(not(feature = "std"))]
    const MACHINE_CONSTANTS: &MachineConsts<Self> = &MACHINE_CONSTANTS_64;

    #[inline]
    fn from_f64(value: f64) -> Self {
        value
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

    const ONE_THIRD: Self = 1.0 / 3.0;
    const TWO_THIRDS: Self = 2.0 / 3.0;
    const HALF: Self = 0.5;
    const TWO: Self = 2.0;

    #[cfg(feature = "std")]
    const MACHINE_CONSTANTS: &'static LazyLock<MachineConsts<Self>> = &MACHINE_CONSTANTS_32;
    #[cfg(not(feature = "std"))]
    const MACHINE_CONSTANTS: &MachineConsts<Self> = &MACHINE_CONSTANTS_32;

    #[inline]
    fn from_f64(value: f64) -> Self {
        value as f32
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

/// Information about a computed Bessel or Hankel sequence.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct SequenceInfo {
    /// The number of components in the destination slice explicitly set to zero due to underflow.
    ///
    /// For $J_\nu$ and $I_\nu$, underflow zeroes occur at the end of the slice (highest orders).
    /// For $Y_\nu$, $K_\nu$, and $H_\nu^{(m)}$, underflow zeroes occur at the start of the slice.
    pub n_zeros: usize,
    /// Whether partial loss of significance occurred during calculation.
    ///
    /// When `true`, results in the destination slice have reduced accuracy
    /// (less than half of machine precision) due to large $|z|$ or `order`.
    pub partial_loss_of_significance: bool,
}

#[allow(type_alias_bounds)]
pub(crate) type BesselValues<FT: BesselFloat = f64, NT = SequenceInfo> = (Vec<Complex<FT>>, NT);

/// A trait for types that can be used as input to Bessel functions.
///
/// This trait is implemented for `f64`, `Complex<f64>`, `f32`, and `Complex<f32>`, allowing
/// the Bessel functions to accept both real and complex arguments.
///
/// This trait is sealed and cannot be implemented outside of `amos-bessel-rs`.
pub trait BesselInput<T: BesselFloat = f64>: Into<Complex<T>> + private::Sealed<T> {}

impl BesselInput<f64> for f64 {}
impl BesselInput<f64> for Complex<f64> {}
impl BesselInput<f32> for f32 {}
impl BesselInput<f32> for Complex<f32> {}

mod private {
    use num::{Complex, complex::ComplexFloat};

    use crate::{BesselError, BesselFloat, amos::MachineConsts};

    /// Private sealing trait used to convert calculation results back into `Self`
    /// (`T` or `Complex<T>`) - used to do a back conversion from Complex<T> to (real) T
    /// if the user has put a real argument into the high-level interface, to allow a real
    /// value to be returned.
    pub trait Sealed<T: BesselFloat>: Sized {
        /// Converts a `Complex<T>` into `Self`,
        /// verifying that the imaginary part is within machine tolerance when `Self = T`.
        fn back_from(val: Complex<T>) -> Result<Self, BesselError<T>>;
    }

    impl<T: BesselFloat> Sealed<T> for T {
        #[inline]
        fn back_from(val: Complex<T>) -> Result<Self, BesselError<T>> {
            let mc: &MachineConsts<T> = T::MACHINE_CONSTANTS;
            let margin = T::from_f64(1000.0);
            let tol = margin * mc.abs_error_tolerance;
            // if the imaginary part is small, pass the value on
            // if the imaginary part is small compared to the real part, pass the value on
            // if the real part is small, the imaginary part is likely inaccurate, so pass the value on
            if val.im().abs() < tol || val.im().abs() < val.re().abs() * tol || val.re().abs() < tol
            {
                Ok(val.re())
            } else {
                Err(BesselError::ComplexOutputForRealInput { output: val })
            }
        }
    }

    impl<T: BesselFloat> Sealed<T> for Complex<T> {
        #[inline]
        fn back_from(val: Self) -> Result<Self, BesselError<T>> {
            Ok(val)
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
/// nature of the error.
///
/// Implements `Copy` and produces zero heap allocations.
#[derive(Error, Debug, PartialEq, Clone, Copy)]
#[repr(i32)]
pub enum BesselError<T: BesselFloat = f64> {
    /// Indicates that the input is invalid (usually out of bounds) in some way.
    /// Documentation for each function lists valid and invalid inputs.
    #[error("Invalid input: {details}")]
    InvalidInput {
        /// Explanation of why the input was invalid.
        details: &'static str,
    } = 1,
    /// Overflow (or underflow) error in calculation: a valid answer cannot be calculated.
    /// Usually caused by a (very) large `order`, or small `z.abs()`.
    #[error("Overflow: order too large or z.abs() too small or both")]
    Overflow = 2,
    /// Complete loss of significance in output. No value could be calculated.
    #[error("Loss of too much significance in output")]
    LossOfSignificance = 4,
    /// Algorithm failed to converge to an answer.
    #[error("Algorithm failed to terminate")]
    DidNotConverge = 5,
    /// Returned only when the input `z` to the `bessel_...` functions is real.
    /// As these functions return a real output for a real input, the output is
    /// only valid if the imaginary part is small. If the imaginary part of the
    /// answer is significant this error is returned. The complex answer is returned
    /// in the output field, if that is wanted.
    #[error("Real input returned complex output. Output value {output}")]
    ComplexOutputForRealInput {
        /// Complex result of the Bessel function calculation.
        output: Complex<T>,
    } = 6,
}

impl<T: BesselFloat> BesselError<T> {
    /// A numeric form of the error equivalent to the error codes returned by the Amos
    /// Fortran code (where equivalence exists).
    pub fn error_code(&self) -> i32 {
        match self {
            BesselError::InvalidInput { .. } => 1,
            BesselError::Overflow => 2,
            BesselError::LossOfSignificance => 4,
            BesselError::DidNotConverge => 5,
            BesselError::ComplexOutputForRealInput { .. } => 6,
        }
    }

    #[doc(hidden)]
    pub fn from_i32(code: i32) -> Option<BesselError<T>> {
        match code {
            1 => Some(BesselError::InvalidInput {
                details: "from i32",
            }),
            2 => Some(BesselError::Overflow),
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
        match *self {
            BesselError::InvalidInput { details } => BesselError::InvalidInput { details },
            BesselError::Overflow => BesselError::Overflow,
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
        $base_func:ident
    ) => {
        paste! {
            $(#[$meta])*
            #[inline]
            fn [<$base_func _single>]<T: BesselFloat>(order: T, z: Complex<T>) -> Result<Complex<T>, BesselError<T>> {
                let mut buf = [T::C_ZERO; 1];
                [<complex_$base_func _into>](z, order, Scaling::Unscaled, &mut buf)?;
                Ok(buf[0])
            }
        }
    };
}

pub(crate) use simple_bessel_wrapper;

pub const DEFAULT_SBO_CAP: usize = 32;

pub enum ScratchBuffer<T: BesselFloat, const CAP: usize = DEFAULT_SBO_CAP> {
    Stack([Complex<T>; CAP], usize),
    Heap(Vec<Complex<T>>),
}

impl<T: BesselFloat> ScratchBuffer<T, DEFAULT_SBO_CAP> {
    #[inline]
    pub fn new(n: usize) -> Self {
        Self::with_capacity(n)
    }
}

impl<T: BesselFloat, const CAP: usize> ScratchBuffer<T, CAP> {
    #[inline]
    pub fn with_capacity(n: usize) -> Self {
        if n <= CAP {
            ScratchBuffer::Stack([T::C_ZERO; CAP], n)
        } else {
            ScratchBuffer::Heap(T::c_zeros(n))
        }
    }
}

impl<T: BesselFloat, const CAP: usize> Deref for ScratchBuffer<T, CAP> {
    type Target = [Complex<T>];
    #[inline]
    fn deref(&self) -> &Self::Target {
        match self {
            Self::Stack(arr, n) => &arr[..*n],
            Self::Heap(vec) => &vec[..],
        }
    }
}

impl<T: BesselFloat, const CAP: usize> DerefMut for ScratchBuffer<T, CAP> {
    #[inline]
    fn deref_mut(&mut self) -> &mut Self::Target {
        match self {
            Self::Stack(arr, n) => &mut arr[..*n],
            Self::Heap(vec) => &mut vec[..],
        }
    }
}
