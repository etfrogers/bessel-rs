use std::f64::consts::PI;

pub use gamma_ln::{GammaError, gamma_ln};
pub use i_power_series::i_power_series;
use num::{Complex, One, Zero, complex::Complex64};
use thiserror::Error;
pub use translator::zbesj;
pub mod bindings;
mod gamma_ln;
mod i_power_series;
mod machine;
mod overflow_checks;
mod translator;
mod utils;
mod z_asymptotic_i;

#[derive(Error, Debug, PartialEq)]
#[repr(i32)]
pub enum BesselError {
    // This correpsonds IERRR
    // 0 = no error
    #[error("Invalid input: {details}")]
    InvalidInput { details: String } = 1,
    #[error("Overflow: order TOO LARGE OR CABS(Z) TOO SMALL OR BOTH")]
    Overflow = 2, //{ too_large: bool },
    #[error("Partial loss of significance in output. Losssy values returned.")]
    PartialLossOfSignificance { y: Vec<Complex64>, nz: usize },
    // IERR = 3 is a warning (and hence some return value, I think) and needs handling elsewhere
    #[error("Loss of too much significance in output")]
    LossOfSignificance = 4,
    #[error("Algorithm failed to terminate")]
    DidNotConverge = 5,
    // Original Docs:
    // IERR   - ERROR FLAG
    //         IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
    //         IERR=1, INPUT ERROR   - NO COMPUTATION
    //         IERR=2, OVERFLOW      - NO COMPUTATION, order TOO
    //                 LARGE OR CABS(Z) TOO SMALL OR BOTH
    //         IERR=3, CABS(Z) OR order+N-1 LARGE - COMPUTATION DONE
    //                 BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
    //                 REDUCTION PRODUCE LESS THAN HALF OF MACHINE
    //                 ACCURACY
    //         IERR=4, CABS(Z) OR order+N-1 TOO LARGE - NO COMPUTA-
    //                 TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
    //                 CANCE BY ARGUMENT REDUCTION
    //         IERR=5, ERROR              - NO COMPUTATION,
    //                 ALGORITHM TERMINATION CONDITION NOT MET
    #[error("not yet implemented")]
    NotYetImplemented = 100,
}

impl BesselError {
    pub fn error_code(&self) -> i32 {
        match self {
            BesselError::InvalidInput { .. } => 1,
            BesselError::Overflow => 2,
            BesselError::PartialLossOfSignificance { .. } => 3,
            BesselError::LossOfSignificance => 4,
            BesselError::DidNotConverge => 5,
            BesselError::NotYetImplemented => 100,
        }
    }
}

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
#[repr(usize)]
pub enum IKType {
    I = 1,
    K = 2,
}

impl From<IKType> for usize {
    fn from(value: IKType) -> Self {
        match value {
            IKType::I => 1,
            IKType::K => 2,
        }
    }
}
impl From<&IKType> for usize {
    fn from(value: &IKType) -> Self {
        (*value).into()
    }
}

impl IKType {
    pub fn index(&self) -> usize {
        self.into()
    }
}

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub enum HankelKind {
    First,
    Second,
}

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
#[repr(i32)]
pub enum Scaling {
    Unscaled = 1,
    Scaled = 2,
}

pub(crate) type BesselValues<T = usize> = (Vec<Complex64>, T);
pub(crate) type BesselResult<T = BesselValues> = Result<T, BesselError>;

/// `elim` is a nunber such that if you take `elim.exp()` or `(-elim).exp()`, you will create a
/// number that is at risk of overflowing `f64`. In this case "at risk of overflow" means
/// within a factor of 1000.0 of `f64::MAX_POSITIVE` for `elim.exp()` or
/// `f64::MIN_POSITIVE` for `(-elim).exp()`.
///
/// As the two conditions (positive and negative) above are not exactly the same, the code chooses
/// the most conservative (which as standard is the MIN_POSITIVE version)
///
/// `approximation_limit` is a number a bit smaller than `elim` above which scaled calculations are used to
/// ensure that over/underflow is mitigated. Here "a bit smaller" means "reduced by a factor of
/// the smaller of 18 digits or the number of (decimal) digits stored in the mantissa of `f64`
/// (which as standard is ~15.6). For `x > approximation_limit`, e^x will start to lose precision.
///
/// # See also
/// tests in `test_machine_consts.rs`
#[derive(Debug, Clone)]
pub(crate) struct MachineConsts {
    // Originally ARM
    pub underflow_limit: f64,
    /// Originally ASCLE
    pub absolute_approximation_limit: f64,
    /// TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
    pub abs_error_tolerance: f64, // TOL
    /// Originally ELIM
    pub exponent_limit: f64,
    /// Originally ALIM
    pub approximation_limit: f64,
    /// DIG NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
    pub _significant_digits: f64, // DIG
    /// RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
    pub asymptotic_z_limit: f64, // RL
    /// FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU.
    pub asymptotic_order_limit: f64, // FNUL
    pub rtol: f64,
}

impl MachineConsts {
    pub fn new() -> Self {
        // Here we use approximate value, rather than calculating `10.0_f64.ln()`, as
        // this matches the the Fotran code, and the exact value causes subtle differences
        // in output (should just be what values are accpeted, but cause tetst to fail)
        let ln_10: f64 = 2.303;

        let underflow_limit = 1.0e+3 * 2.0 * f64::MIN_POSITIVE;
        let abs_error_tolerance = f64::EPSILON.max(1.0e-18);
        // absolute_approximation_limit == (-approximation_limit).exp() -- see test
        let absolute_approximation_limit = underflow_limit / abs_error_tolerance;

        let digits_per_bit = (f64::RADIX as f64).log10();
        let exponent_bit_limit = f64::MIN_EXP.abs().min(f64::MAX_EXP.abs());

        // Subtract 3.0 (digits) to give a number above which 10^decimal_exponent_limit would
        // be close to overflowing (i.e. within 1000 == 10^3.0 of the actual limit)
        let decimal_exponent_limit = (exponent_bit_limit as f64) * digits_per_bit - 3.0;
        // Multiplying by ln_10 converts from 10^x overflowing to e^x overflowing
        let exponent_limit = ln_10 * decimal_exponent_limit;

        let f64_siginficant_digits = digits_per_bit * ((f64::MANTISSA_DIGITS - 1) as f64);
        // siginficant_digits == abs_error_tolerance.log10() -- see test
        let significant_digits = f64_siginficant_digits.min(18.0);
        // Again, multiply number of base 10 digits by ln_10 to convert to e^x
        let approximation_limit = exponent_limit - (significant_digits * ln_10);

        let asymptotic_z_limit = 1.2 * significant_digits + 3.0;
        let asymptotic_order_limit = 10.0 + 6.0 * (significant_digits - 3.0);

        Self {
            underflow_limit,
            absolute_approximation_limit,
            abs_error_tolerance,
            exponent_limit,
            approximation_limit,
            _significant_digits: significant_digits,
            asymptotic_z_limit,
            asymptotic_order_limit,
            rtol: 1.0 / abs_error_tolerance,
        }
    }
}

pub(crate) fn c_one() -> Complex64 {
    Complex64::one()
}

#[inline]
pub(crate) fn c_zero() -> Complex64 {
    Complex64::zero()
}

#[inline]
pub(crate) fn c_zeros(n: usize) -> Vec<Complex64> {
    vec![Complex64::zero(); n]
}

pub(crate) fn max_abs_component(c: Complex64) -> f64 {
    c.re.abs().max(c.im.abs())
}

pub(crate) trait PositiveArg {
    fn parg(&self) -> f64;
}

impl PositiveArg for Complex<f64> {
    fn parg(&self) -> f64 {
        let mut ang = self.arg();
        if ang < 0.0 {
            ang = (PI * 2.0) + ang;
        }
        ang
    }
}
