use num::{
    Complex, One, Zero,
    complex::{Complex64, ComplexFloat},
    traits::Pow,
};
use std::{f64::consts::PI, ops::Neg};
use thiserror::Error;

pub use gamma_ln::{GammaError, gamma_ln};
pub(crate) use i_power_series::i_power_series;
pub(crate) use machine::MACHINE_CONSTANTS;
pub use translator::{complex_airy, complex_bessel_k, complex_bessel_y, complex_bessel_h, complex_bessel_j, complex_bessel_i};

pub(crate) mod bindings;
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
    // This correpsonds to IERR
    // 0 = no error
    #[error("Invalid input: {details}")]
    InvalidInput { details: String } = 1,
    #[error("Overflow: order TOO LARGE OR CABS(Z) TOO SMALL OR BOTH")]
    Overflow = 2, //{ too_large: bool },
    #[error("Partial loss of significance in output. Losssy values returned.")]
    PartialLossOfSignificance { y: Vec<Complex64>, nz: usize } = 3,
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

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
#[repr(usize)]
pub enum HankelKind {
    First = 1,
    Second = 2,
}

impl HankelKind {
    pub fn get_rotation(&self) -> RotationDirection {
        match self {
            HankelKind::First => RotationDirection::Right,
            HankelKind::Second => RotationDirection::Left,
        }
    }
}

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
#[repr(i32)]
pub enum Scaling {
    Unscaled = 1,
    Scaled = 2,
}

impl Scaling {
    pub fn scale_zetas(
        &self,
        z: Complex64,
        modified_order: f64,
        zeta1: Complex64,
        zeta2: Complex64,
    ) -> Complex64 {
        match self {
            Scaling::Unscaled => -zeta1 + zeta2,
            Scaling::Scaled => {
                let mut st = z + zeta2;
                st = st.conj() * (modified_order / st.abs()).pow(2);
                -zeta1 + st
            }
        }
    }
}

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
#[repr(i32)]
pub enum RotationDirection {
    Left = -1,
    None = 0,
    Right = 1,
}

impl RotationDirection {
    pub fn signum(&self) -> f64 {
        (*self as i32 as f64).signum()
    }
}

impl From<RotationDirection> for f64 {
    fn from(value: RotationDirection) -> Self {
        value as i32 as f64
    }
}

impl Neg for RotationDirection {
    type Output = Self;

    fn neg(self) -> Self::Output {
        match self {
            RotationDirection::Left => RotationDirection::Right,
            RotationDirection::None => RotationDirection::None,
            RotationDirection::Right => RotationDirection::Left,
        }
    }
}

pub(crate) type BesselValues<T = usize> = (Vec<Complex64>, T);
pub(crate) type BesselResult<T = BesselValues> = Result<T, BesselError>;

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
            ang += PI * 2.0;
        }
        ang
    }
}
