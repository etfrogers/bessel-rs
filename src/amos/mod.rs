use std::f64::consts::PI;

pub use gamma_ln::{GammaError, gamma_ln};
pub use i_power_series::i_power_series;
use machine::{d1mach, i1mach};
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
// mod zbesh;

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

pub(crate) struct MachineConsts {
    arm: f64,
    ascle: f64,
    tol: f64,
    elim: f64,
    alim: f64,
    dig: f64,
    rl: f64,
    fnul: f64,
    rtol: f64,
}

impl MachineConsts {
    fn new() -> Self {
        let arm = 1.0e+3 * d1mach(1);
        let tol = d1mach(4).max(1.0e-18);
        let ascle = arm / tol;

        let k1 = i1mach(15);
        let k2 = i1mach(16);
        let r1m5 = d1mach(5);
        let k = k1.abs().min(k2.abs());
        let elim = 2.303 * ((k as f64) * r1m5 - 3.0);
        let k1 = i1mach(14) - 1;
        let mut aa = r1m5 * (k1 as f64);
        let dig = aa.min(18.0);
        aa *= 2.303;
        let alim = elim + (-aa).max(-41.45);
        let rl = 1.2 * dig + 3.0;
        let fnul = 10.0 + 6.0 * (dig - 3.0);

        Self {
            arm,
            ascle,
            tol,
            elim,
            alim,
            dig,
            rl,
            fnul,
            rtol: 1.0 / tol,
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
