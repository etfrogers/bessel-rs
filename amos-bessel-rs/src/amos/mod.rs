use num::{
    Complex, Float,
    complex::{Complex64, ComplexFloat},
};
use std::{f64::consts::PI, ops::Neg};

pub use entry_points::*;
pub(crate) use gamma_ln::gamma_ln;
pub(crate) use i_power_series::i_power_series;
pub(crate) use machine::{MACHINE_CONSTANTS_32, MACHINE_CONSTANTS_64, MachineConsts};

#[cfg(test)]
pub(crate) use gamma_ln::GammaError;

use crate::types::BesselFloat;

mod asymptotic_i;
mod entry_points;
mod gamma_ln;
mod i_power_series;
mod machine;
mod overflow_checks;
mod translator;
mod utils;

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
#[repr(usize)]
pub(crate) enum IKType {
    I = 1,
    K = 2,
}

const CIP: [Complex64; 4] = [
    Complex64::new(1.0, 0.0),
    Complex64::new(0.0, 1.0),
    Complex64::new(-1.0, 0.0),
    Complex64::new(0.0, -1.0),
];

/// Used to specify the kind of Hankel function in the [hankel](crate::hankel) and
/// [complex_bessel_h] functions.
#[derive(Debug, PartialEq, Eq, Clone, Copy)]
#[repr(usize)]
pub enum HankelKind {
    /// Hankel function of the first kind, H1v(z) = Jv(z) + i*Yv(z)
    First = 1,
    /// Hankel function of the second kind, H2v(z) = Jv(z) - i*Yv(z)
    Second = 2,
}

impl HankelKind {
    pub(crate) fn get_rotation(&self) -> RotationDirection {
        match self {
            HankelKind::First => RotationDirection::Right,
            HankelKind::Second => RotationDirection::Left,
        }
    }
}

/// Represents the scaling option for Bessel and Airy complex_... functions.
#[derive(Debug, PartialEq, Eq, Clone, Copy)]
#[repr(i32)]
pub enum Scaling {
    /// No scaling is applied.
    Unscaled = 1,
    /// Scaling is applied to remove exponential growth or decay.
    Scaled = 2,
}

impl Scaling {
    pub(crate) fn scale_zetas<T: BesselFloat>(
        &self,
        z: Complex<T>,
        modified_order: T,
        zeta1: Complex<T>,
        zeta2: Complex<T>,
    ) -> Complex<T> {
        match self {
            Scaling::Unscaled => -zeta1 + zeta2,
            Scaling::Scaled => {
                let mut st = z + zeta2;
                st = st.conj() * (modified_order / st.abs()).powi(2);
                -zeta1 + st
            }
        }
    }
}

#[doc(hidden)]
#[derive(Debug, PartialEq, Eq, Clone, Copy)]
#[repr(i32)]
pub(crate) enum RotationDirection {
    Left = -1,
    None = 0,
    Right = 1,
}

impl RotationDirection {
    pub fn signum(&self) -> f64 {
        (*self as i32 as f64).signum()
    }

    #[inline]
    pub fn to_float<T: BesselFloat>(self) -> T {
        T::from_f64(self as i32 as f64)
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

pub(crate) fn max_abs_component<T: Float>(c: Complex<T>) -> T {
    c.re.abs().max(c.im.abs())
}

pub(crate) trait PositiveArg<T> {
    fn parg(&self) -> T;
}

impl<T: BesselFloat> PositiveArg<T> for Complex<T> {
    fn parg(&self) -> T {
        let mut ang = self.arg();
        if ang < T::zero() {
            ang += T::from_f64(PI * 2.0);
        }
        ang
    }
}
