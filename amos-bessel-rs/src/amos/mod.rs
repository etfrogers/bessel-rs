use num::{Complex, Float, complex::ComplexFloat};
use std::{f64::consts::PI, ops::Neg};

pub use entry_points::*;
pub(crate) use gamma_ln::gamma_ln;
pub(crate) use machine::{MACHINE_CONSTANTS_32, MACHINE_CONSTANTS_64, MachineConsts};

#[cfg(test)]
pub(crate) use gamma_ln::GammaError;

use crate::types::BesselFloat;

mod airy;
mod analytic_continuation;
mod asymptotics;
mod entry_points;
mod gamma_ln;
mod limits;
mod machine;
mod power_series;
mod recurrence;
mod right_half_plane;
mod utils;
mod wronskian;

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
#[repr(usize)]
pub(crate) enum IKType {
    I = 1,
    K = 2,
}

pub(crate) fn i_pow<T: BesselFloat>(n: usize) -> Complex<T> {
    match n % 4 {
        0 => Complex::new(T::one(), T::zero()),
        1 => Complex::new(T::zero(), T::one()),
        2 => Complex::new(-T::one(), T::zero()),
        3 => Complex::new(T::zero(), -T::one()),
        _ => unreachable!(),
    }
}

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
    /// Returns the signum of the rotation direction as an `f64`.
    ///
    /// **WARNING:** In Rust, `0.0_f64.signum()` evaluates to `1.0`. So calling this on
    /// `RotationDirection::None` will return `1.0`. Do not "optimize" this to return `0.0`
    /// for `None`! This behavior perfectly mimics the legacy Fortran `SIGN(1.0, 0.0)` which
    /// also evaluates to `1.0` by transferring the positive sign of zero. The codebase
    /// relies on this Fortran quirk.
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
