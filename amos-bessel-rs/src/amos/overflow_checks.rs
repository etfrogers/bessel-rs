use std::ops::Index;

use num::Complex;
use num::complex::ComplexFloat;

use crate::BesselError;
use crate::amos::asymptotics::{hj_uniform_asymp_params, ik_uniform_asymp_params};
use crate::amos::utils::{AIC, imaginary_dominant};
use crate::types::{BesselError::*, BesselFloat};

use super::{IKType, Scaling, utils::will_underflow};

#[allow(unused_imports)]
use super::machine::MachineConsts;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Overflow {
    Over(bool),
    NearOver,
    Under(bool),
    NearUnder,
    None,
}

impl Overflow {
    pub fn find_overflow<T: BesselFloat>(rs1: T, phi: Complex<T>, extra_refinement: T) -> Self {
        //-----------------------------------------------------------------------
        //     TEST FOR UNDERFLOW AND OVERFLOW
        //-----------------------------------------------------------------------
        if rs1.abs() > T::MACHINE_CONSTANTS.exponent_limit {
            return if rs1 > T::zero() {
                Self::Over(false)
            } else {
                Self::Under(false)
            };
        }
        if rs1.abs() < T::MACHINE_CONSTANTS.approximation_limit {
            return Self::None;
        }
        //-----------------------------------------------------------------------
        //     REFINE  TEST AND SCALE
        //-----------------------------------------------------------------------
        let refined_rs1 = rs1 + phi.abs().ln() + extra_refinement;
        if refined_rs1.abs() > T::MACHINE_CONSTANTS.exponent_limit {
            return if refined_rs1 > T::zero() {
                Self::Over(true)
            } else {
                Self::Under(true)
            };
        }
        if refined_rs1 > T::zero() {
            Self::NearOver
        } else {
            Self::NearUnder
        }
    }

    pub fn increment(&mut self) {
        match self {
            Overflow::Over(_) | Overflow::Under(_) => {
                panic!("Overflow and underflow are not valid for incrementation")
            }
            Overflow::NearOver => panic!("NearOver is the largest possible overflow condition"),
            Overflow::NearUnder => *self = Self::None,
            Overflow::None => *self = Self::NearOver,
        }
    }
}

impl<T: BesselFloat> Index<Overflow> for [T] {
    type Output = T;

    fn index(&self, index: Overflow) -> &Self::Output {
        match index {
            Overflow::Over(_) | Overflow::Under(_) => {
                panic!("Overflow and underflow are not valid indices")
            }
            Overflow::NearOver => &self[2],
            Overflow::NearUnder => &self[0],
            Overflow::None => &self[1],
        }
    }
}

/// zuoik computes the leading terms of the uniform asymptotic
/// expansions for the I and K functions and compares them
/// (in logarithmic form) to [MachineConsts::approximation_limit]
/// and [MachineConsts::exponent_limit] for over and underflow
///
/// If the magnitude, based on the leading
/// exponential, is less than alim or greater than -approximation_limit, then
/// the result is on scale. If not, then a refined test using other
/// multipliers (in logarithmic form) is made based on exponent_limit.
///
/// It sets the appropriate input y values to zero on underflow.
///
/// Returns Err(Overflow) on overflow
///
/// ik_type specifies the type of squence to test ([I](IKType::I) or [K](IKType::K))
///
/// Below the returned value is designated n_underflow.
///
/// n_underflow = 0 means the last member of the sequence is on scale
///
/// `ik_type == IKType::I && n_underflow > 0` means the last
///     `n_underflow` y values were set to zero.
///     The first `(n - n_underflow)` values must be set by another routine
///
/// `ik_type == IKType::K && n_underflow == n` means all y values were set to zero
///
/// `ik_type = IKType::K && 0 < n_underflow < n` not considered (it is not
///     a possible return value). y must be set by
///     another routine
///
/// Originally ZUIOK
pub fn check_underflow_uniform_asymp_params<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    ik_type: IKType,
    n_to_test: usize,
    y: &mut [Complex<T>],
) -> Result<usize, BesselError<T>> {
    let mut n_underflow = 0;
    let zr = if z.re < T::zero() { -z } else { z };
    let zn = if z.im <= T::zero() {
        -T::I * zr.conj()
    } else {
        -T::I * zr
    };
    let zb = zr;
    let imaginary_dominant = imaginary_dominant(z);
    //-----------------------------------------------------------------------
    //     only the magnitude of arg and phi are needed along with the
    //     real parts of zeta1, zeta2 and zb. no attempt is made to get
    //     the sign of the imaginary part correct.
    //-----------------------------------------------------------------------

    // This piece of code is used in two places, where there is essentially a switch on which function to use,
    // based on whether z is imaginary dominant or real dominant
    let get_parameters = |modified_order: T| {
        let (mut cz, phi, arg, abs_arg) = if imaginary_dominant {
            let (phi, arg, zeta1, zeta2, _, _) = hj_uniform_asymp_params(zn, modified_order, true);
            (-zeta1 + zeta2, phi, arg, arg.abs())
        } else {
            let (phi, zeta1, zeta2, _) = ik_uniform_asymp_params(zr, modified_order, ik_type, true);
            (-zeta1 + zeta2, phi, T::C_ZERO, T::zero())
        };
        if scaling == Scaling::Scaled {
            cz -= zb;
        }
        let refinement = if imaginary_dominant {
            -abs_arg.ln() * T::from_f64(0.25) - T::from_f64(AIC)
        } else {
            T::zero()
        };
        (cz, phi, arg, refinement)
    };

    // First checks the last element
    let modified_order = match ik_type {
        IKType::K => {
            let float_n = T::from_usize(n_to_test);
            let modified_order = order + float_n - T::one();
            modified_order.max(float_n)
        }
        IKType::I => order.max(T::one()),
    };

    let (mut cz, phi, arg, extra_refinement) = get_parameters(modified_order);
    if ik_type == IKType::K {
        cz = -cz;
    }
    //-----------------------------------------------------------------------
    //     OVERFLOW TEST
    //-----------------------------------------------------------------------
    match Overflow::find_overflow(cz.re, phi, extra_refinement) {
        Overflow::Over(_) => return Err(Overflow),
        Overflow::Under(was_refined) => {
            if !was_refined {
                y[0..n_to_test].iter_mut().for_each(|v| *v = T::C_ZERO);
            }
            return Ok(n_to_test);
        }
        Overflow::NearUnder => {
            cz += phi.ln();
            if imaginary_dominant {
                cz -= T::from_f64(0.25) * arg.ln() + T::from_f64(AIC);
            }
            cz = cz.exp() / T::MACHINE_CONSTANTS.abs_error_tolerance;
            if will_underflow(
                cz,
                T::MACHINE_CONSTANTS.absolute_approximation_limit,
                T::MACHINE_CONSTANTS.abs_error_tolerance,
            ) {
                y[0..n_to_test].iter_mut().for_each(|v| *v = T::C_ZERO);
                return Ok(n_to_test);
            }
        }
        Overflow::None | Overflow::NearOver => (),
    }
    // On K type, we only check the max n value, as per function documentation
    if ik_type == IKType::K || n_to_test == 1 {
        return Ok(n_underflow);
    }
    //-----------------------------------------------------------------------
    //     SET UNDERFLOWS ON I SEQUENCE
    //-----------------------------------------------------------------------
    // Note n_to_test is NOT y.len() in this case.
    for (i, yi) in y.iter_mut().enumerate().take(n_to_test).rev() {
        let modified_order = order + T::from_usize(i);
        let (mut cz, phi, _arg, extra_refinement) = get_parameters(modified_order);
        // Match below says that first time we get here and no underflow is found, we immediately return
        match Overflow::find_overflow(cz.re, phi, extra_refinement) {
            Overflow::Under(was_refined) => {
                if was_refined {
                    // Now do a similar overflow check, but on complex values, rather
                    // than the absolute values used in find_overflow
                    cz += phi.ln();
                    if imaginary_dominant {
                        cz -= arg.ln() * T::from_f64(0.25) + T::from_f64(AIC)
                    }
                    cz = cz.exp() / T::MACHINE_CONSTANTS.abs_error_tolerance;
                    if !will_underflow(
                        cz,
                        T::MACHINE_CONSTANTS.absolute_approximation_limit,
                        T::MACHINE_CONSTANTS.abs_error_tolerance,
                    ) {
                        return Ok(n_underflow);
                    }
                }
            }
            Overflow::NearUnder => (),
            Overflow::None | Overflow::NearOver | Overflow::Over(_) => return Ok(n_underflow),
        }
        *yi = T::C_ZERO;
        n_underflow += 1;
    }
    Ok(n_underflow)
}

/// underflow_add_i_k tests for a possible underflow resulting from the
/// addition of the i and k functions in the analytic con-
/// tinuation formula where s1=k function and s2=i function.
/// on kode=1 the i and k functions are different orders of
/// magnitude, but for kode=2 they can be of the same order
/// of magnitude and the maximum must be at least one
/// precision above the underflow limit.
///
/// Returns 1 if underflow is found, and 0 otherwise, but can also modify
/// the value of s1, s2, and n_underflow.
///
/// If s1 is large enough to be set without underflowing, it is
/// set to (s1.ln() - 2.0 * zr).exp(); That is, 2*zr is subtracted
/// from the exponent.
///
/// If underflow is found, the function sets s1 and s2 to zero.
///
/// Originally ZS1S2
pub fn underflow_add_i_k<T: BesselFloat>(
    zr: Complex<T>,
    s_k: &mut Complex<T>,
    s_i: &mut Complex<T>,
    n_good: &mut isize,
) -> usize {
    let mut abs_s1 = s_k.abs();
    let abs_s2 = s_i.abs();
    if (s_k.re != T::zero() || s_k.im != T::zero()) && (abs_s1 != T::zero()) {
        let test = (-T::two() * zr.re) + abs_s1.ln();
        let s1d = *s_k;
        *s_k = T::C_ZERO;
        abs_s1 = T::zero();
        if test >= (-T::MACHINE_CONSTANTS.approximation_limit) {
            *s_k = (s1d.ln() - zr * T::two()).exp();
            abs_s1 = s_k.abs();
            *n_good += 1;
        }
    }
    if abs_s1.max(abs_s2) > T::MACHINE_CONSTANTS.absolute_approximation_limit {
        0
    } else {
        *s_k = T::C_ZERO;
        *s_i = T::C_ZERO;
        *n_good = 0;
        1
    }
}
