use num::Complex;
use num::complex::ComplexFloat;

use crate::amos::{
    ComplexExt, IKType,
    asymptotics::{AiryGeometry, DebyeGeometry},
    utils::{AIC, imaginary_dominant, will_underflow},
};
use crate::{BesselError, BesselFloat, Scaling};

#[allow(unused_imports)]
use super::machine::MachineConsts;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OverflowState {
    Over { was_refined: bool },
    NearOver,
    Under { was_refined: bool },
    NearUnder,
    None,
}

impl OverflowState {
    pub fn check<T: BesselFloat>(
        rs1: T,
        phi: Complex<T>,
        extra_refinement: T,
        mc: &MachineConsts<T>,
    ) -> Self {
        //-----------------------------------------------------------------------
        //     TEST FOR UNDERFLOW AND OVERFLOW
        //-----------------------------------------------------------------------
        if rs1.abs() > mc.exponent_limit {
            return if rs1 > T::ZERO {
                Self::Over { was_refined: false }
            } else {
                Self::Under { was_refined: false }
            };
        }
        if rs1.abs() < mc.approximation_limit {
            return Self::None;
        }
        //-----------------------------------------------------------------------
        //     REFINE  TEST AND SCALE
        //-----------------------------------------------------------------------
        let refined_rs1 = rs1 + phi.abs().ln() + extra_refinement;
        if refined_rs1.abs() > mc.exponent_limit {
            return if refined_rs1 > T::ZERO {
                Self::Over { was_refined: true }
            } else {
                Self::Under { was_refined: true }
            };
        }
        if refined_rs1 > T::ZERO {
            Self::NearOver
        } else {
            Self::NearUnder
        }
    }

    pub fn increment(&mut self) {
        match self {
            OverflowState::Over { .. } | OverflowState::Under { .. } => {
                panic!("Overflow and underflow are not valid for incrementation")
            }
            OverflowState::NearOver => {
                panic!("NearOver is the largest possible overflow condition")
            }
            OverflowState::NearUnder => *self = Self::None,
            OverflowState::None => *self = Self::NearOver,
        }
    }

    // Originally mc.scaling_factors[overflow_state]
    pub fn scaling_factor<T: BesselFloat>(&self, mc: &MachineConsts<T>) -> T {
        match self {
            OverflowState::NearUnder => mc.rtol,
            OverflowState::None => T::ONE,
            OverflowState::NearOver => mc.abs_error_tolerance,
            _ => panic!("Cannot get scaling factor for fatal overflow/underflow"),
        }
    }

    // Originally mc.reciprocal_scaling_factors[overflow_state]
    pub fn reciprocal_scaling_factor<T: BesselFloat>(&self, mc: &MachineConsts<T>) -> T {
        match self {
            OverflowState::NearUnder => mc.abs_error_tolerance,
            OverflowState::None => T::ONE,
            OverflowState::NearOver => mc.rtol,
            _ => panic!("Cannot get reciprocal scaling factor for fatal overflow/underflow"),
        }
    }

    // Originally mc.overflow_boundary[overflow_state]
    pub fn boundary<T: BesselFloat>(&self, mc: &MachineConsts<T>) -> T {
        match self {
            OverflowState::NearUnder => mc.absolute_approximation_limit,
            OverflowState::None => T::ONE / mc.absolute_approximation_limit,
            OverflowState::NearOver => T::max_value() / T::TWO,
            _ => panic!("Cannot get boundary for fatal overflow/underflow"),
        }
    }

    /// Adjusts the scaling factors during a recurrence step if the magnitude of the
    /// newly computed term exceeds the boundary for the current overflow state.
    ///
    /// The values `s1` and `s2` are the recurrence state variables.
    /// `unscaled_s2` is the newly computed term before scaling.
    ///
    /// Note: `boundary` and `recip_scaling_factor` are explicitly passed in as mutable
    /// references rather than evaluated internally via `self.boundary()` and
    /// `self.reciprocal_scaling_factor()`. This acts as a manual loop-invariant code
    /// motion (hoisting), preventing the compiler from executing the `match` branches
    /// inside those functions on every iteration of the tight innermost recurrence loops.
    pub fn scale_recurrence<T: BesselFloat>(
        &mut self,
        s1: &mut Complex<T>,
        s2: &mut Complex<T>,
        unscaled_s2: Complex<T>,
        boundary: &mut T,
        recip_scaling_factor: &mut T,
        mc: &MachineConsts<T>,
    ) {
        if *self != OverflowState::NearOver && unscaled_s2.linf_norm() > *boundary {
            self.increment();
            *boundary = self.boundary::<T>(mc);
            *s1 *= *recip_scaling_factor;
            *s2 = unscaled_s2;
            let scaling = self.scaling_factor::<T>(mc);
            *s1 *= scaling;
            *s2 *= scaling;
            *recip_scaling_factor = self.reciprocal_scaling_factor::<T>(mc);
        }
    }
}

/// check_underflow_uniform_asymp_params computes the leading terms of the uniform asymptotic
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
///     The first `(n_to_test - n_underflow)` values must be set by another routine
///
/// `ik_type == IKType::K && n_underflow == n_to_test` means all y values were set to zero
///
/// `ik_type = IKType::K && 0 < n_underflow < n_to_test` not considered (it is not
///     a possible return value). y must be set by
///     another routine
///
/// Originally ZUIOK
pub(crate) fn check_underflow_uniform_asymp_params<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    ik_type: IKType,
    n_to_test: usize,
    y: &mut [Complex<T>],
    mc: &MachineConsts<T>,
) -> Result<usize, BesselError<T>> {
    let mut n_underflow = 0;
    let zr = if z.re < T::ZERO { -z } else { z };
    let zn = if z.im <= T::ZERO {
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
            let AiryGeometry {
                phi,
                arg,
                zeta1,
                zeta2,
                ..
            } = AiryGeometry::compute(zn, modified_order);
            (-zeta1 + zeta2, phi, arg, arg.abs())
        } else {
            let DebyeGeometry {
                phi_i,
                phi_k,
                zeta1,
                zeta2,
                ..
            } = DebyeGeometry::compute(zr, modified_order);
            let phi = match ik_type {
                IKType::I => phi_i,
                IKType::K => phi_k,
            };
            (-zeta1 + zeta2, phi, T::C_ZERO, T::ZERO)
        };
        if scaling == Scaling::Scaled {
            cz -= zb;
        }
        let refinement = if imaginary_dominant {
            -abs_arg.ln() * T::from_f64(0.25) - T::from_f64(AIC)
        } else {
            T::ZERO
        };
        (cz, phi, arg, refinement)
    };

    // First checks the last element
    let modified_order = match ik_type {
        IKType::K => {
            let float_n = T::from_usize(n_to_test);
            let modified_order = order + float_n - T::ONE;
            modified_order.max(float_n)
        }
        IKType::I => order.max(T::ONE),
    };

    let (mut cz, phi, arg, extra_refinement) = get_parameters(modified_order);
    if ik_type == IKType::K {
        cz = -cz;
    }
    //-----------------------------------------------------------------------
    //     OVERFLOW TEST
    //-----------------------------------------------------------------------
    match OverflowState::check(cz.re, phi, extra_refinement, mc) {
        OverflowState::Over { .. } => return Err(BesselError::Overflow),
        OverflowState::Under { was_refined } => {
            if !was_refined {
                y[0..n_to_test].fill(T::C_ZERO);
            }
            return Ok(n_to_test);
        }
        OverflowState::NearUnder => {
            cz += phi.ln();
            if imaginary_dominant {
                cz -= T::from_f64(0.25) * arg.ln() + T::from_f64(AIC);
            }
            cz = cz.exp() / mc.abs_error_tolerance;
            if will_underflow(cz, mc) {
                y[0..n_to_test].fill(T::C_ZERO);
                return Ok(n_to_test);
            }
        }
        OverflowState::None | OverflowState::NearOver => (),
    }
    // On K type, we only check the max n_to_test value, as per function documentation
    if ik_type == IKType::K || n_to_test == 1 {
        return Ok(n_underflow);
    }
    //-----------------------------------------------------------------------
    //     SET UNDERFLOWS ON I SEQUENCE
    //-----------------------------------------------------------------------
    // Note n_to_test is NOT y.len() in this case.
    for (i, yi) in y.iter_mut().enumerate().take(n_to_test).rev() {
        let current_order = order + T::from_usize(i);
        let (mut cz, phi, _arg, extra_refinement) = get_parameters(current_order);
        // Match below says that first time we get here and no underflow is found, we immediately return
        match OverflowState::check(cz.re, phi, extra_refinement, mc) {
            OverflowState::Under { was_refined } => {
                if was_refined {
                    // Now do a similar overflow check, but on complex values, rather
                    // than the absolute values used in find_overflow
                    cz += phi.ln();
                    if imaginary_dominant {
                        cz -= arg.ln() * T::from_f64(0.25) + T::from_f64(AIC)
                    }
                    cz = cz.exp() / mc.abs_error_tolerance;
                    if !will_underflow(cz, mc) {
                        return Ok(n_underflow);
                    }
                }
            }
            OverflowState::NearUnder => (),
            OverflowState::None | OverflowState::NearOver | OverflowState::Over { .. } => {
                return Ok(n_underflow);
            }
        }
        *yi = T::C_ZERO;
        n_underflow += 1;
    }
    Ok(n_underflow)
}

/// underflow_add_i_k tests for a possible underflow resulting from the
/// addition of the i and k functions in the analytic
/// continuation formula where s_k = k function and s_k = i function.
/// In unscaled calculations the i and k functions are different orders of
/// magnitude, but for scaled calculations they can be of the same order
/// of magnitude and the maximum must be at least one
/// precision above the underflow limit.
///
/// Returns `true` if underflow is found, and `false` otherwise, but can also modify
/// the value of s1, s2, and n_underflow.
///
/// If s_k is large enough to be set without underflowing, it is
/// set to `(s_k.ln() - 2.0 * zr).exp()`. That is, `2*zr` is subtracted
/// from the exponent.
///
/// If underflow is found, the function sets s1 and s2 to zero.
///
/// Originally ZS1S2
pub(crate) fn underflow_add_i_k<T: BesselFloat>(
    zr: Complex<T>,
    s_k: &mut Complex<T>,
    s_i: &mut Complex<T>,
    mc: &MachineConsts<T>,
) -> bool {
    let abs_s_k = s_k.abs();
    if abs_s_k != T::ZERO {
        let test = (-T::TWO * zr.re) + abs_s_k.ln();
        if test >= (-mc.approximation_limit) {
            *s_k = (s_k.ln() - zr * T::TWO).exp();
        } else {
            *s_k = T::C_ZERO;
        }
    }
    if s_k.abs().max(s_i.abs()) > mc.absolute_approximation_limit {
        false
    } else {
        *s_k = T::C_ZERO;
        *s_i = T::C_ZERO;
        true
    }
}
