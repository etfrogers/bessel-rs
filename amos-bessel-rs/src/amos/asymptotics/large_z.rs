use num::{
    Integer,
    complex::{Complex, ComplexFloat},
};

use crate::{
    BesselError, BesselFloat, Scaling,
    amos::{
        MachineConsts,
        utils::{RECIP_TWO_PI, two_over_z_safe},
    },
    types::BesselResult,
};

/// Computes the sequence $[I_\nu(z), \dots, I_{\nu+n-1}(z)]$ for $\operatorname{Re}(z) \ge 0$
/// using Hankel's asymptotic expansion for large argument $|z| > \max(R_L, \nu^2 / 2)$.
///
/// ### Mathematical Expansion
/// For large $|z|$, $I_\nu(z)$ is given by Hankel's two-series expansion (Abramowitz & Stegun 9.7.1, NIST DLMF 10.40.1):
/// $$I_\nu(z) \sim \frac{e^z}{\sqrt{2\pi z}} \sum_{m=0}^\infty \frac{(-1)^m (\nu, m)}{(8z)^m} + \frac{e^{-z \pm i\pi(\nu + 1/2)}}{\sqrt{2\pi z}} \sum_{m=0}^\infty \frac{(\nu, m)}{(8z)^m}$$
/// where $(\nu, m) / 8^m$ is Hankel's symbol.
///
/// The two highest orders in the sequence ($\nu + n - 1$ and $\nu + n - 2$) are evaluated via the series,
/// and the remaining elements are obtained by backward recurrence:
/// $$I_k(z) = \frac{2(k + 1 + \nu)}{z} I_{k+1}(z) + I_{k+2}(z)$$
///
/// Originally Amos routine `ZASYI`.
///
/// # Arguments
/// * `z` - Complex argument $\operatorname{Re}(z) \ge 0$.
/// * `order` - Starting order $\nu \ge 0$.
/// * `scaling` - `Scaling::Unscaled` ($I_\nu(z)$) or `Scaling::Scaled` ($e^{-|\operatorname{Re}(z)|} I_\nu(z)$).
/// * `n` - Number of terms in the sequence.
pub fn i_asymptotic<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    n: usize,
) -> BesselResult<T, usize> {
    let mc: &MachineConsts<T> = T::MACHINE_CONSTANTS;
    let mut y = T::c_zeros(n);
    let abs_z = z.abs();
    let recip_abs_z = T::one() / abs_z;

    // Overflow / Exponent Range Test
    let exponent_arg = match scaling {
        Scaling::Unscaled => z,
        Scaling::Scaled => Complex::<T>::new(T::ZERO, z.im),
    };
    if exponent_arg.re.abs() > mc.exponent_limit {
        return Err(BesselError::Overflow);
    }
    let scaled_calculations = exponent_arg.re.abs() > mc.approximation_limit;

    let mut prefactor = (T::from_f64(RECIP_TWO_PI) * z.conj() * recip_abs_z.powi(2)).sqrt();
    if !scaled_calculations {
        prefactor *= exponent_arg.exp();
    }
    let eight_z = T::from_f64(8.0) * z;
    let abs_eight_z = T::from_f64(8.0) * abs_z;

    // When z is imaginary, the error test must be made relative to the
    // first reciprocal power since this is the leading term of the
    // expansion for the imaginary part.
    let rel_tol_scale = mc.abs_error_tolerance / abs_eight_z;
    let max_iterations = (mc.asymptotic_z_limit * T::TWO).to_usize().unwrap() + 2;
    let mut phase_factor = if z.im == T::ZERO {
        T::C_ZERO
    } else {
        // Compute the Stokes phase factor exp(i*pi*(0.5 + order + k) * sgn(Im(z)))
        // using fract(order) to prevent loss of precision when order or n is large.
        let arg = order.fract() * T::PI();
        let phase_re = -arg.sin();
        let mut phase_im = arg.cos();
        if z.im < T::ZERO {
            phase_im = -phase_im;
        };
        let phase_factor = Complex::<T>::new(phase_re, phase_im);
        if (order.to_usize().unwrap() + n).is_even() {
            -phase_factor
        } else {
            phase_factor
        }
    };

    for (k, elem) in y.iter_mut().enumerate().rev().take(2.min(n)) {
        let (mut sum_dominant, sum_subdominant) = {
            // this block is just to contain the large number of mutable variables in a small space
            let modified_order = order + T::from_usize(k);
            let four_order_sqr = (T::TWO * modified_order).powf(T::TWO);
            let atol = rel_tol_scale * (four_order_sqr - T::one()).abs();
            let mut sign = T::one();
            let mut sum_alternating = T::C_ONE;
            let mut sum_direct = T::C_ONE;
            let mut term = T::C_ONE;
            let mut term_magnitude = T::one();
            let mut converged = false;
            for i in 0..max_iterations {
                let odd = T::from_usize(2 * i + 1);
                let step = T::from_usize(i + 1);
                let numerator_factor = four_order_sqr - odd.powi(2); // 4ν² - (2i + 1)²
                let denominator_factor = step * eight_z;
                term *= numerator_factor / denominator_factor;
                sum_direct += term;
                sign = -sign;
                sum_alternating += term * sign;
                term_magnitude *= numerator_factor.abs() / (step * abs_eight_z);
                if term_magnitude <= atol {
                    converged = true;
                    break;
                }
            }
            if !converged {
                return Err(BesselError::DidNotConverge);
            }
            (sum_alternating, sum_direct)
        };
        if z.re * T::TWO < mc.exponent_limit {
            sum_dominant += (-z * T::TWO).exp() * phase_factor * sum_subdominant;
        }
        phase_factor = -phase_factor;
        *elem = sum_dominant * prefactor;
    }
    if n > 2 {
        let two_over_z = two_over_z_safe(z);
        // recur downward from the last two elements
        for k in (0..n - 2).rev() {
            y[k] = (two_over_z * y[k + 1]) * (T::from_usize(k + 1) + order) + y[k + 2];
        }
    }
    if scaled_calculations {
        let exp_cz = exponent_arg.exp();
        for yi in y.iter_mut() {
            *yi *= exp_cz;
        }
    }
    Ok((y, 0))
}
