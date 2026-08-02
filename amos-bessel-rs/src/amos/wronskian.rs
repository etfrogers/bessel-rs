use super::Scaling;
use crate::{
    BesselError, BesselFloat,
    amos::{
        IKType, limits::check_underflow_uniform_asymp_params, recurrence::i_ratios,
        right_half_plane::k_right_half_plane,
    },
};

use num::{Complex, complex::ComplexFloat};

/// Computes the $I$ Bessel sequence for $\text{Re}(z) \ge 0$ by
/// normalizing the ratios from [i_ratios] using the Wronskian identity with $K_\nu$ and $K_{\nu+1}$.
///
/// Originally ZWRSK
pub(crate) fn i_wronskian<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    n: usize,
    y: &mut [Complex<T>],
) -> Result<usize, BesselError<T>> {
    match check_underflow_uniform_asymp_params(z, order, scaling, IKType::K, 2, &mut [T::C_ONE; 2])
    {
        Ok(n_underflow) => {
            if n_underflow > 0 {
                return Err(BesselError::Overflow);
            }
        }
        Err(_) => {
            y.fill(T::C_ZERO);
            return Ok(n);
        }
    }

    // 1. Compute K_nu and K_{nu+1} to serve as Wronskian anchors
    let (k_values, _) = k_right_half_plane(z, order, scaling, 2)?;
    // 2. Compute backward recurrence ratios r_{nu+j} = I_{nu+j+1} / I_{nu+j}
    let y_ratios = i_ratios(z, order, n);

    // Initial phase factor for scaled computation (e^{i * Im(z)})
    let mut current_i = if scaling == Scaling::Scaled {
        Complex::<T>::cis(z.im)
    } else {
        T::C_ONE
    };

    // On low-exponent machines, K values can be close to under/overflow limits.
    // Scale K_nu and K_{nu+1} to keep intermediate products well on scale.
    let abs_k_nu_plus_1 = k_values[1].abs();
    let k_scale_factor = if abs_k_nu_plus_1 <= T::MACHINE_CONSTANTS.absolute_approximation_limit {
        T::one() / T::MACHINE_CONSTANTS.abs_error_tolerance
    } else if abs_k_nu_plus_1 >= T::one() / T::MACHINE_CONSTANTS.absolute_approximation_limit {
        T::MACHINE_CONSTANTS.abs_error_tolerance
    } else {
        T::one()
    };
    let scaled_k_nu = k_values[0] * k_scale_factor;
    let scaled_k_nu_plus_1 = k_values[1] * k_scale_factor;

    // Evaluate Wronskian denominator: denom = z * (r_nu * K_nu + K_{nu+1}) = 1 / I_nu
    // Performing division as (current_i / |denom|) * (conj(denom) / |denom|) avoids
    // squaring |denom| which could overflow/underflow prematurely.
    let mut wronskian_denom = z * (y_ratios[0] * scaled_k_nu + scaled_k_nu_plus_1);
    let abs_denom = wronskian_denom.abs();
    wronskian_denom = wronskian_denom.conj() / abs_denom;
    current_i = (current_i / abs_denom) * wronskian_denom;

    // Multiply by k_scale_factor to restore true scale
    y[0] = current_i * k_scale_factor;

    // Step forward: I_{nu+j} = I_{nu+j-1} * r_{nu+j-1}
    for i in 1..n {
        current_i *= y_ratios[i - 1];
        y[i] = current_i * k_scale_factor;
    }
    Ok(0)
}
