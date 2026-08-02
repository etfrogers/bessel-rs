use num::complex::{Complex, ComplexFloat};

use crate::{BesselError, BesselFloat};

/// $1/(2\pi) \approx 0.159154943...$, used in asymptotic prefactors $\sqrt{1/(2\pi z)}$.
pub const RECIP_TWO_PI: f64 = 0.159_154_943_091_895_35;
/// $\sqrt{3} \approx 1.7320508...$, threshold where $|\text{Im}(z)| > \sqrt{3}|\text{Re}(z)| \iff |\arg(z)| > 60^\circ$.
pub const RT_THREE: f64 = 1.73205080757;
/// $\ln(2\sqrt{\pi}) = \ln(2) + \frac{1}{2}\ln(\pi) \approx 1.2655121...$, the Airy function
/// asymptotic normalization constant (also equals $\text{Re}(\ln \Gamma(-0.5))$).
pub const AIC: f64 = 1.265_512_123_484_645_4;

/// Computes $2/z = 2(\bar{z}/|z|) / |z|$ in a numerically safe manner.
///
/// Normalizing by $|z|$ in two stages computes a unit vector $\bar{z}/|z|$ first, avoiding
/// intermediate overflow or underflow of $|z|^2$ when $|z|$ is near machine boundaries.
pub(crate) fn two_over_z_safe<T: BesselFloat>(z: Complex<T>) -> Complex<T> {
    let r_abs_z = T::one() / z.abs();
    let intermediate = z.conj() * r_abs_z;
    (intermediate + intermediate) * r_abs_z
}

/// Tests whether $|\text{Im}(z)| > \sqrt{3}|\text{Re}(z)|$ (i.e. angle is within $30^\circ$ of the imaginary axis).
pub(crate) fn imaginary_dominant<T: BesselFloat>(z: Complex<T>) -> bool {
    z.im.abs() > z.re.abs() * T::from_f64(RT_THREE)
}

/// `y` enters as a scaled quantity whose magnitude is greater than
/// (all names below are from [crate::amos::machine::MachineConsts])
/// `(-approximation_limit).exp() = absolute_approximation_limit
///  = 2.0 * f64::MIN_POSITIVE/abs_error_tolerance`.
/// The test is made to see
/// if the magnitude of the real or imaginary part would underflow
/// when y is scaled (by tol) to its proper value. `y` is accepted
/// if the underflow is at least one precision below the magnitude
/// of the largest component; otherwise the phase angle does not have
/// absolute accuracy and an underflow is assumed
pub(crate) fn will_underflow<T: BesselFloat>(y: Complex<T>) -> bool {
    let re_abs = y.re.abs();
    let im_abs = y.im.abs();
    let min_abs_component = re_abs.min(im_abs);
    if min_abs_component > T::MACHINE_CONSTANTS.absolute_approximation_limit {
        false
    } else {
        let max_abs_component = re_abs.max(im_abs);
        // Accepted if the smaller component is within tolerance of the larger: min / max > tol
        max_abs_component < min_abs_component / T::MACHINE_CONSTANTS.abs_error_tolerance
    }
}

/// Checks whether $|z|$ or $\nu$ exceed numerical precision limits.
///
/// - Returns `Err(LossOfSignificance)` if arguments exceed `upper_size_limit` (total precision loss).
/// - Returns `Ok(true)` if arguments exceed $\sqrt{\text{upper\_size\_limit}}$ (partial loss of significance,
///   exponential scaling recommended).
pub fn is_significance_lost<T: BesselFloat>(
    abs_z: T,
    modified_order: T,
    modify_threshold: bool,
) -> Result<bool, BesselError<T>> {
    let f64_precision_limit = T::half() / T::MACHINE_CONSTANTS.abs_error_tolerance;
    // TODO the below is limited to i32: could push to 64 later, but would change compare to fortran
    let integer_size_limit = T::from_f64((i32::MAX as f64) * 0.5);
    let mut upper_size_limit = f64_precision_limit.min(integer_size_limit);
    if modify_threshold {
        upper_size_limit = upper_size_limit.powf(T::TWO_THIRDS);
    }
    if abs_z > upper_size_limit || modified_order > upper_size_limit {
        return Err(BesselError::LossOfSignificance);
    }
    let scaling_limit = upper_size_limit.sqrt();
    Ok((abs_z > scaling_limit) || (modified_order > scaling_limit))
}

/// Validates Bessel inputs: checks that $z \neq 0$ (if requested), $\text{order} \geq 0$, and sequence length $N \geq 1$.
pub(crate) fn sanitise_inputs<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    n: usize,
    check_z_zero: bool,
) -> Result<(), BesselError<T>> {
    let mut err = None;
    if check_z_zero && z.re == T::zero() && z.im == T::zero() {
        err = Some("z must not be zero");
    }
    if order < T::zero() {
        err = Some("order must be positive");
    };
    if n < 1 {
        err = Some("N must be >= 1");
    };
    if let Some(details) = err {
        Err(BesselError::InvalidInput {
            details: details.to_owned(),
        })
    } else {
        Ok(())
    }
}
