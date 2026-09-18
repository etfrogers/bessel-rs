use num::Complex;
use num_integer::binomial;

use crate::{
    BesselError, BesselFloat, BesselInput, HankelKind, Scaling,
    amos::{
        complex_bessel_i, complex_bessel_j, complex_bessel_k, complex_bessel_y, complex_hankel1,
        complex_hankel2,
    },
    reflections::integer_sign,
    types::BesselValues,
};

#[allow(type_alias_bounds)]
type BesselSig<T: BesselFloat = f64> =
    fn(Complex<T>, T, Scaling, usize) -> Result<BesselValues<T>, BesselError<T>>;

/// Computes the first derivative of the Bessel function of the first kind $J_\nu'(z)$ with respect to $z$.
///
/// # Arguments
/// * `order` - The order $\nu$ of the Bessel function (can be integer or non-integer, positive or negative).
/// * `z` - The complex or real argument.
///
/// # Examples
/// ```
/// use amos_bessel_rs::derivatives::bessel_j_p;
///
/// let dj: f64 = bessel_j_p(0.0, 1.0).unwrap();
/// assert!((dj - (-0.4400505857)).abs() < 1e-6);
/// ```
pub fn bessel_j_p<FT: BesselFloat, ZT: BesselInput<FT>, OT: Into<FT>>(
    order: OT,
    z: ZT,
) -> Result<ZT, BesselError<FT>> {
    bessel_j_derivative(order, z, 1, Scaling::Unscaled)
}

/// Computes the $k$-th derivative of the Bessel function of the first kind $\left(\frac{d}{dz}\right)^k J_\nu(z)$.
///
/// Evaluated via the exact binomial identity from DLMF 10.6.7:
/// $$\left(\frac{d}{dz}\right)^k J_\nu(z) = \frac{1}{2^k} \sum_{n=0}^k (-1)^n \binom{k}{n} J_{\nu - k + 2n}(z)$$
///
/// When `scaling` is [`Scaling::Scaled`], the result is scaled by $e^{-|\mathrm{Im}(z)|}$.
///
/// # Arguments
/// * `order` - The order $\nu$ of the Bessel function.
/// * `z` - The complex or real argument.
/// * `derivative_order` - The order of the derivative $k \ge 0$.
/// * `scaling` - Whether to compute the unscaled or exponentially scaled derivative.
///
/// # Examples
/// ```
/// use amos_bessel_rs::{Scaling, derivatives::bessel_j_derivative};
///
/// // Second derivative of J_0 at z = 1.0:
/// let d2j: f64 = bessel_j_derivative(0.0, 1.0, 2, Scaling::Unscaled).unwrap();
/// assert!((d2j - (-0.3251471008)).abs() < 1e-6);
/// ```
pub fn bessel_j_derivative<FT: BesselFloat, ZT: BesselInput<FT>, OT: Into<FT>>(
    order: OT,
    z: ZT,
    derivative_order: u32,
    scaling: Scaling,
) -> Result<ZT, BesselError<FT>> {
    derivative_internal(
        complex_bessel_j,
        order.into(),
        z.into(),
        derivative_order,
        scaling,
        SignType::Cylinder,
    )
    .and_then(ZT::back_from)
}

/// Computes the first derivative of the Bessel function of the second kind $Y_\nu'(z)$ with respect to $z$.
///
/// # Arguments
/// * `order` - The order $\nu$ of the Bessel function.
/// * `z` - The complex or real argument.
///
/// # Examples
/// ```
/// use amos_bessel_rs::derivatives::bessel_y_p;
///
/// let dy: f64 = bessel_y_p(0.0, 1.0).unwrap();
/// assert!((dy - 0.7812128213).abs() < 1e-6);
/// ```
pub fn bessel_y_p<FT: BesselFloat, ZT: BesselInput<FT>, OT: Into<FT>>(
    order: OT,
    z: ZT,
) -> Result<ZT, BesselError<FT>> {
    bessel_y_derivative(order, z, 1, Scaling::Unscaled)
}

/// Computes the $k$-th derivative of the Bessel function of the second kind $\left(\frac{d}{dz}\right)^k Y_\nu(z)$.
///
/// Evaluated via DLMF 10.6.7:
/// $$\left(\frac{d}{dz}\right)^k Y_\nu(z) = \frac{1}{2^k} \sum_{n=0}^k (-1)^n \binom{k}{n} Y_{\nu - k + 2n}(z)$$
///
/// When `scaling` is [`Scaling::Scaled`], the result is scaled by $e^{-|\mathrm{Im}(z)|}$.
///
/// # Arguments
/// * `order` - The order $\nu$ of the Bessel function.
/// * `z` - The complex or real argument.
/// * `derivative_order` - The order of the derivative $k \ge 0$.
/// * `scaling` - Whether to compute the unscaled or exponentially scaled derivative.
///
/// # Examples
/// ```
/// use amos_bessel_rs::{Scaling, derivatives::bessel_y_derivative};
///
/// // Second derivative of Y_0 at z = 1.0:
/// let d2y: f64 = bessel_y_derivative(0.0, 1.0, 2, Scaling::Unscaled).unwrap();
/// assert!((d2y - (-0.8694697855)).abs() < 1e-6);
/// ```
pub fn bessel_y_derivative<FT: BesselFloat, ZT: BesselInput<FT>, OT: Into<FT>>(
    order: OT,
    z: ZT,
    derivative_order: u32,
    scaling: Scaling,
) -> Result<ZT, BesselError<FT>> {
    derivative_internal(
        complex_bessel_y,
        order.into(),
        z.into(),
        derivative_order,
        scaling,
        SignType::Cylinder,
    )
    .and_then(ZT::back_from)
}

/// Computes the first derivative of the modified Bessel function of the first kind $I_\nu'(z)$ with respect to $z$.
///
/// # Arguments
/// * `order` - The order $\nu$ of the Bessel function.
/// * `z` - The complex or real argument.
///
/// # Examples
/// ```
/// use amos_bessel_rs::derivatives::bessel_i_p;
///
/// let di: f64 = bessel_i_p(0.0, 1.0).unwrap();
/// assert!((di - 0.5651591039).abs() < 1e-6);
/// ```
pub fn bessel_i_p<FT: BesselFloat, ZT: BesselInput<FT>, OT: Into<FT>>(
    order: OT,
    z: ZT,
) -> Result<ZT, BesselError<FT>> {
    bessel_i_derivative(order, z, 1, Scaling::Unscaled)
}

/// Computes the $k$-th derivative of the modified Bessel function of the first kind $\left(\frac{d}{dz}\right)^k I_\nu(z)$.
///
/// Evaluated via DLMF 10.29.5:
/// $$\left(\frac{d}{dz}\right)^k I_\nu(z) = \frac{1}{2^k} \sum_{n=0}^k \binom{k}{n} I_{\nu - k + 2n}(z)$$
///
/// When `scaling` is [`Scaling::Scaled`], the result is scaled by $e^{-|\mathrm{Re}(z)|}$.
///
/// # Arguments
/// * `order` - The order $\nu$ of the Bessel function.
/// * `z` - The complex or real argument.
/// * `derivative_order` - The order of the derivative $k \ge 0$.
/// * `scaling` - Whether to compute the unscaled or exponentially scaled derivative.
///
/// # Examples
/// ```
/// use amos_bessel_rs::{Scaling, derivatives::bessel_i_derivative};
///
/// // Second derivative of I_0 at z = 1.0:
/// let d2i: f64 = bessel_i_derivative(0.0, 1.0, 2, Scaling::Unscaled).unwrap();
/// assert!((d2i - 0.7009067738).abs() < 1e-6);
/// ```
pub fn bessel_i_derivative<FT: BesselFloat, ZT: BesselInput<FT>, OT: Into<FT>>(
    order: OT,
    z: ZT,
    derivative_order: u32,
    scaling: Scaling,
) -> Result<ZT, BesselError<FT>> {
    derivative_internal(
        complex_bessel_i,
        order.into(),
        z.into(),
        derivative_order,
        scaling,
        SignType::I,
    )
    .and_then(ZT::back_from)
}

/// Computes the first derivative of the Hankel function $H_\nu^{(1)\prime}(z)$ or $H_\nu^{(2)\prime}(z)$ with respect to $z$.
///
/// # Arguments
/// * `order` - The order $\nu$ of the Hankel function.
/// * `z` - The complex or real argument.
/// * `kind` - The kind of Hankel function ([`HankelKind::First`] or [`HankelKind::Second`]).
///
/// # Examples
/// ```
/// use amos_bessel_rs::{HankelKind, derivatives::hankel_p};
/// use num::Complex;
///
/// let z = Complex::new(1.0, 0.5);
/// let dh = hankel_p(0.0, z, HankelKind::First).unwrap();
/// assert!(dh.norm() > 0.0);
/// ```
pub fn hankel_p<FT: BesselFloat, ZT: BesselInput<FT>, OT: Into<FT>>(
    order: OT,
    z: ZT,
    kind: HankelKind,
) -> Result<ZT, BesselError<FT>> {
    hankel_derivative(order, z, kind, 1, Scaling::Unscaled)
}

/// Computes the $k$-th derivative of the Hankel function $\left(\frac{d}{dz}\right)^k H_\nu^{(1,2)}(z)$.
///
/// Evaluated via DLMF 10.6.7:
/// $$\left(\frac{d}{dz}\right)^k H_\nu^{(m)}(z) = \frac{1}{2^k} \sum_{n=0}^k (-1)^n \binom{k}{n} H_{\nu - k + 2n}^{(m)}(z)$$
///
/// When `scaling` is [`Scaling::Scaled`], the result is scaled by $e^{\mp i z}$.
///
/// # Arguments
/// * `order` - The order $\nu$ of the Hankel function.
/// * `z` - The complex or real argument.
/// * `kind` - The kind of Hankel function ([`HankelKind::First`] or [`HankelKind::Second`]).
/// * `derivative_order` - The order of the derivative $k \ge 0$.
/// * `scaling` - Whether to compute the unscaled or exponentially scaled derivative.
///
/// # Examples
/// ```
/// use amos_bessel_rs::{HankelKind, Scaling, derivatives::hankel_derivative};
/// use num::Complex;
///
/// let z = Complex::new(1.0, 0.5);
/// let d2h = hankel_derivative(0.0, z, HankelKind::First, 2, Scaling::Unscaled).unwrap();
/// assert!(d2h.norm() > 0.0);
/// ```
pub fn hankel_derivative<FT: BesselFloat, ZT: BesselInput<FT>, OT: Into<FT>>(
    order: OT,
    z: ZT,
    kind: HankelKind,
    derivative_order: u32,
    scaling: Scaling,
) -> Result<ZT, BesselError<FT>> {
    match kind {
        HankelKind::First => hankel1_derivative(order, z, derivative_order, scaling),
        HankelKind::Second => hankel2_derivative(order, z, derivative_order, scaling),
    }
}

/// Computes the first derivative of the Hankel function of the first kind $H_\nu^{(1)\prime}(z)$ with respect to $z$.
///
/// # Arguments
/// * `order` - The order $\nu$ of the Hankel function.
/// * `z` - The complex or real argument.
///
/// # Examples
/// ```
/// use amos_bessel_rs::derivatives::hankel1_p;
/// use num::Complex;
///
/// let z = Complex::new(1.0, 0.5);
/// let dh1 = hankel1_p(0.0, z).unwrap();
/// assert!(dh1.norm() > 0.0);
/// ```
pub fn hankel1_p<FT: BesselFloat, ZT: BesselInput<FT>, OT: Into<FT>>(
    order: OT,
    z: ZT,
) -> Result<ZT, BesselError<FT>> {
    hankel1_derivative(order, z, 1, Scaling::Unscaled)
}

/// Computes the $k$-th derivative of the Hankel function of the first kind $\left(\frac{d}{dz}\right)^k H_\nu^{(1)}(z)$.
///
/// Evaluated via DLMF 10.6.7:
/// $$\left(\frac{d}{dz}\right)^k H_\nu^{(1)}(z) = \frac{1}{2^k} \sum_{n=0}^k (-1)^n \binom{k}{n} H_{\nu - k + 2n}^{(1)}(z)$$
///
/// When `scaling` is [`Scaling::Scaled`], the result is scaled by $e^{-i z}$.
///
/// # Arguments
/// * `order` - The order $\nu$ of the Hankel function.
/// * `z` - The complex or real argument.
/// * `derivative_order` - The order of the derivative $k \ge 0$.
/// * `scaling` - Whether to compute the unscaled or exponentially scaled derivative.
///
/// # Examples
/// ```
/// use amos_bessel_rs::{Scaling, derivatives::hankel1_derivative};
/// use num::Complex;
///
/// let z = Complex::new(1.0, 0.5);
/// let d2h1 = hankel1_derivative(0.0, z, 2, Scaling::Unscaled).unwrap();
/// assert!(d2h1.norm() > 0.0);
/// ```
pub fn hankel1_derivative<FT: BesselFloat, ZT: BesselInput<FT>, OT: Into<FT>>(
    order: OT,
    z: ZT,
    derivative_order: u32,
    scaling: Scaling,
) -> Result<ZT, BesselError<FT>> {
    derivative_internal(
        complex_hankel1,
        order.into(),
        z.into(),
        derivative_order,
        scaling,
        SignType::Cylinder,
    )
    .and_then(ZT::back_from)
}

/// Computes the first derivative of the Hankel function of the second kind $H_\nu^{(2)\prime}(z)$ with respect to $z$.
///
/// # Arguments
/// * `order` - The order $\nu$ of the Hankel function.
/// * `z` - The complex or real argument.
///
/// # Examples
/// ```
/// use amos_bessel_rs::derivatives::hankel2_p;
/// use num::Complex;
///
/// let z = Complex::new(1.0, 0.5);
/// let dh2 = hankel2_p(0.0, z).unwrap();
/// assert!(dh2.norm() > 0.0);
/// ```
pub fn hankel2_p<FT: BesselFloat, ZT: BesselInput<FT>, OT: Into<FT>>(
    order: OT,
    z: ZT,
) -> Result<ZT, BesselError<FT>> {
    hankel2_derivative(order, z, 1, Scaling::Unscaled)
}

/// Computes the $k$-th derivative of the Hankel function of the second kind $\left(\frac{d}{dz}\right)^k H_\nu^{(2)}(z)$.
///
/// Evaluated via DLMF 10.6.7:
/// $$\left(\frac{d}{dz}\right)^k H_\nu^{(2)}(z) = \frac{1}{2^k} \sum_{n=0}^k (-1)^n \binom{k}{n} H_{\nu - k + 2n}^{(2)}(z)$$
///
/// When `scaling` is [`Scaling::Scaled`], the result is scaled by $e^{i z}$.
///
/// # Arguments
/// * `order` - The order $\nu$ of the Hankel function.
/// * `z` - The complex or real argument.
/// * `derivative_order` - The order of the derivative $k \ge 0$.
/// * `scaling` - Whether to compute the unscaled or exponentially scaled derivative.
///
/// # Examples
/// ```
/// use amos_bessel_rs::{Scaling, derivatives::hankel2_derivative};
/// use num::Complex;
///
/// let z = Complex::new(1.0, 0.5);
/// let d2h2 = hankel2_derivative(0.0, z, 2, Scaling::Unscaled).unwrap();
/// assert!(d2h2.norm() > 0.0);
/// ```
pub fn hankel2_derivative<FT: BesselFloat, ZT: BesselInput<FT>, OT: Into<FT>>(
    order: OT,
    z: ZT,
    derivative_order: u32,
    scaling: Scaling,
) -> Result<ZT, BesselError<FT>> {
    derivative_internal(
        complex_hankel2,
        order.into(),
        z.into(),
        derivative_order,
        scaling,
        SignType::Cylinder,
    )
    .and_then(ZT::back_from)
}

/// Computes the first derivative of the modified Bessel function of the second kind $K_\nu'(z)$ with respect to $z$.
///
/// # Arguments
/// * `order` - The order $\nu$ of the Bessel function.
/// * `z` - The complex or real argument.
///
/// # Examples
/// ```
/// use amos_bessel_rs::derivatives::bessel_k_p;
///
/// let dk: f64 = bessel_k_p(0.0, 1.0).unwrap();
/// assert!((dk - (-0.6019072301)).abs() < 1e-6);
/// ```
pub fn bessel_k_p<FT: BesselFloat, ZT: BesselInput<FT>, OT: Into<FT>>(
    order: OT,
    z: ZT,
) -> Result<ZT, BesselError<FT>> {
    bessel_k_derivative(order, z, 1, Scaling::Unscaled)
}

/// Computes the $k$-th derivative of the modified Bessel function of the second kind $\left(\frac{d}{dz}\right)^k K_\nu(z)$.
///
/// Evaluated via DLMF 10.29.5:
/// $$\left(\frac{d}{dz}\right)^k K_\nu(z) = \frac{(-1)^k}{2^k} \sum_{n=0}^k \binom{k}{n} K_{\nu - k + 2n}(z)$$
///
/// When `scaling` is [`Scaling::Scaled`], the result is scaled by $e^{z}$.
///
/// # Arguments
/// * `order` - The order $\nu$ of the Bessel function.
/// * `z` - The complex or real argument.
/// * `derivative_order` - The order of the derivative $k \ge 0$.
/// * `scaling` - Whether to compute the unscaled or exponentially scaled derivative.
///
/// # Examples
/// ```
/// use amos_bessel_rs::{Scaling, derivatives::bessel_k_derivative};
///
/// // Second derivative of K_0 at z = 1.0:
/// let d2k: f64 = bessel_k_derivative(0.0, 1.0, 2, Scaling::Unscaled).unwrap();
/// assert!((d2k - 1.0229316684).abs() < 1e-6);
/// ```
pub fn bessel_k_derivative<FT: BesselFloat, ZT: BesselInput<FT>, OT: Into<FT>>(
    order: OT,
    z: ZT,
    derivative_order: u32,
    scaling: Scaling,
) -> Result<ZT, BesselError<FT>> {
    derivative_internal(
        complex_bessel_k,
        order.into(),
        z.into(),
        derivative_order,
        scaling,
        SignType::K,
    )
    .and_then(ZT::back_from)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum SignType {
    I,
    K,
    Cylinder,
}

/// Implements DLMF 10.6.7 and 10.29.5
fn derivative_internal<T: BesselFloat>(
    func: BesselSig<T>,
    order: T,
    z: Complex<T>,
    derivative_order: u32,
    scaling: Scaling,
    sign_type: SignType,
) -> Result<Complex<T>, BesselError<T>> {
    if derivative_order > 60 {
        return Err(BesselError::InvalidInput {
            details: "Derivative order too large - must be no greater than 60",
        });
    }
    let k = derivative_order as usize;
    let mut prefactor = T::ONE / T::TWO.powi(k as i32);
    if sign_type == SignType::K {
        prefactor *= integer_sign::<T>(k as i64);
    }

    let values = match func(z, order - T::from_usize(k), scaling, 2 * k + 1) {
        Ok((values, _n_zeros)) => values,
        Err(BesselError::PartialLossOfSignificance { y, .. }) => y,
        Err(err) => return Err(err),
    };

    let mut sum = T::C_ZERO;
    for n in 0..=k {
        let n_choose_k = T::from_f64(binomial(k as u64, n as u64) as f64);
        let sign = match sign_type {
            SignType::I => T::ONE,
            SignType::K => T::ONE,
            SignType::Cylinder => integer_sign::<T>(n as i64),
        };
        let v = values[n * 2];
        sum += v * sign * n_choose_k;
    }
    Ok(prefactor * sum)
}
