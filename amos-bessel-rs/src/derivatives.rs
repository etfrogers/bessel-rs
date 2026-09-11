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
pub fn bessel_j_p<FT: BesselFloat, ZT: BesselInput<FT>, OT: Into<FT>>(
    order: OT,
    z: ZT,
) -> Result<ZT, BesselError<FT>> {
    bessel_j_derivative(order, z, 1)
}

/// Computes the $k$-th derivative of the Bessel function of the first kind $\left(\frac{d}{dz}\right)^k J_\nu(z)$.
///
/// Evaluated via the exact binomial identity from DLMF 10.6.7:
/// $$\left(\frac{d}{dz}\right)^k J_\nu(z) = \frac{1}{2^k} \sum_{n=0}^k (-1)^n \binom{k}{n} J_{\nu - k + 2n}(z)$$
///
/// # Arguments
/// * `order` - The order $\nu$ of the Bessel function.
/// * `z` - The complex or real argument.
/// * `derivative_order` - The order of the derivative $k \ge 0$.
pub fn bessel_j_derivative<FT: BesselFloat, ZT: BesselInput<FT>, OT: Into<FT>>(
    order: OT,
    z: ZT,
    derivative_order: u32,
) -> Result<ZT, BesselError<FT>> {
    let order: FT = order.into();
    let z: Complex<FT> = z.into();
    ZT::back_from(&derivative_internal(
        complex_bessel_j,
        order,
        z,
        derivative_order,
        SignType::Cylinder,
    )?)
}

/// Computes the first derivative of the Bessel function of the second kind $Y_\nu'(z)$ with respect to $z$.
///
/// # Arguments
/// * `order` - The order $\nu$ of the Bessel function.
/// * `z` - The complex or real argument.
pub fn bessel_y_p<FT: BesselFloat, ZT: BesselInput<FT>, OT: Into<FT>>(
    order: OT,
    z: ZT,
) -> Result<ZT, BesselError<FT>> {
    bessel_y_derivative(order, z, 1)
}

/// Computes the $k$-th derivative of the Bessel function of the second kind $\left(\frac{d}{dz}\right)^k Y_\nu(z)$.
///
/// Evaluated via DLMF 10.6.7:
/// $$\left(\frac{d}{dz}\right)^k Y_\nu(z) = \frac{1}{2^k} \sum_{n=0}^k (-1)^n \binom{k}{n} Y_{\nu - k + 2n}(z)$$
///
/// # Arguments
/// * `order` - The order $\nu$ of the Bessel function.
/// * `z` - The complex or real argument.
/// * `derivative_order` - The order of the derivative $k \ge 0$.
pub fn bessel_y_derivative<FT: BesselFloat, ZT: BesselInput<FT>, OT: Into<FT>>(
    order: OT,
    z: ZT,
    derivative_order: u32,
) -> Result<ZT, BesselError<FT>> {
    let order: FT = order.into();
    let z: Complex<FT> = z.into();
    ZT::back_from(&derivative_internal(
        complex_bessel_y,
        order,
        z,
        derivative_order,
        SignType::Cylinder,
    )?)
}

/// Computes the first derivative of the modified Bessel function of the first kind $I_\nu'(z)$ with respect to $z$.
///
/// # Arguments
/// * `order` - The order $\nu$ of the Bessel function.
/// * `z` - The complex or real argument.
pub fn bessel_i_p<FT: BesselFloat, ZT: BesselInput<FT>, OT: Into<FT>>(
    order: OT,
    z: ZT,
) -> Result<ZT, BesselError<FT>> {
    bessel_i_derivative(order, z, 1)
}

/// Computes the $k$-th derivative of the modified Bessel function of the first kind $\left(\frac{d}{dz}\right)^k I_\nu(z)$.
///
/// Evaluated via DLMF 10.29.5:
/// $$\left(\frac{d}{dz}\right)^k I_\nu(z) = \frac{1}{2^k} \sum_{n=0}^k \binom{k}{n} I_{\nu - k + 2n}(z)$$
///
/// # Arguments
/// * `order` - The order $\nu$ of the Bessel function.
/// * `z` - The complex or real argument.
/// * `derivative_order` - The order of the derivative $k \ge 0$.
pub fn bessel_i_derivative<FT: BesselFloat, ZT: BesselInput<FT>, OT: Into<FT>>(
    order: OT,
    z: ZT,
    derivative_order: u32,
) -> Result<ZT, BesselError<FT>> {
    let order: FT = order.into();
    let z: Complex<FT> = z.into();
    ZT::back_from(&derivative_internal(
        complex_bessel_i,
        order,
        z,
        derivative_order,
        SignType::I,
    )?)
}

/// Computes the first derivative of the Hankel function $H_\nu^{(1)\prime}(z)$ or $H_\nu^{(2)\prime}(z)$ with respect to $z$.
///
/// # Arguments
/// * `order` - The order $\nu$ of the Hankel function.
/// * `z` - The complex or real argument.
/// * `kind` - The kind of Hankel function ([`HankelKind::First`] or [`HankelKind::Second`]).
pub fn hankel_p<FT: BesselFloat, ZT: BesselInput<FT>, OT: Into<FT>>(
    order: OT,
    z: ZT,
    kind: HankelKind,
) -> Result<ZT, BesselError<FT>> {
    hankel_derivative(order, z, kind, 1)
}

/// Computes the $k$-th derivative of the Hankel function $\left(\frac{d}{dz}\right)^k H_\nu^{(1,2)}(z)$.
///
/// Evaluated via DLMF 10.6.7:
/// $$\left(\frac{d}{dz}\right)^k H_\nu^{(m)}(z) = \frac{1}{2^k} \sum_{n=0}^k (-1)^n \binom{k}{n} H_{\nu - k + 2n}^{(m)}(z)$$
///
/// # Arguments
/// * `order` - The order $\nu$ of the Hankel function.
/// * `z` - The complex or real argument.
/// * `kind` - The kind of Hankel function ([`HankelKind::First`] or [`HankelKind::Second`]).
/// * `derivative_order` - The order of the derivative $k \ge 0$.
pub fn hankel_derivative<FT: BesselFloat, ZT: BesselInput<FT>, OT: Into<FT>>(
    order: OT,
    z: ZT,
    kind: HankelKind,
    derivative_order: u32,
) -> Result<ZT, BesselError<FT>> {
    match kind {
        HankelKind::First => hankel1_derivative(order, z, derivative_order),
        HankelKind::Second => hankel2_derivative(order, z, derivative_order),
    }
}

/// Computes the first derivative of the Hankel function of the first kind $H_\nu^{(1)\prime}(z)$ with respect to $z$.
///
/// # Arguments
/// * `order` - The order $\nu$ of the Hankel function.
/// * `z` - The complex or real argument.
pub fn hankel1_p<FT: BesselFloat, ZT: BesselInput<FT>, OT: Into<FT>>(
    order: OT,
    z: ZT,
) -> Result<ZT, BesselError<FT>> {
    hankel1_derivative(order, z, 1)
}

/// Computes the $k$-th derivative of the Hankel function of the first kind $\left(\frac{d}{dz}\right)^k H_\nu^{(1)}(z)$.
///
/// Evaluated via DLMF 10.6.7:
/// $$\left(\frac{d}{dz}\right)^k H_\nu^{(1)}(z) = \frac{1}{2^k} \sum_{n=0}^k (-1)^n \binom{k}{n} H_{\nu - k + 2n}^{(1)}(z)$$
///
/// # Arguments
/// * `order` - The order $\nu$ of the Hankel function.
/// * `z` - The complex or real argument.
/// * `derivative_order` - The order of the derivative $k \ge 0$.
pub fn hankel1_derivative<FT: BesselFloat, ZT: BesselInput<FT>, OT: Into<FT>>(
    order: OT,
    z: ZT,
    derivative_order: u32,
) -> Result<ZT, BesselError<FT>> {
    let order: FT = order.into();
    let z: Complex<FT> = z.into();
    ZT::back_from(&derivative_internal(
        complex_hankel1,
        order,
        z,
        derivative_order,
        SignType::Cylinder,
    )?)
}

/// Computes the first derivative of the Hankel function of the second kind $H_\nu^{(2)\prime}(z)$ with respect to $z$.
///
/// # Arguments
/// * `order` - The order $\nu$ of the Hankel function.
/// * `z` - The complex or real argument.
pub fn hankel2_p<FT: BesselFloat, ZT: BesselInput<FT>, OT: Into<FT>>(
    order: OT,
    z: ZT,
) -> Result<ZT, BesselError<FT>> {
    hankel2_derivative(order, z, 1)
}

/// Computes the $k$-th derivative of the Hankel function of the second kind $\left(\frac{d}{dz}\right)^k H_\nu^{(2)}(z)$.
///
/// Evaluated via DLMF 10.6.7:
/// $$\left(\frac{d}{dz}\right)^k H_\nu^{(2)}(z) = \frac{1}{2^k} \sum_{n=0}^k (-1)^n \binom{k}{n} H_{\nu - k + 2n}^{(2)}(z)$$
///
/// # Arguments
/// * `order` - The order $\nu$ of the Hankel function.
/// * `z` - The complex or real argument.
/// * `derivative_order` - The order of the derivative $k \ge 0$.
pub fn hankel2_derivative<FT: BesselFloat, ZT: BesselInput<FT>, OT: Into<FT>>(
    order: OT,
    z: ZT,
    derivative_order: u32,
) -> Result<ZT, BesselError<FT>> {
    let order: FT = order.into();
    let z: Complex<FT> = z.into();
    ZT::back_from(&derivative_internal(
        complex_hankel2,
        order,
        z,
        derivative_order,
        SignType::Cylinder,
    )?)
}

/// Computes the first derivative of the modified Bessel function of the second kind $K_\nu'(z)$ with respect to $z$.
///
/// # Arguments
/// * `order` - The order $\nu$ of the Bessel function.
/// * `z` - The complex or real argument.
pub fn bessel_k_p<FT: BesselFloat, ZT: BesselInput<FT>, OT: Into<FT>>(
    order: OT,
    z: ZT,
) -> Result<ZT, BesselError<FT>> {
    bessel_k_derivative(order, z, 1)
}

/// Computes the $k$-th derivative of the modified Bessel function of the second kind $\left(\frac{d}{dz}\right)^k K_\nu(z)$.
///
/// Evaluated via DLMF 10.29.5:
/// $$\left(\frac{d}{dz}\right)^k K_\nu(z) = \frac{(-1)^k}{2^k} \sum_{n=0}^k \binom{k}{n} K_{\nu - k + 2n}(z)$$
///
/// # Arguments
/// * `order` - The order $\nu$ of the Bessel function.
/// * `z` - The complex or real argument.
/// * `derivative_order` - The order of the derivative $k \ge 0$.
pub fn bessel_k_derivative<FT: BesselFloat, ZT: BesselInput<FT>, OT: Into<FT>>(
    order: OT,
    z: ZT,
    derivative_order: u32,
) -> Result<ZT, BesselError<FT>> {
    let order: FT = order.into();
    let z: Complex<FT> = z.into();
    ZT::back_from(&derivative_internal(
        complex_bessel_k,
        order,
        z,
        derivative_order,
        SignType::K,
    )?)
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
    sign_type: SignType,
) -> Result<Complex<T>, BesselError<T>> {
    if derivative_order > 60 {
        return Err(BesselError::InvalidInput {
            details: format!(
                "Derivative order {derivative_order} too large - must be no greater than 60"
            ),
        });
    }
    let k = derivative_order as usize;
    let prefactor = T::ONE / T::TWO.powi(k as i32);
    let (values, _n_zeros) = func(z, order - T::from_usize(k), Scaling::Unscaled, 2 * k + 1)?;

    let mut sum = T::C_ZERO;
    let k_sign = integer_sign::<T>(k as i64);
    for n in 0..=k {
        let n_choose_k = T::from_usize(binomial(k, n));
        let sign = match sign_type {
            SignType::I => T::ONE,
            SignType::K => k_sign,
            SignType::Cylinder => integer_sign::<T>(n as i64),
        };
        let v = values[n * 2];
        sum += v * sign * n_choose_k;
    }
    Ok(prefactor * sum)
}
