use num::Complex;
use num_integer::binomial;

use crate::{
    BesselError, BesselFloat, BesselInput, Scaling, amos::complex_bessel_j,
    reflections::integer_sign, types::BesselValues,
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
    )?)
}

fn derivative_internal<T: BesselFloat>(
    func: BesselSig<T>,
    order: T,
    z: Complex<T>,
    derivative_order: u32,
) -> Result<Complex<T>, BesselError<T>> {
    let k = derivative_order as usize;
    let prefactor = T::ONE / T::TWO.powi(k as i32);
    let (values, _n_zeros) = func(
        z,
        order - T::from_usize(k),
        Scaling::Unscaled,
        2 * k + 1,
    )?;

    let mut sum = T::C_ZERO;
    for n in 0..=k {
        let n_choose_k = T::from_usize(binomial(k, n));
        let sign = integer_sign::<T>(n as i64);
        let v = values[n * 2];
        sum += v * sign * n_choose_k;
    }
    Ok(prefactor * sum)
}
