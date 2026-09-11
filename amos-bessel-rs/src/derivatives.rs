use num::Complex;
use num_integer::binomial;

use crate::{
    BesselError, BesselFloat, BesselInput, Scaling, amos::complex_bessel_j,
    reflections::integer_sign, types::BesselValues,
};

#[allow(type_alias_bounds)]
type BesselSig<T: BesselFloat = f64> =
    fn(Complex<T>, T, Scaling, usize) -> Result<BesselValues<T>, BesselError<T>>;

pub fn bessel_j_p<FT: BesselFloat, ZT: BesselInput<FT>, OT: Into<FT>>(
    order: OT,
    z: ZT,
) -> Result<ZT, BesselError<FT>> {
    bessel_j_derivative(order, z, 1)
}

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
        order - T::from_f64(k as f64),
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
