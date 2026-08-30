use num::Complex;

pub use super::core::{complex_airy, complex_airy_b};
use crate::{
    BesselFloat, Scaling,
    amos::HankelKind,
    reflections::{as_integer, integer_sign, reflect_j_element},
    types::BesselResult,
};

/// Computes the H-Bessel functions (Hankel functions) of a complex argument.
///
/// This function computes a sequence of complex Hankel (Bessel) functions
/// `y[j] = H(order + j - 1, z)` real, non-negative
/// orders `order + j - 1` (`j = 1, ..., n`), and a complex argument `z` which is
/// not equal to `(0.0, 0.0)`. The computation is valid in the cut plane
/// `-PI < z.arg() <= PI`.
///
/// The kind of the Hankel function is specified by the hankel_kind parameter,
/// which can take values [HankelKind::First] or [HankelKind::Second]
///
/// When `scaling` is `Scaling::Scaled`, this function returns the scaled Hankel
/// functions, which remove the exponential behavior in both the upper and
/// lower half-planes.
/// `y(j) = (-(3 - 2 * m)*z*i).exp() * H(order + j - 1, z)` where `m` depends
/// on the kind of Hankel function (1 for First, 2 for second).
///
/// # Arguments
///
/// * `z` - Complex argument `z`, `z != (0.0, 0.0)`, `-PI < z.arg() <= PI`.
/// * `order` - Order of the initial H function, `order >= 0.0`.
/// * `scaling` - A parameter to indicate the scaling option.
///     * `Scaling::Unscaled`: returns `y(j) = H(order + j - 1, z)`.
///     * `Scaling::Scaled`: returns `y(j) = H(m, order + j - 1, z) * (-i * z * (3 - 2*m)).exp()`
///       where `m` is determined by the kind of Hankel function (1 for First, 2 for second).
/// * `hankel_kind` - Kind of Hankel function.
/// * `n` - Number of members in the sequence, `n >= 1`.
///
/// # Returns
///
/// A tuple containing:
/// * `y`: A vector of complex numbers containing the values of the Hankel
///   functions for orders `[order, order + 1, ..., order + n - 1]`.
/// * `n_zeros`: The number of components in `y` set to zero due to underflow.
///   for TODO
pub fn complex_bessel_h<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    hankel_kind: HankelKind,
    n: usize,
) -> BesselResult<T> {
    super::core::complex_bessel_h(z, order, scaling, hankel_kind, n)
}

/// Computes the Hankel function of the first kind H1v(z) for a complex argument.
///
/// An alternative interface to [`complex_bessel_h`]
#[inline]
pub fn complex_hankel1<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    n: usize,
) -> BesselResult<T> {
    complex_bessel_h(z, order, scaling, HankelKind::First, n)
}

/// Computes the Hankel function of the second kind H2v(z) for a complex argument.
///
/// An alternative interface to [`complex_bessel_h`]
#[inline]
pub fn complex_hankel2<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    n: usize,
) -> BesselResult<T> {
    complex_bessel_h(z, order, scaling, HankelKind::Second, n)
}

/// Computes the I-Bessel function of a complex argument.
///
/// This function computes a sequence of complex Bessel functions `y(j) = I(order + j - 1, z)`
/// for real, non-negative orders `order + j - 1` (`j = 1, ..., n`) and a complex argument `z`
/// in the cut plane `-PI < z.arg() <= PI`.
///
/// When `scaling` is `Scaling::Scaled`, this function returns the scaled functions
/// `y(j) = (-(z.re.abs())).exp() * I(order + j - 1, z)` which remove the
/// exponential growth in both the left and right half-planes for `z` to infinity.
///
/// The computation is carried out by the power series for small `z.abs()`,
/// the asymptotic expansion for large `z.abs()`,
/// the Miller algorithm normalized by the Wronskian and a Neumann
/// series for intermediate magnitudes, and the
/// uniform asymptotic expansions for I(order, z) and J(order, z)
/// for large orders. Backward recurrence is used to generate
/// sequences or reduce orders when necessary.
///
/// The calculations above are done in the right half-plane and
/// continued into the left half-plane by the formula
/// `I(order, z * (m * PI).exp()) = (m * PI * order).exp() * I(order, z),   z.re > 0.0`
/// with `m = +i OR -i`,  (`i` is the imaginary unit).
//
///
/// # Arguments
///
/// * `z` - Complex argument `z`, `-PI < z.arg() <= PI`.
/// * `order` - Order of the initial I function, `order >= 0.0`.
/// * `scaling` - A parameter to indicate the scaling option.
///     * `Scaling::Unscaled`: returns `y(j) = I(order + j - 1, z)`.
///     * `Scaling::Scaled`: returns `y(j) = I(order + j - 1, z) * (-z.re().abs()).exp()`.
/// * `n` - Number of members of the sequence, `n >= 1`.
///
/// # Returns
///
/// A tuple containing:
/// * `y`: A vector of complex numbers containing the values of the Bessel
///   functions for orders `[order, order + 1, ..., order + n - 1]`.
/// * `n_zeros`: The number of components in `y` set to zero due to underflow.
///   For `complex_bessel_i`, the first `n_zeros` components are set to zero.
pub fn complex_bessel_i<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    n: usize,
) -> BesselResult<T, usize> {
    super::core::complex_bessel_i(z, order, scaling, n)
}

/// Computes the J-Bessel function of a complex argument.
///
/// This function computes a sequence of complex Bessel functions `y(j) = J(order + j - 1, z)`
/// for real, non-negative orders `order + j - 1` (`j = 1, ..., n`) and a complex argument `z`
/// in the cut plane `-PI < z.arg() <= PI`.
///
/// When `scaling` is `Scaling::Scaled`, this function returns the scaled functions
/// `y(j) = (-(z.im.abs())).exp() * J(order + j - 1, z)`, which removes the
/// exponential growth in both the upper and lower half-planes for `z` to infinity.
///
/// The computation is carried out by the formula
///
/// `J(order, Z) = ( order * PI * i / 2).exp() * I(order, -i*z)`    if `z.im >= 0.0`
///
/// `J(order, Z) = (-order * PI * i / 2).exp() * I(order, i*z)`    if `z.im < 0.0`
///
/// where `i` is the imaginary unit and `I(order, z)` is the I Bessel function.
///
/// # Arguments
///
/// * `z` - Complex argument `z`, `-PI < z.arg() <= PI`.
/// * `order` - Order of the initial J function, `order >= 0.0`.
/// * `scaling` - A parameter to indicate the scaling option.
///     * `Scaling::Unscaled`: returns `y(j) = J(order + j - 1, z)`.
///     * `Scaling::Scaled`: returns `y(j) = J(order + j - 1, z) * (-(z.im.abs())).exp()`.
/// * `n` - Number of members of the sequence, `n >= 1`.
///
/// # Returns
///
/// A tuple containing:
/// * `y`: A vector of complex numbers containing the values of the Bessel
///   functions for orders `[order, order + 1, ..., order + n - 1]`.
/// * `n_zeros`: The number of components in `y` set to zero due to underflow.
///   For `complex_bessel_j`, the last `n_zeros` components are set to zero.
///
pub fn complex_bessel_j<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    n: usize,
) -> BesselResult<T> {
    use super::core::complex_bessel_j;
    if order >= T::ZERO {
        // implies all requested orders are positive
        return complex_bessel_j(z, order, scaling, n);
    }

    // Special case for negative integer order: J(-n, z) = (-1)^n J(n, z)
    let abs_order: T = order.abs();
    if let Some(int_order) = as_integer(abs_order) {
        // if we have a negative integer order, then orders are
        // (if int_order is denoted by o)
        // -o, -o+1, -o+2, ..., 0, 1, 2, ... n-(o+1)
        // e.g for order = -3, int_order = 3, n = 5 orders are -3, -2, -1, 0, 1,
        // of course, if n < int_order, then we never reach order 0
        // -o, -o+1, -o+2, ..., -o+n
        // e.g. if n = 2, int_order = 3, orders are -3, -2,

        // now we need positive forms of all the orders we need in either
        // negative or positive form.

        let n64 = n as i64;
        let min_order = (-int_order + n64).abs().min(0);
        let max_order = (n64 - (int_order + 1)).max(int_order);
        let n_positive = ((max_order - min_order) + 1) as usize;
        let (pos_js, n_zeros_positive) =
            complex_bessel_j(z, T::from_isize(min_order as isize), scaling, n_positive)?;

        let mut answer = Vec::with_capacity(n);

        let mut n_zeros = 0;
        for i in 0..n {
            let current_order = order + T::from_usize(i);
            let abs_order = current_order.abs().to_usize().unwrap();
            if abs_order >= (n_positive - n_zeros_positive) {
                n_zeros += 1;
            }
            if current_order < T::ZERO {
                answer.push(pos_js[abs_order] * integer_sign::<T>(abs_order as i64));
            } else {
                answer.push(pos_js[abs_order]);
            }
        }
        return Ok((answer, n_zeros));
    }

    // General case: need both J and Y at positive |ν|
    //
    // say order = -2.7, then we need
    // -2.7, -1.7, -0.7, 0.3, 1.3, 2.3, ...

    let n_negative = order.abs().ceil().abs().to_usize().unwrap();
    let first_negative = order.fract().abs();
    let (j_negative, n_zeros_j_neg) = complex_bessel_j(z, first_negative, scaling, n_negative)?;
    let (y_negative, n_zeros_y_neg) = complex_bessel_y(z, first_negative, scaling, n_negative)?;

    let n_positive = n.saturating_sub(n_negative);

    let first_positive = T::ONE + order.fract();
    let (j_positive, n_zeros_j_pos) = if n_positive > 0 {
        complex_bessel_j(z, first_positive, scaling, n_positive)?
    } else {
        (Vec::new(), 0)
    };

    let mut answer = Vec::with_capacity(n);
    let mut n_zeros = 0;
    for i in 0..n {
        let current_order = order + T::from_usize(i);
        if current_order < T::ZERO {
            let index = n_negative - 1 - i;
            answer.push(reflect_j_element(
                current_order.abs(),
                j_negative[index],
                y_negative[index],
            ));
            if (index >= n_negative - n_zeros_j_neg) && index >= n_negative - n_zeros_y_neg {
                n_zeros += 1;
            }
        } else {
            let index = i - n_negative;
            answer.push(j_positive[index]);
            if index >= n_positive - n_zeros_j_pos {
                n_zeros += 1;
            }
        }
    }
    Ok((answer, n_zeros))
}

/// Computes the K-Bessel function of a complex argument.
///
/// This function computes a sequence of complex Bessel functions `y(j) = K(order + j - 1, z)`
/// for real, non-negative orders `order + j - 1` (`j = 1, ..., n`) and a complex argument `z`
/// which is not equal to `(0.0, 0.0)`. The computation is valid in the cut plane
/// `-PI < z.arg() <= PI`.
///
/// When `scaling` is `Scaling::Scaled`, this function returns the scaled K functions,
/// `y(j) = z.exp() * K(order + j - 1, z)`, which remove the exponential behavior in both
/// the left and right half-planes for `z` to infinity.
///
/// EQUATIONS ARE IMPLEMENTED FOR SMALL ORDERS
/// order AND order + 1.0 IN THE RIGHT HALF PLANE X >= 0.0. FORWARD
/// RECURRENCE GENERATES HIGHER ORDERS. K IS CONTINUED TO THE LEFT
/// HALF PLANE BY THE RELATION
///
/// `K(order, z * mp.exp()) = (-mp * order).exp() * K(order, z) - mp * I(order, z)`
///
/// where `mp = mr * PI * i`, `mr = +1 OR -1`, `z.re > 0`, `i` is the imaginary unit
/// and `I(order, Z)` is the I Bessel function.
///
/// For large order, `order > MACHINE_CONSTANTS.asymptotic_order_limit`, the K function is computed
/// by means of its uniform asymptotic expansions.
///
/// # Arguments
///
/// * `z` - Complex argument `z`, `z != (0.0, 0.0)`, `-PI < z.arg() <= PI`.
/// * `order` - Order of the initial K function, `order >= 0.0`.
/// * `scaling` - A parameter to indicate the scaling option.
///     * `Scaling::Unscaled`: returns `y(j) = K(order + j - 1, z)`.
///     * `Scaling::Scaled`: returns `y(j) = K(order + j - 1, z) * z.exp()`.
/// * `n` - Number of members of the sequence, `n >= 1`.
///
/// # Returns
///
/// A tuple containing:
/// * `y`: A vector of complex numbers containing the values of the Bessel
///   functions for orders `[order, order + 1, ..., order + n - 1]`.
/// * `n_zeros`: The number of components in `y` set to zero due to underflow.
///   For `complex_bessel_k`, the first `n_zeros` components are set to zero.
pub fn complex_bessel_k<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    n: usize,
) -> BesselResult<T> {
    // K_{-ν}(z) = K_ν(z) (DLMF 10.27.3)
    let abs_order = order.abs();
    super::core::complex_bessel_k(z, abs_order, scaling, n)
}

/// Computes the Y-Bessel function of a complex argument.
///
/// This function computes a sequence of complex Bessel functions `y(j) = Y(order + j - 1, z)`
/// for real, non-negative orders `order + j - 1` (`j = 1, ..., n`) and a complex argument `z`
/// which is not equal to `(0.0, 0.0)`. The computation is valid in the cut plane
/// `-PI < z.arg() <= PI`.
///
/// When `scaling` is `Scaling::Scaled`, this function returns the scaled functions
/// `y(j) = (-(z.im.abs())).exp() * Y(order + j - 1, z)`, which remove the
/// exponential growth in both the upper and lower half-planes for `z` to infinity.
///
/// The computation is carried out in terms of the I(order, z) and
/// K(order, z) Bessel functions in the right half-plane by
///
/// `Y(order, z) = i * cc * I(order, arg) - (2/PI) * cc.conj() * K(order, arg)` if `z.im >= 0`
///
/// `Y(order, z) = Y(order, z.conj()).conj()` if `z.im < 0`
///
/// where
/// `cc = (i* PI * order / 2).exp()`, `arg = z * (-i * PI / 2).exp()` and `i` is the imaginary unit.
///
/// # Arguments
///
/// * `z` - Complex argument `z`, `z != (0.0, 0.0)`, `-PI < z.arg() <= PI`.
/// * `order` - Order of the initial Y function, `order >= 0.0`.
/// * `scaling` - A parameter to indicate the scaling option.
///     * `Scaling::Unscaled`: returns `y(j) = Y(order + j - 1, z)`.
///     * `Scaling::Scaled`: returns `y(j) = Y(order + j - 1, z) * (-(z.im.abs())).exp()`.
/// * `n` - Number of members of the sequence, `n >= 1`.
///
/// # Returns
///
/// A tuple containing:
/// * `y`: A vector of complex numbers containing the values of the Bessel
///   functions for orders `[order, order + 1, ..., order + n - 1]`.
/// * `n_zeros`: The number of components in `y` set to zero due to underflow.
///   For `complex_bessel_y`, the first `n_zeros` components are set to zero.
pub fn complex_bessel_y<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    n: usize,
) -> BesselResult<T> {
    super::core::complex_bessel_y(z, order, scaling, n)
}
