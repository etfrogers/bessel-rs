use num::Complex;

pub use super::core::{complex_airy, complex_airy_b};
use crate::{
    BesselFloat, Scaling,
    amos::HankelKind,
    reflections::{BesselI, BesselJ, BesselK, BesselY, Hankel, reflect_orders},
    types::BesselResult,
};

/// Computes the Hankel function $H_\nu^{(1)}(z)$ or $H_\nu^{(2)}(z)$ for a complex argument.
///
/// Computes a sequence of complex Hankel (Bessel of the third kind) functions
/// `y[j] = H(order + j, z)` for real orders `order + j` (`j = 0, ..., n - 1`) and complex
/// argument `z` ($z \ne 0$) in the cut plane $-\pi < \arg(z) \le \pi$.
///
/// The kind of Hankel function is specified by `hankel_kind` ([HankelKind::First] or [HankelKind::Second]).
/// Negative orders are supported via the reflection relation (DLMF 10.4.6/7):
///
/// $$H_{-\nu}^{(1)}(z) = e^{i\nu\pi} H_\nu^{(1)}(z), \quad H_{-\nu}^{(2)}(z) = e^{-i\nu\pi} H_\nu^{(2)}(z)$$
///
/// # Arguments
///
/// * `z` - Complex argument $z \ne 0$ in the cut plane $-\pi < \arg(z) \le \pi$.
/// * `order` - Order of the initial Hankel function (any real number).
/// * `scaling` - Scaling option:
///     * [Scaling::Unscaled]: returns $H_\nu^{(m)}(z)$.
///     * [Scaling::Scaled]: returns $\exp(-i z (3 - 2m)) H_\nu^{(m)}(z)$ where $m \in \{1, 2\}$,
///       which removes the exponential growth in upper and lower half-planes.
/// * `hankel_kind` - [HankelKind::First] ($H_\nu^{(1)}$) or [HankelKind::Second] ($H_\nu^{(2)}$).
/// * `n` - Number of members in the sequence ($n \ge 1$).
///
/// # Returns
///
/// A `Result` containing a tuple `(y, n_zeros)`:
/// * `y`: A vector of complex values for orders `[order, order + 1, ..., order + n - 1]`.
/// * `n_zeros`: The number of components at the start of `y` set to zero due to underflow.
pub fn complex_bessel_h<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    hankel_kind: HankelKind,
    n: usize,
) -> BesselResult<T> {
    reflect_orders(z, order, scaling, n, Hankel(hankel_kind))
}

/// Computes the Hankel function of the first kind $H_\nu^{(1)}(z)$ for a complex argument.
///
/// Convenience wrapper equivalent to calling [`complex_bessel_h`] with [`HankelKind::First`].
#[inline]
pub fn complex_hankel1<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    n: usize,
) -> BesselResult<T> {
    complex_bessel_h(z, order, scaling, HankelKind::First, n)
}

/// Computes the Hankel function of the second kind $H_\nu^{(2)}(z)$ for a complex argument.
///
/// Convenience wrapper equivalent to calling [`complex_bessel_h`] with [`HankelKind::Second`].
#[inline]
pub fn complex_hankel2<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    n: usize,
) -> BesselResult<T> {
    complex_bessel_h(z, order, scaling, HankelKind::Second, n)
}

/// Computes the modified Bessel function of the first kind $I_\nu(z)$ for a complex argument.
///
/// Computes a sequence of complex modified Bessel functions `y[j] = I(order + j, z)`
/// for real orders `order + j` (`j = 0, ..., n - 1`) and complex argument `z`
/// in the cut plane $-\pi < \arg(z) \le \pi$.
///
/// Negative orders are evaluated via the DLMF reflection formulas (DLMF 10.27.2):
///
/// $$I_{-\nu}(z) = I_\nu(z) + \frac{2}{\pi}\sin(\nu\pi)K_\nu(z), \quad I_{-n}(z) = I_n(z) \ (n \in \mathbb{Z})$$
///
/// # Arguments
///
/// * `z` - Complex argument in the cut plane $-\pi < \arg(z) \le \pi$.
/// * `order` - Order of the initial I function (any real number).
/// * `scaling` - Scaling option:
///     * [Scaling::Unscaled]: returns $I_\nu(z)$.
///     * [Scaling::Scaled]: returns $\exp(-|\text{Re}(z)|) I_\nu(z)$, which removes exponential
///       growth in both left and right half-planes as $|z| \to \infty$.
/// * `n` - Number of members in the sequence ($n \ge 1$).
///
/// # Returns
///
/// A `Result` containing a tuple `(y, n_zeros)`:
/// * `y`: A vector of complex values for orders `[order, order + 1, ..., order + n - 1]`.
/// * `n_zeros`: The number of components at the end of `y` (highest orders) set to zero due to underflow.
pub fn complex_bessel_i<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    n: usize,
) -> BesselResult<T, usize> {
    reflect_orders(z, order, scaling, n, BesselI)
}

/// Computes the Bessel function of the first kind $J_\nu(z)$ for a complex argument.
///
/// Computes a sequence of complex Bessel functions `y[j] = J(order + j, z)`
/// for real orders `order + j` (`j = 0, ..., n - 1`) and complex argument `z`
/// in the cut plane $-\pi < \arg(z) \le \pi$.
///
/// Negative orders are evaluated via the DLMF reflection formulas (DLMF 10.2.3):
///
/// $$J_{-\nu}(z) = \cos(\nu\pi)J_\nu(z) - \sin(\nu\pi)Y_\nu(z), \quad J_{-n}(z) = (-1)^n J_n(z) \ (n \in \mathbb{Z})$$
///
/// # Arguments
///
/// * `z` - Complex argument in the cut plane $-\pi < \arg(z) \le \pi$.
/// * `order` - Order of the initial J function (any real number).
/// * `scaling` - Scaling option:
///     * [Scaling::Unscaled]: returns $J_\nu(z)$.
///     * [Scaling::Scaled]: returns $\exp(-|\text{Im}(z)|) J_\nu(z)$, which removes exponential
///       growth in both upper and lower half-planes as $|z| \to \infty$.
/// * `n` - Number of members in the sequence ($n \ge 1$).
///
/// # Returns
///
/// A `Result` containing a tuple `(y, n_zeros)`:
/// * `y`: A vector of complex values for orders `[order, order + 1, ..., order + n - 1]`.
/// * `n_zeros`: The number of components at the end of `y` (highest orders) set to zero due to underflow.
pub fn complex_bessel_j<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    n: usize,
) -> BesselResult<T> {
    reflect_orders(z, order, scaling, n, BesselJ)
}

/// Computes the modified Bessel function of the second kind $K_\nu(z)$ for a complex argument.
///
/// Computes a sequence of complex modified Bessel functions `y[j] = K(order + j, z)`
/// for real orders `order + j` (`j = 0, ..., n - 1`) and complex argument `z` ($z \ne 0$)
/// in the cut plane $-\pi < \arg(z) \le \pi$.
///
/// Negative orders are evaluated via the reflection identity (DLMF 10.27.3):
///
/// $$K_{-\nu}(z) = K_\nu(z)$$
///
/// # Arguments
///
/// * `z` - Complex argument $z \ne 0$ in the cut plane $-\pi < \arg(z) \le \pi$.
/// * `order` - Order of the initial K function (any real number).
/// * `scaling` - Scaling option:
///     * [Scaling::Unscaled]: returns $K_\nu(z)$.
///     * [Scaling::Scaled]: returns $\exp(z) K_\nu(z)$, which removes exponential decay
///       as $\text{Re}(z) \to +\infty$.
/// * `n` - Number of members in the sequence ($n \ge 1$).
///
/// # Returns
///
/// A `Result` containing a tuple `(y, n_zeros)`:
/// * `y`: A vector of complex values for orders `[order, order + 1, ..., order + n - 1]`.
/// * `n_zeros`: The number of components at the start of `y` set to zero due to underflow.
pub fn complex_bessel_k<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    n: usize,
) -> BesselResult<T> {
    reflect_orders(z, order, scaling, n, BesselK)
}

/// Computes the Bessel function of the second kind $Y_\nu(z)$ for a complex argument.
///
/// Computes a sequence of complex Bessel functions `y[j] = Y(order + j, z)`
/// for real orders `order + j` (`j = 0, ..., n - 1`) and complex argument `z` ($z \ne 0$)
/// in the cut plane $-\pi < \arg(z) \le \pi$.
///
/// Negative orders are evaluated via the DLMF reflection formulas (DLMF 10.2.3):
///
/// $$Y_{-\nu}(z) = \sin(\nu\pi)J_\nu(z) + \cos(\nu\pi)Y_\nu(z), \quad Y_{-n}(z) = (-1)^n Y_n(z) \ (n \in \mathbb{Z})$$
///
/// # Arguments
///
/// * `z` - Complex argument $z \ne 0$ in the cut plane $-\pi < \arg(z) \le \pi$.
/// * `order` - Order of the initial Y function (any real number).
/// * `scaling` - Scaling option:
///     * [Scaling::Unscaled]: returns $Y_\nu(z)$.
///     * [Scaling::Scaled]: returns $\exp(-|\text{Im}(z)|) Y_\nu(z)$, which removes exponential
///       growth in both upper and lower half-planes as $|z| \to \infty$.
/// * `n` - Number of members in the sequence ($n \ge 1$).
///
/// # Returns
///
/// A `Result` containing a tuple `(y, n_zeros)`:
/// * `y`: A vector of complex values for orders `[order, order + 1, ..., order + n - 1]`.
/// * `n_zeros`: The number of components at the start of `y` set to zero due to underflow.
pub fn complex_bessel_y<T: BesselFloat>(
    z: Complex<T>,
    order: T,
    scaling: Scaling,
    n: usize,
) -> BesselResult<T> {
    reflect_orders(z, order, scaling, n, BesselY)
}
