//! Core implementation of the Temme algorithm for finding Bessel zeros.
//!
//! The algorithm is generic over [`BesselBackend`] so that it can be used with
//! either the AMOS or the real-bessel backend without any code duplication.
//!
//! Reference:
//! N. M. Temme, "An Algorithm with ALGOL 60 Program for the Computation of the
//! zeros of the Ordinary Bessel Functions and those of their Derivatives",
//! Journal of Computational Physics, 32, 270-279 (1979).

use std::f64::consts::{FRAC_PI_2, PI};

use crate::BesselFunType;
use crate::backend::BesselBackend;

/// Finds the first `n_zeros` zeros of the specified Bessel function using
/// backend `B` for function evaluations inside the Newton iteration.
pub(crate) fn bessel_zeros_impl<B: BesselBackend>(
    func_type: &BesselFunType,
    order: f64,
    n_zeros: usize,
    precision: f64,
) -> Vec<f64> {
    let order_int = order as i64;
    let mut z = vec![0.0; n_zeros];

    let aa = order.powf(2.0);
    let mu = 4.0 * aa;
    let mu2 = mu.powf(2.0);
    let mu3 = mu.powf(3.0);
    let mu4 = mu.powf(4.0);

    let mut p: f64;
    let p0: f64;
    let p1: f64;
    let q1: f64;
    if func_type.is_non_derivative() {
        p = 7.0 * mu - 31.0;
        p0 = mu - 1.0;

        if (1.0 + p) == p {
            p1 = 0.0;
            q1 = 0.0;
        } else {
            p1 = 4.0 * (253.0 * mu2 - 3722.0 * mu + 17869.0) * p0 / (15.0 * p);
            q1 = 1.6 * (83.0 * mu2 - 982.0 * mu + 3779.0) / p;
        }
    } else {
        p = 7.0 * mu2 + 82.0 * mu - 9.0;
        p0 = mu + 3.0;
        if (p + 1.0) == 1.0 {
            p1 = 0.0;
            q1 = 0.0;
        } else {
            p1 = (4048.0 * mu4 + 131264.0 * mu3 - 221984.0 * mu2 - 417600.0 * mu + 1012176.0)
                / (60.0 * p);
            q1 = 1.6 * (83.0 * mu3 + 2075.0 * mu2 - 3039.0 * mu + 3537.0) / p;
        }
    }

    let t = if (*func_type == BesselFunType::J) || (*func_type == BesselFunType::YP) {
        0.25
    } else {
        0.75
    };

    let tt = 4.0 * t;

    let (pp1, qq1) = if func_type.is_non_derivative() {
        (5. / 48., -5. / 36.)
    } else {
        (-7. / 48., 35. / 288.)
    };

    let y = 0.375 * PI;
    let bb = if order >= 3.0 {
        order.powf(-2.0 / 3.0)
    } else {
        1.0
    };

    let a1 = 3 * order_int - 8;

    for s in 1..=n_zeros {
        let sf = s as f64;
        let mut x: f64;
        let mut w: f64 = 0.0;
        if (order_int == 0) && (s == 1) && (*func_type == BesselFunType::JP) {
            x = 0.0;
        } else {
            if (s as i64) >= a1 {
                let b = (sf + 0.5 * order - t) * PI;
                let c = 0.015625 / (b.powf(2.0));
                x = b - 0.125 * (p0 - p1 * c) / (b * (1.0 - q1 * c));
            } else {
                if s == 1 {
                    x = match func_type {
                        BesselFunType::J => -2.33811,
                        BesselFunType::Y => -1.17371,
                        BesselFunType::JP => -1.01879,
                        BesselFunType::YP => -2.29444,
                    };
                } else {
                    x = y * (4.0 * sf - tt);
                    let v = x.powf(-2.0);
                    x = -(x.powf(2.0 / 3.0)) * (1.0 + v * (pp1 + qq1 * v));
                }
                let u = x * bb;
                let v = fi(2.0 / 3.0 * (-u).powf(1.5));
                w = 1.0 / v.cos();
                let xx = 1.0 - w.powf(2.0);
                let c = (u / xx).sqrt();
                x = if func_type.is_non_derivative() {
                    w * (order + c * (-5.0 / u - c * (6.0 - 10.0 / xx)) / (48.0 * order * u))
                } else {
                    w * (order + c * (7.0 / u + c * (18.0 - 14.0 / xx)) / (48.0 * order * u))
                }
            }

            // Newton refinement: Temme recommends capping at 5 iterations;
            // the initial guess is close enough that convergence is typically
            // reached in 1–2 steps.
            let mut j = 0;
            while (j == 0) || ((j < 5) && ((w / x).abs() > precision)) {
                let xx = x.powf(2.0);
                let x4 = x.powf(4.0);
                let a2 = aa - xx;
                let r0 = bessr::<B>(func_type, order, x);
                j += 1;
                let q: f64;
                let u: f64;
                if func_type.is_non_derivative() {
                    u = r0;
                    w = 6.0 * x * (2.0 * order + 1.0);
                    p = (1.0 - 4.0 * a2) / w;
                    q = (4.0 * (xx - mu) - 2.0 - 12.0 * order) / w;
                } else {
                    u = -xx * r0 / a2;
                    let v = 2.0 * x * a2 / (3.0 * (aa + xx));
                    w = 64.0 * a2.powf(3.0);
                    q = 2.0 * v * (1.0 + mu2 + 32.0 * mu * xx + 48.0 * x4) / w;
                    p = v * (1.0 + (40.0 * mu * xx + 48.0 * x4 - mu2) / w);
                }
                w = u * (1.0 + p * r0) / (1.0 + q * r0);
                x += w;
            }
        }
        z[s - 1] = x;
    }
    z
}

/// Inverse of the function `f(phi) = phi - sin(phi)*cos(phi)`.
/// Used for Airy-function initial guesses in the Temme algorithm.
fn fi(y: f64) -> f64 {
    let c1 = FRAC_PI_2;
    if y == 0.0 {
        0.0
    } else if y > 1e5 {
        c1
    } else {
        let mut p: f64;
        if y < 1.0 {
            p = (3.0 * y).powf(1.0 / 3.0);
            let pp = p.powf(2.0);
            p *= 1.0 + pp * (pp * (27.0 - 2.0 * pp) - 210.0) / 1575.0;
        } else {
            p = 1.0 / (y + c1);
            let pp = p.powf(2.0);
            p = c1
                - p * (1.0
                    + pp * (2310.0
                        + pp * (3003.0 + pp * (4818.0 + pp * (8591.0 + pp * 16328.0))))
                        / 3465.0);
        }
        let pp = (y + p).powf(2.0);
        let r = (p - (p + y).atan()) / pp;
        p - (1.0 + pp) * r * (1.0 + r / (p + y))
    }
}

/// Evaluates the ratio used in the Newton step, dispatching to backend `B`.
fn bessr<B: BesselBackend>(fun_type: &BesselFunType, order: f64, z: f64) -> f64 {
    match fun_type {
        BesselFunType::J => B::j(order, z) / B::j(order + 1.0, z),
        BesselFunType::Y => B::y(order, z) / B::y(order + 1.0, z),
        BesselFunType::JP => order / z - B::j(order + 1.0, z) / B::j(order, z),
        BesselFunType::YP => order / z - B::y(order + 1.0, z) / B::y(order, z),
    }
}
