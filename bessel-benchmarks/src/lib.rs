//! Shared utilities and test point generators for Bessel benchmarks.

use num::Complex;

/// Generates a grid of points (order, z) in the moderate |z| regime (Domain II)
/// that exercise Miller's algorithm (Series Normalization) for I_v(z).
pub fn miller_workload_points() -> Vec<(f64, Complex<f64>)> {
    let orders = [0.0, 0.5, 1.0, 2.0, 5.0];
    let z_values = [
        Complex::new(3.0, 2.0),
        Complex::new(5.0, 4.0),
        Complex::new(8.0, 6.0),
        Complex::new(10.0, 5.0),
        Complex::new(12.0, 8.0),
        Complex::new(2.5, 12.0),
    ];

    let mut points = Vec::with_capacity(orders.len() * z_values.len());
    for &order in &orders {
        for &z in &z_values {
            points.push((order, z));
        }
    }
    points
}

/// Generates a dense grid of 1,000 radial points geometrically spaced
/// across `[0.01, 600.0]`.
///
/// This spans:
/// - Small $r \in [0.01, 2.0]$ (Taylor/power series domain)
/// - Intermediate $r \in [2.0, 25.0]$ (Miller recurrence domain)
/// - Large $r \in [25.0, 600.0]$ (Hankel asymptotic domain)
///
/// $r \le 600.0$ ensures all Bessel families ($J, Y, I, K$) evaluate cleanly
/// without floating-point overflow for unscaled $I_\nu$ ($e^{600} \approx 3.8 \times 10^{260} < 1.8 \times 10^{308}$).
pub fn dense_r_workload_radii() -> Vec<f64> {
    let n = 1000;
    let r_min = 0.01_f64;
    let r_max = 600.0_f64;
    let ratio = (r_max / r_min).powf(1.0 / (n - 1) as f64);
    let mut radii = Vec::with_capacity(n);
    let mut current = r_min;
    for _ in 0..n {
        radii.push(current);
        current *= ratio;
    }
    radii
}

/// Representative (initial_order, z) test cases for sequence recurrence:
/// - Moderate $|z|$: Domain II (Miller recurrence & Wronskians)
/// - Intermediate $|z|$: Domain II transition regime
/// - Large $|z|$: Domain III (Hankel asymptotic expansions)
///
/// Includes both purely real and complex arguments.
pub fn sequence_workload_cases() -> [(f64, Complex<f64>); 4] {
    [
        (0.0, Complex::new(5.0, 0.0)),   // Real radial argument, |z| = 5.0
        (0.0, Complex::new(4.0, 3.0)),   // Complex argument, |z| = 5.0
        (0.5, Complex::new(15.0, 10.0)), // Complex intermediate, |z| ≈ 18.0
        (1.0, Complex::new(60.0, 0.0)),  // Real large argument, |z| = 60.0
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use amos_bessel_rs::{
        Scaling, bessel_i, bessel_j, bessel_k, bessel_y,
    };
    #[cfg(feature = "into")]
    use amos_bessel_rs::amos::{
        complex_bessel_i_into, complex_bessel_j_into, complex_bessel_k_into, complex_bessel_y_into,
    };
    #[cfg(not(feature = "into"))]
    use amos_bessel_rs::amos::{
        complex_bessel_i, complex_bessel_j, complex_bessel_k, complex_bessel_y,
    };

    #[test]
    fn test_dense_r_evaluation() {
        let radii = dense_r_workload_radii();
        assert_eq!(radii.len(), 1000);
        let orders = [0.0, 1.0, 2.5, 10.0];
        for &order in &orders {
            for &r in &radii {
                let j = bessel_j(order, r);
                assert!(j.is_ok(), "bessel_j failed at order {order}, r {r}: {j:?}");
                let y = bessel_y(order, r);
                assert!(y.is_ok(), "bessel_y failed at order {order}, r {r}: {y:?}");
                let i = bessel_i(order, r);
                assert!(i.is_ok(), "bessel_i failed at order {order}, r {r}: {i:?}");
                let k = bessel_k(order, r);
                assert!(k.is_ok(), "bessel_k failed at order {order}, r {r}: {k:?}");
            }
        }
    }

    #[cfg(feature = "into")]
    #[test]
    fn test_sequence_evaluation() {
        let cases = sequence_workload_cases();
        let mut buf_30 = [Complex::new(0.0, 0.0); 30];
        let mut buf_100 = [Complex::new(0.0, 0.0); 100];

        for &(order, z) in &cases {
            for buf in [&mut buf_30[..], &mut buf_100[..]] {
                let j = complex_bessel_j_into(z, order, Scaling::Unscaled, buf);
                assert!(j.is_ok(), "complex_bessel_j_into failed at z={z}, order={order}: {j:?}");
                let y = complex_bessel_y_into(z, order, Scaling::Unscaled, buf);
                assert!(y.is_ok(), "complex_bessel_y_into failed at z={z}, order={order}: {y:?}");
                let i = complex_bessel_i_into(z, order, Scaling::Unscaled, buf);
                assert!(i.is_ok(), "complex_bessel_i_into failed at z={z}, order={order}: {i:?}");
                let k = complex_bessel_k_into(z, order, Scaling::Unscaled, buf);
                assert!(k.is_ok(), "complex_bessel_k_into failed at z={z}, order={order}: {k:?}");
            }
        }
    }

    #[cfg(not(feature = "into"))]
    #[test]
    fn test_sequence_evaluation() {
        let cases = sequence_workload_cases();

        for &(order, z) in &cases {
            for n in [30, 100] {
                let j = complex_bessel_j(z, order, Scaling::Unscaled, n);
                assert!(j.is_ok(), "complex_bessel_j failed at z={z}, order={order}: {j:?}");
                let y = complex_bessel_y(z, order, Scaling::Unscaled, n);
                assert!(y.is_ok(), "complex_bessel_y failed at z={z}, order={order}: {y:?}");
                let i = complex_bessel_i(z, order, Scaling::Unscaled, n);
                assert!(i.is_ok(), "complex_bessel_i failed at z={z}, order={order}: {i:?}");
                let k = complex_bessel_k(z, order, Scaling::Unscaled, n);
                assert!(k.is_ok(), "complex_bessel_k failed at z={z}, order={order}: {k:?}");
            }
        }
    }
}
