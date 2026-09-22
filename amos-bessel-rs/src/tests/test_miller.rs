use num::Complex;

#[test]
fn test_bessel_i_miller_regime() {
    // Tests Bessel I in Domain II (Series-normalized Miller algorithm):
    // max_order <= 1.0 and |z| in [2.5, asymptotic_z_limit]
    let test_points = [
        (0.0, Complex::new(4.0, 3.0)),   // |z| = 5.0
        (0.5, Complex::new(5.0, 4.0)),   // |z| ≈ 6.4
        (1.0, Complex::new(8.0, 6.0)),   // |z| = 10.0
        (0.0, Complex::new(10.0, 5.0)),  // |z| ≈ 11.2
        (0.5, Complex::new(12.0, 8.0)),  // |z| ≈ 14.4
    ];

    for (order, z) in test_points {
        let actual = crate::bessel_i(order, z).expect("actual evaluation failed");
        let expected =
            complex_bessel_rs::bessel_i::bessel_i(order, z).expect("reference evaluation failed");

        let rel_diff = (actual - expected).norm() / expected.norm();
        println!("order={order}, z={z:?}: rel_diff={rel_diff:e}");

        assert!(
            rel_diff < 1e-14,
            "Failed at order={order}, z={z:?}: rel_diff={rel_diff:e}"
        );
    }
}
