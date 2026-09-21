use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicUsize, Ordering};

use amos_bessel_rs::{
    HankelKind, Scaling, airy, airy_b, airy_bp, airyp,
    amos::{
        complex_bessel_i_into, complex_bessel_j_into, complex_bessel_k_into, complex_bessel_y_into,
        complex_hankel1_into, complex_hankel2_into,
    },
    bessel_i, bessel_j, bessel_k, bessel_y,
    derivatives::{
        bessel_i_p, bessel_j_derivative, bessel_j_p, bessel_k_p, bessel_y_p,
    },
    hankel,
};
use num::Complex;

// Tracking allocator to verify zero heap allocations at runtime
struct CountingAlloc;
static ALLOC_COUNT: AtomicUsize = AtomicUsize::new(0);

unsafe impl GlobalAlloc for CountingAlloc {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        ALLOC_COUNT.fetch_add(1, Ordering::SeqCst);
        unsafe { System.alloc(layout) }
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        unsafe { System.dealloc(ptr, layout) }
    }
}

#[global_allocator]
static A: CountingAlloc = CountingAlloc;

fn assert_zero_alloc<R>(f: impl FnOnce() -> R) -> R {
    let before = ALLOC_COUNT.load(Ordering::SeqCst);
    let result = f();
    let after = ALLOC_COUNT.load(Ordering::SeqCst);
    assert_eq!(
        after - before,
        0,
        "Expected zero heap allocations, but {} allocations occurred",
        after - before
    );
    result
}

// =========================================================================
// 1. Single-value convenience APIs (must run with zero heap allocations)
// =========================================================================

fn check_single_value_zero_alloc() {
    assert_zero_alloc(|| {
        // J_0(0.0) = 1.0, J_1(0.0) = 0.0
        let j0: f64 = bessel_j(0.0, 0.0).unwrap();
        assert!((j0 - 1.0).abs() < 1e-14);

        let j1: f64 = bessel_j(1.0, 0.0).unwrap();
        assert!(j1.abs() < 1e-14);

        // I_0(0.0) = 1.0, I_1(0.0) = 0.0
        let i0: f64 = bessel_i(0.0, 0.0).unwrap();
        assert!((i0 - 1.0).abs() < 1e-14);

        // Complex argument
        let z = Complex::new(1.0, 1.0);
        let jz: Complex<f64> = bessel_j(0.0, z).unwrap();
        assert!(jz.re.is_finite() && jz.im.is_finite());

        let yz: Complex<f64> = bessel_y(0.5, z).unwrap();
        assert!(yz.re.is_finite() && yz.im.is_finite());

        let kz: Complex<f64> = bessel_k(1.0, z).unwrap();
        assert!(kz.re.is_finite() && kz.im.is_finite());

        let h1: Complex<f64> = hankel(0.0, z, HankelKind::First).unwrap();
        assert!(h1.re.is_finite() && h1.im.is_finite());

        let h2: Complex<f64> = hankel(0.0, z, HankelKind::Second).unwrap();
        assert!(h2.re.is_finite() && h2.im.is_finite());

        // Airy functions
        let ai: Complex<f64> = airy(z).unwrap();
        let aip: Complex<f64> = airyp(z).unwrap();
        let bi: Complex<f64> = airy_b(z).unwrap();
        let bip: Complex<f64> = airy_bp(z).unwrap();
        assert!(ai.re.is_finite() && aip.re.is_finite() && bi.re.is_finite() && bip.re.is_finite());
    });
}

fn check_negative_order_reflections_zero_alloc() {
    assert_zero_alloc(|| {
        let z = Complex::new(2.0, 1.5);
        // Integer orders: J_{-n}(z) = (-1)^n J_n(z)
        let j_pos: Complex<f64> = bessel_j(2.0, z).unwrap();
        let j_neg: Complex<f64> = bessel_j(-2.0, z).unwrap();
        assert!((j_pos.re - j_neg.re).abs() < 1e-13);
        assert!((j_pos.im - j_neg.im).abs() < 1e-13);

        // Non-integer orders (triggers reflection routine using ScratchBuffer)
        let j_nonint: Complex<f64> = bessel_j(-2.5, z).unwrap();
        assert!(j_nonint.re.is_finite() && j_nonint.im.is_finite());

        let y_neg: Complex<f64> = bessel_y(-1.5, z).unwrap();
        assert!(y_neg.re.is_finite() && y_neg.im.is_finite());

        let i_neg: Complex<f64> = bessel_i(-1.5, z).unwrap();
        assert!(i_neg.re.is_finite() && i_neg.im.is_finite());

        let k_neg: Complex<f64> = bessel_k(-1.5, z).unwrap();
        let k_pos: Complex<f64> = bessel_k(1.5, z).unwrap();
        assert!((k_pos.re - k_neg.re).abs() < 1e-13);
        assert!((k_pos.im - k_neg.im).abs() < 1e-13);
    });
}

// =========================================================================
// 2. In-place sequence APIs (_into) with stack buffers (zero alloc)
// =========================================================================

fn check_into_sequence_stack_buffers_zero_alloc() {
    assert_zero_alloc(|| {
        let z = Complex::new(3.0, 0.5);
        let mut out_j = [Complex::<f64>::ZERO; 5];
        let mut out_y = [Complex::<f64>::ZERO; 5];
        let mut out_i = [Complex::<f64>::ZERO; 5];
        let mut out_k = [Complex::<f64>::ZERO; 5];
        let mut out_h1 = [Complex::<f64>::ZERO; 5];
        let mut out_h2 = [Complex::<f64>::ZERO; 5];

        let info_j = complex_bessel_j_into(z, 0.0, Scaling::Unscaled, &mut out_j).unwrap();
        let info_y = complex_bessel_y_into(z, 0.0, Scaling::Unscaled, &mut out_y).unwrap();
        let info_i = complex_bessel_i_into(z, 0.0, Scaling::Unscaled, &mut out_i).unwrap();
        let info_k = complex_bessel_k_into(z, 0.0, Scaling::Unscaled, &mut out_k).unwrap();
        let info_h1 = complex_hankel1_into(z, 0.0, Scaling::Unscaled, &mut out_h1).unwrap();
        let info_h2 = complex_hankel2_into(z, 0.0, Scaling::Unscaled, &mut out_h2).unwrap();

        assert_eq!(info_j.n_zeros, 0);
        assert_eq!(info_y.n_zeros, 0);
        assert_eq!(info_i.n_zeros, 0);
        assert_eq!(info_k.n_zeros, 0);
        assert_eq!(info_h1.n_zeros, 0);
        assert_eq!(info_h2.n_zeros, 0);

        // Verify Hankel connection: H^(1)_0(z) = J_0(z) + i * Y_0(z)
        let expected_h1 = out_j[0] + Complex::new(0.0, 1.0) * out_y[0];
        assert!((out_h1[0].re - expected_h1.re).abs() < 1e-13);
        assert!((out_h1[0].im - expected_h1.im).abs() < 1e-13);

        // Verify recurrence: J_{n-1}(z) + J_{n+1}(z) = (2n / z) * J_n(z) for n = 1, 2, 3
        for n in 1..=3 {
            let left = out_j[n - 1] + out_j[n + 1];
            let right = (2.0 * n as f64 / z) * out_j[n];
            assert!((left.re - right.re).abs() < 1e-12);
            assert!((left.im - right.im).abs() < 1e-12);
        }
    });
}

// =========================================================================
// 3. Derivatives with internal SBO stack buffers (zero alloc)
// =========================================================================

fn check_derivatives_zero_alloc() {
    assert_zero_alloc(|| {
        let z = Complex::new(2.5, 1.0);
        // Identity: J_0'(z) = -J_1(z)
        let dj0: Complex<f64> = bessel_j_p(0.0, z).unwrap();
        let j1: Complex<f64> = bessel_j(1.0, z).unwrap();
        assert!((dj0.re + j1.re).abs() < 1e-13);
        assert!((dj0.im + j1.im).abs() < 1e-13);

        // Identity: I_0'(z) = I_1(z)
        let di0: Complex<f64> = bessel_i_p(0.0, z).unwrap();
        let i1: Complex<f64> = bessel_i(1.0, z).unwrap();
        assert!((di0.re - i1.re).abs() < 1e-13);
        assert!((di0.im - i1.im).abs() < 1e-13);

        // Higher-order derivative
        let d2j: Complex<f64> = bessel_j_derivative(0.0, z, 2, Scaling::Unscaled).unwrap();
        assert!(d2j.re.is_finite() && d2j.im.is_finite());

        // Other derivatives execute without allocating
        let dy: Complex<f64> = bessel_y_p(1.5, z).unwrap();
        let dk: Complex<f64> = bessel_k_p(1.5, z).unwrap();
        assert!(dy.re.is_finite() && dk.re.is_finite());
    });
}

// =========================================================================
// 4. SBO Capacity Boundary: N <= 32 (Stack) vs N > 32 (Heap fallback / Error)
// =========================================================================

fn check_sbo_capacity_within_limit_succeeds_zero_alloc() {
    assert_zero_alloc(|| {
        let z = Complex::new(5.0, 2.0);
        // N = 32: all 32 orders negative (-32.5 to -1.5) triggers ScratchBuffer::new(32)
        let mut out = [Complex::<f64>::ZERO; 32];
        let res = complex_bessel_j_into(z, -32.5, Scaling::Unscaled, &mut out);
        assert!(res.is_ok(), "N=32 negative orders should succeed on stack without error");
    });
}

#[cfg(not(feature = "alloc"))]
fn check_sbo_capacity_exceeded_without_alloc_returns_invalid_input() {
    use amos_bessel_rs::BesselError;
    let z = Complex::new(5.0, 2.0);
    // N = 33: all 33 orders negative (-33.5 to -1.5) triggers ScratchBuffer::new(33), exceeding CAP=32
    let mut out = [Complex::<f64>::ZERO; 33];
    let res = complex_bessel_j_into(z, -33.5, Scaling::Unscaled, &mut out);
    match res {
        Err(BesselError::InvalidInput { .. }) => {
            // Expected: pure core cannot spill to heap, returns InvalidInput error
        }
        other => panic!("Expected Err(BesselError::InvalidInput) when exceeding stack capacity without alloc, got: {other:?}"),
    }
}

#[cfg(feature = "alloc")]
fn check_sbo_capacity_exceeded_with_alloc_spills_to_heap_and_succeeds() {
    let z = Complex::new(5.0, 2.0);
    // N = 33: exceeds CAP=32, so with feature = "alloc" it should allocate on heap and succeed
    let mut out = [Complex::<f64>::ZERO; 33];
    let before = ALLOC_COUNT.load(Ordering::SeqCst);
    let res = complex_bessel_j_into(z, -33.5, Scaling::Unscaled, &mut out);
    let after = ALLOC_COUNT.load(Ordering::SeqCst);

    assert!(res.is_ok(), "N=33 with alloc should succeed by spilling to heap");
    assert!(
        after > before,
        "Expected heap allocation when spilling beyond CAP=32, but no allocation detected"
    );
}

// =========================================================================
// 5. Alloc-only allocating API tests (under feature = "alloc")
// =========================================================================

#[cfg(feature = "alloc")]
fn check_allocating_apis_work_when_alloc_enabled() {
    use amos_bessel_rs::amos::complex_bessel_j;

    let z = Complex::new(2.0, 1.0);
    let before = ALLOC_COUNT.load(Ordering::SeqCst);
    let (vec_res, info) = complex_bessel_j(z, 1.0, Scaling::Unscaled, 10).unwrap();
    let after = ALLOC_COUNT.load(Ordering::SeqCst);

    assert_eq!(vec_res.len(), 10);
    assert_eq!(info.n_zeros, 0);
    assert!(after > before, "complex_bessel_j should allocate a Vec on the heap");
}

// =========================================================================
// Main Test Entry Point
// =========================================================================

#[test]
fn test_no_std_compliance_suite() {
    check_single_value_zero_alloc();
    check_negative_order_reflections_zero_alloc();
    check_into_sequence_stack_buffers_zero_alloc();
    check_derivatives_zero_alloc();
    check_sbo_capacity_within_limit_succeeds_zero_alloc();

    #[cfg(not(feature = "alloc"))]
    check_sbo_capacity_exceeded_without_alloc_returns_invalid_input();

    #[cfg(feature = "alloc")]
    check_sbo_capacity_exceeded_with_alloc_spills_to_heap_and_succeeds();

    #[cfg(feature = "alloc")]
    check_allocating_apis_work_when_alloc_enabled();
}
