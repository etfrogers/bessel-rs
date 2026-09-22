#[cfg(feature = "into")]
use amos_bessel_rs::amos::{
    complex_bessel_i_into, complex_bessel_j_into, complex_bessel_k_into, complex_bessel_y_into,
};
#[cfg(not(feature = "into"))]
use amos_bessel_rs::amos::{complex_bessel_i, complex_bessel_j, complex_bessel_k, complex_bessel_y};
use amos_bessel_rs::{Scaling, bessel_i, bessel_j, bessel_k, bessel_y};
use bessel_benchmarks::{dense_r_workload_radii, sequence_workload_cases};
use bessel_zeros::{BesselFunType, bessel_zeros};
use criterion::{Criterion, criterion_group, criterion_main};
#[cfg(feature = "into")]
use num::Complex;
use std::sync::LazyLock;

const TYPES: [BesselFunType; 4] = [
    BesselFunType::J,
    BesselFunType::Y,
    BesselFunType::JP,
    BesselFunType::YP,
];

const PRECISIONS: [f64; 2] = [1e-6, 1e-12];
const N_ZEROS: [usize; 2] = [50, 200];
const ORDERS: [f64; 3] = [0.0, 1.0, 75.5];

/// Runs the combined zeros workload across all 4 function types (J, Y, J', Y'),
/// two target precisions (1e-6, 1e-12), and zero counts (50, 200).
#[inline]
pub fn run_zeros_workload() {
    for fun_type in &TYPES {
        for &order in &ORDERS {
            for &n in &N_ZEROS {
                for &prec in &PRECISIONS {
                    #[cfg(bessel_zeros_ref)]
                    let zeros = bessel_zeros(fun_type, order, n, prec);
                    #[cfg(not(bessel_zeros_ref))]
                    let zeros = bessel_zeros(*fun_type, order, n, prec);
                    std::hint::black_box(zeros);
                }
            }
        }
    }
}

static RADII: LazyLock<Vec<f64>> = LazyLock::new(dense_r_workload_radii);
const DENSE_ORDERS: [f64; 4] = [0.0, 1.0, 2.5, 10.0];

/// Runs the dense radial workload across J, Y, I, K for 4 representative orders
/// over 1,000 geometrically spaced radial points spanning the Taylor series,
/// Miller recurrence, and asymptotic expansion domains.
#[inline]
pub fn run_dense_r_workload() {
    let radii = &*RADII;
    for &order in &DENSE_ORDERS {
        for &r in radii {
            let _ = std::hint::black_box(bessel_j(order, r));
            let _ = std::hint::black_box(bessel_y(order, r));
            let _ = std::hint::black_box(bessel_i(order, r));
            let _ = std::hint::black_box(bessel_k(order, r));
        }
    }
}

/// Runs the sequence recurrence workload across J, Y, I, K for sequences
/// of N = 30 and N = 100 orders into pre-allocated buffers across moderate,
/// intermediate, and asymptotic regimes.
#[inline]
pub fn run_sequence_workload() {
    let cases = sequence_workload_cases();

    #[cfg(feature = "into")]
    {
        let mut buf_30 = [Complex::new(0.0, 0.0); 30];
        let mut buf_100 = [Complex::new(0.0, 0.0); 100];

        for &(order, z) in &cases {
            for buf in [&mut buf_30[..], &mut buf_100[..]] {
                let _ = std::hint::black_box(complex_bessel_j_into(z, order, Scaling::Unscaled, buf));
                let _ = std::hint::black_box(complex_bessel_y_into(z, order, Scaling::Unscaled, buf));
                let _ = std::hint::black_box(complex_bessel_i_into(z, order, Scaling::Unscaled, buf));
                let _ = std::hint::black_box(complex_bessel_k_into(z, order, Scaling::Unscaled, buf));
            }
        }
    }

    #[cfg(not(feature = "into"))]
    {
        for &(order, z) in &cases {
            for n in [30, 100] {
                let _ = std::hint::black_box(complex_bessel_j(z, order, Scaling::Unscaled, n));
                let _ = std::hint::black_box(complex_bessel_y(z, order, Scaling::Unscaled, n));
                let _ = std::hint::black_box(complex_bessel_i(z, order, Scaling::Unscaled, n));
                let _ = std::hint::black_box(complex_bessel_k(z, order, Scaling::Unscaled, n));
            }
        }
    }
}

fn bench_workloads(c: &mut Criterion) {
    c.bench_function("zeros_workload", |b| {
        b.iter(run_zeros_workload);
    });

    c.bench_function("dense_r_workload", |b| {
        b.iter(run_dense_r_workload);
    });

    c.bench_function("sequence_workload", |b| {
        b.iter(run_sequence_workload);
    });
}

criterion_group!(benches, bench_workloads);
criterion_main!(benches);
