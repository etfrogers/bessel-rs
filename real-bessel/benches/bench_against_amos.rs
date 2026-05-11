use amos_bessel_rs::bessel_j;
use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use itertools::Itertools;
use real_bessel::{bessel_j0, bessel_j1, bessel_jn};

const ORDERS: [i32; 15] = [0, 1, 2, 5, 10, 25, 50, 75, 85, 90, 100, 150, 200, 500, 1000];

const Z_PARTS: [f64; 37] = [
    -50.0, -40.0, -30.0, -25.0, -20.0, -15.0, -12.0, -10.0, -8.0, -6.0, -4.0, -3.0, -2.0, -1.0,
    -0.5, -0.1, -0.001, -1e-6, 0.0, 1e-6, 0.001, 0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0,
    12.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0,
];

fn bessel_j0_bench(c: &mut Criterion) {
    let mut group = c.benchmark_group("Bessel functions");
    group.bench_with_input(
        BenchmarkId::new("Bessel J0 real", "Fixed cases"),
        &Z_PARTS,
        |b, z_parts| {
            b.iter(|| {
                z_parts.iter().for_each(|z| {
                    let _ = bessel_j0(*z);
                })
            })
        },
    );

    group.bench_with_input(
        BenchmarkId::new("Bessel J0 complex", "Fixed cases"),
        &Z_PARTS,
        |b, z_parts| {
            b.iter(|| {
                z_parts.iter().for_each(|z| {
                    let _ = bessel_j(0.0, *z);
                })
            })
        },
    );

    group.bench_with_input(
        BenchmarkId::new("Bessel J1 real", "Fixed cases"),
        &Z_PARTS,
        |b, z_parts| {
            b.iter(|| {
                z_parts.iter().for_each(|z| {
                    let _ = bessel_j1(*z);
                })
            })
        },
    );

    group.bench_with_input(
        BenchmarkId::new("Bessel J1 complex", "Fixed cases"),
        &Z_PARTS,
        |b, z_parts| {
            b.iter(|| {
                z_parts.iter().for_each(|z| {
                    let _ = bessel_j(1.0, *z);
                })
            })
        },
    );

    let cases: Vec<_> = ORDERS.into_iter().cartesian_product(Z_PARTS).collect();

    group.bench_with_input(
        BenchmarkId::new("Bessel Jn real", "Fixed cases"),
        &cases,
        |b, cases| {
            b.iter(|| {
                cases.iter().for_each(|(order, z)| {
                    let _ = bessel_jn(*order, *z);
                })
            })
        },
    );

    group.bench_with_input(
        BenchmarkId::new("Bessel Jn complex", "Fixed cases"),
        &cases,
        |b, cases| {
            b.iter(|| {
                cases.iter().for_each(|(order, z)| {
                    let _ = bessel_j(*order, *z);
                })
            })
        },
    );
}

criterion_group!(benches, bessel_j0_bench);
criterion_main!(benches);
