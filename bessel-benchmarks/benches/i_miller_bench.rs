use amos_bessel_rs::bessel_i;
use bessel_benchmarks::miller_workload_points;
use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};

fn bench_i_miller(c: &mut Criterion) {
    let mut group = c.benchmark_group("Bessel I Miller Regime");

    let points = miller_workload_points();

    group.bench_function(BenchmarkId::new("bessel_i", "miller_grid"), |b| {
        b.iter(|| {
            for &(order, z) in &points {
                let _ = bessel_i(order, z);
            }
        });
    });

    for &(order, z) in points.iter().take(4) {
        let label = format!("order={order}, z={:.1}+{:.1}i", z.re, z.im);
        group.bench_with_input(BenchmarkId::new("single_eval", &label), &(order, z), |b, &(o, z)| {
            b.iter(|| bessel_i(o, z));
        });
    }

    group.finish();
}

criterion_group!(benches, bench_i_miller);
criterion_main!(benches);
