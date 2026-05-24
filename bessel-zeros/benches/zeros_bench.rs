use bessel_zeros::{BesselFunType, bessel_zeros};
use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};

fn bench_zeros(c: &mut Criterion) {
    let mut group = c.benchmark_group("Bessel Zeros Runtime");

    let types = [
        BesselFunType::J,
        BesselFunType::Y,
        BesselFunType::JP,
        BesselFunType::YP,
    ];

    let precisions: Vec<f64> = (4..=16)
        .step_by(4)
        .map(|i| 10.0f64.powi(-(i as i32)))
        .collect();

    let order = 1.0;

    for fun_type in &types {
        let type_name = match fun_type {
            BesselFunType::J => "J",
            BesselFunType::Y => "Y",
            BesselFunType::JP => "JP",
            BesselFunType::YP => "YP",
        };

        for &precision in &precisions {
            for n_zeros in [50, 100] {
                group.bench_with_input(
                    BenchmarkId::new(
                        format!("Type_{}", type_name),
                        format!("prec_{:.0e}, n={}", precision, n_zeros),
                    ),
                    &precision,
                    |b, &p| b.iter(|| bessel_zeros(fun_type, order, n_zeros, p)),
                );
            }
        }
    }
    group.finish();
}

criterion_group!(benches, bench_zeros);
criterion_main!(benches);
