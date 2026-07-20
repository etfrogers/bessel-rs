use bessel_zeros::{BesselFunType, bessel_zeros, fast};
use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};

fn bench_backends(c: &mut Criterion) {
    let mut group = c.benchmark_group("AMOS vs fast backend");

    let types = [
        ("J", BesselFunType::J),
        ("Y", BesselFunType::Y),
        ("JP", BesselFunType::JP),
        ("YP", BesselFunType::YP),
    ];

    let orders: &[i32] = &[0, 1, 5, 20];
    let n_zeros_values = [10, 100];

    for (type_name, fun_type) in types {
        for &order in orders {
            for &n in &n_zeros_values {
                let id = format!("order={order}, n={n}");

                group.bench_with_input(
                    BenchmarkId::new(format!("amos/{type_name}"), &id),
                    &(order, n),
                    |b, &(ord, nz)| b.iter(|| bessel_zeros(fun_type, ord, nz, 1e-14)),
                );

                group.bench_with_input(
                    BenchmarkId::new(format!("fast/{type_name}"), &id),
                    &(order, n),
                    |b, &(ord, nz)| b.iter(|| fast::bessel_zeros(fun_type, ord, nz, 1e-14)),
                );
            }
        }
    }

    group.finish();
}

criterion_group!(benches, bench_backends);
criterion_main!(benches);
