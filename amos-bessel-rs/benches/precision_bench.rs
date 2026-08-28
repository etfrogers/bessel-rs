use amos_bessel_rs::Scaling;
use amos_bessel_rs::amos::{
    complex_airy, complex_airy_b, complex_bessel_i, complex_bessel_j, complex_bessel_k,
    complex_bessel_y, complex_hankel1, complex_hankel2,
};
use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use num::Complex;

const CASES: [(f64, f64, f64); 46] = [
    (4.0, 2.1, 0.0),
    (5.0, 2.0001, 0.0),
    (340.0, 35.0001, 0.0),
    (407.3, -325.1, 635.2),
    (465.0, -867.0, -448.0),
    (10.711871220659752, -6.89931119845653, -9.408182256887017),
    (8.544950848779838, -8.645033163963603, 18.976439189605003),
    (21.04, 53.19, -40.77),
    (4.0, 2.1, 0.0),
    (5.0, 2.0001, 0.0),
    (340.0, 35.0001, 0.0),
    (899.6, -35.7, 317.8),
    (531.0, -106.7, -16.0),
    (531.0, -106.0, -16.0),
    (433.0, 16.874, -38.17),
    (433.0, 16.8, -38.17),
    (311.2078694557452, -10.990141170168044, -25.70154097357056),
    (8.544950848779838, -8.645033163963603, 18.976439189605003),
    (17.5, 70.3, 37.4),
    (13.337522865795481, -29.8266399174247, 17.66323218839807),
    (5423.24, -7915.11, -3113.95),
    (2213.0, -1813.0, -1033.0),
    (5514.86274463943, -9489.650336481069, 4951.6909981261),
    (2.74e-288, 6.33e-166, 7.53e-275),
    (1.51e-150, -3.07e-118, 3.51e-42),
    (2.637e-27, -4.01e-50, 0.0),
    (4.0e-132, 0.0, 445.0),
    (8714.0, 8904.0, -10.0),
    (60.9, 246.2, -982.5),
    (40.5, 1673.3, -4.0),
    (2634.5, -2634.5, 14.1),
    (5.007e-14, 4.401331657952316e-5, -3.6e-6),
    (1719.3, 920.1, 0.0),
    (
        3.5695132850479827e3,
        -2.2313404290100934e3,
        8.646324128723001e3,
    ),
    (0.28008208034835413, -2435.84398720043, -9106.813568430613),
    (35.42423142304685, 2689.1019240048972, -688.7899868054337),
    (1.0111752223029848, 7037.518427975952, -685.0803465010631),
    (9491.159287083694, -2404.8869667701747, -6391.664651975572),
    (
        3.468367867017804e0,
        -1.8067397106295227e-254,
        -3.0255676077184667e-21,
    ),
    (6.946702885186345e-149, 0.0, -6.691424259254966e2),
    (
        3.684122892548987e3,
        -5.107972475729046e3,
        5.916387337090975e3,
    ),
    (
        7.107636998006379e3,
        -1.867258055869096e3,
        4.865284129480511e3,
    ),
    (
        172302836.50840142,
        1.2494954195932068e-254,
        -981457506.3179193,
    ),
    (645.0, -736006017.5, 0.0),
    (1253.5, 0.0, 2102.4),
    (1.0, -4816.864663442315, 9.992997770079455),
];

const ORDERS: [f64; 21] = [
    0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 5.0, 10.0, 25.0, 50.0, 75.0, 85.0, 90.0, 100.0, 150.0, 200.0,
    500.0, 1000.0, -0.5, -1.5, -2.0,
];

const Z_PARTS: [f64; 37] = [
    -50.0, -40.0, -30.0, -25.0, -20.0, -15.0, -12.0, -10.0, -8.0, -6.0, -4.0, -3.0, -2.0, -1.0,
    -0.5, -0.1, -0.001, -1e-6, 0.0, 1e-6, 0.001, 0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0,
    12.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0,
];

fn precision_bench(c: &mut Criterion) {
    let mut group = c.benchmark_group("f32 vs f64 Performance");

    let scattered_f64: Vec<(f64, Complex<f64>)> = CASES
        .iter()
        .map(|(o, re, im)| (*o, Complex::new(*re, *im)))
        .collect();

    let scattered_f32: Vec<(f32, Complex<f32>)> = CASES
        .iter()
        .map(|(o, re, im)| (*o as f32, Complex::new(*re as f32, *im as f32)))
        .collect();

    // To prevent the grid from taking forever in benchmarks, we sample a subset or just run the whole thing.
    // 21 * 37 * 37 = 28,749 cases. This is very fast in a compiled bench, but we might want to scale it down slightly
    // for Criterion iterations. We will use the full grid.
    let mut grid_f64 = Vec::with_capacity(ORDERS.len() * Z_PARTS.len() * Z_PARTS.len());
    let mut grid_f32 = Vec::with_capacity(grid_f64.capacity());
    for &o in &ORDERS {
        for &re in &Z_PARTS {
            for &im in &Z_PARTS {
                grid_f64.push((o, Complex::new(re, im)));
                grid_f32.push((o as f32, Complex::new(re as f32, im as f32)));
            }
        }
    }

    // We only benchmark with Unscaled scaling and n=1 to measure the core function speed

    // --- Bessel J ---
    group.bench_function(BenchmarkId::new("J (f64)", "Scattered"), |b| {
        b.iter(|| {
            scattered_f64.iter().for_each(|(o, z)| {
                let _ = complex_bessel_j(*z, *o, Scaling::Unscaled, 1);
            })
        })
    });
    group.bench_function(BenchmarkId::new("J (f32)", "Scattered"), |b| {
        b.iter(|| {
            scattered_f32.iter().for_each(|(o, z)| {
                let _ = complex_bessel_j(*z, *o, Scaling::Unscaled, 1);
            })
        })
    });
    group.bench_function(BenchmarkId::new("J (f64)", "Grid"), |b| {
        b.iter(|| {
            grid_f64.iter().for_each(|(o, z)| {
                let _ = complex_bessel_j(*z, *o, Scaling::Unscaled, 1);
            })
        })
    });
    group.bench_function(BenchmarkId::new("J (f32)", "Grid"), |b| {
        b.iter(|| {
            grid_f32.iter().for_each(|(o, z)| {
                let _ = complex_bessel_j(*z, *o, Scaling::Unscaled, 1);
            })
        })
    });

    // --- Bessel Y ---
    group.bench_function(BenchmarkId::new("Y (f64)", "Scattered"), |b| {
        b.iter(|| {
            scattered_f64.iter().for_each(|(o, z)| {
                let _ = complex_bessel_y(*z, *o, Scaling::Unscaled, 1);
            })
        })
    });
    group.bench_function(BenchmarkId::new("Y (f32)", "Scattered"), |b| {
        b.iter(|| {
            scattered_f32.iter().for_each(|(o, z)| {
                let _ = complex_bessel_y(*z, *o, Scaling::Unscaled, 1);
            })
        })
    });
    group.bench_function(BenchmarkId::new("Y (f64)", "Grid"), |b| {
        b.iter(|| {
            grid_f64.iter().for_each(|(o, z)| {
                let _ = complex_bessel_y(*z, *o, Scaling::Unscaled, 1);
            })
        })
    });
    group.bench_function(BenchmarkId::new("Y (f32)", "Grid"), |b| {
        b.iter(|| {
            grid_f32.iter().for_each(|(o, z)| {
                let _ = complex_bessel_y(*z, *o, Scaling::Unscaled, 1);
            })
        })
    });

    // --- Bessel I ---
    group.bench_function(BenchmarkId::new("I (f64)", "Scattered"), |b| {
        b.iter(|| {
            scattered_f64.iter().for_each(|(o, z)| {
                let _ = complex_bessel_i(*z, *o, Scaling::Unscaled, 1);
            })
        })
    });
    group.bench_function(BenchmarkId::new("I (f32)", "Scattered"), |b| {
        b.iter(|| {
            scattered_f32.iter().for_each(|(o, z)| {
                let _ = complex_bessel_i(*z, *o, Scaling::Unscaled, 1);
            })
        })
    });
    group.bench_function(BenchmarkId::new("I (f64)", "Grid"), |b| {
        b.iter(|| {
            grid_f64.iter().for_each(|(o, z)| {
                let _ = complex_bessel_i(*z, *o, Scaling::Unscaled, 1);
            })
        })
    });
    group.bench_function(BenchmarkId::new("I (f32)", "Grid"), |b| {
        b.iter(|| {
            grid_f32.iter().for_each(|(o, z)| {
                let _ = complex_bessel_i(*z, *o, Scaling::Unscaled, 1);
            })
        })
    });

    // --- Bessel K ---
    group.bench_function(BenchmarkId::new("K (f64)", "Scattered"), |b| {
        b.iter(|| {
            scattered_f64.iter().for_each(|(o, z)| {
                let _ = complex_bessel_k(*z, *o, Scaling::Unscaled, 1);
            })
        })
    });
    group.bench_function(BenchmarkId::new("K (f32)", "Scattered"), |b| {
        b.iter(|| {
            scattered_f32.iter().for_each(|(o, z)| {
                let _ = complex_bessel_k(*z, *o, Scaling::Unscaled, 1);
            })
        })
    });
    group.bench_function(BenchmarkId::new("K (f64)", "Grid"), |b| {
        b.iter(|| {
            grid_f64.iter().for_each(|(o, z)| {
                let _ = complex_bessel_k(*z, *o, Scaling::Unscaled, 1);
            })
        })
    });
    group.bench_function(BenchmarkId::new("K (f32)", "Grid"), |b| {
        b.iter(|| {
            grid_f32.iter().for_each(|(o, z)| {
                let _ = complex_bessel_k(*z, *o, Scaling::Unscaled, 1);
            })
        })
    });

    // --- Hankel H1 ---
    group.bench_function(BenchmarkId::new("H1 (f64)", "Scattered"), |b| {
        b.iter(|| {
            scattered_f64.iter().for_each(|(o, z)| {
                let _ = complex_hankel1(*z, *o, Scaling::Unscaled, 1);
            })
        })
    });
    group.bench_function(BenchmarkId::new("H1 (f32)", "Scattered"), |b| {
        b.iter(|| {
            scattered_f32.iter().for_each(|(o, z)| {
                let _ = complex_hankel1(*z, *o, Scaling::Unscaled, 1);
            })
        })
    });
    group.bench_function(BenchmarkId::new("H1 (f64)", "Grid"), |b| {
        b.iter(|| {
            grid_f64.iter().for_each(|(o, z)| {
                let _ = complex_hankel1(*z, *o, Scaling::Unscaled, 1);
            })
        })
    });
    group.bench_function(BenchmarkId::new("H1 (f32)", "Grid"), |b| {
        b.iter(|| {
            grid_f32.iter().for_each(|(o, z)| {
                let _ = complex_hankel1(*z, *o, Scaling::Unscaled, 1);
            })
        })
    });

    // --- Hankel H2 ---
    group.bench_function(BenchmarkId::new("H2 (f64)", "Scattered"), |b| {
        b.iter(|| {
            scattered_f64.iter().for_each(|(o, z)| {
                let _ = complex_hankel2(*z, *o, Scaling::Unscaled, 1);
            })
        })
    });
    group.bench_function(BenchmarkId::new("H2 (f32)", "Scattered"), |b| {
        b.iter(|| {
            scattered_f32.iter().for_each(|(o, z)| {
                let _ = complex_hankel2(*z, *o, Scaling::Unscaled, 1);
            })
        })
    });
    group.bench_function(BenchmarkId::new("H2 (f64)", "Grid"), |b| {
        b.iter(|| {
            grid_f64.iter().for_each(|(o, z)| {
                let _ = complex_hankel2(*z, *o, Scaling::Unscaled, 1);
            })
        })
    });
    group.bench_function(BenchmarkId::new("H2 (f32)", "Grid"), |b| {
        b.iter(|| {
            grid_f32.iter().for_each(|(o, z)| {
                let _ = complex_hankel2(*z, *o, Scaling::Unscaled, 1);
            })
        })
    });

    // --- Airy Ai (Return Derivative = False) ---
    // Note: Airy functions do not take an 'order' argument. We use the scattered and grid Z points directly.
    let z_scattered_f64: Vec<Complex<f64>> = scattered_f64.iter().map(|(_, z)| *z).collect();
    let z_scattered_f32: Vec<Complex<f32>> = scattered_f32.iter().map(|(_, z)| *z).collect();
    // For grid, we just need Z_PARTS x Z_PARTS (37 x 37 = 1369 cases)
    let mut z_grid_f64 = Vec::with_capacity(Z_PARTS.len() * Z_PARTS.len());
    let mut z_grid_f32 = Vec::with_capacity(Z_PARTS.len() * Z_PARTS.len());
    for &re in &Z_PARTS {
        for &im in &Z_PARTS {
            z_grid_f64.push(Complex::new(re, im));
            z_grid_f32.push(Complex::new(re as f32, im as f32));
        }
    }

    group.bench_function(BenchmarkId::new("Ai (f64)", "Scattered"), |b| {
        b.iter(|| {
            z_scattered_f64.iter().for_each(|z| {
                let _ = complex_airy(*z, false, Scaling::Unscaled);
            })
        })
    });
    group.bench_function(BenchmarkId::new("Ai (f32)", "Scattered"), |b| {
        b.iter(|| {
            z_scattered_f32.iter().for_each(|z| {
                let _ = complex_airy(*z, false, Scaling::Unscaled);
            })
        })
    });
    group.bench_function(BenchmarkId::new("Ai (f64)", "Grid"), |b| {
        b.iter(|| {
            z_grid_f64.iter().for_each(|z| {
                let _ = complex_airy(*z, false, Scaling::Unscaled);
            })
        })
    });
    group.bench_function(BenchmarkId::new("Ai (f32)", "Grid"), |b| {
        b.iter(|| {
            z_grid_f32.iter().for_each(|z| {
                let _ = complex_airy(*z, false, Scaling::Unscaled);
            })
        })
    });

    // --- Airy Bi (Return Derivative = False) ---
    group.bench_function(BenchmarkId::new("Bi (f64)", "Scattered"), |b| {
        b.iter(|| {
            z_scattered_f64.iter().for_each(|z| {
                let _ = complex_airy_b(*z, false, Scaling::Unscaled);
            })
        })
    });
    group.bench_function(BenchmarkId::new("Bi (f32)", "Scattered"), |b| {
        b.iter(|| {
            z_scattered_f32.iter().for_each(|z| {
                let _ = complex_airy_b(*z, false, Scaling::Unscaled);
            })
        })
    });
    group.bench_function(BenchmarkId::new("Bi (f64)", "Grid"), |b| {
        b.iter(|| {
            z_grid_f64.iter().for_each(|z| {
                let _ = complex_airy_b(*z, false, Scaling::Unscaled);
            })
        })
    });
    group.bench_function(BenchmarkId::new("Bi (f32)", "Grid"), |b| {
        b.iter(|| {
            z_grid_f32.iter().for_each(|z| {
                let _ = complex_airy_b(*z, false, Scaling::Unscaled);
            })
        })
    });
}

criterion_group!(benches, precision_bench);
criterion_main!(benches);
