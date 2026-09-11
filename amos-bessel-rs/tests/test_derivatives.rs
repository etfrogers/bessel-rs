use amos_bessel_rs::{BesselError, Scaling, derivatives::bessel_j_derivative};
use approx::assert_relative_eq;
use num::Complex;
use rstest::rstest;
mod common;
use amos_bessel_rs::amos::complex_bessel_j;
use common::{BesselSig, ORDERS, Z_PARTS};

#[rstest]
#[case(0, -0.44005)]
#[case(1, 0.32515)]
#[case(2, 0.21033)]
fn test_bessel_j_hardcoded(#[case] order: u32, #[case] expected: f64) {
    let dz = bessel_j_derivative(order, 1.0, 1).unwrap();
    assert_relative_eq!(dz, expected, epsilon = 1e-4);
}

#[rstest]
fn test_against_numerical_derivative(
    #[values(1, 2, 3, 4, 5)] derivative_order: u32,
    // #[values(Scaling::Unscaled, Scaling::Scaled)] scaling: Scaling,
) {
    for order in ORDERS {
        for zr in Z_PARTS {
            for zi in Z_PARTS {
                let z = Complex::new(zr, zi);

                if let Ok(dz) = bessel_j_derivative(order, z, derivative_order)
                    && let Ok(dz_num) =
                        calculate_numerical_derivative(complex_bessel_j, derivative_order, order, z)
                {
                    assert_relative_eq!(dz, dz_num, max_relative = 1e-7);
                }
            }
        }
    }
}

fn calculate_numerical_derivative(
    func: BesselSig,
    derivative_order: u32,
    order: f64,
    z: Complex<f64>,
) -> Result<Complex<f64>, BesselError> {
    const EPSILON: f64 = 1e-6;
    let n_points = (derivative_order + 1) as i32;
    let mut points = Vec::with_capacity(n_points as usize);
    for i in 0..n_points {
        let offset = EPSILON * ((i - (derivative_order as i32) / 2) as f64);
        let z_loop = z + offset;
        let result = func(z_loop, order, Scaling::Unscaled, 1)?;
        points.push(result.0[0]);
    }
    let mut ans = points;
    for _ in 0..derivative_order {
        ans = diff(&ans, EPSILON);
    }
    Ok(ans[0])
}

// fn calculate_numerical_derivative(order: f64, z: Complex<f64>) -> Complex<f64> {
//     const EPSILON: f64 = 1e-6;
//     let n_points = 10;
//     letpoints, EPSILON);
//     }
// }

fn diff(points: &[Complex<f64>], eps: f64) -> Vec<Complex<f64>> {
    let mut result = Vec::with_capacity(points.len() - 1);
    for i in 0..points.len() - 1 {
        let diff = (points[i + 1] - points[i]) / eps;
        result.push(diff);
    }
    result
}
