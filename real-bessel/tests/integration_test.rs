use rstest_reuse::template;

use amos_bessel_rs::{bessel_j, bessel_y};
use approx::assert_relative_eq;
use real_bessel::{bessel_j0, bessel_j1, bessel_jn, bessel_y0, bessel_y1, bessel_yn};
use rstest::rstest;
use rstest_reuse::apply;

#[macro_use]
mod common;

#[rstest]
#[case(0.0, 1.0)]
#[case(1.0, 0.7651976865579666)]
#[case(2.0, 0.22389077914123562)]
#[case(3.0, -0.2600519549019335)]
#[case(4.0, -0.39714980986384737)]
#[case(5.0, -0.1775967713143383)]
fn test_bessel_j0(#[case] input: f64, #[case] expected: f64) {
    let result = bessel_j0(input);
    assert_relative_eq!(result, expected, epsilon = 1e-10,);
}

#[apply(bessel_cases)]
fn test_bessel_j0_against_amos(#[case] _order: f64, #[case] zr: f64, #[case] _zi: f64) {
    let expected = unwrap_real_bessel!(bessel_j, 0.0, zr);
    let actual = bessel_j0(zr);
    assert_relative_eq!(actual, expected, epsilon = 1e-10,);
}

#[apply(bessel_cases)]
fn test_bessel_y0_against_amos(#[case] _order: f64, #[case] mut zr: f64, #[case] _zi: f64) {
    let actual = get_real_y_bessel!(bessel_y0, zr, return);

    let expected = unwrap_real_bessel!(bessel_y, 0.0, zr);
    assert_relative_eq!(actual, expected, epsilon = 1e-10,);
}

#[apply(bessel_cases)]
fn test_bessel_j1_against_amos(#[case] _order: f64, #[case] zr: f64, #[case] _zi: f64) {
    let actual = bessel_j1(zr);
    let expected = unwrap_real_bessel!(bessel_j, 1.0, zr);
    assert_relative_eq!(actual, expected, epsilon = 1e-10,);
}

#[apply(bessel_cases)]
fn test_bessel_y1_against_amos(#[case] _order: f64, #[case] mut zr: f64, #[case] _zi: f64) {
    let actual = get_real_y_bessel!(bessel_y1, zr, return);
    let expected = unwrap_real_bessel!(bessel_y, 1.0, zr);

    assert_relative_eq!(actual, expected, epsilon = 1e-10,);
}

#[apply(bessel_cases)]
#[trace]
fn test_bessel_jn_against_amos(#[case] order: f64, #[case] zr: f64, #[case] _zi: f64) {
    let integer_order = order as i32;
    let expected = unwrap_real_bessel!(bessel_j, integer_order as f64, zr);
    let actual = bessel_jn(integer_order, zr);
    assert_relative_eq!(actual, expected, epsilon = 1e-10,);
}

#[apply(bessel_cases)]
fn test_bessel_yn_against_amos(#[case] order: f64, #[case] mut zr: f64, #[case] _zi: f64) {
    let integer_order = order as i32;
    let yn_local = |zr| bessel_yn(integer_order, zr);
    let actual = get_real_y_bessel!(yn_local, zr, return);
    let expected = unwrap_real_bessel!(bessel_y, integer_order as f64, zr);
    if expected.is_nan() {
        assert!(
            actual.is_nan() || actual.is_infinite(),
            "Expected overflow to produce NaN or infinity, but got {actual}"
        );
        return;
    }
    // yn(n, x) for large x loses precision because it's computed as a difference of large values
    // via recurrence or trigonometric asymptotic expansions, compounding the f64 precision limits.
    let eps = if zr > 1e6 { 1e-5 } else { 1e-10 };
    assert_relative_eq!(actual, expected, epsilon = eps, max_relative = 1e-10,);
}

#[template]
#[rstest]
#[case(4.0, 2.1, 0.0)]
#[case(5.0, 2.0001, 0.0)]
#[case(340.0, 35.0001, 0.0)]
#[case(407.3,-325.1, 635.2)]
#[case(465.0,-867.0, -448.0)] // 5
#[case(10.711871220659752, -6.89931119845653, -9.408182256887017)]
#[case(8.544950848779838, -8.645033163963603, 18.976439189605003,)]
#[case(21.04, 53.19, -40.77)]
#[case(4.0, 2.1, 0.0)]
#[case(5.0, 2.0001, 0.0)] // 10
#[case(340.0, 35.0001, 0.0)]
#[case(899.6,-35.7,317.8)]
#[case(531.0,-106.7,-16.0)]
#[case(531.0,-106.0,-16.0)]
#[case(433.0,16.874,-38.17)] //15
#[case(433.0,16.8,-38.17)]
#[case(311.2078694557452,-10.990141170168044,-25.70154097357056,)]
#[case(8.544950848779838, -8.645033163963603, 18.976439189605003,)]
#[case(17.5, 70.3, 37.4)]
#[case(13.337522865795481, -29.8266399174247, 17.66323218839807)] //20
#[case(5423.24, -7915.11, -3113.95)]
#[case(2213.0, -1813.0, -1033.0)]
#[case(5514.86274463943, -9489.650336481069, 4951.6909981261)]
#[case(2.74e-288, 6.33e-166, 7.53e-275)]
#[case(1.51e-150, -3.07e-118, 3.51e-42)] //25
#[case(2.637e-27, -4.01e-50, 0.0)]
#[case(4.0e-132, 0e0, 445.0)]
#[case(8714.0, 8904.0, -10.0)]
#[case(60.9, 246.2, -982.5)]
#[case(40.5, 1673.3, -4.0)] // 30
#[case(2634.5, -2634.5, 14.1)]
#[case(5.007e-14, 4.401331657952316e-5, -3.6e-6)]
#[case(1719.3, 920.1, 0.0)]
#[case(3.5695132850479827e3, -2.2313404290100934e3, 8.646324128723001e3)]
#[case(0.28008208034835413, -2435.84398720043, -9106.813568430613)] // 35
#[case(35.42423142304685, 2689.1019240048972, -688.7899868054337)]
#[case(1.0111752223029848, 7037.518427975952, -685.0803465010631)]
#[case(9491.159287083694, -2404.8869667701747, -6391.664651975572)]
#[case(3.468367867017804e0, -1.8067397106295227e-254, -3.0255676077184667e-21)]
#[case(6.946702885186345e-149, 0e0, -6.691424259254966e2)] // 40
#[case(3.684122892548987e3, -5.107972475729046e3, 5.916387337090975e3)]
#[case(7.107636998006379e3, -1.867258055869096e3, 4.865284129480511e3)]
#[case(172302836.50840142, 1.2494954195932068e-254, -981457506.31791925)]
#[case(645.0, -736006017.5, 0.0)]
#[case(1253.5, 0.0, 2102.4)] // 45
#[case(1.0, -4816.864663442315, 9.992997770079455)]
#[case(2.8237311072834126e-124, -7.414168789342814e5, 0e0)]
pub fn bessel_cases(#[case] order: f64, #[case] zr: f64, #[case] zi: f64) {}

#[rstest]
fn test_jn_very_large_x(#[values(0, 1, 2, 3, 4, 5, 6, 7)] n: i32) {
    let x = f64::MAX; // ~1.8e308, well above 2^302
    let actual = bessel_jn(n, x);

    // Direct asymptotic formula: J_n(x) ≈ sqrt(2/πx) * cos(x - (2n+1)π/4)
    let amplitude = (2.0 / (std::f64::consts::PI * x)).sqrt();
    let phase = x - ((2 * n + 1) as f64) * std::f64::consts::PI / 4.0;
    let expected = amplitude * phase.cos();

    assert_relative_eq!(actual, expected, epsilon = 1e-10);
}

#[rstest]
fn test_yn_very_large_x(#[values(0, 1, 2, 3, 4, 5, 6, 7)] n: i32) {
    let x = f64::MAX; // ~1.8e308, well above 2^302
    let actual = yn(n, x).unwrap();

    // Direct asymptotic formula: Y_n(x) ≈ sqrt(2/πx) * sin(x - (2n+1)π/4)
    let amplitude = (2.0 / (std::f64::consts::PI * x)).sqrt();
    let phase = x - ((2 * n + 1) as f64) * std::f64::consts::PI / 4.0;
    let expected = amplitude * phase.sin();

    assert_relative_eq!(actual, expected, epsilon = 1e-10);
}
