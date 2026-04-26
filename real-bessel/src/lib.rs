// mod types;
// pub use types::BesselError;

mod j0;
mod j1;
mod jn;
pub use j0::j0 as bessel_j0;
pub use j0::y0 as bessel_y0;
pub use j1::j1 as bessel_j1;
pub use j1::y1 as bessel_y1;
pub use jn::jn as bessel_jn;
use std::sync::LazyLock;

const _HUGE: f64 = 1e300;
const TWO_M13: f64 = 1.0 / (1 << 13) as f64; // 2**-13 0x3f20000000000000
const TWO_M27: f64 = 1.0 / (1 << 27) as f64; // 2**-27 0x3e40000000000000
const TWO_M29: f64 = 1.0 / (1 << 29) as f64; // 2**-29
const TWO_M54: f64 = 1.0 / (1_u64 << 54) as f64; // 2**-54 0x3c90000000000000
static TWO_129: LazyLock<f64> = LazyLock::new(|| 2.0_f64.powi(129)); // 2**129 0x4800000000000000
static TWO_302: LazyLock<f64> = LazyLock::new(|| 2.0_f64.powi(302)); // 2**302
