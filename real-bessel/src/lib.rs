mod types;
pub use types::BesselError;

mod j0;
mod j1;
mod jn;
pub use j0::j0 as bessel_j0;
pub use j0::y0 as bessel_y0;
pub use j1::j1 as bessel_j1;
pub use j1::y1 as bessel_y1;
pub use jn::jn as bessel_jn;
pub use jn::yn as bessel_yn;

pub(crate) const TWO_M13: f64 = 1.0 / (1 << 13) as f64; // 2**-13 0x3f20000000000000
pub(crate) const TWO_M27: f64 = 1.0 / (1 << 27) as f64; // 2**-27 0x3e40000000000000
pub(crate) const TWO_M29: f64 = 1.0 / (1 << 29) as f64; // 2**-29
pub(crate) const TWO_M54: f64 = 1.0 / (1_u64 << 54) as f64; // 2**-54 0x3c90000000000000
pub(crate) const TWO_129: f64 = 6.805_647_338_418_77e38; // 2**129 0x4800000000000000
pub(crate) const TWO_302: f64 = 8.148_143_905_337_944e90; // 2**302 0x52D0000000000000
