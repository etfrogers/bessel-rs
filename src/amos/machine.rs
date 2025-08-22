use std::sync::LazyLock;

/// `elim` is a nunber such that if you take `elim.exp()` or `(-elim).exp()`, you will create a
/// number that is at risk of overflowing `f64`. In this case "at risk of overflow" means
/// within a factor of 1000.0 of `f64::MAX_POSITIVE` for `elim.exp()` or
/// `f64::MIN_POSITIVE` for `(-elim).exp()`.
///
/// As the two conditions (positive and negative) above are not exactly the same, the code chooses
/// the most conservative (which as standard is the MIN_POSITIVE version)
///
/// `approximation_limit` is a number a bit smaller than `elim` above which scaled calculations are used to
/// ensure that over/underflow is mitigated. Here "a bit smaller" means "reduced by a factor of
/// the smaller of 18 digits or the number of (decimal) digits stored in the mantissa of `f64`
/// (which as standard is ~15.6). For `x > approximation_limit`, e^x will start to lose precision.
///
/// # See also
/// tests in `test_machine_consts.rs`
#[derive(Debug, Clone)]
pub(crate) struct MachineConsts {
    /// Originally ARM
    pub underflow_limit: f64,
    /// Originally ASCLE
    pub absolute_approximation_limit: f64,
    /// TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
    pub abs_error_tolerance: f64, // TOL
    /// ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT
    pub exponent_limit: f64,
    /// Originally ALIM
    pub approximation_limit: f64,
    /// DIG NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
    pub _significant_digits: f64, // DIG
    /// RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
    pub asymptotic_z_limit: f64, // RL
    /// FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU.
    pub asymptotic_order_limit: f64, // FNUL
    pub rtol: f64,
    /// CSSR
    pub scaling_factors: [f64; 3],
    /// CSRR
    pub reciprocal_scaling_factors: [f64; 3],
    pub bry: [f64; 3],
}

impl MachineConsts {
    fn new() -> Self {
        // Here we use approximate value, rather than calculating `10.0_f64.ln()`, as
        // this matches the the Fotran code, and the exact value causes subtle differences
        // in output (should just be what values are accpeted, but cause tetst to fail)
        let ln_10: f64 = 2.303;

        let underflow_limit = 1.0e+3 * 2.0 * f64::MIN_POSITIVE;
        let abs_error_tolerance = f64::EPSILON.max(1.0e-18);
        // absolute_approximation_limit == (-approximation_limit).exp() -- see test
        let absolute_approximation_limit = underflow_limit / abs_error_tolerance;

        let digits_per_bit = (f64::RADIX as f64).log10();
        let exponent_bit_limit = f64::MIN_EXP.abs().min(f64::MAX_EXP.abs());

        // Subtract 3.0 (digits) to give a number above which 10^decimal_exponent_limit would
        // be close to overflowing (i.e. within 1000 == 10^3.0 of the actual limit)
        let decimal_exponent_limit = (exponent_bit_limit as f64) * digits_per_bit - 3.0;
        // Multiplying by ln_10 converts from 10^x overflowing to e^x overflowing
        let exponent_limit = ln_10 * decimal_exponent_limit;

        let f64_siginficant_digits = digits_per_bit * ((f64::MANTISSA_DIGITS - 1) as f64);
        // siginficant_digits == abs_error_tolerance.log10() -- see test
        let significant_digits = f64_siginficant_digits.min(18.0);
        // Again, multiply number of base 10 digits by ln_10 to convert to e^x
        let approximation_limit = exponent_limit - (significant_digits * ln_10);

        let asymptotic_z_limit = 1.2 * significant_digits + 3.0;
        let asymptotic_order_limit = 10.0 + 6.0 * (significant_digits - 3.0);
        let rtol = 1.0 / abs_error_tolerance;
        //-----------------------------------------------------------------------
        //     COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG-
        //     NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE,
        //     EXP(ALIM)=EXP(ELIM)*TOL
        //-----------------------------------------------------------------------
        // let CSCL = machine_consts.rtol;
        // let CRSC = machine_consts.abs_error_tolerance;
        let scaling_factors = [rtol, 1.0, abs_error_tolerance];
        let reciprocal_scaling_factors = [abs_error_tolerance, 1.0, rtol];
        let bry = [
            absolute_approximation_limit,
            1.0 / absolute_approximation_limit,
            f64::MAX / 2.0,
        ];
        Self {
            underflow_limit,
            absolute_approximation_limit,
            abs_error_tolerance,
            exponent_limit,
            approximation_limit,
            _significant_digits: significant_digits,
            asymptotic_z_limit,
            asymptotic_order_limit,
            rtol,
            scaling_factors,
            reciprocal_scaling_factors,
            bry,
        }
    }
}

pub(crate) static MACHINE_CONSTANTS: LazyLock<MachineConsts> = LazyLock::new(MachineConsts::new);
