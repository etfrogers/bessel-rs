use std::sync::LazyLock;

use crate::types::BesselFloat;

/// `exponent_limit` is a number such that if you take `exponent_limit.exp()` or `(-exponent_limit).exp()`
/// , you will create a number that is at risk of overflowing `T`. In this case "at risk of overflow" means
/// within a factor of 1000.0 of `T::MAX_POSITIVE` for `exponent_limit.exp()` or
/// `T::MIN_POSITIVE` for `(-exponent_limit).exp()`.
///
/// As the two conditions (positive and negative) above are not exactly the same, the code chooses
/// the most conservative (which as standard is the MIN_POSITIVE version)
///
/// `approximation_limit` is a number a bit smaller than `exponent_limit` above which scaled calculations are used to
/// ensure that over/underflow is mitigated. Here "a bit smaller" means "reduced by a factor of
/// the smaller of 18 digits or the number of (decimal) digits stored in the mantissa of `T`
/// (which as standard is ~15.6). For `x > approximation_limit`, e^x will start to lose precision.
///
/// # See also
/// tests in `test_machine_consts.rs`
#[derive(Debug, Clone)]
pub(crate) struct MachineConsts<T: BesselFloat> {
    /// Originally ARM
    pub underflow_limit: T,
    /// Originally ASCLE
    pub absolute_approximation_limit: T,
    /// Originally TOL. The approximate unit roundoff limited to 1.0e-18.
    pub abs_error_tolerance: T,
    /// Originally ELIM. The approximate exponential over- and under-flow limit
    pub exponent_limit: T,
    /// Originally ALIM
    pub approximation_limit: T,
    /// Originally DIG. Number of base 10 digits in abs_error_tolerance = 10.pow(-dig).
    pub _significant_digits: T,
    /// Originally RL. The lower boundary of the asymptotic expansion for large z.
    pub asymptotic_z_limit: T,
    /// Originally FNUL.  The lower boundary of the asymptotic series for large order.
    pub asymptotic_order_limit: T,
    pub rtol: T,
    /// Originally CSSR
    pub scaling_factors: [T; 3],
    /// Originally CSRR
    pub reciprocal_scaling_factors: [T; 3],
    /// Originally BRY
    pub overflow_boundary: [T; 3],
}

impl<T: BesselFloat> MachineConsts<T> {
    fn new() -> Self {
        // Here we use approximate value, rather than calculating `10.0_f64.ln()`, as
        // this matches the the Fotran code, and the exact value causes subtle differences
        // in output (should just be what values are accepted, but cause tests to fail)
        let ln_10: T = T::from_f64(2.303);

        let underflow_limit = T::from_f64(1.0e+3 * 2.0) * T::min_positive_value();
        let abs_error_tolerance = T::epsilon().max(T::from_f64(1.0e-18));
        // absolute_approximation_limit == (-approximation_limit).exp() -- see test
        let absolute_approximation_limit = underflow_limit / abs_error_tolerance;

        let digits_per_bit: T = (T::from_f64(T::RADIX as f64)).log10();
        let exponent_bit_limit: T = T::from_f64((T::MIN_EXP.abs().min(T::MAX_EXP.abs())) as f64);

        // Subtract 3.0 (digits) to give a number above which 10^decimal_exponent_limit would
        // be close to overflowing (i.e. within 1000 == 10^3.0 of the actual limit)
        let decimal_exponent_limit = exponent_bit_limit * digits_per_bit - T::from_f64(3.0);
        // Multiplying by ln_10 converts from 10^x overflowing to e^x overflowing
        let exponent_limit = ln_10 * decimal_exponent_limit;

        let base_type_siginficant_digits: T =
            digits_per_bit * (T::from_f64(T::MANTISSA_DIGITS as f64) - T::one());
        // siginficant_digits == abs_error_tolerance.log10() -- see test
        let significant_digits = base_type_siginficant_digits.min(T::from_f64(18.0));
        // Again, multiply number of base 10 digits by ln_10 to convert to e^x
        let approximation_limit = exponent_limit - (significant_digits * ln_10);

        let asymptotic_z_limit = T::from_f64(1.2) * significant_digits + T::from_f64(3.0);
        let asymptotic_order_limit =
            T::from_f64(10.0) + T::from_f64(6.0) * (significant_digits - T::from_f64(3.0));
        let rtol = T::from_f64(1.0) / abs_error_tolerance;
        //-----------------------------------------------------------------------
        //     COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG-
        //     NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE,
        //     EXP(ALIM)=EXP(ELIM)*TOL
        //-----------------------------------------------------------------------
        let scaling_factors = [rtol, T::from_f64(1.0), abs_error_tolerance];
        let reciprocal_scaling_factors = [abs_error_tolerance, T::from_f64(1.0), rtol];
        let overflow_boundary = [
            absolute_approximation_limit,
            T::from_f64(1.0) / absolute_approximation_limit,
            T::max_value() / T::from_f64(2.0),
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
            overflow_boundary,
        }
    }
}

pub(crate) static MACHINE_CONSTANTS: LazyLock<MachineConsts<f64>> =
    LazyLock::new(MachineConsts::new);

pub(crate) static MACHINE_CONSTANTS_64: LazyLock<MachineConsts<f64>> =
    LazyLock::new(MachineConsts::new);
