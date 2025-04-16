use approx::assert_relative_eq;

use crate::amos::MachineConsts;

#[test]
fn test_machine_consts() {
    let mc = MachineConsts::new();
    dbg!((-mc.exponent_limit).exp(), 2.0 * f64::MIN_POSITIVE * 1000.0,);
    assert_relative_eq!((-mc.exponent_limit).exp(), 2.0 * f64::MIN_POSITIVE * 1000.0,);
    dbg!(mc._significant_digits, -mc.abs_error_tolerance.log10());
    assert_relative_eq!(mc._significant_digits, -mc.abs_error_tolerance.log10());
    dbg!(
        mc.absolute_approximation_limit,
        (-mc.approximation_limit).exp()
    );
    assert_relative_eq!(
        mc.absolute_approximation_limit,
        (-mc.approximation_limit).exp()
    );
}
