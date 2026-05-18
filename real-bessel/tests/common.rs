#[macro_export]
macro_rules! unwrap_real_bessel {
    ($fun: ident, $order:expr, $zr: expr) => {
        match $fun::<f64, f64>($order, $zr) {
            Ok(v) => v,
            Err(::amos_bessel_rs::BesselError::Overflow) => f64::NAN,
            Err(e) => panic!("Unexpected error from {}: {:?}", stringify!($fun), e),
        }
    };
}

#[macro_export]
macro_rules! get_real_y_bessel {
    ($fun: ident, $zr: expr, $action: expr) => {
        if $zr == 0.0 {
            // bessel_y0 returns -inf for 0, while bessel_y returns an error.
            $action;
        } else {
            // bessel_yX does not support negative inputs, but bessel_y does.
            // therefore if we get an invalid input error, we retry with the absolute value.
            // we could ignore this case, but it's better to test something.
            match $fun($zr) {
                Ok(v) => v,
                Err(details) => {
                    if details.contains("z < 0") {
                        $zr = $zr.abs();
                        $fun($zr).unwrap()
                    } else {
                        panic!("Unexpected error from {}: {:?}", stringify!($fun), details);
                    }
                }
            }
        }
    };
}
