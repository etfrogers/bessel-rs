// Translation of Go's math.Jn, Yn function
// from
// https://cs.opensource.google/go/go/+/master:src/math/j1.go
// That code has this notice:
// Copyright 2010 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

/*
    Bessel function of the first and second kinds of order n.
*/

// The original C code and the long comment below are
// from FreeBSD's /usr/src/lib/msun/src/e_jn.c and
// came with this notice.
// The go code is a simplified
// version of the original C.
//
// ====================================================
// Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
//
// Developed at SunPro, a Sun Microsystems, Inc. business.
// Permission to use, copy, modify, and distribute this
// software is freely granted, provided that this notice
// is preserved.
// ====================================================
//
// __ieee754_jn(n, x), __ieee754_yn(n, x)
// floating point Bessel's function of the 1st and 2nd kind
// of order n
//
// Special cases:
//      y0(0)=y1(0)=yn(n,0) = -inf with division by zero signal;
//      y0(-ve)=y1(-ve)=yn(n,-ve) are NaN with invalid signal.
// Note 2. About jn(n,x), yn(n,x)
//      For n=0, j0(x) is called,
//      for n=1, j1(x) is called,
//      for n<x, forward recursion is used starting
//      from values of j0(x) and j1(x).
//      for n>x, a continued fraction approximation to
//      J(n,x)/j(n-1,x) is evaluated and then backward
//      recursion is used starting from a supposed value
//      for J(n,x). The resulting value of J(0,x) is
//      compared with the actual value to correct the
//      supposed value of J(n,x).
//
//      yn(n,x) is similar in all respects, except
//      that forward recursion is used for all
//      values of n>1.
#![allow(clippy::excessive_precision)]
use num::Integer;

use super::{
    TWO_302, TWO_M29,
    j0::{j0, y0},
    j1::{j1, y1},
};
use std::f64::{self, consts::PI};

// Jn returns the order-n Bessel function of the first kind.
//
// Special cases are:
//
//	Jn(n, ±Inf) = 0
//	Jn(n, NaN) = NaN
pub fn jn(n: i32, x: f64) -> f64 {
    // special cases
    if x.is_nan() {
        return f64::NAN;
    }
    if x.is_infinite() {
        return 0.0;
    }
    // J(-n, x) = (-1)**n * J(n, x), J(n, -x) = (-1)**n * J(n, x)
    // Thus, J(-n, x) = J(n, -x)

    if n == 0 {
        return j0(x);
    }
    if x == 0.0 {
        return 0.0;
    }
    let (n, x) = if n < 0 { (-n, -x) } else { (n, x) };
    if n == 1 {
        return j1(x);
    }
    let (x, sign) = if x < 0.0 {
        (-x, n.is_odd()) // odd n and negative x
    } else {
        (x, false)
    };

    let b = if n as f64 <= x {
        // Safe to use J(n+1,x)=2n/x *J(n,x)-J(n-1,x)
        if x >= *TWO_302 {
            // x > 2**302

            // (x >> n**2)
            //          Jn(x) = cos(x-(2n+1)*pi/4)*sqrt(2/x*pi)
            //          Yn(x) = sin(x-(2n+1)*pi/4)*sqrt(2/x*pi)
            //          Let s=sin(x), c=cos(x),
            //              xn=x-(2n+1)*pi/4, sqt2 = sqrt(2),then
            //
            //                 n    sin(xn)*sqt2    cos(xn)*sqt2
            //              ----------------------------------
            //                 0     s-c             c+s
            //                 1    -s-c            -c+s
            //                 2    -s+c            -c-s
            //                 3     s+c             c-s

            let (s, c) = x.sin_cos();
            let temp = match n % 3 {
                0 => c + s,
                1 => -c + s,
                2 => -c - s,
                3 => c - s,
                _ => panic!("Not reachable"),
            };
            (1.0 / PI.sqrt()) * temp / x.sqrt()
        } else {
            let mut b = j1(x);
            let mut a = j0(x);
            for i in 1..n {
                let two_i: f64 = (i + i).into();
                (a, b) = (b, b * (two_i / x) - a) // avoid underflow
            }
            b
        }
    } else {
        if x < TWO_M29 {
            // x < 2**-29
            // x is tiny, return the first Taylor expansion of J(n,x)
            // J(n,x) = 1/n!*(x/2)**n  - ...
            if n > 33 {
                // underflow
                0.0
            } else {
                let n_factorial = (1..=n).fold(1.0, |acc, i| acc * i as f64);
                (x / 2.0).powi(n) / n_factorial
            }
            // b
        } else {
            // use backward recurrence
            //                      x      x**2      x**2
            //  J(n,x)/J(n-1,x) =  ----   ------   ------   .....
            //                      2n  - 2(n+1) - 2(n+2)
            //
            //                      1      1        1
            //  (for large x)   =  ----  ------   ------   .....
            //                      2n   2(n+1)   2(n+2)
            //                      -- - ------ - ------ -
            //                       x     x         x
            //
            // Let w = 2n/x and h=2/x, then the above quotient
            // is equal to the continued fraction:
            //                  1
            //      = -----------------------
            //                     1
            //         w - -----------------
            //                        1
            //              w+h - ---------
            //                     w+2h - ...
            //
            // To determine how many terms needed, let
            // Q(0) = w, Q(1) = w(w+h) - 1,
            // Q(k) = (w+k*h)*Q(k-1) - Q(k-2),
            // When Q(k) > 1e4	good for single
            // When Q(k) > 1e9	good for double
            // When Q(k) > 1e17	good for quadruple

            // determine k
            let two_n: f64 = (n + n).into();
            let w = two_n / x;
            let h = 2.0 / x;
            let mut q0 = w;
            let mut z = w + h;
            let mut q1 = w * z - 1.0;
            let mut k = 1;
            while q1 < 1e9 {
                k += 1;
                z += h;
                (q0, q1) = (q1, z * q1 - q0);
            }
            let m = n + n;
            let mut t = 0.0;
            for i in (m..(2 * (n + k))).step_by(2).rev() {
                t = 1.0 / (((i as f64) / x) - t);
            }
            let mut a = t;
            let mut b = 1.0;
            //  estimate log((2/x)**n*n!) = n*log(2/x)+n*ln(n)
            //  Hence, if n*(log(2n/x)) > ...
            //  single 8.8722839355e+01
            //  double 7.09782712893383973096e+02
            //  long double 1.1356523406294143949491931077970765006170e+04
            //  then recurrent value may overflow and the result is
            //  likely underflow to zero

            let mut tmp: f64 = n.into();
            let v = 2.0 / x;
            tmp = tmp * (v * tmp).abs().ln();
            if tmp < 7.09782712893383973096e+02 {
                // for i := n - 1; i > 0; i-- {
                for i in (1..n).rev() {
                    let di: f64 = (i + i).into();
                    (a, b) = (b, b * di / x - a);
                }
            } else {
                // for i := n - 1; i > 0; i-- {
                for i in (1..n).rev() {
                    let di: f64 = (i + i).into();
                    (a, b) = (b, b * di / x - a);
                    // scale b to avoid spurious overflow
                    if b > 1e100 {
                        a /= b;
                        t /= b;
                        b = 1.0;
                    }
                }
            }
            b = t * j0(x) / b;
            b
        }
    };
    if sign { -b } else { b }
}

// Yn returns the order-n Bessel function of the second kind.
//
// Special cases are:
//
//	Yn(n, +Inf) = 0
//	Yn(n ≥ 0, 0) = -Inf
//	Yn(n < 0, 0) = +Inf if n is odd, -Inf if n is even
//	Yn(n, x < 0) = NaN
//	Yn(n, NaN) = NaN
pub fn yn(n: i32, x: f64) -> Result<f64, String> {
    // special cases
    if x < 0.0 {
        return Err("yn is complex for z < 0, and this function is soley real".to_string());
    }
    if x.is_nan() {
        return Ok(f64::NAN);
    }
    if x.is_infinite() && x.is_sign_positive() {
        return Ok(0.0);
    }

    if n == 0 {
        return y0(x);
    }
    if x == 0.0 {
        return Ok(if n < 0 && n & 1 == 1 {
            f64::INFINITY
        } else {
            f64::NEG_INFINITY
        });
    }
    let (n, sign) = if n < 0 {
        (-n, n.is_odd()) // sign true if n < 0 && |n| odd
    } else {
        (n, false)
    };
    if n == 1 {
        return Ok(if sign { -y1(x)? } else { y1(x)? });
    }
    let mut b = if x >= *TWO_302 {
        // x > 2**302
        // (x >> n**2)
        //	    Jn(x) = cos(x-(2n+1)*pi/4)*sqrt(2/x*pi)
        //	    Yn(x) = sin(x-(2n+1)*pi/4)*sqrt(2/x*pi)
        //	    Let s=sin(x), c=cos(x),
        //		xn=x-(2n+1)*pi/4, sqt2 = sqrt(2),then
        //
        //		   n	sin(xn)*sqt2	cos(xn)*sqt2
        //		----------------------------------
        //		   0	 s-c		 c+s
        //		   1	-s-c 		-c+s
        //		   2	-s+c		-c-s
        //		   3	 s+c		 c-s
        let (s, c) = x.sin_cos();
        let temp = match n % 3 {
            0 => s - c,
            1 => -s - c,
            2 => -s + c,
            3 => s + c,
            _ => panic!("Not reachable"),
        };
        (1.0 / PI.sqrt()) * temp / x.sqrt()
    } else {
        let mut a = y0(x)?;
        let mut b = y1(x)?;
        for i in 1..n {
            let two_i: f64 = (i + i).into();
            (a, b) = (b, b * (two_i / x) - a) // avoid underflow
        }
        b
    };
    if sign {
        b = -b;
    }
    Ok(b)
}
