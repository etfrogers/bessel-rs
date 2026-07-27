#![allow(non_snake_case)]
use super::Scaling;
use crate::{
    BesselError, BesselFloat,
    amos::{recurrence::i_ratios, right_half_plane::k_right_half_plane},
};

use num::{Complex, complex::ComplexFloat};

// i_wronksian computes the i bessel function for re(z) >= 0.0 by
// normalizing the i function ratios from zrati by the Wronskian
// Originally ZWRSK
pub(crate) fn i_wronksian<T: BesselFloat>(
    zr: Complex<T>,
    order: T,
    scaling: Scaling,
    n: usize,
    y: &mut [Complex<T>],
) -> Result<usize, BesselError<T>> {
    //-----------------------------------------------------------------------
    //     I(FNU+I-1,Z) BY BACKWARD RECURRENCE FOR RATIOS
    //     Y(I)=I(FNU+I,Z)/I(FNU+I-1,Z) FROM CRATI NORMALIZED BY THE
    //     WRONSKIAN WITH K(FNU,Z) AND K(FNU+1,Z) FROM CBKNU.
    //-----------------------------------------------------------------------

    let (cw, _) = k_right_half_plane(zr, order, scaling, 2)?;
    let y_ratios = i_ratios(zr, order, n);
    //-----------------------------------------------------------------------
    //     RECUR FORWARD ON I(FNU+1,Z) = R(FNU,Z)*I(FNU,Z),
    //     R(FNU+J-1,Z)=Y(J),  J=1,...,N
    //-----------------------------------------------------------------------
    let mut cinu = T::C_ONE;
    if scaling == Scaling::Scaled {
        cinu = Complex::<T>::cis(zr.im);
    }
    //-----------------------------------------------------------------------
    //     ON LOW EXPONENT MACHINES THE K FUNCTIONS CAN BE CLOSE TO BOTH
    //     THE UNDER AND OVERFLOW LIMITS AND THE NORMALIZATION MUST BE
    //     SCALED TO PREVENT OVER OR UNDERFLOW. CUOIK HAS DETERMINED THAT
    //     THE RESULT IS ON SCALE.
    //-----------------------------------------------------------------------
    let acw = cw[1].abs();
    let CSCLR = if acw <= T::MACHINE_CONSTANTS.absolute_approximation_limit {
        T::one() / T::MACHINE_CONSTANTS.abs_error_tolerance
    } else if acw >= T::one() / T::MACHINE_CONSTANTS.absolute_approximation_limit {
        T::MACHINE_CONSTANTS.abs_error_tolerance
    } else {
        T::one()
    };
    let c1 = cw[0] * CSCLR;
    let c2 = cw[1] * CSCLR;
    //-----------------------------------------------------------------------
    //     CINU=CINU*(CONJG(CT)/CABS(CT))*(1.0/CABS(CT) PREVENTS
    //     UNDER- OR OVERFLOW PREMATURELY BY SQUARING CABS(CT)
    //-----------------------------------------------------------------------
    let mut ct = zr * (y_ratios[0] * c1 + c2);
    let abs_ct = ct.abs();
    ct = ct.conj() / abs_ct;
    cinu = (cinu / abs_ct) * ct;
    y[0] = cinu * CSCLR;
    for i in 1..n {
        cinu *= y_ratios[i - 1];
        y[i] = cinu * CSCLR;
    }
    Ok(0)
}
