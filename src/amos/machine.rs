pub fn d1mach(i: u8) -> f64 {
    //*********************************************************************72
    //
    //  D1MACH returns double precision real machine-dependent constants.
    //
    //  Discussion:
    //
    //    D1MACH can be used to obtain machine-dependent parameters
    //    for the local machine environment.  It is a function
    //    with one input argument, and can be called as follows:
    //
    //      D = D1MACH ( I )
    //
    //    where I=1,...,5.  The output value of D above is
    //    determined by the input value of I:.
    //
    //    D1MACH ( 1) = B**(EMIN-1), the smallest positive magnitude.
    //    D1MACH ( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
    //    D1MACH ( 3) = B**(-T), the smallest relative spacing.
    //    D1MACH ( 4) = B**(1-T), the largest relative spacing.
    //    D1MACH ( 5) = LOG10(B)
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license.
    //
    //  Modified:
    //
    //    25 April 2007
    //
    //  Author:
    //
    //    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
    //    This FORTRAN77 version by John Burkardt.
    //
    //  Reference:
    //
    //    Phyllis Fox, Andrew Hall, Norman Schryer,
    //    Algorithm 528:
    //    Framework for a Portable Library,
    //    ACM Transactions on Mathematical Software,
    //    Volume 4, Number 2, June 1978, page 176-188.
    //
    //  Parameters:
    //
    //    Input, integer I, the index of the desired constant.
    //
    //    Output, double precision D1MACH, the value of the constant.
    //
    // implicit none

    // double precision d1mach
    // integer i

    match i {
        1 => 4.450147717014403e-308,
        2 => 8.988465674311579e+307,
        3 => 1.110223024625157e-016,
        4 => 2.220446049250313e-016,
        5 => 0.301029995663981e+000,
        _ => panic!(
            "D1MACH - Fatal error!
  The input argument I is out of bounds.
  Legal values satisfy 1 <= I <= 5.
  I = {i}"
        ),
    }
}

pub fn i1mach(i: u8) -> i32 {
    // *********************************************************************72
    //
    //c I1MACH returns integer machine dependent constants.
    //
    //  Discussion:
    //
    //    Input/output unit numbers.
    //
    //      I1MACH(1) = the standard input unit.
    //      I1MACH(2) = the standard output unit.
    //      I1MACH(3) = the standard punch unit.
    //      I1MACH(4) = the standard error message unit.
    //
    //    Words.
    //
    //      I1MACH(5) = the number of bits per integer storage unit.
    //      I1MACH(6) = the number of characters per integer storage unit.
    //
    //    Integers.
    //
    //    Assume integers are represented in the S digit base A form:
    //
    //      Sign * (X(S-1)*A**(S-1) + ... + X(1)*A + X(0))
    //
    //    where 0 <= X(1:S-1) < A.
    //
    //      I1MACH(7) = A, the base.
    //      I1MACH(8) = S, the number of base A digits.
    //      I1MACH(9) = A**S-1, the largest integer.
    //
    //    Floating point numbers
    //
    //    Assume floating point numbers are represented in the T digit
    //    base B form:
    //
    //      Sign * (B**E) * ((X(1)/B) + ... + (X(T)/B**T) )
    //
    //    where 0 <= X(I) < B for I=1 to T, 0 < X(1) and EMIN <= E <= EMAX.
    //
    //      I1MACH(10) = B, the base.
    //
    //    Single precision
    //
    //      I1MACH(11) = T, the number of base B digits.
    //      I1MACH(12) = EMIN, the smallest exponent E.
    //      I1MACH(13) = EMAX, the largest exponent E.
    //
    //    Double precision
    //
    //      I1MACH(14) = T, the number of base B digits.
    //      I1MACH(15) = EMIN, the smallest exponent E.
    //      I1MACH(16) = EMAX, the largest exponent E.
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license.
    //
    //  Modified:
    //
    //    25 April 2007
    //
    //  Author:
    //
    //    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
    //    This FORTRAN77 version by John Burkardt.
    //
    //  Reference:
    //
    //    Phyllis Fox, Andrew Hall, Norman Schryer,
    //    Algorithm 528,
    //    Framework for a Portable Library,
    //    ACM Transactions on Mathematical Software,
    //    Volume 4, Number 2, June 1978, page 176-188.
    //
    //  Parameters:
    //
    //    Input, integer I, chooses the parameter to be returned.
    //    1 <= I <= 16.
    //
    //    Output, integer I1MACH, the value of the chosen parameter.
    //
    // implicit none

    // integer i
    // integer i1mach

    match i {
        1 => 5,
        2 => 6,
        3 => 7,
        4 => 6,
        5 => 32,
        6 => 4,
        7 => 2,
        8 => 31,
        9 => 2147483647,
        10 => 2,
        11 => 24,
        12 => -125,
        13 => 128,
        14 => 53,
        15 => -1021,
        16 => 1024,
        _ => panic!(
            "I1MACH - Fatal error!
The input argument I is out of bounds.
Legal values satisfy 1 <= I <= 16.
I = {i}"
        ),
    }
}
/*      function r1mach ( i )

// *********************************************************************72
//
//c R1MACH returns single precision real machine constants.
//
//  Discussion:
//
//    Assume that single precision real numbers are stored with a mantissa
//    of T digits in base B, with an exponent whose value must lie
//    between EMIN and EMAX.  Then for values of I between 1 and 5,
//    R1MACH will return the following values:
//
//      R1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
//      R1MACH(2) = B**EMAX*(1-B**(-T)), the largest magnitude.
//      R1MACH(3) = B**(-T), the smallest relative spacing.
//      R1MACH(4) = B**(1-T), the largest relative spacing.
//      R1MACH(5) = log10(B)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 April 2007
//
//  Author:
//
//    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
//    This FORTRAN77 version by John Burkardt.
//
//  Reference:
//
//    Phyllis Fox, Andrew Hall, Norman Schryer,
//    Algorithm 528,
//    Framework for a Portable Library,
//    ACM Transactions on Mathematical Software,
//    Volume 4, Number 2, June 1978, page 176-188.
//
//  Parameters:
//
//    Input, integer I, chooses the parameter to be returned.
//    1 <= I <= 5.
//
//    Output, real R1MACH, the value of the chosen parameter.
//
      implicit none

      integer i
      real r1mach

      if ( i < 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R1MACH - Fatal error!'
        write ( *, '(a)' ) '  The input argument I is out of bounds.'
        write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
        write ( *, '(a,i12)' ) '  I = ', i
        r1mach = 0.0E+00
        stop
      else if ( i == 1 ) then
        r1mach = 1.1754944E-38
      else if ( i == 2 ) then
        r1mach = 3.4028235E+38
      else if ( i == 3 ) then
        r1mach = 5.9604645E-08
      else if ( i == 4 ) then
        r1mach = 1.1920929E-07
      else if ( i == 5 ) then
        r1mach = 0.3010300E+00
      else if ( 5 < i ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R1MACH - Fatal error!'
        write ( *, '(a)' ) '  The input argument I is out of bounds.'
        write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
        write ( *, '(a,i12)' ) '  I = ', i
        r1mach = 0.0E+00
        stop
      end if

      return
      end
      subroutine timestamp ( )

// *********************************************************************72
//
//c TIMESTAMP prints out the current YMDHMS date as a timestamp.
//
//  Discussion:
//
//    This FORTRAN77 version is made available for cases where the
//    FORTRAN90 version cannot be used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
      implicit none

      character * ( 8 ) ampm
      integer d
      character * ( 8 ) date
      integer h
      integer m
      integer mm
      character * ( 9 ) month(12)
      integer n
      integer s
      character * ( 10 ) time
      integer y

      save month

      data month /
     &  'January  ', 'February ', 'March    ', 'April    ',
     &  'May      ', 'June     ', 'July     ', 'August   ',
     &  'September', 'October  ', 'November ', 'December ' /

      call date_and_time ( date, time )

      read ( date, '(i4,i2,i2)' ) y, m, d
      read ( time, '(i2,i2,i2,1x,i3)' ) h, n, s, mm

      if ( h .lt. 12 ) then
        ampm = 'AM'
      else if ( h .eq. 12 ) then
        if ( n .eq. 0 .and. s .eq. 0 ) then
          ampm = 'Noon'
        else
          ampm = 'PM'
        end if
      else
        h = h - 12
        if ( h .lt. 12 ) then
          ampm = 'PM'
        else if ( h .eq. 12 ) then
          if ( n .eq. 0 .and. s .eq. 0 ) then
            ampm = 'Midnight'
          else
            ampm = 'AM'
          end if
        end if
      end if

      write ( *,
     &  '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' )
     &  d, month(m), y, h, ':', n, ':', s, '.', mm, ampm

      return
      end
*/
