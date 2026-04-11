

program test

   implicit none
   double precision :: zr, zi, nu!, cyr, cyi, cwrkr, cwrki
   integer :: kode, n, nz, ierr, kind, id
   double precision, DIMENSION (1:1) :: cyr, cyi, cwrkr, cwrki
   ! x = 10.0
   ! print*, x

   nu =  2.74d-288
   zr = 6.33d-166
   zi = 7.53d-275
   ! nu = 7107.636998006379
   ! zr = -1867.258055869096
   ! zi = 4865.284129480511
   ! kind = 2
   n = 1

   kode = 1
   id = 0
   call zbesy(zr, zi, nu, kode, n, cyr, cyi, nz, cwrkr, cwrki, ierr)
   print*, cyr
   print*, cyi
   print*, nz
   print*, ierr

end program test
