

program test

   implicit none
   double precision :: zr, zi, nu, cyr, cyi, cwrkr, cwrki
   integer :: kode, n, nz, ierr, kind, id
   ! double precision, DIMENSION (1:9) :: cyr, cyi, cwrkr, cwrki
   ! x = 10.0
   ! print*, x

   nu = 0e0
   zr = 2.220446049250313e-16
   zi = 1e0
   kind = 2
   n = 1

   kode = 2
   id = 0
   call zbiry(zr, zi, id, kode, cyr, cyi, ierr)
   print*, cyr
   print*, cyi
   print*, nz
   print*, ierr

end program test
