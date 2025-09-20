

program test

   implicit none
   double precision :: zr, zi, nu!, cyr, cyi, cwrkr, cwrki
   integer :: kode, n, nz, ierr, kind, id
   double precision, DIMENSION (1:9) :: cyr, cyi, cwrkr, cwrki
   ! x = 10.0
   ! print*, x

   nu = 1.0
   zr = 1.0
   zi = 1.0
   kind = 2
   n = 9

   kode = 2
   id = 0
   call zbesy(zr, zi, nu, kode, n, cyr, cyi, nz, cwrkr, cwrki, ierr)
   print*, cyr
   print*, cyi
   print*, nz
   print*, ierr

end program test
