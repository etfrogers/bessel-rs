

program test

   implicit none
   double precision :: zr, zi, nu!, cyr, cyi, cwrkr, cwrki
   integer :: kode, n, nz, ierr, kind, id
   double precision, DIMENSION (1:15) :: cyr, cyi, cwrkr, cwrki
   ! x = 10.0
   ! print*, x

   nu = 7107.636998006379
   zr = -1867.258055869096
   zi = 4865.284129480511
   kind = 2
   n = 15

   kode = 1
   id = 0
   call zbesh(zr, zi, nu, kode, kind, n, cyr, cyi, nz, ierr)
   print*, cyr
   print*, cyi
   print*, nz
   print*, ierr

end program test
