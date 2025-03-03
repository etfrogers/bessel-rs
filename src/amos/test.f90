

program test

   implicit none
   double precision :: zr, zi, nu, cyr, cyi
   integer :: kode, n, nz, ierr

   ! x = 10.0
   ! print*, x

   zr = 2.1
   zi = 0.0
   nu = 4.0
   n = 1
   kode = 1
   call zbesj(zr, zi, nu, kode, n, cyr, cyi, nz, ierr)
   print*, cyr
   print*, cyi
   print*, nz
   print*, ierr

end program test
