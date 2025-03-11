

program test

   implicit none
   double precision :: zr, zi, nu!, cyr, cyi
   integer :: kode, n, nz, ierr
   double precision, DIMENSION (1:9) :: cyr, cyi
   ! x = 10.0
   ! print*, x

   zr = 35.0001
   zi = 0.0
   nu = 340.0
   n = 9
   
   kode = 1
   call zbesj(zr, zi, nu, kode, n, cyr, cyi, nz, ierr)
   print*, cyr
   print*, cyi
   print*, nz
   print*, ierr

end program test
