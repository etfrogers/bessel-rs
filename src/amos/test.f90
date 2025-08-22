

program test

   implicit none
   double precision :: zr, zi, nu!, cyr, cyi
   integer :: kode, n, nz, ierr, kind
   double precision, DIMENSION (1:9) :: cyr, cyi
   ! x = 10.0
   ! print*, x

   nu = 172302836.50840142
   zr = 1.2494954195932068e-254
   zi = -981457506.3179193
   kind = 2
   n = 9

   kode = 1
   call zbesi(zr, zi, nu, kode, n, cyr, cyi, nz, ierr)
   print*, cyr
   print*, cyi
   print*, nz
   print*, ierr

end program test
