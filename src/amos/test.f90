

program test

   implicit none
   double precision :: zr, zi, nu, cyr, cyi
   integer :: kode, n, nz, ierr, kind, id
   ! double precision, DIMENSION (1:9) :: cyr, cyi
   ! x = 10.0
   ! print*, x

   nu = 172302836.50840142
   zr = -4816.864663442315
   zi =  9.992997770079455
   kind = 2
   n = 9

   kode = 1
   id = 0
   call zairy(zr, zi, id, kode, cyr, cyi, nz, ierr)
   print*, cyr
   print*, cyi
   print*, nz
   print*, ierr

end program test
