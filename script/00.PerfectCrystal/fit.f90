program fit
   implicit none
   real*8 :: emin, v0, k0, xsi, birch, demin, dv0, dk0, dxsi, chi, chi_new, &
             v0_new, k0_new, emin_new
   integer :: nseek, i, j, n, icycle, ncycle
   real*8, allocatable :: x(:), e(:)

   print *, "Enter: emin, v0, k0, xsi, demin, dv0, dk0, dxsi, nseek, ncycle"
   read (*, *) emin, v0, k0, xsi, demin, dv0, dk0, dxsi, nseek, ncycle
   open (1, file="tmp", status="old")
   n = 0
   do
      read (1, *, end=10)
      n = n + 1
   end do
10 continue
   allocate (x(n), e(n))
   rewind (1)
   chi = 0.0
   print *, "chi x(i) e(i) i"
   do i = 1, n
      read (1, *) x(i), e(i)
      chi = chi + (e(i) - birch(emin, v0, k0, xsi, x(i)))**2
      print *, chi, x(i), e(i), i
   end do

   do icycle = 1, ncycle
      do i = 1, nseek
         emin_new = emin + demin*(rand() - 0.5)
         v0_new = v0 + dv0*(rand() - 0.5)
         k0_new = k0 + dk0*(rand() - 0.5)
         chi_new = 0.0
         do j = 1, n
            chi_new = chi_new + (e(j) - birch(emin_new, v0_new, k0_new, xsi, x(j)))**2
         end do
         if (chi_new < chi) then
            emin = emin_new
            v0 = v0_new
            k0 = k0_new
            chi = chi_new
         end if
      end do
      demin = demin/10
      dv0 = dv0/10
      dk0 = dk0/10
   end do
   print *, "Final set of parameters:"
   print *, "emin = ", emin
   print *, "v0 = ", v0
   print *, "k0 = ", k0
   print *, "chi = ", chi
end program

function birch(emin, v0, k0, xsi, x)
   implicit none
   real*8 :: emin, v0, k0, xsi, x, birch
   birch = emin + 1.5*v0*k0/1602*(0.75*(1 + 2*xsi)*(v0/x)**(4.0/3.0) &
                                  - xsi/2*(v0/x)**2 - 3.0/2.0*(1 + xsi)*(v0/x)**(2.0/3.0) + 0.5*(xsi + 1.5))
end
