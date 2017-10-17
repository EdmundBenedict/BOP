
   function kentropy(ef)
      use mod_precision
      use mod_all_scalar, only : mag, nsp, kt
      use mod_const
      use mod_kspace
      use topologia, only : iproc, mpmap

!
!    This is a function to evaluate the entropic contribution to
!     the free energy of the electrons.
!

      implicit none
      
      real(dp), intent(in) :: ef
      
      real(dp) :: kentropy
      real(dp) :: entden
      real(dp) :: lkent

      integer ik,in,isp



      kentropy = 0.0_dp
      
      do isp = 1, nsp
         do ik = mpmap(iproc)+1, mpmap(iproc+1)
            lkent = 0.0_dp
            do in = 1, knh
               lkent = lkent + entden( enk(in, ik, isp) - ef, kt)
            end do
            kentropy = kentropy + wtk(ik) * lkent
         end do
      end do

      kentropy = -kt*kentropy
      
      if (.not. mag) kentropy = 2*kentropy

   end function kentropy

