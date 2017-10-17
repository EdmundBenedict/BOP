 
      subroutine dradf(r,dvsss,dvsps,dvpss,dvpps,dvppp, &
     &                   dvsds,dvdss,dvpds,dvdps,dvpdp,dvdpp, & 
     &                   dvdds,dvddp,dvddd,bnddat,bndscl)

!
!    This is a subroutine to evaluate the derivative of radial functions.
!    The energies are in eV, and the distances are in Angstroms.
!
          use mod_precision
      implicit none

!
!    Declare the simple variables.
!

      real(dp) :: r,dscale
      real(dp) :: dvsss,dvsps,dvpss,dvpps,dvppp
      real(dp) :: dvsds,dvdss,dvpds,dvdps,dvpdp,dvdpp
      real(dp) :: dvdds,dvddp,dvddd
      integer i

!
!    Declare the arrays.
!

      real(dp) :: bnddat(15)
      real(dp) :: bndscl(14,15)
      real(dp) :: f(14)

!
!    Evaluate the derivative of the radial function.
!

      do i=1,14,1

         if (bndscl(1,i) == 0.0) then
            f(i) = 0.0_dp
         else

         f(i) = dscale(r,bndscl(1,i),bndscl(2,i),bndscl(5,i), & 
     &                 bndscl(6,i),bndscl(3,i),bndscl(4,i), & 
     &                 bndscl(7,i),bndscl(8,i),bndscl(9,i), & 
     &                 bndscl(10,i),bndscl(11,i),bndscl(12,i))
         endif

      enddo

      dvsss = bnddat(1)*f(1)
      dvsps = bnddat(2)*f(2)
      dvpss = bnddat(3)*f(3)
      dvpps = bnddat(4)*f(4)
      dvppp = bnddat(5)*f(5)
      dvsds = bnddat(6)*f(6)
      dvdss = bnddat(7)*f(7)
      dvpds = bnddat(8)*f(8)
      dvdps = bnddat(9)*f(9)
      dvpdp = bnddat(10)*f(10)
      dvdpp = bnddat(11)*f(11)
      dvdds = bnddat(12)*f(12)
      dvddp = bnddat(13)*f(13)
      dvddd = bnddat(14)*f(14)

      end
