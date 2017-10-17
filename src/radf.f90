 
      subroutine radf(r,vsss,vsps,vpss,vpps,vppp,vsds,vdss,vpds,vdps,vpdp,vdpp,vdds,vddp,vddd,bnddatin,bndsclin)

!
!    This is a subroutine to evaluate the radial functions. The energies are
!     in eV, and the distances are in Angstroms.
!
          use mod_precision


      implicit none

!
!    Declare the simple variables.
!

      real(dp) :: r,vscale
      real(dp) :: vsss,vsps,vpss,vpps,vppp
      real(dp) :: vsds,vdss,vpds,vdps,vpdp,vdpp
      real(dp) :: vdds,vddp,vddd

      integer i, j

!
!    Declare the arrays.
!

      real(dp) :: bnddatin(15)
      real(dp) :: bndsclin(14,15)
      real(dp) :: f(14)
!
!      EXTERNAL SCALE
!
!    Evaluate the radial function.
!

      do i=1,14,1
!         write(*,*) bndsclin(1,I)
         if (bndsclin(1,i) /= 0.0_dp) then 
!
         f(i) = vscale(r,bndsclin(1,i),bndsclin(2,i),bndsclin(5,i), & 
     &                bndsclin(6,i),bndsclin(3,i),bndsclin(4,i), & 
     &                bndsclin(7,i),bndsclin(8,i),bndsclin(9,i), & 
     &                bndsclin(10,i),bndsclin(11,i),bndsclin(12,i))
!
         else
            f(i) = 0.0_dp
         endif
      enddo

      vsss = bnddatin(1)*f(1)
      vsps = bnddatin(2)*f(2)
      vpss = bnddatin(3)*f(3)
      vpps = bnddatin(4)*f(4)
      vppp = bnddatin(5)*f(5)
      vsds = bnddatin(6)*f(6)
      vdss = bnddatin(7)*f(7)
      vpds = bnddatin(8)*f(8)
      vdps = bnddatin(9)*f(9)
      vpdp = bnddatin(10)*f(10)
      vdpp = bnddatin(11)*f(11)
      vdds = bnddatin(12)*f(12)
      vddp = bnddatin(13)*f(13)
      vddd = bnddatin(14)*f(14)

      end
