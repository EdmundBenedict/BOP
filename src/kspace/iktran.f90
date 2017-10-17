 
      subroutine iktran(kk,nk,a)
          use mod_precision
         use mod_const

!
!    This is a routine to transform the k points from
!     the Cartesian basis to the reciprocal lattice basis.
!

      implicit none

!
!    Declare the simple variables.
!

      integer nk,ik

!
!    Define the parameters.
!


!
!    Declare the arrays.
!

      real(dp) :: kk(3,nk)
      real(dp) :: tmp(3)
      real(dp) :: a(3,3)

!
!    Transform the K vectors.
!
!       print '(/)'
      do ik = 1,nk,1
!          print '("o ",3(x,f8.5))', kk(:,ik)
         tmp(1) = kk(1,ik)/(2.0_dp*pi)
         tmp(2) = kk(2,ik)/(2.0_dp*pi)
         tmp(3) = kk(3,ik)/(2.0_dp*pi)
         kk(1,ik) = tmp(1)*a(1,1)+tmp(2)*a(1,2)+tmp(3)*a(1,3)
         kk(2,ik) = tmp(1)*a(2,1)+tmp(2)*a(2,2)+tmp(3)*a(2,3)
         kk(3,ik) = tmp(1)*a(3,1)+tmp(2)*a(3,2)+tmp(3)*a(3,3)
!          print '("n ",3(x,f8.5))', kk(:,ik)
      enddo
      
!       print *,a

      end

