 
      subroutine rotbo(dr,bo,rbo)
          use mod_precision


!
!    This is a routine to evaluate the bond order matrix in
!     the rotated frame in which the bond lies along the Z axis.
!

      implicit none

!
!    Declare the arrays.
!

      real(dp) :: orm(9,9)
      real(dp) :: hcorm(9,9)
      real(dp) :: tmp(9,9)
      real(dp) :: dr(3)
      real(dp) :: bo(9,9)
      real(dp) :: rbo(9,9)

!
!    Evaluate the rotation matrix.
!

      call orbrot(dr,orm,hcorm)

!
!    Carry out the rotation.
!

      call dmatml(tmp,9,bo,9,9,9,orm,9,9)
      call dmatml(rbo,9,hcorm,9,9,9,tmp,9,9)

      end

