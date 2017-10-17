 
      subroutine rescale(v,n,x)
          use mod_precision


!
!    This is a routine to rescale a list of vectors by
!     a constant factor.
!

      implicit none

!
!    Declare the simple variables.
!

      real(dp) :: x

      integer n
      integer i

!
!    Declare the arrays.
!

      real(dp) :: v(3,n)

!
!    Rescale the vectors.
!

      do i = 1,n,1
         v(1,i) = v(1,i)*x
         v(2,i) = v(2,i)*x
         v(3,i) = v(3,i)*x
      end do

      end

