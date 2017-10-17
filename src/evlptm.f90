 
      subroutine evlptm(ai,bi,i,a,b,g,phi,w)
          use mod_precision


!
!    This is a routine to evaluate the periodic terminator.
!

      implicit none

!
!    Declare the simple variables.
!

      real(dp) :: ai,bi
      real(dp) :: a,b,g,phi,w
      real(dp) :: x

      integer i

!
!    Evaluate the function.
!

      x = 2.0_dp*phi*real(i, dp)-w
      ai = a + g*cos(x)
      bi = b + 0.5_dp*g*cos(x-phi)

      end

