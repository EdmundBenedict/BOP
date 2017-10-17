 
       function ddelta(x,kt)
          use mod_precision


!
!    This is a function to evaluate the second
!     derivative of the fermi function.
!

      implicit none
      real(dp) :: ddelta

!
!    Declare the simple variables.
!

      real(dp) :: x,y,z,kt

!
!    Evaluate the function.
!

      y = x/kt
      if (y < -20.0_dp) then
         ddelta = 0.0_dp
      elseif (y > 20.0_dp) then
         ddelta = 0.0_dp
      else
         z = exp(y)
         ddelta = (z*z-z)/(((1.0_dp+z)**3)*kt*kt)
      endif

      end

