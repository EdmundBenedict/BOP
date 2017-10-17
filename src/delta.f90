 
       function delta(x,kt)
          use mod_precision


!
!    This is a function to evaluate the derivative of the fermi function.
!     NOTE: X = E-Ef
!

      implicit none
      real(dp) :: delta

!
!    Declare the simple variables.
!

      real(dp) :: x,y,z,kt,z1
      real(8), save :: old_kt = 0.3_dp, inv_kt = 1.0_dp/0.3_dp
      
      if (kt /= old_kt) inv_kt = 1.0_dp/kt
!
!    Evaluate the function.
!

      y = x*inv_kt
      if (y < -32.0_dp .or. y > 32.0_dp ) then
         delta = 0.0_dp
      else
         z = exp(y)
         z1 = z+1.0_dp
         delta = -inv_kt*z/(z1*z1)
      endif

      end

