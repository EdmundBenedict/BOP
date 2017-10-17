 
       function entden(x,kt)
          use mod_precision


!
!    This is a function to evluate the entropy density function.
!

      implicit none
      real(dp) :: entden

!
!    Declare the simple variables.
!

      real(dp) :: x,kt,theta,n,m

!
!    Evaluate the function.
!

      n = theta(x,kt)
      m = 1.0_dp - n

      if ((n < 1.0e-10_dp).or.(m < 1.0e-10_dp)) then
         entden = 0.0_dp
      else
         entden = n*log(n) + m*log(m)
      endif

      end

