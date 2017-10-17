 
       function parfrac(c,z0,n,z)
          use mod_precision


!
!    This is a function to evaluate partial fractions.
!

      implicit none
      complex(dp) :: parfrac

!
!    Declare the simple variables.
!

      complex(dp) :: z

      integer n,i

!
!    Declare the arrays.
!

      complex(dp) :: c(n)
      complex(dp) :: z0(n)

!
!    Evaluate the partial fraction
!

      parfrac = cmplx(0.0_dp, kind=dp)
      do i = 1,n,1
         parfrac = parfrac + cmplx(1.0_dp, kind=dp)/(c(i)*(z-z0(i)))
      enddo

      end

