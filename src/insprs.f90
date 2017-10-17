 
       function insprs(nd,t,vol,sum)
          use mod_precision


!
!    This is a function to evaluate the instantaneous pressure.
!

      implicit none
      real(dp) :: insprs

!
!    Declare the simple variables.
!

      real(dp) :: t,vol
      real(dp) :: sum
      real(dp) :: fac,ev,ang,kb

      integer nd

!
!    Define the parameters.
!

      parameter(kb = 1.38066e-23_dp)
      parameter(ev = 1.602189e-19_dp)
      parameter(ang = 1.0e-10_dp)
      parameter(fac = 1.0e-3_dp/(ang*ang*ang))

!
!    Evaluate the pressure.
!

      insprs = fac*(3.0_dp*real(nd, dp)*kb*t-ev*sum)/(3.0_dp*vol)

      end

