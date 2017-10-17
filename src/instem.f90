 
       function instem(vel,nd,mass)
          use mod_precision


!
!    This is a routine to evaluate the instantaneous
!     temperature of the system.
!

      implicit none
      real(dp) :: instem

!
!    Declare the simple variables.
!

      real(dp) :: sum
      real(dp) :: fac,ang,tscale,amu
      real(dp) :: kb

      integer nd
      integer i,j

!
!    Define the parameters.
!

      parameter (kb = 1.38066e-23_dp)
      parameter (ang = 1.0e-10_dp)
      parameter (tscale = 1.0e-15_dp)
      parameter (amu = 1.660566e-27_dp)
      parameter (fac = amu*(ang/tscale)*(ang/tscale))

!
!    Declare the arrays.
!

      real(dp) :: vel(3,nd)
      real(dp) :: mass(nd)

!
!    Evaluate the temperature.
!

      sum = 0.0_dp

      do 1 i = 1,nd,1
         do 2 j = 1,3,1
            sum = sum + mass(i)*(vel(j,i)**2)
 2       continue
 1    continue

      instem = fac*sum/(3.0_dp*kb*real(nd-1, dp))

      end

