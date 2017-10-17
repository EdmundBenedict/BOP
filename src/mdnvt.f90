 
      subroutine mdnvt(f,dr,vel,nd,mass,dt,temp)
          use mod_precision


!
!    This is a routine to perform a molecular dynamics movement of
!     the ions under the conditions of constant NVT. A Gaussian
!     thermostat and the Verlet algorithm are used here.
!

      implicit none

!
!    Declare the simple variables.
!

      real(dp) :: temp
      real(dp) :: dt
      real(dp) :: xfac,fac
      real(dp) :: amu,ang,ev,kb
      real(dp) :: tscale
      real(dp) :: beta,sum

      integer nd
      integer i,j

!
!    Define the parameters.
!

      parameter (amu = 1.660566e-27_dp)
      parameter (ang = 1.0e-10_dp)
      parameter (ev = 1.602189e-19_dp)
      parameter (kb = 1.3807e-23_dp)
      parameter (tscale = 1.0e-15_dp)
      parameter (xfac = (tscale/ang)*(tscale/ang)*(ev/amu))

!
!    Declare the arrays.
!

      real(dp) :: f(3,nd)
      real(dp) :: dr(3,nd)
      real(dp) :: vel(3,nd)
      real(dp) :: mass(nd)

!
!    Calculate the thermostat factor.
!

      sum = 0.0_dp
      do 3 i = 1,nd,1
         fac = 0.5_dp*dt*xfac/mass(i)
         sum = sum + mass(i)*((vel(1,i)+fac*f(1,i))**2 & 
     &                      + (vel(2,i)+fac*f(2,i))**2 & 
     &                      + (vel(3,i)+fac*f(3,i))**2)
 3    continue
      sum = sum*amu*ang*ang/(tscale*tscale)
      beta = sqrt(3.0_dp*real(nd-1, dp)*kb*temp/sum)

!
!    Update the positions and the velocities.
!

      do 1 i = 1,nd,1
         fac = dt*xfac/mass(i)
         do 2 j = 1,3,1
            vel(j,i) = vel(j,i)*(2.0_dp*beta-1.0_dp)+beta*f(j,i)*fac
            dr(j,i) = vel(j,i)*dt
 2       continue
 1    continue

      end

