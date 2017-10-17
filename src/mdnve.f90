 
      subroutine mdnve(oldf,f,dr,vel,nd,mass,dt,step)
          use mod_precision


!
!    This is a routine to perform a molecular dynamics movement of
!     the ions under the conditions of constant NVE.
!     The Beeman algorithm is used here.
!     (D.Beeman, J.Comput.Phys., vol 20, pp130-9 (1976))
!

      implicit none

!
!    Declare the simple variables.
!

      real(dp) :: dt
      real(dp) :: xfac,yfac
      real(dp) :: amu,ang,ev
      real(dp) :: tscale

      integer nd
      integer i,j
      integer step

!
!    Define the parameters.
!

      parameter (amu = 1.660566e-27_dp)
      parameter (ang = 1.0e-10_dp)
      parameter (ev = 1.602189e-19_dp)
      parameter (tscale = 1.0e-15_dp)
      parameter (xfac = (tscale/ang)*(tscale/ang)*(ev/amu))

!
!    Declare the arrays.
!

      real(dp) :: f(3,nd)
      real(dp) :: oldf(3,nd)
      real(dp) :: dr(3,nd)
      real(dp) :: vel(3,nd)
      real(dp) :: mass(nd)

!
!    Update the positions and the velocities.
!

      if (step == 1) then
         do i = 1,nd,1
            yfac = xfac/mass(i)
            do j = 1,3,1
               dr(j,i) = vel(j,i)*dt & 
     &                 + (2.0_dp/3.0_dp)*f(j,i)*yfac*dt*dt & 
     &                 - (1.0_dp/6.0_dp)*oldf(j,i)*yfac*dt*dt
               vel(j,i) = vel(j,i) & 
     &                  + (5.0_dp/6.0_dp)*f(j,i)*yfac*dt & 
     &                  - (1.0_dp/6.0_dp)*oldf(j,i)*yfac*dt
               oldf(j,i) = f(j,i)
            enddo
         enddo
      else
         do i = 1,nd,1
            yfac = xfac/mass(i)
            do j = 1,3,1
               vel(j,i) = vel(j,i) & 
     &                  + (1.0_dp/3.0_dp)*f(j,i)*yfac*dt
            enddo
         enddo
      endif

      end

