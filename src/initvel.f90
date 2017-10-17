 
      subroutine initvel(vel,nd,temp,mass)
          use mod_precision


!
!    This is a routine to assign random velocities
!     the ions consistent with a given temperature.
!

      implicit none

!
!    Declare the simple variables.
!

      real(dp) :: temp
      real(dp) :: fac,kb,amu,ang
      real(dp) :: sd
      real(dp) :: gasdev
      real(dp) :: tscale
      real(dp) :: cmmass

      integer nd
      integer i,j

!
!    Define the parameters.
!

      parameter(kb = 1.38066e-23_dp)
      parameter(amu = 1.660566e-27_dp)
      parameter(ang = 1.0e-10_dp)
      parameter(tscale = 1.0e-15_dp)
      parameter(fac = (kb/amu)*(tscale/ang)*(tscale/ang))

!
!    Declare the arrays.
!

      real(dp) :: vel(3,nd)
      real(dp) :: pcm(3)
      real(dp) :: vcm(3)
      real(dp) :: mass(nd)

!
!    Initialise the velocities.
!

      do 1 i = 1,nd,1
         sd = sqrt(fac*temp/mass(i))
         do 2 j = 1,3,1
            vel(j,i) = gasdev(100)*sd
 2       continue
 1    continue

!
!    Remove the CM motion.
!

      pcm(1) = 0.0_dp
      pcm(2) = 0.0_dp
      pcm(3) = 0.0_dp
      cmmass = 0.0_dp

      do i = 1,nd,1
         pcm(1) = pcm(1) + mass(i)*vel(1,i)
         pcm(2) = pcm(2) + mass(i)*vel(2,i)
         pcm(3) = pcm(3) + mass(i)*vel(3,i)
         cmmass = cmmass + mass(i)
      enddo

      vcm(1) = pcm(1)/cmmass
      vcm(2) = pcm(2)/cmmass
      vcm(3) = pcm(3)/cmmass

      do i = 1,nd,1
         vel(1,i) = vel(1,i) - vcm(1)
         vel(2,i) = vel(2,i) - vcm(2)
         vel(3,i) = vel(3,i) - vcm(3)
      enddo

      end

!----------------------------------------------------------------------

       function gasdev(i)
          use mod_precision


!
!    This is a numerical recipes routine to produce random numbers
!     with a Gaussian distribution with mean 0 and standard deviation
!     of 1.
!

      implicit none
      real(dp) :: gasdev

!
!    Declare the simple variables.
!

      real(dp) :: v1,v2
      real(dp) :: r,fac
      real(dp) :: gset
      real(dp) :: ran1

      integer i
      integer iset

!
!    Assign a data value.
!

      data iset/0/

      save gset

!
!    Evaluate the random numbers.
!

      if (iset == 0) then
 1       v1 = 2.0_dp*ran1(i)-1.0_dp
         v2 = 2.0_dp*ran1(i)-1.0_dp
         r = v1**2 + v2**2
         if (r > 1.0_dp) goto 1
         fac = sqrt(-2.0_dp*log(r)/r)
         gset = v1*fac
         gasdev = v2*fac
         iset = 1
      else
         gasdev = gset
         iset = 0
      endif

      end

