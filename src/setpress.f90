 
      subroutine setpress()
          use mod_precision


          use mod_all_scalar

          use mod_const

!
!    This is a routine to find the lattice constant that gives the correct
!     pressure.
!

      implicit none

!
!    Include constants.
!

!      include "Include/ALL.const"

!
!    Include scalars.
!

!      include "Include/ALL.scalar"

!
!    Include arrays.
!

      include "Include/PosVel.array"

!
!    Declare the simple variables.
!

      real(dp) :: insprs,p,perr
      real(dp) :: x,dx,dxmax,p1,delta,deltamx

      integer flag

!
!    Define the parameters.
!

      parameter (deltamx = 1.0e-3_dp)
      parameter (dxmax = 1.0e-2_dp)
      parameter (perr = 0.01_dp)

!
!    Calculate the present pressure.
!

      flag = 1
      call getetot(flag)
      p = insprs(nd,temp,vol,rdfbs+rdfpp)

!
!    Increase the volume a little, and calculate the new pressure.
!

      delta = deltamx

 98   if (p > press) then
         dx = delta
      else
         dx = -delta
      endif

      x = 1.0_dp + dx
      call rescale(ad,nd,x)
      call rescale(adinert,ninert,x)
      call rescale(a,3,x)
      vol = vol*(1.0_dp+dx*(3.0_dp+dx*(3.0_dp+dx)))

      flag = 1
      call getetot(flag)
      p1 = insprs(nd,temp,vol,rdfbs+rdfpp)

!
!    Calculate the change in volume needed to reach desired
!     pressure.
!

      delta = min(deltamx,abs(dx*(p-press)/(p-p1)))
      dx = (press-p1)*dx/(p1-p)
      if (abs(dx) > dxmax) then
         if (dx < 0.0_dp) then
            dx = -dxmax
         else
            dx = dxmax
         endif
      endif

      x = 1.0_dp + dx
      call rescale(ad,nd,x)
      call rescale(adinert,ninert,x)
      call rescale(a,3,x)
      vol = vol*(1.0_dp+dx*(3.0_dp+dx*(3.0_dp+dx)))

!
!    Calculate new pressure.
!

      flag = 1
      call getetot(flag)
      p = insprs(nd,temp,vol,rdfbs+rdfpp)

      write(6,'(''Pressure = '',G12.5)') p

!
!    If it is too far from the desired value, rescale again.
!

      if ((abs((press-p)/(press+p)) > perr).and. & 
     &    (abs(delta) > 1.0e-6_dp)) goto 98

      end

