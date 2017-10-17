 
      subroutine getptm(a0,a1,b1,a,b,g,phi,w)
          use mod_precision


!
!    This is a routine to evaluate the coefficients for
!     the periodic terminator. The terminator has the
!     form:
!          a(i) = a +     G*cos(2*i*phi - w)
!          b(i) = b + 0.5*G*cos((2*i-1)*phi - w)
!    (See Turchi, Ducastelle, Treglia,
!     J.Phys.C:Sol.State Phys. vol 15, p2891 (1982))
!

      implicit none

!
!    Declare the simple variables.
!

      real(dp) :: a0,a1,b1
      real(dp) :: a,b,g,phi,w
      real(dp) :: x,y

!
!    Evaluate the coefficients.
!

      if (abs(b1-b) > 1.0e-6_dp) then
         x = (a0+a1-2.0_dp*a)/(4.0_dp*(b1-b))
         if (abs(x) <= 1.0_dp) then
            phi = acos(x)
            y = (a0-a)*sqrt(1.0_dp-x*x)
            if (abs(y) > 1.0e-6_dp) then
               w = atan((2.0_dp*(b1-b)-(a0-a)*x)/y)
               g = (a0-a)/cos(w)
            else
               g = 0.0_dp
               phi = 0.0_dp
               w = 0.0_dp
            endif
         else
            g = 0.0_dp
            phi = 0.0_dp
            w = 0.0_dp
         endif
      else
         g = 0.0_dp
         phi = 0.0_dp
         w = 0.0_dp
      endif

      end

