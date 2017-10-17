 
       function integ(a,b,emax,f)
          use mod_precision


!
!    This is a function to use Simpson's rule to integrate a
!     given function over a given interval.
!

      implicit none
      real(dp) :: integ

!
!    Declare the simple variables.
!

      real(dp) :: h
      real(dp) :: fxm2h,fxmh,fx,fxph,fxp2h
      real(dp) :: f
      real(dp) :: x
      real(dp) :: emax
      real(dp) :: es
      real(dp) :: esmax
      real(dp) :: a,b
      real(dp) :: sum

      integer stkmax
      integer stkptr

!
!    Define the parameters.
!

      parameter(stkmax = 42)

!
!    Declare the arrays.
!

      real(dp) :: stka(stkmax)
      real(dp) :: stkfa(stkmax)
      real(dp) :: stkb(stkmax)
      real(dp) :: stkfb(stkmax)
      real(dp) :: stkm(stkmax)
      real(dp) :: stkfm(stkmax)

!
!    Declare external functions.
!

      external f

!
!    Initialise some variables.
!

      sum = 0.0_dp
      stkptr = 0
      esmax = 2.0_dp*emax/(b-a)

!
!    Push the limits of the integral onto the stack.
!

      stkptr = stkptr + 1

      stka(stkptr) = a
      stkfa(stkptr) = f(a)

      stkb(stkptr) = b
      stkfb(stkptr) = f(b)

      stkm(stkptr) = (a+b)/2.0_dp
      stkfm(stkptr) = f((a+b)/2.0_dp)

!
!    Evaluate the integral, bisecting the intervals until the
!     error from the interval is small enough.
!

!
!    WHILE the stack is not empty DO ....
!

 1    if (stkptr > 0) then

!
!       POP the limits from the stack.
!

         x = stkm(stkptr)
         h = (stkb(stkptr) - stka(stkptr))/4.0_dp

!
!       POP the function values needed for the error estimate.
!

         fxm2h = stkfa(stkptr)
         fxp2h = stkfb(stkptr)
         fx = stkfm(stkptr)

         stkptr = stkptr - 1

!
!       Evaluate the remaining function values needed for the
!        error estimate.
!

         fxmh = f(x-h)
         fxph = f(x+h)

!
!       Evaluate the error estimate.
!

         es = abs(fxm2h-4.0_dp*fxmh+6.0_dp*fx-4.0_dp*fxph+fxp2h)/90.0_dp

!
!       If the error is too large then bisect the interval, otherwise
!        add the contribution from the interval to the integral.
!

         if (es > esmax) then

!
!          PUSH the lower limit to the mid-point onto the stack.
!

            stkptr = stkptr + 1

            stka(stkptr) = x - 2.0_dp*h
            stkfa(stkptr) = fxm2h

            stkb(stkptr) = x
            stkfb(stkptr) = fx

            stkm(stkptr) = x - h
            stkfm(stkptr) = fxmh

!
!          PUSH the mid-point to the upper limit onto the stack.
!

            stkptr = stkptr + 1

            stka(stkptr) = x
            stkfa(stkptr) = fx

            stkb(stkptr) = x + 2.0_dp*h
            stkfb(stkptr) = fxp2h

            stkm(stkptr) = x + h
            stkfm(stkptr) = fxph

         else

!
!          Add the contribution to the total integral.
!

            sum = sum + 2.0_dp*h*(fxm2h + 4.0_dp*fx + fxp2h)/3.0_dp

         endif

         if (stkptr <= stkmax-2) then
            goto 1
         else
            print *,'Stack too full. Abandoning integral.'
            sum = -9.999999e30_dp
         endif

      endif

      integ = sum

      end

