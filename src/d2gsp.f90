 
      function d2gsp(r,r0,rc,n,nc)
          use mod_precision


!
!    This is a function that evaluates the second derivative of the
!     Goodwin-Skinner-Pettifor scaling function.
!    NOTE: N=0 implies just an exponential decay
!          NC = 0 implies just a power law decay
!

      implicit none
      real(dp) :: d2gsp

!
!    Declare the simple variables.
!

      real(dp), intent(in) :: r,r0,rc,n,nc
      real(dp) :: gsp
      real(dp) :: f1,f2

!
!    Evaluate the function.
!

      f1 = gsp(r,r0,rc,n,nc)

      if (n == 0.0_dp) then
         if (nc == 0.0_dp) then
            f2 = 0.0_dp
         elseif (nc == 1.0_dp) then
            f2 = 1.0_dp/rc**2
         else
            f2 = (nc/r**2)*(r/rc)**nc*(nc*(r/rc)**nc-nc+1)
         endif
      else
         if (nc == 0.0_dp) then
            f2 = n*(n+1)/r**2
         else

            f2 = n/r**2 * ( 1 + n*(1+nc*(r/rc)**nc)**2 + &
     &                     nc*(r/rc)**nc - nc**2*(r/rc)**nc ) 

         endif
      endif

      d2gsp = f1*f2

      end

