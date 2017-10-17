 
    module mod_gsp
          use mod_precision
    implicit none
    
    contains
 
       pure function gsp(r,r0,rc,n,nc)
    


!
!    This is a function that evaluates the Goodwin-Skinner-Pettifor
!     scaling function.
!    NOTE: N=0 implies just an exponential decay
!          NC = 0 implies just a power law decay
!

      
      real(dp) :: gsp

!
!    Declare the simple variables.
!

      real(dp), intent(in) :: r,r0,rc,n,nc
      real(dp) :: f1,f2

!
!    Evaluate the function.
!

      if (n == 0.0_dp) then
         f1 = 1.0_dp
         if (nc == 0.0_dp) then
            f2 = 1.0_dp
         elseif (nc == 1.0_dp) then
            f2 = exp((r0/rc)-(r/rc))
         else
            f2 = exp((r0/rc)**nc-(r/rc)**nc)
         endif
      else
!          f1 = exp(n*log(r0/r))
         f1 = (r0/r)**n
         if (nc == 0.0_dp) then
            f2 = 1.0_dp
         elseif (nc == 1.0_dp) then
            f2 = exp(n*((r0/rc)-(r/rc)))
         else
            f2 = exp(n*((r0/rc)**nc-(r/rc)**nc))
         endif
      endif

      gsp = f1*f2
!       
!       gsp = exp(n*log(r0/r))
!       gsp = (r0/r)**n * exp(n*((r0/rc)**nc - (r/rc)**nc))
        
        
!        gsp = exp(n*(log(r0/r) + exp(nc*log(r0/rc)) - exp(nc*log(r/rc)) ))        
        
        
!       gsp = exp(n*log(r0/r)) * exp( n*( exp(nc*log(r0/rc)) - exp(nc*log(r/rc)) )  )
      

      end function gsp

 
       pure function dgsp(r,r0,rc,n,nc)


!
!    This is a function that evaluates the derivative of the
!     Goodwin-Skinner-Pettifor scaling function.
!    NOTE: N=0 implies just an exponential decay
!          NC = 0 implies just a power law decay
!

      real(dp) :: dgsp

!
!    Declare the simple variables.
!

      real(dp), intent(in) :: r,r0,rc,n,nc

      real(dp) :: f1,f2

!
!    Evaluate the function.
!

      f1 = gsp(r,r0,rc,n,nc)

      if (n == 0.0_dp) then
         if (nc == 0.0_dp) then
            f2 = 0.0_dp
         elseif (nc == 1.0_dp) then
            f2 = -1.0_dp/rc
         else
            f2 = -(nc/r)*((r/rc)**nc)
         endif
      else
         if (nc == 0.0_dp) then
            f2 = -n/r
         elseif (nc == 1.0_dp) then
            f2 = -n*(1.0_dp/r + 1.0_dp/rc)
         else
            f2 = -(n/r)*(1.0_dp+nc*(r/rc)**nc)
         endif
      endif

      dgsp = f1*f2

     end function dgsp
        
        
        
       pure function d2gsp(r,r0,rc,n,nc)


!
!    This is a function that evaluates the second derivative of the
!     Goodwin-Skinner-Pettifor scaling function.
!    NOTE: N=0 implies just an exponential decay
!          NC = 0 implies just a power law decay
!

      real(dp) :: d2gsp

!
!    Declare the simple variables.
!

      real(dp), intent(in) :: r,r0,rc,n,nc

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

      end function d2gsp



      

      
      
      
      






    end module mod_gsp









