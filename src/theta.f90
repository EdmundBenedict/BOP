 
       function theta(x,kt)
          use mod_precision


!
!    This is a function to evaluate the fermi function.
!

      implicit none
      real(dp) :: theta

!
!    Declare the simple variables.
!

      real(dp) :: x,y,kt

!
!    Evaluate the function.
!

      y = x/kt
!       print '("y: ",f18.5)',y
! 
!       if (y.lt.-20.0d0) then
!          theta = 1.0d0
!       elseif (y.gt.20.0d0) then
!          theta = 0.0d0
!       else
!          theta = 1.0d0/(1.0d0 + exp(y))
!       endif
    
      if (y < -32.0_dp) then
         theta = 1.0_dp
      elseif (y > 32.0_dp) then
         theta = 0.0_dp
      else
         theta = 1.0_dp/(1.0_dp + exp(y))
      endif

      

!       
!       if (y.gt.20.0d0) then
!          theta = 0.0d0
!       elseif (y.ge.-20.0d0) then
!          theta = 1.0d0/(1.0d0 + exp(y))
!       else
!          theta = 1.0d0
!       endif
      
! 
!     if (.not. (abs(y) > 20.0d0)) then
!         theta = 1.0d0/(1.0d0 + exp(y))
!     else if (y.gt.20.0d0) then
!         theta = 0.0d0
!     else
!         theta = 1.0d0
!     end if

!       if (-20<=y .and. y>=20) then 
!           theta = 1.0d0/(1.0d0 + exp(y))
!       else if (y>20) then
!           theta = 0.0d0
!       else 
!           theta = 1.0d0
!       end if

      end

