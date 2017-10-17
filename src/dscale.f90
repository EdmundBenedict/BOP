 
       function dscale(r,r0,rc,n,nc,r1,rcut,c0,c1,c2, c3,c4,c5)

!
!    This is a function that evaluates the scaling function.
!    NOTE: N=0 implies just an exponential decay
!          NC = 0 implies just a power law decay
!
          use mod_precision
          use mod_gsp, only : gsp, dgsp
      implicit none
      real(dp) :: dscale

!
!    Declare the simple variables.
!

      real(dp) :: r,r0,rc,n,nc
      real(dp) :: r1,rcut,c0,c1,c2,c3,c4,c5
      real(dp) :: y

!
!    Evaluate the function.
!

      if (r < r1) then
         dscale = dgsp(r,r0,rc,n,nc)
      elseif (r < rcut) then
         y = (r-r1)/(rcut-r1)
!          dscale = dgsp(r,r0,rc,n,nc) * (((-6.0_dp*y + 15.0_dp)*y -10.0_dp)*y*y*y + 1.0_dp  ) &
         dscale = dgsp(r,r0,rc,n,nc) * ((6.0_dp*y + 3.0_dp)*y + 1.0_dp) * (1.0_dp - y)**3 &
               & + gsp(r,r0,rc,n,nc) * (-30.0_dp*(y*( y  -  1.0_dp ))**2)/(rcut-r1)
         
!          dscale = c1+y*(2.0_dp*c2+y*(3.0_dp*c3+y*(4.0_dp*c4+y*5.0_dp*c5)))
         
      else
         dscale = 0.0_dp
      endif

      end

