 
       function vscale(r,r0,rc,n,nc,r1,rcut,c0,c1,c2,c3,c4,c5)

!
!    This is a function that evaluates the scaling function.
!    NOTE: N=0 implies just an exponential decay
!          NC = 0 implies just a power law decay
!
          use mod_precision
          use mod_gsp, only : gsp
          
      implicit none
      real(dp) :: vscale

!
!    Declare the simple variables.
!

      real(dp) :: r,r0,rc,n,nc
      real(dp) :: r1,rcut,c0,c1,c2,c3,c4,c5
      real(dp) :: y
!
!    Evaluate the function.
!

!      write(*,*) r0
      if (r0 == 0) then
         write(6,'(''Undefined bond type requested.'')')
         call panic()
      endif


!       print *, r

      if (r < r1) then
         vscale = gsp(r,r0,rc,n,nc)
      elseif (r < rcut) then
!    print *, 'hop in tail'
         y = (r-r1)/(rcut-r1)
         vscale = gsp(r,r0,rc,n,nc) * ((6.0_dp*y + 3.0_dp)*y + 1.0_dp) * (1.0_dp - y)**3
                                    ! (((-6.0_dp*y + 15.0_dp)*y -10.0_dp)*y*y*y + 1.0_dp  )
!          vscale = c0 + y*(c1 + y*(c2 + y*(c3+y*(c4+y*c5))))
      else
         vscale = 0.0_dp
      endif

      end

