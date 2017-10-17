 
       function orderp()
          use mod_precision


          use mod_all_scalar

          use mod_const

!
!    This is a function to evaluate the
!     translational order parameter (cf Allen & Tildesley p.171)
!

      implicit none
      real(dp) :: orderp

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
      include "Include/Misc.array"

      integer i
      real(dp) :: kx,ky,kz,kdotr

!
!    Evaluate the order parameter.
!

      kx = (2.0_dp*pi/latparam)*kcpt(1)
      ky = (2.0_dp*pi/latparam)*kcpt(2)
      kz = (2.0_dp*pi/latparam)*kcpt(3)

      orderp = 0.0_dp

      do i = 1, nd

         kdotr  = kx*ad(1,i) + ky*ad(2,i) + kz*ad(3,i)
         orderp = orderp + cos(kdotr)

      enddo

      orderp = orderp/nd

      end

