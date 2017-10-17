 
       function msd()
          use mod_precision


          use mod_all_scalar

          use mod_const

!
!    This is a function to evaluate the mean squared displacement
!

      implicit none
      real(dp) :: msd

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
!      include "Include/Misc.array"

      integer i
      real(dp) :: rt

!
!    Evaluate the order parameter.
!

      msd = 0.0_dp

      do i = 1, nd, 1

         rt = sqrt(adnopb(1,i)**2 + adnopb(2,i)**2 + adnopb(3,i)**2)
         msd = msd + (rt-dispo(i))**2

      enddo

      msd = msd/nd

      end

