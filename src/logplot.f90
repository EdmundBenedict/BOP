 
!****************************************************************************
!                                                                           *
!                                                                           *
!                             Alexey Girshick                               *
!              Department of Materials Science and Engineering              *
!                        University of Pennsylvania                         *
!                                                                           *
!                                                                           *
!****************************************************************************
!                                                                           *
!                                                                           *
!      This subroutine calculate the logariphmic radial dependence of       *
!                        the dislocation self-energy.                       *
!                                                                           *
!                                                                           *
!****************************************************************************


      subroutine logplot( )
          use mod_precision


          use mod_all_scalar

          use mod_const


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
 
!      include "Include/Atom.array"
!      include "Include/Misc.array"
!      include "Include/PosVel.array"

!
!    Common blocks introduced by A.Girshick.
!    
 
      include "Include/ag.conn"


      real(dp) ::  eout
      integer flag

      flag = 1
      call getetot(flag)
      eout = eprom + ebond + epair - eatom - uent
      write ( 84, '(F10.6, 4X, E14.7)') &
     &          log_rad, eout + ( nd * e_coh )
      write ( 85, '(F10.6, 4X, E14.7)')  log_rad, nd*1.d0
      call flush( 84 )
      call flush( 85 )

      return
      end
