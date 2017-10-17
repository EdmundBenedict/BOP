 
!****************************************************************************
!                                                                           *
!                                                                           *
!                             Alexey Girshick                               *
!              Department of Materials Science and Engineering              *
!                        University of Pennsylvania                         *
!                                                                           *
!                           +DLP, QUB                                       *
!****************************************************************************
!                                                                           *
!                                                                           *
!          This subroutine calculates the maximum force acting on           *
!                   atom for the dislocation relaxation.                    *
!                                                                           *
!                                                                           *
!****************************************************************************


   subroutine maxf_ds( maxf )
      use mod_precision
      use mod_all_scalar
      use mod_const
      use topologia

      implicit none


      include "Include/Force.array"
      include "Include/BEC.array"
      include "Include/PosVel.array"
      include "Include/ag.conn"

      real(dp) ::  flayer( 3 ),f_sq
      real(dp) ::  maxf
      real(dp) ::  disorg2
! keep is array because maxloc insists on it
      integer           i, j, k, keep(1), stat
!      INTEGER GFBCON


      maxf = 0.d0
      keep = 0
!
!      CALL FNDREC(8,'GFBCON',STAT)
!      READ(8,*) GFBCON
!      
!

!  Search for the maximum force among the individual atom forces.
      
!       if (gfbcon  ==  0) then
!          do i = 1, nd
!             f_sq = 0.d0
!             do j = 1, 3
!                f_sq = f_sq + ftot( j, i )*ftot( j, i )
!             enddo
!             
!             if ( sqrt( f_sq )  >  maxf ) then
!                maxf = sqrt( f_sq )
!                keep = i
!             endif
!          enddo
!       endif
!       
      if (gfbcon  ==  0) then
         
         
         keep = maxloc(sum(ftot*ftot, dim = 1))
         maxf = sum(ftot(:,keep(1))*ftot(:,keep(1)))
         
         
      else if (gfbcon  ==  1) then
      
!
!     If the GFBCs have been switched on, we only look for the maximum
!     force in the region R1
!
         do i = 1, nd
            disorg2 = sqrt(ad(1,i)*ad(1,i)+ad(2,i)*ad(2,i))
            if (disorg2  <=  radiusi) then
               f_sq = sum(ftot( 1:3, i )*ftot( 1:3, i ))
               
               if ( f_sq  >  maxf ) then
                  maxf = f_sq 
                  keep(1) = i
               endif
            endif
         enddo
      endif

      maxf = sqrt(maxf)
      
      if (iproc == master) write( 72, '(" Maximum force is ",E14.7," acting on atom ",I5/)')   maxf, keep(1)

      return
      end

