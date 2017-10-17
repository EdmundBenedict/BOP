 
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
!         This subroutine relaxes an atomic configuration for the           * 
!           grain boundary, stacking fault or one point on the              *
!               gamma-surfaces (using the sub-block forces).                *
!                                                                           *
!                                                                           *
!****************************************************************************


      subroutine relax_sb( ezero, surf )
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
 
      include "Include/Atom.array"
!      include "Include/PosVel.array"

!
!    Common blocks introduced by A.Girshick.
!    
 
      include "Include/ag.conn"

!
!    Declare the simple variables.
!

      real(dp) :: eout, u0
      real(dp) :: maxf
      real(dp) :: ezero, surf

      integer usrexit
      integer i, j, k, flag
      character*80 filename

!      FILENAME = GENFILE(1:LENGFN)//'.72'
!      OPEN(UNIT = 72,FILE = FILENAME,STATUS = 'NEW')

      write( 6, '("*** RELAXATON ***")' )

      eflo = 0.d0
      do i = 1, nd, 1
          de( i ) = 0.d0
      enddo

      flag = 1
      call getetot(flag)

      call erasab()
      eout = eprom + ebond + epair - eatom - uent
      u0   = eout + 0.5_dp*uent

      if ( mvflag  ==  15  .or.  pointsgs  ==  1 ) then
         write( 72, '("Energy: ",F12.6)' ) ((eout/surf) - ezero)
      endif
      call flush( 72 ) 

      if      ( mvflag  ==  14 ) then
         call maxf_gs( maxf )
      else if ( mvflag  ==  15  .and.  gbmoveopt  ==  1 ) then
         call maxf_gs( maxf )
      else if ( mvflag  ==  15  .and.  gbmoveopt  ==  2 ) then
         call maxf_gb( maxf )
      endif


!  Iteration loop.

      do while ((iter <= mxiter).and.(maxf > ftol) &
     &                          .and.(usrexit() == 0) )
       
!        Move the atoms.

         call move_sb( ) 

!        Recalculate the energy.

         if (dndm > 1.0e-6_dp) lef = lef + (locc-totnia)/dndm      
         flag = 1

         call getetot(flag)

         eout = eprom + ebond + epair - eatom - uent
         u0   = eout  + 0.5_dp*uent
         if ( mvflag  ==  15  .or.  pointsgs  ==  1 ) then
            write( 72, '("Energy: ",F12.6)' ) ((eout/surf)-ezero)
         endif
         call flush( 72 )

         if      ( mvflag  ==  14 ) then
            call maxf_gs( maxf )
         else if ( mvflag  ==  15  .and.  gbmoveopt  ==  1 ) then
            call maxf_gs( maxf )
         else if ( mvflag  ==  15  .and.  gbmoveopt  ==  2 ) then
            call maxf_gb( maxf ) 
         endif

         iter = iter + 1

      enddo

      if     ( maxf  <  ftol )  then
         write( 72, '(''Relaxation successfully completed.'')')
      elseif ( iter  >  mxiter )  then
         write( 72, '(''Maximum number of iterations exceeded.'')')
      endif


!  Write down the displacements of the active atoms (if we were
!  doing grain boundary or considering only one point in the 
!  gamma-surface).

      if ( mvflag  ==  15  .or.  pointsgs  ==  1 )   then
          open (unit = 82,file = 'displ.out')
          write ( 82, '(//10X,"Displacements of the active atoms")' )
              
          do  k = 1, num_of_layers
              write( 82, '(//" LAYER ",I2,":"/)' ) k
              do  i = finlay(k), finlay(k) + ninlay(k) - 1
                  write( 82, '("Atom",I2,":",3(6X,F8.4))') &
     &              i + 1 - finlay( k ),  &
     &              ( shcrd( 1, i ) - unshad( 1, i ) ) / latpar ,  &
     &              ( shcrd( 2, i ) - unshad( 2, i ) ) / latpar , &
     &              ( shcrd( 3, i ) - unshad( 3, i ) ) / latpar
              enddo      
          enddo

          close ( 82 )
      endif

      close ( 72 ) 

      return
      end

