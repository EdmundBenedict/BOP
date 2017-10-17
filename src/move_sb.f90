 
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
!         This is an atom moving subroutine for the grain boundary,         *
!         stacking fault or gamma-surface relaxation - it involves          *
!                         the sub-block forces.                             *
!                                                                           *
!                                                                           *
!****************************************************************************


      subroutine move_sb()
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

      include "Include/PosVel.array"
      include "Include/Force.array"

!
!    Common blocks introduced by A.Girshick.
!    
 
      include "Include/ag.conn"

!
!    Declare the simple variables and arrays.
!

      real(dp) ::  maxdad
      real(dp) ::  dis( 3 ), totdis( 3 ), displ
      integer           i, j, k

!
!    Define parameters.
!

      parameter( maxdad = 0.20_dp )


      do j = 1, 3
          totdis( j ) = 0.d0
      enddo


!  Displace the atoms according to the sub-block forces.

      do k = 1, num_of_layers
          do j = 1, 3
              dis( j ) = fblock( j, finlay( k ) ) * ok
          enddo 

          displ = sqrt(dis(1)*dis(1) + dis(2)*dis(2) + dis(3)*dis(3)) 
          if ( displ  >  maxdad ) then
              write(6,'("MAXIMUM DISPLACEMENT EXCEEDED: ",F16.8)') displ
              do j = 1, 3
                  dis( j ) = maxdad * dis( j ) / displ 
              enddo
          endif
 
          do i = finlay( k ), nd
              do j = 1, 3
                  ad( j, i ) = ad( j, i ) + dis( j )                 
              enddo
              lef = lef + dis( 1 ) * dmdl( 1, i )    &
     &                  + dis( 2 ) * dmdl( 2, i )    &
     &                  + dis( 3 ) * dmdl( 3, i )   
          enddo

          if ( mvflag  ==  15  .or.  pointsgs  ==  1 )  then
              do i = finlay( k ), nd
                  do j = 1, 3
                      shcrd( j, i ) = shcrd( j, i ) + dis( j )
                  enddo
              enddo
          endif

          do j = 1, 3
              totdis( j ) = totdis( j ) + dis( j )
          enddo
      enddo

 
!  To calculate the total displacement of atoms above the relaxed region,
!  we have to add the displacement of the first unrelaxed layer above
!  the relaxed region.

      do j = 1, 3
          dis( j ) = fblock( j, nd+1 ) * ok  
      enddo
      
      displ = sqrt(dis(1)*dis(1) + dis(2)*dis(2) + dis(3)*dis(3))
      if ( displ  >  maxdad )  then
          do j = 1, 3
              dis( j ) = maxdad * dis( j ) / displ
          enddo
      endif
     
      do j = 1, 3
          totdis( j ) = totdis( j ) + dis( j )
      enddo


!  Displace all the atoms above the relaxed region according to the
!  total displacement.

      do i = 1, ninert
          if ( adinert( 3, i )  >  zeroheight ) then
              do j = 1, 3
                  adinert( j, i ) = adinert( j, i ) + totdis( j )
              enddo
          endif
      enddo

 
!  Displace the atoms according to the individual atomic forces.

!      DO  K = 1, NUM_OF_LAYERS
!          DO I = FINLAY( K ), FINLAY( K ) + NINLAY( K ) - 1
!              DO J = 1, 3
!                  DIS( J ) = FATOM( J, I ) * OK  
!              ENDDO
!
!              DISPL = DSQRT(DIS(1)*DIS(1)+DIS(2)*DIS(2)+DIS(3)*DIS(3))
!              IF ( DISPL .GT. MAXDAD )  THEN
!                  DO J = 1, 3
!                      DIS( J ) = MAXDAD * DIS( J ) / DISPL
!                  ENDDO
!              ENDIF
!
!              DO J = 1, 3 
!                  AD( J, I ) = AD( J, I ) + DIS( J )
!              ENDDO 
!              LEF = LEF + DIS( 1 ) * DMDL( 1, I )   
!     &                  + DIS( 2 ) * DMDL( 2, I )   
!     &                  + DIS( 3 ) * DMDL( 3, I )   
!
!              IF ( MVFLAG .EQ. 15  .OR.  POINTSGS .EQ. 1 )  THEN
!                  DO J = 1, 3
!                      SHCRD( J, I ) = SHCRD( J, I ) + DIS( J )
!                  ENDDO
!              ENDIF
!
!          ENDDO
!      ENDDO


      return
      end

	
