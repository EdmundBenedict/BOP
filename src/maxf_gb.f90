 
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
!        This subrroutine calculates the individual atomic forces           *
!        (adjusted for the sub-block movement) and finds the                *
!        maximum force acting on atom among both sub-block and              *
!        individual atomic forces for the grain boundary relaxation.        *
!                                                                           *
!                                                                           *
!****************************************************************************


      subroutine maxf_gb( maxf )
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

      include "Include/Force.array"

!
!    Common blocks introduced by A.Girshick.
!    
 
      include "Include/ag.conn"

!
!    Declare the simple variables.
!

      real(dp) ::  flayer( 3 ),f_sq
      real(dp) ::  fmax1, fmax2, maxf
      integer           i, j, k, keep1, keep2


      fmax1 = 0.d0
      fmax2 = 0.d0
      keep1 = 0
      keep2 = 0


!  If we perform the relaxation at constant volume, then we shouldn't
!  account for the sub-block forces in Z direction.

      if ( const_volume  ==  1 ) then
          do i = 1, nd+1           
              fblock( 3, i ) = 0.d0
          enddo
      endif      


!  Calculate forces on individual atoms.
              
      do  k = 1, num_of_layers

          if ( ninlay( k )  ==  1 ) then

              if ( const_volume  ==  1 ) then
                  fatom( 1, finlay(k) ) = 0.d0
                  fatom( 2, finlay(k) ) = 0.d0
                  fatom( 3, finlay(k) ) = ftot( 3, finlay(k) ) 
              else
                  fatom( 1, finlay(k) ) = 0.d0
                  fatom( 2, finlay(k) ) = 0.d0
                  fatom( 3, finlay(k) ) = 0.d0
              endif

          else

!             Average forces acting on the atoms in layer K
!             except the first atom in the layer.
          
              do  j = 1, 3
                  flayer( j ) = ( fblock( j, finlay( k ) + 1 ) - &
     &                 fblock( j, finlay( k ) + ninlay( k ) ) ) / &
     &                 ( ninlay( k ) - 1 )
              enddo

!             Substract the average force acting on layer
!             from the force acting on atom.
 
              do  i = finlay(k) + 1, finlay(k) + ninlay(k) - 1
                  do  j = 1, 3
                      fatom( j, i ) = ftot( j, i ) - flayer( j )
                  enddo      
              enddo

          endif
      enddo


!  Search for the maximum force among the sub-block forces.

      do k = 1, num_of_layers
          i = finlay( k )
          f_sq = 0.d0
          do j = 1, 3
              f_sq = f_sq + fblock( j, i )*fblock( j, i )
          enddo
          f_sq = f_sq / ninlay( k ) / ninlay( k )
          
          if ( sqrt( f_sq )  >  fmax1 ) then
              fmax1 = sqrt( f_sq )
              keep1 = k
          endif
      enddo


!  Search for the maximum force among the individual atom forces.

      do i = 1, nd
          f_sq = 0.d0
          do j = 1, 3
              f_sq = f_sq + fatom( j, i )*fatom( j, i )
          enddo
          
          if ( sqrt( f_sq )  >  fmax2 ) then
              fmax2 = sqrt( f_sq )
              keep2 = i
          endif
      enddo


!  Compare those two maximum forces.

      if ( fmax1  >  fmax2 )  then
          maxf = fmax1
          write( 72, '(" Maximum force is ",E14.7," acting on " &
     &           "layer ",I2/ )') maxf, keep1
      else
          maxf = fmax2
          write( 72, '(" Maximum force is ",E14.7," acting on "  &
     &           "atom ",I2," out of ",I2," in layer ",I2/)') &
     &           maxf, keep2-finlay(lay(keep2))+1, &
     &           ninlay(lay(keep2)), lay(keep2)
      endif

      call flush( 72 )


      return
      end

