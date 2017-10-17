 
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
!            This subroutine sorts out the atoms in every layer.            *
!            (Z=const) according to their X and Y coordinates -             *
!            - it is necessary for using the sub-block forces.              *
!                                                                           *
!                                                                           *
!****************************************************************************


      subroutine sort_out()
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
      include "Include/Misc.array"
      include "Include/PosVel.array"

!
!    Common blocks introduced by A.Girshick.
!    
 
      include "Include/ag.conn"



      real(dp), parameter ::  delta = 1.0e-6_dp
      integer           i, j, k, ii, check, layer

!       parameter ( delta = 1.0e-6_dp )


!       do i = 1, 100
!           finlay( i ) = 0
!           ninlay( i ) = 0
!       enddo
!       
      finlay( 1:100 ) = 0
      ninlay( 1:100 ) = 0


!  Rearrange atoms according to the increase of Z.

      check = 1
      do while ( check  ==  1)
          check=0
          do i = 1, nd-1
              if ( d( 3, i )  >   d( 3, i+1 ) ) then
                  call swap( i, i+1 )
                  check = 1
              endif
          enddo
      enddo


!  Calculate the number of layers and fill out the layer arrays

      layer = 1
      finlay( 1 ) = 1
      ninlay( 1 ) = ninlay( 1 ) + 1
      lay( 1 ) = 1 
      do i = 2, nd
          if ( d( 3, i )  >  ( d( 3, i-1 ) + delta ) ) then
              layer = layer + 1
              finlay( layer ) = i
          endif
          ninlay( layer ) = ninlay( layer ) + 1
          lay( i ) = layer
      enddo

      num_of_layers = layer
      if ( num_of_layers  >  100 ) then
          write( 6, '(/"Increase the allowed number of layers."/)' )
          stop
      endif


!  Sort our atoms within the layers in the order of increasing X and Y.

      do k = 1, num_of_layers
          check = 1
          do while ( check  ==  1) 
              check=0
              do i = finlay( k ), finlay( k ) + ninlay( k ) - 2
                  if ( d( 1, i )  >   (d( 1, i+1 ) + delta) ) then
                      call swap( i, i+1 )
                      check = 1
                  else if ( (abs( d(1,i)-d(1,i+1) )  <  delta) &
     &               .and.(d( 2, i )  >  (d( 2, i+1 ) + delta))) then
                      call swap( i, i+1 )
                      check = 1
                  endif
              enddo
          enddo
      enddo


      return
      end

!_________________________________________________________________________
!

      subroutine swap( atom_1, atom_2 )
          use mod_precision


          use mod_all_scalar

          use mod_const

!
!    This is a routine to change the numbers in the list for 2 atoms.
!

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
      include "Include/Misc.array"
      include "Include/PosVel.array"


      integer           atom_1, atom_2, zdummy, j
      real(dp) ::  dummy
      character*2       symdummy


      do j = 1, 3
          dummy = d( j, atom_1 ) 
          d( j, atom_1 ) = d( j, atom_2 )
          d( j, atom_2 ) = dummy
      enddo

      zdummy = z( atom_1 )
      z( atom_1 ) = z( atom_2 )
      z( atom_2 ) = zdummy

      symdummy = symb( atom_1 )
      symb( atom_1 ) = symb( atom_2 )
      symb( atom_2 ) = symdummy


      return
      end
