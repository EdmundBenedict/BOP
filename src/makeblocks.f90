 
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
!         This subroutine creats atomic blocks corresponding to the         *
!         different size of the active part, needed to calculate the        * 
!       logariphmic radial dependence of the dislocation's self-energy.     *
!                                                                           *
!                                                                           *
!****************************************************************************
!
!     Modified by M. Cawkwell as we now have the dislocation line
!     parallel to z and XXMIN etc. are not used anymore. 29th July 2003
!

      subroutine makeblocks( )
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

!
!    Declare the variables.
!

      real(dp) ::  ad_0 ( 3, mxnd )
      real(dp) ::  adinert_0( 3, minert ) 
      real(dp) ::  unrld_0( 3, mxnd )
      real(dp) ::  de_0(mxnd)
      real(dp) ::  eflo_0
      character*2       symb_0( mxnd )
      character*2       symi_0( minert )
      integer           z_0( mxnd )
      integer           zinert_0( minert )
      integer           nd_0, ninert_0

      real(dp) ::  distance
      integer           i, j, ii, jj
!
      x_center = ( x_center ) * latpar
      y_center = ( y_center ) * latpar
      r_min = r_min * latpar
      r_max = r_max * latpar

!  Record the complete relaxed block information.

      do i = 1, nd
          do j = 1, 3
              ad_0( j, i )    = ad( j, i )
              unrld_0( j, i ) = unrld( j, i )
          enddo 
          de_0( i )   = de( i )
          symb_0( i ) = symb( i )
          z_0( i ) = z( i )
      enddo

      do i = 1, ninert
          do j = 1, 3
              adinert_0( j, i )    = adinert( j, i )
          enddo
          symi_0( i ) = symi( i )
          zinert_0( i ) = zinert( i )
      enddo

      nd_0     = nd
      ninert_0 = ninert
      eflo_0   = eflo


!  Main cycle: make blocks for all the required radii.

      do klp = 1, n_mesh

!         Renull the arrays we have to fill out.

          do i = 1, nd_0
              do j = 1, 3
                  ad( j, i )    = 0.0_dp
                  unrld( j, i ) = 0.0_dp
              enddo
              de( i ) = 0.0_dp
              z( i )  = 0
          enddo
          do i = 1, (nd_0 + ninert_0)
              do j = 1, 3
                  adinert( j, i ) = 0.0_dp
              enddo
              zinert( i ) = 0
          enddo
          nd     = 0
          ninert = 0
          ii = 0
          jj = 0

          radius = r_min + ( r_max - r_min )*(klp-1)/(n_mesh-1)

!         Fill in the arrays going through all the active atoms.

          do i = 1, nd_0
              distance = sqrt( ( ad_0(1,i) - x_center )**2 + &
     &                          ( ad_0(2,i) - y_center )**2 )
              if ( distance  <=  radius ) then
                  ii = ii + 1
                  do j = 1, 3
                      ad( j, ii )    = ad_0 ( j, i )
                      unrld( j, ii ) = unrld_0 ( j, i )
                  enddo
                  z( ii )    = z_0( i )
                  de( ii )   = de_0( i )
                  symb( ii ) = symb_0( i )
              else
                  jj = jj + 1
                  do j = 1, 3 
                      adinert( j, jj )    = ad_0 ( j, i )
                  enddo
                  symi( jj )   = symb_0( i )
                  zinert( jj ) = z_0( i )
              endif
          enddo

!         Add the inert atoms to the lit.

          do i = 1, ninert_0
              jj = jj + 1
              do j = 1, 3
                  adinert( j, jj ) = adinert_0 ( j, i )     
              enddo
              symi( jj ) = symi_0( i )   
              zinert( jj ) = zinert_0( i )
          enddo

          nd     = ii
          ninert = jj
          if( (nd+ninert)  /=  (nd_0+ninert_0) ) then
              write( 6, '(//"Error in the atom rearrangement " &
     &             "routine."//)' )
              stop
          endif
         
          eflo = eflo_0

!         Write down the block.

          call outblock( 2 )

      enddo


      return
      end

