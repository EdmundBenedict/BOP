 
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
!         This subroutine minimizes the energy with respect to the          *
!         two independent lattice parameters A and C (L10, HCP...)          *
!                                                                           *
!                                                                           *
!****************************************************************************



      subroutine  findminhcp()
          use mod_precision


          use mod_all_scalar


      implicit  none

!
!    Include constants.
!

!      include "Include/ALL.const"

!
!    Include scalars.
!

!      include "Include/ALL.scalar"


      integer           i,  imin
      real(dp) ::  energy( 9 )
      real(dp) ::  steph, min_steph, max_steph
      real(dp) ::  alpha, beta, energy_min


      open( unit = 82, file='min')
      open( unit = 80, file='min.l')
      
      max_steph = 1.0e-2_dp
      min_steph = 1.0e-8_dp

      alpha = 1.d0
      beta  = 1.d0
      steph = max_steph

 801  format( 2x, f12.10, 5x, f12.10, 5x, f17.12 )
 802  format( / 5x, f12.10, 5x, f12.10, 5x, f17.12 / )


 100  continue
    
      call get_en( alpha - steph,  beta + steph, energy( 1 ) )
      call get_en( alpha       ,  beta + steph, energy( 2 ) )
      call get_en( alpha + steph,  beta + steph, energy( 3 ) )

      call get_en( alpha - steph,  beta       , energy( 4 ) )
      call get_en( alpha       ,  beta       , energy( 5 ) )
      call get_en( alpha + steph,  beta       , energy( 6 ) )

      call get_en( alpha - steph,  beta - steph, energy( 7 ) )
      call get_en( alpha       ,  beta - steph, energy( 8 ) )
      call get_en( alpha + steph,  beta - steph, energy( 9 ) )

      write( 80, 801 ) alpha - steph,  beta + steph, energy( 1 )
      write( 80, 801 ) alpha       ,  beta + steph, energy( 2 ) 
      write( 80, 801 ) alpha + steph,  beta + steph, energy( 3 )

      write( 80, 801 ) alpha - steph,  beta       , energy( 4 )
      write( 80, 801 ) alpha       ,  beta       , energy( 5 )
      write( 80, 801 ) alpha + steph,  beta       , energy( 6 )

      write( 80, 801 ) alpha - steph,  beta - steph, energy( 7 )
      write( 80, 801 ) alpha       ,  beta - steph, energy( 8 ) 
      write( 80, 801 ) alpha + steph,  beta - steph, energy( 9 )


 200  continue

!  Determine which of the nine points has the lowest energy.

      energy_min = energy( 1 )
      imin = 1
      do i = 1, 9
          if( energy( i )  <  energy_min ) then
              energy_min = energy( i )
              imin = i
          end if
      enddo


!  According to which point has the lowest energy, recalculate
!  only the points which were not calculated before.

      if      ( imin  ==  1 )  then

          alpha = alpha - steph
          beta  = beta  + steph
          
          write( 82, 801 )   alpha, beta, energy_min
          write( 80, 802 )   alpha, beta, energy_min
          call flush( 80 )

          energy( 9 ) = energy( 5 )
          energy( 5 ) = energy( 1 )
          energy( 6 ) = energy( 2 )
          energy( 8 ) = energy( 4 )
          call get_en( alpha - steph,  beta + steph, energy( 1 ) )
          call get_en( alpha       ,  beta + steph, energy( 2 ) )
          call get_en( alpha + steph,  beta + steph, energy( 3 ) )
          call get_en( alpha - steph,  beta       , energy( 4 ) )
          call get_en( alpha - steph,  beta - steph, energy( 7 ) )

          write( 80, 801 ) alpha-steph,  beta+steph, energy( 1 )
          write( 80, 801 ) alpha     ,  beta+steph, energy( 2 )
          write( 80, 801 ) alpha+steph,  beta+steph, energy( 3 )
          write( 80, 801 ) alpha-steph,  beta     , energy( 4 )
          write( 80, 801 ) alpha     ,  beta     , energy( 5 )
          write( 80, 801 ) alpha+steph,  beta     , energy( 6 )
          write( 80, 801 ) alpha-steph,  beta-steph, energy( 7 )
          write( 80, 801 ) alpha     ,  beta-steph, energy( 8 )
          write( 80, 801 ) alpha+steph,  beta-steph, energy( 9 )

          go to 200

      else if ( imin  ==  2 )  then
          
          alpha = alpha
          beta  = beta  + steph

          write( 82, 801 )   alpha, beta, energy_min
          write( 80, 802 )   alpha, beta, energy_min
          call flush( 80 )      

          energy( 7 ) = energy( 4 )
          energy( 8 ) = energy( 5 )
          energy( 9 ) = energy( 6 )
          energy( 4 ) = energy( 1 )
          energy( 5 ) = energy( 2 )
          energy( 6 ) = energy( 3 )
          call get_en( alpha - steph,  beta + steph, energy( 1 ) )  
          call get_en( alpha       ,  beta + steph, energy( 2 ) )  
          call get_en( alpha + steph,  beta + steph, energy( 3 ) )

          write( 80, 801 ) alpha-steph,  beta+steph, energy( 1 )    
          write( 80, 801 ) alpha     ,  beta+steph, energy( 2 )    
          write( 80, 801 ) alpha+steph,  beta+steph, energy( 3 )    
          write( 80, 801 ) alpha-steph,  beta     , energy( 4 )
          write( 80, 801 ) alpha     ,  beta     , energy( 5 )    
          write( 80, 801 ) alpha+steph,  beta     , energy( 6 )    
          write( 80, 801 ) alpha-steph,  beta-steph, energy( 7 )    
          write( 80, 801 ) alpha     ,  beta-steph, energy( 8 )
          write( 80, 801 ) alpha+steph,  beta-steph, energy( 9 )    
          
          go to 200

      else if ( imin  ==  3 )  then

          alpha = alpha + steph
          beta  = beta  + steph
          
          write( 82, 801 )   alpha, beta, energy_min
          write( 80, 802 )   alpha, beta, energy_min
          call flush( 80 )

          energy( 7 ) = energy( 5 )
          energy( 5 ) = energy( 3 )
          energy( 4 ) = energy( 2 )
          energy( 8 ) = energy( 6 )
          call get_en( alpha - steph,  beta + steph, energy( 1 ) )
          call get_en( alpha       ,  beta + steph, energy( 2 ) )
          call get_en( alpha + steph,  beta + steph, energy( 3 ) )
          call get_en( alpha + steph,  beta       , energy( 6 ) )
          call get_en( alpha + steph,  beta - steph, energy( 9 ) )
          
          write( 80, 801 ) alpha-steph,  beta+steph, energy( 1 )    
          write( 80, 801 ) alpha     ,  beta+steph, energy( 2 )    
          write( 80, 801 ) alpha+steph,  beta+steph, energy( 3 )    
          write( 80, 801 ) alpha-steph,  beta     , energy( 4 )
          write( 80, 801 ) alpha     ,  beta     , energy( 5 )    
          write( 80, 801 ) alpha+steph,  beta     , energy( 6 )    
          write( 80, 801 ) alpha-steph,  beta-steph, energy( 7 )    
          write( 80, 801 ) alpha     ,  beta-steph, energy( 8 )
          write( 80, 801 ) alpha+steph,  beta-steph, energy( 9 )    
          
          go to 200

      else if ( imin  ==  4 )  then

          alpha = alpha - steph
          beta  = beta
          
          write( 82, 801 )   alpha, beta, energy_min
          write( 80, 802 )   alpha, beta, energy_min
          call flush( 80 )

          energy( 3 ) = energy( 2 )
          energy( 2 ) = energy( 1 )
          energy( 6 ) = energy( 5 )
          energy( 5 ) = energy( 4 )
          energy( 9 ) = energy( 8 )
          energy( 8 ) = energy( 7 )
          call get_en( alpha - steph,  beta + steph, energy( 1 ) )
          call get_en( alpha - steph,  beta       , energy( 4 ) )
          call get_en( alpha - steph,  beta - steph, energy( 7 ) )

          write( 80, 801 ) alpha-steph,  beta+steph, energy( 1 )    
          write( 80, 801 ) alpha     ,  beta+steph, energy( 2 )    
          write( 80, 801 ) alpha+steph,  beta+steph, energy( 3 )    
          write( 80, 801 ) alpha-steph,  beta     , energy( 4 )
          write( 80, 801 ) alpha     ,  beta     , energy( 5 )    
          write( 80, 801 ) alpha+steph,  beta     , energy( 6 )    
          write( 80, 801 ) alpha-steph,  beta-steph, energy( 7 )    
          write( 80, 801 ) alpha     ,  beta-steph, energy( 8 )
          write( 80, 801 ) alpha+steph,  beta-steph, energy( 9 )    
          
          go to 200

      else if ( imin  ==  5 )  then

          alpha = alpha
          beta  = beta
          
          write( 82, 801 )   alpha, beta, energy_min
          write( 80, 802 )   alpha, beta, energy_min
          call flush( 80 )


          steph = steph / 10.d0
          if( steph  <  min_steph )  then
              write(  6, 801 )   alpha, beta, energy_min
              write( 82, 801 )   alpha, beta, energy_min
              write( 80, 802 )   alpha, beta, energy_min

              write( 80, '( "EBOND = ", F9.6 )' )    ebond / nd
              write( 80, '( "EPAIR = ", F9.6 )' )    epair / nd            
              write( 80, '( "EPROM = ", F9.6 )' )    eprom / nd        
              write( 80, '( "EATOM = ", F9.6 )' )  - eatom / nd         
              write( 80, '( "UENT  = ", F9.6 )' )  - uent  / nd        

              call tag( 80 )

              close( 82 )
              close( 80 )

              return

          end if

          go to 100

      else if ( imin  ==  6 )  then

          alpha = alpha + steph
          beta  = beta
          
          write( 82, 801 )   alpha, beta, energy_min
          write( 80, 802 )   alpha, beta, energy_min
          call flush( 80 )

          energy( 1 ) = energy( 2 )
          energy( 2 ) = energy( 3 )
          energy( 4 ) = energy( 5 )
          energy( 5 ) = energy( 6 )
          energy( 7 ) = energy( 8 )
          energy( 8 ) = energy( 9 )
          call get_en( alpha + steph,  beta + steph, energy( 3 ) )
          call get_en( alpha + steph,  beta       , energy( 6 ) )
          call get_en( alpha + steph,  beta - steph, energy( 9 ) ) 

          write( 80, 801 ) alpha-steph,  beta+steph, energy( 1 )    
          write( 80, 801 ) alpha     ,  beta+steph, energy( 2 )    
          write( 80, 801 ) alpha+steph,  beta+steph, energy( 3 )    
          write( 80, 801 ) alpha-steph,  beta     , energy( 4 )
          write( 80, 801 ) alpha     ,  beta     , energy( 5 )    
          write( 80, 801 ) alpha+steph,  beta     , energy( 6 )    
          write( 80, 801 ) alpha-steph,  beta-steph, energy( 7 )    
          write( 80, 801 ) alpha     ,  beta-steph, energy( 8 )
          write( 80, 801 ) alpha+steph,  beta-steph, energy( 9 )    
          
          go to 200

      else if ( imin  ==  7 )  then

          alpha = alpha - steph
          beta  = beta  - steph
          
          write( 82, 801 )   alpha, beta, energy_min
          write( 80, 802 )   alpha, beta, energy_min
          call flush( 80 )

          energy( 2 ) = energy( 4 )
          energy( 3 ) = energy( 5 )
          energy( 5 ) = energy( 7 )
          energy( 6 ) = energy( 8 )
          call get_en( alpha - steph,  beta + steph, energy( 1 ) )
          call get_en( alpha - steph,  beta       , energy( 4 ) )
          call get_en( alpha - steph,  beta - steph, energy( 7 ) )
          call get_en( alpha       ,  beta - steph, energy( 8 ) )
          call get_en( alpha + steph,  beta - steph, energy( 9 ) )

          write( 80, 801 ) alpha-steph,  beta+steph, energy( 1 )    
          write( 80, 801 ) alpha     ,  beta+steph, energy( 2 )    
          write( 80, 801 ) alpha+steph,  beta+steph, energy( 3 )    
          write( 80, 801 ) alpha-steph,  beta     , energy( 4 )
          write( 80, 801 ) alpha     ,  beta     , energy( 5 )    
          write( 80, 801 ) alpha+steph,  beta     , energy( 6 )    
          write( 80, 801 ) alpha-steph,  beta-steph, energy( 7 )    
          write( 80, 801 ) alpha     ,  beta-steph, energy( 8 )
          write( 80, 801 ) alpha+steph,  beta-steph, energy( 9 )    
          
          go to 200

      else if ( imin  ==  8 )  then

          alpha = alpha
          beta  = beta  - steph
          
          write( 82, 801 )   alpha, beta, energy_min
          write( 80, 802 )   alpha, beta, energy_min
          call flush( 80 )

          energy( 1 ) = energy( 4 )
          energy( 2 ) = energy( 5 )
          energy( 3 ) = energy( 6 )
          energy( 4 ) = energy( 7 )
          energy( 5 ) = energy( 8 )
          energy( 6 ) = energy( 9 )
          call get_en( alpha - steph,  beta - steph, energy( 7 ) )
          call get_en( alpha       ,  beta - steph, energy( 8 ) )  
          call get_en( alpha + steph,  beta - steph, energy( 9 ) )

          write( 80, 801 ) alpha-steph,  beta+steph, energy( 1 )    
          write( 80, 801 ) alpha     ,  beta+steph, energy( 2 )    
          write( 80, 801 ) alpha+steph,  beta+steph, energy( 3 )    
          write( 80, 801 ) alpha-steph,  beta     , energy( 4 )
          write( 80, 801 ) alpha     ,  beta     , energy( 5 )    
          write( 80, 801 ) alpha+steph,  beta     , energy( 6 )    
          write( 80, 801 ) alpha-steph,  beta-steph, energy( 7 )    
          write( 80, 801 ) alpha     ,  beta-steph, energy( 8 )
          write( 80, 801 ) alpha+steph,  beta-steph, energy( 9 )    

          go to 200          

      else if ( imin  ==  9 )  then

          alpha = alpha + steph
          beta  = beta  - steph
          
          write( 82, 801 )   alpha, beta, energy_min
          write( 80, 802 )   alpha, beta, energy_min
          call flush( 80 )

          energy( 1 ) = energy( 5 )
          energy( 2 ) = energy( 6 )
          energy( 4 ) = energy( 8 )
          energy( 5 ) = energy( 9 )
          call get_en( alpha + steph,  beta + steph, energy( 3 ) )
          call get_en( alpha + steph,  beta       , energy( 6 ) )          
          call get_en( alpha - steph,  beta - steph, energy( 7 ) )
          call get_en( alpha       ,  beta - steph, energy( 8 ) )
          call get_en( alpha + steph,  beta - steph, energy( 9 ) )

          write( 80, 801 ) alpha-steph,  beta+steph, energy( 1 )    
          write( 80, 801 ) alpha     ,  beta+steph, energy( 2 )    
          write( 80, 801 ) alpha+steph,  beta+steph, energy( 3 )    
          write( 80, 801 ) alpha-steph,  beta     , energy( 4 )
          write( 80, 801 ) alpha     ,  beta     , energy( 5 )    
          write( 80, 801 ) alpha+steph,  beta     , energy( 6 )    
          write( 80, 801 ) alpha-steph,  beta-steph, energy( 7 )    
          write( 80, 801 ) alpha     ,  beta-steph, energy( 8 )
          write( 80, 801 ) alpha+steph,  beta-steph, energy( 9 )    
          
          go to 200

      end if


      go to 100


      return
      end
