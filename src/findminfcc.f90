 
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
!           This subroutine minimizes the energy with respect to            *
!              one independent lattice parameter (FCC, BCC...)              *
!                                                                           *
!                                                                           *
!****************************************************************************


   subroutine  findminfcc()
      use mod_precision
      use mod_all_scalar

      implicit  none


      integer  ::  i,  imin
      real(dp) ::  energy( 3 )
      real(dp) ::  stepf, min_stepf, max_stepf
      real(dp) ::  alpha, energy_min

      character(len=40) :: fmt1, fmt2

      open( unit = 82, file='min')
      open( unit = 80, file='min.l')
      
      max_stepf = 1.0e-2_dp
      min_stepf = 1.0e-8_dp

      alpha = 1.d0
      stepf = max_stepf

      fmt1 = '( 2x, f8.6, 5x, f16.12 )'
      fmt2 = '( /,5x,f8.6,5x,f16.12,/)'


 100  continue

      call get_en_2( alpha - stepf,  energy( 1 ) )
      call get_en_2( alpha        ,  energy( 2 ) )
      call get_en_2( alpha + stepf,  energy( 3 ) )

      write( 80, fmt1 ) alpha - stepf,  energy( 1 )
      write( 80, fmt1 ) alpha       ,  energy( 2 ) 
      write( 80, fmt1 ) alpha + stepf,  energy( 3 )


 200  continue

!  Determine which of the nine points has the lowest energy.

      energy_min = energy( 1 )
      imin = 1
      do i = 1, 3
          if( energy( i )  <  energy_min ) then
              energy_min = energy( i )
              imin = i
          end if
      enddo


!  According to which point has the lowest energy, recalculate
!  only the points which were not calculated before.

      if      ( imin  ==  1 )  then

          alpha = alpha - stepf
          
          write( 82, fmt1 )   alpha, energy_min
          write( 80, fmt2 )   alpha, energy_min
          call flush( 80 )

          energy( 3 ) = energy( 2 )
          energy( 2 ) = energy( 1 )
          call get_en_2( alpha - stepf,  energy( 1 ) )

          write( 80, fmt1 ) alpha-stepf,  energy( 1 )
          write( 80, fmt1 ) alpha     ,  energy( 2 )
          write( 80, fmt1 ) alpha+stepf,  energy( 3 )

          go to 200


      else if ( imin  ==  2 )  then

          alpha = alpha  

          write( 82, fmt1 )   alpha, energy_min
          write( 80, fmt2 )   alpha, energy_min
          call flush( 80 )


          stepf = stepf / 10.d0
          if( stepf  <  min_stepf )  then
              write(  6, fmt1 )   alpha, energy_min
              write( 82, fmt1 )   alpha, energy_min
              write( 80, fmt2 )   alpha, energy_min

              write( 80, '( "EBOND = ", F12.6 )' )    ebond / nd
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


      else if ( imin  ==  3 )  then

          alpha = alpha + stepf
          
          write( 82, fmt1 )   alpha, energy_min
          write( 80, fmt2 )   alpha, energy_min
          call flush( 80 )

          energy( 1 ) = energy( 2 )
          energy( 2 ) = energy( 3 )
          call get_en_2( alpha + stepf,  energy( 3 ) )
          
          write( 80, fmt1 ) alpha-stepf,  energy( 1 )    
          write( 80, fmt1 ) alpha     ,  energy( 2 )    
          write( 80, fmt1 ) alpha+stepf,  energy( 3 )    
          
          go to 200

      end if


      go to 100


      return
      end










