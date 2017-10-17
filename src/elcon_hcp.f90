 
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
!        This subroutine  calculates all the elastic constants and          *
!          stresses for the HCP unit cell. The end of the output            *
!             file may be later used in the fitting procedure.              *
!                                                                           *
!                                                                           *
!****************************************************************************


      subroutine elcon_hcp()
          use mod_precision


          use mod_const


      implicit none

      integer    strains   
      parameter( strains = 20 )

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
      include "Include/BEC.array"
!      include "Include/BondOrder.array"
!      include "Include/Force.array"
!       include "Include/Hamilt.array"
!      include "Include/Misc.array"
!      include "Include/Moment.array"
!      include "Include/NebList.array"
!      include "Include/NotSRT.array"
      include "Include/PosVel.array"
!      include "Include/RecurCoef.array"
!      include "Include/Relax.array"
!      include "Include/SRT.array"

!
!    Common blocks introduced by A.Girshick.
!

!      include "Include/ag.conn"

!
!    Declare the simple variables.
!

      real(dp) :: strain( strains, 3, 3 )
      character*20  fileext( strains )
      integer  number_of_strains
      integer  i, j, k
      real(dp) :: volume, curvature( strains ), slope( strains )
      real(dp) :: c11exp, c12exp, c13exp, c33exp, c44exp, c66exp 
      real(dp) :: c11,    c12,    c13,    c33,    c44,    c66
      real(dp) :: kkexp, rrexp, hhexp
      real(dp) :: kk,    rr,    hh
      real(dp) :: sigma11,  sigma33
!Dimitar Pashov
        real(8) time_1, time_2

!Dimitar Pashov
        call cpu_time(time_1)
      c11exp = 1.099
      c12exp = 0.542
      c13exp = 0.426
      c33exp = 1.189
      c44exp = 0.317
      c66exp = 0.281
      kkexp  = 0.687
      rrexp  = 0.386
      hhexp  = 0.305


      strain(  1, 1, 1 ) = 1.d0
      strain(  1, 1, 2 ) = 0.d0
      strain(  1, 1, 3 ) = 0.d0

      strain(  1, 2, 1 ) = 0.d0
      strain(  1, 2, 2 ) = 0.d0
      strain(  1, 2, 3 ) = 0.d0

      strain(  1, 3, 1 ) = 0.d0
      strain(  1, 3, 2 ) = 0.d0
      strain(  1, 3, 3 ) = 0.d0

      fileext(  1 ) = 'C11'


      strain(  2, 1, 1 ) = 0.d0
      strain(  2, 1, 2 ) = 0.d0
      strain(  2, 1, 3 ) = 0.d0
      
      strain(  2, 2, 1 ) = 0.d0
      strain(  2, 2, 2 ) = 1.d0
      strain(  2, 2, 3 ) = 0.d0
 
      strain(  2, 3, 1 ) = 0.d0
      strain(  2, 3, 2 ) = 0.d0
      strain(  2, 3, 3 ) = 0.d0

      fileext(  2 ) = 'C11'
      
      
      strain(  3, 1, 1 ) = 0.d0
      strain(  3, 1, 2 ) = 0.d0
      strain(  3, 1, 3 ) = 0.d0
     
      strain(  3, 2, 1 ) = 0.d0
      strain(  3, 2, 2 ) = 0.d0
      strain(  3, 2, 3 ) = 0.d0
 
      strain(  3, 3, 1 ) = 0.d0
      strain(  3, 3, 2 ) = 0.d0
      strain(  3, 3, 3 ) = 1.d0

      fileext(  3 ) = 'C33'
      
      
      strain(  4, 1, 1 ) = 1.d0
      strain(  4, 1, 2 ) = 0.d0
      strain(  4, 1, 3 ) = 0.d0
      
      strain(  4, 2, 1 ) = 0.d0
      strain(  4, 2, 2 ) = 1.d0
      strain(  4, 2, 3 ) = 0.d0
 
      strain(  4, 3, 1 ) = 0.d0
      strain(  4, 3, 2 ) = 0.d0
      strain(  4, 3, 3 ) = 0.d0

      fileext(  4 ) = '2C11+2C12'
      
      
      strain(  5, 1, 1 ) = 0.5_dp
      strain(  5, 1, 2 ) = 0.d0
      strain(  5, 1, 3 ) = 0.d0
      
      strain(  5, 2, 1 ) = 0.d0
      strain(  5, 2, 2 ) = 1.d0
      strain(  5, 2, 3 ) = 0.d0
 
      strain(  5, 3, 1 ) = 0.d0
      strain(  5, 3, 2 ) = 0.d0
      strain(  5, 3, 3 ) = 0.d0

      fileext(  5 ) = '1.25C11+C12'
      
      
      strain(  6, 1, 1 ) = 1.d0
      strain(  6, 1, 2 ) = 0.d0
      strain(  6, 1, 3 ) = 0.d0
      
      strain(  6, 2, 1 ) = 0.d0
      strain(  6, 2, 2 ) = 0.d0
      strain(  6, 2, 3 ) = 0.d0
 
      strain(  6, 3, 1 ) = 0.d0
      strain(  6, 3, 2 ) = 0.d0
      strain(  6, 3, 3 ) = 1.d0

      fileext(  6 ) = 'C11+C33+2C13'
      
      
      strain(  7, 1, 1 ) = 0.d0
      strain(  7, 1, 2 ) = 0.d0
      strain(  7, 1, 3 ) = 0.d0
      
      strain(  7, 2, 1 ) = 0.d0
      strain(  7, 2, 2 ) = 1.d0
      strain(  7, 2, 3 ) = 0.d0
 
      strain(  7, 3, 1 ) = 0.d0
      strain(  7, 3, 2 ) = 0.d0
      strain(  7, 3, 3 ) = 1.d0

      fileext(  7 ) = 'C11+C33+2C13'
      
      
      strain(  8, 1, 1 ) = 1.d0
      strain(  8, 1, 2 ) = 0.d0
      strain(  8, 1, 3 ) = 0.d0
      
      strain(  8, 2, 1 ) = 0.d0
      strain(  8, 2, 2 ) = 1.d0
      strain(  8, 2, 3 ) = 0.d0
 
      strain(  8, 3, 1 ) = 0.d0
      strain(  8, 3, 2 ) = 0.d0
      strain(  8, 3, 3 ) = 1.d0

      fileext(  8 ) = '2C11+2C12+C33+4C13'
      
      
      strain(  9, 1, 1 ) = 0.d0
      strain(  9, 1, 2 ) = 0.d0
      strain(  9, 1, 3 ) = 1.d0
      
      strain(  9, 2, 1 ) = 0.d0
      strain(  9, 2, 2 ) = 0.d0
      strain(  9, 2, 3 ) = 0.d0
 
      strain(  9, 3, 1 ) = 1.d0
      strain(  9, 3, 2 ) = 0.d0
      strain(  9, 3, 3 ) = 0.d0

      fileext(  9 ) = '4C44'
      
      
      strain( 10, 1, 1 ) = 0.d0
      strain( 10, 1, 2 ) = 0.d0
      strain( 10, 1, 3 ) = 0.d0
      
      strain( 10, 2, 1 ) = 0.d0
      strain( 10, 2, 2 ) = 0.d0
      strain( 10, 2, 3 ) = 1.d0
 
      strain( 10, 3, 1 ) = 0.d0
      strain( 10, 3, 2 ) = 1.d0
      strain( 10, 3, 3 ) = 0.d0

      fileext( 10 ) = '4C44'
      
      
      strain( 11, 1, 1 ) = 0.d0
      strain( 11, 1, 2 ) = 0.d0
      strain( 11, 1, 3 ) = 1.d0
      
      strain( 11, 2, 1 ) = 0.d0
      strain( 11, 2, 2 ) = 0.d0
      strain( 11, 2, 3 ) = 1.d0
 
      strain( 11, 3, 1 ) = 1.d0
      strain( 11, 3, 2 ) = 1.d0
      strain( 11, 3, 3 ) = 0.d0

      fileext( 11 ) = '8C44'
      
      
      strain( 12, 1, 1 ) = 1.d0
      strain( 12, 1, 2 ) = 0.d0
      strain( 12, 1, 3 ) = 0.d0
      
      strain( 12, 2, 1 ) = 0.d0
      strain( 12, 2, 2 ) = -1.d0
      strain( 12, 2, 3 ) = 0.d0
 
      strain( 12, 3, 1 ) = 0.d0
      strain( 12, 3, 2 ) = 0.d0
      strain( 12, 3, 3 ) = 0.d0

      fileext(  12 ) = '4C66'
      
      
      strain( 13, 1, 1 ) = 0.d0
      strain( 13, 1, 2 ) = 1.d0
      strain( 13, 1, 3 ) = 0.d0
      
      strain( 13, 2, 1 ) = 1.d0
      strain( 13, 2, 2 ) = 0.d0
      strain( 13, 2, 3 ) = 0.d0
 
      strain( 13, 3, 1 ) = 0.d0
      strain( 13, 3, 2 ) = 0.d0
      strain( 13, 3, 3 ) = 0.d0

      fileext( 13 ) = '4C66'
      
      
      strain( 14, 1, 1 ) = -0.5_dp
      strain( 14, 1, 2 ) = 0.d0
      strain( 14, 1, 3 ) = 0.d0
      
      strain( 14, 2, 1 ) = 0.d0
      strain( 14, 2, 2 ) = -0.5_dp
      strain( 14, 2, 3 ) = 0.d0
 
      strain( 14, 3, 1 ) = 0.d0
      strain( 14, 3, 2 ) = 0.d0
      strain( 14, 3, 3 ) = 1.d0

      fileext( 14 ) = '3R'
      
      
      strain( 15, 1, 1 ) = 1.d0
      strain( 15, 1, 2 ) = 0.d0
      strain( 15, 1, 3 ) = 0.d0
      
      strain( 15, 2, 1 ) = 0.d0
      strain( 15, 2, 2 ) = -0.5_dp
      strain( 15, 2, 3 ) = 0.d0
 
      strain( 15, 3, 1 ) = 0.d0
      strain( 15, 3, 2 ) = 0.d0
      strain( 15, 3, 3 ) = -0.5_dp

      fileext( 15 ) = '3H'
      
 
      number_of_strains = 15


      open( unit = 80, file = 'ec.out' ) 

      do i = 1, number_of_strains

          do j = 1, 3
              do k = 1, 3
                  tmatrix( j, k ) = strain( i, j, k ) 
              enddo
          enddo

          call getelast( i, curvature( i ), slope( i ) )

          write( 80,'(5X, F8.2, 8X, A20)' ) curvature( i ), fileext( i )

      enddo

      c11 = 0.5_dp   * ( curvature(1) + curvature(2) )
      c33 = 1.0_dp   *   curvature(3)
      c12 = 1.d0/3.d0  * ( curvature(4) + curvature(5) - 3.25*c11 )
      c13 = 0.25_dp * ( curvature(6) + curvature(7) - 2*c11 - 2*c33 )
      c44 = 0.0625_dp* ( curvature(9) + curvature(10) + curvature(11) )
      c66 = 0.125_dp * ( curvature( 12 ) + curvature( 13 ) )
      kk  = curvature(  8 ) / 9.d0
      rr  = curvature( 14 ) / 3.d0
      hh  = curvature( 15 ) / 3.d0 

!Dimitar Pashov
        call cpu_time(time_2)
        write(6,'("time: ",f12.6," s")') time_2-time_1

      write( 80,'(//"COMPARE THE SOLUTIONS"/)' )

      write( 80,'(5X, F8.2, 8X, F8.2)' ) curvature( 1 ), c11
      write( 80,'(5X, F8.2, 8X, F8.2)' ) curvature( 2 ), c11         
      write( 80,'(5X, F8.2, 8X, F8.2)' ) curvature( 3 ), c33         
      write( 80,'(5X, F8.2, 8X, F8.2)' ) curvature( 4 ),  &
     &                                   2.d0*c11 + 2.d0*c12 
      write( 80,'(5X, F8.2, 8X, F8.2)' ) curvature( 5 ), &
     &                                   1.25_dp*c11 + c12
      write( 80,'(5X, F8.2, 8X, F8.2)' ) curvature( 6 ), &
     &                                   c11 + c33 +2.d0*c13
      write( 80,'(5X, F8.2, 8X, F8.2)' ) curvature( 7 ), &
     &                                   c11 + c33 +2.d0*c13  
      write( 80,'(5X, F8.2, 8X, F8.2)' ) curvature( 8 ), &
     &                                   2*c11+2*c12+c33+4*c13
      write( 80,'(5X, F8.2, 8X, F8.2)' ) curvature( 9 ), 4*c44
      write( 80,'(5X, F8.2, 8X, F8.2)' ) curvature( 10 ), 4*c44
      write( 80,'(5X, F8.2, 8X, F8.2)' ) curvature( 11 ), 8*c44
      write( 80,'(5X, F8.2, 8X, F8.2)' ) curvature( 12 ), 4*c66
      write( 80,'(5X, F8.2, 8X, F8.2)' ) curvature( 13 ), 4*c66
      write( 80,'(5X, F8.2, 8X, F8.2)' ) curvature( 14 ),  &
     &                        0.5*c11 + c33 + 0.5*c12 - 2*c13
      write( 80,'(5X, F8.2, 8X, F8.2)' ) curvature( 15 ),       &
     &                       1.25*c11 + 0.25*c33 - c12 - 0.5*c13

         
      volume = a(1,1)*a(2,2)*a(3,3)
      write( 80,'(/5X,"VOLUME = ", F10.6/)' )  volume 

      write( 80,'( 62("_")//19X,"in eV/(A**3)",12X,"in 10**11 Pa"// &
     &             17X,"calc.",6X,"exp.",10X,"calc.",6X,"exp."/ &
     &             62("_")/ )' )

      write( 80,'(5X,"C11",8X, F6.3, 5X, F6.3, 8X, F6.3, 5X, F6.3)' ) &
     &        c11/volume, c11exp, c11*1.602/volume, c11exp*1.602

      write( 80,'(5X,"C12",8X, F6.3, 5X, F6.3, 8X, F6.3, 5X, F6.3)' ) &
     &        c12/volume, c12exp, c12*1.602/volume, c12exp*1.602

      write( 80,'(5X,"C13",8X, F6.3, 5X, F6.3, 8X, F6.3, 5X, F6.3)' ) &
     &        c13/volume, c13exp, c13*1.602/volume, c13exp*1.602

      write( 80,'(5X,"C33",8X, F6.3, 5X, F6.3, 8X, F6.3, 5X, F6.3)' ) &
     &        c33/volume, c33exp, c33*1.602/volume, c33exp*1.602

      write( 80,'(5X,"C44",8X, F6.3, 5X, F6.3, 8X, F6.3, 5X, F6.3)' ) &
     &        c44/volume, c44exp, c44*1.602/volume, c44exp*1.602

      write( 80,'(5X,"C66",8X, F6.3, 5X, F6.3, 8X, F6.3, 5X, F6.3)' ) &
     &        c66/volume, c66exp, c66*1.602/volume, c66exp*1.602

      write( 80,'(5X," K ",8X, F6.3, 5X, F6.3, 8X, F6.3, 5X, F6.3)' ) &
     &        kk/volume, kkexp, kk*1.602/volume, kkexp*1.602

      write( 80,'(5X," R ",8X, F6.3, 5X, F6.3, 8X, F6.3, 5X, F6.3)' ) &
     &        rr/volume, rrexp, rr*1.602/volume, rrexp*1.602

      write( 80,'(5X," H ",8X, F6.3, 5X, F6.3, 8X, F6.3, 5X, F6.3)' ) &
     &        hh/volume, hhexp, hh*1.602/volume, hhexp*1.602


!  This part is for writing out an input file for the fitting program.

      sigma11 = slope( 4 ) / volume / 2.d0
      sigma33 = slope( 3 ) / volume
      write ( 80,'(//"Sigma11 =",F16.12/)' ) sigma11
      write ( 80,'("Sigma33 =",F16.12/)' ) sigma33 
      write ( 80,'("dE/dA =",F16.12/)' )  &
     &               slope( 4 ) / a( 1, 1 ) / 4.d0
      write ( 80,'("dE/dc =",F16.12/)' )  &
     &               slope( 3 ) / a( 3, 3 ) / 4.d0

      write ( 80,'(65("_"))' )
      write ( 80,'("ECOH_BOND"/)' )
      write ( 80,'("DEDA_BOND")' )
      write ( 80,'(F18.12)' )  slope( 4 ) / a( 1, 1 ) / 4.d0
      write ( 80,'("DEDC_BOND")' )
      write ( 80,'(F18.12)' )  slope( 3 ) / a( 3, 3 ) / 4.d0
      write ( 80,'("SIG11_BOND")' )
      write ( 80,'(F18.12)' )  sigma11
      write ( 80,'("SIG33_BOND")' )
      write ( 80,'(F18.12)' )  sigma33                      
      write ( 80,'("C11_BOND")' )
      write ( 80,'(F18.12)' )  c11 / volume - sigma11
      write ( 80,'("C12_BOND")' )
      write ( 80,'(F18.12)' )  c12 / volume
      write ( 80,'("C13_BOND")' )
      write ( 80,'(F18.12)' )  c13 / volume
      write ( 80,'("C33_BOND")' )
      write ( 80,'(F18.12)' )  c33 / volume - sigma33
      write ( 80,'("C44_BOND")' )
      write ( 80,'(F18.12)' )  c44 / volume - ( sigma11+sigma33 )/4.d0
      write ( 80,'("C66_BOND")' )
      write ( 80,'(F18.12)' )  c66 / volume - sigma11 / 2.d0

      write ( 80,'(///16X,"CAUCHY  PRESSURES",8X,"CALC.",8X,"EXP.")' )
      write ( 80,'(/20X,"C13 - C44",10X,F8.4,7X,F6.3)' ) &
     &   (c13 / volume) - (c44 / volume - ( sigma11+sigma33 )/4.d0), &
     &    c13exp - c44exp
      write ( 80,'(/20X,"C12 - C66",10X,F8.4,7X,F6.3/)' ) &
     &   (c12 / volume) - (c66 / volume - sigma11 / 2.d0), &
     &    c12exp - c66exp
     
      write ( 6,'(///16X,"CAUCHY  PRESSURES",8X,"CALC.",8X,"EXP.")' )
      write ( 6,'(/20X,"C13 - C44",10X,F8.4,7X,F6.3)' ) &
     &   (c13 / volume) - (c44 / volume - ( sigma11+sigma33 )/4.d0), &
     &    c13exp - c44exp
      write ( 6,'(/20X,"C12 - C66",10X,F8.4,7X,F6.3/)' ) &
     &   (c12 / volume) - (c66 / volume - sigma11 / 2.d0), &
     &    c12exp - c66exp

! bases
! +0.5852
! +0.7755

      call tag( 80 )

      close( 80 )


      return
      end
