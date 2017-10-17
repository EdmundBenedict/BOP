 
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
!          stresses for the FCC unit cell. The end of the output            *
!             file may be later used in the fitting procedure.              *
!                                                                           *
!                                                                           *
!****************************************************************************


      subroutine elcon_fcc()
          use mod_precision


          use mod_all_scalar

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
!      include "Include/Hamilt.array"
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
      real(dp) :: c11exp, c12exp, c44exp 
      real(dp) :: c11,    c12,    c44
      real(dp) :: kkexp, rrexp
      real(dp) :: kk,    rr
      real(dp) :: sigma11
!
!     Experimental data for Ir taken from Simmons & Wang, ref. 159
!
      c11exp = 3.620_dp
      c12exp = 1.510_dp
      c44exp = 1.598_dp
      kkexp  = 2.214_dp
      rrexp  = 0.000


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


      strain(  2, 1, 1 ) = 1.d0
      strain(  2, 1, 2 ) = 0.d0
      strain(  2, 1, 3 ) = 0.d0
      
      strain(  2, 2, 1 ) = 0.d0
      strain(  2, 2, 2 ) = 1.d0
      strain(  2, 2, 3 ) = 0.d0
 
      strain(  2, 3, 1 ) = 0.d0
      strain(  2, 3, 2 ) = 0.d0
      strain(  2, 3, 3 ) = 0.d0

      fileext(  2 ) = '2C11+2C12'
      
      
      strain(  3, 1, 1 ) = 1.d0
      strain(  3, 1, 2 ) = 0.d0
      strain(  3, 1, 3 ) = 0.d0
     
      strain(  3, 2, 1 ) = 0.d0
      strain(  3, 2, 2 ) = 1.d0
      strain(  3, 2, 3 ) = 0.d0
 
      strain(  3, 3, 1 ) = 0.d0
      strain(  3, 3, 2 ) = 0.d0
      strain(  3, 3, 3 ) = 1.d0

      fileext(  3 ) = '3C11+6C12'
      
      
      strain(  4, 1, 1 ) = 1.d0
      strain(  4, 1, 2 ) = 0.d0
      strain(  4, 1, 3 ) = 0.d0
      
      strain(  4, 2, 1 ) = 0.d0
      strain(  4, 2, 2 ) =-1.d0
      strain(  4, 2, 3 ) = 0.d0
 
      strain(  4, 3, 1 ) = 0.d0
      strain(  4, 3, 2 ) = 0.d0
      strain(  4, 3, 3 ) = 0.d0

      fileext(  4 ) = '2C11-2C12'
      
      
      strain(  5, 1, 1 ) = 0.d0
      strain(  5, 1, 2 ) = 0.d0
      strain(  5, 1, 3 ) = 1.d0
      
      strain(  5, 2, 1 ) = 0.d0
      strain(  5, 2, 2 ) = 0.d0
      strain(  5, 2, 3 ) = 0.d0
 
      strain(  5, 3, 1 ) = 1.d0
      strain(  5, 3, 2 ) = 0.d0
      strain(  5, 3, 3 ) = 0.d0

      fileext(  5 ) = '4C44'
      
      
      strain(  6, 1, 1 ) = 0.d0
      strain(  6, 1, 2 ) = 0.d0
      strain(  6, 1, 3 ) = 1.d0
      
      strain(  6, 2, 1 ) = 0.d0
      strain(  6, 2, 2 ) = 0.d0
      strain(  6, 2, 3 ) = 1.d0
 
      strain(  6, 3, 1 ) = 1.d0
      strain(  6, 3, 2 ) = 1.d0
      strain(  6, 3, 3 ) = 0.d0

      fileext(  6 ) = '8C44'
      
 
      number_of_strains = 6


      open( unit=80, file = 'ec.out' ) 
      do i = 1, number_of_strains

          do j = 1, 3
              do k = 1, 3
                  tmatrix( j, k ) = strain( i, j, k ) 
              enddo
          enddo

          call getelast( i, curvature( i ), slope( i ) )

          write( 80,'(5X, F8.2, 8X, A20)' ) curvature( i ), fileext( i )

      enddo

      c11 = curvature( 1 )
      c12 = 0.5_dp * ( curvature( 2 ) - 2.d0 * curvature( 1 ) )
      c44 = 1.0_dp * ( curvature( 5 ) +  curvature( 6 ) ) / 12.d0
      kk  = curvature( 3 ) / 9.d0
      rr  = curvature( 4 ) / 4.d0

      write( 80,'(//"COMPARE THE SOLUTIONS"/)' )

      write( 80,'(5X, F8.2, 8X, F8.2)' ) curvature( 1 ), c11
      write( 80,'(5X, F8.2, 8X, F8.2)' ) curvature( 2 ),  &
     &                                   2.d0*c11 + 2.d0*c12
      write( 80,'(5X, F8.2, 8X, F8.2)' ) curvature( 3 ),           &
     &                                   3.d0*c11 + 6.d0*c12
      write( 80,'(5X, F8.2, 8X, F8.2)' ) curvature( 4 ),  &
     &                                   2.d0*c11 - 2.d0*c12 
      write( 80,'(5X, F8.2, 8X, F8.2)' ) curvature( 5 ), 4.d0 * c44
      write( 80,'(5X, F8.2, 8X, F8.2)' ) curvature( 6 ), 8.d0 * c44

         
      volume = a(1,1)*a(2,2)*a(3,3)
      write( 80,'(/5X,"VOLUME = ", F10.6/)' )  volume 

      write( 80,'( 62("_")//19X,"in eV/(A**3)",12X,"in 10**11 Pa"// &
     &             17X,"calc.",6X,"exp.",10X,"calc.",6X,"exp."/ &
     &             62("_")/ )' )

      write( 80,'(5X,"C11",8X, F6.3, 5X, F6.3, 8X, F6.3, 5X, F6.3)' ) &
     &        c11/volume, c11exp, c11*1.602/volume, c11exp*1.602

      write( 80,'(5X,"C12",8X, F6.3, 5X, F6.3, 8X, F6.3, 5X, F6.3)' ) &
     &        c12/volume, c12exp, c12*1.602/volume, c12exp*1.602

      write( 80,'(5X,"C44",8X, F6.3, 5X, F6.3, 8X, F6.3, 5X, F6.3)' ) &
     &        c44/volume, c44exp, c44*1.602/volume, c44exp*1.602

      write( 80,'(5X," K ",8X, F6.3, 5X, F6.3, 8X, F6.3, 5X, F6.3)' ) &
     &        kk/volume, kkexp, kk*1.602/volume, kkexp*1.602

      write( 80,'(5X," R ",8X, F6.3, 5X, F6.3, 8X, F6.3, 5X, F6.3)' ) &
     &        rr/volume, rrexp, rr*1.602/volume, rrexp*1.602



!  This part is for writing out an input file for the fitting program.

      sigma11 = slope( 3 ) / volume / 3.d0
      write ( 80,'(//"Sigma11 =",F16.12/)' ) sigma11
      write ( 80,'("dE/dA =",F16.12/)' )  &
     &               slope( 3 ) / a( 1, 1 ) / 4.d0

      write ( 80,'(65("_"))' )
      write ( 80,'("ECOH_BOND"/)' )
      write ( 80,'("DEDA_BOND")' )
      write ( 80,'(F18.12)' )  slope( 3 ) / a( 1, 1 ) / 4.d0
      write ( 80,'("SIG11_BOND")' )
      write ( 80,'(F18.12)' )  sigma11
      write ( 80,'("C11_BOND")' )
      write ( 80,'(F18.12)' )  c11 / volume - sigma11
      write ( 80,'("C12_BOND")' )
      write ( 80,'(F18.12)' )  c12 / volume
      write ( 80,'("C44_BOND")' )
      write ( 80,'(F18.12)' )  c44 / volume - sigma11 / 2.d0
!
!     Added by Marc Cawkwell, 22 September 2002
!
      write ( 80,'("Experimental Cauchy Pressure = -0.044")' )      
      write ( 80,'("Cauchy pressure, 0.5*(c12 - c44)")' )
      write ( 80,'(F18.12)' ) 0.5_dp * ((c12 / volume) -  &
     &                         (c44 / volume - sigma11 / 2.d0))
      write ( 80,'("Env. term must contribute")' )
      write ( 80,'(F18.12)' ) -0.044_dp - 0.5_dp * ((c12 / volume) -  &
     &                        (c44 / volume - sigma11 / 2.d0))
!
      print *, "CAUCHY PRESSURE"
      print *, 0.5_dp * ((c12 / volume) -  & 
     &                  (c44 / volume - sigma11 / 2.d0))
! 
      call tag( 80 )

      close( 80 )


      return
      end

