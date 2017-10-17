 
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
!    This subroutine calculates energies at close atomic separations and    *
!           compares it with the Rose universal equation of state.          *
!                                                                           *
!                                                                           *
!****************************************************************************


      subroutine rose_curve( )
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
 
!      include "Include/Atom.array"
      include "Include/BEC.array"
      include "Include/PosVel.array"

!
!    Common blocks introduced by A.Girshick.
!    
 
      include "Include/ag.conn"

!
!    Declare the simple variables.
!

      real(dp) ::  eout
      real(dp) ::  alpha_lo, alpha_hi, alpha_step, alpha
      real(dp) ::  ecoh, scaling, ws_rad, bulk_mod
      real(dp) ::  a_rose, e_rose
      real(dp) ::  curvature, slope
      real(dp) ::  dum1, dum2, frombond
      integer           i, flag
      integer           points


!  Define the mesh for the Rose curve.

      alpha_lo = -0.10_dp
      alpha_hi =  0.00_dp
      points = 25
      alpha_step = ( alpha_hi - alpha_lo ) / ( points - 1 )


!  Calculate the parameters for the Rose curve
!
!     DATA FOR IRIDIUM - BULK MODULUS IN UNITS OF eV / A^3
!
      ecoh = -5.71333333333_dp
! Ir     ECOH = -6.93D0
!      ECOH = -4.85D0
! Ti      ECOH = -4.85D0
! Mo      ECOH = -6.82D0
! Nb      ECOH = -7.57D0
! Ta      ECOH = -8.10D0
! W       ECOH = -8.90D0

      write( 6, '(/"Cohesive energy: ",F16.12)' ) ecoh


      ws_rad = exp(log((a(1,1)*a(2,2)*a(3,3)*3)/(4.d0*pi*nd))/3.d0 )
      write( 6, '(/"Wigner-Seitz radius: ",F16.12)' ) ws_rad
!
      bulk_mod = 1.327_dp
! Ir    BULK_MOD = 2.2140D0
!
!      BULK_MOD = 0.719555D0
! Ti      BULK_MOD = 0.719555D0
! Mo      BULK_MOD = 1.639D0
! Nb      BULK_MOD = 1.067D0
! Ta      BULK_MOD = 1.224D0
! W       BULK_MOD = 1.9376D0

      scaling = sqrt( -ecoh/(12.d0*pi*ws_rad*bulk_mod) )
      write( 6, '(/"Distance scaling: ",F16.12)' ) scaling


      tlo = alpha_lo
      thi = alpha_hi
      nt = points
      polyord = 3

      tmatrix( 1, 1 ) = 1.d0
      tmatrix( 1, 2 ) = 0.d0
      tmatrix( 1, 3 ) = 0.d0

      tmatrix( 2, 1 ) = 0.d0
      tmatrix( 2, 2 ) = 1.d0
      tmatrix( 2, 3 ) = 0.d0

      tmatrix( 3, 1 ) = 0.d0
      tmatrix( 3, 2 ) = 0.d0
      tmatrix( 3, 3 ) = 1.d0

      call getelast( 1, curvature, slope )


!      OPEN ( UNIT = 90, FILE = "bond")
!      DO I = 1, POINTS
!
!          READ ( 90, * ) DUM1, FROMBOND, DUM2
!          BEC(I) = BEC(I)  +  FROMBOND * ND
!
!      ENDDO
!      CLOSE( 90 )


      open ( unit = 90, file = "rose")
      do i = 1, points

          alpha = alpha_lo + ( i - 1 ) * alpha_step
          a_rose = alpha * ws_rad / scaling
          e_rose = ecoh * exp( - a_rose ) * &
     &              (1.d0 + a_rose + 0.05_dp*( a_rose**3 ) )
          write( 90, '( 3(F12.8, 4X) )' ) alpha, e_rose, bec(i)/nd      

      enddo
      close( 90 )


      return
      end

