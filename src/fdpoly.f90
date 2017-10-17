 
      subroutine fdpoly()
          use mod_precision


          use mod_const

          use mod_srt

!
!    This is a routine to evaluate the coefficients of the Fermi-Dirac Polynomial
!     for finite temperature work with the square root terminator.
!

      implicit none

!
!    Declare the simple variables.
!

      real(dp) :: x

!
!    Include constants.
!

!      include "Include/ALL.const"

!
!    Include arrays.
!

!      include "Include/SRT.array"

!
!    Evaluate the parameters.
!

      x = fdx0

      fdpol(0) = -0.25_dp

      fdpol(1) = 1.0_dp/108.0_dp

!     FDPOL(1) = -0.07291666666666666D0*(-180.D0*(0.5D0-0.25D0*X)
!    +                                   -21.D0*X)/(X*X*X) + 
!    +        3.333333333333333D0*(0.05D0/(X*X)
!    +                           + 0.02916666666666666D0*
!    +             (-180.D0*(0.5D0-0.25D0*X)-21.*X)/(X*X*X))

!     FDPOL(2) = -0.05D0/(X*X*X*X)
!    +           - 0.02916666666666666D0*(-180.D0*(0.5D0-0.25D0*X)
!    +                              -21.D0*X)/(X**5)

!     FDPOL(3) = 0.01041666666666666D0*(-180.D0*(0.5D0-0.25D0*X)
!    +                                  -21.D0*X)/(X**7)

      end

