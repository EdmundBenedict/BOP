 
      subroutine scrcut(sfin,sfout)
          use mod_precision


          use mod_const

!
!    This is a subroutine which adds various cut-offs
!    for the screening function.
!    (c) Matous Mrovec, MPI, August 2002
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
!      INCLUDE "Include/Hamilt.array"
!      INCLUDE "Include/Moment.array"
!      INCLUDE "Include/NebList.array"
!      INCLUDE "Include/PosVel.array"

!
!    Declare the simple variables.
!
      integer i

      real(dp) :: ff
      real(dp) :: sfin(14), sfout(14)


! Various cut-offs for the screening function.
! If the SCF_CUT=0 there is no cut-off.

      if ( scf_cut == 0 ) then
         do i=1,14
            sfout(i) = ( 1.0_dp - sfin(i) )
         enddo

! Step cut-off: if SF>1 => SF=1

      elseif ( scf_cut == 1 ) then
         do i=1,14
            if (sfin(i) > 1.0_dp) then
               sfout(i) = 0.0_dp
            else
               sfout(i) = ( 1.0_dp - sfin(i) )
            endif
         enddo

! Addition of Fermi-like function to smooth-out the cut-off
! of bond integrals when the screening function is too big.

      elseif ( scf_cut == 2 ) then
         do i=1,14
            if (sfin(i) > 20.0_dp) then
               ff = 0.0_dp
            elseif (sfin(i) < -20.0_dp) then
               ff = 1.0_dp
               write(6,'(/" *** WARNING ***")')
               write(6,'("Screening function has a large negative value")')
               write(6,'(" SF = ",F7.3/)') sfin(i)
            else
               ff = 1.0_dp/(1.0_dp+exp(10.0_dp*(sfin(i)-1.5_dp))) 
            endif
            sfout(i) = ( 1.0_dp - sfin(i) ) * ff
         enddo

! Hyperbolic tangent cut-off.

      elseif ( scf_cut == 3 ) then
         do i=1,14
            sfout(i) = 1.0_dp
            sfout(i) = 1.0_dp - tanh(sfin(i))
         enddo
      else
         write(6,*) "Wrong cut-off flag for the screening function !"
         call panic()
      endif

      end
