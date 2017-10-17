 
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
!         This subroutine caclulates the energy of the unit cell            *
!              with one independent lattice parameter given                 *
!                     the values of these parameter.                        *
!                                                                           *
!                                                                           *
!****************************************************************************


      subroutine get_en_2( lambda, energy )
          use mod_precision


          use mod_all_scalar

          use mod_const


      implicit none

      real(dp) ::  lambda, energy

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

      include "Include/ag.conn"

!
!    Declare the simple variables.
!

      integer  i, flag
      real(dp) :: a1(3,3)


!  Save the original values of A matrix.

      a1(1,1) = a(1,1)
      a1(1,2) = a(1,2)
      a1(1,3) = a(1,3)

      a1(2,1) = a(2,1)
      a1(2,2) = a(2,2)
      a1(2,3) = a(2,3)

      a1(3,1) = a(3,1)
      a1(3,2) = a(3,2)
      a1(3,3) = a(3,3)


!  Creating the strained block.

      a( 1, 1 ) = a1( 1, 1 ) * lambda
      a( 2, 2 ) = a1( 2, 2 ) * lambda
      a( 3, 3 ) = a1( 3, 3 ) * lambda

      do i = 1,nd,1
         ad(1,i) = a(1,1)*d(1,i)+a(2,1)*d(2,i)+a(3,1)*d(3,i)
         ad(2,i) = a(1,2)*d(1,i)+a(2,2)*d(2,i)+a(3,2)*d(3,i)
         ad(3,i) = a(1,3)*d(1,i)+a(2,3)*d(2,i)+a(3,3)*d(3,i)
      enddo

      do i = 1,ninert,1
         adinert(1,i) = a(1,1)*dinert(1,i)+a(2,1)*dinert(2,i)+ & 
     &                  a(3,1)*dinert(3,i)
         adinert(2,i) = a(1,2)*dinert(1,i)+a(2,2)*dinert(2,i)+ & 
     &                  a(3,2)*dinert(3,i)
         adinert(3,i) = a(1,3)*dinert(1,i)+a(2,3)*dinert(2,i)+ & 
     &                  a(3,3)*dinert(3,i)
         enddo


!  Evaluate the enrgy of the strained configuration.

      eflo = 0.d0
      do i = 1, nd, 1
          de( i ) = 0.d0
      enddo

      flag = 1
      call getetot(flag)
 
      energy = (eprom+ebond+epair-eatom-uent) / nd


!  Reassign the original values of A matrix.

      a(1,1) = a1(1,1)
      a(1,2) = a1(1,2)
      a(1,3) = a1(1,3)

      a(2,1) = a1(2,1)
      a(2,2) = a1(2,2)
      a(2,3) = a1(2,3)

      a(3,1) = a1(3,1)
      a(3,2) = a1(3,2)
      a(3,3) = a1(3,3)


      return
      end

