 
!****************************************************************************
!                                                                           *
!                                                                           *
!                               Matous Mrovec                               *
!              Department of Materials Science and Engineering              *
!                        University of Pennsylvania                         *
!                                                                           *
!                                                                           *
!****************************************************************************
!                                                                           *
!                                                                           *
! This subroutine relaxes block for different volumes.                      *
!                                                                           *
!                                                                           *
!****************************************************************************


      subroutine volrelax( )
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
!      include "Include/Hamilt.array"
!      include "Include/Moment.array"
!      include "Include/NebList.array"
      include "Include/PosVel.array"
      include "Include/BEC.array"
!      include "Include/Relax.array"


!
!    Common blocks introduced by A.Girshick.
!    
 
      include "Include/ag.conn"

!
!    Declare the simple variables.
!

      real(dp) :: rxab,ryab,rzab,rab,dalpha,vid
      real(dp) ::  alpha_lo, alpha_hi, alpha_step, alpha
      real(dp) ::  a1a2,a1a3,a2a3
      integer           i,j,k,it,bt
      integer           ia,ib
      integer           points, flag


!
!    Declare the arrays.
!
      real(dp) :: ta(3,3)
      real(dp) :: a1(3,3)
      real(dp) :: e0(30),e1(30)
      real(dp) :: v(10)
      real(dp) :: d1(3,100)

!  Define the mesh for scaling (usually 3 points).
      write(*,'(A)')'VOLRELAX'

      alpha_hi = thi
      alpha_lo = tlo
      points = nt

      write(6,*) tlo,thi,nt

      if (points < 2) then
         alpha_hi = 0.0
         alpha_lo = 0.0
         alpha_step = 0.0
      else
         alpha_step = ( alpha_hi - alpha_lo ) / ( points - 1 )
      endif

      do i = 1,nd,1  
         d1(1,i)=d(1,i)
         d1(2,i)=d(2,i)
         d1(3,i)=d(3,i)
      enddo

      a1(1,1) = a(1,1)
      a1(1,2) = a(1,2)
      a1(1,3) = a(1,3)

      a1(2,1) = a(2,1)
      a1(2,2) = a(2,2)
      a1(2,3) = a(2,3)

      a1(3,1) = a(3,1)
      a1(3,2) = a(3,2)
      a1(3,3) = a(3,3)

      vid = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) & 
     &    + a(1,2)*(a(2,3)*a(3,1)-a(2,1)*a(3,3)) & 
     &    + a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))

      do i = 1,3,1
         do j = 1,3,1
            ta(i,j) = 0.0_dp
            do k = 1,3,1
               ta(i,j) = ta(i,j) + a1(i,k)*tmatrix(k,j)
            enddo
         enddo
      enddo

      dalpha = (thi - tlo)/real(nt-1, dp)
      talpha(1) = tlo

      idebug=0

      do it = 1,nt,1

         eflo = 0.d0
         do i = 1, nd, 1
             de( i ) = 0.d0
         enddo

!  Creating the strained block.

         do i = 1,3,1
            do j = 1,3,1
               a(j,i) = a1(j,i) + talpha(it)*ta(j,i)
            enddo
         enddo

         v(it) = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) & 
     &        + a(1,2)*(a(2,3)*a(3,1)-a(2,1)*a(3,3)) & 
     &        + a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))

         write(9,'(''A = '',3F10.5)') a 
         do i = 1,nd,1
            ad(1,i) = a(1,1)*d1(1,i)+a(2,1)*d1(2,i)+a(3,1)*d1(3,i)
            ad(2,i) = a(1,2)*d1(1,i)+a(2,2)*d1(2,i)+a(3,2)*d1(3,i)
            ad(3,i) = a(1,3)*d1(1,i)+a(2,3)*d1(2,i)+a(3,3)*d1(3,i)
            if (idebug == 1) write(9    ,'(''AD('',I4,'') = '', & 
     &           3G13.5)') i,ad(1,i),ad(2,i),ad(3,i)
         enddo
         do i = 1,ninert,1
            adinert(1,i) = a(1,1)*dinert(1,i)+a(2,1)*dinert(2,i)+ & 
     &                     a(3,1)*dinert(3,i)
            adinert(2,i) = a(1,2)*dinert(1,i)+a(2,2)*dinert(2,i)+ & 
     &                     a(3,2)*dinert(3,i)
            adinert(3,i) = a(1,3)*dinert(1,i)+a(2,3)*dinert(2,i)+ & 
     &                     a(3,3)*dinert(3,i)
         enddo


         flag = 1
         call getetot(flag)
         e0(it) = eprom+ebond+epair-eatom-uent

!     Reset chain characteristics and occupancies.

         call erasab()
         call getocc(1,locc)
!         EATOM = GETEATOM(Z,ND)

         call relax()
         e1(it) = eprom+ebond+epair-eatom-uent

         call erasab()      

         talpha(it+1) = talpha(it) + dalpha

      enddo


      open ( unit = 90, file = "volrel.dat")
      do i = 1, points

          alpha = alpha_lo + ( i - 1 ) * alpha_step
          write( 90, '( 3(F12.8, 4X) )' ) alpha, e0(i)/nd, e1(i)/nd

      enddo
      close( 90 )

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

