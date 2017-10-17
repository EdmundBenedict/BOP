 
      subroutine get_quick_elast( number, curvature, slope )
          use mod_precision


          use mod_all_scalar

          use mod_const

!
!    This is a routine to evaluate elastic constants.
!
!       use mod_precision
      implicit none
      

      integer  :: number
      real(dp) :: curvature, slope
      
      include "Include/Atom.array"
      include "Include/BEC.array"
      include "Include/PosVel.array"
      include "Include/Relax.array"
      include "Include/ag.conn"

      integer :: minflg

      real(dp) :: polyt,dx
      real(dp) :: x
      real(dp) :: dalpha

      integer :: it
      integer :: i,j,k
      integer :: flag

      character(len=80) filename, filename2

!*** Lattice vector stretch.
      real(dp) :: ta(3,3)
      real(dp) :: a1(3,3)


!
!    Evaluate the elastic constants.
!

      a1(:,:) = a(:,:)

      if( number  <=  9 ) then
          write(filename,'("bec.",i1)') number
      else
          write(filename,'("bec.",i2)') number
      endif
      open(unit = 2,file = filename,action='write')


      
      ta = matmul(a1,tmatrix)

      nt = 3

!      if (nt > mxnt) then
!         write(6,'(''too many points on binding energy curve.'')')
!         write(6,'(''==> setting nt = '',i2)') mxnt
!         nt = mxnt
!      elseif (nt < 2) then
!         write(6,'(''too few points on binding energy curve.'')')
!         write(6,'(''==> setting nt = '',i2)') mxnt
!         nt = mxnt
!      endif

      
      dalpha = (thi - tlo)/real(nt-1, dp)
      talpha(1) = tlo

!  cycle over all the points in the binding energy curve. 
      

      do it = 1,nt,1

!  creating the strained block.

         a(:,:) = a1(:,:) + talpha(it)*ta(:,:)
         ad(:,:nd) = matmul(transpose(a(:,:)), d(:,:nd))
         adinert(:,:ninert) = matmul(transpose(a(:,:)), dinert(:,:ninert))


         eflo = 0.0_dp
         de(:nd ) = 0.0_dp
         flag = 1
         call getetot(flag)
         
!          write(6,'(6(a,1x,f18.8,/))') 'eprom =', eprom, 'ebond =', ebond, 'epair =', epair, 'eatom =', eatom,&
!           &'uent =', uent, 'tot =',  eprom+ebond+epair-eatom-uent
         
         bec(it) = eprom+ebond+epair-eatom-uent
         talpha(it+1) = talpha(it) + dalpha

      enddo


      do i = 1,nt,1
         write(2,'(2g23.15)') talpha(i),bec(i)/nd
      enddo
      
      slope = 0.5_dp * (bec(3) - bec(1)) / dalpha
      curvature = (bec(3) - 2*bec(2) + bec(1))/(dalpha*dalpha)
      
      
      
      write(2,'(''# position of minimum  = '',g12.5)') x
      write(2,'(''# value at minimum     = '',g12.5)') bec(2)
      write(2,'(''# slope at minimum     = '',g12.5)') slope
      write(2,'(''# curvature at minimum = '',g12.5)') curvature
      
      a(:,:) = a1(:,:)
      
      end subroutine get_quick_elast

