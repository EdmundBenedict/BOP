 
   subroutine getelast( number, curvature, slope )
      use mod_precision
      use mod_all_scalar
      use mod_const
      use mod_conf, only : rtc
      use mod_io, only : print_c
!
!    This is a routine to evaluate elastic constants.
!

      implicit none

      integer      number
      real(dp) ::  curvature, slope


      include "Include/Atom.array"
      include "Include/BEC.array"
      include "Include/PosVel.array"
      include "Include/Relax.array"

!
!    Common blocks introduced by A.Girshick.
!
      integer minflg
      include "Include/ag.conn"



      real(dp) :: polyt,dx
      real(dp) :: x
      real(dp) :: dalpha

      integer it
      integer i,j,k
      integer flag

      character*80 filename, filename2

!*** Lattice vector stretch.
      real(dp) :: ta(3,3)
      real(dp) :: a1(3,3)
!*** Arrays for fitting polynomials.
      real(dp) :: y(0:2)

!
!    Evaluate the elastic constants.
!

      a1 = a
      ta = matmul(a,tmatrix)
      
      if( number  <=  9 ) then
          write(filename,'("bec.",I1)') number
      else
          write(filename,'("bec.",I2)') number
      endif
      open(unit = 2,file = filename)

      

      if (nt > mxnt) then
         write(6,'(''Too many points on binding energy curve.'')')
         write(6,'(''==> Setting NT = '',I2)') mxnt
         nt = mxnt
      elseif (nt < 2) then
         write(6,'(''Too few points on binding energy curve.'')')
         write(6,'(''==> Setting NT = '',I2)') mxnt
         nt = mxnt
      endif

      dalpha = (thi - tlo)/real(nt-1, dp)
      talpha(1) = tlo

!  Cycle over all the points in the binding energy curve. 

      do it = 1,nt,1

         eflo = 0.d0
         de( 1:nd ) = 0.d0

!  Creating the strained block.

         a = a1 + talpha(it)*ta
         
         if (idebug == 1) write(9,'(''A = '',9G13.5)') a 
         
         ad(:,:nd) = matmul(transpose(a), d(:,:nd))
         
         if (idebug == 1) call print_c(ad,6,'ad:')
         
         adinert(:,:ninert) = matmul(transpose(a), dinert(:,:ninert))

         if ( elplusrel  ==  1) then
!          Mark Cawkwell
            print *, "STARTING RELAX"
            call relax()
         endif
!
!  Evaluate the enrgy of the strained configuration.

         flag = 1
         call getetot(flag)
         
!          write(6,'(6(a,1x,f18.8,/))') 'eprom =', eprom, 'ebond =', ebond, 'epair =', epair, 'eatom =', eatom,&
!           &'uent =', uent, 'tot =',  EPROM+EBOND+EPAIR-EATOM-UENT
         
!          bec(it) = eprom+ebond+epair-eatom-uent
         bec(it) = rtc % tote % bind
         talpha(it+1) = talpha(it) + dalpha

      enddo


!  Interpolating the binding energy curve.

      if (polyord > nt-1) then
         write(6,'(''Too high an order for polynomial.'')')
         write(6,'(''==> Setting POLYORD = '',I2)') nt-1
         polyord = nt-1
      elseif (polyord < 1)then
         write(6,'(''Too few points on binding energy curve.'')')
         write(6,'(''==> Setting POLYORD = '',I2)') nt-1
         polyord = nt-1
      endif

      polyb = 0.0_dp
      polys = 0.0_dp

      do i = 1,nt
         polyt = 1.0_dp
         do j = 0,polyord
            polys(j) = polys(j) + polyt
            polyb(j) = polyb(j) + polyt*bec(i)
            polyt = polyt*talpha(i)
         enddo
         do j = polyord+1,2*polyord
            polys(j) = polys(j) + polyt
            polyt = polyt*talpha(i)
         enddo
      enddo

      do i = 0,polyord
         do j = 0,polyord
            polya(i,j) = polys(i+j)
         enddo
      enddo

      call ludcmp(polya,polyord,pmax+1,indx,polyt)
      call lubksb(polya,polyord,pmax+1,indx,polyb)

      do i = 1,nt
         call rpoly(polyb,polyord,talpha(i),y,2)
         write(2,'(5G23.15)') talpha(i),bec(i)/nd,y(0)/nd
      enddo

      write(2,'(''# A['',I2,''] = '',G14.7)')  (i,polyb(i),i=0,polyord,1)

      x = 0.d0
      
      call rpoly(polyb,polyord,x,y,2)
      
      write(2,'(''# Position of minimum  = '',G12.5)') x
      write(2,'(''# Value at minimum     = '',G12.5)') y(0)
      write(2,'(''# Slope at minimum     = '',G12.5)') y(1)
      write(2,'(''# Curvature at minimum = '',G12.5)') y(2)
      
      
      slope = y(1)
      curvature = y(2)

      a = a1

      end

