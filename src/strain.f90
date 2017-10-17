
   subroutine strain()
      
      use mod_precision
      use mod_const
      use mod_all_scalar
      use topologia
      
      implicit none

      include "Include/Atom.array"
      include "Include/PosVel.array"

      integer  i, j, ndis, num, nni06
      character(len=80) :: filename
      character(len=80) :: info

      real(dp) ::  c11, c22, c12, c13, c33, c66, c44, c55, c45

      common/elc/c11, c22, c12, c13, c33, c66, c44, c55, c45

      
      if (.not. quiet) call writexyz('xyz/prestrain.xyz','prestrain')
      
      read ( 21, '(i1)' ) ndis
      if (.not. quiet) write(  9, '("number of dislocations: ",i2)' ) ndis
      if (.not. quiet) write(  *, '("number of dislocations: ",i2)' ) ndis

      read ( 21, '(1x,a80)') info
      if (.not. quiet) write(  9, '(1x,a80)') info
      if (.not. quiet) write(  *, '(1x,a80)') info
      
      read(21,*)c11, c22, c12, c13, c33, c66, c44, c55, c45
      if (.not. quiet) write(9,*)c11, c22, c12, c13, c33, c66, c44, c55, c45
      if (.not. quiet) write(*,*)c11, c22, c12, c13, c33, c66, c44, c55, c45
      
      do i = 1, ndis
         call impose_strain()
      enddo
      
      iter = 0
      
      
      if (.not. quiet) then
         
         call writeddp('ddplot/strained.dat')
         call writexyz('xyz/strained.xyz','strained')
         call writecfg("cfg/strained.cfg")
         call writexbs("bs/strained.bs")
   
      end if
   end subroutine strain


!_________________________________________________________________________
!


   subroutine impose_strain()
! dlp: this routine has plenty of potentially dangerous usage of ints and real(4) in place of real(8)
! I have fixed a few but not all
      use mod_precision
      use mod_all_scalar
      use mod_const
      use topologia
      
      implicit none

      include "Include/PosVel.array"
      include "Include/ag.conn"


!       real(dp) ::  pr(3), pp(3), f(3,3), g(3,3)
!       real(dp) ::  u(3), phicut(3), x, y, shift(3)
      real(dp) ::  f(3,2), g(3,2)
      real(dp) ::  u(3), x, y, shift(3)
      real(dp) ::  xcore, ycore, ang
!       real(dp) ::  etar, etai, arctg, phcut

! c      real(dp) :: ug,z,arctg0

      real(dp) :: unew(3), uold(3)
      real(dp) :: xs, ys, dist2

      integer           i, k, n, nang, it
      character*80      info

      real(dp) :: c11, c22, c12, c13, c33, c66, c44, c55, c45
      real(dp) :: cc11, lam, fi
      real(dp) :: q, t, x0, y0, x1, y1, u1, u2, u3
      real(dp) :: gat, dat
      real(dp) :: rlz(3), ilz(3), rp(3), ip(3)
      real(dp) :: bx0, by0, bz0
      real(dp) :: bx, by, bz

      
      complex(dp) :: akd(3,3), uc(3), lz(3), p(3)
      
      common/elc/c11, c22, c12, c13, c33, c66, c44, c55, c45


      
      
      
      

! THIS SUBROUTINE IMPOSES A DISLOCATION STRAINFIELD ACCORDING TO THE
! PARAMETERS P,F,G.
! <ANG> IS THE ANGLE (DEGREES, LATER CONVERTED TO RADIANS) THAT THE
! SLIT IN THE ARCTAN PART OF THE DISPLACEMENT FIELD MAKES WITH THE
! POS. X-AXIS THIS IS THE ANGLE IN REAL X,Y-SPACE , WHICH IS OF COURSE
! DIFFERENT FROM ETA (OR X+P(N)*Y, P(N) COMPLEX)-SPACE
! <XCORE,YCORE> ARE THE POSITION OF THE CORE OF THE DISLOCATION
! <pi> IS JUST USED AS THE CONSTANT PI

!  All that this subroutine does with respect to other parts of the 
!  program is that it modifies ("slightly") the values of coordinates.


!  Reading data from the channel 21 and echo-print to channel 9

      if (.not. quiet) write (*,*) "READING FROM CHANNEL 21"
      do i = 1, 2
          read ( 21, '(1X,A80)') info
          if (.not. quiet) write(  9, '(1X,A80)') info
          if (.not. quiet) write(  6, '(1X,A80)') info
      enddo
      
      read ( 21, '(1x,a80)') info
      if (.not. quiet) print *, info
      read ( 21, * )  akd
      if (.not. quiet) print *, akd
      
!       akd = transpose(akd)
      
      read ( 21, '(1x,a80)') info
      if (.not. quiet) print *, info
      read ( 21, * )  p
      if (.not. quiet) print *, p
      
      read ( 21, '(1x,a80)') info
      if (.not. quiet) print *, info
      read ( 21, * ) xcore, ycore, ang
      if (.not. quiet) write(  9, * ) xcore, ycore, ang
      if (.not. quiet) write(  6, * ) xcore, ycore, ang

!     Transfer ANG to radians (to be an argument for COS and SIN)

      nang = ang / 360.0_dp
      if(ang >= 360.0_dp) ang = ang - nang * 360.0_dp
      if(ang < 0.0_dp)   ang = ang - (nang - 1.0_dp) * 360.0_dp
      ang = ang * pi / 180.0_dp

!     Read the displacements parameters
! C           Modified by I. Katzarov



      read ( 21, '(1x,a80)') info
      if (.not. quiet) write(  9, '(1x,a80)') info
      if (.not. quiet) write(  *, '(1x,a80)') info

! c         Read Burger's vector

      read(21,*) bx0, by0, bz0
      if (.not. quiet) write(9,*) bx0, by0, bz0
      if (.not. quiet) write(6,*) bx0, by0, bz0

      call strain_eventually(    nd,      ad, xcore, ycore, ang, akd, p )
      call strain_eventually(ninert, adinert, xcore, ycore, ang, akd, p )
      
      
   end subroutine impose_strain







   
   subroutine strain_eventually(nd, ad, xcore, ycore, ang, akd, p )
      use mod_precision
      use mod_const
      use mod_all_scalar, only : quiet
      
      implicit none

      integer, intent(in) :: nd
      real(dp), intent(inout) :: ad(3,nd)
      real(dp), intent(in) :: xcore, ycore, ang
      complex(dp), intent(in) :: akd(3,3), p(3)
      
      real(dp) :: rp(3), ip(3), x, y, xs, ys, gat, dat, rlz(3), ilz(3), u(3), u1, u2, u3
      
      complex(dp) :: lz(3), uc(3)
      
      integer :: i, n, k

      
      
      rp = real(p)
      ip = aimag(p)
         
      do i= 1, nd
         x =  ad( 1, i )  - xcore
         y =  ad( 2, i )  - ycore
         
! C     -----------------------------------------------------------------
! C     The CLEANING part 1 from Genesis code (Baskes?)

         xs =  x * cos(ang) + y * sin(ang)
         ys = -x * sin(ang) + y * cos(ang)  


         do n = 1, 3
            gat = ip(n) * ys
            dat = xs + rp(n) * ys 

!    c         gat=(ip(n)*ys)/1.d+10
!    c         dat=(xs+rp(n)*ys)/1.d+10

!    c         rlz(n)=sqrt((xs+rp(n)*ys)**2.+(ip(n)*ys)**2.)
!    c         rlz(n)=log(sqrt(dat**2.+gat**2.))

            rlz(n) = log(sqrt(dat*dat + gat*gat)) ! shall be the same as 0.5*log(dat*dat + gat*gat) really            
            ilz(n) = atan2(gat, dat)
            if (ilz(n) < 0.0_dp) ilz(n) = ilz(n) + 2*pi
            
!             if (dat == 0.0_dp) dat = dat + 0.1e-06_dp
!             ilz(n) = atan(gat/dat)
!             if(dat < 0.0_dp) ilz(n) = ilz(n)+pi
!             if((dat > 0.0_dp).and.(gat < 0.0_dp) ) ilz(n) = ilz(n)+2*pi


            
!    c         if(ilz(n) < 0.0) ilz(n)=ilz(n)+2.*pi

            lz(n) = cmplx(rlz(n), ilz(n),dp)
   
         enddo  
! 
! 
!          do k = 1,3
!             uc(k) = 0.0_dp
!             
!             do n = 1,3
!             
! ! cc         uc(k)=uc(k)+(-1./2.*pi*sqrt(-1.))*akd(k,n)*
! ! cc     &   clog(xs+p(n)*ys)
! ! 
! ! c         uc(k)=uc(k)+(0.,1.)*sign(1.,ip(n))*
! ! c     &   akd(k,n)*clog((xs+p(n)*ys)/1.d+0)
! 
!                uc(k) = uc(k)+cmplx(0.0_dp,1.0_dp,kind=dp)*akd(n,k)*lz(n)
!             enddo
! !    c         write(*,*)k, uc(k)
!          enddo
!          
         do k = 1,3
            uc(k) = sum(akd(1:3,k)*lz(1:3))
         end do
!          uc = uc * cmplx(0.0_dp,1.0_dp)
         uc = cmplx(-aimag(uc),real(uc))
         
! !          (a+ib)(c+id) = ac -bd + i(ad + bc)
! !           (0+i)(c+id) = -d + ic
         

         u1 = (0.5_dp/pi)*real(uc(1))
         u2 = (0.5_dp/pi)*real(uc(2))
         u3 = (0.5_dp/pi)*real(uc(3))

         u(1) = u1*cos(ang) - u2*sin(ang)
         u(2) = u1*sin(ang) + u2*cos(ang)
         u(3) = u3            

! c  Now transfer the strain coordinates from the Finnis-Sin        
! c  to the BOP scale.

!          u = u * latpar


! c  Correcting the coordinates.

         ad( 1:3, i ) = ad( 1:3, i ) + u


         if (.not. quiet) write(  6, '(3(2x,f14.7))' )  u
         if (.not. quiet) write(  9, '(3(2x,f14.7))' )  u
         

      end do

      

   
   end subroutine strain_eventually










   ! 
   ! 
   !       do n = 1, 3
   !           read ( 21, '(2(2X,E14.7))' )  pr( n ), pp( n )
   !           if (.not. quiet) write(  9, '(2(2X,E14.7))' )  pr( n ), pp( n )
   !       enddo
   !       do k = 1, 3
   !           read ( 21, '(3(2X,E14.7))' )  f(k,1), f(k,2), f(k,3)
   !           if (.not. quiet) write(  9, '(3(2X,E14.7))' )  f(k,1), f(k,2), f(k,3)
   !       enddo
   !       do k = 1, 3
   !           read ( 21, '(3(2X,E14.7))' )  g(k,1), g(k,2), g(k,3)
   !           if (.not. quiet) write(  9, '(3(2X,E14.7))' )  g(k,1), g(k,2), g(k,3)
   !       enddo
   !       if (.not. quiet) write( 9, '(/)' )
   ! 
   ! 
   ! ! CALCULATE FOR THE DIFFERENT VALUES OF PR,PP THE ANGLE OF THE CUT
   ! ! IN THE ETA (OR X+PR*Y,PP*Y)-PLANE THIS IS <PHICUT>
   ! ! ( COS(ANG),SIN(ANG) ) IS A POINT ON THE CUT IN THE X,Y-PLANE; TRANSFORM
   ! ! THIS POINT TO THE ETA-PLANE AND CALCULATE THE ANGLE.
   ! 
   ! !  Fill out the PHICUT array - needed to get the correct ARCTG
   ! 
   !       do n = 1, 3
   !           etar=cos(ang)+pr(n)*sin(ang)
   !           etai=pp(n)*sin(ang)
   !           if (etar == 0.0) etar=etar+0.1e-06_dp
   !           phcut=atan(etai/etar)
   !           if(etar < 0.0) phcut=phcut+pi
   !           if( (etar > 0.0).and.(etai < 0.0) ) &
   !      &        phcut=phcut+2*pi
   !           phicut(n)=phcut
   !       enddo
   ! 
   ! ! NOW, PHICUT IS THE ANGLE OF THE CUT IN THE ETA-PLANE. NOTE THAT PHICUT IS
   ! ! IN [0,2*pi]
   ! 
   ! !
   ! !  Now we are ready to impose the  dislocation strainfield.
   ! !
   ! 
   ! !  First modify the positions of the active atoms.
   ! 
   !       do i= 1, nd
   ! 
   !           x =  ad( 1, i ) / latpar - xcore
   !           y =  ad( 2, i ) / latpar - ycore
   ! 
   !           do k = 1, 3
   !               u(k)=0
   !          
   !               do n = 1, 3
   !                   if((x+pr(n)*y)  ==  0.0)  x=x+0.1e-06
   ! 
   ! !  FIRST DETERMINE THE ARGUMENT OF THE COMPLEX NUMBER X+P*Y, BUT BE
   ! !  CAREFUL TO HAVE ONLY ONE DISCONTINUITY OF 2*pi ALONG THE POSITIVE
   ! !  X-AXIS.  THE ARGUMENT WILL BE IN THE INTERVAL [0,2*pi].
   ! 
   ! !  Get the correct value of ARCTG
   ! 
   !                   arctg=atan( (pp(n)*y)/(x+pr(n)*y) )
   !                   if((x+pr(n)*y) < 0) arctg=arctg+pi
   !                   if( ((x+pr(n)*y) > 0).and.((pp(n)*y) < 0) ) &
   !      &            arctg=arctg+2*pi
   ! 
   ! ! NOW THE SLIT IS ALONG THE POSITIVE X-AXIS IN THE ETA (OR X+PR*Y,PP*Y)
   ! ! -PLANE.
   ! 
   ! ! IF THE SLIT SHOULD BE SOMEWHERE ELSE WE MAKE A CORRECTION,
   ! ! ACCORDING TO PHICUT
   ! 
   !                   if (arctg < phicut(n)) arctg=arctg+2*pi
   ! 
   ! !  Imposing the displacements !
   ! 
   !                   u(k) = u(k)+f(k,n)*log((x+pr(n)*y)**2+(pp(n)*y)**2)+ g(k,n)*arctg
   ! 
   !               enddo
   ! 
   ! !  Normalizing the displacements.
   ! 
   ! ! Could be dangerous not to specify the precision as single is the default
   ! ! dlp
   ! !               u(k) = u(k) * (-1./(2.*pi))
   !               u(k) = u(k) * (-1.0_dp/(2*pi))
   ! 
   ! 
   ! !  Now transfer the strain coordinates from the Finnis-Sinclair scale
   ! !  to the BOP scale.
   ! 
   !               u(k) = u(k) * latpar
   ! 
   ! !  Correcting the coordinates.
   ! 
   !           enddo
   ! 
   !           ad( 1:3, i ) = ad( 1:3, i ) + u
   ! 
   !           if (.not. quiet) write(  9, '(3(2X,E14.7))' )  u
   ! 
   !       enddo
   ! 
   ! 
   ! !  Now do exactly the same thing for inert atoms.
   ! 
   !       do i= 1, ninert
   ! 
   !           x = adinert( 1, i ) / latpar - xcore
   !           y = adinert( 2, i ) / latpar - ycore
   ! 
   !           do k = 1, 3
   !               u(k)=0
   !          
   !               do n = 1, 3
   !                   if((x+pr(n)*y)  ==  0.0_dp)  x=x+0.1e-06_dp
   ! 
   !                   arctg = atan( (pp(n)*y)/(x+pr(n)*y) )
   !                   if((x+pr(n)*y) < 0) arctg = arctg + pi
   !                   if( ((x+pr(n)*y) > 0).and.((pp(n)*y) < 0) ) arctg=arctg+2*pi
   ! 
   !                   if (arctg < phicut(n)) arctg = arctg + 2*pi
   ! 
   !                   u(k)=u(k)+f(k,n)*log((x+pr(n)*y)**2+(pp(n)*y)**2)+ g(k,n)*arctg
   !               enddo
   ! 
   ! !  dlp             
   ! !               u(k) = u(k) * (-1./(2.*pi))
   ! !               u(k) = u(k) * latpar
   !                u(k) = u(k) * (-1.0_dp/(2*pi)) * latpar
   ! 
   !           enddo
   ! 
   !           adinert( 1:3, i ) = adinert( 1:3, i ) + u
   ! 
   !       enddo
   ! 
   ! 
   !       end
