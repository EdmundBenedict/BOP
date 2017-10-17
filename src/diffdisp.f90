 
      subroutine diffdisp(ad,munrl,mcad,dispor,nd,mfcomp,lena, ftot)
          use mod_precision

!
!     This subroutine calculates the difference between the lattice and
!     elastic Greens function displacments as a function of distance
!     from the origin. It is called from latgfbc.f.
!
!     This is specific for the fcc lattice. Sorry about that...
!
!     M. Cawkwell 24/6/2003
!
      implicit none
!
      integer i, j, nd, mfcomp, stat, nn, shlno, incld
      parameter (nn = 41)
      integer count(nn)
      real(dp) :: ellatd(nn), rfromo(nn), latpar, lena(3), mdscut
      real(dp) :: ad(3,nd), munrl(3,nd), mcad(3,nd), dispor(nd)
      real(dp) :: latdsp(nd), eldisp(nd), ltmpd(3,nd), etmpd(3,nd)
      real(dp) :: lgfdisp(nd), ftot(3,nd)
      character*16 disfile
      parameter (mdscut = 20.0_dp)
!
      call fndrec(31,'LATPAR',stat)
      read(31,*) latpar
!
      do i=1,nn
         ellatd(i) = 0.0_dp
         count(i) = 0
      enddo
!
      do i = 1, nd
         do j = 1,3
            ltmpd(j,i) = ad(j,i) - mcad(j,i)
            etmpd(j,i) = munrl(j,i) - mcad(j,i)
         enddo
!
         if (ltmpd(3,i) > 0.5_dp*lena(3)) ad(3,i)=ad(3,i) & 
     &        -lena(3)
         if (ltmpd(3,i) < -0.5_dp*lena(3)) ad(3,i)=ad(3,i) & 
     &        +lena(3)
         if (etmpd(3,i) > 0.5_dp*lena(3)) munrl(3,i)=munrl(3,i) & 
     &        -lena(3)
         if (etmpd(3,i) < -0.5_dp*lena(3)) munrl(3,i)=munrl(3,i) & 
     &        +lena(3)
!
      enddo
!
      do i = 1,nd
         lgfdisp(i) = sqrt((ad(1,i) - munrl(1,i))**2.0_dp + & 
     &   (ad(2,i) - munrl(2,i))**2.0_dp + (ad(3,i) - munrl(3,i))**2.0_dp)
!         LGFDISP(I) = DSQRT((MCAD(1,I) - AD(1,I))**2.0D0 +
!     +   (MCAD(2,I) - AD(2,I))**2.0D0 + 
!     +   (MCAD(3,I) - AD(3,I))**2.0D0)
      enddo
!
      incld = 0
      do 10 i = 1,nd
!
         if (dispor(i) > mdscut) goto 10
!
         if (abs(dispor(i)) < 1.0e-4_dp) then
            incld = 1
            shlno=41
            goto 20
         endif         
         if (abs(dispor(i)-latpar/sqrt(2.0_dp)) < 1.0e-4_dp) then
            incld = 1
            shlno=1
            goto 20
         endif
         if (abs(dispor(i)-latpar) < 1.0e-4_dp) then
            incld = 1
            shlno=2
            goto 20
         endif
         if (abs(dispor(i)-latpar*sqrt(1.5_dp)) < 1.0e-4_dp) then
            incld = 1
            shlno=3
            goto 20
         endif
         if (abs(dispor(i)-latpar*sqrt(2.0_dp)) < 1.0e-4_dp) then
            incld = 1
            shlno=4
            goto 20
         endif
         if (abs(dispor(i)-latpar*sqrt(2.5_dp)) < 1.0e-4_dp) then
            incld = 1
            shlno=5
            goto 20
         endif
         if (abs(dispor(i)-latpar*sqrt(3.5_dp)) < 1.0e-4_dp) then
            incld = 1
            shlno=6
            goto 20
         endif
         if (abs(dispor(i)-latpar*2.0_dp) < 1.0e-4_dp) then
            incld = 1
            shlno=7
            goto 20
         endif         
         if (abs(dispor(i)-latpar*sqrt(4.5_dp)) < 1.0e-4_dp) then
            incld = 1
            shlno=8
            goto 20
         endif
         if (abs(dispor(i)-latpar*sqrt(5.5_dp)) < 1.0e-4_dp) then
            incld = 1
            shlno=9
            goto 20
         endif
         if (abs(dispor(i)-latpar*sqrt(6.0_dp)) < 1.0e-4_dp) then
            incld = 1
            shlno=10
            goto 20
         endif
         if (abs(dispor(i)-latpar*sqrt(6.5_dp)) < 1.0e-4_dp) then
            incld = 1
            shlno=11
            goto 20
         endif
         if (abs(dispor(i)-latpar*sqrt(7.5_dp)) < 1.0e-4_dp) then
            incld = 1
            shlno=12
            goto 20
         endif
         if (abs(dispor(i)-latpar*sqrt(8.0_dp)) < 1.0e-4_dp) then
            incld = 1
            shlno=13
            goto 20
         endif     
         if (abs(dispor(i)-latpar*sqrt(8.5_dp)) < 1.0e-4_dp) then
            incld = 1
            shlno=14
            goto 20
         endif    
         if (abs(dispor(i)-latpar*3.0_dp) < 1.0e-4_dp) then
            incld = 1
            shlno=15
            goto 20
         endif
         if (abs(dispor(i)-latpar*sqrt(9.5_dp)) < 1.0e-4_dp) then
            incld = 1
            shlno=16
            goto 20
         endif
         if (abs(dispor(i)-latpar*sqrt(10.5_dp)) < 1.0e-4_dp) then
            incld = 1
            shlno=17
            goto 20
         endif
         if (abs(dispor(i)-latpar*sqrt(11.0_dp)) < 1.0e-4_dp) then
            incld = 1
            shlno=18
            goto 20
         endif
         if (abs(dispor(i)-latpar*sqrt(12.0_dp)) < 1.0e-4_dp) then
            incld = 1
            shlno=19
            goto 20
         endif
         if (abs(dispor(i)-latpar*sqrt(12.5_dp)) < 1.0e-4_dp) then
            incld = 1
            shlno=20
            goto 20
         endif
         if (abs(dispor(i)-latpar*sqrt(13.5_dp)) < 1.0e-4_dp) then
            incld = 1
            shlno=21
            goto 20
         endif
         if (abs(dispor(i)-latpar*sqrt(15.5_dp)) < 1.0e-4_dp) then
            incld = 1
            shlno=22
            goto 20
         endif
         if (abs(dispor(i)-latpar*4.0_dp) < 1.0e-4_dp) then
            incld = 1
            shlno=23
            goto 20
         endif
         if (abs(dispor(i)-latpar*sqrt(16.5_dp)) < 1.0e-4_dp) then
            incld = 1
            shlno=24
            goto 20
         endif
         if (abs(dispor(i)-latpar*sqrt(17.0_dp)) < 1.0e-4_dp) then
            incld = 1
            shlno=25
            goto 20
         endif
         if (abs(dispor(i)-latpar*sqrt(17.5_dp)) < 1.0e-4_dp) then
            incld = 1
            shlno=26
            goto 20
         endif
         if (abs(dispor(i)-latpar*sqrt(18.0_dp)) < 1.0e-4_dp) then
            incld = 1
            shlno=27
            goto 20
         endif
         if (abs(dispor(i)-latpar*sqrt(18.5_dp)) < 1.0e-4_dp) then
            incld = 1
            shlno=28
            goto 20
         endif
         if (abs(dispor(i)-latpar*sqrt(19.0_dp)) < 1.0e-4_dp) then
            incld = 1
            shlno=29
            goto 20
         endif
         if (abs(dispor(i)-latpar*sqrt(20.5_dp)) < 1.0e-4_dp) then
            incld = 1
            shlno=30
            goto 20
         endif
         if (abs(dispor(i)-latpar*sqrt(21.5_dp)) < 1.0e-4_dp) then
            incld = 1
            shlno=31
            goto 20
         endif
         if (abs(dispor(i)-latpar*sqrt(22.0_dp)) < 1.0e-4_dp) then
            incld = 1
            shlno=32
            goto 20
         endif
         if (abs(dispor(i)-latpar*sqrt(22.5_dp)) < 1.0e-4_dp) then
            incld = 1
            shlno=33
            goto 20
         endif
         if (abs(dispor(i)-latpar*sqrt(23.5_dp)) < 1.0e-4_dp) then
            incld = 1
            shlno=34
            goto 20
         endif
         if (abs(dispor(i)-latpar*sqrt(24.0_dp)) < 1.0e-4_dp) then
            incld = 1
            shlno=35
            goto 20
         endif
         if (abs(dispor(i)-latpar*sqrt(24.5_dp)) < 1.0e-4_dp) then
            incld = 1
            shlno=36
            goto 20
         endif
         if (abs(dispor(i)-latpar*5.0_dp) < 1.0e-4_dp) then
            incld = 1
            shlno=37
            goto 20
         endif
         if (abs(dispor(i)-latpar*sqrt(25.5_dp)) < 1.0e-4_dp) then
            incld = 1
            shlno=38
            goto 20
         endif
         if (abs(dispor(i)-latpar*sqrt(26.5_dp)) < 1.0e-4_dp) then
            incld = 1
            shlno=39
            goto 20
         endif
         if (abs(dispor(i)-latpar*sqrt(27.0_dp)) < 1.0e-4_dp) then
            incld = 1
            shlno=40
            goto 20
         endif
!
 20      if (incld == 1) then
!            RFROMO(SHLNO) = DISPOR(I)
!            COUNT(SHLNO) = COUNT(SHLNO) + 1
!            ELLATD(SHLNO) = ELLATD(SHLNO) + LGFDISP(I)
!            INCLD = 0
         endif
!
 10   continue
!

      write(disfile,'("ellatdiffs.",I1)') mfcomp
      open (unit = 98, status="UNKNOWN", file=disfile)
!
      do i=1,nd
!         IF (LGFDISP(I).GT.0.0035D0) THEN
!            WRITE(98,*) "ND =", I
!            WRITE(98,*) "FTOT =", FTOT(1,I), FTOT(2,I), FTOT(3,I)
!            WRITE(98,*) "AD =", AD(1,I), AD(2,I), AD(3,I)
!            WRITE(98,*) "MUNRL =", MUNRL(1,I), MUNRL(2,I), MUNRL(3,I)
!            WRITE(98,*) "MCAD =", MCAD(1,I), MCAD(2,I), MCAD(3,I)
!            WRITE(98,*) "DISPOR=", DISPOR(I), "  LGFDISP =", LGFDISP(I)
!            WRITE(98,*) "***********************************"
         write(98,*) dispor(i), lgfdisp(i)
!         ENDIF
      enddo
!
!      WRITE (98,*) RFROMO(41), (ELLATD(41)/DFLOAT(COUNT(41)))
!      DO 30 I = 1, (NN-1)
!         IF (I.EQ.23) GOTO 30
!         WRITE (98,*) RFROMO(I), (ELLATD(I)/DFLOAT(COUNT(I)))
! 30   CONTINUE
!     
      return
!     
      end
