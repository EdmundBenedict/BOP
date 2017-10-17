 
      subroutine getroots(a,b,nrec,nmax,ainf,binf, &
     &                    root,f1,f2,forder, & 
     &                    alpha,beta,p,f,fnag,work, & 
     &                    mrec)

!
!    This is a routine to find the roots of the polynomial f_N(E), and
!     the first and second derivatives of f_N at these roots. This is
!     used for the analytic integration scheme of Masato Aoki.
!
          use mod_precision


      implicit none

!
!    Declare the simple variables.
!

!     DOUBLE COMPLEX PARFRAC,Y

!      DOUBLE PRECISION X,ERRMAX
!      DOUBLE COMPLEX ROOMAX

      real(dp) :: ainf,binf

      integer nrec,mrec
      integer i,n,nmax
      integer forder
      integer ifail

!
!    Declare the arrays.
!

!     DOUBLE COMPLEX DF1(0:2)
!     DOUBLE COMPLEX DF2(0:2)
      complex(dp) :: f1(2*(mrec+3))
      complex(dp) :: f2(2*(mrec+3))
      complex(dp) :: root(2*(mrec+3))
      complex(dp) :: df(0:2)

      real(dp) :: a(0:mrec)
      real(dp) :: b(0:mrec+1)
      real(dp) :: alpha(0:mrec+3)
      real(dp) :: beta(0:mrec+3)
      real(dp) :: p(-1:mrec+4,-1:mrec+4)
      real(dp) :: f(0:2*(mrec+3))
      real(dp) :: fnag(0:2*(mrec+3))
      real(dp) :: work(4*(mrec+3))

!
!    Evaluate the normalised coefficients.
!

      alpha(0) = -(a(0)-ainf)/(2.0_dp*binf)
      beta(0) = 0.0_dp

      do i = 1,nrec,1
         alpha(i) = -(a(i)-ainf)/(2.0_dp*binf)
         beta(i) = b(i)/(2.0_dp*binf)
      enddo

      do i = nrec+1,nmax+1,1
         alpha(i) = 0.0_dp
         beta(i) = 0.5_dp
      enddo

!     WRITE(9,'(''ALPHA('',I2,'') = '',G12.5,'' BETA('',I2,'') = '',
!    +          G12.5)') (I,ALPHA(I),I,BETA(I),I=0,NMAX+1,1)

!
!    Evaluate the polynomial coefficients for P_n.
!

      do n = -1,nmax+2,1
         do i = -1,nmax+2,1
            p(i,n) = 0.0_dp
         enddo
         p(n,n) = 1.0_dp
      enddo
      p(-1,-1) = 0.0_dp
      p(0,1) = -alpha(0)

      do n = 2,nmax+2,1
         do i = 0,n,1
            p(i,n) = p(i-1,n-1)-alpha(n-1)*p(i,n-1)- & 
     &               (beta(n-1)**2)*p(i,n-2)
         enddo
      enddo

!     DO N = 0,NMAX+2,1
!        WRITE(9,'(''N = '',I2)') N
!        WRITE(9,'(''P('',I2,'') = '',G12.5)')
!    +        (I,P(I,N),I=0,N,1)
!     ENDDO

!
!    Evaluate the polynomial coefficients for f_N.
!

      do i = 0,2*nmax+2,1
         f(i) = 0.0_dp
      enddo

      call prodcoeff(p(0,nmax+1),nmax+1,p(0,nmax+1),nmax+1,f,1.0_dp)
      call prodcoeff(p(0,nmax+2),nmax+2,p(0,nmax),nmax,f,-1.0_dp)

!
!    Find the order of the polynomial.
!

      forder = 2*nmax+2
      do while ((abs(f(forder)) <= 1.0e-10_dp).and.(forder > 0))
         forder = forder - 1
      enddo
!     WRITE(9,'(''Order of the polynomial FN = '',I2)')FORDER 
!     DO I = 0,2*NMAX+2,1
!        WRITE(9,'(''F('',I2,'') = '',G22.15)') I,F(I)
!     ENDDO


!
!    Find the roots of the polynomial.
!
!     NOTE: NAG uses the following stupid convention:
!              f(x) = a0 x^n + a1 x^(n-1) + ... + an
!           whereas I use the following completely logical convention:
!              f(x) = a0 + a1 x + ..... + an x^n
!

      if (forder == 0) then

         f1(1) = cmplx(f(0), kind=dp)

      else

!        DO I = 0,FORDER,1
!           FNAG(I) = F(FORDER-I)
!        ENDDO

!  If you uncomment the call of LNAG routine C02AGF, make sure to
!  uncomment three lines above!

!        CALL C02AGF(FNAG,FORDER,.TRUE.,ROOT,WORK,IFAIL)
         call findroots(f,forder,root,ifail)

         if (ifail /= 0) then
            write(6,'(''Unable to find the roots of the polynomial.'')')
            call panic()
         endif

!         DO I = 1,FORDER,1
!            CALL POLY(F,FORDER,ROOT(I),X,0)
!            IF (ABS(DREAL(X)).GT.1.0D-4) THEN
!               WRITE(9,'(''F('',I3,'') = '',2G23.15)') I,X
!               WRITE(6,'(''F('',I3,'') = '',2G23.15)') I,X
!            ENDIF
!         ENDDO

!          ERRMAX = -10.0D0
!          DO I = 1,FORDER,1
!             CALL POLY(F,FORDER,ROOT(I),X,0)
!             IF (ABS(DREAL(X)).GT.ERRMAX) THEN
!                ERRMAX = ABS(DREAL(X))
!                ROOMAX = ROOT(I)
!             ENDIF
!          ENDDO
!          PRINT *,'MAX ERR = ',ERRMAX,'ROO = ',ROOMAX

!
!       Evaluate the first and second derivatives of f_N at the roots.
!

         do i = 1,forder,1

            call poly(f,forder,root(i),df,2)
            f1(i) = df(1)
            f2(i) = df(2)

!           CALL POLY(F,FORDER,ROOT(I)+DCMPLX(1.0D-6),DF1,1)
!           CALL POLY(F,FORDER,ROOT(I)-DCMPLX(1.0D-6),DF2,1)
!           WRITE(9,'(''DF(Analytic) = '',2G23.15)') DF(1)
!           WRITE(9,'(''DF(Numeric ) = '',2G23.15)')
!    +               (DF1(0)-DF2(0))/DCMPLX(2.0D-6)
!           WRITE(9,'(''D2F(Analytic) = '',2G23.15)') DF(2)
!           WRITE(9,'(''D2F(Numeric ) = '',2G23.15)')
!    +               (DF1(1)-DF2(1))/DCMPLX(2.0D-6)
!           WRITE(9,'(''ROOT('',I2,'') = ('',G12.5,'','',G12.5,'')'')')
!    +           I,ROOT(I)
!           WRITE(9,'(''F1('',I2,'') = ('',G12.5,'','',G12.5,'')'')')
!    +           I,F1(I)
!           WRITE(9,'(''F2('',I2,'') = ('',G22.15,'','',G22.15,'')'')')
!    +           I,F2(I)
!           WRITE(9,'(''F2('',I2,'')/F1('',I2,'') = '',2G13.5)')
!    +           I,I,F2(I)/F1(I)

         enddo

      endif

!***********************************************************************
!    Test results.
!
!     DO X = -1.0D0,1.0D0,0.1D0
!        CALL POLY(F,FORDER,DCMPLX(X),Y,0)
!        Y = Y*PARFRAC(F1,ROOT,FORDER,DCMPLX(X))
!        WRITE(9,'(''X = '',G22.15,'' Y = '',2G23.15)') X,Y
!     ENDDO
!
!***********************************************************************

      end

