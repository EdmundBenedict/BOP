 
      subroutine evalfn(nmax,la)
          use mod_precision


          use mod_all_scalar

          use mod_srt

          use ab_io

!
!    This is a routine to evaluate the coefficients used for
!     susceptibility calculation.
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

!      include "Include/RecurCoef.array"
!      include "Include/SRT.array"

!
!    Declare the simple variables.
!

      integer i,n,nmax,la

!
!    Evaluate the normalised coefficients.
!

      salpha(0) = -(arec(0,la)-lainf(la))/(2.0_dp*lbinf(la))
      sbeta(0) = 0.0_dp

      do i = 1,nrec,1
         salpha(i) = -(arec(i,la)-lainf(la))/(2.0_dp*lbinf(la))
         sbeta(i) = brec(i,la)/(2.0_dp*lbinf(la))
      enddo

      do i = nrec+1,nmax+1,1
         salpha(i) = 0.0_dp
         sbeta(i) = 0.5_dp
      enddo

!     WRITE(9,'(''ALPHA('',I2,'') = '',G12.5,'' BETA('',I2,'') = '',
!    +          G12.5)') (I,SALPHA(I),I,SBETA(I),I=0,NMAX+1,1)

!
!    Find the coefficients for the polynomial w_lm.
!

      do n = 0,nmax,1
         do i = -1,nmax,1
            wn(i,n) = 0.0_dp
         enddo
      enddo
      wn(0,nmax-1) = 1.0_dp

      do n = nmax-2,0,-1
!        WRITE(9,'(''N = '',I3)') N
         do i = 0,nmax-n-1,1
            wn(i,n) = wn(i-1,n+1) - salpha(n+1)*wn(i,n+1) & 
     &              - (sbeta(n+2)**2)*wn(i,n+2)
!           WRITE(9,'(''WN('',I2,'') = '',G22.15)') I,WN(I,N)
         enddo
      enddo

      do n = 0,nmax+1,1
         do i = -1,nmax+1,1
            wn1(i,n) = 0.0_dp
         enddo
      enddo
      wn1(0,nmax) = 1.0_dp

      do n = nmax-1,0,-1
!        WRITE(9,'(''N = '',I3)') N
         do i = 0,nmax-n,1
            wn1(i,n) = wn1(i-1,n+1) - salpha(n+1)*wn1(i,n+1) & 
     &              - (sbeta(n+2)**2)*wn1(i,n+2)
!           WRITE(9,'(''WN1('',I2,'') = '',G22.15)') I,WN1(I,N)
         enddo
      enddo

!
!    Evaluate the polynomial coefficients for P_n.
!

      do n = -1,nmax+2,1
         do i = -1,nmax+2,1
            sp(i,n) = 0.0_dp
         enddo
         sp(n,n) = 1.0_dp
      enddo
      sp(-1,-1) = 0.0_dp
      sp(0,1) = -salpha(0)

      do n = 2,nmax+2,1
         do i = 0,n,1
            sp(i,n) = sp(i-1,n-1)-salpha(n-1)*sp(i,n-1)- & 
     &               (sbeta(n-1)**2)*sp(i,n-2)
         enddo
      enddo

!
!    Find the coefficients for the polynomials D_nN.
!

      spa(0) = sp(0,nmax+1)
      do i = 1,nmax+1,1
         spa(i) = sp(i,nmax+1)-0.5_dp*sp(i-1,nmax)
!        WRITE(9,'(''PA('',I2,'') = '',G22.15)') I,SPA(I)
      enddo

      spb(0) = -0.25_dp*sp(0,nmax)
      do i = 1,nmax+2,1
         spb(i) = 0.5_dp*sp(i-1,nmax+1) - 0.25_dp*sp(i,nmax)
!        WRITE(9,'(''PB('',I2,'') = '',G22.15)') I,SPB(I)
      enddo

      do n = 0,nmax,1
         do i = 0,2*nmax+1,1
            sd(i,n) = 0.0_dp
         enddo
         call prodcoeff(wn1(0,n),nmax-n,spa,nmax+1,sd(0,n),1.0_dp)
         call prodcoeff(wn(0,n),nmax-n-1,spb,nmax+2,sd(0,n),-1.0_dp)
      enddo

!
!    Find the coefficients for the polynomials P_nD_nN and P_(n-1)D_nN.
!

      pdord = 2*nmax+1
      pd1ord = 2*nmax
      do n = 0,nmax,1
         if (idebug == 1) write(9,'(''N = '',I3)') n
         do i = 0,2*nmax+1,1
            pd(i,n) = 0.0_dp
            pd1(i,n) = 0.0_dp
         enddo
         call prodcoeff(sp(0,n),n,sd(0,n),2*nmax-n+1,pd(0,n),1.0_dp)
         if (n > 0) call prodcoeff(sp(0,n-1),n-1,sd(0,n),2*nmax-n+1, & 
     &                              pd1(0,n),1.0_dp)

         if (idebug == 1) then
            do i = 0,2*nmax-n+1,1
               write(9,'(''D('',I2,'') = '',G22.15)') i,sd(i,n)
            enddo
            do i = 0,pdord,1
               write(9,'(''PD('',I2,'') = '',G22.15)') i,pd(i,n)
            enddo
            do i = 0,pd1ord,1
               write(9,'(''PD1('',I2,'') = '',G22.15)') i,pd1(i,n)
            enddo
        endif

      enddo

      end

