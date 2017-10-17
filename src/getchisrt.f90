 
      subroutine getchisrt_a(ef,nmax,la)
          use mod_precision
          use mod_all_scalar
          use mod_const
          use mod_srt
          use ab_io
          use mod_chi
!
!    This is a routine to evaluate the susceptibilities.
!

      implicit none

      real(dp), intent(in) :: ef
      integer, intent(in) :: nmax, la
      
      real(dp) :: prod,sum
      real(dp) :: u0

      integer n
      integer i

!
!    Evaluate the susceptibilities.
!

      call evalfn(nmax,la)
      call get_srtvn(pdord,ef,la)
      call get_srtun(0,ef,la)
      u0 = srtun(0)

      prod = 1.0_dp
      if (nmax > nrec) then
         do i = 1,nrec,1
            prod = prod*brec(i,la)/(2.0_dp*lbinf(la))
         enddo
         do i = nrec+1,nmax,1
            prod = prod*0.5_dp
         enddo
      else
         do i = 1,nmax,1
            prod = prod*brec(i,la)/(2.0_dp*lbinf(la))
         enddo
      endif
      prod = prod*prod/(2.0_dp*lbinf(la)*pi)

!
!    Chi_nn
!

      do n = 0,nmax,1
         sum = 0.0_dp
         do i = 0,pdord,1
            sum = sum + pd(i,n)*real(srtvn(i))
         enddo
         chia(n,la) = -prod*sum
      enddo

!
!    Chi_n(n-1)
!

      chib(0,la) = 0.0_dp
      do n = 1,nmax,1
         sum = 0.0_dp
         do i = 0,pd1ord,1
            sum = sum + pd1(i,n)*real(srtvn(i))
         enddo
         if (n <= nrec) then
            chib(n,la) = prod*((lbinf(la)/brec(n,la))*u0 & 
     &              +       (brec(n,la)/(2.0_dp*lbinf(la)))*sum)
         else
            chib(n,la) = prod*(u0+0.5_dp*sum)
         endif
      enddo

      end

!---------------------------------------------------------------

      subroutine getchisrt_c(ef,nmax,la)
          use mod_precision
          use mod_all_scalar
          use mod_const
          use mod_srt
          use ab_io
          use mod_g0n, only : getg0n
          use mod_chi
!     
!    This is a routine to evaluate the susceptibilities
!    AMB : via fast integration in complex plane
!

      implicit none

      real(dp), intent(in) :: ef
      integer, intent(in) :: nmax, la

      include "Include/Atom.array"
!       include "Include/Misc.array"

!
!    Declare the simple variables.
!

      complex(dp) :: zp,zfac,zep, g0n(0:nmax)
      real(dp) :: phase,w0
      integer m,p,n

!
!    Evaluate the susceptibilities.
!

!      M     = NINT(MFAC*0.25D0*(EF-LAINF(LA)+2.0D0*LBINF(LA))/KT)
!      WRITE(6,'(/"M = ",I4)') M
      m     = nint(mfac*0.25_dp*(ef - lainf(la) + 2*lbinf(la))/kt)

!       m     = nint(mfac*0.25_dp*1.8_dp/kt)
      w0    = kt*(2*m)
      phase = pi/real(2*m, dp)
      zp    = cmplx(cos(phase),sin(phase), kind=dp)
      zfac  = zp*zp

      chia(:,la) = 0.0_dp
      chib(:,la) = 0.0_dp

      do p = 0,m-1
!          zep = dcmplx(ef)+w0*(zp-dcmplx(1.0d0))
         zep = ef+w0*(zp-1)
         call getg0n(zep,arec(0:lchain(la),la),brec(0:lchain(la)+1,la),g0n,nmax, & 
     &               lchain(la),lainf(la),lbinf(la))

         chia(0:nmax,la) = chia(0:nmax,la) + real(zp*g0n(0:nmax  )*g0n(0:nmax))
         chib(1:nmax,la) = chib(1:nmax,la) + real(zp*g0n(0:nmax-1)*g0n(1:nmax))

         zp = zfac*zp
      enddo

      chia(0:nmax,la) = -2*kt*chia(0:nmax,la)
      chib(1:nmax,la) = -2*kt*chib(1:nmax,la)
      
!      print "('chia(',i0,',',i0,'):',g0)", n,la
!      print *,  chia(:n,la)

      end

