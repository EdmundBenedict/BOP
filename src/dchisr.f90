 
      subroutine dchisr_a(ef,nmax,la)
          use mod_precision
          use mod_all_scalar
          use mod_const
          use mod_srt
          use ab_io
          use mod_chi


!
!    This is a routine to evaluate the derivatives of the susceptibilities.
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
      call get_dsrtvn(pdord,ef,la)
      call get_dsrtun(0,ef,la)
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
      prod = prod*prod/(2.0_dp*lbinf(la)*pi*kt)

!
!    dChi_nn
!

      do n = 0,nmax,1
         sum = 0.0_dp
         do i = 0,pdord,1
            sum = sum + pd(i,n)*real(srtvn(i))
         enddo
         dchia(n,la) = -prod*sum
      enddo

!
!    dChi_n(n-1)
!

      chib(0,la) = 0.0_dp
      do n = 1,nmax,1
         sum = 0.0_dp
         do i = 0,pd1ord,1
            sum = sum + pd1(i,n)*real(srtvn(i))
         enddo
         if (n <= nrec) then
            dchib(n,la) = prod*((lbinf(la)/brec(n,la))*u0 & 
     &              +       (brec(n,la)/(2.0_dp*lbinf(la)))*sum)
         else
            dchib(n,la) = prod*(u0+0.5_dp*sum)
         endif
      enddo

      end

!-------------------------------------------------------------------------

      subroutine dchisr_c(ef,nmax,la)
          use mod_precision
          use mod_all_scalar
          use mod_const
          use mod_srt
          use ab_io
          use mod_g0n, only : getg0n
          use mod_chi
!
!    This is a routine to evaluate the derivatives of the susceptibilities.
!

      implicit none

      
      real(dp), intent(in) :: ef
      integer, intent(in) :: nmax, la
      
      include "Include/Atom.array"
      include "Include/Misc.array"

      complex(dp) :: e

      real(dp) :: fac,delta,simp
      real(dp) :: x,x0,x1,dx, dlt

      integer :: n
      integer :: p,m
      integer, parameter :: m0=10,mx=1000
!       integer :: nrec

      real(dp) :: ya(0:nmax,mx)
      real(dp) :: yb(1:nmax,mx)


!
!    Evaluate the susceptibilities.
!

      dchia(:,la) = 0.0_dp
      dchib(:,la) = 0.0_dp

      x0 = max((lainf(la)-2*lbinf(la)-ef)/kt,-10.0_dp)
      x1 = min((lainf(la)+2*lbinf(la)-ef)/kt,10.0_dp)
      m = min(mx,max(m0,nint(100*real(nmax, dp)*kt/lbinf(la))))
      m = 2*(m/2)+1

      x = x0
      dx = (x1-x0)/real(m-1, dp)

      do p = 1,m
         e = cmplx(ef + kt*x, kind=dp)
         call getg0n(e, arec(0:nrec,la), brec(0:nrec+1,la), g0n, nmax, nrec, lainf(la), lbinf(la))
         dlt = delta(kt*x,kt)
         ya(0:nmax,p) = aimag(g0n(0:nmax  )*g0n(0:nmax))*dlt
         yb(1:nmax,p) = aimag(g0n(0:nmax-1)*g0n(1:nmax))*dlt
         x = x + dx
      enddo


      fac = (2*kt)/pi
      dchia(0:nmax,la) = fac*simp(dx,ya(0:nmax,1),m)
      dchib(1:nmax,la) = fac*simp(dx,yb(1:nmax,1),m)
      
!       print *, la
!       print *, dchia(0:nmax,la)
!       print *, dchib(1:nmax,la)
!       
      end subroutine dchisr_c

