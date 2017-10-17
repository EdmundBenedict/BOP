 
       function dndmfn(ef,la)
          use mod_precision
          use mod_all_scalar
          use mod_const
          use mod_srt
          use ab_io
          use mod_g0n, only : g00
!
!    This is a routine to evaluate the derivative of the number of
!     electrons with respect to the chemical potential.
!
!   the integration probably has to be synced with the others in ebsos etc...
!

      implicit none
      real(dp) :: dndmfn

      include "Include/Atom.array"

      complex(dp) :: zi

      real(dp) :: ef, delta, wi,ei, prod, x0,x1,x,dx,simp

      integer :: i,la,m
      integer, parameter :: mx = 111

      real(dp) :: y(mx)


!    Evaluate the derivative.


      if (term == 1) then

         if (chi_meth == 1) then
            call get_dsrtun(0,ef,la)
            prod = 0.5_dp
            do i = 1,lchain(la)
               prod = prod*(brec(i,la)/(2*lbinf(la)))
            enddo
            prod = prod*prod/(pi*kt)
            dndmfn = prod*srtun(0)*wt(la)
         else
            x0 = max((lainf(la)-2*lbinf(la)-ef)/kt,-10.0_dp)
            x1 = min((lainf(la)+2*lbinf(la)-ef)/kt, 10.0_dp)
            m = min(mx,max(11,nrec*nint(50*kt)))
            m = 2*(m/2)+1
            x = x0
            dx = (x1-x0)/real(m-1, dp)
            do i = 1,m
               zi = cmplx(ef + kt*x, kind=dp)
               y(i) = aimag(g00(zi,arec(0:lchain(la),la), &
                                 & brec(0:lchain(la)+1,la),lchain(la), & 
                                 & lainf(la),lbinf(la)))*delta(kt*x,kt)
               x = x + dx
            enddo
            dndmfn = (2*wt(la)*kt/pi)*simp(dx,y,m)
         endif

      elseif ((term == 2).or.(term == 3)) then

         dndmfn = 0.0_dp
         do i = 1,nrec+1
            wi = eigvec(1,i,la)
            ei = diag(i,la)
            dndmfn = dndmfn - 2*wi*wi*delta(ei-ef,kt)*wt(la)
         enddo

      endif

!       print *,'dndmfn:',dndmfn

      end

!-----------------------------------------------------------------------------

       function if1(x)
          use mod_precision


          use mod_const

          use mod_srt

          use ab_io
        use mod_g0n, only : g00
!
!    This is a function to evaluate the integrand for performing the integral
!     for dmu/dN.
!

      implicit none
      real(dp) :: if1


      include "Include/Atom.array"

      complex(dp) :: e
      real(dp) :: x,ef,kt,delta

      integer la

!
!    Declare the common block.
!

      common/ifcom/ef,kt,la

!
!    Evaluate the integrand.
!

      e = cmplx(ef + kt*x, kind=dp)

      if1 = aimag(g00(e,arec(0:lchain(la),la),brec(0:lchain(la)+1,la),lchain(la), & 
     &                lainf(la),lbinf(la)))*delta(kt*x,kt)

      end

!-----------------------------------------------------------------------------

       function if2(x)
          use mod_precision


          use mod_const

          use mod_srt

          use ab_io
          use mod_g0n, only : g00
!
!    This is a function to evaluate the integrand for performing the integral
!     for the entropy.
!

      implicit none
      real(dp) :: if2

!
!    Include constants.
!

!      include "Include/ALL.const"

!
!    Include arrays.
!

      include "Include/Atom.array"
!      include "Include/SRT.array"
!      include "Include/RecurCoef.array"

!
!    Declare the simple variables.
!

      complex(dp) :: e

      real(dp) :: x,ef,kt,entden

      integer la

!
!    Declare the common block.
!

      common/ifcom/ef,kt,la

!
!    Evaluate the integrand.
!

      e = cmplx(ef + kt*x, kind=dp)

      if2 = aimag(g00(e,arec(0:lchain(la),la),brec(0:lchain(la)+1,la),lchain(la), & 
     &                lainf(la),lbinf(la)))*entden(kt*x,kt)

      end

