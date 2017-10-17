 
      subroutine get_srtun(m,mu,la)
          use mod_precision


          use mod_all_scalar

          use mod_const

          use mod_srt

          use ab_io

!
!    This is a function to evaluate the finite temperature Aoki integrals Un.
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

      real(dp) :: nu,stau,mu
      real(dp) :: taur,nur,nurs,sum
      real(dp) :: e0,e1

      integer la,r,s,n,m

!
!    Evaluate the integrals.
!

      nu = -cmplx((mu-lainf(la))/(2.0_dp*lbinf(la)), kind=dp)
      stau = -cmplx(kt/(2.0_dp*lbinf(la)), kind=dp)

      e0 = nu - fdx0*stau
      e1 = nu + fdx0*stau

      call get_ai_un(e0,mrec,root(1,la),f1(1,la),f2(1,la), & 
     &               forder(la),2*mfdpol+1+m,srtin,srtrn,srtun0)

      call get_ai_un(e1,mrec,root(1,la),f1(1,la),f2(1,la), & 
     &               forder(la),2*mfdpol+1+m,srtin,srtrn,srtun1)

      do n = 0,m,1
         srtun(n) = 0.5_dp*(srtun0(n)+srtun1(n))
         taur = stau
         nur = -nu
         do r = 0,mfdpol,1
            nurs = nur
            sum = 0.0_dp
            do s = 0,2*r+1,1
               sum = sum & 
     &             + triang(2*r+2,s+1)*nurs*(srtun1(s+n)-srtun0(s+n))
               nurs = -nurs/nu
            enddo
            srtun(n) = srtun(n) + fdpol(r)*sum/taur
            taur = taur*stau*stau
            nur = nur*nu*nu
         enddo
      enddo

      end

!----------------------------------------------------------------

      subroutine get_srtvn(m,mu,la)
          use mod_precision


          use mod_all_scalar

          use mod_const

          use mod_srt

          use ab_io

!
!    This is a function to evaluate the finite temperature Aoki integrals Vn.
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

      complex(dp) :: sum

      real(dp) :: nu,stau,mu
      real(dp) :: taur,nur,nurs
      real(dp) :: e0,e1

      integer la,r,s,n,m

!
!    Evaluate the integrals.
!

      nu = -cmplx((mu-lainf(la))/(2.0_dp*lbinf(la)), kind=dp)
      stau = -cmplx(kt/(2.0_dp*lbinf(la)), kind=dp)

      e0 = nu - fdx0*stau
      e1 = nu + fdx0*stau

      call get_ai_vn(2*mfdpol+1+m,e0,srtvn0,mrec,root(1,la), & 
     &               f1(1,la),f2(1,la),forder(la),srtin,srtjn,srtrn)

      call get_ai_vn(2*mfdpol+1+m,e1,srtvn1,mrec,root(1,la), & 
     &               f1(1,la),f2(1,la),forder(la),srtin,srtjn,srtrn)

      do n = 0,m,1
         srtvn(n) = 0.5_dp*(srtvn0(n)+srtvn1(n))
         taur = stau
         nur = -nu
         do r = 0,mfdpol,1
            nurs = nur
            sum = cmplx(0.0_dp, kind=dp)
            do s = 0,2*r+1,1
               sum = sum + cmplx(triang(2*r+2,s+1)*nurs, kind=dp)* & 
     &                           (srtvn1(s+n)-srtvn0(s+n))
               nurs = -nurs/nu
            enddo
            srtvn(n) = srtvn(n) + cmplx(fdpol(r)/taur, kind=dp)*sum
            taur = taur*stau*stau
            nur = nur*nu*nu
         enddo
      enddo

      end

!----------------------------------------------------------------

      subroutine get_dsrtun(m,mu,la)
          use mod_precision


          use mod_all_scalar

          use mod_const

          use mod_srt

          use ab_io

!
!    This is a function to evaluate the derivatives of the finite
!     temperature Aoki integrals Un.
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

      real(dp) :: nu,stau,mu
      real(dp) :: taur,nur,nurs,sum
      real(dp) :: e0,e1

      integer la,r,s,n,m

!
!    Evaluate the integrals.
!

      nu = -cmplx((mu-lainf(la))/(2.0_dp*lbinf(la)), kind=dp)
      stau = -cmplx(kt/(2.0_dp*lbinf(la)), kind=dp)

      e0 = nu - fdx0*stau
      e1 = nu + fdx0*stau

      call get_ai_un(e0,mrec,root(1,la),f1(1,la),f2(1,la), & 
     &               forder(la),2*mfdpol+m,srtin,srtrn,srtun0)

      call get_ai_un(e1,mrec,root(1,la),f1(1,la),f2(1,la), & 
     &               forder(la),2*mfdpol+m,srtin,srtrn,srtun1)

      do n = 0,m,1
         srtun(n) = 0.0_dp
         taur = 1.0_dp
         nur = 1.0_dp
         do r = 0,mfdpol,1
            nurs = nur
            sum = 0.0_dp
            do s = 0,2*r,1
               sum = sum + triang(2*r+1,s+1)*nurs* & 
     &                     (srtun1(s+n)-srtun0(s+n))
               nurs = -nurs/nu
            enddo
            srtun(n) = srtun(n) + fdpol(r)*real(2*r+1, dp)*sum/taur
            taur = taur*stau*stau
            nur = nur*nu*nu
         enddo
      enddo

      end

!----------------------------------------------------------------

      subroutine get_dsrtvn(m,mu,la)
          use mod_precision


          use mod_all_scalar

          use mod_const

          use mod_srt

          use ab_io

!
!    This is a function to evaluate the derivative of the finite
!     temperature Aoki integrals Vn.
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

      complex(dp) :: sum

      real(dp) :: nu,stau,mu
      real(dp) :: taur,nur,nurs
      real(dp) :: e0,e1

      integer la,r,s,n,m

!
!    Evaluate the integrals.
!

      nu = -cmplx((mu-lainf(la))/(2.0_dp*lbinf(la)), kind=dp)
      stau = -cmplx(kt/(2.0_dp*lbinf(la)), kind=dp)

      e0 = nu - fdx0*stau
      e1 = nu + fdx0*stau

      call get_ai_vn(2*mfdpol+m,e0,srtvn0,mrec,root(1,la), & 
     &               f1(1,la),f2(1,la),forder(la),srtin,srtjn,srtrn)

      call get_ai_vn(2*mfdpol+m,e1,srtvn1,mrec,root(1,la), & 
     &               f1(1,la),f2(1,la),forder(la),srtin,srtjn,srtrn)

      do n = 0,m,1
         srtvn(n) = 0.0_dp
         taur = 1.0_dp
         nur = 1.0_dp
         do r = 0,mfdpol,1
            nurs = nur
            sum = cmplx(0.0_dp, kind=dp)
            do s = 0,2*r,1
               sum = sum + cmplx(triang(2*r+1,s+1)*nurs, kind=dp)* & 
     &                           (srtvn1(s+n)-srtvn0(s+n))
               nurs = -nurs/nu
            enddo
            srtvn(n) = srtvn(n) + cmplx(real(2*r+1, dp)*fdpol(r)/taur, kind=dp)*sum
            taur = taur*stau*stau
            nur = nur*nu*nu
         enddo
      enddo

      end

