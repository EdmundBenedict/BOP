 
      subroutine get_ai_un(ef,mrec,z,f1,f2,forder,m,ai_in,ai_rn,ai_un)

!
!    This is a function to evaluate the Aoki integrals U0 through Un.
!
          use mod_precision


      implicit none

!
!    Declare the simple variables.
!

      complex(dp) :: ai_i0

      real(dp) :: ef,phif,pi

      integer forder
      integer mrec
      integer i,n,m

!
!    Define parameters.
!

      parameter (pi = 3.1415926535_dp)

!
!    Declare the arrays.
!

      complex(dp) :: z(2*(mrec+3))
      complex(dp) :: f1(2*(mrec+3))
      complex(dp) :: f2(2*(mrec+3))
      complex(dp) :: ai_in(0:m)

      real(dp) :: ai_rn(0:m)
      real(dp) :: ai_un(0:m)


!
!    Evaluate the number of electrons.
!

      if (ef > 1.0_dp-1.0e-6_dp) then
          phif = 0.0_dp
      elseif (ef < -1.0_dp+1.0e-6_dp) then
          phif = pi
      else
         phif = acos(ef)
      endif

      call get_ai_rn(m,phif,ai_rn)

      if (forder > 0) then
         do n = 0,m,1
            ai_un(n) = 0.0_dp
         enddo
         do i = 1,forder,1
            ai_in(0) = ai_i0(z(i),phif)
            ai_un(0) = ai_un(0) + real(ai_in(0)/f1(i))
            do n = 1,m,1
               ai_in(n) = cmplx(ai_rn(n-1), kind=dp) + z(i)*ai_in(n-1)
               ai_un(n) = ai_un(n) + real(ai_in(n)/f1(i))
            enddo
         enddo
      else
         do n = 0,m,1
            ai_un(n) = ai_rn(n)/f1(1)
         enddo
      endif

      end

!----------------------------------------------------------------

      subroutine get_ai_vn(n,ef,ai_vn,mrec,z,f1,f2,forder,ai_in,ai_jn,ai_rn)

!
!    This is a routine to evaluate the Aoki integrals Vn.
!
          use mod_precision

     
      implicit none

!
!    Declare the simple variables.
!

      real(dp) :: ef,phif,pi

      integer n,mrec
      integer forder,i,j

!
!    Define parameters.
!

      parameter (pi = 3.1415926535_dp)

!
!    Declare the arrays.
!

      complex(dp) :: z(2*(mrec+3))
      complex(dp) :: f1(2*(mrec+3))
      complex(dp) :: f2(2*(mrec+3))
      complex(dp) :: ai_vn(0:n)
      complex(dp) :: ai_in(0:2*(mrec+3))
      complex(dp) :: ai_jn(0:2*(mrec+3))

      real(dp) :: ai_rn(0:2*(mrec+3))

!
!    Evaluate the integral.
!

      if (ef > 1.0_dp-1.0e-6_dp) then
          phif = 0.0_dp
      elseif (ef < -1.0_dp+1.0e-6_dp) then
          phif = pi
      else
         phif = acos(ef)
      endif

      do i = 0,n,1
         ai_vn(i) = cmplx(0.0_dp, kind=dp)
      enddo

      if (forder > 0) then
         do j = 1,forder,1
            call get_ai_jn(n,z(j),phif,ai_jn,ai_in,ai_rn,mrec)
            do i = 0,n,1
               ai_vn(i) = ai_vn(i) & 
     &                  +(ai_jn(i)-ai_in(i)*f2(j)/f1(j))/f1(j)/f1(j)
            enddo
         enddo
      else
         call get_ai_rn(n,phif,ai_rn)
         do i = 0,n,1
            ai_vn(i) = cmplx(ai_rn(i), kind=dp)/(f1(1)**2)
         enddo
      endif

      end

!----------------------------------------------------------------

      subroutine get_ai_in(n,z,phif,ai_in,ai_rn,mrec)
          use mod_precision


!
!    This is a routine to evaluate the Aoki integrals In.
!

      implicit none

!
!    Declare the simple variables.
!

      complex(dp) :: z,ai_i0

      real(dp) :: phif

      integer n,i
      integer mrec

!
!    Declare the arrays.
!

      complex(dp) :: ai_in(0:n)

      real(dp) :: ai_rn(0:n-1)

!
!    Evaluate the function.
!

      if (n == 0) then
         ai_in(0) = ai_i0(z,phif)
      else
         ai_in(0) = ai_i0(z,phif)
         call get_ai_rn(n-1,phif,ai_rn)
         do i = 1,n,1
            ai_in(i) = cmplx(ai_rn(i-1), kind=dp) + z*ai_in(i-1)
         enddo
      endif

      end

!----------------------------------------------------------------

      subroutine get_ai_rn(n,phif,ai_rn)
          use mod_precision


!
!    This is a routine to evaluate the Aoki integrals Rn.
!

      implicit none

!
!    Declare the simple variables.
!

      real(dp) :: phif
      real(dp) :: c,s,s3,s3cn

      integer n,i

!
!    Declare the arrays.
!

      real(dp) :: ai_rn(0:n)

!
!    Evaluate the integral.
!

      s = sin(phif)
      c = cos(phif)
      if (n == 0) then
         ai_rn(0) = 0.5_dp*(s*c-phif)
      elseif (n == 1) then
         ai_rn(0) = 0.5_dp*(s*c-phif)
         ai_rn(1) = -(s*s*s)/3.0_dp
      else
         ai_rn(0) = 0.5_dp*(s*c-phif)
         s3 = s*s*s
         ai_rn(1) = -s3/3.0_dp
         s3cn = s3
         do i = 2,n,1
            s3cn = s3cn*c
            ai_rn(i) = (real(i-1, dp)*ai_rn(i-2)-s3cn)/real(i+2, dp)
         enddo
      endif

      end

!----------------------------------------------------------------

      subroutine get_ai_jn(n,z,phif,ai_jn,ai_in,ai_rn,mrec)
          use mod_precision


!
!    This is a routine to evaluate the Aoki integrals JN.
!

      implicit none

!
!    Declare the simple variables.
!

      complex(dp) :: z,ai_j0

      real(dp) :: phif

      integer n,i,mrec

!
!    Declare the arrays.
!

      complex(dp) :: ai_jn(0:n)
      complex(dp) :: ai_in(0:n)

      real(dp) :: ai_rn(0:2*(mrec+3))

!
!    Evaluate the integrals.
!

      call get_ai_in(n,z,phif,ai_in,ai_rn,mrec)
      if (n == 0) then
         ai_jn(0) = ai_j0(z,phif)
      else
         ai_jn(0) = ai_j0(z,phif)
         do i = 1,n,1
            ai_jn(i) = ai_in(i-1) + z*ai_jn(i-1)
         enddo
      endif

      end

!----------------------------------------------------------------

       function ai_i0(z,phif)
          use mod_precision


!
!    This is a function to evaluate the Aoki integral I0.
!

      implicit none
      complex(dp) :: ai_i0

!
!    Declare the simple variables.
!

      complex(dp) :: z,ai_l

      real(dp) :: phif

!
!    Evaluate the function.
!

      ai_i0 = cmplx(sin(phif), kind=dp) + z*cmplx(phif, kind=dp) & 
     &      - sqrt(cmplx(1.0_dp, kind=dp)-z)*sqrt(cmplx(1.0_dp, kind=dp)+z)* & 
     &        ai_l(z,phif)

      end

!----------------------------------------------------------------

       function ai_j0(z,phif)
          use mod_precision


!
!    This is a function to evaluate the Aoki integral J0.
!

      implicit none
      complex(dp) :: ai_j0

!
!    Declare the simple variables.
!

      complex(dp) :: z,ai_l

      real(dp) :: phif

!
!    Evaluate the function.
!

      ai_j0 = cmplx(phif, kind=dp) - cmplx(sin(phif), kind=dp)/(cmplx(cos(phif), kind=dp)-z) & 
     &  + z*ai_l(z,phif)/sqrt(cmplx(1.0_dp, kind=dp)-z)/sqrt(cmplx(1.0_dp, kind=dp)+z)

      end

!----------------------------------------------------------------

       function ai_l(z,phif)
          use mod_precision


!
!    This is a function to evaluate the Aoki integral L.
!

      implicit none
      complex(dp) :: ai_l

!
!    Declare the simple variables.
!

      complex(dp) :: z,c,s

      real(dp) :: phif

!
!    Evaluate the function.
!

      c = sqrt(cmplx(1.0_dp, kind=dp) - z)*cmplx(cos(0.5_dp*phif), kind=dp)
      s = sqrt(cmplx(1.0_dp, kind=dp) + z)*cmplx(sin(0.5_dp*phif), kind=dp)

      ai_l = log(c+s)-log(c-s)

      end

