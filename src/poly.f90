 
      subroutine poly(c,n,z,f,nd)
          use mod_precision


!
!    This is a function to evaluate a polynomial with real coefficients
!     but for complex values fo the argument.
!

      implicit none

!
!    Declare the simple variables.
!

      complex(dp) :: z

      real(dp) :: fac

      integer n,i,nd,id,nnd

!
!    Declare the arrays.
!

      real(dp) :: c(0:n)
      complex(dp) :: f(0:nd)

!
!    Evaluate the polynomial, and its derivatives.
!

      f(0) = cmplx(c(n), kind=dp)
      do id = 1,nd,1
         f(id) = cmplx(0.0_dp, kind=dp)
      enddo

      do i = n-1,0,-1
         nnd = min(nd,n-i)
         do id = nnd,1,-1
           f(id) = f(id)*z + f(id-1)
         enddo
         f(0) = f(0)*z+cmplx(c(i), kind=dp)
      enddo

      fac = 1.0_dp
      do id = 2,nd,1
         fac   = fac*real(id, dp)
         f(id) = cmplx(fac, kind=dp)*f(id)
      enddo

      end
!
!
!**********************************************************************
!
!
      subroutine poly1(c,n,z,f)
          use mod_precision


!
!    This is a function to evaluate a polynomial with real coefficients
!     but for complex values fo the argument.
!

      implicit none

!
!    Declare the simple variables.
!

      complex(dp) :: f
      complex(dp) :: z

      real(dp) :: fac

      integer n,i,nd,id

!
!    Declare the arrays.
!

      real(dp) :: c(0:n)

!
!    Evaluate the polynomial, and its derivatives.
!

      f = cmplx(c(n), kind=dp)

      do i = n-1,0,-1
         f = f*z+cmplx(c(i), kind=dp)
      enddo

      end
