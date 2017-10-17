 
      subroutine prodcoeff(a,na,b,nb,c,s)
          use mod_precision


!
!    This is a routine to return the coefficients of a polynomial formed
!     from the product of two other polynomials.
!

      implicit none

!
!    Declare the simple variables.
!

      real(dp) :: s

      integer na,nb
      integer i,j

!
!    Declare the arrays.
!

      real(dp) :: a(0:na)
      real(dp) :: b(0:nb)
      real(dp) :: c(0:na+nb)

!
!    Evaluate the coefficients.
!    NOTE: C is not initialised. This is to allow the accumulation
!          of products.
!

!D
!       print *, 'C from prodcoeff', c
      do i = 0,na,1
         do j = 0,nb,1
            c(i+j) = c(i+j) + s*a(i)*b(j)
         enddo
      enddo

      end
