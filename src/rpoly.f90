 
      subroutine rpoly(c,n,z,f,nd)
          use mod_precision


!
!    This is a function to evaluate a polynomial and its derivatives.
!

      implicit none

!
!    Declare the simple variables.
!

      real(dp) :: z
      real(dp) :: fac

      integer n,i,nd,id,nnd

!
!    Declare the arrays.
!

      real(dp) :: c(0:n)
      real(dp) :: f(0:nd)

!
!    Evaluate the polynomial, and its derivatives.
!

      f(0) = c(n)
      do id = 1,nd,1
         f(id) = 0.0_dp
      enddo

      do i = n-1,0,-1
         nnd = min(nd,n-i)
         do id = nnd,1,-1
           f(id) = f(id)*z + f(id-1)
         enddo
         f(0) = f(0)*z+c(i)
      enddo

      fac = 1.0_dp
      do id = 2,nd,1
         fac = fac*real(id, dp)
         f(id) = fac*f(id)
      enddo

      end

