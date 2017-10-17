 
      subroutine fft2d(f,m,n1,n2,wk,mwk,isign)
          use mod_precision


!
!    This is a routine to evaluate the Fourier transform of a two dimensional function.
!

      implicit none

!
!    Declare the simple variables.
!

      integer m,n1,n2,mwk,isign
      integer i1,i2

!
!    Declare the arrays.
!

      complex(dp) :: f(m,n2)
      complex(dp) :: wk(mwk)

!
!    Evaluate the transform.
!

      if (mwk < n2) then
         write(6,'(''The work space is too small.'')')
         stop
      endif

      do i2 = 1,n2,1
         call four1(f(1,i2),n1,isign)
      enddo

      do i1 = 1,n1,1
         do i2 = 1,n2,1
            wk(i2) = f(i1,i2)
         enddo
         call four1(wk,n2,isign)
         do i2 = 1,n2,1
            f(i1,i2) = wk(i2)
         enddo
      enddo

      end

