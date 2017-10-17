 
      subroutine utranv(m1,h,m2,work,mxn,n1,n2)
          use mod_precision


!
!    This is a routine to carry out a unitary transformation
!     on a matrix of vectors.
!

      implicit none

!
!    Declare the simple variables.
!

      integer mxn,n1,n2
      integer i1,j1,i2,j2
      integer k

!
!    Declare the arrays.
!

      real(dp) :: m1(mxn,mxn)
      real(dp) :: m2(mxn,mxn)
      real(dp) :: h(3,mxn,mxn)
      real(dp) :: work(mxn,mxn)

!
!    Carry out the transformation.
!

      do k = 1,3,1

         do i2 = 1,n2,1
            do j1 = 1,n1,1
               work(j1,i2) = 0.0_dp
               do i1 = 1,n1,1
                  work(j1,i2) = work(j1,i2) + m1(i1,j1)*h(k,i1,i2)
               enddo
            enddo
         enddo

         do j2 = 1,n2,1
            do j1 = 1,n1,1
               h(k,j1,j2) = 0.0_dp
               do i2 = 1,n2,1
                  h(k,j1,j2) = h(k,j1,j2) + work(j1,i2)*m2(i2,j2)
               enddo
            enddo
         enddo

      enddo

      end

