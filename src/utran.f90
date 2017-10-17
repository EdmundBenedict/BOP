 
      subroutine utran(m1,h,m2,work,mxn,n1,n2,nlo,nhi,mrec,mxnnb,ja)

!
!    This is a routine to carry out a unitary transformation
!     on a matrix.
!
          use mod_precision

      implicit none

!
!    Declare the simple variables.


      integer mxn,n1,n2,mrec,nlo,nhi
      integer i1,i2,j1,j2,n,mxnnb,ja

!
!    Declare the arrays.
!

      real(dp) :: m1(mxn,mxn)
      real(dp) :: m2(mxn,mxn)
      real(dp) :: h(2*mrec+2,mxn,mxnnb,mxn)
      real(dp) :: work(mxn,mxn)

!
!    Carry out the transformation.
!

      do n = nlo,nhi,1

         do i2 = 1,n2,1
            do j1 = 1,n1,1
               work(j1,i2) = 0.0_dp
               do i1 = 1,n1,1
                  work(j1,i2) = work(j1,i2) + m1(i1,j1)*h(n,i1,ja,i2)
               enddo
            enddo
         enddo

         do j2 = 1,n2,1
            do j1 = 1,n1,1
               h(n,j1,ja,j2) = 0.0_dp
               do i2 = 1,n2,1
                  h(n,j1,ja,j2) = h(n,j1,ja,j2) + work(j1,i2)*m2(i2,j2)
               enddo
            enddo
         enddo

      enddo

      end

