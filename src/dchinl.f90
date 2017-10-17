 
      subroutine dchinl(efermi,nrec,mrec,dchia,dchib,diag,eigvec,kt)

!
!    This is a routine to evaluate the derivatives of the suscetibilities.
!
          use mod_precision

      implicit none

!
!    Declare the simple variables.
!

      real(dp) :: efermi,kt
      real(dp) :: delta,ddelta
      real(dp) :: sum,wi,wj,ei,ej

      integer nrec,mrec
      integer n
      integer i,j

!
!    Declare the arrays.
!

      real(dp) :: dchia(0:mrec+2)
      real(dp) :: dchib(1:mrec+2)
      real(dp) :: diag(mrec+2)
      real(dp) :: eigvec(mrec+2,mrec+2)

!
!    Evaluate the derivatives of the susceptibilities.
!

      do n = 0,nrec

         sum = 0.0_dp
         do i = 1,nrec+1
            wi = eigvec(1,i)*eigvec(n+1,i)
            ei = diag(i)
            do j = 1,nrec+1
               wj = eigvec(1,j)*eigvec(n+1,j)
               ej = diag(j)
               if (abs(ei-ej) > 1.0e-5_dp*kt) then
                  sum = sum - wi*wj/(ei-ej)* & 
     &                  (delta(ei-efermi,kt)-delta(ej-efermi,kt))
               else
                  sum = sum - wi*wj*ddelta(ei-efermi,kt)
               endif
            enddo
         enddo

         dchia(n) = sum

      enddo


      do n = 1,nrec

         sum = 0.0_dp
         do i = 1,nrec+1
            wi = eigvec(1,i)*eigvec(n+1,i)
            ei = diag(i)
            do j = 1,nrec+1
               wj = eigvec(1,j)*eigvec(n,j)
               ej = diag(j)
               if (abs(ei-ej) > 1.0e-5_dp*kt) then
                  sum = sum - wi*wj/(ei-ej)* & 
     &                 (delta(ei-efermi,kt)-delta(ej-efermi,kt))
               else
                  sum = sum - wi*wj*ddelta(ei-efermi,kt)
               endif
            enddo
         enddo

         dchib(n) = sum

      enddo

      end

