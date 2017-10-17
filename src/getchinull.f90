 
      subroutine getchinull(efermi,nrec,mrec,chia,chib, diag,eigvec,kt)

!
!    This is a routine to evaluate the susceptibilities.
!
          use mod_precision

      implicit none

!
!    Declare the simple variables.
!

      real(dp) :: efermi,kt,theta
      real(dp) :: delta
      real(dp) :: sum,wi,wj,ei,ej

      integer nrec,mrec
      integer n
      integer i,j

!
!    Declare the arrays.
!

      real(dp) :: chia(0:mrec+2)
      real(dp) :: chib(1:mrec+2)
      real(dp) :: diag(mrec+2)
      real(dp) :: eigvec(mrec+2,mrec+2)

!
!    Evaluate the susceptibilities.
!

      do n = 0,nrec

         sum = 0.0_dp
!        DO I = 1,NREC,1
         do i = 1,nrec+1
            wi = eigvec(1,i)*eigvec(n+1,i)
            ei = diag(i)
!           DO J = I+1,NREC+1
            do j = 1,nrec+1,1
               wj = eigvec(1,j)*eigvec(n+1,j)
               ej = diag(j)
               if (abs(ei-ej) > 1.0e-9_dp*kt) then
                  sum = sum - wi*wj/(ei-ej)* & 
     &                  (theta(ei-efermi,kt)-theta(ej-efermi,kt))
!                 WRITE(9,'(I3,I3,G14.5,G14.5,G14.5)') I,J,WI*WJ,
!    +             (THETA(EI-EFERMI,KT)-THETA(EJ-EFERMI,KT))/
!    +             (EI-EJ),SUM
               else
                  sum = sum - wi*wj*delta(ei-efermi,kt)
!                 WRITE(9,'(I3,I3,G14.5,G14.5,G14.5)') I,J,WI*WJ,
!    +                  DELTA(EI-EFERMI,KT),SUM
               endif
            enddo
         enddo

!        CHIA(N) = 2.0D0*SUM
         chia(n) = sum

      enddo


      do n = 1,nrec

         sum = 0.0_dp
         do i = 1,nrec+1
            wi = eigvec(1,i)*eigvec(n+1,i)
            ei = diag(i)
            do j = 1,nrec+1
               wj = eigvec(1,j)*eigvec(n,j)
               ej = diag(j)
               if (abs(ei-ej) > 1.0e-9_dp*kt) then
                  sum = sum - wi*wj/(ei-ej)* & 
     &                 (theta(ei-efermi,kt)-theta(ej-efermi,kt))
               else
                  sum = sum - wi*wj*delta(ei-efermi,kt)
               endif
            enddo
         enddo

         chib(n) = sum

      enddo

      end

