 
       function numelnull(efermi,nrec,mrec, diag,eigvec,kt)
!
!    This is a function to evaluate the number of electrons in a given state.
!
          use mod_precision

      implicit none
      real(dp) :: numelnull

!
!    Declare the simple variables.
!

      real(dp) :: efermi,theta,kt

      integer mrec,nrec
      integer i

!
!    Declare the arrays.
!

      real(dp) :: diag(mrec+2)
      real(dp) :: eigvec(mrec+2,mrec+2)

!
!    Evaluate the number of electrons.
!

      numelnull = 0.0_dp
      do i = 1,nrec+1,1
         numelnull = numelnull + (eigvec(1,i)*eigvec(1,i))*theta(diag(i)-efermi, kt)
      enddo
      numelnull = 2.0_dp*numelnull

      end

