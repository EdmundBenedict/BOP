 
       function getebsnull(efermi,nrec,mrec, diag,eigvec,kt)

!
!    This is a function to evaluate the band structure energy.
!
          use mod_precision

      implicit none
      real(dp) :: getebsnull

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

      getebsnull = 0.0_dp
      do i = 1,nrec+1,1
         getebsnull = getebsnull + (eigvec(1,i)**2)*diag(i)*theta(diag(i)-efermi,kt)
      enddo
      getebsnull = 2.0_dp*getebsnull

      end

