 
      subroutine inv3x3(r)
          use mod_precision


          use mod_all_scalar

!
!    This is a subroutine that inverts a 3x3 matrix.
!

      implicit none

!      include 'Include/ALL.scalar'

!
!    Declare the simple variables.
!

      real(dp) :: det, inv_det

!
!    Declare the arrays.
!

      real(dp) :: r(3,3)
      real(dp) :: cof(3,3)
!       integer i,j
      real(dp) :: tin

!
!    Evalauate the determinant.
!


         cof(1,1) =  r(2,2)*r(3,3)-r(2,3)*r(3,2)
         cof(1,2) =  r(2,3)*r(3,1)-r(2,1)*r(3,3)
         cof(1,3) =  r(2,1)*r(3,2)-r(2,2)*r(3,1)

      det = r(1,1)*cof(1,1) + r(1,2)*cof(1,2) + r(1,3)*cof(1,3)

      inv_det = 1.0_dp/det
!
!    Evaluate the cofactors.
!
      if (abs(det) > 1.0e-10_dp) then ! tiny(0.0_dp)
!             write(6,'(3(3(x,f8.5),/))') ((r(i,j), j=1,3), i=1,3) 
!          cof(1,1) =  r(2,2)*r(3,3)-r(2,3)*r(3,2)
!          cof(1,2) =  r(2,3)*r(3,1)-r(2,1)*r(3,3)
!          cof(1,3) =  r(2,1)*r(3,2)-r(2,2)*r(3,1)

         cof(2,1) =  r(1,3)*r(3,2)-r(1,2)*r(3,3)
         cof(2,2) =  r(1,1)*r(3,3)-r(1,3)*r(3,1)
         cof(2,3) =  r(1,2)*r(3,1)-r(1,1)*r(3,2)

         cof(3,1) =  r(1,2)*r(2,3)-r(1,3)*r(2,2)
         cof(3,2) =  r(1,3)*r(2,1)-r(1,1)*r(2,3)
         cof(3,3) =  r(1,1)*r(2,2)-r(1,2)*r(2,1)

!
!       Evaluate the inverse.
!

         r(1,1) = cof(1,1)*inv_det
         r(1,2) = cof(2,1)*inv_det
         r(1,3) = cof(3,1)*inv_det
         r(2,1) = cof(1,2)*inv_det
         r(2,2) = cof(2,2)*inv_det
         r(2,3) = cof(3,2)*inv_det
         r(3,1) = cof(1,3)*inv_det
         r(3,2) = cof(2,3)*inv_det
         r(3,3) = cof(3,3)*inv_det


      else

         if (.not. quiet) write(6,'(''Unable to invert matrix. Setting it to I'')')
!          write(6,'(''Unable to invert matrix. Setting it to I'')')
!          write(6,'(3(3(x,f8.5),/))') ((r(i,j), j=1,3), i=1,3) 
         
         r(1,1) = 1.0_dp
         r(1,2) = 0.0_dp
         r(1,3) = 0.0_dp
         r(2,1) = 0.0_dp
         r(2,2) = 1.0_dp
         r(2,3) = 0.0_dp
         r(3,1) = 0.0_dp
         r(3,2) = 0.0_dp
         r(3,3) = 1.0_dp

      endif

      end

