 
      subroutine pascal()
          use mod_precision


          use mod_const

          use mod_srt

!
!    This is a routine for building Pascal's triangle.
!

      implicit none

!
!    Include constants.
!

!      include "Include/ALL.const"

!
!    Include arrays.
!

!      include "Include/SRT.array"

!
!    Declare the simple variables.
!

      integer m,i,j

!
!    Initialise the triangle.
!

      m = 2*mfdpol+3

      do i = 1,m,1
         do j = 1,m,1
            triang(j,i) = 0.0_dp
         enddo
      enddo
      triang(1,1) = 1.0_dp
      triang(2,1) = 1.0_dp
      triang(2,2) = 1.0_dp

!
!    Build the triangle.
!

      do i = 3,m,1
         triang(i,1) = 1.0_dp
         do j = 2,i,1
            triang(i,j) = triang(i-1,j-1) + triang(i-1,j)
         enddo
      enddo

      end

