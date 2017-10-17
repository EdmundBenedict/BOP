 
      subroutine orbrot(dr,orm,hcorm)
          use mod_precision


!
!    This is a routine to make the full orbital rotation matrix.
!

      implicit none

!
!    Declare the simple variables.
!

      complex(dp) :: ca,cb,sa,sb

      real(dp) :: alpha,beta,lng

      integer i1,i2

!
!    Declare the arrays.
!

      real(dp) :: orm(9,9)
      real(dp) :: hcorm(9,9)
      real(dp) :: r(5,5)
      real(dp) :: dr(3)
      real(dp) :: n(3)

!
!    Evaluate the rotation angles.
!

      n = -dr/sqrt(sum(dr*dr))

      lng = sqrt(sum(n(:2)*n(:2)))
      if (lng > 1.0e-9_dp) then
         alpha = atan2(n(2),n(1))
         beta = atan2(lng,n(3))
      else
         alpha = 0.0_dp
         beta = 0.0_dp
      endif

!     WRITE(6,'(''ALPHA = '',G12.5)') ALPHA*180.0D0/3.14159D0
!     WRITE(6,'(''BETA = '',G12.5)') BETA*180.0D0/3.14159D0

!
!    Clear the transformation matrix.
!

      do i1 = 1,9,1
         do i2 = 1,9,1
            orm(i1,i2) = 0.0_dp
         enddo
      enddo

!
!    s-s block.
!

      orm(1,1) = 1.0_dp

!
!    p-p block.
!

      ca = cos(alpha)
      sa = sin(alpha)
      cb = cos(beta)
      sb = sin(beta)

      orm(2,2) =  ca*cb
      orm(2,3) = -sa
      orm(2,4) =  ca*sb

      orm(3,2) =  sa*cb
      orm(3,3) =  ca
      orm(3,4) =  sa*sb

      orm(4,2) = -sb
      orm(4,3) =  0.0_dp
      orm(4,4) =  cb

!
!    d-d block.
!

      call maker(alpha,beta,r)

      do i1 = 1,5,1
         do i2 = 1,5,1
            orm(4+i1,4+i2) = r(i1,i2)
         enddo
      enddo

!
!    Evaluate the Hermitian conjugate of the transformation matrix.
!

      do i1 = 1,9,1
         do i2 = 1,9,1
            hcorm(i1,i2) = orm(i2,i1)
         enddo
      enddo

      end

