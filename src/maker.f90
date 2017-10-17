 
      subroutine maker(alpha,beta,r)
          use mod_precision


!
!    This is a routine to make the rotation matrix for d orbitals.
!

      implicit none

!
!    Declare the simple variables.
!

      real(dp) :: alpha,beta,root
      real(dp) :: sb,cb,s2b,c2b

!
!    Declare the arrays.
!

      real(dp) :: z(5,5)
      real(dp) :: y(5,5)
      real(dp) :: r(5,5)

!
!    Make the transformation matrix.
!

      z(1,1) =  1.0_dp
      z(1,2) =  0.0_dp
      z(1,3) =  0.0_dp
      z(1,4) =  0.0_dp
      z(1,5) =  0.0_dp

      z(2,1) =  0.0_dp
      z(2,2) =  cos(2.0_dp*alpha)
      z(2,3) = -sin(2.0_dp*alpha)
      z(2,4) =  0.0_dp
      z(2,5) =  0.0_dp

      z(3,1) =  0.0_dp
      z(3,2) =  sin(2.0_dp*alpha)
      z(3,3) =  cos(2.0_dp*alpha)
      z(3,4) =  0.0_dp
      z(3,5) =  0.0_dp

      z(4,1) =  0.0_dp
      z(4,2) =  0.0_dp
      z(4,3) =  0.0_dp
      z(4,4) =  cos(alpha)
      z(4,5) = -sin(alpha)

      z(5,1) =  0.0_dp
      z(5,2) =  0.0_dp
      z(5,3) =  0.0_dp
      z(5,4) =  sin(alpha)
      z(5,5) =  cos(alpha)

      sb = sin(beta)
      cb = cos(beta)
      s2b = 2.0_dp*sb*cb
      c2b = cb*cb - sb*sb
      root = sqrt(3.0_dp)*0.5_dp

      y(1,1) = 0.5*(3.0_dp*cb*cb-1.0_dp)
      y(1,2) = root*sb*sb
      y(1,3) = 0.0_dp
      y(1,4) =-root*s2b
      y(1,5) = 0.0_dp

      y(2,1) = root*sb*sb
      y(2,2) = 1.0_dp-0.5_dp*sb*sb
      y(2,3) = 0.0_dp
      y(2,4) = 0.5_dp*s2b
      y(2,5) = 0.0_dp

      y(3,1) = 0.0_dp
      y(3,2) = 0.0_dp
      y(3,3) = cb
      y(3,4) = 0.0_dp
      y(3,5) = sb

      y(4,1) = root*s2b
      y(4,2) =-0.5_dp*s2b
      y(4,3) = 0.0_dp
      y(4,4) = c2b
      y(4,5) = 0.0_dp

      y(5,1) = 0.0_dp
      y(5,2) = 0.0_dp
      y(5,3) =-sb
      y(5,4) = 0.0_dp
      y(5,5) = cb

      call dmatml(r,5,z,5,5,5,y,5,5)

      end

