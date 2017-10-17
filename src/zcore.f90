 
       function zcore(z)
          use mod_precision


!
!    This is a routine to return the core charge for an element.
!

      implicit none
      real(dp) :: zcore

!
!    Declare the simple variables.
!

      integer z,mxz

!
!    Define the parameters.
!

      parameter (mxz = 103)

!
!    Declare the arrays.
!

      real(dp) :: zc(0:mxz)

!
!    Declare the common blocks.
!

      common /zcom/zc

!
!    Find the core charge.
!

      zcore = zc(z)

      end

