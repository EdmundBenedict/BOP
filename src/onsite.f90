 
         subroutine onsite(z,esz,epz,edz)
             use mod_precision


!
!    This is a routine to assign the on site energies. The values are in eV.
!

      implicit none

!
!    Declare the simple variables.
!

      real(dp) :: esz,epz,edz

      integer z,mxz

!
!    Define parameters.
!

      parameter(mxz = 103)

!
!    Declare the arrays.
!

      real(dp) :: es(0:mxz)
      real(dp) :: ep(0:mxz)
      real(dp) :: ed(0:mxz)

!
!    Declare the common block.
!

      common /ecom/es,ep,ed

!
!    Assign the values.
!

      esz = es(z)
      epz = ep(z)
      edz = ed(z)

      end

