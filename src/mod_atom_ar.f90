

   module mod_atom_ar
   
      use mod_precision
      use mod_const, only : mxz, matype, mbtype
!       use mod_all_scalar
!       use topologia
      
      implicit none
      
      
! zatype(z) -> type
! atype(type) -> z
      integer :: zatype(0:mxz)

!*** Sort codes for finding the input atom/bond type index from the internal z sorted atom/bond type index
      integer :: asort(0:matype-1), bsort(mbtype)

!*** The bond type pointers.
      integer btype(0:mxz,0:mxz)

      
!*** Errors in onsite charges.
      real(dp), allocatable :: dq(:)

!*** The masses of the atoms.
      real(dp), allocatable :: mass(:)
      
!*** Magnetic moment
      real(dp), allocatable :: mg(:),  demi(:) !mginert(:),
          
      
   
   end module mod_atom_ar