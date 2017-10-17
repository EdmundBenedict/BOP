 
      subroutine pbc(nd,a,rlim, ad)
          use mod_precision


!
!    This is a subroutine to impose periodic boundary conditions.
!

      implicit none

      integer, intent(in) :: nd,rlim(3)

      real(dp), intent(inout) :: ad(3,nd)
      real(dp), intent(in) :: a(3,3)
      
      real(dp) :: b(3,3), d(3)
      integer i,j

!
!    Evaluate the transformation that takes cartesian
!     coordinates into reduced coordinates.
!


      b = a

      call inv3x3(b)

!
!       Impose periodic boundary conditions on each ion.
!

      do i = 1,nd
         call mul3x3(b,ad(1,i),d)
         where (rlim > 0) d = modulo(d,1.0_dp)               
         call mul3x3(a,d,ad(1,i))
      enddo


      end subroutine pbc

