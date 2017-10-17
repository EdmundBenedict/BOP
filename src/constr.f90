 
      subroutine constr(f,nd,cnst_a,cnst_v,cnst_n,mcnst_n)
          use mod_precision


!
!    This is a routine to constrain the forces on atoms.
!

      implicit none
      real(dp) :: dot

      integer cnst_n,mcnst_n,ia
      integer nd
      integer i


      real(dp) :: f(3,nd)
      real(dp) :: cnst_v(3,mcnst_n)

      integer :: cnst_a(mcnst_n)

!
!    Constrain atoms that need to be constrained.
!

      do i = 1,cnst_n
         ia = cnst_a(i)
         f(1:3,ia) = f(1:3,ia) - cnst_v(1:3,i)*sum(f(1:3,ia)*cnst_v(1:3,i)) 
      enddo

      end

