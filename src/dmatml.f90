 
      subroutine dmatml(ab,lxab,a,lxa,l1a,l2a,b,lxb,l2b)
          use mod_precision


!
!    This is a routine to multiply together two real matrices.
!

      implicit none

!
!    Declare the simple variables.
!

      integer lxab
      integer lxa,l1a,l2a
      integer lxb,l2b
      integer i1,i2,i3

!
!    Declare the arrays.
!

      real(dp) :: ab(lxab,l2b)
      real(dp) :: a(lxa,l2a)
      real(dp) :: b(lxb,l2b)

!
!    Carry out the multiplication.
!

      do i1 = 1,l1a,1
         do i2 = 1,l2b,1
            ab(i1,i2) = 0.0_dp
            do i3 = 1,l2a,1
               ab(i1,i2) = ab(i1,i2) + a(i1,i3)*b(i3,i2)
            enddo
         enddo
      enddo

      end

