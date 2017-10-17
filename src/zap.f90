 
      subroutine zap(za,z,ad,d,nd)
          use mod_precision


!
!    This is a routine to remove an atom of type ZA from the list
!     of atoms. The atom is cosen at random.
!

      implicit none

!
!    Declare the simple variables.
!

      real(dp) :: ran1

      integer za,nd
      integer na,i,j,nzap
      integer idum

!
!    Define the parameters.
!

!
!    Declare the arrays.
!

      real(dp) :: d(3,nd)
      real(dp) :: ad(3,nd)

      integer z(nd)

!
!    Assign data values.
!

      data idum/1/

!
!    Save variables.
!

      save idum

!
!    Count the number of atoms of type ZA in the unit cell.
!

      na = 0
      do i = 1,nd,1
         if (za == z(i)) na = na + 1
      enddo

!
!    Zap an atom.
!

      if (na > 0) then
         nzap = 1 + nint(real(na-1, dp)*ran1(idum))
         na = 0
         i = 1
         do while((i <= nd).and.(na < nzap))
            if (za == z(i)) na = na + 1
            i = i + 1
         enddo
         do j = i-1,nd-1,1
            z(j) = z(j+1)
            d(1,j) = d(1,j+1)
            d(2,j) = d(2,j+1)
            d(3,j) = d(3,j+1)
            ad(1,j) = ad(1,j+1)
            ad(2,j) = ad(2,j+1)
            ad(3,j) = ad(3,j+1)
         enddo
         nd = nd - 1
      endif

      end

