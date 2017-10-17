 
      subroutine getocc(n,occ)
          use mod_precision


          use mod_all_scalar

          use mod_const

!
!    This is a routine to find the occupancies.
!

      implicit none

!
!    Include constants.
!

!      include "Include/ALL.const"

!
!    Include scalars.
!

!      include "Include/ALL.scalar"

!
!    Include arrays.
!

      include "Include/Atom.array"

!
!    Declare the simple variables.
!

      real(dp) :: nelec
      real(dp) :: occmin,occmax,docc

      integer ia
      integer n,i

!
!    Declare the local arrays.
!

      real(dp) :: occ(n)

!
!    Find the occupancy.
!

!     WRITE(6,'(''Finding occupancies.'')')

      if (n == 1) then
         nelec = 0.0_dp
         do ia = 1,nd,1
            
!             write(6,*) z(ia), zc(z(ia))
            
            nelec = nelec + zc(z(ia))
            
         enddo
         occ(1) = nelec
      else
         occmax = -1.0e-2_dp*real(nd, dp)
         do ia = 1,nd,1
            occmax = occmax + real(2*nstt(z(ia)), dp)
         enddo
         occmin = 1.0e-2_dp*real(nd, dp)
         docc = (occmax-occmin)/real(n-1, dp)
         occ(1) = occmin
         do i = 2,n,1
            occ(i) = occ(i-1) + docc
         enddo
      endif

!     write(6,*) occ(:), size(occ)

      end

