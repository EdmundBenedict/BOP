 
      function eprommix(ia,ef,isp)
      
      use mod_precision
      use mod_all_scalar
      use mod_const
      use ab_io
      use mod_funptr
      use mod_chi
!
!    This is a subroutine to evaluate the promotion energy.
!

      implicit none
      
      integer, intent(in) :: ia, isp ! isp is only needed to keep to the feprom interface
      real(dp), intent(in) :: ef

      include "Include/Atom.array"
      include "Include/NebList.array"

      real(dp) :: eprommix
      procedure(real(dp)) :: numelnull

      integer nmax, j
      integer nstta,nla,la,lb
      integer za

!
!    Evaluate the on site energies.
!

      eprommix = 0.0_dp

      za = z(ia)
      call states(za,nla,nstta,llista)

      do la = 1,nstta
         nmax = lchain(la)

!       Find the diagonal components of the on site density matrix.!
         if (term == 1) then
            bo(la,la) = numelsrt(ef,la)
         elseif ((term == 2).or.(term == 3)) then
            bo(la,la) = numelnull(ef,nmax,mrec, diag(1,la),eigvec(1,1,la),kt)
         endif
         
!       Find the off diagonal elements of the on site density matrix.
         do lb = 1,nstta  
            if (la /= lb) then
               bo(la,lb) = chia(0,la) * darec(0,la,lb,1)
               do j = 1,nmax
                  bo(la,lb) = bo(la,lb) + chia(j,la)*darec(j,la,lb,1) + 2*chib(j,la)*dbrec(j,la,lb,1)
               enddo
               if (term == 1) bo(la,lb) = bo(la,lb) + 2*chib(nmax+1,la) * dbrec(nmax+1,la,lb,1)
               bo(la,lb) = -2*bo(la,lb)
            endif
         enddo
      enddo

      do la = 1,nstta
         do lb = 1,nstta
            eprommix = eprommix + darec(0,la,lb,1)*bo(lb,la)
         enddo
      enddo

      end function eprommix

