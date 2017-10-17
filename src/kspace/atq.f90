
   subroutine atqsym(nd, knh, krhod,  dq)
      use mod_precision
      use mod_const
      use mod_kspace, only : khpos, ateqv, symon
         
!    This is a routine to evaluate the charge on each
!     atomic site.
      implicit none
      
      integer, intent(in) ::  nd, knh
      real(dp), intent(in) :: krhod(knh)
      real(dp), intent(out) :: dq(nd)
      
      
      
      include "../Include/Atom.array"

      real(dp) :: zcore

      integer :: i,j,ip,ik,ierr, awt(nd)

      logical usepot

      
      
      dq = 0.0_dp
      awt = 0
      
      
      
      do j = 1,nd
         dq(j) = dq(j) + sum(krhod(khpos(j)+1:khpos(j+1)))
!          print *, 'ksum', j, sum(krhod(khpos(j)+1:khpos(j+1)))
      enddo

      
      if (symon) then
!  symmetrisation
!    gather on the unique atoms and generate weights
         do i = 1, nd
            if (ateqv(i) /= i) dq(ateqv(i)) = dq(ateqv(i)) + dq(i)
            awt(ateqv(i)) = awt(ateqv(i)) + 1
         end do
         
!    average   
         do i = 1, nd
            if (ateqv(i) == i) dq(i) = dq(i)/real(awt(i), dp)
         end do
      
      
!    scatter from the unique to the rest   
!      This loop is separate from the above in case the index of the unique atom is higher than the current
!       this is probably never true as I think spglib always makes the first encounter unique, but still...      
         do i = 1, nd
            if (ateqv(i) /= i) dq(i) = dq(ateqv(i))
         end do
      
      end if
      
   end subroutine atqsym

