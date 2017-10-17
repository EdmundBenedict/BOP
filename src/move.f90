 
      subroutine move(ad,nd,del,alpha,dmdl,dmu)
          use mod_precision


!
!    This is an ion mover routine.
!

      implicit none

      integer , intent(in) :: nd
      real(dp), intent(inout) :: ad(3,nd), dmu
      real(dp), intent(in) :: del(3,nd), dmdl(3,nd), alpha
      
      
      real(dp) :: dad
      real(dp), parameter :: maxdad = 0.1_dp
      integer :: i,j
      

      do i = 1,nd
         do j = 1,3
            dad = alpha*del(j,i)
            if (dad < 0.0_dp) then
               dad = max(dad,-maxdad)
            else
               dad = min(dad,maxdad)
            endif
            ad(j,i) = ad(j,i) + dad           
!             dmu = dmu + dmdl(j,i)*dad
         enddo
      enddo

      end subroutine move

