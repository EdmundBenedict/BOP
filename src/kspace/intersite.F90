

   subroutine bsinter()
      use mod_precision
      use mod_all_scalar
      use mod_const
      use mod_kspace
      use mod_atom_ar, only : btype
      use topologia
   
      implicit none

      integer :: isp, ik, i, j, k
      real(dp) :: ro, bsi, bso
      complex(dp) :: t
      complex(dp), allocatable :: h(:,:)
      
      allocate(h(knh,knh))
      
      bsi = 0.0_dp
      bso = 0.0_dp
      
      do isp = 1, nsp
         do ik = 1, nk
            call shiftons(isp, knh, kham(1,1,ik),  h)
            do i = 1, knh
               do j = 1, knh
                  ro = 0.0_dp
                  
                  do k = 1, knh
                     ro = ro + occ(k,ik,isp) &
                              & *( real(kpsi(k,i,ik,isp))* real(kpsi(k,j,ik,isp)) &
                              & + aimag(kpsi(k,i,ik,isp))*aimag(kpsi(k,j,ik,isp)))
                  end do

                  
                  if (i /= j) then
                     bsi = bsi + h(i,j)*ro
                  else
!                      print *, isp, ik, i, j, ro, h(i,j)
                     bso = bso + h(i,j)*ro
                  end if
                  
               end do
            end do
         end do
      end do
      
      deallocate(h)
      
      print *, 'bsi:', bsi
      print *, 'bso:', bso
      
   end subroutine bsinter