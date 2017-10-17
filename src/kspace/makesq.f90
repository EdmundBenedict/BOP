   subroutine makesq()
   use mod_precision
   use mod_kspace, only : nk, knh, kpsi, ccw, occ

   implicit none

   integer :: i, ik, in, j
   complex(dp) :: alpha, beta
   complex(dp), allocatable :: opsi(:,:)
   real(dp) :: q
   alpha = 1.0_dp
   beta = 1.0_dp

   ccw = 0.0_dp

   allocate(opsi(knh,knh))

   

   do ik = 1, nk
      do in = 1, knh
         opsi(:,in) = kpsi(:,in,ik)*occ(in,ik) 
      end do
      call zgemm('n', 'c', knh, knh, knh, &
               & alpha, opsi, knh, &
               &        kpsi(1,1,ik), knh, &
               & beta,  ccw,  knh)
            
   end do

   write(6,'("ccwd: ",4("|",9(x,f11.8)))') (real(ccw(i,i),dp), i = 1,36)

   do i = 0, 3
      q = 0.0_dp
      do j = i*9+1, (i+1)*9
         q = q + real(ccw(j,j),dp)
      end do
      print*,i, q
   end do
!    write(6,'("ccwds: ",4(x,f11.8))') (sum(((/real(ccw(i*j:(i+1)*j,i*j:(i+1)*j),dp), j=1,9/))), i = 0,3)


   deallocate(opsi)



   end subroutine makesq