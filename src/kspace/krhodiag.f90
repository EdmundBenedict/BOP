   subroutine krhodiag(knh, nk, occ, kpsi,  krhod,kovl)
   
      use mod_precision
!       use mod_all_scalar
      use mod_const
      use mod_conf
!       use topologia, only : iproc, mpmap

      implicit none

      integer, intent(in) :: knh, nk
      real(dp), intent(in) :: occ(knh,nk)
      complex(dp), intent(in) :: kpsi(knh,knh,nk),kovl(knh,knh)
      real(dp), intent(out) :: krhod(knh)
      
      integer :: ik, i,l
      logical :: ovl
      
      ovl = rtc % ham_conf % ovl
      
   ! 
   !     forall (ik=mpmap(iproc)+1:mpmap(iproc+1), i=1:knh, ip=1:knh)
   !         krho(i,ip,ik) = occ(i,ik) *(real(kpsi(ip,i,ik))*real(kpsi(ip,i,ik)) + aimag(kpsi(ip,i,ik))*aimag(kpsi(ip,i,ik)))
   !     end forall
   !     
      krhod(:) = 0.0_dp
   
      do ik = 1,nk
         do i = 1, knh
            if (.not. ovl) then 
                krhod(:) = krhod(:) + occ(i,ik)*(real(kpsi(:,i,ik))* real(kpsi(:,i,ik)) &
                                        & + aimag(kpsi(:,i,ik))*aimag(kpsi(:,i,ik)))
            else
                do l =  1, knh                           
                    krhod(l) = krhod(l) + occ(i,ik)*(conjg(kpsi(l,i,ik))*sum(kovl(l,:)*kpsi(:,i,ik)))
                enddo
            endif
         end do
      end do


   end subroutine krhodiag