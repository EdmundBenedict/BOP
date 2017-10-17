
   subroutine get_dq2chia(efermi)
      use mod_precision
      use mod_const
      use mod_all_scalar, only : term, nsp, kt
      use mod_chi
      use ab_io
      use mod_funptr
      use topologia, only : iproc, mpmap, aproc, master
      use mod_atom_ar, only : dq
      
      implicit none
   
      include "Include/Atom.array"

      real(dp), intent(in) :: efermi
      integer :: ia, la, nmax, isp
      real(dp) :: chiaia

      
      do isp = 1, nsp  
         do ia = mpmap(iproc)+1, mpmap(iproc+1)

            call assoc_ab(ia,isp)
            call assoc_chi(ia,isp)
            
   ! Evaluate the susceptibilities
            chiaia = 0.0_dp
            if (term == 1) then
                  do la = 1,nchain
                     nmax = lchain(la)
                     call getchisrt(efermi,nmax+1,la)
! dchi needed for dmdl                     
!                      call dchisr(efermi,nmax+1,la) 
                     chiaia = chiaia + chia(0,la)*wt(la)
                  end do
            elseif ((term == 2).or.(term == 3)) then
               do la = 1,nchain
                  nmax = lchain(la)
                  call getchinull(efermi,nmax,mrec,chia(0,la), & 
      &                         chib(1,la),diag(1,la), & 
      &                         eigvec(1,1,la),kt)
   ! dchi needed dmdl      
!                   call dchinl(efermi,nmax,mrec,dchia(0,la), & 
!          &                     dchib(0,la),diag(1,la), & 
!          &                     eigvec(1,1,la),kt)
                  chiaia = chiaia + chia(0,la)*wt(la)
               enddo
            endif
            
            if (isp == 1) then
                dq2chia(ia) = chiaia
            else ! isp == 2 and mag == true
                dq2chia(ia) = dq2chia(ia) + chiaia
            end if
            
            if (isp == nsp) then
                if (isp == 1) dq2chia(ia) = 2*dq2chia(ia)
                
                dq2chia(ia) = dq(ia)/dq2chia(ia)
            end if
            
         end do
      end do   
      
   end subroutine get_dq2chia

