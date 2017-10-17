
      subroutine makechi(efermi)
         use mod_precision
         use mod_const
         use mod_all_scalar
         use mod_chi
         use ab_io
         use mod_funptr
         use topologia, only : iproc, mpmap, aproc, master


         implicit none
      
         include "Include/Atom.array"

         real(dp), intent(inout) :: efermi
         integer :: ia, la, nmax
         
         real(dp) :: chia_loc





         do ia = mpmap(iproc)+1, mpmap(iproc+1)

            call assoc_ab(ia)
            call assoc_chi(ia)

 ! Evaluate the susceptibilities
            chiaia(ia) = 0.0_dp
            if (term == 1) then
                do la = 1,nchain
                    nmax = lchain(la)
                    call getchisrt(efermi,nmax+1,la)
                    call dchisr(efermi,nmax+1,la) ! needed for dmdl
                    chiaia(ia) = chiaia(ia) + chia(0,la)*wt(la)
                end do
            elseif ((term == 2).or.(term == 3)) then
               do la = 1,nchain
                  nmax = lchain(la)
                  call getchinull(efermi,nmax,mrec,chia(0,la), & 
      &                         chib(1,la),diag(1,la), & 
      &                         eigvec(1,1,la),kt)
! dchi needed for dmdl     
                  call dchinl(efermi,nmax,mrec,dchia(0,la), & 
      &                     dchib(0,la),diag(1,la), & 
      &                     eigvec(1,1,la),kt)
                  chiaia(ia) = chiaia(ia) + chia(0,la)*wt(la)
               enddo
            endif
         end do


      end subroutine makechi

