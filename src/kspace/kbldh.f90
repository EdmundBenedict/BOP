 
   subroutine kbldh(addons)

      use mod_precision
      use mod_all_scalar
      use mod_const
      use mod_kspace
      use mod_conf
      use topologia, only : iproc, mpmap
      use mod_atom_ar, only : btype, bsort
      use mod_io, only : print_c
      
      implicit none

      logical, intent(in) :: addons ! whether to add onsite energies
      
      include "../Include/Atom.array"
      include "../Include/NebList.array"
      include "../Include/PosVel.array"

      
      
      real(dp) :: dea,phase
      
      integer :: i,j,ik
      integer :: za,zb,mapb
      integer :: ja,ja0,ia,ib,la,lb
      integer :: nstta,nsttb
      integer :: llist1(ml)
      integer :: llist2(ml)
      integer :: nl1,nl2

      real(dp) :: dr(3)
      real(dp) :: ksubh(mxnstat,mxnstat)
      complex(dp) :: fac
      type(ham_conf_t),pointer :: hamtb
      
      real(dp) :: scfcut(14)
      
      logical :: ovl

      scfcut = 1.0_dp

      dea = 0.0_dp
      ovl = rtc % ham_conf % ovl

      do ik = mpmap(iproc)+1, mpmap(iproc+1)                     
         kham(:,:,ik) = cmplx(0.0_dp,0.0_dp,kind=dp)
!         write(6,'("ik,k:",i2,3(x,f6.3))') ik, kk(:,ik)
         do ia = 1,nd
!             print *,'    ia:',ia
            za = z(ia)
            nstta = nstt(za)
            nl1 = nl(za)
            llist1 = llist(:,za)
            if (addons) dea = de(ia)
            
!
!          For each neighboring ion (B) do ...
!

            ja = aptr(ia)
            ja0 = ja - 1
            do while (bptr(ja) /= eol)               
               ib = bptr(ja)
               mapb = map(ib)
               zb = z(mapb)
               nsttb = nstt(zb)
               nl2 = nl(zb)
               llist2 = llist(:,zb)
               
               dr = ad(:,ia)-ad(:,ib)
               
               phase = sum(kk(:,ik)*dr)
               
               fac = cmplx(cos(phase),sin(phase),kind=dp)

               hamtb => rtc % ham_conf 

               call matel(za,nstta,zb,nsttb,dr,ksubh,dea,scfcut,hamtb,.false.)

               do la = 1,nstta
                  do lb = 1,nsttb
                     kham(khpos(ia)+la, khpos(mapb)+lb, ik) = kham(khpos(ia)+la, khpos(mapb)+lb, ik) + fac*ksubh(la,lb)
                  end do
               end do
               
               
               ja = ja+1
            end do

         end do
      end do
      

      if (ovl) then
       
        do ik = mpmap(iproc)+1, mpmap(iproc+1)                     
          kovl(:,:,ik) = cmplx(0.0_dp,0.0_dp,kind=dp)
         
  !         write(6,'("ik,k:",i2,3(x,f6.3))') ik, kk(:,ik)
          do ia = 1,nd
  !             print *,'    ia:',ia
              za = z(ia)
              nstta = nstt(za)
              nl1 = nl(za)
              llist1 = llist(:,za)
              dea = 0.0_dp
              
  !
  !          For each neighboring ion (B) do ...
  !

              ja = aptr(ia)
              ja0 = ja - 1
              do while (bptr(ja) /= eol)               
                ib = bptr(ja)
                mapb = map(ib)
                zb = z(mapb)
                nsttb = nstt(zb)
                nl2 = nl(zb)
                llist2 = llist(:,zb)
                
                dr = ad(:,ia)-ad(:,ib)
                
                phase = sum(kk(:,ik)*dr)
                
                fac = cmplx(cos(phase),sin(phase),kind=dp)

                
                call matel(za,nstta,zb,nsttb,dr,ksubh,dea,scfcut,hamtb,.true.)

                do la = 1,nstta
                    do lb = 1,nsttb
                      kovl(khpos(ia)+la, khpos(mapb)+lb, ik) = kovl(khpos(ia)+la, khpos(mapb)+lb, ik) + fac*ksubh(la,lb)
                    end do
                end do
                
                
                ja = ja+1
              end do

          end do
        end do
      
      endif

   end subroutine kbldh
