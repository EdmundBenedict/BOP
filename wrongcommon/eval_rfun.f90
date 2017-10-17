

   subroutine eval_rfun(z1, z2, r, d, rf,tb)
   
   !relies on correct repeating info in b%h for ss, pp, dd, ...
   
      use mod_precision
      use mod_const
      use mod_conf
      use mod_atom_ar, only : btype, bsort
      use mod_pft, only : tailed_fun
      
      implicit none
      
      integer, intent(in) :: z1, z2, d
      real(dp), intent(in) :: r
      real(dp), intent(out) :: rf(0:d,14)

      logical :: rev
      integer :: oi, oj, li, lj, bi, bt, zbt, sli, slj,tb
      
      type(bond_conf_t), pointer :: b
      type(hop_t), pointer :: h
      
      real(dp) :: fv(0:d),a2,v
      
      rev = .false.
   
      zbt = btype(z1,z2)
      bt = bsort(zbt)
      
!       print *, 'zbt,bt',bt,zbt,z1,z2
      
      if (bt == -1) then
         rev = .true.
         bt = bsort(btype(z1,z2))
      end if
   
      if (tb == 0) then
          b => rtc % ham_conf % b(bt)
      else
           b => rtc % tbham_conf % b(bt)
      endif
!       print *, 'bt', bt,btype(z1,z2)
      rf = 0.0_dp
      
!       print *, 'dtau',rtc % ham_conf % b(bt)%h(7)%dtau
      
      if (.not.rev) then
         do oj = 1, b % atom1 % b % norbs
            lj = b % atom1 % b % orblist(oj)
            do oi = 1, b % atom0 % b % norbs
               li = b % atom0 % b % orblist(oi)
               if (li <= lj) then
                  do bi = 0,li
                     h => b % h(hsidx(li,lj,bi))
!                       print *, b % name, ' ', hnam(hsidx(sli,slj,bi)), ' ', li, lj, bi                       
                     if (tb==2) then
                        a2 = h%sc%a(2)
                        v = h % v
                        h%v = h%oscr
                        h%sc%a(2) = h%dtau
                     endif
                  
                     call tailed_fun(r, h % sc, h % tl, d, fv)
                     rf(0:d,hlidx(li,lj,bi)) = h % v * fv(0:d)
                     
                     if (tb == 2) then
                        h%sc%a(2) = a2
                        h%v = v
                     endif

!                      print *, 'ENDS',r,fv, h % v 
                  end do
               else
                  do bi = 0,lj
                     h => b % rev % h(hsidx(lj,li,bi))
                        
                     if (tb==2) then
                        a2 = h%sc%a(2)
                        v = h % v
                        h%v = h%oscr
                        h%sc%a(2) = h%dtau
                     endif
                  
                     call tailed_fun(r, h % sc, h % tl, d, fv)
                     rf(0:d,hlidx(li,lj,bi)) = h % v * fv(0:d)
                     
                     if (tb == 2) then
                        h%sc%a(2) = a2
                        h%v = v
          
                     endif
                  end do
               end if 
            end do
         end do 
      else
         stop 'The intelligent handling of nonrepeating hamiltonian input '&
         & //'data is not completed yet. Please provide data for all bond types'
      
      end if
   
   end subroutine eval_rfun