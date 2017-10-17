
   subroutine classic ( )
      use mod_precision
      use mod_const
      use mod_all_scalar      
      use mod_clock
      use mod_gsp, only : gsp
      use mod_conf, only : rtc, atom_conf_t, bond_conf_t, pwp_t
      use mod_atom_ar, only : asort, bsort, btype
      use mod_pft, only : tailed_fun
      
      implicit none

      include 'Include/Atom.array'
      include "Include/NebList.array"
      include "Include/NebListL.array"
      include "Include/ag.conn"
      include "Include/PosVel.array"
      include 'Include/Force.array'

      integer :: i, j, ja, bidx, bt, df
      real(dp) :: r, rs, vr(3), lam_j, lam_i, eenv_ij, epair_ij, lam(totnd), fv(0:1), f_i(3)
      type(atom_conf_t), pointer :: atm(:)
      type(bond_conf_t), pointer :: bnd(:)
      type(pwp_t), pointer :: pwp
      
      integer(8) :: c1,c2

      
      df = 0
      if (forces) df = 1
      
      epair = 0.0_dp
      eenv = 0.0_dp
      eclas = 0.0_dp

      atm => rtc % ham_conf % a
      bnd => rtc % ham_conf % b


      
!       if ( rtc % env /= 0) then
!          call system_clock(c1)
! 
!          
! !          bnd(1) % a = sqrt(bnd(0) % a * bnd(2) % a)
! 
! 
!          do i = 1, nne
!             lam_i = 0.0_dp
!             ja = aptrl(i)
!             j  = bptrl(ja)
!             do while (j /= eol)
!                r = sqrt(sum((ad(:,j) - ad(:,i))**2))
! !                   bidx = ind(i) + ind(j)
!                bidx = 2*ind(j)
!                if (r < bnd(bidx) % rcut) then
!                   lam_j = bnd(bidx) % c * exp(-bnd(bidx) % nu * (r-bnd(bidx)%rnu))
!                   if (.not. quiet) &
!                 & print '(12x,i3,x,a,2(x,f12.6))', map(j), atm(ind(j)) % symb, r, lam_j
!                   if (r > bnd(bidx) % r1 ) then
!                       print *, 'in tail'
!                      rs = (r - bnd(bidx) % r1)/(bnd(bidx) % rcut - bnd(bidx) % r1)
!                      lam_j = lam_j * ((6.0_dp*rs + 3.0_dp)*rs + 1.0_dp) * (1.0_dp - rs)**3
!                   end if
!                else
!                   lam_j = 0.0_dp
!                end if
!                lam_i = lam_i + lam_j
!                ja = ja + 1
!                j = bptrl(ja)
!             end do
!             lam(i) = atm(ind(i)) % e % lam_0 + lam_i ** (1.0_dp/atm(ind(i)) % e % m)
!             if (.not. quiet) print '(8x,i3,x,a,2(x,f12.6))', i, atm(ind(i)) % symb, lam(i), lam_i ** (1.0_dp/atm(ind(i)) % e % m)
!          end do
!          do  i = nne+1, totnd
!             lam( i )  = lam( map( i ) )
!          enddo
!          call system_clock(c2)
!          print *,'lam:',real(c2-c1,dp)/real(cr,dp),'s'
!       end if



      call system_clock(c1)
      
      do i = 1, nd
!          if (.not. quiet) write(6,'(12x,i3,x,a)') i, atm(ind(i)) % symb

!          write(6,'(a,2(x,i0))') 'i:', i, z(i)

         f_i = 0.0_dp
         
         ja = aptrl(i)
         j  = bptrl(ja)
         do while (j /= eol)
            
            bt = bsort(btype(z(i),z(j)))
            if (bt == -1) bt = bsort(btype(z(j),z(i)))
            pwp => bnd(bt) % pwp 
            
            
            
            vr = ad(:,j) - ad(:,i)
            r = sqrt(sum(vr*vr))

!             write(6,'(a,2(x,i0),4(x,f20.12))') '  j:', j, z(j), vr, r
            
            if (rtc % vpair == 1) then
               
               call tailed_fun(r, pwp % fp, pwp % tl, df, fv)
               
               epair = epair + fv(0)
               
!                if (forces) f_i = f_i + fv(1)*vr/r
               if (forces) fpp(:,i) = fpp(:,i) + fv(1)*vr/r
            end if
            
!             write(6,'(a,2(x,i0),3(x,f20.12))') '   j:', j, z(j), fv, epair
            ja = ja + 1
            j = bptrl( ja )
         end do
         
!          fpp(:,i) = fpp(:,i) + f_i
         
      end do

      


!       eenv  = 0.5_dp * eenv
      epair = 0.5_dp * epair
      eclas = eenv + epair
      
      call system_clock(c2)
      
      if (.not. quiet) print "('t(class@root):', g0)", real(c2-c1,dp)/real(cr,dp)

      
   end subroutine classic
   
   
   
!             if (.not. quiet) write(6,'(16x,i3,x,a,4(x,f6.3))') j, atm(ind(j)) % symb, ad(:,j) - ad(:,i), r
!             write(200+i,'(a,4(x,f12.6))') atm(ind(j)) % symb, (ad(:,j) - ad(:,i))/r, r
            
!             bidx = ind(i) + ind(j)
!             if (rtc % env == 1) then
!                if (r < bnd(bidx) % rcut) then
!                   eenv_ij = bnd(bidx)%a &
!                            & * exp(-0.5_dp * (lam(i) + lam(j)) * (r - (atm(ind(i))%e%r_core + atm(ind(j))%e%r_core)))/r
! !                   eenv_ij = exp(bnd(bidx)%a -0.5_dp * (lam(i) + lam(j)) * (r - (atm(ind(i))%e%r_core + atm(ind(j))%e%r_core)))/r
! 
! !                            print *, 'r - R',  (r - (atm(ind(i))%e%r_core + atm(ind(j))%e%r_core))
!                   if (r > bnd(bidx) % r1 ) then
!                      rs = (r - bnd(bidx) % r1)/(bnd(bidx) % rcut - bnd(bidx) % r1)
!                      eenv_ij = eenv_ij * ((6.0_dp*rs + 3.0_dp)*rs + 1.0_dp) * (1.0_dp - rs)**3
!                   end if
!                else
!                   eenv_ij = 0.0_dp
!                end if
!                eenv = eenv + eenv_ij
!             end if

