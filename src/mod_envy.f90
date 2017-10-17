
   module mod_envy
!    envy : env yukawa
   

      use mod_precision
      use mod_const
      use mod_all_scalar      
      use mod_clock
      use mod_conf, only : rtc, atom_conf_t, bond_conf_t
      
      implicit none
      
   contains
   
   subroutine envy()
   
   
      include "Include/NebList.array"
      include "Include/NebListL.array"
      include "Include/ag.conn"
      include "Include/PosVel.array"

      integer :: i, j, ja, bidx
      real(dp) :: r, rs, lam_j, lam_i, eenv_ij, epair_ij, lam(totnd)
      type(atom_conf_t), pointer :: atm(:)
      type(bond_conf_t), pointer :: bnd(:)
      integer(8) :: c1,c2


      epair = 0.0_dp
      eenv = 0.0_dp
      eclas = 0.0_dp

      atm => rtc % ham_conf % a
      bnd => rtc % ham_conf % b


      
      if ( rtc % env /= 0) then
         call system_clock(c1)

         
!          bnd(1) % a = sqrt(bnd(0) % a * bnd(2) % a)


         do i = 1, nne
            lam_i = 0.0_dp
            ja = aptrl(i)
            j  = bptrl(ja)
            do while (j /= eol)
               r = sqrt(sum((ad(:,j) - ad(:,i))**2))
!                   bidx = ind(i) + ind(j)
               bidx = 2*ind(j)
               if (r < bnd(bidx) % rcut) then
                  lam_j = bnd(bidx) % c * exp(-bnd(bidx) % nu * (r-bnd(bidx)%rnu))
                  if (.not. quiet) &
                & print '(12x,i3,x,a,2(x,f12.6))', map(j), atm(ind(j)) % symb, r, lam_j
                  if (r > bnd(bidx) % r1 ) then
                      print *, 'in tail'
                     rs = (r - bnd(bidx) % r1)/(bnd(bidx) % rcut - bnd(bidx) % r1)
                     lam_j = lam_j * ((6.0_dp*rs + 3.0_dp)*rs + 1.0_dp) * (1.0_dp - rs)**3
                  end if
               else
                  lam_j = 0.0_dp
               end if
               lam_i = lam_i + lam_j
               ja = ja + 1
               j = bptrl(ja)
            end do
            lam(i) = atm(ind(i)) % e % lam_0 + lam_i ** (1.0_dp/atm(ind(i)) % e % m)
            if (.not. quiet) print '(8x,i3,x,a,2(x,f12.6))', i, atm(ind(i)) % symb, lam(i), lam_i ** (1.0_dp/atm(ind(i)) % e % m)
         end do
         do  i = nne+1, totnd
            lam( i )  = lam( map( i ) )
         enddo
         call system_clock(c2)
         print *,'lam:',real(c2-c1,dp)/real(cr,dp),'s'
      end if



      call system_clock(c1)
      do i = 1, nd
!          if (.not. quiet) write(6,'(12x,i3,x,a)') i, atm(ind(i)) % symb

!          write(200+i,'(a)') atm(ind(i)) % symb
         ja = aptrl(i)
         j  = bptrl(ja)
         do while (j /= eol)
            
            r = sqrt(sum((ad(:,j) - ad(:,i))**2))
!             if (.not. quiet) write(6,'(16x,i3,x,a,4(x,f6.3))') j, atm(ind(j)) % symb, ad(:,j) - ad(:,i), r
!             write(200+i,'(a,4(x,f12.6))') atm(ind(j)) % symb, (ad(:,j) - ad(:,i))/r, r
            
            bidx = ind(i) + ind(j)
            if (rtc % env == 1) then
               if (r < bnd(bidx) % rcut) then
                  eenv_ij = bnd(bidx)%a &
                           & * exp(-0.5_dp * (lam(i) + lam(j)) * (r - (atm(ind(i))%e%r_core + atm(ind(j))%e%r_core)))/r
!                   eenv_ij = exp(bnd(bidx)%a -0.5_dp * (lam(i) + lam(j)) * (r - (atm(ind(i))%e%r_core + atm(ind(j))%e%r_core)))/r

!                            print *, 'r - R',  (r - (atm(ind(i))%e%r_core + atm(ind(j))%e%r_core))
                  if (r > bnd(bidx) % r1 ) then
                     rs = (r - bnd(bidx) % r1)/(bnd(bidx) % rcut - bnd(bidx) % r1)
                     eenv_ij = eenv_ij * ((6.0_dp*rs + 3.0_dp)*rs + 1.0_dp) * (1.0_dp - rs)**3
                  end if
               else
                  eenv_ij = 0.0_dp
               end if
               eenv = eenv + eenv_ij
            end if

            if (rtc % vpair == 1) then
               if (r < bnd(bidx) % pair % rcut) then
                  epair_ij = bnd(bidx) % pair % v &
                           & * gsp(r, bnd(bidx)%pair%r0, bnd(bidx)%pair%rc, bnd(bidx)%pair%n, bnd(bidx)%pair%nc)
                  if (r > bnd(bidx) % pair % r1 ) then
                     rs = (r - bnd(bidx) % pair % r1)/(bnd(bidx) % pair % rcut - bnd(bidx) % pair % r1)
                     epair_ij = epair_ij * ((6.0_dp*rs + 3.0_dp)*rs + 1.0_dp) * (1.0_dp - rs)**3
                  end if
               else
                  epair_ij = 0.0_dp
               end if
               epair = epair + epair_ij
            end if
            ja = ja + 1
            j = bptrl( ja )
         end do
      end do

      call system_clock(c2)
!       print *,'class',real(c2-c1,dp)/real(cr,dp)


      eenv  = 0.5_dp * eenv
      epair = 0.5_dp * epair
      eclas = eenv + epair

   end subroutine envy
   
   
   
   
   
   
   
   
   
   
   
   end module mod_envy