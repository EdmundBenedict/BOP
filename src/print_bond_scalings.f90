

   subroutine print_bond_scalings()
      use mod_precision
      use mod_const
      use mod_conf
      use mod_io  
      use mod_pft
      
      implicit none
      
      type(ham_conf_t), pointer :: hc
      type(bond_conf_t), pointer :: b
      type(hop_t), pointer :: h
      
      integer :: i, j, oi, oj, li, lj, bi, u, n
      character(len=3) :: hsymb
      character(len=2) :: si, sj
      
      real(dp) :: r, s, d, e, fv(0:2)
      
      n = 500
      s = 1.5_dp 
      e = maxval( rtc % rcut) + 0.2_dp
      d = (e - s) / real(n,dp)
      
      u = get_new_unit()
      
      open(u, file= 'bond_scls.out', action='write')
      
      write(u, '("#",7x,"r",12x)',advance = 'no' )
      
      hc => rtc % ham_conf
      
      do i = 0, hc % nb - 1
         b => hc % b(i)
         do oj = 1, b % atom1 % b % norbs
            lj = b % atom1 % b % orblist(oj)
            sj = b % atom1 % symb
            do oi = 1, b % atom0 % b % norbs
               li = b % atom0 % b % orblist(oi)
               si = b % atom1 % symb
               if (li > lj) cycle
               do bi = 0, li 
                  hsymb = hnam(hsidx(li,lj,bi))
                  write(u, "(3(9x,a2,'-',a2,x,a3,x, a))",advance = 'no') &
                     & si, sj, hsymb, 'v ', si, sj, hsymb, "v'", si, sj, hsymb, 'v"'
               end do
            end do
         end do 
      end do   
      
      hsymb = 'pwp'
      write(u, "(3(8x,a2,'-',a2,x,a3,x, a))") &
         &  si, sj, hsymb, 'v', si, sj, hsymb,"v'", si, sj, hsymb, 'v"'
      
      
      do j = 0, n - 1
         r = s + j * d
         write(u,'(x,f20.12,x)', advance='no') r
         do i = 0, hc % nb - 1
            b => hc % b(i)
            do oj = 1, b % atom1 % b % norbs
               lj = b % atom1 % b % orblist(oj)
               do oi = 1, b % atom0 % b % norbs
                  li = b % atom0 % b % orblist(oi)
                  if (li > lj) cycle
                  do bi = 0, li 
                     h => b % h(hsidx(li,lj,bi))
                     call tailed_fun(r, h % sc, h % tl, 2, fv)
                     write(u, '(3(x,f20.12))',advance = 'no') fv * h % v
                  end do
               end do
            end do 
            
            call tailed_fun(r, b % pwp % fp, b % pwp % tl, 2, fv)
            write(u, '(3(x,f20.12))',advance = 'no') fv
            
         end do
         write(u,'()')
      end do
      
      
      
      call close_unit(u)

   end subroutine print_bond_scalings
   