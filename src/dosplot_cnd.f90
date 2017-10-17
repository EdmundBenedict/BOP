 
   subroutine dosplot()

      use mod_precision
      use mod_all_scalar
      use mod_const
      use ab_io          
      use mod_io, only : get_new_unit, close_unit
      use mod_g0n, only : g00
!       use mod_clock
      use mod_conf, only : rtc
      implicit none


      include "Include/Atom.array"


      

      real(dp) :: e, estep, ildos, ldos, emin, emax, tdos, itdos
      real(dp), allocatable, dimension(:,:) :: dos_tbl, idos_tbl, sdos_tbl, sidos_tbl

      real(dp) :: numelsrt, numelnull

      integer :: i, j, ia, la, ne, dos_u, idos_u, nld, zmap(mxz)

      integer, allocatable :: offsets(:)

      write(6,'("number of divisions?: ")',advance='no')
      read(5,*) ne

      write(6,'("adjusting energy window")')
!
!    Find the limits of the density of states.
!

      emin = huge(0.0_dp)
      emax = -emin
      if (term == 1) then
         do ia = 1,nd
            call assoc_ab(ia)
            print *,'ia,nchain:',ia,nchain
            do la = 1,nchain
               emin = min((lainf(la)-2.0_dp*lbinf(la)),emin)
               emax = max((lainf(la)+2.0_dp*lbinf(la)),emax)
            end do
         end do
      elseif ((term == 2).or.(term == 3)) then
         do ia = 1,nd
            call assoc_ab(ia)
            do la = 1,nchain
               emin = min(diag(1,la),emin)
               emax = max(diag(lchain(la)+1,la),emax)
            end do
         end do
         emin = emin - etail
         emax = emax + etail
      endif

      estep = (emax-emin)/real(ne-1, dp)



      


      nld = 0
      do ia = 1, nd
         call assoc_ab(ia)
         nld = nld + nchain
      end do
      
      allocate(dos_tbl(ne,nld), idos_tbl(ne,nld))


      
      j = 1
      if (term == 1) then
         do ia = 1, nd
            call assoc_ab(ia)
            do la = 1, nchain   
               write(6,'(2(x,i0))') ia, la               
               flush(6)            
               do i = 1, ne
                  e = emin + (i-1)*estep
                               
                  ldos = -aimag(g00(cmplx(e, kind=dp), & 
      &               arec(0:lchain(la),la),brec(0:lchain(la)+1,la), & 
      &               lchain(la),lainf(la),lbinf(la)))/pi

                  ildos = numelsrt(e,la)
                  
                  ldos = wt(la)*ldos
                  ildos = wt(la)*ildos

                  dos_tbl(i,j) = ldos
                  idos_tbl(i,j) = ildos

               end do
               j = j + 1               
            end do
         end do
      elseif ((term == 2).or.(term == 3)) then
         do ia = 1, nd
            call assoc_ab(ia)
            do la = 1, nchain
               do i = 1, ne
                  e = emin + (i-1)*estep

                  ldos = -aimag(g00(cmplx(e,kt, kind=dp), & 
      &               arec(0:lchain(la),la),brec(0:lchain(la)+1,la), & 
      &               lchain(la),lainf(la),lbinf(la)))/pi

                  ildos = numelnull(e,lchain(la),mrec,diag(1,la), eigvec(1,1,la),kt)

                  ldos = wt(la)*ldos
                  ildos = wt(la)*ildos

                  dos_tbl(i,j) = ldos
                  idos_tbl(i,j) = ildos

               end do
               j = j + 1               
            end do
         end do
      end if


! printout full



      dos_u  = get_new_unit()
      idos_u = get_new_unit()


      open( dos_u, file='dos', action='write' )
      open(idos_u, file='idos', action='write')

      write( dos_u,'("ef: ",e13.6)', advance='no') lef
      write(idos_u,'("ef: ",e13.6)', advance='no') lef

      do ia = 1, nd
         call assoc_ab(ia)
         do la = 1, nchain
            write( dos_u,'(" a:",i0," l:",i0,6x)', advance='no') ia, la
            write(idos_u,'(" a:",i0," l:",i0,6x)', advance='no') ia, la
         end do
      end do
      write( dos_u,'(" total")')
      write(idos_u,'(" total")')


      


      do i = 1, ne
!          call system_clock(c)
!          write(6,'(x,i0)', advance='no') i
!          write(6,'(x,i0)') i
!          
!          flush(6)

         e = emin + (i-1)*estep

         write( dos_u,'(x,e13.6,x)',advance='no') e
         write(idos_u,'(x,e13.6,x)',advance='no') e

         tdos = 0.0_dp
         itdos = 0.0_dp
         
         do j = 1, nld
            write( dos_u,'(x,e13.6)',advance='no')  dos_tbl(i,j)
            write(idos_u,'(x,e13.6)',advance='no') idos_tbl(i,j)
            tdos  = tdos  + dos_tbl(i,j)
            itdos = itdos + idos_tbl(i,j)
         end do
!          tdos = tdos/real(nd,kind=dp)
!          itdos = itdos/real(nd,kind=dp)

         write( dos_u,'(2x,e13.6)')  tdos
         write(idos_u,'(2x,e13.6)') itdos
         
      end do

!       call system_clock(c2)
      write(6,'(/)')

      call close_unit( dos_u)
      call close_unit(idos_u)


! compact

      allocate(offsets(0:natype))

      
      nld = 0
      do i = 0, natype
         offsets(i) = nld
         nld = nld + nl(atype(i))
         zmap(atype(i)) = i         
      end do

      allocate(sdos_tbl(ne,nld), sidos_tbl(ne,nld))
      
      sdos_tbl = 0.0_dp
      sidos_tbl = 0.0_dp

      
      j = 1
      do ia = 1, nd
         call assoc_ab(ia)
         i = offsets(zmap(z(ia)))
         sdos_tbl(:,i+1:i+nchain) = sdos_tbl(:,i+1:i+nchain) + dos_tbl(:,j:j+nchain-1)
         sidos_tbl(:,i+1:i+nchain) = sidos_tbl(:,i+1:i+nchain) + idos_tbl(:,j:j+nchain-1)
         j = j + nchain         
      end do




! printout compact



       dos_u = get_new_unit()
      idos_u = get_new_unit()


      open( dos_u, file= 'dos_c', action='write')
      open(idos_u, file='idos_c', action='write')

      write( dos_u,'("ef: ",e13.6)', advance='no') lef
      write(idos_u,'("ef: ",e13.6)', advance='no') lef

      do i = 0, natype
         do la = 1, nl(atype(i))
            write( dos_u,'(x,a,x,a,9x)', advance='no') listsymb(atype(i)),orbs(llist(la,atype(i)))
            write(idos_u,'(x,a,x,a,9x)', advance='no') listsymb(atype(i)),orbs(llist(la,atype(i)))
         end do
      end do

      write( dos_u,'(" total")')
      write(idos_u,'(" total")')



      do i = 1, ne
!          call system_clock(c)
!          write(6,'(x,i0)', advance='no') i
!          write(6,'(x,i0)') i
!          
!          flush(6)

         e = emin + (i-1)*estep

         write( dos_u,'(x,e13.6,x)',advance='no') e
         write(idos_u,'(x,e13.6,x)',advance='no') e

         tdos = 0.0_dp
         itdos = 0.0_dp
         
         do j = 1, nld
            write( dos_u,'(x,e13.6)',advance='no')  sdos_tbl(i,j)
            write(idos_u,'(x,e13.6)',advance='no') sidos_tbl(i,j)
            tdos  =  tdos +  sdos_tbl(i,j)
            itdos = itdos + sidos_tbl(i,j)
         end do
!          tdos = tdos/real(nd,kind=dp)
!          itdos = itdos/real(nd,kind=dp)

         write( dos_u,'(2x,e13.6)')  tdos
         write(idos_u,'(2x,e13.6)') itdos
         
      end do

!       call system_clock(c2)
      write(6,'(/)')

      call close_unit( dos_u)
      call close_unit(idos_u)


      deallocate(dos_tbl, idos_tbl,sdos_tbl, sidos_tbl, offsets)

      end subroutine dosplot

