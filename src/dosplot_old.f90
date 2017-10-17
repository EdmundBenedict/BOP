 
   subroutine dosplot()

      use mod_precision
      use mod_all_scalar
      use mod_const
      use ab_io          
      use mod_io, only : get_new_unit, close_unit
      use mod_g0n, only : g00
!       use mod_clock

      implicit none


      include "Include/Atom.array"


      

      real(dp) :: e, estep, ildos, ldos, emin, emax, tdos, itdos
      real(dp) :: numelsrt, numelnull

      integer :: i, ia, la, ne, dos_u, idos_u

      

      write(6,'("number of divisions?: ")',advance='no')
      read(5,*) ne

      write(6,'("adjusting energy window")')
!
!    Find the limits of the density of states.
!

      emin = 1.0e30_dp
      emax = -emin
      if (term == 1) then
         do ia = 1,nd
            call assoc_ab(ia)
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


      

      
!       write (6,'("generating dos ...")',advance='no')
      write (6,'("generating dos ...")')
      flush(6)

!       call system_clock(c1)

      do i = 0, ne-1
!          call system_clock(c)
!          write(6,'(x,i0)', advance='no') i+1
         write(6,'(x,i0)') i+1
         
         flush(6)

         e = emin + i*estep

         write( dos_u,'(x,e13.6,x)',advance='no') e
         write(idos_u,'(x,e13.6,x)',advance='no') e

         tdos = 0.0_dp
         itdos = 0.0_dp

         do ia = 1, nd
            call assoc_ab(ia)
            do la = 1, nchain
               if (term == 1) then
                  ldos = -aimag(g00(cmplx(e, kind=dp), & 
      &               arec(0:lchain(la),la),brec(0:lchain(la)+1,la), & 
      &               lchain(la),lainf(la),lbinf(la)))/pi
                  ildos = numelsrt(e,la)
               elseif ((term == 2).or.(term == 3)) then
                  ldos = -aimag(g00(cmplx(e,kt, kind=dp), & 
      &               arec(0:lchain(la),la),brec(0:lchain(la)+1,la), & 
      &               lchain(la),lainf(la),lbinf(la)))/pi
                  ildos = numelnull(e,lchain(la),mrec,diag(1,la), eigvec(1,1,la),kt)
               end if
               ldos = wt(la)*ldos
               ildos = wt(la)*ildos
               write( dos_u,'(x,e13.6)',advance='no')  ldos
               write(idos_u,'(x,e13.6)',advance='no') ildos
               tdos  = tdos  + ldos
               itdos = itdos + ildos
            end do
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



      end subroutine dosplot

