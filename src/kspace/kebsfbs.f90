
   subroutine kebsfbs()
      use mod_precision
      use mod_all_scalar
      use mod_const
      use topologia
!           use spread_stuff, only : spread_aveclusiz
      use mod_clock
      use mod_conf, only : rtc
      use mod_kspace
      use mod_atom_ar, only : dq, mg

!       use mpi

      implicit none

      include "../Include/Atom.array"



      integer :: it, mxit, ia, mgit, mit
      real(dp) :: etott
      real(qp) :: avdde, sumz, merr
      real(qp), allocatable :: de_prev(:), dq_prev(:), dec(:),  td(:), dql(:)
      real(dp), allocatable :: mg_o(:)

      character(len=25) :: ffmt

      integer,parameter :: mixing_type = 2
      integer :: loc_clusiz
      integer(8) :: c1,c2

      integer :: i
      real(dp) :: ef,nf

      integer :: ib, ierr
      real(dp) :: dqmax_old, aa, aamin, aamax

      ef = 0.0_dp
      !       mxit = 500
      mxit = rtc % bop_conf % mxit
      if (mag) then
        mgit = rtc % bop_conf % mgit
        merr = rtc % bop_conf % merr
        mit = 1
        if (merr < 0.0_dp) mgit = 1
      end if

      allocate(de_prev(nd), dq_prev(nd), td(nd), dec(nd), dql(nd))
      dec = real(de(:nd),qp)
      sumz = sum( (/ (real(zc(z(ia)), qp), ia = 1,nd) /) )



      

      call kbldh(.false.)
      if ((qerr < 0.0_dp).or.(mixing_type == 0)) mxit = 1
      if (.not. quiet) then
         write(9,"('Imposing LCN ...')")
         write(ffmt,"(a,i0,a)") '(a,":",', nd, '(x,f22.14))'
      end if

      if (mag) allocate(mg_o(nd))

      do
         if (mag) then
            if (.not. quiet) write(6,'(/,"mit:",i3)') mit
            mg_o = mg

            do ia = 1, nd
               dem(ia) = 0.5_dp*istn(z(ia))*mg(ia)
            end do

         end if

!          if linmix
         dqmax_old = 0.0_dp
         aamax=1.0_dp
         aamin=1.0e-4_dp
         aa=10.0_dp
!          endif

         it = 1
         do
      !           print "('sumzde,de:',2(x,f22.14))", sum( (/ (de(ia)*zc(z(ia)), ia = 1,nd) /) ), de(9)
            de(:nd) = real(dec - sum( (/ (dec(ia)*zc(z(ia)), ia = 1,nd) /) ) / sumz, dp)
      !           print "('bde:',x,f22.14)", de(9)
            if (.not. quiet) then
               write(6,'(/,"it:",i3)') it
               write(6, ffmt) 'de ', de(:nd)
               write(6, ffmt) 'dem', dem(:nd)
            end if

            call spnkspc()

      !                      dec(:nd) = de(:nd)
            dql = real(dq,qp)


      !
            if (abs(dqmax) <= qerr .or. it == mxit) exit

!             if (it > 1) then
!                do ia = 1, nd
!       !     if (.not. quiet) print *, 'dif:',   (dec(ia) - de_prev(ia))/ (dq(ia) - dq_prev(ia) ), abs(dq(ia) - dq_prev(ia))
!                   if ( abs(dec(ia) - de_prev(ia)) > 6.0_qp*abs(dql(ia) - dq_prev(ia))) then
!                      td(ia) = -dql(ia)
!       !    if (.not. quiet) print *,red,'potential daner', (dec(ia) - de_prev(ia)), (dq(ia) - dq_prev(ia)),endc
!                   else
!                      td(ia) = dql(ia) * (dec(ia) - de_prev(ia)) / ( dql(ia) - dq_prev(ia) )
!                   endif
!                end do
!                de_prev = dec
!                dq_prev = dql
!                dec = dec - td
!       !               write(6,"('td,nde:',2(x,f22.14))") td(9),de(9)
!             else
!                de_prev = dec
!                dq_prev = dql
!                dec = dec + dql
!             end if

            if (     (dqmax_old*dqmax > 0.0_dp &
            &  .and. abs(dqmax_old) > abs(dqmax)) &
            &  .or. (dqmax_old*dqmax < 0.0_dp)) then
               aa=aa*abs(dqmax_old)/abs(dqmax_old-dqmax)
            end if

            aa=max(aamin,aa)
            aa=min(aamax,aa)

            if (abs(dqmax) <= qerr) aa=aamin

            dec(:nd) = de(:nd) + aa*dq(:nd)


            dqmax_old = dqmax

            it = it + 1
         end do

         if (.not. mag) exit
         mgmax = maxval(abs(mg - mg_o))

         if (.not. quiet) then
            write(6,ffmt) 'mgo', mg_o
            write(6,ffmt) 'mg ', mg
            write(6,*) 'mgmax:', mgmax
         end if
         if (mgmax <= merr .or. mit == mgit) exit
         mit = mit + 1
      end do

      if (mag) deallocate(mg_o)
      deallocate(de_prev, dq_prev, td, dec, dql)


      if (forces) call bsf()

   end subroutine kebsfbs








!    subroutine kebsfbs()
!       use mod_precision
!       use mod_all_scalar
!       use mod_const
!       use topologia
! !           use spread_stuff, only : spread_aveclusiz
!       use mod_clock
!       use mod_conf, only : rtc
!       use mod_kspace
!
!       implicit none
!
!       include "../Include/Atom.array"
!
!
!
!       integer :: it, mxit, ia
!       real(dp) :: etott
!       real(dp) :: sumz
!
!       character(len=20) :: ffmt
!
!       integer,parameter :: mixing_type = 2
!
!       integer :: i
!       real(dp) :: ef, dqmax_old, aa, aamin, aamax
!
!
!
! !       mxit = 500
!       mxit = rtc % bop_conf % mxit
!       it = 1
!       sumz = sum( (/ (zc(z(ia)), ia = 1,nd) /) )
!
!
!       dqmax=0.0_dp
!       aamax=1.0_dp
!       aamin=1.0e-4_dp
!       aa=10.0_dp
!
!
!
!       call kbldh()
!
!       if ((qerr < 0.0_dp).or.(mixing_type == 0)) mxit = 1
!       if (.not. quiet) then
!          write(9,"('Imposing LCN ...')")
!          write(ffmt,"(a,i0,a)") '(a,":",', nd, '(x,f22.14))'
!       end if
!
!
!
!
!       do
!          if (it > 1) then
!             if (     (dqmax_old*dqmax > 0.0_dp &
!               &  .and. abs(dqmax_old) > abs(dqmax)) &
!               &  .or. (dqmax_old*dqmax < 0.0_dp)) then
!                aa=aa*abs(dqmax_old)/abs(dqmax_old-dqmax)
!             end if
!
!             aa=max(aamin,aa)
!             aa=min(aamax,aa)
!
!             if (abs(dqmax) <= qerr) aa=aamin
!
!             de(:nd) = de(:nd) + aa*dq(:nd)
!          end if
!
!          de(:nd) = de(:nd) - sum( (/ (de(ia)*zc(z(ia)), ia = 1,nd) /) ) / sumz
!
!          dqmax_old = dqmax
!
!          if (.not. quiet) then
!             write(6,'(/,"it:",i3)') it
!             write(6, ffmt) 'de', de(:nd)
!          end if
!
!          call kbldhdo()
!          call kdiag()
!          call fndocc(ef,totnia)  !    Find the Fermi energy.
!          call atq() ! find atom charges
!
!          if (.not. quiet) then
!             write(6, ffmt) 'dq', dq(:nd)
!             write(6,'("dqmax:",x,f22.14,2x,"ef:",x,f22.14,2x,"numel:",x,f22.14,x,f22.14)') dqmax, ef, totnia, sum(dq(:nd))
!          end if
!
!          call eval_kens(ef)
!
!
!          if (abs(dqmax) <= qerr .or. it == mxit) exit
!
!          it = it + 1
!       end do
!
!
!    end subroutine kebsfbs
!
!
!
! !             if (.not. quiet) then
! !                etott = eprom + ebond + epair - eatom - uent
! !                write(6,'(/,"it:",i3)') it
! !                write(6, ffmt) 'de', de(:nd)
! !                write(6, ffmt) 'dec', dec
! !                write(6, ffmt) 'dq', dq(:nd)
! !    !             write(6, ffmt) 'ch', chian(:nd)
! !    !             write(6,"('de,dq:',2(x,f22.14))") de(9),dq(9)
! !                write(6, '("ef:",2(x,f22.16))') ef, totnia, dqmax
! !                write(6, '("etot:",x,f22.16)') etott
! !                write(6, '("eprom:",x,f22.16)') eprom
! !                write(6, '("eatom:",x,f22.16)') eatom
! !                if (abs(dqmax) <= qerr) then
! !                   write(9,"('c')",advance='no')
! !                else if (it == 1) then
! !                   write(9,"('s')",advance='no')
! !                else if (it == mxit) then
! !                   write(9,"('x')",advance='no')
! !                else
! !                   write(9,"('i')",advance='no')
! !                end if
! !
! !                write(9,"(2(x,f22.14),x,i3)") etott, dqmax, it
! !
! !             end if
