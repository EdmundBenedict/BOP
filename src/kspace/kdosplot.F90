
   subroutine kdosplot()

      use mod_precision
      use mod_all_scalar
      use mod_const
      use mod_io, only : get_new_unit, close_unit
!       use mod_clock
      use mod_kspace
      use topologia, only : iproc, master, mpmap
#ifdef MPI
      use mpi, only : mpi_in_place, mpi_real8, mpi_integer, mpi_sum, mpi_min, mpi_max, mpi_comm_world
#endif

      implicit none


      include "../Include/Atom.array"
!       include "../Include/KHamilt.array"



      real(dp) :: e, estep, ildos, ldos, emin, emax, tdos, itdos
      real(dp), allocatable, dimension(:,:) :: dos_tbl, idos_tbl, sdos_tbl, sidos_tbl
      real(dp), pointer :: nullp => null()

      real(dp) :: numelsrt, numelnull

      integer :: i, j, ia, la, ne, dos_u, idos_u, nld, zmap(mxz)

      integer, allocatable :: offsets(:)
      integer :: ik, nma, nla, im, in, jd, jh, nstta,ib, ierr, isp


      if (iproc==master) then
        write(6,'("number of divisions?: ")',advance='no')
        read(5,*) ne

        write(6,'("adjusting energy window")')

      end if

#ifdef MPI
      call mpi_bcast(ne, 1, mpi_integer, master, mpi_comm_world, ierr)
#endif

      emin = huge(0.0_dp)
      emax = -emin

      isp = 1

      do ik = mpmap(iproc)+1, mpmap(iproc+1)
         emin = min(emin,enk(1,ik,isp))
         emax = max(emax,enk(knh,ik,isp))
      end do


#ifdef MPI
      call mpi_allreduce(mpi_in_place, emin, 1, mpi_real8, mpi_min, mpi_comm_world, ierr)
      call mpi_allreduce(mpi_in_place, emax, 1, mpi_real8, mpi_max, mpi_comm_world, ierr)
#endif


      estep = (emax-emin)/real(ne-1, dp)


      nld = 0
      do ia = 1, nd
         nld = nld + nl(z(ia))
      end do

      allocate(dos_tbl(ne,nld), idos_tbl(ne,nld))

      dos_tbl(:,:) = 0.0_dp
      idos_tbl(:,:) = 0.0_dp

      print *,'kt:', kt
!       stop
      print *,'calc all pdos'

      do ik = mpmap(iproc)+1, mpmap(iproc+1)
        jd = 1
        jh = 0
        do ia = 1, nd
    !          print *, 'ia:', ia
            call states(z(ia),nla,nstta,llista)
            do la = 1,nla
    !             print *, '    la:', ia
                nma = 2*llista(la)+1
                do i = 1, ne
    !                print *, '        i:', i
                    e = emin + (i-1)*estep
                    do in = 1,khpos(nd+1)
                        do im = 1, nma
                            dos_tbl(i,jd) = dos_tbl(i,jd) &
                                        & -aimag(1.0_dp/cmplx(e-enk(in,ik,isp),kt,kind=dp))/pi &
                                        & * wtk(ik) &
                                        & *( real(kpsi(jh+im,in,ik,isp))* real(kpsi(jh+im,in,ik,isp)) &
                                        & + aimag(kpsi(jh+im,in,ik,isp))*aimag(kpsi(jh+im,in,ik,isp)))
    !                         print *, jh, khpos(ia), kpsi(jh+im,in,ik)
                        end do
                    end do
                end do
                jh = jh + nma
                jd = jd + 1
            end do
        end do
     end do



#ifdef MPI
      if (iproc /= master) then
        call mpi_reduce(dos_tbl, nullp, ne*nld , mpi_real8, mpi_sum, master, mpi_comm_world, ierr)
      else
        call mpi_reduce(mpi_in_place, dos_tbl, ne*nld , mpi_real8, mpi_sum, master, mpi_comm_world, ierr)
      endif
#endif

!
!    do j = 1, nld
!       do i = 1, ne
!
!       end do
!    end do


! print *, 'z(ia),nla,nstta,llista',z(ia),nla,nstta,llista,khpos(ia)
!
!
!

! printout full

     if (iproc /= master) then
        deallocate(dos_tbl,idos_tbl)
        return
     end if

      print *,'print all pdos'

      dos_u  = get_new_unit()
      idos_u = get_new_unit()


      open( dos_u, file='kdos', action='write' )
      open(idos_u, file='kidos', action='write')

      write( dos_u,'("ef: ",e13.6)', advance='no') lef
      write(idos_u,'("ef: ",e13.6)', advance='no') lef

      do ia = 1, nd
         do la = 1, nl(z(ia))
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


      print *,'compact dos'

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
         nla = nl(z(ia))
         i = offsets(zmap(z(ia)))
         sdos_tbl(:,i+1:i+nla) = sdos_tbl(:,i+1:i+nla) + dos_tbl(:,j:j+nla-1)
         sidos_tbl(:,i+1:i+nla) = sidos_tbl(:,i+1:i+nla) + idos_tbl(:,j:j+nla-1)
         j = j + nla
      end do




! printout compact


      print *,'print compact dos'

       dos_u = get_new_unit()
      idos_u = get_new_unit()


      open( dos_u, file= 'kdos_c', action='write')
      open(idos_u, file='kidos_c', action='write')

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

      end subroutine kdosplot

