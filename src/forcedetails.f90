
   module mod_forcedetails

      use mod_precision
      use mod_const
      use mod_all_scalar
      use mod_conf
      use ab_io
      use mod_ham
      use mod_chi
      use topologia, only : iproc, mpmap, aproc, master, nproc
      use mod_clock
      use mod_atom_ar, only : asort, bsort, btype
      use mod_pft, only : tailed_fun
      use mod_io

      implicit none

      private

      include "Include/Atom.array"
      include "Include/NebList.array"
      include "Include/PosVel.array"





      logical :: inited = .false.


      real(dp) :: disl_core(2), strl = huge(0.0_dp)
!       integer, allocatable :: disl_atoms(:)



      real(dp), allocatable, dimension(:,:,:,:,:) :: rho, gdh, frl
      real(dp), allocatable, dimension(:,:,:) :: fbb, fbp, fbt
      real(dp), allocatable :: fat(:,:)

      integer, allocatable :: cls(:,:), pcl(:)

      public :: print_forcedetails

   contains

!
!    subroutine select_disl_atoms()
!
!       integer :: ia, j
!       integer, allocatable :: tmp(:)
!
!
!       allocate(tmp(nd))
!
!       j = 0
!       do ia = 1, nd
!          if (sum((ad(1:2,i) - disl_core)**2) < 2*rcutl) then
!             j = j + 1
!             tmp(j) = i
!          end if
!       end do
!
!       allocate(disl_atoms(j))
!       disl_atoms = tmp(:j)
!       deallocate(tmp)
!
!
!    end subroutine select_disl_atoms
!

   subroutine print_forcedetails()

!       will have to reload the neighlists cause the current after getetot is the pairpot

!       a logical :: fs ardument shall be passedd and there shalt be 2 calls from getetto 1 after getebsfs and one after cassic doing diff things of cource
!       then gatherv mozhe da e polezno
      use mpi

      include "Include/Moment.array"
      include "Include/Force.array"
      include "Include/ag.conn" ! for totstr


      integer :: ia, ib, jb, za, zb, la, lb, lla0, nla, nlb, nstta, nsttb, isp, ia1, ia2, nma, nmax, ma, mib, jd
      integer :: bidx, bt, na, nb
      real(dp) :: r, dr(3), fv(0:1), ferr, dda, ddb
      type(atom_conf_t), pointer :: atm(:)
      type(bond_conf_t), pointer :: bnd(:)
      type(pwp_t), pointer :: pwp
      integer(8) :: c1,c2, c1_m, c2_m
      real(dp) :: rhol(mxnstat,mxnstat), frll(4,mxnstat,mxnstat), fab(3), fabl(3), fatl(3)

      real(dp) :: grad(3,mxnstat,mxnstat)

      character(len=80) :: filename
      integer :: u

      integer :: usz, ierr
      integer, allocatable :: mpcnt(:), szs(:), displs(:)


      real(dp) :: scfcut(14),dsf(14,3)


      if (.not. inited) call init()

!       Let it overwrite the old file on every step. In this way only the relaxed state will be left
!       if (totstr == strl)
!       strl = totstr


      call system_clock(c1)

      scfcut = 1.0_dp
      dsf = 0.0_dp

      atm => rtc % ham_conf % a
      bnd => rtc % ham_conf % b
           


      if (iproc /= master) then
         ia1 = mpmap(iproc)+1   ! do not use these ia1 ia2 as loop range. the master will be wrong and calc all atoms
         ia2 = mpmap(iproc+1)
      else
         ia1 = 1
         ia2 = nd
      end if

      allocate(   rho(mxnstat, mxnstat, nsp, mxnnb, ia1:ia2), &
               &  gdh(3, mxnstat, mxnstat, mxnnb, ia1:ia2), &
               &  frl(4, mxnstat, mxnstat, mxnnb, ia1:ia2), &
               &  fbb(                  4, mxnnb, ia1:ia2), &
               &  fbp(                  4, mxnnb, ia1:ia2), &
               &  fbt(                  4, mxnnb, ia1:ia2), &
               &  cls(                     mxnnb, ia1:ia2), &
               &  pcl(                            ia1:ia2)) ! for reference with the proper force from forcebop.


!       real(dp) :: rho(mxnstat,mxnstat), grdmat(3,mxnstat,mxnstat), frl(mxnstat,mxnstat,3), fbb(3), fbp(3), fbt(3)


      rho = 0.0_dp
      gdh = 0.0_dp
      frl = 0.0_dp
      fbb = 0.0_dp
      fbp = 0.0_dp
      fbt = 0.0_dp

      do ia = mpmap(iproc)+1, mpmap(iproc+1)

         za = z(ia)
         call states(za,nla,nstta,llista)

         if (momflg == 1) then
            lla0 = 0
            do la = 1,nla
               nma = 2*llista(la) + 1
               mlista(lla0+1:lla0+nma) = la
               lla0 = lla0 + nma
            enddo
         else
            do la = 1, nstta
               mlista(la) = la
            end do
         endif

         call assoc_ham(ia)
         pcl(ia) = pcluster(2) - pcluster(1) ! avoid onsite by starting at pcuster(1)
         cls(1:pcl(ia),ia) = cluster(pcluster(1):pcluster(2)-1)
         jd = pcluster(1) - 1
         fatl = 0.0_dp

      !    For each neighbor ....
         do jb = 1, pcl(ia)
            ib = cls(jb, ia)

            mib = map(ib)
            zb = z(ib)
            call states(zb,nlb,nsttb,llistb)

            bt = bsort(btype(za,zb))
            if (bt == -1) bt = bsort(btype(zb,za))
            pwp => bnd(bt) % pwp

            dr = ad(1:3,ia) - ad(1:3,ib)



!     it may be possible to avoid computing grad 2 times in the magnetic case
            call grdmat(grad,dr,za,nla,nstta,llista,zb,nlb,nsttb,llistb,scfcut,dsf)

            if (momflg == 2) call utranv(transm(1,1,ia),grad,transm(1,1,mib),trwork,mxnstat,nstta,nsttb)

            gdh(1:3, 1:nstta, 1:nsttb, jb, ia) = grad(1:3, 1:nstta, 1:nsttb)

            frll = 0.0_dp
            fab = 0.0_dp


            do isp = 1, nsp
               call assoc_ab(ia,isp)
               call assoc_chi(ia,isp)

               rhol = 0.0_dp

               do la = 1,nstta
                  ma = mlista(la)
                  nmax = lchain(ma)
                  do lb = 1, nsttb
! may be done with dot_product instead sum, not sure which is faster
                     rhol(la,lb) =  sum(chia(0:nmax,ma) * darec(0:nmax,la,lb,jb+jd)) &
                              & + 2*sum(chib(1:nmax,ma) * dbrec(1:nmax,la,lb,jb+jd))
                     if (term == 1) rhol(la,lb) = rhol(la,lb) + 2*chib(nmax+1,ma)*dbrec(nmax+1,la,lb,jb+jd)
                     fabl = grad(1:3, la, lb) * rhol(la, lb)
                     frll(1:3, la, lb) = frll(1:3, la, lb) +  fabl
                     fab = fab + fabl
                  end do
               end do
               rho(1:nstta, 1:nsttb, isp, jb, ia) = rhol(1:nstta, 1:nsttb)
            end do ! isp

            if (.not. mag) then
               fab = 2*fab
               frll = 2*frll
            end if

            do lb = 1, nsttb
               do la = 1, nstta
                  frll(4, la, lb) = sqrt(sum(frll(1:3, la, lb)*frll(1:3, la, lb)))
               end do
            end do

            frl(1:4, 1:nstta, 1:nsttb, jb, ia) = frll

            fbb(1:3, jb, ia) = 2*fab
            fbb(4, jb, ia) = sqrt(sum(fbb(1:3, jb, ia)*fbb(1:3, jb, ia)))

            r = sqrt(sum(dr*dr))
            call tailed_fun(r, pwp % fp, pwp % tl, 1, fv)


            fbp(1:3, jb, ia) = -fv(1)*dr/r
            fbp(4, jb, ia) = abs(fv(1))

            fbt(1:3, jb, ia) = fbp(1:3, jb, ia) + fbb(1:3, jb, ia)
!             fatl = fatl +  fbb (1:3, jb, ia) !!!! only the band forces will be checked !+ fbt(1:3, jb, ia)
            fbt(4, jb, ia) = sqrt(sum(fbt(1:3, jb, ia)*fbt(1:3, jb, ia)))

         enddo ! ib

!          fat(1:3, ia) = fatl
!          fat(4, ia) = sqrt(sum(fatl*fatl))
      end do ! ia

!       the sharing and mapping shalt happen here and not later when the map is rewritten

      if (nproc > 1) then
          call system_clock(c1_m)

          allocate(mpcnt(0:nproc-1), szs(0:nproc-1), displs(0:nproc))

          mpcnt = mpmap(1:nproc)-mpmap(0:nproc-1)

          usz = 1
          szs = mpcnt*usz
          displs = mpmap*usz

    !       print '(a,i0)', 'mpi sharing from :', iproc


          call mpi_gatherv(pcl(mpmap(iproc)+1), szs(iproc), mpi_integer, &
                         & pcl, szs, displs, mpi_integer, master, mpi_comm_world, ierr)
    !       print '(a,i0)', 'pcl shared from :', iproc


          usz = mxnnb
          szs = mpcnt*usz
          displs = mpmap*usz


          call mpi_gatherv(cls(1, mpmap(iproc)+1), szs(iproc), mpi_integer, &
                         & cls, szs, displs, mpi_integer, master, mpi_comm_world, ierr)
    !       print '(a,i0)', 'cls shared from :', iproc


          usz = mxnstat*mxnstat*nsp*mxnnb
          szs = mpcnt*usz
          displs = mpmap*usz

          call mpi_gatherv(rho(1, 1, 1, 1, mpmap(iproc)+1), szs(iproc), mpi_real8, &
                         & rho, szs, displs, mpi_real8, master, mpi_comm_world, ierr)

    !       print '(a,i0)', 'rho shared from :', iproc

          usz = 3*mxnstat*mxnstat*mxnnb
          szs = mpcnt*usz
          displs = mpmap*usz

          call mpi_gatherv(gdh(1, 1, 1, 1, mpmap(iproc)+1), szs(iproc), mpi_real8, &
                         & gdh, szs, displs, mpi_real8, master, mpi_comm_world, ierr)

    !       print '(a,i0)', 'gdh shared from :', iproc

          usz = 4*mxnstat*mxnstat*mxnnb
          szs = mpcnt*usz
          displs = mpmap*usz

          call mpi_gatherv(frl(1, 1, 1, 1, mpmap(iproc)+1), szs(iproc), mpi_real8, &
                         & frl, szs, displs, mpi_real8, master, mpi_comm_world, ierr)
    !       print '(a,i0)', 'frl shared from :', iproc

          usz = 4*mxnnb
          szs = mpcnt*usz
          displs = mpmap*usz

          call mpi_gatherv(fbb(1, 1, mpmap(iproc)+1), szs(iproc), mpi_real8, &
                         & fbb, szs, displs, mpi_real8, master, mpi_comm_world, ierr)

          call mpi_gatherv(fbp(1, 1, mpmap(iproc)+1), szs(iproc), mpi_real8, &
                         & fbp, szs, displs, mpi_real8, master, mpi_comm_world, ierr)

          call mpi_gatherv(fbt(1, 1, mpmap(iproc)+1), szs(iproc), mpi_real8, &
                         & fbt, szs, displs, mpi_real8, master, mpi_comm_world, ierr)

    !       print '(a,i0)', 'f** shared from :', iproc

    !       usz = 4
    !       szs = mpcnt*usz
    !       displs = mpmap*usz
    !
    !       call mpi_gatherv(fat(1, mpmap(iproc)+1), szs(iproc), mpi_real8, &
    !                      & fat, szs, displs, mpi_real8, master, mpi_comm_world, ierr)

          deallocate(mpcnt, szs, displs)

          call system_clock(c2_m)
          if (iproc == master) print *, 't(fdnet@root):', real(c2_m-c1_m,dp)/real(cr,dp),'s'
      end if

      if (iproc == master) then

         allocate(fat( 4, 1:nd))
         fat = 0.0_dp

         do ia = 1, nd
            do jb = 1, pcl(ia) ! avoid onsite by starting at pcuster(1)
               ib = cls(jb,ia)
               mib = map(ib)
               if (mib /= 0) then
                  fat(1:3,  ia) = fat(1:3,  ia) + 0.5_dp*fbb(1:3, jb, ia)
                  fat(1:3, mib) = fat(1:3, mib) - 0.5_dp*fbb(1:3, jb, ia)
               else
                  fat(1:3,  ia) = fat(1:3,  ia) + fbb(1:3, jb, ia)
               endif
            end do
         end do

         ferr = sum(abs(fbs(1:3, 1:nd) - fat(1:3, 1:nd)))

         if (ferr > 0.0001_dp) then
            print *, redb//'ferr:', ferr,endc
         else
            print *, 'ferr:', ferr
         endif

!          do ia = 1, nd
!             print '(2(x,3(x,f10.5)))', fbs(1:3, ia), fat(1:3, ia)
!          end do

         deallocate (fat)

         na = 0
         nb = 0

         do ia = 1, nd
            if (sqrt(sum((ad(1:2,ia) - disl_core)**2)) > 2*rcut) cycle
            na = na + 1
            do jb = 1, pcl(ia)
               if (fbt(4,jb,ia) == 0.0_dp) cycle
               nb = nb + 1
            end do
         end do


         write(filename,'("fdetails/fd.st.rlx.",f5.3,".dat")') totstr

         u = get_new_unit()
         open(u, file=trim(filename), action='write')

         write(u, "('na: ', i0, '; nb: ', i0)") na, nb

         do ia = 1, nd
            dda = sqrt(sum((ad(1:2,ia) - disl_core)**2))
            if (dda > 3*rcut) cycle

            za = z(ia)
            call states(za,nla,nstta,llista)

            do jb = 1, pcl(ia)
               if (fbt(4,jb,ia) == 0.0_dp) cycle

               ib = cls(jb,ia)
               zb = z(ib)
               call states(zb,nlb,nsttb,llistb)

               ddb = sqrt(sum((ad(1:2,ib) - disl_core)**2))
               dr = ad(1:3,ib) - ad(1:3,ia)

               write(u, "('b: ',i5,'-',i5,';  z: ',i3,x,i3,';  nL: ',i2,x,i2,';  cd: ',f6.3,x,f6.3)") &
                                                               & ia, ib, za, zb, nstta, nsttb, dda, ddb
               write(u, "(2x,'r:',2(x,3(x,f8.3),x,a))") ad(1:3,ia), '|', ad(1:3,ib), ''
               write(u, "(2x,'dr:',x,4(x,f8.3))") sqrt(sum(dr*dr)), dr
               write(u, "(2x,'ft:',x,4(x,f8.3))") fbt(4,jb,ia), fbt(1:3,jb,ia)
               write(u, "(2x,'fb:',x,4(x,f8.3))") fbb(4,jb,ia), fbb(1:3,jb,ia)
               write(u, "(2x,'fp:',x,4(x,f8.3))") fbp(4,jb,ia), fbp(1:3,jb,ia)

               write(u, "(2x,'frl: |f| x y z')")

               do la = 1, nstta
                  write(u, "(4x)", advance='no')
                  do lb = 1, nsttb
                     write(u, "(x,f8.3)", advance='no') frl(4,la,lb,jb,ia)
                  end do
                  write(u, "(x,'|',x)", advance='no')
                  do lb = 1, nsttb
                     write(u, "(x,f8.3)", advance='no') frl(1,la,lb,jb,ia)
                  end do
                  write(u, "(x,'|',x)", advance='no')
                  do lb = 1, nsttb
                     write(u, "(x,f8.3)", advance='no') frl(2,la,lb,jb,ia)
                  end do
                  write(u, "(x,'|',x)", advance='no')
                  do lb = 1, nsttb
                     write(u, "(x,f8.3)", advance='no') frl(3,la,lb,jb,ia)
                  end do
                  write(u,'(x)')
               end do

               write(u, "(2x,'rho: nsp: ', i1)") nsp

               do la = 1, nstta
                  write(u, "(4x)", advance='no')
                  do isp = 1, nsp
                     do lb = 1, nsttb
                        write(u, "(x,f8.3)", advance='no') rho(la,lb,isp,jb,ia)
                     end do
                     if (isp /= 2 .and. nsp == 2) write(u, "(x,'|',x)", advance='no')
                  end do
                  write(u,'(x)')
               end do

               write(u, "(2x,'grad H: x y z')")

               do la = 1, nstta
                  write(u, "(4x)", advance='no')
                  do lb = 1, nsttb
                     write(u, "(x,f8.3)", advance='no') gdh(1,la,lb,jb,ia)
                  end do
                  write(u, "(x,'|',x)", advance='no')
                  do lb = 1, nsttb
                     write(u, "(x,f8.3)", advance='no') gdh(2,la,lb,jb,ia)
                  end do
                  write(u, "(x,'|',x)", advance='no')
                  do lb = 1, nsttb
                     write(u, "(x,f8.3)", advance='no') gdh(3,la,lb,jb,ia)
                  end do
                  write(u,'(x)')
               end do
               write(u, '(/)')
            end do
            write(u, '(/)')
         end do


         call close_unit(u)


      end if

      deallocate(rho, gdh, frl, fbb, fbp, fbt, pcl, cls)

      call system_clock(c2)

      if (iproc == master) print *, 't(fd):', real(c2-c1,dp)/real(cr,dp), 's'

   end subroutine print_forcedetails



!    subroutine print_forcedetails()
!
!       include "Include/ag.conn"
!       include 'Include/Force.array'
!
!
!             if (mib /= 0) then
! ! Symmetrisation
!                fat(1:3,  ia) = fat(1:3,  ia) + fab
!                fat(1:3, mib) = fat(1:3, mib) - fab
!             else
!                fat(1:3,  ia) = fat(1:3,  ia) + fab*2
!             endif
!
!    end subroutine print_forcedetails


    subroutine init()

      include "Include/ag.conn"

      character(len=80) :: info
      complex(dp) :: akd(3,3), p(3)
      real(dp) :: ang
      real(dp) :: bx0, by0, bz0
      real(dp) :: c11, c22, c12, c13, c33, c66, c44, c55, c45
      integer :: ndis, i, j, u

      if (rcutl /= rcut) stop 'forcedetails: only valid for rcutl == rcut'

      u = get_new_unit()
      open(u, file='fort.21', action='read')

      read (u, '(i1)' ) ndis

      if (ndis /= 1) stop 'handling more than 1 dislocations is not implemented in fdetails yet.'
      
      read ( u, '(1x,a80)') info
      read ( u, *) c11, c22, c12, c13, c33, c66, c44, c55, c45

      do j = 1, ndis
         do i = 1, 2
            read ( u, '(1X,A80)') info
         end do

         read ( u, '(1X,A80)') info
         read ( u, * )  akd

         read ( u, '(1x,a80)') info
         read ( u, * )  p

         read ( u, '(1x,a80)') info
         read ( u, * ) disl_core, ang

         read ( u, '(1x,a80)') info
         read ( u,*) bx0, by0, bz0
      end do

      call close_unit(u)

      inited = .true.

    end subroutine init

   end module mod_forcedetails
