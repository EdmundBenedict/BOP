

   subroutine bldnebt(fs, rlim, lena, a, rprune, rcut, stressed, et, &
                  & nd, ad, z, ninert, adinert, zinert, &
                  & map, aptr, bptr, maxnumnb, totnd, nne, mapi,maxnb)
      use mod_precision
      use mod_all_scalar, only : quiet, mag, mvflag
      use mod_const
      use mod_clock, only : cr

      use topologia, only : sdistrib, nproc, iproc, master ! truly no more than these shall be taken cause the rest will be local

      use mpi

!    This is a routine to build the neighbor tables for atoms.
!
! fs needs proper map() while the rest of the code uses map=0 for inert atoms.
! fs also sets the self atom at first position rather than wherever it comes.

      implicit none

      logical, intent(in) :: fs, stressed

!*** The number of unit cells to include in each dimension.
      integer, intent(in) ::  nd, rlim(3), ninert,maxnb

      real(dp), intent(in) :: a(3, 3), rprune, rcut, et(3, 3)
      real(dp), intent(inout) :: lena(3), ad(3, mxtotnd), adinert(3, minert)

!*** MAP is an array linking the central cell atoms to their periodically
!    repeated neighbors.
      integer, intent(inout) :: map(mxtotnd), z(mxtotnd), zinert(minert), mapi(mxtotnd)

!*** APTR is an array of pointers to the neighbor list.
!*** BPTR is the neighbor list.
      integer, intent(out) :: aptr(mxtotnd), bptr(mxnbl), maxnumnb, totnd, nne




!       include "Include/PosVel.array"
!       include "Include/ag.conn"





      real(dp) :: rcut2

      integer inda,indb,totnbox
      integer ia1,ia2,ia3,ib1,ib2,ib3,i
      integer ia,ib,ic
      integer ja,jb,p
      integer nb
      integer ns
      integer :: maxnnbia
!     INTEGER SUM


!*** Size of periodically repeated unit cell.
      real(dp) :: span(3)
!*** Origin for total cell of atoms.
      real(dp) :: org(3)
!*** Size of one of partitions.
      real(dp) :: dl(3)
!*** DLIM sets the limits on the number of neighboring cells to be scanned.
      real(dp) :: dlim(2,3)
!*** Extremities of atomic coordinates.
      real(dp) :: rmax(3)
      real(dp) :: rmin(3)

!*** NNA is an array of the number of neighbors for each atom
!     INTEGER NNA(MXTOTND)
!*** The number of boxes in a given direction used to partition space.
      integer nbox(3)


!*** The number of atoms in a box.
      integer boxcount(mbox)

!*** The list of atoms for each box.
      integer box(mcount,mbox)

!*** The type of neighbour list builder to be used. 1 => O(N^2)
!***                                                2 => O(N)
      integer :: nbtflg


      real(dp) :: a1a2,a1a3,a2a3

      logical :: lmapi

      real(dp) :: strad(3)


      integer(8) :: c1, c2, lc1, lc2, nbc, nwc


      integer :: aproc, maxnn(2), ierr
      integer, allocatable :: mpmap(:), nbs(:), lbuf(:), szs(:)

      call system_clock(c1)

      lmapi = (.not. fs) .and. mag .and. (ninert > 0)

      maxnnbia = 0

      a1a2 = abs(a(1,1)*a(2,1)+a(1,2)*a(2,2)+a(1,3)*a(2,3))
      a1a3 = abs(a(1,1)*a(3,1)+a(1,2)*a(3,2)+a(1,3)*a(3,3))
      a2a3 = abs(a(2,1)*a(3,1)+a(2,2)*a(3,2)+a(2,3)*a(3,3))

      if ((a1a2 < 1.0e-8_dp) .and. (a1a3 < 1.0e-8_dp) .and. (a2a3 < 1.0e-8_dp) .and. (.not. stressed)) then
         nbtflg = 2
      else
         nbtflg = 1
      endif

      lena = sqrt(sum(a*a, dim=1))

      dlim(1,:) = -rprune/lena
      dlim(2,:) = 1.0_dp + rprune/lena

      do ia = 1,nd
         map(ia) = ia
      enddo

      nne = nd + ninert

      if (ninert > 0) then
         if (fs) then
            do ia = nd+1,nne
               map(ia) = ia
            end do
         else
            map(nd+1:nne) = 0
            if (lmapi) then ! or just check mag
               do ia = 1,ninert
                  mapi(nd+ia) = ia
               end do
            end if
         end if
         z(nd+1:nne) = zinert(:ninert)
         ad(1:3,nd+1:nne) = adinert(1:3,1:ninert)
      end if

      if (mvflag /= 10 .and. any(rlim > 0)) call pbc(nne, a, rlim, ad)

      if (ninert > 0) adinert(1:3, 1:ninert) = ad(1:3, nd+1:nne)



      call bldlist(nne, rlim, a, dlim, ad, totnd, map, z, lmapi, nd, mapi)

      if (stressed) then
         do  i = 1, totnd
            strad = et(1:3, 1) * ad(1, i) + et(1:3, 2) * ad(2, i) + et(1:3, 3) * ad(3, i)
            ad(1:3, i) = ad(1:3, i) + strad
         enddo
      end if

      ns = totnd

!      WRITE(6,'("TOTND = ",I5)') TOTND

      maxnumnb = 0

      rcut2 = rcut*rcut

      if (nbtflg == 1) then
!    Build the neighbor list - this is simple O(N^2) method.

! This is parallelised in the hope that mpi collectives will be faster that the actual neighbour finding.
! Even if the above condition is true the idea is still not great. A much better approach will be to make/take O(N)
! neighour finder for nonorthogonal cell.


         call system_clock(lc1)
         if (.not. quiet) write(6,'(/''Using simple neighbor table builder.'')')
         nb = 1

         allocate(mpmap(0:nproc))
         call sdistrib(ns, nproc, mpmap, aproc)
         if (.not. quiet) print *, 'nbmpmap:', mpmap

         do ia = mpmap(iproc)+1, mpmap(iproc+1) ! do not change to 1, nsia because the real ia is needed
!          do ia = 1, ns
            aptr(ia) = nb

            if (.not. fs) then ! the pair potential bptr contains the self atom in first position while the tb list does not.
               bptr(nb) = ia
               nb = nb + 1
            end if

            do ib = 1, ns
               if (sum((ad(1:3,ia) - ad(1:3,ib))**2) <= rcut2 .and. (ib /= ia)) then
                  bptr(nb) = ib
                  nb = nb + 1
               endif
            enddo


            if (nb-aptr(ia) > maxnb) then
               write(6,'(''Number of neighbors for atom # '',i5,'' = '',i3)') ia, nb - aptr(ia)
               write(6,'(''This is more than the maximum allowed of '', i3)') maxnb
               write(6,'(''=> Increase MXNNB.'')')
               call panic()
            endif
!           NNA(IA) = NB-APTR(IA)-1
            if ( maxnumnb  <  ( nb-aptr(ia)-1 ) ) then
               maxnumnb = nb-aptr(ia)-1
               maxnnbia = ia
            end if

            bptr(nb) = eol
            nb = nb + 1
         enddo

         nb = nb - 1

         call system_clock(lc2)
         nbc = lc2 - lc1
         if (.not. quiet) print *, 't(bldn@root):', real(nbc,dp)/real(cr,dp)

         if (nproc > 1) then
             call system_clock(lc1)

             allocate(nbs(0:nproc-1), szs(0:nproc-1))
             szs = mpmap(1:nproc) - mpmap(0:nproc-1)
    !          print *,'szs',szs
             allocate(lbuf(szs(iproc)))

             call mpi_allgather(nb, 1, mpi_integer, nbs, 1, mpi_integer, mpi_comm_world, ierr)
    !          print *, 'nbs:', iproc, ':', nbs

             if (iproc == 0) then  ! the condition is to avoid illegal index in nbs
                lbuf = aptr(mpmap(iproc)+1 : mpmap(iproc+1))
             else
                lbuf = aptr(mpmap(iproc)+1 : mpmap(iproc+1)) + sum(nbs(0:iproc-1))
             end if

             call mpi_allgatherv(lbuf, szs(iproc), mpi_integer, aptr, szs, mpmap, mpi_integer, mpi_comm_world, ierr)

             deallocate(lbuf)

             allocate(lbuf(nb))

             lbuf = bptr(:nb)

    !          Size is now reused as displacements/offsets array to bptr as mpmap is used for aptr
             szs(0) = 0
             do i = 1, nproc-1
                szs(i) = szs(i-1) + nbs(i-1)
             end do

             call mpi_allgatherv(lbuf, nb, mpi_integer, bptr, nbs, szs, mpi_integer, mpi_comm_world, ierr)

             nb = sum(nbs)

             deallocate(szs, lbuf, nbs, mpmap)

             maxnn = [maxnumnb, maxnnbia]
             call mpi_allreduce(maxnn, mpi_in_place, 1, mpi_2integer, mpi_maxloc, mpi_comm_world, ierr)
             maxnumnb = maxnn(1)
             maxnnbia = maxnn(2)

    !          call mpi_barrier(mpi_comm_world, ierr)

             call system_clock(lc2)
             nwc = lc2 - lc1
             if (.not. quiet) then
                print *, 't(netw@root):', real(nwc,dp)/real(cr,dp)
                if (((nproc-1)*nbc) <= nwc) then
                   print *, 'You might be better off using the serial neighbour finder.'
                   print *, 'Try adding #undef MPI in the beginning of the first line of bldnebt.F90'
                endif
             end if

    !          call mpi_finalize(ierr)
    !          stop
        end if

        if (.not.quiet) write(6,'(''Total number of atoms = '', i6)') ns
        if (.not.quiet) write(6,'(''Size of neighbor list = '', i6)') nb
!        WRITE(10,'(15I5)') (BPTR(IA),IA=1,NB,1)



      elseif (nbtflg == 2) then
!    Build the neighbor list - this is the linked list O(N) method.
!     Note: this requires a tetrahedral unit cell.

         if (any(rlim == 0)) then
            rmin = min(minval(ad, dim=2), minval(adinert, dim=2))
            rmax = max(maxval(ad, dim=2), maxval(adinert, dim=2))
         endif

         do i = 1, 3
            if (rlim(i) == 0) then
               span(i) = rmax(i)-rmin(i)+0.001_dp
               org(i) = rmin(i)
            else
               span(i) = lena(i) + 2*rprune
               org(i) = dlim(1,1)*a(1,i)+dlim(1,2)*a(2,i)+dlim(1,3)*a(3,i)
            endif
         end do

         nbox = max(1,int(span/rcut))


!          write(6,'(''SPAN = '',3G13.5)') span
!          write(6,'(''ORG = '',3G13.5)') org
!          write(6,'(''NBOX = '',3I5)') nbox

         totnbox = product(nbox)
         if (totnbox > mbox) then
            write(6,'(''Too many boxes for neighbor list calculation.'')')
            write(6,'(''Increase MBOX to at least '',I5)') totnbox
            call panic()
         endif

         dl(:) = span(:)/real(nbox(:), dp)
!          write(6,'(''DL = '',3G13.5)') dl

         boxcount(1:totnbox) = 0

         do ia = 1,ns

            ia1 = int((ad(1,ia)-org(1))/dl(1))
            ia2 = int((ad(2,ia)-org(2))/dl(2))
            ia3 = int((ad(3,ia)-org(3))/dl(3))

            inda = 1 + ia1 + nbox(1)*(ia2 + nbox(2)*ia3)

!             if ((inda < 1).or.(inda > totnbox)) then
!                write(6,'(''Attempting to store atom in illegal box.'')')
!                write(6,'(''INDA = '',I5)') inda
!                call panic()
!             endif
!

            boxcount(inda) = boxcount(inda) + 1
            if (boxcount(inda) > mcount) then
               write(6,'(''Too many atoms in one box.'')')
               write(6,'(''Increase MCOUNT.'')')
               call panic()
            endif

            box(boxcount(inda),inda) = ia

         enddo

!        WRITE(6,'(''Boxcount:'')')
!        WRITE(6,'(10I5)') (BOXCOUNT(INDA),INDA=1,
!     +                      NBOX(1)*NBOX(2)*NBOX(3),1)

         nb = 1
         do ia = 1,ns
            aptr(ia) = nb

            if (.not. fs) then  ! the pair potential bptr contains the self atom in first position while the tb list does not.
               bptr(nb) = ia
               nb = nb + 1
            end if

            ia1 = int((ad(1,ia)-org(1))/dl(1))
            ia2 = int((ad(2,ia)-org(2))/dl(2))
            ia3 = int((ad(3,ia)-org(3))/dl(3))

            do ib1 = max(0,ia1-1),min(nbox(1)-1,ia1+1)
               do ib2 = max(0,ia2-1),min(nbox(2)-1,ia2+1)
                  do ib3 = max(0,ia3-1),min(nbox(3)-1,ia3+1)
                     indb = 1 + ib1 + nbox(1)*(ib2 + nbox(2)*ib3)
                     do jb = 1,boxcount(indb)
                        ib = box(jb,indb)
                        if (sum((ad(:,ia) - ad(:,ib))**2) <= rcut2 .and. (ia /= ib)) then
!                            deleted the ordering by ib, if needed better use a better algorithm.
                           bptr(nb) = ib
                           nb = nb + 1
                        endif
                     enddo
                  enddo
               enddo
            enddo

            if (nb-aptr(ia) > maxnb) then
               write(6,'(''Number of neighbors for atom # '',I5,'' = '',I3)') ia,nb-aptr(ia)
               write(6,'(''This is more than the maximum allowed of '',I3)') maxnb
               write(6,'(''=> Increase MXNNB.'')')
               call panic()
            endif

!           NNA(IA) = NB-APTR(IA)-1

            if ( maxnumnb  <  ( nb-aptr(ia)-1 ) )   then
                maxnumnb = nb-aptr(ia)-1
                maxnnbia = ia
            endif

            bptr(nb) = eol
            nb = nb + 1

         end do



!        WRITE(6,'(''Size of neighbor list = '',I5)') NB
!        WRITE(11,'(15I5)') (APTR(IA),IA=1,TOTND,1)
!        WRITE(11,'(15I5)') (BPTR(IA),IA=1,NB,1)

      endif

      aptr(ns+1) = nb



      call system_clock(c2)

      if (.not. quiet) print '("maxnumnb ",i0," of atom ",i0)',  maxnumnb, maxnnbia
      if (.not. quiet) print *, 't(bldnebt):', real(c2-c1,dp)/real(cr,dp)


      !
!       do ia = 1, nd
!          print *, 'ia:', ia
!          p = 0
!          do while (bptr(aptr(ia)+p)/=eol)
!             print *, '  ib:', p+1, bptr(aptr(ia)+p)
!             p = p + 1
!          end do
!       end do

      end subroutine bldnebt


