   subroutine getebsfbs(flag)
      use mod_precision
      use mod_all_scalar
      use mod_const
      use topologia
      use mod_ham
      use ab_io
      use mod_chi
      use mod_clock
      use mod_conf, only : rtc
      use mod_atom_ar
      use mpi
      use mod_par

      implicit none

      include "Include/Atom.array"
      include "Include/ag.conn"

      integer, intent(inout) :: flag

      integer :: it, mxit, ia,  inflag, ierr, mit, mgit
      real(dp) :: etott, merr

      character(len=30) :: ffmt, ffmtp

      procedure(real(dp)) :: getefsrt, getefnull, numel, eff

      integer :: loc_clusiz
      integer(8) :: c1,c2

      integer :: i
      real(dp) :: ef,nf,sumz,sumzde

      real(dp), allocatable :: mg_o(:), mg_nd(:)

      type(trans_pack_t), allocatable :: tpack(:)

      write(ffmt,"(a,i0,a)") '(a,":",', nd, '(x,f22.14))'
      write(ffmtp,"(a,i0,a)") '(a,x,i0,x,":",', nd, '(x,f22.14))'

      sumz = 0.0_dp
      sumzde = 0.0_dp

      loc_clusiz = 0

      call system_clock(c1)

      do ia = mpmap(iproc)+1, mpmap(iproc+1)
         call assoc_ab(ia,1)
         call assoc_ham(ia)
!          print *, 'ia:', ia
      !          call system_clock(c1)
         call getnch(ia)   !Find number of linear chains per site.
         call bldclus(ia) !Build cluster.
         loc_clusiz = loc_clusiz + pcluster(nbase+1)-1
         call bldh() !Build Hamiltonian.
      !          call system_clock(c2)
   !          print *,'hamrel',real(c2-c1,dp)/real(cr,dp)
      end do


      call system_clock(c2)

      if (.not. quiet) print *, 't(bldh@root):', real(c2-c1,dp)/real(cr,dp)

      if (nproc > 1) call mpi_allreduce(mpi_in_place, loc_clusiz, 1, mpi_integer, mpi_sum, mpi_comm_world,ierr)

      aveclusiz = loc_clusiz/real(nd,kind=dp)

      if (.not. quiet) print *, 'aveclusiz:', aveclusiz

      inflag = flag
      mxit = rtc % bop_conf % mxit
      if (mag) then
        mgit = rtc % bop_conf % mgit
        merr = rtc % bop_conf % merr
        mit = 1
        if (merr < 0.0_dp) mgit = 1
      end if


      sumz = sum( (/ (zc(z(ia)), ia = mpmap(iproc)+1,mpmap(iproc+1)) /) )
      if (nproc > 1) call mpi_allreduce(mpi_in_place, sumz, 1, mpi_real8, mpi_sum, mpi_comm_world,ierr)

      if (qerr < 0.0_dp) mxit = 1
      if (.not. quiet) write(6,"('Imposing LCN ...')")

      if (mag) allocate(mg_o(mpmap(iproc)+1:mpmap(iproc+1)))

      do
         if (mag) then

            if (.not. quiet) write(6,'(/,"mit:",i3)') mit

            mg_o = mg(mpmap(iproc)+1:mpmap(iproc+1))
            if (nproc > 1) call prep_vlvec_ir_allgather(tpack,iproc,aproc,master,mpmap,ierr)

            do ia = mpmap(iproc)+1, mpmap(iproc+1)
               dem(ia) =  0.5_dp*istn(z(ia))*mg(ia)
            end do
            if (nproc > 1) call exec_vlvec_ir_allgather(tpack,iproc,aproc,master,mpmap,nd,dem,ierr)
         end if
         it = 1
         do
            if (.not. quiet) write(6,'(/,"it:",i3)') it
            if (nproc > 1) call prep_vlvec_ir_allgather(tpack,iproc,aproc,master,mpmap,ierr)
            sumzde = sum((/ (de(ia)*zc(z(ia)), ia = mpmap(iproc)+1,mpmap(iproc+1)) /))
            if (nproc > 1) call mpi_allreduce(mpi_in_place, sumzde, 1, mpi_real8, mpi_sum, mpi_comm_world,ierr)
            de(mpmap(iproc)+1:mpmap(iproc+1)) = de(mpmap(iproc)+1:mpmap(iproc+1)) - sumzde/sumz
            
            if (.not. quiet) then
                  print *, 'de',de(:nd)
                    print *, 'mg',mg(:nd)
            endif
            
            if (nproc > 1) call exec_vlvec_ir_allgather(tpack,iproc,aproc,master,mpmap,nd,de,ierr)
            call spnbop(2) ! set to 1 to avoid energy calculation
            if (nproc > 1) call mpi_allreduce(mpi_in_place, dqmax, 1, mpi_real8, mpi_max, mpi_comm_world,ierr)
            if (dqmax <= qerr .or. it == mxit) exit
            de(mpmap(iproc)+1:mpmap(iproc+1)) = de(mpmap(iproc)+1:mpmap(iproc+1)) + dq2chia(mpmap(iproc)+1:mpmap(iproc+1))
            it = it + 1
         end do
         if (.not. mag) exit
         mgmax = maxval(abs(mg(mpmap(iproc)+1:mpmap(iproc+1)) - mg_o))
         
         if (nproc > 1) call mpi_allreduce(mpi_in_place, mgmax, 1, mpi_real8, mpi_max, mpi_comm_world,ierr)
         if (iproc == master) print *, 'mgmax:', mgmax
         
         if (mgmax <= merr .or. mit == mgit) exit
         mit = mit + 1
      end do

      call system_clock(c1)
      call eval_bsens(lef)
      if (forces) call forcebop()
      call system_clock(c2)
      if (.not. quiet) print *, 't(enforce@root):', real(c2-c1,dp)/real(cr,dp)

      if (mag) then
         if (nproc > 1) call prep_vlvec_ir_allgather(tpack,iproc,aproc,master,mpmap,ierr)
!          mg(mpmap(iproc)+1:mpmap(iproc+1)) = mg_o(mpmap(iproc)+1:mpmap(iproc+1))
         if (nproc > 1) call exec_vlvec_ir_allgather(tpack,iproc,aproc,master,mpmap,nd,mg,ierr)
!          if (.not. quiet) write(6, ffmt) 'mg', mg
      end if


      if (mag) deallocate(mg_o)

!     may be write de and mg to the atoms in the cell object


   end subroutine getebsfbs

