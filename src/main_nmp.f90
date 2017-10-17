! 
   program main

!       An atom in the parallel version should ideally broadcast only to
!       the atoms in its cluster up to the chosen recursion level.
!       This is not what is happening yet!!! To be fixed when very many
!       processors are used on large systems

      use mod_const
      use mod_conf
      use mod_io
      use mod_clock
      use topologia, only : iproc, nproc, master, mpmap
      use mod_misc, only : get_urandom, get_stamp, in_argv
!       use ssa, only : ssiman
      use mpi

      implicit none

      type(run_conf_t) :: rtconf
      type(cell_t), target :: ucell
      type(ifl_t) :: conf_fl, cell_fl
      character(len=1), allocatable :: tmp_str(:)
      integer :: flsz

      character(len=200) :: cmd
      integer :: cmd_len

      integer :: namelen, version, subversion, ierror
      character(len=100) :: fname
      character(len=50) :: stamp


      character(len=100) :: processor_name
      real(8) :: start_time, end_time

      integer :: seed_size, seed_io
      integer, allocatable :: seed(:)

      integer :: u,conf_io


      call mpi_init(ierror)
      call init_messy_consts()
      call get_clock_prec()
      call init_io()

      call mpi_comm_rank(mpi_comm_world, iproc,ierror);
      call mpi_comm_size(mpi_comm_world, nproc, ierror);
      if (nproc > 1) then
          call mpi_get_processor_name(processor_name, namelen, ierror);
          call mpi_get_version(version, subversion,ierror);
          write(*,'(a,i0,a,i0,a,a,a,i0,a,i0)') &
             & "process ", iproc, " out of ", nproc, &
             & ", name: ", trim(processor_name), &
             & ", mpi version: ", version,".",subversion
      end if
      master = 0
      allocate(mpmap(0:nproc))




   !Distribute the seed so every process has the same.
      if ( iproc == master ) call random_seed(size=seed_size)
      if (nproc > 1) call mpi_bcast(seed_size, 1, mpi_integer, master,  mpi_comm_world, ierror )
      allocate(seed(seed_size))
      if (iproc == master) call get_urandom(seed)
      if (nproc > 1) call mpi_bcast(seed, seed_size, mpi_integer, master,  mpi_comm_world, ierror )
      call random_seed(put=seed)


! could be useful later
!       if (iproc == master) then
!          seed_io = get_new_unit()
!          open(unit=seed_io, file='seed.out', action='write')
!          call print_c(seed,seed_io,'seed:')
!          call close_unit(seed_io)
!       end if

      if (iproc == master) then
         call fl_to_s('conf.in',tmp_str)
         flsz = size(tmp_str)
      end if

      if (nproc > 1) call mpi_bcast(flsz,1, mpi_integer, master, mpi_comm_world, ierror)
      if (iproc /= master) allocate(tmp_str(flsz))
      if (nproc > 1) call mpi_bcast(tmp_str, flsz, mpi_character, master, mpi_comm_world, ierror)
      call set_s(conf_fl,tmp_str)
      deallocate(tmp_str)

      call read_run_conf(rtconf, conf_fl)

      if (in_argv('-fsa',1)) then
         stop 'remove this stop when the fit mods are brought up to the changes in mod_conf and uncomment "call ssiman" below'
         if (iproc == master) then
               call get_stamp(stamp)
               print *, 'stamp: ', stamp
               u = get_new_unit()
               write(fname, '("out-",a,".out")') trim(stamp)
               !             write(fname, '("out-p",i0,"-.out")') iproc
               open(u, file=trim(fname),action='write')
               start_time = mpi_wtime()
         end if
         rtconf % quiet = .true.
!          call ssiman(rtconf,u)
         if (iproc == master) call close_unit(u)
      else
         if (iproc == master) then


!    This is a program to evaluate the Bond Order Potentials using the
!     square root terminator.

!     This call starts up either the scalar or parallel version depending
!     on which library has been linked in.
!
!        open(unit = 8, file = 'conf.in', action = 'read' )
!        open(unit = 31, file = 'ucell.block', action = 'read' )
!
            open( unit = 9, file = 'out', action = 'write' )
            open( unit = 32, file = 'block.out', action = 'write' )
!        open( unit = 72, file = 'rl', action = 'write' )

         else
            rtconf % quiet = .true.
         end if
            conf_io = get_new_unit()
            open(unit=conf_io, file='conf.out', action='write')
            call print_run_conf(rtconf, '', conf_io)
            call close_unit(conf_io)
!


         if (iproc == master) then
            call fl_to_s('cell.in',tmp_str)
            flsz = size(tmp_str)
         end if


         if (nproc > 1) call mpi_bcast(flsz,1, mpi_integer, master, mpi_comm_world, ierror)
         if (iproc /= master) allocate(tmp_str(flsz))
         if (nproc > 1) call mpi_bcast(tmp_str, flsz, mpi_character, master, mpi_comm_world, ierror)
         call set_s(cell_fl,tmp_str)
         deallocate(tmp_str)

         
         if (rtconf % pot_flg == 9 .or. rtconf % pot_flg == 10) then
             call read_cell(ucell, cell_fl,.true.)
         else
             call read_cell(ucell, cell_fl,.false.)
         endif


         start_time = mpi_wtime()
         rtconf % cell => ucell

         call get_command_argument(0,cmd,cmd_len)
         if (cmd_len > 7) then
            if (cmd(cmd_len-7:cmd_len) == 'bop_elas') rtconf % mvflag = 12
         end if

         if (cmd_len > 6) then
            if (cmd(cmd_len-6:cmd_len) == 'bop_dos') then
                rtconf % mvflag = -5
            end if
         end if

         if (cmd_len > 7) then
            if (cmd(cmd_len-7:cmd_len) == 'bop_kdos') then
               rtconf % mvflag = -5
               rtconf % pot_flg = 4
            end if
         end if

         call bop(rtconf)
      endif

!        call print_run_conf(rtconf,'',800+iproc)
!        call print_cell(rtconf % cell,'',900+iproc)


      if (iproc == master) then
         end_time = mpi_wtime()


         write(6, '(/,14x,"energies/cell",8x,"energies/atom")')
         call print_ens(rtconf % tote, rtconf % avre, '', 6)

!       write(6, '(/,"energies/cell")')
!       call print_ens(rtconf % tote, '', 6)
!
!       write(6, '(/,"energies/atom")')
!       call print_ens(rtconf % avre, '', 6)

         if ( rtconf % mvflag == 12 .or. rtconf % mvflag == 36) then
         write(6, '(/,"elastics")')
         call print_elastic_consts(rtconf % elastic_const_conf, '', 6)
         end if

!           call print_elastic_const_conf(rtconf % elastic_const_conf, '', 6)
         write(6,'("Total wallclock time: ",g30.10,"s")') end_time - start_time
      end if

      deallocate(mpmap, seed)
      call close_io()
      call close_all_files()

      call mpi_finalize(ierror)

   end program main
