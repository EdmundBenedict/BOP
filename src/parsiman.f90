   
   program parsiman
      use mod_const
      use mpi
      use mod_misc
      use sa
      use mod_conf, only : run_conf_t
!         use data_t
      !       use specific
      use topologia, only: iproc, nproc, master, mpmap
      implicit none



      integer :: id, comm_size, namelen, version, subversion, ierror
      character(len=100) :: processor_name, fname
      character(len=50) :: stamp
      integer,allocatable :: seed(:), allseeds(:,:)
      integer :: u, seed_size, i, su
      real(8) :: start_time, end_time
      type(run_conf_t) :: rtconf
      type(ifl_t) :: conf_fl
      character(len=1), allocatable :: tmp_str(:)
      integer :: flsz
      
      
      call mpi_init(ierror)
      call mpi_comm_rank(mpi_comm_world, id,ierror);
      call mpi_comm_size(mpi_comm_world, comm_size, ierror);
      call mpi_get_processor_name(processor_name, namelen, ierror);
      call mpi_get_version(version, subversion,ierror);
      write(*,'(a,i0,a,i0,a,a,a,i0,a,i0)')"process ", id, " out of ", comm_size, ", name: ", trim(processor_name), &
      ", mpi version: ", version,".",subversion
      
      master = 0
      iproc = 0
      nproc = 1
      allocate(mpmap(0:1))
      
      call init_messy_consts()
      
      call init_io()
      
      if (id == master) then
         start_time = mpi_wtime()
         call get_stamp(stamp)

         print *, 'stamp: ', stamp

         u = get_new_unit()
         write(fname, '("out-",a,".out")') trim(stamp)
!             write(fname, '("out-p",i0,"-.out")') id
         open(u, file=trim(fname),action='write')
         
      end if
      
      call random_seed(size=seed_size)
      
      allocate(seed(seed_size), allseeds(seed_size, comm_size))
      
      call get_urandom(seed)
      
!         call mpi_bcast(seed, 8, mpi_integer, master,  mpi_comm_world, ierror)
      call mpi_gather( seed, seed_size, mpi_integer, &
                     & allseeds, seed_size, mpi_integer, master, mpi_comm_world, ierror)
      
      call random_seed(put=seed)
      
      if (id == master) then
         write(*,'("seed: ")')
         write(u,'("seed: ")')
         su = get_new_unit()
         open(unit=su,file='seed.out',action='write')
         write(su,*) comm_size, seed_size
         write(u,*) comm_size, seed_size
         do i = 1, comm_size
            print *, allseeds(:,i)
            write(su,*) allseeds(:,i)
            write(u,*) allseeds(:,i)
         end do
         call close_unit(su)
      end if

      if (iproc == master) then
         call fl_to_s('conf.in',tmp_str)   
         flsz = size(tmp_str)
      end if
      
      call mpi_bcast(flsz,1, mpi_integer, master, mpi_comm_world, ierror)
      if (iproc /= master) allocate(tmp_str(flsz))
      call mpi_bcast(tmp_str, flsz, mpi_character, master, mpi_comm_world, ierror)
      call set_s(conf_fl,tmp_str)
      deallocate(tmp_str)
      
      call read_run_conf(rtconf, conf_fl)
            
      call siman(u,rtconf)
      
      if (id == master) then
         call close_unit(u)
         end_time = mpi_wtime()
      end if
      
      call close_io()
      
      deallocate(mpmap)
      
      if (id == master) write(*,"('Total time: ', f30.8,'s' )") end_time - start_time
      
      call MPI_Finalize(ierror)
      
      
      
   end program parsiman
   
   