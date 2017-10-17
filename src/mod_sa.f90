

    module sa
        use mpi
        use distrib
        use mod_precision
        use mod_io
        use mod_fit_model
!         use specific
        
        implicit none
        
        
        
        real(dp) :: init_t, min_t
        integer :: max_isoterm_steps, sampling_points, max_t_steps
        
        
        
    contains
    
        
    subroutine load_sa_conf()
        integer :: u
        
        u = get_new_unit()
        open(u, file='sac.in', action='read')
        
        call get_var(init_t, 'init_t', u)
        call get_var(min_t, 'min_t', u)
        call get_var(max_t_steps, 'max_t_steps', u)
        call get_var(max_isoterm_steps, 'max_isoterm_steps', u)
        call get_var(sampling_points, 'sampling_points', u)
        
        call close_unit(u)
        
        print *, 'sa_conf: ', init_t, min_t, max_t_steps, max_isoterm_steps, sampling_points
        
    end subroutine load_sa_conf
    
    subroutine neighbour(neigh, base, t)
        type(conf_t), intent(inout) :: neigh
        type(conf_t), intent(in) :: base
        real(dp), intent(in) :: t
        
        real(dp), allocatable :: rnd(:)
        
        allocate(rnd(base%npar))
        
        
        call random_number(rnd)
        

        where (base % par % on)
            neigh % par % val = base % par % val + base % par % scl * (rnd - 0.5_dp) !* t
        end where
        
        where (base % par % on .and. base % par % mn % on)
            neigh % par % val = base % par % mn % val + abs( base % par % mn % val - neigh % par % val )
        end where 
        
        where (base % par % on .and. base % par % mx % on)
            neigh % par % val = base % par % mx % val - abs( base % par % mx % val - neigh % par % val )
        end where 
        
        deallocate(rnd)
        
    end subroutine neighbour
    
    
    subroutine siman(u,rtconf)
        
        type(run_conf_t), intent(inout) :: rtconf
        integer :: ierror,task_increment,win,task_counter,end_task, master, rank, smp_counter
        integer,allocatable :: req(:), stat(:,:)
        integer(kind=mpi_address_kind) :: win_size
        real(8) :: start_time, end_time, t_start_time
        integer, intent(in) :: u
        
        real(dp) :: t
        real(dp) ::  prnd, prob, reach
        
        integer(8) :: counter, max_counter
        integer :: iso_step, t_step
        
        
        real(dp),allocatable :: prop_vals(:,:), par_vals(:,:)
        
        type(conf_t), allocatable :: smpl(:)
        type(conf_t) :: init, pres, new, opt

!         real(dp), allocatable :: fvals(:)

        integer :: i, j, k, n, smp, min_ind(1), min_proc
        real(dp) :: time
        
        call mpi_comm_rank(mpi_comm_world, rank, ierror)
        
        
        
        master = 0
        
        reach = 1.0_dp

        call load_sa_conf()
        

        if (rank == master) then
            win_size = 1
        else
            win_size = 0
        end if
!         write(0,*) 'id, size: ', rank, win_size
        call mpi_win_create(smp_counter, win_size, 1, mpi_info_null, mpi_comm_world, win,ierror)
!         write(0,*) 'id, err:', rank, ierror


!         call CreateCounter( rank, master, smp_counter, 1, win, mpi_comm_world)


        max_counter = max_t_steps*max_isoterm_steps*sampling_points


        call read_conf(init)
        init % val = huge(1.0_dp)

        if (rank == master) call print_conf(init, '', 6)
        
        call conf_eval(init,rtconf)
!         call print_conf(init, '', 6)
       
!         call mpi_bcast(init % val, 1, mpi_real8, master,  mpi_comm_world, ierror )
!         call mpi_bcast(init % prop % val, init%nprop, mpi_real8, master,  mpi_comm_world, ierror)
!         call mpi_barrier(mpi_comm_world, ierror)
!         print *,'after sampl:', rank
!         stop 
        
        allocate(pres%par(init%npar),   new%par(init%npar),   opt%par(init%npar))
        allocate(pres%prop(init%nprop), new%prop(init%nprop), opt%prop(init%nprop))

        allocate(smpl(sampling_points))
        allocate(req(4*sampling_points), stat(mpi_status_size, 4*sampling_points))
        allocate(prop_vals(init%nprop, sampling_points), par_vals(init%npar, sampling_points))

        do smp=1,sampling_points
            allocate(smpl(smp)%par(init%npar), smpl(smp)%prop(init%nprop))
        end do


        
        t = init_t
        pres = init
        opt = init
        new = init
        do smp=1,sampling_points
            smpl(smp) = init
        end do
!         print *, 'allocated pres: ', allocated(pres%par)
!         call write_conf_params(pres)
        
        counter = 0
        task_increment=1
!         call plot_conf(init,'init_plot')
!         call CreateCounter(rank,winRank,smp_counter,1,win,MPI_COMM_WORLD)
        t_start_time = mpi_wtime()
        do t_step = 1,max_t_steps
            do iso_step = 1,max_isoterm_steps
                
                
                
                if (rank == master) then 
                    write(6,'("t_step, iso_step, counter, wtime :",3(x,i0),x,f30.6)') &
                               & t_step, iso_step, counter, mpi_wtime() - t_start_time
                    call resetcounter(master, 1 , win)
                    do smp = 1,sampling_points
                        call mpi_irecv(smpl(smp) % proc,1,mpi_integer,mpi_any_source,smp,mpi_comm_world,req(smp),ierror)
!                         
                        call mpi_irecv(smpl(smp) % val,1,mpi_real8, &
                                        & mpi_any_source, smp+sampling_points, mpi_comm_world, &
                                        & req(smp+sampling_points),ierror)
                       
                        
                        call mpi_irecv(prop_vals(1,smp), smpl(smp) % nprop, mpi_real8, &
                                        & mpi_any_source, smp+2*sampling_points, mpi_comm_world, &
                                        & req(smp+2*sampling_points),ierror)
                        
                        call mpi_irecv(par_vals(1,smp), smpl(smp) % npar, mpi_real8, &
                                        & mpi_any_source, smp+3*sampling_points, mpi_comm_world, &
                                        & req(smp+3*sampling_points),ierror)
                        
                    end do
                end if
                
                
                
                call mpi_barrier(mpi_comm_world, ierror)
                
                if (rank /= master ) then
!                 slaves'
                    do
                        
                        start_time = mpi_wtime()
                        task_counter= GetCounter(master, task_increment , win)
                        end_time = mpi_wtime()
                        if ((end_time - start_time) > 1. ) print *, 'rank, getcounter_time: ', rank,  end_time - start_time
                        if (task_counter > sampling_points) exit;
                        
                        end_task = min( task_counter + task_increment - 1, sampling_points )
                        
!                         write(*,'("Process ",i0," obtained the following task(s): ",i0," - ", i0," t, isot ",i0,x,i0)')&
!                         & rank, task_counter, end_task, t_step, iso_step
                        
                        do smp = task_counter, end_task
    !                             write(*,'("Process ", i0, " executes the following task: ", i0)') rank, smp
                            
                            smpl(smp) % proc = rank
    !                             print *, 'send_proc', rank
                            call mpi_irsend(smpl(smp) % proc, 1, mpi_integer, master, smp, mpi_comm_world, req(smp), ierror )
    !                             print *, 'after send_proc', rank
!                             smpl(smp) % val = opt % val
                            smpl(smp) % counter = counter
                            call neighbour(smpl(smp), pres, t)
                            
                            call conf_eval(smpl(smp),rtconf)
                            
                            prop_vals(:,smp) = smpl(smp) % prop % val
                            par_vals(:,smp) = smpl(smp) % par % val
                            
                            call mpi_irsend(smpl(smp) % val, 1, mpi_real8, &
                                            & master, smp+sampling_points, mpi_comm_world, &
                                            & req(smp+sampling_points), ierror )
                          
                            
                            call mpi_irsend(prop_vals(1,smp), smpl(smp) % nprop, mpi_real8, &
                                            & master, smp+2*sampling_points, mpi_comm_world, &
                                            & req(smp+2*sampling_points), ierror )
                            
                            call mpi_irsend(par_vals(:,smp), smpl(smp) % npar, mpi_real8, &
                                            & master, smp+3*sampling_points, mpi_comm_world, &
                                            & req(smp+3*sampling_points), ierror )
                                            
                            
                        enddo
                        
    !                     task_increment=2
                    enddo
                else
!                 master's    
                    
                    call mpi_waitall(sampling_points*4, req, stat, ierror)                    
                    
                    min_ind = minloc(smpl % val)
                    
                    do smp=1,sampling_points
                        smpl(smp) % par % val = par_vals(:,smp)
                        smpl(smp) % prop % val = prop_vals(:,smp)
                    end do
                    
                    new = smpl(min_ind(1))
                    
                    if (new % val < opt % val) then
                        opt = new
                        write(u,'("move:")')
                        call print_conf(opt, '', u)
!                         write(u,*)'F: ', opt % val
!                         write(u,*)'par: ', opt % par % val
                        pres = new
                    else
                        call random_number(prnd)
    
                        prob = exp((pres % val - new % val)/t)
                        
                        if (prob > prnd) then
                            pres = new
                            write(u,'("explore:")')
                            call print_conf(pres, '', u)
                        end if    
                    end if
                end if
                
                call mpi_barrier(mpi_comm_world, ierror)
                start_time = mpi_wtime()
                call mpi_bcast(pres % val, 1, mpi_real8, master,  mpi_comm_world, ierror )
                call mpi_bcast(pres % par % val, pres % npar, mpi_real8, master,  mpi_comm_world, ierror)
                call mpi_bcast(pres % prop % val, pres % nprop, mpi_real8, master,  mpi_comm_world, ierror)
                end_time = mpi_wtime()
                call mpi_barrier(mpi_comm_world, ierror)
                if ((end_time - start_time) > 1. ) print *,"bcasts' time: ", end_time - start_time
                
                counter = counter + 1
            end do
!             write(10,*) t
            t = init_t*exp(log(min_t/init_t)*real(t_step,dp)/real(max_t_steps,dp))
            
        end do
        
!         call DestroyCounter(rank,winRank,smp_counter,win)
        
        call MPI_Win_free(win,ierror)
        
        if (rank == master) then
            call print_conf(opt, '', u)
            write(u, '("c:",i0)') counter
!             print *, 'c,best: ',counter, opt % par % val
!             call plot_conf(opt,'opt_plot')
        end if
        
    end subroutine siman
    
    
    
    
    end module sa
    
    
