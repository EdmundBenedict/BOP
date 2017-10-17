

    module ssa
        use topologia, only : master, iproc
        use mpi, only : mpi_wtime
        use mod_precision
        use mod_io
        use mod_fit_model
        use mod_all_scalar, only : quiet
        
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
        
        call close_unit(u)
        
        if (iproc == master) print *, 'sa_conf: ', init_t, min_t, max_t_steps, max_isoterm_steps, sampling_points
        
    end subroutine load_sa_conf
    
    subroutine neighbour(neigh, base, t)
        type(conf_t), intent(inout) :: neigh
        type(conf_t), intent(in) :: base
        real(dp), intent(in) :: t
        
        real(dp), allocatable :: rnd(:)
        
        allocate(rnd(base%npar))
        
        
        call random_number(rnd)
        

        where (base % par % on)
            neigh % par % val = base % par % val + base % par % scl * (rnd - 0.5_dp) * t
        end where
        
        where (base % par % on .and. base % par % mn % on)
            neigh % par % val = base % par % mn % val + abs( base % par % mn % val - neigh % par % val )
        end where 
        
        where (base % par % on .and. base % par % mx % on)
            neigh % par % val = base % par % mx % val - abs( base % par % mx % val - neigh % par % val )
        end where 
        
        deallocate(rnd)
        
    end subroutine neighbour
    
    
    subroutine ssiman(rtconf,u)
        use mod_conf, only : run_conf_t

        type(run_conf_t), intent(inout) :: rtconf        
        integer, intent(in) :: u


        integer :: ierror
        integer,allocatable :: req(:), stat(:,:)
        real(8) :: start_time, end_time, t_start_time
        
        
        real(dp) :: t
        real(dp) ::  prnd, prob
        
        integer(8) :: counter, max_counter
        integer :: iso_step, t_step
        
        
        type(conf_t) :: pres, new, opt

!         real(dp), allocatable :: fvals(:)

        integer :: i, j, k, n, min_ind(1), min_proc
        real(dp) :: time
        

        
        call load_sa_conf()

        max_counter = max_t_steps * max_isoterm_steps


        call read_conf(pres)
        pres % val = huge(1.0_dp)

!         if (iproc == master) call print_conf(pres, '', 6)
        
        call conf_eval(pres,rtconf)
        if (iproc == master) call print_conf(pres, '', 6)
       
        
        allocate(new%par(pres%npar),   opt%par(pres%npar))
        allocate(new%prop(pres%nprop), opt%prop(pres%nprop))


        
        t = init_t
        opt = pres
        new = pres
        
        counter = 0

        t_start_time = mpi_wtime()
        do t_step = 1,max_t_steps
            do iso_step = 1,max_isoterm_steps
               if (iproc == master .and. mod(counter,1000) == 0 ) &
                  & write(6,'("t_step, iso_step, counter, wtime :",3(x,i0),x,f20.5)') &
                  & t_step, iso_step, counter, mpi_wtime() - t_start_time
               new % counter = counter
               call neighbour(new, pres, t)
               
               call conf_eval(new,rtconf)
               
               if (new % val < opt % val) then
                  opt = new
                  if (iproc == master) then
                     write(u,'("move:")')
                     call print_conf(opt, '', u)
                  end if
!                         write(u,*)'F: ', opt % val
!                         write(u,*)'par: ', opt % par % val
                  pres = new
               else
                  call random_number(prnd)

                  prob = exp((pres % val - new % val)/t)
                  
                  if (prob > prnd) then
                     pres = new
!                     if (iproc == master) then
!                        write(u,'("explore:")')
!                        call print_conf(pres, '', u)
!                     end if
                  end if    
               end if
               
               counter = counter + 1
            end do
!             write(10,*) t
!            t = init_t*exp(log(min_t/init_t)*real(t_step,dp)/real(max_t_steps,dp))
            t = (min_t-init_t)*real(t_step,dp)/real(max_t_steps,dp) + init_t
            
        end do
        
        if (iproc == master) then
            call print_conf(opt, '', u)
            write(u, '("c:",i0)') counter
!             print *, 'c,best: ',counter, opt % par % val
!             call plot_conf(opt,'opt_plot')
        end if
        

        deallocate(new%par , opt%par, new%prop, opt%prop)


    end subroutine ssiman
    
    
    
    
    end module ssa
    
    
