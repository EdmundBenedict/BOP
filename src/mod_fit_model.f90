
    module mod_fit_model
        
        use mod_precision
!         use math
        use mod_io
        use mod_conf
        use sa_link
        
        implicit none
        
        public
        
!         integer, parameter :: n_prop = 101, n_par = 6
    
        type :: lim_t
            character(len=10) :: name
            real(dp) :: val = 0.0_dp
            logical :: on = .false.
        end type lim_t
        
        type :: par_t
            character(len=15) :: name
            real(dp) :: val = 0.0_dp, scl = 0.0_dp
            type(lim_t) :: mx, mn
            logical :: on = .false.
        end type par_t
        
        type :: prop_t
            character(len=15) :: name
            real(dp) :: val = 0.0_dp, trg = 0.0_dp, w = 0.0_dp, x = 0.0_dp
            logical :: on = .false.
!         contains
!             procedure :: load_data => prop_data
        end type prop_t
        
        
        type :: conf_t
            integer :: npar = 0, nprop = 0, counter = 0
            character(len=10) :: name
            real(dp) :: val = 0.0_dp
            type(par_t), allocatable :: par(:) 
            type(prop_t), allocatable :: prop(:) 
            integer :: proc = 99999
            character(len=200) :: rund
!         contains
!             procedure :: eval => conf_eval
!             procedure :: load_data => conf_prop_data
!             procedure :: load_params => load_conf_params
!             procedure :: write_params => write_conf_params
        end type conf_t
        

    contains
        
        
        subroutine print_par(par, l, u)
            type(par_t), intent(in) :: par
            character(len=*), intent(in) :: l
            integer, intent(in) :: u
            
            write(u,'(a,x,a,x,l,x,g20.12,x,f12.6,2(x,l,x,f12.6))') l, par % name, par % on, &
            & par % val, par % scl, par % mn % on, par % mn % val ,  par % mx % on, par % mx % val  
            
            
        end subroutine print_par
        
        
        
        subroutine print_prop(prop, l, u)
            type(prop_t), intent(in) :: prop
            character(len=*), intent(in) :: l
            integer, intent(in) :: u
            
            write(u,'(a,x,a,x,l,3(x,g20.12))') l, prop%name, prop % on, &
            & prop % trg, prop % val, prop % w
            
        
        end subroutine print_prop
        
        
        
        
        subroutine print_conf(conf, l, u)
            type(conf_t), intent(in) :: conf
            character(len=*), intent(in) :: l
            integer, intent(in) :: u
            
            integer :: i
            
            write(u,'(a,x,i0,x,a)') l, conf % npar, &
            &' name             sw    val                  scl       mnsw  min       mxsw  max'
            do i=1,conf % npar
                call print_par(conf % par(i), l//'    ', u)
            end do
            
            write(u,'(a,x,i0,x,a)') l, conf % nprop, &
            &' name             sw    target               val                  weight'
            do i=1,conf % nprop
                call print_prop(conf % prop(i), l//'    ', u)
            end do
            
            write(u,'(a,x,a,x,g20.12,2x,a,x,i0)') l, conf % name, conf % val, 'c:', conf%counter
            
        
        end subroutine print_conf
        
        
        
        
        
        
        subroutine read_par(par, u)
            type(par_t), intent(out) :: par
            integer, intent(in) :: u
            
            
            
            read(u, * ) par % name, par % on, &
            & par % val, par % scl, par % mn % on, par % mn % val ,  par % mx % on, par % mx % val  
            
            
        end subroutine read_par
        
        
        
        subroutine read_prop(prop, u)
            type(prop_t), intent(out) :: prop
            integer, intent(in) :: u
            
            read(u,*) prop%name, prop % on, &
            & prop % trg, prop % val, prop % w
            
        
        end subroutine read_prop
        
        
        
        
        subroutine read_conf(conf, iu)
            type(conf_t), intent(inout) :: conf
            integer, intent(in), optional :: iu
            
            integer :: i, u
            
            if (.not. present(iu)) then
                u = get_new_unit()
                open(unit=u, file='model.in', action='read')
            else
                u = iu
            end if
            
            read(u,*) conf % npar
            
            if (.not. allocated(conf%par)) then
                allocate(conf%par(conf%npar))
            else
                if (size(conf%par) /= conf % npar) then
                    deallocate(conf%par)
                    allocate(conf%par(conf%npar))
                end if
            end if
            
            do i=1,conf % npar
                call read_par(conf % par(i), u)
            end do
            
            read(u,*) conf % nprop
            
            if (.not. allocated(conf%prop)) then
                allocate(conf%prop(conf%nprop))
            else
                if (size(conf%prop) /= conf % nprop) then
                    deallocate(conf%prop)
                    allocate(conf%prop(conf%nprop))
                end if
            end if
            
            do i=1,conf % nprop
                call read_prop(conf % prop(i), u)
            end do
            
            read(u,*) conf % name, conf % val
            
            if (.not. present(iu)) then
                call close_unit(u)
            end if
            
        
        end subroutine read_conf
        
        
        
        
        
        
        subroutine decode(rtconf, fitconf)
            type(run_conf_t), intent(inout) :: rtconf
            type(conf_t), intent(in) :: fitconf
            
!             Al_cc, Ti_cc, 
!             pps_v, pps_n, pps_rc, pps_nc, 
!             ppp_v, ppp_n, ppp_rc, ppp_nc,
!             pds_v, pds_n, pds_rc, pds_nc,
!             pdp_v, pdp_n, pdp_rc, pdp_nc,
!             dds_v, dds_n, dds_rc, dds_nc,
!             ddp_v, ddp_n, ddp_rc, ddp_nc,
!             ddd_v, ddd_n, ddd_rc, ddd_nc
            
            
!             rtconf % ham_conf % a(0) % b % cc  = fitconf % par( 1) % val
!             rtconf % ham_conf % a(1) % b % cc  = fitconf % par( 2) % val
!             
!             rtconf % ham_conf % b(0) % pps % v  = fitconf % par(3) % val            
!             rtconf % ham_conf % b(0) % ppp % v  = fitconf % par(4) % val            
!            
!             rtconf % ham_conf % b(1) % pds % v  = fitconf % par(5) % val            
!             rtconf % ham_conf % b(1) % pdp % v  = fitconf % par(6) % val
!                         
!             rtconf % ham_conf % b(2) % dds % v  = fitconf % par(7) % val            
!             rtconf % ham_conf % b(2) % ddp % v  = fitconf % par(8) % val            
!             rtconf % ham_conf % b(2) % ddd % v  = fitconf % par(9) % val
            
           rtconf % ham_conf % a(0) % e % lam_0  = fitconf % par(1) % val           
           rtconf % ham_conf % a(0) % e % r_core = fitconf % par(2) % val           

           rtconf % ham_conf % a(1) % e % lam_0  = fitconf % par(3) % val
           rtconf % ham_conf % a(1) % e % r_core = fitconf % par(4) % val           

           rtconf % ham_conf % b(0) % a   = fitconf % par(5) % val
           rtconf % ham_conf % b(0) % c   = fitconf % par(6) % val
           rtconf % ham_conf % b(0) % nu  = fitconf % par(7) % val
          
           rtconf % ham_conf % b(2) % a   = fitconf % par(8) % val
           rtconf % ham_conf % b(2) % c   = fitconf % par(9) % val
           rtconf % ham_conf % b(2) % nu  = fitconf % par(10) % val                                                                                      
                    
        end subroutine decode
        
        
        
        
        subroutine conf_eval(confin,rtconf)
            use topologia, only : iproc,master

            type(conf_t), intent(inout) :: confin
            type(run_conf_t), intent(inout) :: rtconf
            

            type(cell_patch), save, target :: d19, d24, l12, d22, l10, b19
            type(elastic_const_t), save :: d19_bel, l10_bel
            real(dp), save :: d19_ben, d24_ben, l12_ben, d22_ben, l10_ben, b19_ben

            real(dp) :: f
            integer, save :: c = 1
            integer :: i, rt
            real(dp) :: ti_al_diff
            
!             integer(8) :: seed(2)
            character(len=20) :: fnfmt
            real(dp) :: tmp

            type(ifl_t) :: cell_fl
            
            f = 0.0_dp

            call decode(rtconf, confin)

!             if (.not. (abs(rtconf % ham_conf % b(0) % pps % v) > abs(rtconf % ham_conf % b(0) % ppp % v) &
!               & .and.  abs(rtconf % ham_conf % b(1) % pds % v) > abs(rtconf % ham_conf % b(1) % pdp % v) )) then
!                 confin % val = 10000.0_dp
!                 c = c + 1
!                 return
!             end if

            if (c == 1) then
            
                call fl_to_ifl('ucell.d19',cell_fl)
                call read_cell( d19 % cell,cell_fl)
                
                call fl_to_ifl('ucell.l10',cell_fl)
                call read_cell( l10 % cell,cell_fl)

                call fl_to_ifl('ucell.b19',cell_fl)
                call read_cell( b19 % cell,cell_fl)
                
                call fl_to_ifl('ucell.l12',cell_fl)                
                call read_cell( l12 % cell,cell_fl)
                
                call fl_to_ifl('ucell.d22',cell_fl)                
                call read_cell( d22 % cell,cell_fl)
                
                call fl_to_ifl('ucell.d24',cell_fl)
                call read_cell( d24 % cell,cell_fl)

!                print *, 'start bond'
               
               rtconf % fs_only = 0
               rtconf % fs = 0
               rtconf % env = 0
               rtconf % vpair = 0


                call bop_el(rtconf, d19)               
                call bop_el(rtconf, l10)

               d19_ben = d19 % en % elct
               d19_bel  = d19 % elcon 

               l10_ben = l10 % en % elct
               l10_bel = l10 % elcon 

                call bop_en(rtconf, b19)               
                call bop_en(rtconf, l12)
                call bop_en(rtconf, d22)               
                call bop_en(rtconf, d24)
               
               b19_ben = b19 % en % elct
               l12_ben = l12 % en % elct
               d22_ben = d22 % en % elct
               d24_ben = d24 % en % elct

               rtconf % fs_only = 1
               rtconf % fs = 1
               rtconf % env = 1
               rtconf % vpair = 0

!                 print *, 'done bond'
                

            end if


      print *, rtconf % fs_only, rtconf % fs, rtconf % env, rtconf % vpair

!             if ( iproc == master ) then
!                iounits(1) = get_new_unit()
!                write(fnfmt,'("model2conf.out_",i0,"_",i0)') ,iproc,c
!                open(unit=iounits(1), file=trim(fnfmt), action='write')
!                call print_run_conf(rtconf, '', iounits(1))
!                call close_unit(iounits(1))
!                print  *,'wrote '//trim(fnfmt)//' exitting..'
!                stop 0
!             end if
            


            call bop_el(rtconf, d19)
            call bop_el(rtconf, l10)

            d19 % en % bind = d19 % en % bind + d19_ben

            l10 % en % bind = l10 % en % bind + l10_ben            

            d19 % elcon % s11 % tot = d19 % elcon % s11 % tot + d19_bel % s11 % elc
            d19 % elcon % s33 % tot = d19 % elcon % s33 % tot + d19_bel % s33 % elc

            l10 % elcon % s11 % tot = l10 % elcon % s11 % tot + l10_bel % s11 % elc
            l10 % elcon % s33 % tot = l10 % elcon % s33 % tot + l10_bel % s33 % elc

            d19 % elcon % c11 % tot = d19 % elcon % c11 % tot + d19_bel % c11 % elc
            d19 % elcon % c33 % tot = d19 % elcon % c33 % tot + d19_bel % c33 % elc
            d19 % elcon % c12 % tot = d19 % elcon % c12 % tot + d19_bel % c12 % elc
            d19 % elcon % c13 % tot = d19 % elcon % c13 % tot + d19_bel % c13 % elc
            d19 % elcon % c44 % tot = d19 % elcon % c44 % tot + d19_bel % c44 % elc
            d19 % elcon % c66 % tot = d19 % elcon % c66 % tot + d19_bel % c66 % elc
            
            l10 % elcon % c11 % tot = l10 % elcon % c11 % tot + l10_bel % c11 % elc
            l10 % elcon % c33 % tot = l10 % elcon % c33 % tot + l10_bel % c33 % elc
            l10 % elcon % c12 % tot = l10 % elcon % c12 % tot + l10_bel % c12 % elc
            l10 % elcon % c13 % tot = l10 % elcon % c13 % tot + l10_bel % c13 % elc
            l10 % elcon % c44 % tot = l10 % elcon % c44 % tot + l10_bel % c44 % elc
            l10 % elcon % c66 % tot = l10 % elcon % c66 % tot + l10_bel % c66 % elc

            d19 % elcon % c1344 % tot = d19 % elcon % c1344 % tot + d19_bel % c1344 % elc
            d19 % elcon % c1266 % tot = d19 % elcon % c1266 % tot + d19_bel % c1266 % elc

            l10 % elcon % c1344 % tot = l10 % elcon % c1344 % tot + l10_bel % c1344 % elc
            l10 % elcon % c1266 % tot = l10 % elcon % c1266 % tot + l10_bel % c1266 % elc


            call bop_en(rtconf, b19)               
            b19 % en % bind = b19 % en % bind + b19_ben

            call bop_en(rtconf, l12)
            l12 % en % bind = l12 % en % bind + l12_ben

            call bop_en(rtconf, d22)               
            d22 % en % bind = d22 % en % bind + d22_ben

            call bop_en(rtconf, d24)
            d24 % en % bind = d24 % en % bind + d24_ben


            confin % prop( 1) % val = d19 % en % bind

            confin % prop( 2) % val = d24 % en % bind

            confin % prop( 3) % val = l12 % en % bind

            confin % prop( 4) % val = d22 % en % bind

            confin % prop( 5) % val = l10 % en % bind

            confin % prop( 6) % val = b19 % en % bind

            confin % prop( 7) % val = d19 % elcon % s11 % tot
            confin % prop( 8) % val = d19 % elcon % s33 % tot

            confin % prop( 9) % val = l10 % elcon % s11 % tot
            confin % prop(10) % val = l10 % elcon % s33 % tot

            confin % prop(11) % val = d19 % elcon % c11 % tot
            confin % prop(12) % val = d19 % elcon % c33 % tot
            confin % prop(13) % val = d19 % elcon % c12 % tot
            confin % prop(14) % val = d19 % elcon % c13 % tot
            confin % prop(15) % val = d19 % elcon % c44 % tot
            confin % prop(16) % val = d19 % elcon % c66 % tot
            
            confin % prop(17) % val = l10 % elcon % c11 % tot
            confin % prop(18) % val = l10 % elcon % c33 % tot
            confin % prop(19) % val = l10 % elcon % c12 % tot
            confin % prop(20) % val = l10 % elcon % c13 % tot
            confin % prop(21) % val = l10 % elcon % c44 % tot
            confin % prop(22) % val = l10 % elcon % c66 % tot
                          
            confin % prop(23) % val = d19 % elcon % c1344 % tot
            confin % prop(24) % val = d19 % elcon % c1266 % tot

            confin % prop(25) % val = l10 % elcon % c1344 % tot
            confin % prop(26) % val = l10 % elcon % c1266 % tot


      


!            f = f+  abs(d19 % elcon % c1344 % tot - confin % prop(23) % trg)*confin % prop(23) % w
!            f = f+  abs(d19 % elcon % c1266 % tot - confin % prop(24) % trg)*confin % prop(24) % w
! 
!            f = f+  abs(l10 % elcon % c1344 % tot - confin % prop(25) % trg)*confin % prop(25) % w
!            f = f+  abs(l10 % elcon % c1266 % tot - confin % prop(26) % trg)*confin % prop(26) % w


           do i=7,10
               if (confin % prop(i) % val < 0.0_dp .and. confin % prop(i) % on) f = f+  abs(confin % prop(i) % val)*100.0_dp
           end do 


! 
           do i=23,26 !confin % nprop
!               tmp = confin % prop(i) % val-3.0_dp
!               if (tmp > 0) f = f+ 100.0_dp*tmp
              if (confin % prop(i) % on) f = f+ abs(confin % prop(i) % val - confin % prop(i) % trg)*confin % prop(i) % w
           end do
        
            confin % val = f
        
            c = c + 1
            
        
        end subroutine conf_eval
        
    end module mod_fit_model










