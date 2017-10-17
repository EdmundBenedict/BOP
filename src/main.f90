
   program main
      use mod_const
      use mod_conf
      use mod_io
      use mod_clock
      use topologia, only : iproc,nproc,mpmap,master
      implicit none
      
      integer(8) :: c1, c2
      
      type(run_conf_t) :: rtconf
      type(cell_t), target :: ucell
      
      integer :: conf_io, cell_io
      character(len=200) :: cmd 
      integer :: cmd_len
      
      type(ifl_t) :: conf_fl, cell_fl


      call init_messy_consts()
      call get_clock_prec()

      
!
!    This is a program to evaluate the Bond Order Potentials using the
!     square root terminator.

!
!     This call starts up either the scalar or parallel version depending
!     on which library has been linked in.
!
! 
!        open(unit = 8, file = 'conf.in', action = 'read' )
!        open(unit = 31, file = 'ucell.block', action = 'read' )
!        
      open( unit = 9, file = 'out', action = 'write' )
      open( unit = 32, file = 'block.out', action = 'write' )
!        open( unit = 72, file = 'rl', action = 'write' )

      
      
            
      call init_io()


      
      call fl_to_ifl('conf.in',conf_fl)   
      call read_run_conf(rtconf, conf_fl)
      
      conf_io = get_new_unit()
      open(unit=conf_io, file='conf.out', action='write')
      call print_run_conf(rtconf, '', conf_io)
      call close_unit(conf_io)
      
!       stop
      
      
!        call print_run_conf(rtconf,'',865)


      call get_command_argument(0,cmd,cmd_len)

      if (cmd_len > 7) then
         if (cmd(cmd_len-7:cmd_len) == 'bop_elas') rtconf % mvflag = 12
      end if

      if (cmd_len > 6) then
         if (cmd(cmd_len-6:cmd_len) == 'bop_dos') rtconf % mvflag = -5
      end if

      if (cmd_len > 7) then
         if (cmd(cmd_len-7:cmd_len) == 'bop_kdos') then
            rtconf % mvflag = -5
            rtconf % pot_flg = 4
         end if
      end if
   

   
      call fl_to_ifl('cell.in',cell_fl)
      
      if (rtconf % pot_flg == 9 .or. rtconf % pot_flg == 10) then
          call read_cell(ucell, cell_fl,.true.)
      else
          call read_cell(ucell, cell_fl,.false.)
      endif

      rtconf % cell => ucell
      

!        call print_cell(rtconf % cell,'',866)


      
      master = 0
      iproc = 0
      nproc = 1
      allocate(mpmap(0:1))

      call system_clock(c1)
      
      call bop(rtconf)

      call system_clock(c2)       

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
      
      
      
      write(6,'("Total wallclock time: ",g30.10,"s")'), real(c2 - c1,8)/real(cr,8)
      
      deallocate(mpmap)
      call close_io()
      call close_all_files()
      
      
   end program main
