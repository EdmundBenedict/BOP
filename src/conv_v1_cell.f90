
   program conv_v1_cell

   use mod_precision
   use mod_io
   use mod_conf
   
   implicit none
   
   character(len=1000) :: ifln
   type(cell_t) :: c
   type(ifl_t) :: ifl
   
   integer :: i, u, argc, stat
   
   
   argc = command_argument_count()
   
   if ( argc <= 0 ) then
      print *, 'Error: file not supplied'
      print *, ''
      print *, 'Usage: conv_v1_cell oldcellfile'
      print *, ''
      print *, 'The new cell will be written to stdout'
      stop
   end if
   
   call get_command_argument(1, ifln, argc, stat)
   
   if (stat /= 0) stop 'some error during retriaval of the command line arguments'
   
   
   call init_io()
      
      
   call fl_to_ifl(ifln(1:argc),ifl)
      
   c % name = 'name'
   c % str = 'str'
   
   if (.not. goto_next_uncmt(ifl,'A')) stop 'A not found'
!    print *, ifl % p
!    print *, ifl % s (ifl % p: ifl % p+10)
!    stop 'fvaf'
   call read_m(ifl, c % a)

!    print *, c % a
   c % a = transpose(c % a)
   
   if (.not. goto_next_uncmt(ifl,'LEN')) stop 'LEN not found'
   call read_m(ifl, c % lena)
   
!    print *, c % lena
   
   if (.not. goto_next_uncmt(ifl,'LATPAR')) stop 'LATPAR not found'
   call read_m(ifl, c % latpar)
   
!    print *, c % latpar
   
   if (.not. goto_next_uncmt(ifl,'ND')) stop 'ND not found'
   call read_m(ifl, c % nd)
!    print *, c % nd
   
   if (.not. goto_next_uncmt(ifl,'D')) stop 'D not found'
   allocate( c % d ( c % nd))
   
   do i = 1, c % nd
      call read_m(ifl, c % d(i) % crds)
      call read_m(ifl, c % d(i) % symb)
      call read_m(ifl, c % d(i) % mg)
      c % d(i) % de = 0.0_dp
   end do
   
   if (goto_next_uncmt(ifl,'NINERT')) then
      call read_m(ifl, c % ninert)
      if (c%ninert>0) then
         if (.not. goto_next_uncmt(ifl,'DINERT')) stop 'DINERT not found'
         allocate( c % dinert ( c % ninert))   
         do i = 1, c % ninert
            call read_m(ifl, c % dinert(i) % crds)
            call read_m(ifl, c % dinert(i) % symb)
            call read_m(ifl, c % dinert(i) % mg)
            c % dinert(i) % de = 0.0_dp
         end do      
      end if   
   end if
   
   if (goto_next_uncmt(ifl,'UNRLD')) then
      allocate( c % unrld (c % nd))   
      do i = 1, c % nd
         call read_m(ifl, c % unrld(i) % crds)
         call read_m(ifl, c % unrld(i) % symb)
         c % unrld(i) % de = 0.0_dp
         c % unrld(i) % mg = 0.0_dp
      end do         
   end if
   
   if (goto_next_uncmt(ifl,'DE')) then
      do i = 1, c % nd
         call read_m(ifl, c % d(i) % de)
      end do
   end if
   
   
   call print_cell(c,'',6)
   
   call close_io()
   
   
   end program conv_v1_cell