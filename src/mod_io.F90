
   module mod_io
   
      use mod_precision
   !        use mod_clock

      implicit none

      private
      
#ifndef no_protecteds      
      integer, protected :: fin = 5, fout = 6, ferr = 0
#else
      integer :: fin = 5, fout = 6, ferr = 0
#endif

      logical, allocatable :: io_units(:)

      public :: init_io, get_new_unit, close_unit, set_in, set_out, set_err, close_io, fl_to_s
      public :: fin, fout, ferr
            
      type ifl_t
         integer :: p, n  ! p: position; n: size of s
         character(len=2**20) :: s
      end type
      
      interface read_m
#ifndef no_quad       
         module procedure    read_i1, read_i2, read_i4, read_i8, &
                           & read_r4, read_r8, read_r16, &
                           & read_z4, read_z8, read_z16, &
                           & read_c , read_l, &
                           & read_i1a, read_i2a, read_i4a, read_i8a, &
                           & read_r4a, read_r8a, read_r16a, &
                           & read_z4a, read_z8a, read_z16a, &
                           & read_ca , read_la, &
                           & read_r8a2
#else
         module procedure    read_i1, read_i2, read_i4, read_i8, &
                           & read_r4, read_r8, &
                           & read_z4, read_z8, &
                           & read_c , read_l, &
                           & read_i1a, read_i2a, read_i4a, read_i8a, &
                           & read_r4a, read_r8a, &
                           & read_z4a, read_z8a, &
                           & read_ca , read_la, &
                           & read_r8a2
#endif                           
                           
      end interface read_m
      
      public :: set_s, read_m, ifl_t, fl_to_ifl, goto_next_uncmt
      
      
      
      interface print_c
         module procedure print_v, print_mr, print_mz, print_lv, print_iv
      end interface print_c

      interface get_var
         module procedure  get_var_r, get_var_i, get_var_ra, get_var_ia, get_var_a, get_var_aa, get_var_l, get_var_la
      end interface get_var
      
      public :: print_c, print_v, print_mr, print_mz, print_lv, print_iv
      public :: get_var, get_var_r, get_var_i, get_var_ra, get_var_ia, get_var_a, get_var_aa, get_var_l, get_var_la
      
   contains

   
   subroutine init_io()
      allocate(io_units(100:299))
      io_units = .false.
   end subroutine init_io

   function get_new_unit()
      integer :: get_new_unit
      integer :: i
      i = 100
      do while (io_units(i) .and. (i<300))
            i = i + 1
      end do
      if (i == 300) write(ferr, '("Run out of free units")')
      io_units(i) = .true.
      get_new_unit =  i
      
   end function get_new_unit
   
   subroutine close_unit(u)
      integer, intent(in) :: u
      close(u)
      io_units(u) = .false.
   end subroutine close_unit

   subroutine close_io()
      integer :: u
      do u=100,299
         if (io_units(u)) close(u)
      end do
      deallocate(io_units)
   end subroutine close_io


   subroutine set_in(u)
      integer, intent(in) :: u
      fin = u
   end subroutine set_in

   subroutine set_out(u)
      integer, intent(in) :: u
      fout = u
   end subroutine set_out
   
   subroutine set_err(u)
      integer, intent(in) :: u
      ferr = u
   end subroutine set_err
   
   
   subroutine fl_to_s(fln,s)
      use iso_c_binding, only : c_new_line, c_carriage_return 
      character(len=*), intent(in) :: fln
      character(len=1), intent(inout), allocatable :: s(:)
      integer :: u, stat, sz
      
#ifdef no_f03_filefuns   
      character(len=1), allocatable :: stmp(:)
      integer :: i
#endif
      
      u = get_new_unit()

#ifndef no_f03_filefuns
      open(unit=u, form='unformatted',file=trim(fln), action='read',access='stream')      
      inquire(unit=u,size=sz)
      sz = sz + 1
      allocate(s(sz))
      s(1) = ' '
      read(u,iostat=stat) s(2:)
      
#else
! old f is nasty
      sz = 32*2**20
      open(unit=u, form='unformatted',file=trim(fln), action='read',recl=sz-1)      
      allocate(stmp(sz)) ! 32 MiB should be enough for input text file, 1 MiB reserved for ending, see a few lines lower
      stmp = 'e'
      stmp(1) = ' '
      read(u,iostat=stat) stmp(2:)
      stat = 0
      
!       find the end marker, 1MiB of 'e'
      do i=1,size(stmp)
         if (stmp(i) == 'e') then
            stat = stat + 1
            if (stat > 2**20) then 
               sz = i - stat
               exit
            end if
         else
            stat = 0
         end if
      end do
      allocate(s(sz))
      s = stmp(:sz)
      deallocate(stmp)
#endif
            
      call close_unit(u)
      
!       dump dumb line ends
      where (s == c_carriage_return) s = c_new_line

#ifdef ifort_read_bug      
      call please_ifort(s)  ! hopefully not for long    
#endif

   end subroutine fl_to_s
   
   
#ifdef ifort_read_bug
   subroutine please_ifort(s)
!  ifort compiled executables issue error and exit 
!   if there is no empty space (' ') after a string being read during read(u,..) operation.
!   This is probably very old. The first ref I have seen is in Mark Cawkwell's documentation for BOP. He suspected it is peculiarity of BOP itself.
!   Someone should file a bug report to intel sometime. 
!  Until then, this works around the problem

      use iso_c_binding, only : c_new_line
      character(len=1), allocatable, intent(inout) :: s(:)
      character(len=1), allocatable :: t(:)
            
      integer :: i, j, n
      
!       print *, 'allocated(s)', allocated(s)
      n = size(s)
      
      j = n + count(s == c_new_line) + 1
      
      allocate(t(j))
      
      j = 1
      do i = 1, n         
         if (s(i) == c_new_line) then
            t(j) = ' '
            j = j + 1
         end if         
         t(j) = s(i)
         j = j + 1
      end do
      t(j) = ' '
      
      deallocate(s)
      allocate(s(j))
      
      s = t(:j)
      
      deallocate(t)
   end subroutine please_ifort
#endif      
   
   
   subroutine fl_to_ifl(fln,ifl)
      character(len=*), intent(in) :: fln
      type(ifl_t), intent(inout) :: ifl
      character(len=1), allocatable :: s(:)
      
      
      call fl_to_s(trim(fln),s)
      call set_s(ifl,s)
      deallocate(s)
      
   end subroutine fl_to_ifl
   
   
! Hopefully is_blank and is_comment will be inlined with with -O2/3.
! The probability of inlining may decrease if they are exported
   
   function is_blank(c)
!             This can probably be replaced with the intrinsic 'scan' at later point
      use iso_c_binding, only :  &
         &   c_null_char      , & ! \0
         &   c_alert          , & ! \a
         &   c_backspace      , & ! \b
         &   c_form_feed      , & ! \f
         &   c_new_line       , & ! \n
         &   c_carriage_return, & ! \r
         &   c_horizontal_tab , & ! \t
         &   c_vertical_tab       ! \v

      character(len=1), intent(in) :: c
      logical :: is_blank
      is_blank =  ((c == ' '              ) &
            & .or. (c == c_null_char      ) &                      
            & .or. (c == c_alert          ) &                    
            & .or. (c == c_backspace      ) &     
            & .or. (c == c_form_feed      ) &     
            & .or. (c == c_new_line       ) &     
            & .or. (c == c_carriage_return) &     
            & .or. (c == c_horizontal_tab ) &     
            & .or. (c == c_vertical_tab   )) 
   end function is_blank     
   
   function is_comment(c)
!             This can probably be replaced with the intrinsic 'scan' at later point
      character(len=1), intent(in) :: c
      logical :: is_comment
      is_comment =  ((c == '#') .or. (c == '!')) 
   end function is_comment
   
   
   subroutine set_s(f,cs)
! !  copy the char array cs to the internal string and build map pointing to all noncommented words
      type(ifl_t), intent(inout) :: f        ! may be not a bad idea to turn this to class and add the procedure to the type def.
      character(len=1), intent(in) :: cs(:)
      integer :: i
   
      f % n = size(cs)
      if (f%n > len(f%s)) then
         write(ferr,"('Increase the static size in ifl_t%s to at least ',i0'. Exitting ..')") f%n
         stop 255
      end if
      do i = 1, f % n
         f % s(i:i) = cs(i)
      end do
      
      f % p = 1
   end subroutine set_s

   function next_wi(f,n) result(i)
      use iso_c_binding, only : c_new_line
      type(ifl_t), intent(inout) :: f
      integer, intent(in) :: n
      integer :: i, skips
      
      character :: c
      logical :: cmt, pcws, ccws ! comment; previous cchar was white space; current char is white space
      
      cmt  = .false.
      i = f%p
      pcws = i == 1
      skips = 0
      do while (i < f % n .and. skips < n)
         i = i + 1
         c = f%s(i:i)
         
         if (is_comment(c)) then
            cmt = .true.
            cycle
         else if (c == c_new_line) then
            if (cmt) cmt = .false.
            cycle
         end if
         
         ccws = is_blank(c)
         if ((.not. cmt) .and. pcws .and. (.not. ccws)) then
            skips = skips + 1
         end if
         
         pcws = ccws
      end do
      
      if (skips < n) i = -1
      
   end function next_wi
   
   function next_uncmt(f,w,ini) result(i)
      use iso_c_binding, only : c_new_line
      type(ifl_t), intent(inout) :: f
      character(len=*), intent(in) :: w
      integer, intent(in), optional :: ini
      integer :: i,j,n,o
      logical :: c
      character :: predi,sled
!       Return index of the next occurence of word w or 0 if not found
!       This is horribly messy but right now i cannot improve it, so i am just getting it to work somehow
      i = 0
      n = len(trim(w))
      if (n == 0) return

      if (present(ini)) then
         i = ini
      else
         i = f%p
      end if
      
      c = .false.
      outer: do
         o = index(f%s(i:),w)
         i = i + o
!          print *, 'i,o', i, o
         
         if (o == 0) exit ! not found
         
         predi = f%s(i-2:i-2)
         sled = f%s(i+n-1:i+n-1)
        
         if (.not. (is_blank(predi) .and. is_blank(sled))) then
            print *, 'cycling', 's:'//f%s(i-2:i+n-1)//':e'
            cycle
         end if
         
         if (o == 1) exit ! no place for comment so this must be it
         
         if (o >= 2) then
            do j = i-1, 1, -1 ! make sure it is not commented out
               c = is_comment(f%s(j:j))
   !             print *, 'j,c', j, c
               if (c) exit
               if (f%s(j:j) == c_new_line) exit outer
            end do
         end if
         
         if (.not. c) exit
      end do outer      
      
      i = i-1
      if (o /= 0) then
         i = i - 1
      else
         i = 0
      end if
      
   end function next_uncmt
 
   function goto_next_uncmt(f,w) result (exst)
      type(ifl_t), intent(inout) :: f
      character(len=*), intent(in) :: w
      logical :: exst
      integer :: i
      exst = .false.
      i = next_uncmt(f,w)
      if (i == 0) i = next_uncmt(f,w,1)
      if (i /= 0) then
         f % p = i
         exst = .true.
!          print *, 'i,ex', i, exst
      end if      
   end function goto_next_uncmt
 
   subroutine goto_next_wi(f,n)
      type(ifl_t), intent(inout) :: f
      integer, intent(in) :: n
      f%p = next_wi(f,n)
   end subroutine goto_next_wi

   subroutine handle_rd_err(ierr,f)
      integer, intent(in) :: ierr
      type(ifl_t), intent(in) :: f
      
      if (ierr>0) then 
         write(ferr,"('Error ocured during reading of: ',a)") f%s(f%p:f%p+16)
      else if (ierr==-1) then
         write(ferr,"('End of file encountered during reading')")
      end if
      stop 255
      
   end subroutine handle_rd_err
 
   subroutine read_i1(f, r)
      type(ifl_t), intent(inout) :: f
      integer(1), intent(out) :: r
      integer :: ierr
      call goto_next_wi(f,1)
      read(f%s(f%p:),*,iostat=ierr) r
      if (ierr/=0) call handle_rd_err(ierr,f)
#ifdef dbg
      print *, 'read_i1'
#endif
   end subroutine read_i1     
   
   subroutine read_i2(f, r)
      type(ifl_t), intent(inout) :: f
      integer(2), intent(out) :: r
      integer :: ierr
      call goto_next_wi(f,1)
      read(f%s(f%p:),*,iostat=ierr) r
      if (ierr/=0) call handle_rd_err(ierr,f)
#ifdef dbg      
      print *, 'read_i2'
#endif
   end subroutine read_i2

   subroutine read_i4(f, r)
      type(ifl_t), intent(inout) :: f
      integer(4), intent(out) :: r
      integer :: ierr
      call goto_next_wi(f,1)
      read(f%s(f%p:),*,iostat=ierr) r
      if (ierr/=0) call handle_rd_err(ierr,f)
#ifdef dbg
      print *, 'read_i4'
#endif
   end subroutine read_i4
         
   subroutine read_i8(f, r)
      type(ifl_t), intent(inout) :: f
      integer(8), intent(out) :: r
      integer :: ierr
      call goto_next_wi(f,1)
      read(f%s(f%p:),*,iostat=ierr) r
      if (ierr/=0) call handle_rd_err(ierr,f)
#ifdef dbg
      print *, 'read_i8'
#endif
   end subroutine read_i8
   
   
   
   subroutine read_r4(f, r)
      type(ifl_t), intent(inout) :: f
      real(4), intent(out) :: r
      integer :: ierr
      call goto_next_wi(f,1)
      read(f%s(f%p:),*,iostat=ierr) r
      if (ierr/=0) call handle_rd_err(ierr,f)
#ifdef dbg
      print *, 'read_r4'
#endif
   end subroutine read_r4
   
   subroutine read_r8(f, r)
      type(ifl_t), intent(inout) :: f
      real(8), intent(out) :: r
      integer :: ierr
      call goto_next_wi(f,1)
      read(f%s(f%p:f%n),*,iostat=ierr) r
      if (ierr/=0) call handle_rd_err(ierr,f)
#ifdef dbg
      print *, 'read_r8'
#endif
   end subroutine read_r8

#ifndef no_quad   
   subroutine read_r16(f, r)
      type(ifl_t), intent(inout) :: f
      real(16), intent(out) :: r
      integer :: ierr
      call goto_next_wi(f,1)
      read(f%s(f%p:),*,iostat=ierr) r
      if (ierr/=0) call handle_rd_err(ierr,f)
#ifdef dbg
      print *, 'read_r16'
#endif
   end subroutine read_r16   
#endif

   subroutine read_z4(f, r)
      type(ifl_t), intent(inout) :: f
      complex(4), intent(out) :: r
      integer :: ierr
      call goto_next_wi(f,2)
      read(f%s(f%p:),*,iostat=ierr) r
      if (ierr/=0) call handle_rd_err(ierr,f)
#ifdef dbg
      print *, 'read_z4'
#endif
   end subroutine read_z4

   subroutine read_z8(f, r)
      type(ifl_t), intent(inout) :: f
      complex(8), intent(out) :: r
      integer :: ierr
      call goto_next_wi(f,2)
      read(f%s(f%p:),*,iostat=ierr) r
      if (ierr/=0) call handle_rd_err(ierr,f)
#ifdef dbg
      print *, 'read_z8'
#endif
   end subroutine read_z8

#ifndef no_quad    
   subroutine read_z16(f, r)
      type(ifl_t), intent(inout) :: f
      complex(16), intent(out) :: r
      integer :: ierr
      call goto_next_wi(f,2)
      read(f%s(f%p:),*,iostat=ierr) r
      if (ierr/=0) call handle_rd_err(ierr,f)
#ifdef dbg
      print *, 'read_z16'
#endif
   end subroutine read_z16
#endif   

   subroutine read_c(f, r)
      type(ifl_t), intent(inout) :: f
      character(len=*), intent(out) :: r
      integer :: ierr
      call goto_next_wi(f,1)
      read(f%s(f%p:),*,iostat=ierr) r
      if (ierr/=0) call handle_rd_err(ierr,f)
#ifdef dbg
      print *, 'read_c'
#endif
   end subroutine read_c
   
   subroutine read_l(f, r)
      type(ifl_t), intent(inout) :: f
      logical, intent(out) :: r
      integer :: ierr
      call goto_next_wi(f,1)
      read(f%s(f%p:),*,iostat=ierr) r
      if (ierr/=0) call handle_rd_err(ierr,f)
#ifdef dbg
      print *, 'read_l'
#endif
   end subroutine read_l

   
  
   subroutine read_i1a(f, r)
      type(ifl_t), intent(inout) :: f
      integer(1), intent(out) :: r(1:)
      integer :: i
      do i = 1, size(r)
         call read_m(f,r(i))
      end do
#ifdef dbg
      print *, 'read_i1a'
#endif
   end subroutine read_i1a     
   
   subroutine read_i2a(f, r)
      type(ifl_t), intent(inout) :: f
      integer(2), intent(out) :: r(1:)
      integer :: i
      do i = 1, size(r)
         call read_m(f,r(i))
      end do
#ifdef dbg      
      print *, 'read_i2a'
#endif
   end subroutine read_i2a

   subroutine read_i4a(f, r)
      type(ifl_t), intent(inout) :: f
      integer(4), intent(out) :: r(1:)
      integer :: i
      do i = 1, size(r)
         call read_m(f,r(i))
      end do      
#ifdef dbg
      print *, 'read_i4a'
#endif
   end subroutine read_i4a
         
   subroutine read_i8a(f, r)
      type(ifl_t), intent(inout) :: f
      integer(8), intent(out) :: r(1:)
      integer :: i
      do i = 1, size(r)
         call read_m(f,r(i))
      end do
#ifdef dbg
      print *, 'read_i8a'
#endif
   end subroutine read_i8a
   
   
   subroutine read_r4a(f, r)
      type(ifl_t), intent(inout) :: f
      real(4), intent(out) :: r(1:)
      integer :: i
      do i = 1, size(r)
         call read_m(f,r(i))
      end do      
#ifdef dbg
      print *, 'read_r4a'
#endif
   end subroutine read_r4a
   
   subroutine read_r8a(f, r)
      type(ifl_t), intent(inout) :: f
      real(8), intent(out) :: r(1:)
      integer :: i
      do i = 1, size(r)
         call read_m(f,r(i))
      end do
#ifdef dbg
      print *, 'read_r8a'
#endif
   end subroutine read_r8a

#ifndef no_quad       
   subroutine read_r16a(f, r)
      type(ifl_t), intent(inout) :: f
      real(16), intent(out) :: r(1:)
      integer :: i
      do i = 1, size(r)
         call read_m(f,r(i))
      end do      
#ifdef dbg
      print *, 'read_r16a'
#endif
   end subroutine read_r16a
#endif   

   subroutine read_z4a(f, r)
      type(ifl_t), intent(inout) :: f
      complex(4), intent(out) :: r(1:)
      integer :: i
      do i = 1, size(r)
         call read_m(f,r(i))
      end do
#ifdef dbg
      print *, 'read_z4a'
#endif
   end subroutine read_z4a

   subroutine read_z8a(f, r)
      type(ifl_t), intent(inout) :: f
      complex(8), intent(out) :: r(1:)
      integer :: i
      do i = 1, size(r)
         call read_m(f,r(i))
      end do      
#ifdef dbg
      print *, 'read_z8a'
#endif
   end subroutine read_z8a

#ifndef no_quad    
   subroutine read_z16a(f, r)
      type(ifl_t), intent(inout) :: f
      complex(16), intent(out) :: r(1:)
      integer :: i
      do i = 1, size(r)
         call read_m(f,r(i))
      end do
#ifdef dbg
      print *, 'read_z16a'
#endif
   end subroutine read_z16a
#endif   

   subroutine read_ca(f, r)
      type(ifl_t), intent(inout) :: f
      character(len=*), intent(out) :: r(1:)
      integer :: i
      do i = 1, size(r)
         call read_m(f,r(i))
      end do      
#ifdef dbg
      print *, 'read_ca'
#endif
   end subroutine read_ca
   
   
   subroutine read_la(f, r)
      type(ifl_t), intent(inout) :: f
      logical, intent(out) :: r(1:)
      integer :: i
      do i = 1, size(r)
         call read_m(f,r(i))
      end do
#ifdef dbg
      print *, 'read_la'
#endif
   end subroutine read_la
   
   
   
   subroutine read_r8a2(f, r)
      type(ifl_t), intent(inout) :: f
      real(8), intent(out) :: r(:,:)
      integer :: i, j
      do j = 1, size(r,2) 
         do i = 1, size(r,1)
            call read_m(f,r(i,j))
         end do
      end do   
#ifdef dbg
      print *, 'read_r8a'
#endif
   end subroutine read_r8a2
   
   

   subroutine print_iv(v,u,l,f,s)

      integer, intent(in) :: v(1:)
      integer,intent(in) :: u
      character(len=*), intent(in) :: l
      character(len=*), intent(in), optional :: f, s
      character(len=10) :: fm, sm
      integer :: i
      fm = 'i0'
      sm = ''
      if (present(f)) fm = f
      if (present(s)) sm = s
      
      write(u,'(a,a)',advance='no') l, trim(sm)

      do i = 1, size(v)
         write(u,'(x,'//trim(fm)//')',advance='no') v(i)
      end do
      write(u,'()')

   end subroutine print_iv

   subroutine print_lv(v,u,l,f,s)
      logical, intent(in) :: v(:)
      integer,intent(in) :: u
      character(len=*), intent(in) :: l
      character(len=*), intent(in),optional :: f, s
      character(len=10) :: fm, sm
      sm = ''
      if (present(s)) sm = s
      write(u,'(a,a)',advance='no') l, trim(sm)
      write(u, *) v
      
   end subroutine  print_lv
   


   subroutine print_v(v,u,l,f,s)

      real(dp), intent(in) :: v(1:)
      integer,intent(in) :: u
      character(len=*), intent(in) :: l
      character(len=*), intent(in), optional :: f, s
      character(len=10) :: fm, sm
      integer :: i
      fm = 'f16.10'
      sm = ''
      if (present(f)) fm = f
      if (present(s)) sm = s
      
      write(u,'(a,a)',advance='no') l, trim(sm)

      do i = 1, size(v)
         write(u,'(x,'//trim(fm)//')',advance='no') v(i)
      end do
      write(u,'()')

   end subroutine print_v


   subroutine print_mr(a,u,l,f,s)

      real(dp), intent(in) :: a(1:,1:)
      integer,intent(in) :: u
      character(len=*), intent(in) :: l
      character(len=*), intent(in), optional :: f, s
      character(len=10) :: fm, sm
      integer :: i,j
      
      fm = 'f16.10'
      sm = ''
      if (present(f)) fm = f
      if (present(s)) sm = s
      
      write(u,'(a,a)',advance='no') l, trim(sm)

      do i = 1, size(a,1)
         do j = 1, size(a,2)
            write(u,'(x,'//trim(fm)//')',advance='no'), a(i,j)
         end do
         write(u,'()')
      end do

   end subroutine print_mr

   subroutine print_mz(a,u,l,f,s)

      complex(dp), intent(in) :: a(1:,1:)
      integer,intent(in) :: u
      character(len=*), intent(in) :: l
      character(len=*), intent(in), optional :: f, s
      character(len=10) :: fm, sm            
      integer :: i,j
      
      fm = 'f16.10'
      sm = ''
      if (present(f)) fm = f
      if (present(s)) sm = s


      write(u,'(a,a)',advance='no') l, trim(sm)

      do i = 1, size(a,1)
         do j = 1, size(a,2)
            write(u,'(x,("(",'//trim(fm)//',",",'//trim(fm)//',")"))',advance='no'), a(i,j)
         end do
         write(u,'()', advance='yes')
      end do

   end subroutine print_mz

        

      
! The following routines are old and possibly not clean of bugs      
      

   subroutine get_var_la(var,name,u)
      logical, intent(out) :: var(:)
      character(len=*), intent(in) :: name
      integer, intent(in) :: u
      character(len=20) :: sname
      integer :: i,r ,stat, cp
      character(len=300) :: s
      
      rewind(u)
      do 
         read(u,"(a)",iostat=stat) s
         if (stat == 0) then
            cp = index(s,'#')
            if (cp > 0) s = trim(s(1:cp-1))
!                 print *,s
            if (len(trim(s)) < len(trim(name))+1) cycle
            i = index(s, trim(name)//' ', back=.true.)
            if ((i == 1) .or. (s(max(i-1,1):max(i-1,1)) == ' ')) then
                  
               read(s(i+len(trim(name))+1:), fmt=*,  iostat=stat) var 
               if (stat == 0) then
!                         print *,var
               else if (stat < 0) then
!                         print *,'cycling'
                  cycle
               else if (stat > 0) then
                  write(0,'("conversion error")')
                  exit
               end if
                  
            end if
            
         else if (stat < 0 ) then
            exit
         else if (stat > 0) then
            write(0,'("error reading the file")')
            exit
         end if
      end do
      
         
   end subroutine get_var_la
   
   
   
   
   subroutine get_var_l(var,name,u)
      logical, intent(out) :: var
      character(len=*), intent(in) :: name
      integer, intent(in) :: u
      character(len=20) :: sname
      integer :: i,r ,stat, cp
      character(len=300) :: s
      
      rewind(u)
      do 
         read(u,"(a)",iostat=stat) s
         if (stat == 0) then
            cp = index(s,'#')
            if (cp > 0) s = trim(s(1:cp-1))
!                 print *,s
            if (len(trim(s)) < len(trim(name))+1) cycle
            i = index(s, trim(name)//' ', back=.true.)
            if ((i == 1) .or. (s(max(i-1,1):max(i-1,1)) == ' ')) then
               
               read(s(i+len(trim(name))+1:), fmt=*,  iostat=stat) var 
               if (stat == 0) then
!                         print *,var
               else if (stat < 0) then
!                         print *,'cycling'
                  cycle
               else if (stat > 0) then
                  write(0,'("conversion error")')
                  exit
               end if
               
            end if
            
            
            
         else if (stat < 0 ) then
            exit
         else if (stat > 0) then
            write(0,'("error reading the file")')
            exit
         end if
      end do
      
      
   end subroutine get_var_l
   


   subroutine get_var_aa(var,name,u)
      character(len=*), intent(out) :: var(:)
      character(len=*), intent(in) :: name
      integer, intent(in) :: u
      character(len=20) :: sname
      integer :: i,r ,stat, cp
      character(len=300) :: s
      
      rewind(u)
      do 
         read(u,"(a)",iostat=stat) s
         if (stat == 0) then
            cp = index(s,'#')
            if (cp > 0) s = trim(s(1:cp-1))
!                 print *,s
            if (len(trim(s)) < len(trim(name))+1) cycle
            i = index(s, trim(name)//' ', back=.true.)
            if ((i == 1) .or. (s(max(i-1,1):max(i-1,1)) == ' ')) then
                  
               read(s(i+len(trim(name))+1:), fmt=*,  iostat=stat) var 
               if (stat == 0) then
!                         print *,var
               else if (stat < 0) then
!                         print *,'cycling'
                  cycle
               else if (stat > 0) then
                  write(0,'("conversion error")')
                  exit
               end if
               
            end if
            
            
            
         else if (stat < 0 ) then
            exit
         else if (stat > 0) then
            write(0,'("error reading the file")')
            exit
         end if
      end do
      
         
   end subroutine get_var_aa
   
   
      
      
   subroutine get_var_a(var,name,u)
      character(len=*), intent(out) :: var
      character(len=*), intent(in) :: name
      integer, intent(in) :: u
      character(len=20) :: sname
      integer :: i,r ,stat, cp
      character(len=300) :: s
      
      rewind(u)
      do 
         read(u,"(a)",iostat=stat) s
         if (stat == 0) then
            cp = index(s,'#')
            if (cp > 0) s = trim(s(1:cp-1))
!                 print *,s
            if (len(trim(s)) < len(trim(name))+1) cycle
            i = index(s, trim(name)//' ', back=.true.)
            if ((i == 1) .or. (s(max(i-1,1):max(i-1,1)) == ' ')) then
                  
               read(s(i+len(trim(name))+1:), fmt=*,  iostat=stat) var 
               if (stat == 0) then
!                         print *,var
               else if (stat < 0) then
!                         print *,'cycling'
                  cycle
               else if (stat > 0) then
                  write(0,'("conversion error")')
                  exit
               end if
                  
            end if
            
            
            
         else if (stat < 0 ) then
            exit
         else if (stat > 0) then
            write(0,'("error reading the file")')
            exit
         end if
      end do
      
         
   end subroutine get_var_a
   
   
   subroutine get_var_ra(var,name,u)
      real(dp), intent(out) :: var(:)
      character(len=*), intent(in) :: name
      integer, intent(in) :: u
      character(len=20) :: sname
      integer :: i,r ,stat, cp
      character(len=300) :: s
      
      rewind(u)
      do 
         read(u,"(a)",iostat=stat) s
         if (stat == 0) then
!                 print *,s
            cp = index(s,'#')
            if (cp > 0) s = trim(s(1:cp-1))
!                 print *,s
            if (len(trim(s)) < len(trim(name))+1) cycle
            i = index(s, trim(name)//' ', back=.true.)
            if ((i == 1) .or. (s(max(i-1,1):max(i-1,1)) == ' ')) then
                  
                  read(s(i+len(trim(name))+1:), fmt=*,  iostat=stat) var 
                  if (stat == 0) then
!                         print *,var
                  else if (stat < 0) then
!                         print *,'cycling'
                     cycle
                  else if (stat > 0) then
                     write(0,'("conversion error")')
                     exit
                  end if
                  
            end if
            
            
            
         else if (stat < 0 ) then
            exit
         else if (stat > 0) then
            write(0,'("error reading the file")')
            exit
         end if
      end do
      
      
   end subroutine get_var_ra
   
   subroutine get_var_ia(var,name,u)
      integer, intent(out) :: var(:)
      character(len=*), intent(in) :: name
      integer, intent(in) :: u
      character(len=20) :: sname
      integer :: i,r ,stat, cp
      character(len=300) :: s
      
      rewind(u)
      do 
         read(u,"(a)",iostat=stat) s
         if (stat == 0) then
            cp = index(s,'#')
            if (cp > 0) s = trim(s(1:cp-1))
!                 print *,s
            if (len(trim(s)) < len(trim(name))+1) cycle
            i = index(s, trim(name)//' ', back=.true.)
            if ((i == 1) .or. (s(max(i-1,1):max(i-1,1)) == ' ')) then
               
               read(s(i+len(trim(name))+1:), fmt=*,  iostat=stat) var 
               if (stat == 0) then
!                         print *,var
               else if (stat < 0) then
!                         print *,'cycling'
                  cycle
               else if (stat > 0) then
                  write(0,'("conversion error")')
                  exit
               end if
               
            end if
            
            
            
         else if (stat < 0 ) then
            exit
         else if (stat > 0) then
            write(0,'("error reading the file")')
            exit
         end if
      end do
      
      
   end subroutine get_var_ia
   
      
      
   subroutine get_var_r(var,name,u)
      real(dp), intent(out) :: var
      character(len=*), intent(in) :: name
      integer, intent(in) :: u
      character(len=20) :: sname
      integer :: i,r ,stat, cp
      character(len=300) :: s
      
      rewind(u)
      do 
         read(u,"(a)",iostat=stat) s
         if (stat == 0) then
            cp = index(s,'#')
            if (cp > 0) s = trim(s(1:cp-1))
!                 print *,s
            if (len(trim(s)) < len(trim(name))+1) cycle
            i = index(s, trim(name)//' ', back=.true.)
            if (i<1) cycle
            if ((i == 1) .or. (s(max(i-1,1):max(i-1,1)) == ' ')) then
               
               read(s(i+len(trim(name))+1:), fmt=*,  iostat=stat) var 
               if (stat == 0) then
!                         print *,var
               else if (stat < 0) then
!                         print *,'cycling'
                  cycle
               else if (stat > 0) then
                  write(0,'("conversion error")')
                  exit
               end if
                  
            end if
            
            
            
         else if (stat < 0 ) then
            exit
         else if (stat > 0) then
            write(0,'("error reading the file")')
            exit
         end if
      end do
      
      
   end subroutine get_var_r
      
   subroutine get_var_i(var,name,u)
      integer, intent(out) :: var
      character(len=*), intent(in) :: name
      integer, intent(in) :: u
      character(len=20) :: sname
      integer :: i,r ,stat, cp
      character(len=300) :: s
      
      rewind(u)
      do 
         read(u,"(a)",iostat=stat) s
         if (stat == 0) then
            cp = index(s,'#')
            if (cp > 0) s = trim(s(1:cp-1))
!                 print *,s
            if (len(trim(s)) < len(trim(name))+1) cycle
            i = index(s, trim(name)//' ', back=.true.)
            if ((i == 1) .or. (s(max(i-1,1):max(i-1,1)) == ' ')) then
                  
               read(s(i+len(trim(name))+1:), fmt=*,  iostat=stat) var 
               if (stat == 0) then
!                         print *,var
               else if (stat < 0) then
!                         print *,'cycling'
                  cycle
               else if (stat > 0) then
                  write(0,'("conversion error")')
                  exit
               end if
                  
            end if
         else if (stat < 0 ) then
            exit
         else if (stat > 0) then
            write(0,'("error reading the file")')
            exit
         end if
      end do
   
   
   end subroutine get_var_i



   end module mod_io

