    module mod_misc
      
      use mod_precision
      use mod_io
      implicit none
      
      
      interface get_urandom
        module procedure get_urandom_s, get_urandom_v, get_urandom_m
      end interface get_urandom
      
      contains
      
      subroutine get_stamp(stamp)
        character(len=*), intent(out) :: stamp
        character(len=8) :: stamp_date
        character(len=10) :: stamp_time
        
        call date_and_time(date=stamp_date, time=stamp_time)
      
        stamp = stamp_date(7:8)//'.'//stamp_date(5:6)//'.'//stamp_date(1:4)//'-'// &
                  & stamp_time(1:2)//'.'//stamp_time(3:4)//'.'//stamp_time(5:)
                  
      end subroutine get_stamp
      
      
      subroutine get_urandom_s(a)
        integer, intent(out) :: a
        integer :: u
        u = get_new_unit()
        open(unit=u, file='/dev/urandom', access='stream', action='read')
        read(unit=u) a
        call close_unit(u)
      end subroutine get_urandom_s

      subroutine get_urandom_v(a)
        integer, intent(out) :: a(:)
        integer :: u
        u = get_new_unit()
        open(unit=u, file='/dev/urandom', access='stream', action='read')
        read(unit=u) a
        call close_unit(u)
      end subroutine get_urandom_v
      
      subroutine get_urandom_m(a)
        integer, intent(out) :: a(:,:)
        integer :: u
        u = get_new_unit()
        open(unit=u, file='/dev/urandom', access='stream', action='read')
        read(unit=u) a
        call close_unit(u)
      end subroutine get_urandom_m



      function in_argv(a,p)
         character(len=*), intent(in) :: a
         integer, intent(in), optional :: p
         logical :: in_argv
         integer :: argc, i
         character(len=100) :: iarg

         in_argv = .false.

         argc = command_argument_count()

         if (argc > 0) then 
            if (present(p)) then
               if (argc >= p) then
                  call get_command_argument(p,iarg)
                  if (iarg(1:len(a)) == a) in_argv = .true.
               end if
            else
               do i=1,argc
                  call get_command_argument(i,iarg)
                  if (iarg(1:len(a)) == a) then
                     in_argv = .true.
                     exit
                  endif
               end do
            end if
         end if

      end function in_argv

    end module mod_misc