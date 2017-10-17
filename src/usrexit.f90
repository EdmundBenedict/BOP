 
      integer function usrexit()
          use mod_precision


!
!    This is a function that tests to see if the usr has requested the
!     smooth termination of the calling program.
!    A value of 1 is returned if the file "exit.d" exists, otherwise
!     a value of 0 is returned.
!

      implicit none

!
!    Declare the simple variables.
!

      logical ex

!
!    Test for the existence of "exit.d".
!

      inquire(file = 'BOP.exit', exist = ex)

      if (ex) then
         write(6,'(''User request to exit program found.'')')
         usrexit = 1
         open(unit=51,file='BOP.exit',status='OLD')
         close(51,status='DELETE')
      else
         usrexit = 0
      endif

      end

