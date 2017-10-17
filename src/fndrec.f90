 

      subroutine fndrec(unt,str,stat)
          use mod_precision


!
!    This is a routine to find the position of a record in a file.
!

!
!    Declare the simple variables.
!

      integer, intent(in) :: unt
      integer, intent(inout) :: stat
      integer strlen

      character(len=*), intent(in) :: str
      character*100 instr

!
!    Scan though the file looking for the required string.
!

      strlen = len(str)

      rewind(unt)

 2    read(unt,'(A)',end = 1) instr

      if (instr(1:strlen) == str) then
         stat = 0
!         BACKSPACE(UNT)
         goto 3
      else
         goto 2
      endif

 1    write(6,'('' Unable to find record '',A)')str
      stat = 1

 3    end

