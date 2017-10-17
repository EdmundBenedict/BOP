 
      subroutine rdtokn(record,reclen,tokbgn,tokend,error)

!
!    This is a routine to find the beginning and end of a token in a
!     record. The separators between tokens are taken to be spaces.
!

      implicit none

!
!    Declare the simple variables.
!

      integer reclen,tokbgn,tokend,error

!
!    Declare the arrays.
!

      character record(reclen)

!
!    Pass over leading spaces.
!

 1    if (tokbgn <= reclen) then
         if (record(tokbgn) == ' ') then
            tokbgn = tokbgn + 1
            goto 1
         endif
      endif

      if (tokbgn > reclen) then
         write(6,'(''RDTOKN: No token found.'')')
         error = 1
         return
      endif

!
!    Find the end of the token.
!

      tokend = tokbgn + 1
 2    if (tokend <= reclen) then
         if (record(tokend) /= ' ') then
            tokend = tokend + 1
            goto 2
         endif
      endif
      tokend = tokend - 1

      error = 0

      end

