 

      integer function findj(bptr,mxnnb,id)
          use mod_precision
          use mod_const, only : eol


!
!    This is a function to find the position in a neighbor list of a particular
!     atom. A return value of EOL means that the atom could not be found.
!

      implicit none

!
!    Declare the simple variables.
!

      integer id,mxnnb

!
!    Declare the arrays.
!

      integer bptr(mxnnb)

!
!    Search for the atom.
!

      findj = 1
 1    if (bptr(findj) == eol) then
         findj = eol
      elseif (bptr(findj) /= id) then
         findj = findj+1
         goto 1
      endif

      end

