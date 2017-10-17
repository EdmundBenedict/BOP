 
      subroutine states(z,nlz,nsttz,llistz)
          use mod_precision
          use mod_const, only : ml, mxz


!
!    This is a routine that returns the number of angular momentum states
!     used for the tight-binding parameterization of an element.
!

      implicit none

!
!    Declare the simple variables.
!

      integer z,nlz,nsttz

!    Declare the arrays.
!

      integer nstt(0:mxz)
      integer llist(ml,0:mxz)
      integer llistz(ml)
      integer nl(0:mxz)

!
!    Declare the common block.
!

      common/lcom/nl,nstt,llist

!
!    Assign the state parameters.
!

      nlz = nl(z)
      nsttz = nstt(z)
!       llistz(1) = llist(1,z)
!       llistz(2) = llist(2,z)
!       llistz(3) = llist(3,z)
      
      llistz = llist(:,z)

      end

 
  subroutine tbstates(z,nlz,nsttz,llistz)
          use mod_precision
          use mod_const, only : ml, mxz


!
!    This is a routine that returns the number of angular momentum states
!     used for the tight-binding parameterization of an element.
!

      implicit none

!
!    Declare the simple variables.
!

      integer z,nlz,nsttz

!    Declare the arrays.
!

      integer tbnstt(0:mxz)
      integer tbllist(ml,0:mxz)
      integer llistz(ml)
      integer tbnl(0:mxz)

!
!    Declare the common block.
!

      common/lcom/tbnl,tbnstt,tbllist

!
!    Assign the state parameters.
!

      nlz = tbnl(z)
      nsttz = tbnstt(z)
      
      llistz = tbllist(:,z)

      end