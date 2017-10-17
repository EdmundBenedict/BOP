 
      subroutine locabinf(ain,bin,nrecin,ainf,binf)
          use mod_precision


!
!    This is a routine to find the parameters for the square root
!     terminator using the method of Beer, Pettifor and Aoki.
!

      implicit none

!
!    Declare the simple variables.
!

      real(dp) :: ainf,binf
      real(dp) :: astep,alower,aupper
      real(dp) :: bupper,blower
      real(dp) :: erreta,eps,getb

      integer locmrec,nrec,nrecin
      integer i,ifail

!
!    Define the parameters.
!

      parameter (locmrec = 100)
      parameter (eps = 1.0e-4_dp)
      parameter (erreta = 0.0_dp)
      parameter (astep = 0.1_dp)

!
!    Declare the arrays.
!

      real(dp) :: ain(0:nrecin)
      real(dp) :: bin(0:nrecin+1)
      real(dp) :: a(0:locmrec)
      real(dp) :: b(0:locmrec+1)

!
!    Declare the common blocks.
!

      common /params/a,b,nrec
      common /term/bupper,blower

!
!    Declare the external functions.
!

      external getb

!
!    Copy values of A and B over, ready for passing to GETB.
!

      if (nrecin > locmrec) then
         write(6,'(''Too many recursion levels for limits finder.'')')
         write(6,'(''Increase LOCMREC in LOCABINF and GETB.'')')
         call panic()
      endif

      nrec = nrecin

     a(0:nrec) = ain(0:nrec)
     b(0:nrec) = bin(0:nrec)
     b(0) = 0.0_dp
!
!    Estimate starting values of a and b.
!

      ainf = a(nrec)
      binf = b(nrec)

!
!    Find the final values for a and b.
!

!     CALL C05AGF(AINF,ASTEP,EPS,ERRETA,GETB,AUPPER,ALOWER,IFAIL)
      alower = ainf - 0.1_dp
      aupper = ainf + 0.1_dp
!       print *, ainf, binf, alower, aupper
      call findzero(getb,alower,aupper,eps,ainf,ifail)
      if (ifail /= 0) then
         write(6,'('' Unable to find zero of function.'')')
         call panic()
      endif

!
!    Take the final b value to be the mean of the two calculated values.
!

      binf = 0.5_dp*(abs(bupper)+abs(blower))

!      WRITE(6,'(''AINF = '',G22.15,'' BINF = '',G22.15)') AINF,BINF

      end

