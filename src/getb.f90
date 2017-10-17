 
       function getb(ainf)
          use mod_precision


!
!    This is the function to evaluate the upper and lower eigenvalues of
!     the Beer-Pettifor-Aoki matrix, and to return a value equal to the
!     difference in magnitudes between the upper and lower limits.
!

      implicit none
      real(dp) :: getb

!
!    Declare the simple variables.
!

      real(dp) :: ainf
      real(dp) :: bupper,blower
      real(dp), parameter :: sqrt2 = 1.4142135623730950_dp !1.4142136_dp !

      integer locmrec,nrec
      integer i,ifail

!
!    Define the parameters.
!

      parameter (locmrec = 100)
!       parameter (sqrt2 = 1.4142136_dp)
      

!
!    Declare the arrays.
!

      real(dp) :: a(0:locmrec)
      real(dp) :: b(0:locmrec+1)
      real(dp) :: d(locmrec+2)
      real(dp) :: e(locmrec+2)

!
!    Declare the common blocks.
!

      common /params/a,b,nrec
      common /term/bupper,blower

!
!    Set up the Beer-Pettifor-Aoki matrix.
!

      if (nrec == 0) then

         ainf = a(0)
         bupper = 0.0_dp
         blower = -bupper
         getb = 0.0_dp

      else

         if (nrec == 1) then

            d(1) = a(0) - ainf
            d(2) = a(1) - ainf
            e(2) = b(1)/sqrt2

         else

            d(1) = a(0) - ainf
            e(2) = b(1)/sqrt2
            d(2) = 0.5_dp*(a(1)-ainf)
            
            do i = 2,nrec-1,1
               d(i+1) = 0.5_dp*(a(i)-ainf)
               e(i+1) = 0.5_dp*b(i)
            enddo
            
            d(nrec+1) = a(nrec) - ainf
            e(nrec+1) = b(nrec)/sqrt2

         endif

!
!       Diagonalise the matrix.
!

!        CALL F02AVF(NREC+1,1.0D-15,D,E,IFAIL)
         call tql1(nrec+1,d,e,ifail)
         if (ifail /= 0) then
            write(6,'('' Unable to diagonalise matrix.'')')
            call panic()
         endif

!
!       Find the difference between the magnitudes of the limiting
!        eigenvalues.
!

         blower = d(1)
         bupper = d(nrec+1)

         getb = bupper - abs(blower)

      endif

      end

