 
      subroutine four1(data,nn,isign)
          use mod_precision


!
!    This is a routine to carry out a discrete fast fourier transform.
!    It is taken from Numerical Recipes by Press et al.
!    ISIGN = 1   ==>    Forward transform
!    ISIGN =-1   ==>    Inverse transform
!

      implicit none

!
!    Declare the simple variables.
!

      real(dp) :: wr,wi,wpr,wpi,wtemp,theta
      real(dp) :: tempr,tempi

      integer nn,isign
      integer m,n,i,j
      integer mmax,istep

!
!    Declare the arrays.
!

      real(dp) :: data(2*nn)

!
!    Evaluate the transform.
!

      n = 2*nn
      j = 1
      do i = 1,n,2
         if (j > i) then
            tempr = data(j)
            tempi = data(j+1)
            data(j) = data(i)
            data(j+1) = data(i+1)
            data(i) = tempr
            data(i+1) = tempi
         endif
         m = n/2
         do while ((m >= 2).and.(j > m))
            j = j - m
            m = m/2
         enddo
         j = j + m
      enddo

      mmax = 2
      do while (n > mmax)
         istep = 2*mmax
         theta = 6.28318530717959_dp/real(isign*mmax, dp)
         wpr = -2.0_dp*(sin(0.5_dp*theta)**2)
         wpi = sin(theta)
         wr = 1.0_dp
         wi = 0.0_dp
         do m = 1,mmax,2
            do i = m,n,istep
               j = i + mmax
               tempr = wr*data(j) - wi*data(j+1)
               tempi = wr*data(j+1) + wi*data(j)
               data(j) = data(i) - tempr
               data(j+1) = data(i+1) - tempi
               data(i) = data(i) + tempr
               data(i+1) = data(i+1) + tempi
            enddo
            wtemp = wr
            wr = wr*wpr - wi*wpi + wr
            wi = wi*wpr + wtemp*wpi + wi
         enddo
         mmax = istep
      enddo

      if (isign == -1) then
         do i = 1,2*nn,1
            data(i) = data(i)/real(nn, dp)
         enddo
      endif

      end

