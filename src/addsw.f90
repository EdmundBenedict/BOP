 
      subroutine addsw(ad,totnd,nd,bin,nbins,aptr,bptr,mxnb,de)
          use mod_precision


!
!    This is a routine to evaluate the Stillinger Weber 3 body potential.
!

      implicit none

!
!    Declare the simple variables.
!

      real(dp) :: v3b
      real(dp) :: mr12,mr23,mr31
      real(dp) :: g12,g23,g31
      real(dp) :: gamma,efac,rcut
      real(dp) :: ct123,ct312,ct231
      real(dp) :: de

      integer nd,totnd
      integer i1
      integer i1a,i1b,i2a,i2b
      integer nbins,ibin
      integer mxnb,eol

!
!    Define the parameters.
!

      parameter (gamma = 2.5141_dp)
      parameter (rcut = 3.7712_dp)
      parameter (efac = 2.98103_dp)
      parameter (eol = -1000000)

!
!    Declare the arrays.
!

      real(dp) :: ad(3,totnd)
      real(dp) :: r12(3)
      real(dp) :: r23(3)
      real(dp) :: r31(3)

      integer bin(nbins)
      integer aptr(nd)
      integer bptr(mxnb)

!
!    Evaluate the SW potential and add the contributions to
!     the histogram.
!

      do 1 i1 = 1,nd,1

         v3b = 0.0_dp

         i1a = aptr(i1)
 2       i1b = bptr(i1a)

         if (i1b /= eol) then

            if (i1b /= i1) then

               r12(1) = ad(1,i1) - ad(1,i1b)
               r12(2) = ad(2,i1) - ad(2,i1b)
               r12(3) = ad(3,i1) - ad(3,i1b)

               mr12 = sqrt(r12(1)**2+r12(2)**2+r12(3)**2)
               if (mr12 < rcut-0.1_dp) then
                  g12 = exp(gamma/(mr12-rcut))
               else
                  g12 = 0.0_dp
               endif

               i2a = i1a + 1
 3             i2b = bptr(i2a)

               if (i2b /= eol) then

                  if ((i2b /= i1b).and.(i2b /= i1)) then

                     r23(1) = ad(1,i1b) - ad(1,i2b)
                     r23(2) = ad(2,i1b) - ad(2,i2b)
                     r23(3) = ad(3,i1b) - ad(3,i2b)

                     mr23 = sqrt(r23(1)**2+r23(2)**2+r23(3)**2)
                     if (mr23 < rcut-0.1_dp) then
                        g23 = exp(gamma/(mr23-rcut))
                        ct123 = -(r12(1)*r23(1)+r12(2)*r23(2)+ & 
     &                            r12(3)*r23(3))/(mr12*mr23)
                        v3b = v3b + g12*g23*(ct123+1.0_dp/3.0_dp)**2
                     else
                        g23 = 0.0_dp
                     endif

                     r31(1) = ad(1,i2b) - ad(1,i1)
                     r31(2) = ad(2,i2b) - ad(2,i1)
                     r31(3) = ad(3,i2b) - ad(3,i1)

                     mr31 = sqrt(r31(1)**2+r31(2)**2+r31(3)**2)
                     if (mr31 < rcut-0.1_dp) then
                        g31 = exp(gamma/(mr31-rcut))
                        ct312 = -(r31(1)*r12(1)+r31(2)*r12(2)+ & 
     &                            r31(3)*r12(3))/(mr31*mr12)
                        ct231 = -(r23(1)*r31(1)+r23(2)*r31(2)+ & 
     &                            r23(3)*r31(3))/(mr23*mr31)
                        v3b = v3b + g31*g12*(ct312+1.0_dp/3.0_dp)**2 & 
     &                            + g23*g31*(ct231+1.0_dp/3.0_dp)**2
                     endif

                  endif

                  i2a = i2a + 1
                  goto 3

               endif

            endif

            i1a = i1a + 1
            goto 2

         endif

         v3b = v3b*efac

         ibin = int(v3b/de+0.5_dp)+1
         if (ibin < 1) then
            ibin = 1
         elseif (ibin > nbins) then
            ibin = nbins
         endif
         bin(ibin) = bin(ibin) + 1

 1    continue

      end
