 
      subroutine addg3(ad,totnd,bins,nbins,aptr,nd,bptr,mxnb,rmax)
          use mod_precision


!
!    This procedure updates the g3(theta) summation.
!

      implicit none

!
!    Declare the simple variables.
!

      real(dp) :: mrab,mrac,rmax
      real(dp) :: costht,theta

      integer ind
      integer mxnb
      integer ia,ib,ic
      integer ja,jb
      integer nd,nbins,totnd
      integer eol

!
!    Define parameters.
!

      parameter(eol = -1000000)

!
!    Declare the arrays.
!

      real(dp) :: rab(3)
      real(dp) :: rac(3)
      real(dp) :: ad(3,totnd)

      integer bins(nbins)
      integer aptr(nd)
      integer bptr(mxnb)

!
!    Update the histogram.
!

      do ia = 1,nd,1
         ja = aptr(ia)
         do while (bptr(ja) /= eol)
            ib = bptr(ja)
            if (ia /= ib) then
               rab(1) = ad(1,ib) - ad(1,ia)
               rab(2) = ad(2,ib) - ad(2,ia)
               rab(3) = ad(3,ib) - ad(3,ia)
               mrab = sqrt(rab(1)**2+rab(2)**2+rab(3)**2)
               if (mrab <= rmax) then
                  rab(1) = rab(1)/mrab
                  rab(2) = rab(2)/mrab
                  rab(3) = rab(3)/mrab
                  jb = ja+1
                  do while (bptr(jb) /= eol)
                     ic = bptr(jb)
                     rac(1) = ad(1,ic) - ad(1,ia)
                     rac(2) = ad(2,ic) - ad(2,ia)
                     rac(3) = ad(3,ic) - ad(3,ia)
                     mrac = sqrt(rac(1)**2+rac(2)**2+rac(3)**2)
                     if (mrac <= rmax) then
                        rac(1) = rac(1)/mrac
                        rac(2) = rac(2)/mrac
                        rac(3) = rac(3)/mrac
                        costht = rab(1)*rac(1)+rab(2)*rac(2) & 
     &                           +rab(3)*rac(3)
                        costht = max(min(1.0_dp,costht),-1.0_dp)
                        theta = acos(costht)* & 
     &                          180.0_dp/3.1415926535_dp
                        ind = int(theta*real(nbins, dp)/180.0_dp)+1
                        if ((ind >= 1).and.(ind <= nbins)) & 
     &                      bins(ind) = bins(ind) + 1
                     endif
                     jb = jb + 1
                  enddo
               endif
            endif
            ja = ja + 1
         enddo
      enddo

      end
