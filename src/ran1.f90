 
       function ran1(idum)
          use mod_precision

      implicit double precision  (a-h,o-z)
      real(dp) :: ran1
      dimension r(97)
      parameter (m1=259200,ia1=7141,ic1=54773,rm1=3.8580247e-6)
      parameter (m2=134456,ia2=8121,ic2=28411,rm2=7.4373773e-6)
      parameter (m3=243000,ia3=4561,ic3=51349)
      data iff /0/
      save ix1,ix2,ix3,r
      if (idum < 0.or.iff == 0) then
        iff=1
        ix1=mod(ic1-idum,m1)
        ix1=mod(ia1*ix1+ic1,m1)
        ix2=mod(ix1,m2)
        ix1=mod(ia1*ix1+ic1,m1)
        ix3=mod(ix1,m3)
        do j=1,97
          ix1=mod(ia1*ix1+ic1,m1)
          ix2=mod(ia2*ix2+ic2,m2)
          r(j)=(real(ix1, dp)+real(ix2, dp)*rm2)*rm1
        end do
        idum=1
      endif
      ix1=mod(ia1*ix1+ic1,m1)
      ix2=mod(ia2*ix2+ic2,m2)
      ix3=mod(ia3*ix3+ic3,m3)
      j=1+(97*ix3)/m3
      if(j > 97.or.j < 1) print *,'"pause" was here if it tells you anything.'
      ran1=r(j)
      r(j)=(real(ix1, dp)+real(ix2, dp)*rm2)*rm1
      return
      end

