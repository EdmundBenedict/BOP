 
      integer function iyamax(n,zxr,zxi,incx)
          use mod_precision

!     implicit none
!
!     finds the index of element having max. absolute value.
!     Adapted from blas, but avoids separates real, imaginary parts
!     jack dongarra, 1/15/85.
!     modified 3/93 to return if incx .le. 0.
!
      real(dp) :: zxr(1),zxi(1)
      real(dp) :: smax
      integer i,incx,ix,n
      real(dp) :: dcabs1
!
      iyamax = 0
      if( n < 1 .or. incx <= 0 )return
      iyamax = 1
      if(n == 1)return
      if(incx == 1)go to 20
!
!        code for increment not equal to 1
!
      ix = 1
      smax = abs(zxr(1)) + abs(zxi(1))
      ix = ix + incx
      do 10 i = 2,n
         if(abs(zxr(ix))+abs(zxi(ix)) <= smax) go to 5
         iyamax = i
         smax = abs(zxr(ix))+abs(zxi(ix))
    5    ix = ix + incx
   10 continue
      return
!
!        code for increment equal to 1
!
   20 smax = abs(zxr(1)) + abs(zxi(1))
      do 30 i = 2,n
         if(abs(zxr(i))+abs(zxi(i)) <= smax) go to 30
         iyamax = i
         smax = abs(zxr(i))+abs(zxi(i))
   30 continue
      return
      end
