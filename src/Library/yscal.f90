 
      subroutine  yscal(n,zar,zai,zxr,zxi,incx)
          use mod_precision

!
!     scales a complex vector by a constant, using real arithmetic
!     jack dongarra, 3/11/78.
!
      real(dp) :: zar,zai,zxr(1),zxi(1),tmp
      if( n <= 0 .or. incx <= 0 )return
      if(incx == 1)go to 20
!
!        code for increment not equal to 1
!
      ix = 1
!     if(incx.lt.0)ix = (-n+1)*incx + 1
      do 10 i = 1,n
        tmp      = zar*zxr(ix) - zai*zxi(ix)
        zxi(ix) = zar*zxi(ix) + zai*zxr(ix)
        zxr(ix) = tmp
        ix = ix + incx
   10 continue
      return
!
!        code for increment equal to 1
!
   20 do 30 i = 1,n
        tmp    = zar*zxr(i) - zai*zxi(i)
        zxi(i) = zar*zxi(i) + zai*zxr(i)
        zxr(i) = tmp
   30 continue
      return
      end
