 
      subroutine dswap (n,dx,incx,dy,incy)
          use mod_precision

!
!     interchanges two vectors. Adapted from:
!     jack dongarra, linpack, 3/11/78.
!
      real(dp) :: dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,n
!
!#ifdefC APOLLO
!      if (n .le. 0) return
!      if (incx .eq. 1 .and. incy .eq. 1) then
!        call vec_$dswap(dx,dy,n)
!      else
!        call vec_$dswap_i(dx,incx,dy,incy,n)
!      endif
!#else
      ix = 1
      iy = 1
      if (incx < 0) ix = (1-n)*incx + 1
      if (incy < 0) iy = (1-n)*incy + 1
      do 10 i = 1,n
        dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 continue
!#endif
      end
