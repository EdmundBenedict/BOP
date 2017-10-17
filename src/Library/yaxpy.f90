 
!#define HP
      subroutine yaxpy(n,ar,ai,xr,xi,incx,yr,yi,incy,sw)
          use mod_precision

!- Complex daxpy, using real arithmetic
      implicit none
! Passed parameters
      logical sw
      integer n,incx,incy
      real(dp) :: ar,ai,xr(1),xi(1),yr(1),yi(1)
! Local variables:
      integer ix,iy,i

      if (n  <=  0) return

!#ifdefC BLAS
!      call daxpy(n, ar,xr,incx,yr,incy)
!      if (sw) then
!        call daxpy(n,-ai,xi,incx,yr,incy)
!        call daxpy(n, ar,xi,incx,yi,incy)
!        call daxpy(n, ai,xr,incx,yi,incy)
!      endif
!#elseif APOLLO | HP
! --- Do real * real ---
!     if (incx .eq. 1 .and. incy .eq. 1) then
!       if (ar .ne. 0) then
!         call vec_$dmult_add(yr,xr,n, ar,yr)
!         if (sw) call vec_$dmult_add(yi,xi,n, ar,yi)
!       endif
!       if (ai .ne. 0 .and. sw) then
!         call vec_$dmult_add(yi,xr,n, ai,yi)
!         call vec_$dmult_add(yr,xi,n,-ai,yr)
!       endif
!     else
!       if (ar .ne. 0) then
!         call vec_$dmult_add_i(yr,incy,xr,incx,n, ar,yr,incy)
!         if (sw) call vec_$dmult_add_i(yi,incy,xi,incx,n, ar,yi,incy)
!       endif
!       if (ai .ne. 0 .and. sw) then
!         call vec_$dmult_add_i(yi,incy,xr,incx,n, ai,yi,incy)
!         call vec_$dmult_add_i(yr,incy,xi,incx,n,-ai,yr,incy)
!       endif
!     endif
!#elseC
      ix = 1
      iy = 1
      if (incx  <  0) ix = (1-n)*incx + 1
      if (incy  <  0) iy = (1-n)*incy + 1
      if (sw .and. ai  /=  0) then
        if (ar  /=  0) then
!     --- Case ar != 0 && ai != 0 ---
          do  10  i = 1, n
            yr(iy) = yr(iy) + ar*xr(ix) - ai*xi(ix)
            yi(iy) = yi(iy) + ar*xi(ix) + ai*xr(ix)
            ix = ix + incx
            iy = iy + incy
   10     continue
        else
!     --- Case ar == 0 && ai != 0 ---
          do  20  i = 1, n
            yr(iy) = yr(iy) - ai*xi(ix)
            yi(iy) = yi(iy) + ai*xr(ix)
            ix = ix + incx
            iy = iy + incy
   20     continue
        endif
      else
        if (ar  ==  0) return
!     --- Case ar != 0 && ai == 0 ---
        if (sw) then
          do  30  i = 1, n
            yr(iy) = yr(iy) + ar*xr(ix)
            yi(iy) = yi(iy) + ar*xi(ix)
            ix = ix + incx
            iy = iy + incy
   30     continue
        else
          do  40  i = 1, n
            yr(iy) = yr(iy) + ar*xr(ix)
            ix = ix + incx
            iy = iy + incy
   40     continue
        endif        
      endif
!#endif
      end
