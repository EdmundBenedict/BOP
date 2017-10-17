 
      subroutine cdiv(tr,ti,dr,di,t1,t2)
          use mod_precision

!- complex divide (t1,t2) = (tr,ti) / (dr,di) 
!r Remarks
!r   Adapted from eispack.
!r   It is permissible for (t1,t2) to occupy the same address space
!r   as either (tr,ti) or (dr,di)
      real(dp) :: tr,ti,dr,di,t1,t2
      real(dp) :: rr,d,tmp

      if (abs(di)  >  abs(dr)) then
        rr = dr / di
        d = dr * rr + di
        tmp = (tr * rr + ti) / d
        t2 = (ti * rr - tr) / d
        t1 = tmp
      else
        rr = di / dr
        d = dr + di * rr
        tmp = (tr + ti * rr) / d
        t2 = (ti - tr * rr) / d
        t1 = tmp
      endif
      end
!      subroutine cdiv(ar,ai,br,bi,cr,ci)
!      double precision ar,ai,br,bi,cr,ci
!c
!c     complex division, (cr,ci) = (ar,ai)/(br,bi)
!c
!      double precision s,ars,ais,brs,bis
!      s = abs(br) + abs(bi)
!      ars = ar/s
!      ais = ai/s
!      brs = br/s
!      bis = bi/s
!      s = brs**2 + bis**2
!      cr = (ars*brs + ais*bis)/s
!      ci = (ais*brs - ars*bis)/s
!      return
!      end
