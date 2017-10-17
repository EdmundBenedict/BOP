

   subroutine vbdist(n, idx, np, amap)
!- Distribute n variable size elements in np boxes as evenly as you can.
!-------------------------------------------------------------------------------
!i Inputs:
!i   n    : number of elements
!i   idx  : cumulative array , ie. idx(1) == size(1)
!i                                 idx(i) == idx(i-1) + size(i)
!i                                 idx(n) == the sum of the sizes of all elements
!i   np    : number of boxes (processes in the currently intended usage)
!
!o Output:
!o   amap  : cumulative array of evenly mapped elements
!o           box i is assigned the elements from (amap(i)+1) to (amap(i+1)) inclusive
!
!ir Remarks
!ir This may benefit from some improvements and sophistication.


      integer, intent(in) :: n, idx(n), np
      integer, intent(out) :: amap(0:np)

!       extra +1 for forces..

      real :: av, tg, dfp, df !, dfs(0:np), tgs(np)
      integer :: a, p

      av = real(idx(n))/real(np)

!       dfs(0) = 0
      amap(0) = 0
      a = 1
!       print *, 'p,a,idx,tg:', amap(0),0,0,0
      do p = 1, np
         tg = p*av

         dfp = huge(0.0)
         df = abs(tg - idx(a))
!          print *, 'p,tg,a,idx,dfp:',p,tg,a,idx(a),dfp
         do while (df < dfp .and. a<n)
            dfp = df
            a = a + 1
            df = abs(tg - idx(a))
!             print *, '     a,idx,df:',a,idx(a),dfp,df
         end do
         if (a < n .or. df > dfp) a = a-1
!          if (a == amap(p-1)) then
!             if (a == n)
!          end if
         amap(p) = a
!          dfs(p) = abs(tg - idx(a))
!          print *, 'p,a,idx,tg:', p,a,idx(a),tg
      end do

!       tgs = (/(p*av, p=1,np)/)

!       print *, 'idx:             ', idx
!       print *, 'amap:', amap
!       print '(" dfs:    ",5(8x,f4.2))', dfs
!       print '(" tgs:             ",4(6x,f6.2))', tgs
!       print *, 'rdx:             ', idx(amap(1:))
!       print *, 'asum1:', sum(abs(idx(amap(1:)) - tgs)), sum(dfs)
! !       print *, 'asum2:', sum(abs(idx([1,2]) - tgs))

   end subroutine vbdist



