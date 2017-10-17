 
      subroutine tql1(n,d,e,ierr)
          use mod_precision

!
      implicit none
      integer i,j,k,l,m,n,ii,l1,l2,mml,ierr
      real(dp) :: d(n),e(n)
      real(dp) :: c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
!
!     this subroutine is a translation of the algol procedure tql2,
!     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
!     wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
!
!     this subroutine finds the eigenvalues [and eigenvectors]
!     of a symmetric tridiagonal matrix by the ql method.
!     the eigenvectors of a full symmetric matrix can also
!     be found if  tred2  has been used to reduce this
!     full matrix to tridiagonal form.
!
!     on input
!
!        n is the order of the matrix.
!
!        d contains the diagonal elements of the input matrix.
!
!        e contains the subdiagonal elements of the input matrix
!          in its last n-1 positions.  e(1) is arbitrary.
!
!      on output
!
!        d contains the eigenvalues in ascending order.  if an
!          error exit is made, the eigenvalues are correct but
!          unordered for indices 1,2,...,ierr-1.
!
!        e has been destroyed.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the j-th eigenvalue has not been
!                     determined after 30 iterations.
!
!     calls pythag for  sqrt(a*a + b*b) .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!

!        integer(8) :: clocks_1, clocks_2
!        integer(8), save :: cmltv_clocks=0,  count_rate=0, count_max = 0
!        
!        if (cmltv_clocks == 0 ) call system_clock(count_rate = count_rate, count_max = count_max )
!        call system_clock(clocks_1)

!        print *, 'called tql1'

      ierr = 0
      if (n  ==  1) go to 1001
!
      do 100 i = 2, n
  100 e(i-1) = e(i)
!
      f = 0.0_dp
      tst1 = 0.0_dp
      e(n) = 0.0_dp
!
      do 240 l = 1, n
         j = 0
         h = abs(d(l)) + abs(e(l))
         if (tst1  <  h) tst1 = h
!     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + abs(e(m))
            if (tst2  ==  tst1) go to 120
!     .......... e(n) is always zero, so there is no exit
!                through the bottom of the loop ..........
  110    continue
!
  120    if (m  ==  l) go to 220
  130    if (j  ==  30) go to 1000
         j = j + 1
!     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0_dp * e(l))
         r = pythag(p,1.0_dp)
         d(l) = e(l) / (p + sign(r,p))
         d(l1) = e(l) * (p + sign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2  >  n) go to 145
!
         do 140 i = l2, n
  140    d(i) = d(i) - h
!
  145    f = f + h
!     .......... ql transformation ..........
         p = d(m)
         c = 1.0_dp
         c2 = c
         el1 = e(l1)
         s = 0.0_dp
         mml = m - l
!     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
!
  200    continue
!
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + abs(e(l))
         if (tst2  >  tst1) go to 130
  220    d(l) = d(l) + f
  240 continue
!     .......... order eigenvalues and eigenvectors ..........
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
!
         do 260 j = ii, n
            if (d(j)  >=  p) go to 260
            k = j
            p = d(j)
  260    continue
!
         if (k  ==  i) go to 300
         d(k) = d(i)
         d(i) = p
!
  300 continue
!
      go to 1001
!     .......... set error -- no convergence to an
!                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
!       call system_clock(clocks_2)
!        cmltv_clocks = cmltv_clocks + clocks_2 - clocks_1
!        write(6,'("tql1 wallclock time: ",g,"s cc:",g,"s")') &
!        &      real(clocks_2 - clocks_1,8)/real(count_rate,8), &
!        &      cmltv_clocks/real(count_rate,8)
       
      end
