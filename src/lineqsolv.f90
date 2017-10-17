 
      subroutine ludcmp(a,n,np,indx,d)
          use mod_precision


      implicit none

!
!    Declare the simple variables.
!

      real(dp) :: tiny,d,aamax,sum,dum

      integer nmax,n,np,i,j,k,imax

!
!    Define the parameters.
!

      parameter (nmax=100,tiny=1.0e-20_dp)

!
!    Declare the arrays.
!

      real(dp) :: a(np,np)
      real(dp) :: vv(nmax)

      integer indx(n)

!
!    Carry out the decomposition.
!

      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)) > aamax) aamax=abs(a(i,j))
11      continue
        if (aamax == 0.) print *, "pause 'Singular matrix.'"
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        if (j > 1) then
          do 14 i=1,j-1
            sum=a(i,j)
            if (i > 1)then
              do 13 k=1,i-1
                sum=sum-a(i,k)*a(k,j)
13            continue
              a(i,j)=sum
            endif
14        continue
        endif
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          if (j > 1)then
            do 15 k=1,j-1
              sum=sum-a(i,k)*a(k,j)
15          continue
            a(i,j)=sum
          endif
          dum=vv(i)*abs(sum)
          if (dum >= aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j /= imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(j /= n)then
          if(a(j,j) == 0.)a(j,j)=tiny
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      if(a(n,n) == 0.)a(n,n)=tiny
      return
      end

!--------------------------------------------------------------------------

      subroutine lubksb(a,n,np,indx,b)
          use mod_precision


      implicit none

!
!    Declare the simple variables.
!

      real(dp) :: sum

      integer n,np,i,ii,ll,j

!
!    Declare the arrays.
!

      real(dp) :: a(np,np)
      real(dp) :: b(n)

      integer indx(n)

!
!    Carry out the back substitution.
!

      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii /= 0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum /= 0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        if(i < n)then
          do 13 j=i+1,n
            sum=sum-a(i,j)*b(j)
13        continue
        endif
        b(i)=sum/a(i,i)
14    continue
      return
      end
