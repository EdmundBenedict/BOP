 
      subroutine ygefa(ar,ai,lda,n,ipvt,info)
          use mod_precision

! Adapted from linpack, but uses no complex arithmetic
!     implicit none
      integer lda,n,ipvt(1),info
      real(dp) :: ar(lda,1),ai(lda,1)
!
!     ygefa factors a complex*16 matrix by gaussian elimination.
!
!     ygefa is an adaptation of zgefa, using real arithmetic
!     on entry
!
!        ar,ai   double precision (lda, n)
!                the matrix to be factored.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!     on return
!
!        a       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                the factorization can be written  a = l*u  where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.
!
!        ipvt    integer(n)
!                an integer vector of pivot indices.
!
!        info    integer
!                = 0  normal value.
!                = k  if  u(k,k) .eq. 0.0 .  this is not an error
!                     condition for this subroutine, but it does
!                     indicate that ygedi will divide by zero
!                     if called.
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas-adapted yaxpy,yscal,iyamax
!     fortran abs
!
!     internal variables
!
      integer iyamax,j,k,kp1,l,nm1
!#ifdefC DCMPLX
!      complex*16 t,zdum,zdumr,zdumi
!      double precision cabs1,dreal,dimag
!      dreal(zdumr) = zdumr
!      dimag(zdumi) = (0.0d0,-1.0d0)*zdumi
!      cabs1(zdum) = abs(dreal(zdum)) + abs(dimag(zdum))
!#else
      real(dp) :: tmp,t(2)
!#endif
!
!     gaussian elimination with partial pivoting
!
      info = 0
      nm1 = n - 1
      if (nm1  <  1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
!
!        find l = pivot index
!
         l = iyamax(n-k+1,ar(k,k),ai(k,k),1) + k - 1
         ipvt(k) = l
!
!        zero pivot implies this column already triangularized
!
         if (abs(ar(l,k))+abs(ai(l,k))  ==  0.0_dp) go to 40
!
!           interchange if necessary
!
            if (l  ==  k) go to 10
               tmp = ar(l,k)
               ar(l,k) = ar(k,k)
               ar(k,k) = tmp
               tmp = ai(l,k)
               ai(l,k) = ai(k,k)
               ai(k,k) = tmp
   10       continue
!
!           compute multipliers
!
            call cdiv(-1e0_dp,0e0_dp,ar(k,k),ai(k,k),t(1),t(2))
            call yscal(n-k,t(1),t(2),ar(k+1,k),ai(k+1,k),1)
!
!           row elimination with column indexing
!
            do 30 j = kp1, n
               t(1) = ar(l,j)
               t(2) = ai(l,j)
               if (l  ==  k) go to 20
                  ar(l,j) = ar(k,j)
                  ai(l,j) = ai(k,j)
                  ar(k,j) = t(1)
                  ai(k,j) = t(2)
   20          continue
               call yaxpy(n-k,t(1),t(2),ar(k+1,k),ai(k+1,k), & 
     &           1,ar(k+1,j),ai(k+1,j),1,.true.)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (abs(ar(n,n))+abs(ai(n,n))  ==  0.0_dp) info = n
      end

      subroutine ygedi(ar,ai,lda,n,ipvt,det,work,job)
          use mod_precision

!     implicit none
      integer lda,n,ipvt(1),job
!r  This version is a real adaptatin of zgedi
      real(dp) :: ar(lda,1),ai(lda,1),work(2,1),det(2,2)
!
!     ygedi computes the determinant and inverse of a matrix
!     using the factors computed by ygefa.
!
!     on entry
!
!        ar,ai   double precision (lda, n)
!                the output from ygefa.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!        ipvt    integer(n)
!                the pivot vector from ygefa.
!
!        work    complex*16(n)
!                work vector.  contents destroyed.
!
!        job     integer
!                = 11   both determinant and inverse.
!                = 01   inverse only.
!                = 10   determinant only.
!
!     on return
!
!        a       inverse of original matrix if requested.
!                otherwise unchanged.
!
!        det     complex*16(2)
!                determinant of original matrix if requested.
!                otherwise not referenced.
!                determinant = det(1) * 10.0**det(2)
!                with  1.0 .le. dcabs1<(det(1)) .lt. 10.0
!                or  det(1) .eq. 0.0 .
!
!     error condition
!
!        a division by zero will occur if the input factor contains
!        a zero on the diagonal and the inverse is requested.
!        it will not occur if the subroutines are called correctly
!        and if zgeco has set rcond .gt. 0.0 or ygefa has set
!        info .eq. 0 .
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas-adapted yaxpy,yscal,yswap
!
!     internal variables
!
      integer i,j,k,kb,kp1,l,nm1
      real(dp) :: ten
      real(dp) :: t(2),tmp(2)
!
!     compute determinant
!
      if (job/10  ==  0) go to 70
         det(1,1) = 1
         det(2,1) = 0
         det(1,2) = 0
         det(2,2) = 0
         ten = 10e0_dp
         do 50 i = 1, n
            if (ipvt(i)  /=  i) det(1,1) = -det(1,1)
            if (ipvt(i)  /=  i) det(2,1) = -det(2,1)
            call cpy(ar(i,i),ai(i,i),det,det(2,1),det,det(2,1))
!        ...exit
            if (abs(det(1,1))+abs(det(2,1))  ==  0.0_dp) go to 60
   10       if (abs(det(1,1))+abs(det(2,1))  >=  1.0_dp) go to 20
!        ...   det(1) = dcmplx(ten,0d0)*det(1)
               call cpy(ten,0e0_dp,det(1,1),det(2,1),det(1,1),det(2,1))
               det(1,2) = det(1,2) - 1
            go to 10
   20       continue
   30       if (abs(det(1,1))+abs(det(2,1))  <  ten) go to 40
!        ... complex divide ... det(1) = det(1)/dcmplx(ten,0.0d0)
               call cdiv(det(1,1),det(2,1),ten,0e0_dp,det(1,1),det(2,1))
               det(1,2) = det(1,2) + 1
            go to 30
   40       continue
   50    continue
   60    continue
   70 continue
!
!     compute inverse(u)
!
      if (mod(job,10)  ==  0) go to 150
         do 100 k = 1, n
            call cdiv(1e0_dp,0e0_dp,ar(k,k),ai(k,k),ar(k,k),ai(k,k))
            t(1) = -ar(k,k)
            t(2) = -ai(k,k)
            call yscal(k-1,t(1),t(2),ar(1,k),ai(1,k),1)
            kp1 = k + 1
            if (n  <  kp1) go to 90
            do 80 j = kp1, n
               t(1) = ar(k,j)
               t(2) = ai(k,j)
               ar(k,j) = 0
               ai(k,j) = 0
               call yaxpy(k,t(1),t(2),ar(1,k),ai(1,k),1, & 
     &           ar(1,j),ai(1,j),1,.true.)
   80       continue
   90       continue
  100    continue
!
!        form inverse(u)*inverse(l)
!
         nm1 = n - 1
         if (nm1  <  1) go to 140
         do 130 kb = 1, nm1
            k = n - kb
            kp1 = k + 1
            do 110 i = kp1, n
               work(1,i) = ar(i,k)
               work(2,i) = ai(i,k)
               ar(i,k) = 0
               ai(i,k) = 0
  110       continue
            do 120 j = kp1, n
               t(1) = work(1,j)
               t(2) = work(2,j)
               call yaxpy(n,t(1),t(2),ar(1,j),ai(1,j),1, & 
     &           ar(1,k),ai(1,k),1,.true.)
  120       continue
            l = ipvt(k)
            if (l  /=  k) then
              call dswap(n,ar(1,k),1,ar(1,l),1)
              call dswap(n,ai(1,k),1,ai(1,l),1)
            endif
  130    continue
  140    continue
  150 continue
      return
      end
