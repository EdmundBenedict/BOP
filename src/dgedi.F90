      subroutine dgedi(a,lda,n,ipvt,det,work,job)
      integer lda,n,ipvt(*),job
      double precision a(lda,n),det(2),work(n)
! c
! c     dgedi computes the determinant and inverse of a matrix
! c     using the factors computed by dgeco or dgefa.
! c
! c     on entry
! c
! c        a       double precision(lda, n)
! c                the output from dgeco or dgefa.
! c
! c        lda     integer
! c                the leading dimension of the array  a .
! c
! c        n       integer
! c                the order of the matrix  a .
! c
! c        ipvt    integer(n)
! c                the pivot vector from dgeco or dgefa.
! c
! c        work    double precision(n)
! c                work vector.  contents destroyed.
! c
! c        job     integer
! c                = 11   both determinant and inverse.
! c                = 01   inverse only.
! c                = 10   determinant only.
! c
! c     on return
! c
! c        a       inverse of original matrix if requested.
! c                otherwise unchanged.
! c
! c        det     double precision(2)
! c                determinant of original matrix if requested.
! c                otherwise not referenced.
! c                determinant = det(1) * 10.0**det(2)
! c                with  1.0 .le. dabs(det(1)) .lt. 10.0
! c                or  det(1) .eq. 0.0 .
! c
! c     error condition
! c
! c        a division by zero will occur if the input factor contains
! c        a zero on the diagonal and the inverse is requested.
! c        it will not occur if the subroutines are called correctly
! c        and if dgeco has set rcond .gt. 0.0 or dgefa has set
! c        info .eq. 0 .
! c
! c     linpack. this version dated 08/14/78 .
! c     cleve moler, university of new mexico, argonne national lab.
! c
! c     subroutines and functions
! c
! c     blas daxpy,dscal,dswap
! c     fortran dabs,mod
! c
! c     internal variables
! c
      double precision t
      double precision ten
      integer i,j,k,kb,kp1,l,nm1
! c
! c
! c     compute determinant
! c
! C#ifdefC CRAY
! C      call sgedi(a,lda,n,ipvt,det,work,job)
! C#else
!       print *, 'dgedi ipv',ipvt(:n)

      if (job/10 .eq. 0) goto 70
         det(1) = 1.0d0
         det(2) = 0.0d0
         ten = 10.0d0
         do  50  i = 1, n
            if (ipvt(i) .ne. i) det(1) = -det(1)
            det(1) = a(i,i)*det(1)
! c        ...exit
            if (det(1) .eq. 0.0d0) goto 60
   10       if (dabs(det(1)) .ge. 1.0d0) goto 20
               det(1) = ten*det(1)
               det(2) = det(2) - 1.0d0
            goto 10
   20       continue
   30       if (dabs(det(1)) .lt. ten) goto 40
               det(1) = det(1)/ten
               det(2) = det(2) + 1.0d0
            goto 30
   40       continue
   50    continue
   60    continue
   70 continue
! c
! c     compute inverse(u)
! c
!        do i = 1,lda
!           do k = 1,n
!               print *, 'a',i,k,a(i,k)
!           enddo
!        enddo
       
      if (mod(job,10) .eq. 0) goto 150
         do  100  k = 1, n
!             print *, 'akk',a(k,k)
            a(k,k) = 1.0d0/a(k,k)
!             print *, '1 k,1/ak,a1k',k,a(k,k),a(1,k)
            t = -a(k,k)
            call dscal(k-1,t,a(1,k),1)
            
!             print *, '2 k-1,t,a(1,k),a(k,k)',k-1,t,a(1,k)
            kp1 = k + 1
            if (n .lt. kp1) goto 90
            do  80  j = kp1, n
               t = a(k,j)
               a(k,j) = 0.0d0
               call daxpy(k,t,a(1,k),1,a(1,j),1)
!                 print *, '3 t,a(1,k),ak,k',t,a(1,k),a(k,k)
   80       continue
   90       continue
  100    continue
  
!          print *, 'a',a(:,:)
! c
! c        form inverse(u)*inverse(l)
! c
         nm1 = n - 1
         if (nm1 .lt. 1) goto 140
         do  130  kb = 1, nm1
            k = n - kb
            kp1 = k + 1
            do  110  i = kp1, n
               work(i) = a(i,k)
               a(i,k) = 0.0d0
  110       continue
            do  120  j = kp1, n
               t = work(j)
               call daxpy(n,t,a(1,j),1,a(1,k),1)
  120       continue
            l = ipvt(k)
!             print *, 'dgedi'
!             print *, ipvt(:n)
!             print *, l,k
            if (l .ne. k) call dswap(n,a(1,k),1,a(1,l),1)
            
!             print *, 'a',a(:,:)
  130    continue
  140    continue
  150 continue
      return
! C#endif
      end
