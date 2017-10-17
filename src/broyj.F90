! C#define AWRITE
      integer function broyj(n,xin,gin,ir,isw,beta,dxmx,xtol,gtol, &
     &  wc,wk,ndw,xnew)
     use mod_all_scalar, only: quiet
!      use mod_precision
! C- One Broyden step in finding gin = f[xin]-xin = 0
! C ----------------------------------------------------------------------
! Ci Inputs
! Ci   n:     number of variables
! Ci   ir:    Number of iterations of x and g.
! Ci          1 initiates a new sequence of mixing;
! Ci          broyj uses linear mixing for this iteration.
! Ci   isw    1s digit (not implemented)
! Ci          0  find minimum
! Ci          1  find maximum
! Ci         10s digit not used
! Ci        100s digit not used
! Ci       1000s digit governs convergence criterion:
! Ci          1 return when |grad| < gtol
! Ci          2 return when max dx < xtol
! Ci          3 return when either (1) or (2) satisfied
! Ci          4 return when both (1) and (2) are satisfied
! Ci   beta:  linear mixing parameter (ir=1 only)
! Ci   xin:   input vector, this iteration
! Ci   gin:   output-input vector, f[xin]-xin, this iteration
! Ci   wc:    weighting for this iteration
! Cio  wk     workspace of 2*ndw*(ir+2), ndw>=n
! Cio         wk must be preserved between calls to broyj.
! Cio         (*,1,0) x of the prior iteration.
! Cio         (*,2,0) g of the prior iteration.
! Cio         (*,1..2,1..ir-1) u and vt of this and prior iterations
! Cio         (*,1,ir) g(this iter) - g (prior iter).
! Co Outputs
! Co   xnew   estimate of x
! Co   broyj
! Cr Remarks
! Cr   2016.04 (DMT) fix aliasing issue leading to subtle failure on wide vector registers
! Cr   Adapted from Duane Johnson
! C ----------------------------------------------------------------------
      implicit none
      integer isw,ir,n,ndw
      double precision beta,dxmx,wc,xin(n),gin(n),xnew(n),xtol,gtol,wk(ndw,2,0:ir)
! C Local variables
      integer i,ip,j,k,irm1,irm2,lm,ln,nn,dinv,i1mach,isw1,isw2,isw3
! C     integer ierr
      double precision omb,dfnorm,fac1,fac2,gmi,one,zero,ddot,w0
      parameter (zero=0d0,one=1d0,nn=20)
      double precision a(nn,nn),cm(nn),w(nn),d(nn,nn)
      double precision betx,diff,gmax,xmax
! C     double precision wl(nn,3),u(nn,nn),v(nn,nn)
      save w,cm,a,w0
!       ,d
! NOTE: tbe seems to retain it's old values of d, although not for any reason I could see

!       d= 0.0_dp
      isw1 = mod(isw/10,10)
      isw2 = mod(isw/100,10)
      isw3 = mod(isw/1000,10)
      if (ir .gt. nn) then
           write(6,'(/,"broyj: increase nn, need",i3)') ir

          call panic()
!       rxi('broyj: increase nn, need',ir)
      endif
     
     
!      print *, 'XINPUT' 
!      print *, 'n, ir, isw,ndw',n, ir, isw,ndw
!      print *, xin 
!      print *, ' '
!      print *, gin
!      
!      
!       print *, 'first d',d
!       
!       print *, 'first a', a

! C --- First iteration: simple mixing ---
      if (ir .eq. 1) then
!       print *, ' ir .eq. 1'
        w0 = wc
        betx = beta
        gmax = 0
        do k = 1, n
            gmax = max(gmax,abs(gin(k)))
        end do

! C#ifdef AWRITE
!         if (ipr .ge. 30) then
!           j = isw3
!           call awrit8(' broyj:  start'//
!      .      '%?#(n==2|n==3|n==4)#  xtol=%1;2g#%j#'//
!      .      '%?#(n==1|n==3|n==4)#  gtol=%1;2g#%j#  beta=%1;2g'//
!      .      '  w0=%1;2g  isw=%i  gmax=%1;2g',' ',80,i1mach(2),
!      .      j,xtol,j,gtol,beta,w0,isw,gmax)
!         endif
! C#endif

        if (dxmx .gt. 0d0 .and. gmax .gt. dxmx) then
          betx = beta*dxmx/gmax
! C#ifdef AWRITE
!           call awrit3(' broyj:  max shift = %1;3g'//
!      .      ' is larger than dxmx = %1;3g.  Scale by %1;3g',
!      .      ' ',80,i1mach(2),gmax,dxmx,dxmx/gmax)
! C#endif
        endif
        do k = 1, n
            xnew(k) = xin(k) + betx*gin(k)
        end do

! C --- Subsequent iterations: Broyden mixing ---
      else
!       print *,'ir greater'

! C   ... Make xold, gold
        do k = 1, n
          wk(k,1,0) = xin(k) - wk(k,1,0)
          wk(k,1,ir) = gin(k) - wk(k,2,0)
        end do

! C   --- Coefficient matrices and the sum for corrections ---
! C   ... dfnorm = |g(i)-g(i-1)|, used for normalization
        dfnorm = dsqrt(ddot(n,wk(1,1,ir),1,wk(1,1,ir),1))
        fac2 = one/dfnorm
        fac1 = beta*fac2
! C   ... Shuffle each prior u,vt to prior+1 iteration
        irm1 = ir-1
        irm2 = ir-2
        do j = irm2, 1, -1
          call dcopy(n,wk(1,1,j),1,wk(1,1,j+1),1)
          call dcopy(n,wk(1,2,j),1,wk(1,2,j+1),1)
        end do
! C   ... Make u,vt for this iteration
        do k = 1, n
            wk(k,1,1) = fac1*wk(k,1,ir) + fac2*wk(k,1,0)
        end do

        do k = 1, n
            wk(k,2,1) = fac2*wk(k,1,ir)
        end do

! C   --- Make  a and b = ( w0**2 I + a )^-1 (symmetric) ---
        do  j = 1, irm2
          a(irm1,j) = ddot(n,wk(1,2,ir-j),1,wk(1,2,1),1)
          a(j,irm1) = a(irm1,j)
          cm(j) = ddot(n,wk(1,2,ir-j),1,gin(1),1)
        end do
        a(irm1,irm1) = ddot(n,wk(1,2,1),1,wk(1,2,1),1)
        cm(irm1) = ddot(n,wk(1,2,1),1,gin(1),1)
        w(irm1) = wc

! C   ... Set up and calculate beta matrix
!         print *, 'w0',w0
        do lm = 1, irm1
            do ln = lm+1, irm1
                d(ln,lm) = a(ln,lm)*w(ln)*w(lm)
                d(lm,ln) = d(ln,lm)
!                 print *, 'ln,lm'
!                 print *, 'anm,wn,wm',a(ln,lm),w(ln),w(lm)
            end do
            d(lm,lm) = w0*w0 + a(lm,lm)*w(lm)*w(lm)
!             print *, 'amm',a(lm,lm)
        end do

! C   --- Invert to make d ---
!         print *, 'd',d
        
        
        if (dinv(' ',irm1,nn,d) .ne. 0) then
!           call rx('broyj: matrix singular')
            write(6,'(/"broyj: matrix singular")')
            call panic()
        endif
!         print *, 'd',d
! C   ... Invert with singular value decomposition
! C        call svd(nn,irm1,irm1,d,wl(1,2),.true.,u,.true.,v,ierr,wl)
! C        call dpzero(d,nn**2)
! C        do  60  ln = 1, irm1
! C   60   d(ln,ln) = 1
! C        call svbksb(nn,irm1,irm1,irm1,wl(1,2),u,v,d,d,wl)
! C    ... This one sometimes hangs up
! C        call rs(nn,irm1,d,wl(1,3),1,v,wl,wl(1,2),ierr)
! C        call dpzero(d,nn**2)
! C        do  60  ln = 1, irm1
! C          print *, 'evl',ln, wl(ln,3)
! C   60   d(ln,ln) = 1
! C        call svbksb(nn,irm1,irm1,irm1,wl(1,3),v,v,d,d,wl)
! 
! C   --- xnew <- vector for the new iteration ---
        do k = 1, n
            xnew(k) = xin(k) + beta*gin(k)
        end do
        do i = 1, irm1
            gmi = zero
            do ip = 1, irm1
                gmi = gmi + cm(ip)*d(ip,i)*w(ip)
            end do
            call daxpy(n,-gmi*w(i),wk(1,1,ir-i),1,xnew(1),1)
        end do

! C   ... Cap to maximum allowed shift xnew-xin
        if (dxmx .gt. 0d0) then
          diff = 0
          do k = 1, n
            diff = max(diff,abs(xnew(k)-xin(k)))
          end do
          if (diff .gt. dxmx) then
            betx = dxmx/diff
! C#ifdef AWRITE
          if (.not. quiet) then
            write(6,'(" broyj:  max shift = ",f10.8)')diff
            write(6,'(" is larger than dxmx = ",f10.8)')dxmx
            write(6,'(" Scale by ",f10.5)') dxmx/diff
          endif
!             call awrit3(' broyj:  max shift = %1;3g'//
!      .        ' is larger than dxmx = %1;3g.  Scale by %1;3g',
!      .        ' ',80,i1mach(2),diff,dxmx,dxmx/diff)
! C#endif
            omb = 1d0-betx
            do k = 1, n
                xnew(k) = betx*xnew(k) + omb*xin(k)
            end do
          endif
        endif

      endif

! C --- Cleanup, setup for next call ---
      xmax = 0
      gmax = 0
      diff = 0
      do  k = 1, n
        xmax = max(xmax,abs(xnew(k)-xin(k)))
        gmax = max(gmax,dabs(gin(k)))
        diff = diff + (xnew(k)-xin(k))**2
        wk(k,2,0) = gin(k)
        wk(k,1,0) = xin(k)
      end do
      diff = dsqrt(diff/n)

      j = ir+1
      if (isw3 .ne. 0 .and. (gmax .eq. 0 .or. &
     &  gmax .lt. gtol .and. xmax.lt.xtol .and. isw3.eq.4 .or. &
     &  gmax .lt. gtol .and. (isw3.eq.1 .or. isw3.eq.3)  .or. &
     &  xmax .lt. xtol .and. (isw3.eq.2 .or. isw3.eq.3)  .or. &
     &  gmax .lt. gtol .and. (isw3.eq.1 .or. isw3.eq.3))) j = 0
! C#ifdef AWRITE
      if (j .eq. 0) then
!         call awrit3(' broyj: converged to max dx'//
!      .    '=%1;2g, gmax=%1;2g using %i iterations',' ',80,
!      .    i1mach(2),xmax,gmax,ir)
          if (.not. quiet) then
              write(6,'("broyj: converged to max dx ",f10.8)') xmax
              write(6,'("gmax=",f10.8," using ",i3," iterations")') gmax,ir
          endif
      elseif (ir .ne. 1) then
!         call awrit4(' broyj:  ir=%i  dxmax=%1;2g  gmax=%1;2g'//
!      .    '  wc=%1;2g',' ',80,i1mach(2),ir,xmax,gmax,wc)
          if (.not. quiet) write(6,'(" broyj:  ir=",i3,"  dxmax=",f10.8," gmax=",f10.8," wc=",f5.2)') ir,xmax, gmax,wc
     
      endif
! C#endif
      broyj = j
      end
      
   subroutine dgefa(a,lda,n,ipvt,info)
      integer lda,n,info
      double precision a(lda,n),ipvt(*)
! c     dgefa factors a double precision matrix by gaussian elimination.
! c
! c     dgefa is usually called by dgeco, but it can be called
! c     directly with a saving in time if  rcond  is not needed.
! c     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
! c
! c     on entry
! c
! c        a       double precision(lda, n)
! c                the matrix to be factored.
! c
! c        lda     integer
! c                the leading dimension of the array  a .
! c
! c        n       integer
! c                the order of the matrix  a .
! c
! c     on return
! c
! c        a       an upper triangular matrix and the multipliers
! c                which were used to obtain it.
! c                the factorization can be written  a = l*u  where
! c                l  is a product of permutation and unit lower
! c                triangular matrices and  u  is upper triangular.
! c
! c        ipvt    integer(n)
! c                an integer vector of pivot indices.
! c
! c        info    integer
! c                = 0  normal value.
! c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
! c                     condition for this subroutine, but it does
! c                     indicate that dgesl or dgedi will divide by zero
! c                     if called.  use  rcond  in dgeco for a reliable
! c                     indication of singularity.
! c
! c     linpack. this version dated 08/14/78 .
! c     cleve moler, university of new mexico, argonne national lab.
! c
! c     subroutines and functions
! c
! c     blas daxpy,dscal,idamax
! c
! c     internal variables
! c
      double precision t
      integer idamax,j,k,kp1,l,nm1
! c
! c
! c     gaussian elimination with partial pivoting
! c
! C#ifdefC CRAY
! C      call sgefa(a,lda,n,ipvt,info)
! C#else
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) goto 70
      do  60  k = 1, nm1
         kp1 = k + 1
! c
! c        find l = pivot index
! c
         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
!          print *, '1 ipvt',ipvt(:n)
         
! c
! c        zero pivot implies this column already triangularized
! c
         if (a(l,k) .eq. 0.0d0) goto 40
! c
! c           interchange if necessary
! c
            if (l .eq. k) goto 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
! c
! c           compute multipliers
! c
            t = -1.0d0/a(k,k)
            call dscal(n-k,t,a(k+1,k),1)
! c
! c           row elimination with column indexing
! c
            do  30  j = kp1, n
               t = a(l,j)
               if (l .eq. k) goto 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         goto 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      
!       print *, '2 ipvt',ipvt(:n)
      if (a(n,n) .eq. 0.0d0) info = n
      return
! C#endif
      end
      
      
    double precision function dsum(n,dx,incx)
! c
! c     takes the sum of the values.  Adapted from:
! c     jack dongarra, linpack, 3/11/78.
! c
      double precision dx(*)
      integer i,incx,n,nincx

      dsum = 0d0
      if (n .le. 0) return

      nincx = n*incx
      do  10  i = 1, nincx,incx
        dsum = dsum + dx(i)
   10 continue
      end
      
      
      integer function dinv(cs,n,lda,a)
! C- Inversion of a double precision matrix
      implicit none
      character*1 cs
      integer n,lda
      double precision a(lda,lda)
      integer ldw,i
      real(8), allocatable :: w(:)

      ldw = n
      allocate(w(ldw*(n+1)))
      call dqinv(cs,a,lda,2,n,w,ldw,i)
      deallocate(w)
      dinv = i
      end

      
      
