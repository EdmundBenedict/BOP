 
      subroutine findroots_qr(a_arg,m,root_arg,ifail)
          use mod_precision


!
!    This is a routine to find the roots of a real polynomial.
!     This is based on Numerical Recipe routines.
!

!
!  Attempted conversion to quadruple precision.
!

      implicit none

!
!    Declare the simple variables.
!

      real*16 x(2),b(2),c(2), b2(2)
      real*16 eps

      real(dp) :: realp, imp
      integer maxm,m,j,jj,ifail

!
!    Define the parameters.
!

      parameter(eps  = 1.0q-9)
      parameter(maxm = 101)

!
!    Define the arrays.
!

      complex(dp) :: root_arg(m)
      real(dp) :: a_arg(m+1)

      real*16 root(m,2)
      real*16 a(m+1)
      real*16 copy_root(2)

      real*16 ad(maxm,2)

!
!  Copy polynomial coefficients into real*16 arrays.
!   We can only get precision in A to 15 sig. digits from the REAL*8 argument.
!
      
      do j = 1, m+1, 1
         a(j) = qextd(a_arg(j))
      enddo

!
!    Find first estimate of roots.
!

      do j = 1,m+1,1
         ad(j,1) = a(j)
         ad(j,2) = 0.0q0
      enddo

      do j = m,1,-1

         x(1) = 0.0q0
         x(2) = 0.0q0

         call laguer_qr(ad,j,x,eps,.false.,ifail)

         if (ifail /= 0) return

         if (qabs(x(2)) <= 2.0q0*(eps**2)*qabs(x(1))) then
            x(1) = x(1)
            x(2) = 0.0_dp
         endif

         root(j,1) = x(1)
         root(j,2) = x(2)

         b(1) = ad(j+1,1)
         b(2) = ad(j+1,2)

         do jj = j,1,-1

            c(1) = ad(jj,1)
            c(2) = ad(jj,2)
            ad(jj,1) = b(1)
            ad(jj,2) = b(2)

            b2(1) = x(1)*b(1) - x(2)*b(2) + c(1)
            b2(2) = x(2)*b(1) + x(1)*b(2) + c(2)

            b(1) = b2(1)
            b(2) = b2(2)

         enddo


      enddo


!
!    Polish roots.
!

      do j = 1,m+1,1

         ad(j,1) = a(j)
         ad(j,2) = 0.0q0

      enddo

      do j = 1,m,1
         copy_root(1) = root(j,1)
         copy_root(2) = root(j,2)
         call laguer_qr(ad,m,copy_root,eps,.true.,ifail)
         if (ifail /= 0) return
         root(j,1) = copy_root(1)
         root(j,2) = copy_root(2)
      enddo

!
! Copy roots back into REAL*8 FORMAT.
!
      do j = 1, m

         realp = dbleq(root(j,1))         
         imp   = dbleq(root(j,2))
         root_arg(j) = cmplx(realp,imp, kind=dp)

      enddo


      end

!-----------------------------------------------------------------

      subroutine laguer_qr(a,m,x,eps,polish,ifail)
          use mod_precision


!
!    Finds one root of a polynomial by Laguerre's method.
!

!
!    Quadruple precision version. Can hopefully push EPSS
!    down below 1.0Q-15 which is machine precision with
!    double precision to polish off roots more accurately.
!

      implicit none

!
!    Declare the simple variables.
!

      real*16 x(2),zero(2),dx(2),x1(2),b(2), fdb(2), arg(2)
      real*16 d(2),f(2),g(2),h(2),sq(2),gp(2),gm(2),g2(2), marg(2)

      real*16 f2(2),d2(2),b2(2)

      real*16 eps,epss,dxold,err,abx,cdx

      integer m,maxit,iter,j,ifail,maxm

      logical polish

!
!    Define parameters.
!

      parameter(epss = 1.0q-20)  

      parameter(maxit = 100)    
      parameter(maxm = 101)
!
!    Declare the arrays.
!

      real*16 a(maxm,2)

      zero(1) = 0.0q0
      zero(2) = 0.0q0

!
!    Find the root.
!

      ifail = 0

      dxold = qsqrt(x(1)**2 + x(2)**2)

      do iter = 1,maxit,1

         b(1) = a(m+1,1)
         b(2) = a(m+1,2)

         err = qsqrt(b(1)**2 + b(2)**2)

         d(1) = zero(1)
         d(2) = zero(2)

         f(1) = zero(1)
         f(2) = zero(2)

         abx = qsqrt(x(1)**2 + x(2)**2)

         do j = m,1,-1

            f2(1) = x(1)*f(1) - x(2)*f(2) + d(1)
            f2(2) = x(2)*f(1) + x(1)*f(2) + d(2)
            
            f(1) = f2(1)
            f(2) = f2(2)

            d2(1) = x(1)*d(1) - x(2)*d(2) + b(1)
            d2(2) = x(2)*d(1) + x(1)*d(2) + b(2)
            
            d(1) = d2(1)
            d(2) = d2(2)

            b2(1) = x(1)*b(1) - x(2)*b(2) + a(j,1)
            b2(2) = x(2)*b(1) + x(1)*b(2) + a(j,2)

            b(1) = b2(1)
            b(2) = b2(2)

            err = qsqrt(b(1)**2 + b(2)**2) + abx*err

         enddo

         err = epss*err

          if (qsqrt(b(1)**2 + b(2)**2) <= err) then

             return

          else

             call comdiv(d,b,g)

             g2(1) = g(1)**2 - g(2)**2
             g2(2) = 2.0q0*g(1)*g(2)

             call comdiv(f,b,fdb)
             fdb(1) = 2.0q0*fdb(1)
             fdb(2) = 2.0q0*fdb(2)

             h(1) = g2(1) - fdb(1)
             h(2) = g2(2) - fdb(2)

             arg(1) = qflotj(m-1) * ((qflotj(m)*h(1))-g2(1))
             arg(2) = qflotj(m-1) * ((qflotj(m)*h(2))-g2(2))
             call comsqrt(arg,sq)

             gp(1) = g(1) + sq(1)
             gp(2) = g(2) + sq(2)
             gm(1) = g(1) - sq(1)
             gm(2) = g(2) - sq(2)
             
             if ((qsqrt(gp(1)**2 + gp(2)**2)) <  & 
     &          (qsqrt(gm(1)**2 + gm(2)**2))) then
                gp(1) = gm(1)
                gp(2) = gm(2)
             endif

             marg(1) = qflotj(m)
             marg(2) = 0.0q0
             call comdiv(marg,gp,dx)

          endif

          x1(1) = x(1) - dx(1)
          x1(2) = x(2) - dx(2)

          if ((x(1) == x1(1)).and.(x(2) == x1(2))) return

          x(1) = x1(1)
          x(2) = x1(2)

          cdx = qsqrt(dx(1)**2 + dx(2)**2)
          dxold = cdx

          if ((.not.polish).and. & 
     &         (cdx <= eps*qsqrt(x(1)**2 + x(2)**2))) then
             return
          endif

      enddo

      write(6,'(''LAGUER: Too many iterations.'')')
      ifail = 1

      end


      subroutine comdiv(x,y,z)
          use mod_precision


      implicit none
!
!   Computes Z = X/Y where X and Y are complex numbers represented
!    in 2 element arrays.

      real*16 x(2),y(2),z(2)
      real*16 theta

      if (qabs(y(1)) > qabs(y(2))) then
         theta = y(2)/y(1)
         z(1) = (x(1) + theta*x(2)) / (theta*y(2) + y(1))
         z(2) = (x(2) - theta*x(1)) / (theta*y(2) + y(1))
      else
         theta = y(1)/y(2)
         z(1) = (theta*x(1) + x(2)) / (theta*y(1) + y(2))
         z(2) = (theta*x(2) - x(1)) / (theta*y(1) + y(2))
      endif

      return
      end


      subroutine comsqrt(x,y)
          use mod_precision


!  Computes complex Y = SQRT(X) where X and Y are represented
!  as two element REAL*16 arrays.

      implicit none

      real*16 x(2),y(2), sig

      if ((x(1) == 0.0q0).and.(x(2) == 0.0q0)) then
         y(1) = 0.0q0
         y(2) = 0.0q0
         return
      endif

      if (x(1) >= 0.0q0) then
         y(1) = qsqrt( (x(1) + qsqrt(x(1)**2 + x(2)**2))/2.0q0 )
         y(2) = x(2) / (2.0q0*y(1))
      else
         if (x(2) >= 0.0q0) then
            sig = 1.0q0
         else
            sig = -1.0q0
         endif
         y(2) = sig * qsqrt( ( qabs(x(1)) + & 
     &        qsqrt( x(1)**2 + x(2)**2) ) / 2.0q0)
         y(1) = x(2) / (2.0q0*y(2))
      endif

      return
      end
      
