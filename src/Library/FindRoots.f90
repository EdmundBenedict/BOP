 
      subroutine findroots(a,m,root,ifail)
          use mod_precision


!
!    This is a routine to find the roots of a real polynomial.
!     This is based on Numerical Recipe routines.
!

      implicit none

!
!    Declare the simple variables.
!

      complex(dp) :: x,b,c

      real(dp) :: eps

      integer maxm,m,j,jj,ifail

!
!    Define the parameters.
!

      parameter(eps = 1.0e-9_dp)
      parameter(maxm = 101)

!
!    Define the arrays.
!

      complex(dp) :: root(m)
      complex(dp) :: ad(maxm)

      real(dp) :: a(m+1)

!
!    Find first estimate of roots.
!

      do j = 1,m+1,1
         ad(j) = cmplx(a(j), kind=dp)
      enddo

      do j = m,1,-1
         x = cmplx(0.0_dp, kind=dp)
         call laguer(ad,j,x,eps,.false.,ifail)
         if (ifail /= 0) return
         if (abs(aimag(x)) <= 2.0_dp*(eps**2)*abs(real(x))) & 
     &       x = cmplx(real(x), kind=dp)
         root(j) = x
         b = ad(j+1)
         do jj = j,1,-1
            c = ad(jj)
            ad(jj) = b
            b = x*b+c
         enddo
      enddo

!
!    Polish roots.
!

      do j = 1,m+1,1
         ad(j) = cmplx(a(j), kind=dp)
      enddo

      do j = 1,m,1
         call laguer(ad,m,root(j),eps,.true.,ifail)
         if (ifail /= 0) return
      enddo

      end

!-----------------------------------------------------------------

      subroutine laguer(a,m,x,eps,polish,ifail)
          use mod_precision


!
!    Finds one root of a polynomial by Laguerre's method.
!

      implicit none

!
!    Declare the simple variables.
!

      complex(dp) :: x,zero,dx,x1,b,d,f,g,h,sq,gp,gm,g2

      real(dp) :: eps,epss,dxold,err,abx,cdx

      integer m,maxit,iter,j,ifail

      logical polish

!
!    Define parameters.
!

      parameter(zero = (0.0_dp,0.0_dp))
      parameter(epss = 1.0e-14_dp)
      parameter(maxit = 100)

!
!    Declare the arrays.
!

      complex(dp) :: a(m+1)

!
!    Find the root.
!

      ifail = 0

      dxold = abs(x)
      do iter = 1,maxit,1
         b = a(m+1)
         err = abs(b)
         d = zero
         f = zero
         abx = abs(x)
         do j = m,1,-1
            f = x*f+d
            d = x*d+b
            b = x*b+a(j)
            err = abs(b)+abx*err
         enddo
         err = epss*err
         if (abs(b) <= err) then
            return
         else
            g = d/b
            g2 = g*g
            h = g2 - 2.0_dp*f/b
            sq = sqrt(cmplx(m-1, kind=dp)*(cmplx(m, kind=dp)*h-g2))
            gp = g + sq
            gm = g - sq
            if (abs(gp) < abs(gm)) gp = gm
            dx = cmplx(m, kind=dp)/gp
         endif
         x1 = x - dx
         if (x == x1) return
         x = x1
         cdx = abs(dx)
         dxold = cdx
         if ((.not.polish).and.(cdx <= eps*abs(x))) return
      enddo

      write(6,'(''LAGUER: Too many iterations.'')')
      ifail = 1

      end

