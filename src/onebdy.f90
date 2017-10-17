 
      subroutine onebdy(u,f)
          use mod_precision


          use mod_all_scalar

          use mod_const

!
!    This is a routine to evaluate a one body time dependent potential energy,
!     and the resulting force contribution on the ions.
!
!    V1FLAG = 1  ==>  Single C(100) plane
!    V1FLAG = 2  ==>  Two sliding C(100) planes
!    V1FLAG = 3  ==>  Spherical bomb
!

      implicit none

!
!    Include constants.
!

!      include "Include/ALL.const"

!
!    Include scalars.
!

!      include "Include/ALL.scalar"

!
!    Include arrays.
!

      include "Include/Atom.array"
      include "Include/PosVel.array"

!
!    Declare the simple variables.
!

      real(dp) :: u,csrf,csrf2,lj,lj2
      real(dp) :: ftlpln,alpha
      real(dp) :: dv,bombv,bombpt
      real(dp) :: dxpln,dypln

      integer i

!
!    Declare the arrays.
!

      real(dp) :: f(3,nd)
      real(dp) :: du(3)
      real(dp) :: duup(3)
      real(dp) :: dudn(3)

!
!    Evaluate the energy and forces.
!

      if (v1flag == 1)then
         u = 0.0_dp
         do i = 1,nd,1
!           U = U + CSRF(AD(1,I),Z(I))
            u = u + lj(ad(1,i),z(i))
!           CALL DCSRF(AD(1,I),Z(I),DU)
            call dlj(ad(1,i),z(i),du)
            f(1,i) = f(1,i) - du(1)
            f(2,i) = f(2,i) - du(2)
            f(3,i) = f(3,i) - du(3)
         enddo
      elseif (v1flag == 2) then
         u = 0.0_dp
         fplnx = 0.0_dp
         fplny = 0.0_dp
         fplnz = 0.0_dp
         dxpln = vxpln*dt*real(iter, dp)
         dypln = vypln*dt*real(iter, dp)
         do i = 1,nd,1
!           U = U + CSRF2(AD(1,I),Z(I),DXPLN,DYPLN,DZPLN)
            u = u + lj2(ad(1,i),z(i),dxpln,dypln,dzpln)
!           CALL DCSRF2(AD(1,I),Z(I),DXPLN,DYPLN,DZPLN,
!    +                  DUDN,DUUP)
            call dlj2(ad(1,i),z(i),dxpln,dypln,dzpln, & 
     &                dudn,duup)
            f(1,i) = f(1,i) - dudn(1) - duup(1)
            f(2,i) = f(2,i) - dudn(2) - duup(2)
            f(3,i) = f(3,i) - dudn(3) - duup(3)
            fplnx = fplnx + duup(1)
            fplny = fplny + duup(2)
            fplnz = fplnz + duup(3)
         enddo
         ftlpln = fplnz+fxtpln
         alpha = min(1.0_dp,0.1_dp/abs(ftlpln))
         dzpln = dzpln + alpha*ftlpln
      elseif (v1flag == 3) then
         u = 0.0_dp
         bombp = 0.0_dp
         do i = 1,nd,1
            u = u + bombv(ad(1,i),bombr0)
            call dbombv(ad(1,i),bombr0,du,dv)
            f(1,i) = f(1,i) - du(1)
            f(2,i) = f(2,i) - du(2)
            f(3,i) = f(3,i) - du(3)
            bombp = bombp + dv
         enddo
         bombp = bombp/(4.0_dp*pi*bombr0*bombr0)
         bombpt = bombp + bombpx
         alpha = min(1.0_dp,0.1_dp/abs(bombpt))
         bombr0 = bombr0 + alpha*bombpt
      endif

      end

!---------------------------------------------------------------------------

       function csrf(r,za)
          use mod_precision


!
!    This is a function to evaluate the potential due to a carbon (100) surface.
!

      implicit none
      real(dp) :: csrf

!
!    Declare the simple variables.
!


      real(dp) :: x,y,z
      real(dp) :: twopi,fg1,fg2

      integer za

!
!    Define parameters.
!

      parameter (twopi = 2.0_dp*3.141592653589_dp)

!
!    Declare the arrays.
!

      real(dp) :: r(3)
      real(dp) :: p00(4)
      real(dp) :: p01(2)
      real(dp) :: p10(2)
      real(dp) :: p11(2)

!
!    Evaluate the repulsive potential.
!

      x = r(1)/2.52_dp
      y = r(2)/2.52_dp
      z = r(3)

      p00(1) =   .5845088_dp
      p00(2) =   .5750731_dp
      p00(3) =   3.239172_dp
      p00(4) =   1.625920_dp

      p01(1) =   .5580577_dp
      p01(2) =   3.668953_dp

      p10(1) =  -.5707599_dp
      p10(2) =   3.668051_dp

      p11(1) =  -.1250139_dp
      p11(2) =   3.940773_dp

      csrf = fg2(z,p00) & 
     &     + 2.0_dp*(fg1(z,p01)*cos(twopi*y)  & 
     &            + fg1(z,p10)*cos(twopi*x) & 
     &            + fg1(z,p11)*cos(twopi*(x+y)) & 
     &            + fg1(z,p11)*cos(twopi*(x-y)))

      if (za == 1) then
         csrf = csrf*0.54_dp*0.15_dp
      elseif (za == 6) then
         csrf = csrf*0.54_dp*0.40_dp
      endif

      end

!---------------------------------------------------------------------------

       function fg1(x,p)
          use mod_precision


!
!    This is a function to evaluate the variation of the fourier coefficients with Z.
!

      implicit none
      real(dp) :: fg1

!
!    Declare the simple variables.
!

      real(dp) :: x,y
      real(dp) :: a,n

      integer npar

!
!    Define the parameters.
!

      parameter(npar = 2)

!
!    Declare the arrays.
!

      real(dp) :: p(npar)

!
!    Evaluate the function
!

      a = p(1)
      n = p(2)**2
      y = 2.52_dp/x
      fg1 = a*(y**n)

      end

!---------------------------------------------------------------------------

       function fg2(x,p)
          use mod_precision


!
!    This is a function to evaluate the variation of the fourier coefficients with Z.
!

      implicit none
      real(dp) :: fg2

!
!    Declare the simple variables.
!

      real(dp) :: x,y
      real(dp) :: a,k,n1,n2

      integer npar

!
!    Define the parameters.
!

      parameter(npar = 4)

!
!    Declare the arrays.
!

      real(dp) :: p(npar)

!
!    Evaluate the function
!

      a = p(1)
      k = p(2)**2
      y = 1.0_dp/(k*x)
      n1 = p(3)**2
      n2 = p(4)**2
      fg2 = a*(y**n1 - y**n2)

      end

!---------------------------------------------------------------------------

      subroutine dcsrf(r,za,f)
          use mod_precision


!
!    This is a function to evaluate the derivative of the potential due to a
!     carbon (100) surface.
!

      implicit none

!
!    Declare the simple variables.
!

      real(dp) :: x,y,z
      real(dp) :: twopi,fg1,dfg1,dfg2

      integer za

!
!    Define parameters.
!

      parameter (twopi = 2.0_dp*3.141592653589_dp)

!
!    Declare the arrays.
!

      real(dp) :: r(3)
      real(dp) :: f(3)
      real(dp) :: p00(4)
      real(dp) :: p01(2)
      real(dp) :: p10(2)
      real(dp) :: p11(2)

!
!    Evaluate the repulsive potential.
!

      x = r(1)/2.52_dp
      y = r(2)/2.52_dp
      z = r(3)

      p00(1) =   .5845088_dp
      p00(2) =   .5750731_dp
      p00(3) =   3.239172_dp
      p00(4) =   1.625920_dp

      p01(1) =   .5580577_dp
      p01(2) =   3.668953_dp

      p10(1) =  -.5707599_dp
      p10(2) =   3.668051_dp

      p11(1) =  -.1250139_dp
      p11(2) =   3.940773_dp

      f(1) = -(2.0_dp*twopi/2.52_dp)* & 
     &             (fg1(z,p10)*sin(twopi*x) & 
     &            + fg1(z,p11)*sin(twopi*(x+y)) & 
     &            + fg1(z,p11)*sin(twopi*(x-y)))

      f(2) = -(2.0_dp*twopi/2.52_dp)* & 
     &             (fg1(z,p01)*sin(twopi*(y))  & 
     &            + fg1(z,p11)*sin(twopi*(x+y)) & 
     &            - fg1(z,p11)*sin(twopi*(x-y)))

      f(3) = dfg2(z,p00) & 
     &     + 2.0_dp*(dfg1(z,p01)*cos(twopi*y)  & 
     &            + dfg1(z,p10)*cos(twopi*x) & 
     &            + dfg1(z,p11)*cos(twopi*(x+y)) & 
     &            + dfg1(z,p11)*cos(twopi*(x-y)))

      if (za == 1) then
         f(1) = f(1)*0.54_dp*0.15_dp
         f(2) = f(2)*0.54_dp*0.15_dp
         f(3) = f(3)*0.54_dp*0.15_dp
      elseif (za == 6) then
         f(1) = f(1)*0.54_dp*0.40_dp
         f(2) = f(2)*0.54_dp*0.40_dp
         f(3) = f(3)*0.54_dp*0.40_dp
      endif

      end

!---------------------------------------------------------------------------

       function dfg1(x,p)
          use mod_precision


!
!    This is a function to evaluate the derivative of the variation of the
!     fourier coefficients with Z.
!

      implicit none
      real(dp) :: dfg1

!
!    Declare the simple variables.
!

      real(dp) :: x,y
      real(dp) :: a,n

      integer npar

!
!    Define the parameters.
!

      parameter(npar = 2)

!
!    Declare the arrays.
!

      real(dp) :: p(npar)

!
!    Evaluate the function
!

      a = p(1)
      n = p(2)**2
      y = 2.52_dp/x
      dfg1 = -(a*n/x)*(y**n)

      end

!---------------------------------------------------------------------------

       function dfg2(x,p)
          use mod_precision


!
!    This is a function to evaluate the derivative of the variation of the
!     fourier coefficients with Z.
!

      implicit none
      real(dp) :: dfg2

!
!    Declare the simple variables.
!

      real(dp) :: x,y
      real(dp) :: a,k,n1,n2

      integer npar

!
!    Define the parameters.
!

      parameter(npar = 4)

!
!    Declare the arrays.
!

      real(dp) :: p(npar)

!
!    Evaluate the function
!

      a = p(1)
      k = p(2)**2
      y = 1.0_dp/(k*x)
      n1 = p(3)**2
      n2 = p(4)**2
      dfg2 = -(a/x)*(n1*(y**n1) - n2*(y**n2))

      end

!---------------------------------------------------------------------------

       function lj(rin,za)
          use mod_precision


!
!    This is a function to evaluate the Lennard-Jones energy.
!

      implicit none
      real(dp) :: lj

!
!    Declare the simple variables.
!

      real(dp) :: a,b,magdr2,fac,sigma,epsilon

      integer l,lmax,m,mmax,n,nmax
      integer za,nbasis,i

!
!    Define parameters.
!

      parameter (lmax = 3)
      parameter (mmax = lmax)
      parameter (nmax = lmax)
      parameter (nbasis = 4)
      parameter (a = 2.52_dp)
      parameter (b = 3.57_dp)

!
!    Declare the arrays.
!

      real(dp) :: r(3)
      real(dp) :: rin(3)
      real(dp) :: dr(3)
      real(dp) :: d(3,nbasis)

!
!    Evaluate the potential.
!

      if (za == 1) then
         sigma = 3.783_dp
         epsilon = 0.013_dp*0.54_dp*0.15_dp
      elseif (za == 6) then
         sigma = 3.783_dp
         epsilon = 0.013_dp*0.54_dp*0.40_dp
      endif

      d(1,1) =  0.0_dp
      d(2,1) =  0.0_dp
      d(3,1) = -0.75_dp*b

      d(1,2) =  0.0_dp
      d(2,2) =  0.5_dp*a
      d(3,2) = -0.5_dp*b

      d(1,3) =  0.5_dp*a
      d(2,3) =  0.5_dp*a
      d(3,3) = -0.25_dp*b

      d(1,4) =  0.5_dp*a
      d(2,4) =  0.0_dp
      d(3,4) =  0.0_dp

      r(1) = mod(rin(1),a)
      if (r(1) < 0.0_dp) r(1) = r(1) + a
      r(2) = mod(rin(2),a)
      if (r(2) < 0.0_dp) r(2) = r(2) + a
      r(3) = rin(3)

      lj = 0.0_dp
      do l = -lmax,lmax,1
         do m = -mmax,mmax,1
            do n = -nmax,0,1
               do i = 1,nbasis,1
                  dr(1) = r(1) - real(l, dp)*a - d(1,i)
                  dr(2) = r(2) - real(m, dp)*a - d(2,i)
                  dr(3) = r(3) - real(n, dp)*b - d(3,i)
                  magdr2 = dr(1)**2 + dr(2)**2 + dr(3)**2
                  fac = (sigma*sigma/magdr2)**3
                  lj = lj + fac*fac - fac
               enddo
            enddo
         enddo
      enddo
      lj = 4.0_dp*epsilon*lj

      end

!---------------------------------------------------------------------------

      subroutine dlj(rin,za,f)
          use mod_precision


!
!    This is a function to evaluate the Lennard-Jones forces.
!

      implicit none

!
!    Declare the simple variables.
!

      real(dp) :: a,b,magdr2,fac,sigma,epsilon,lj

      integer l,lmax,m,mmax,n,nmax
      integer za,nbasis,i

!
!    Define parameters.
!

      parameter (lmax = 3)
      parameter (mmax = lmax)
      parameter (nmax = lmax)
      parameter (nbasis = 4)
      parameter (a = 2.52_dp)
      parameter (b = 3.57_dp)

!
!    Declare the arrays.
!

      real(dp) :: f(3)
      real(dp) :: r(3)
      real(dp) :: rin(3)
      real(dp) :: dr(3)
      real(dp) :: d(3,nbasis)

!
!    Evaluate the potential.
!

      if (za == 1) then
         sigma = 3.783_dp
         epsilon = 0.013_dp*0.54_dp*0.15_dp
      elseif (za == 6) then
         sigma = 3.783_dp
         epsilon = 0.013_dp*0.54_dp*0.40_dp
      endif

      d(1,1) =  0.0_dp
      d(2,1) =  0.0_dp
      d(3,1) = -0.75_dp*b

      d(1,2) =  0.0_dp
      d(2,2) =  0.5_dp*a
      d(3,2) = -0.5_dp*b

      d(1,3) =  0.5_dp*a
      d(2,3) =  0.5_dp*a
      d(3,3) = -0.25_dp*b

      d(1,4) =  0.5_dp*a
      d(2,4) =  0.0_dp
      d(3,4) =  0.0_dp

      r(1) = mod(rin(1),a)
      if (r(1) < 0.0_dp) r(1) = r(1) + a
      r(2) = mod(rin(2),a)
      if (r(2) < 0.0_dp) r(2) = r(2) + a
      r(3) = rin(3)

      f(1) = 0.0_dp
      f(2) = 0.0_dp
      f(3) = 0.0_dp

      do l = -lmax,lmax,1
         do m = -mmax,mmax,1
            do n = -nmax,0,1
               do i = 1,nbasis,1
                  dr(1) = r(1) - real(l, dp)*a - d(1,i)
                  dr(2) = r(2) - real(m, dp)*a - d(2,i)
                  dr(3) = r(3) - real(n, dp)*b - d(3,i)
                  magdr2 = dr(1)**2 + dr(2)**2 + dr(3)**2
                  fac = (sigma*sigma/magdr2)**3
                  lj = (2.0_dp*fac*fac - fac)/magdr2
                  f(1) = f(1) + lj*dr(1)
                  f(2) = f(2) + lj*dr(2)
                  f(3) = f(3) + lj*dr(3)
               enddo
            enddo
         enddo
      enddo
      f(1) = -24.0_dp*epsilon*f(1)
      f(2) = -24.0_dp*epsilon*f(2)
      f(3) = -24.0_dp*epsilon*f(3)

      end

!---------------------------------------------------------------------------

       function lj2(rdn,za,dxpln,dypln,dzpln)
          use mod_precision


!
!    This is a function to evaluate the potential for sliding planes
!

      implicit none
      real(dp) :: lj2

!
!    Declare simple variables.
!

      real(dp) :: dxpln,dypln,dzpln,lj

      integer za

!
!    Declare the arrays.
!

      real(dp) :: rdn(3)
      real(dp) :: rup(3)

!
!    Evaluate the energy
!

      rup(1) = rdn(1) - dxpln
      rup(2) = rdn(2) - dypln
      rup(3) = dzpln - rdn(3)
      lj2 = lj(rdn,za) + lj(rup,za)

      end

!---------------------------------------------------------------------------

      subroutine dlj2(rdn,za,dxpln,dypln,dzpln,dudn,duup)
          use mod_precision


!
!    This is a function to evaluate the potential for sliding planes
!

      implicit none

!
!    Declare simple variables.
!

      real(dp) :: dxpln,dypln,dzpln

      integer za

!
!    Declare the arrays.
!

      real(dp) :: rdn(3)
      real(dp) :: rup(3)
      real(dp) :: dudn(3)
      real(dp) :: duup(3)

!
!    Evaluate the force
!

      rup(1) = rdn(1) - dxpln
      rup(2) = rdn(2) - dypln
      rup(3) = dzpln - rdn(3)
      call dlj(rdn,za,dudn)
      call dlj(rup,za,duup)
      duup(3) = - duup(3)

      end

!---------------------------------------------------------------------------

       function bombv(r,r0)
          use mod_precision


!
!    This is a function to evaluate the spherical bomb potential.
!

      implicit none
      real(dp) :: bombv

!
!    Declare simple variables.
!

      real(dp) :: magr,delta,r0

!
!    Define the parameters.
!

      parameter (delta = 0.3_dp)

!
!    Declare the arrays.
!

      real(dp) :: r(3)

!
!    Evaluate the potential
!

      magr = sqrt(r(1)**2+r(2)**2+r(3)**2)
      if (magr <= r0) then
         bombv = 0.0_dp
      else
         bombv = ((magr-r0)/delta)**2
      endif

      end

!---------------------------------------------------------------------------

      subroutine dbombv(r,r0,f,dv)
          use mod_precision


!
!    This is a function to evaluate the derivative of the spherical
!     bomb potential.
!

      implicit none

!
!    Declare simple variables.
!

      real(dp) :: magr,delta,dv,r0

!
!    Define the parameters.
!

      parameter (delta = 0.3_dp)

!
!    Declare the arrays.
!

      real(dp) :: r(3)
      real(dp) :: f(3)

!
!    Evaluate the derivative of the potential
!

      magr = sqrt(r(1)**2+r(2)**2+r(3)**2)
      if (magr <= r0) then
         dv = 0.0_dp
      else
         dv = (2.0_dp/delta)*(magr-r0)/delta
      endif
      f(1) = dv*r(1)/magr
      f(2) = dv*r(2)/magr
      f(3) = dv*r(3)/magr

      end

!---------------------------------------------------------------------------

       function csrf2(rdn,za,dxpln,dypln,dzpln)
          use mod_precision


!
!    This is a function to evaluate the potential for sliding planes
!

      implicit none
      real(dp) :: csrf2

!
!    Declare simple variables.
!

      real(dp) :: dxpln,dypln,dzpln,csrf

      integer za

!
!    Declare the arrays.
!

      real(dp) :: rdn(3)
      real(dp) :: rup(3)

!
!    Evaluate the energy
!

      rup(1) = rdn(1) - dxpln
      rup(2) = rdn(2) - dypln
      rup(3) = dzpln - rdn(3)
      csrf2 = csrf(rdn,za) + csrf(rup,za)

      end

!---------------------------------------------------------------------------

      subroutine dcsrf2(rdn,za,dxpln,dypln,dzpln,dudn,duup)
          use mod_precision


!
!    This is a function to evaluate the potential for sliding planes
!

      implicit none

!
!    Declare simple variables.
!

      real(dp) :: dxpln,dypln,dzpln

      integer za

!
!    Declare the arrays.
!

      real(dp) :: rdn(3)
      real(dp) :: rup(3)
      real(dp) :: dudn(3)
      real(dp) :: duup(3)

!
!    Evaluate the force
!

      rup(1) = rdn(1) - dxpln
      rup(2) = rdn(2) - dypln
      rup(3) = dzpln - rdn(3)
      call dcsrf(rdn,za,dudn)
      call dcsrf(rup,za,duup)
      duup(3) = - duup(3)

      end

