 
      subroutine findzero(func,x1,x2,xacc,x,ifail)
          use mod_precision


!
!    This is a subroutine to find the zero of a smooth function.
!     The routines are based on the Numerical Recipes routines.
!

      implicit none

!
!    Declare the simple variables.
!

      real(dp) :: x1,x2,x,xacc
      real(dp) :: func,bisect,linint

      integer ifail

!
!    Declare the external functions.
!

      external func

!
!    Bracket the zero.
!
!     print *, x1, x2
      call bracket(func,x1,x2,ifail)

!
!    Refine the root using a binary section method.
!

      x = linint(func,x1,x2,xacc,ifail)
      if (ifail > 0) x = bisect(func,x1,x2,xacc,ifail)

      end

!-------------------------------------------------------------------------

      subroutine bracket(func,x1,x2,ifail)
          use mod_precision


      implicit none

      real(dp) :: x1,x2,factor,f1,f2,func

      integer ntry,j,ifail

      external func

      parameter (factor=1.6_dp)
      parameter (ntry=50)

      ifail = 0
      if(x1 == x2) then
         write(6,'(''BRACKET: You have to guess an initial range'')')
         ifail = 1
         return
      endif

      f1=func(x1)
      f2=func(x2)
      ifail = 0
      do j=1,ntry,1
         if(f1*f2 < 0.)return
         if(abs(f1) < abs(f2))then
           x1=x1+factor*(x1-x2)
           f1=func(x1)
         else
           x2=x2+factor*(x2-x1)
           f2=func(x2)
         endif
      enddo

      write(6,'(''BRACKET: Unable to bracket minimum.'')')
      write(6,'("X1 = ",F18.14,"       X2 = ",G18.14)') x1, x2
      ifail = 1

      end

!-------------------------------------------------------------------------

       function bisect(func,xlo,xhi,err,ifail)
          use mod_precision


      implicit none
      real(dp) :: bisect

      real(dp) :: xlo,xhi,xmid,err,func
      real(dp) :: fmid,flo,fhi

      integer jmax,j,ifail

      external func

      parameter (jmax=400)

      ifail = 0

      if (xlo > xhi) then
         xmid = xlo
         xlo = xhi
         xhi = xmid
      endif

      flo = func(xlo)
      fhi = func(xhi)
      if (flo*fhi >= 0.) then
         write(6,'(''BISECT: Root must be bracketed for bisection.'')')
         ifail = 1
         return
      endif

      fmid = flo
      j = 1
      do while((abs(xhi-xlo) > err).and.(abs(fmid) > 1.0e-10_dp) & 
     &         .and.(j <= jmax))
         xmid = (xlo+xhi)*0.5_dp
         fmid = func(xmid)
         if (fmid*flo > 0.0_dp) then
            xlo = xmid
            flo = fmid
         else
            xhi = xmid
            fhi = fmid
         endif
         j = j + 1
      enddo

      if (j > jmax) then
         write(6,'(''BISECT: Too many bisections'')')
         ifail = 1
      else
         bisect = xmid
      endif

      end

!-------------------------------------------------------------------------

       function linint(func,xlo,xhi,err,ifail)
          use mod_precision


      implicit none
      real(dp) :: linint

      real(dp) :: xlo,xhi,xmid,err,func
      real(dp) :: fmid,flo,fhi

      integer jmax,j,ifail

      external func

      parameter (jmax=30)

      ifail = 0

      if (xlo > xhi) then
         xmid = xlo
         xlo = xhi
         xhi = xmid
      endif

      flo = func(xlo)
      fhi = func(xhi)
      if (flo*fhi >= 0.) then
         write(6,'(''LININT: Root must be bracketed for'', & 
     &             '' interpolation.'')')
         ifail = 1
         return
      endif

      fmid = flo
      j = 1
      do while((abs(xhi-xlo) > err).and.(abs(fmid) > 1.0e-10_dp) & 
     &         .and.(j <= jmax))
         xmid = (xlo*fhi-xhi*flo)/(fhi-flo)
         fmid = func(xmid)
         if (fmid*flo > 0.0_dp) then
            xlo = xmid
            flo = fmid
         else
            xhi = xmid
            fhi = fmid
         endif
         j = j + 1
      enddo

      if (j > jmax) then
         write(6,'(''LININT: Too many bisections'')')
         ifail = 1
      else
         linint = xmid
      endif

      end
