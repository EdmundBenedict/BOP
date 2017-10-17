 
      subroutine mstmin(n,nmax,w,g1,g2,del,wrk,dspmax,tol,alpha,icom, &
     &                  gd1,alpha1,alpha2,alphau,alphal,nsrch,iprint)

     use mod_precision
      implicit double precision (a-h,o-z)
!  minimisation of a function of n coordinate variables. 
!  M.J.Norgett & R.Fletcher, J.Phys.C, 3, L190 (1972)
!  w      = inverse 2nd deriv. matrix supplied on first entry to mstmin
!           and updated at each subsequent iteration. On first entry
!           you may pass a guess (e.g. a unit matrix)  
!  g1     = first derivative at current iteration (first iteration only)
!           first derivative of previous iteration (subsequent iter)
!  g2     = first derivative at current iteration (after first iter)
!  del    = coordinate displacements calculated by mstmin
!  wrk    = work array
!  dspmax = step length in search direction
!  tol    = tolerance in minimisation (try 0.0000001, say)
!  alpha  = scaling factor to be applied in calling routine:
!           on each exit, the variables of the function being minimised
!           are updated by adding alpha * del(j) to the j'th variable
!  n      = first and second dimension of w (number of variables)
!  icom   = call indicator
!         = 0 for first entry (set in calling routine)
!         = 1 when relaxation is complete (returned by mstmin)
!  gd1,alpha1,alpha2,alphau,alphal,nsrch: parameters generated by mstmin
!  iprint = 0 for suppression of all printing
!         = 1 for some printing to logical unit iprint
      dimension w(nmax,n),g1(n),g2(n),del(n),wrk(n)
      save delmax
      zero=0e0_dp
      one=1e0_dp
      two=2e0_dp
      half=0.5_dp
      qrter=0.25_dp

      if (icom > 0) goto 10

      if (iprint > 0) write(iprint,100)
  100 format(/6x,'minimisation set-up.')
      gd1 = zero
      gd2 = zero
      delmax = zero
      do 2 i = 1,n
      del(i) = zero
      do 1 j = 1,n
      del(i) = del(i) - w(i,j) * g1(j)
    1 continue
      delmax = max(delmax,abs(del(i)))
      gd1 = gd1 + g1(i) * del(i)
    2 continue
      if (iprint > 0) write(iprint,101) gd1,gd2
  101 format(6x,'gd1 = ',f10.5,6x,'gd2 = ',f10.5)

      factor = min(one,dspmax/delmax)
      delmax = factor * delmax
      if (gd1 > zero) factor = - factor
      gd1 = factor * gd1
      do 3 i = 1,n
      del(i) = factor * del(i)
    3 continue
      if (iprint > 0) write(iprint,101) gd1,gd2
      alpha = one
      alpha1 = zero
      alpha2 = one
      if (iprint > 0) write(iprint,102) alpha,alpha1,alpha2
  102 format(6x,'alpha = ',f10.5,6x,'alpha1= ',f10.5,6x,'alpha2 = ', &
     &     f10.5)
      if(delmax < tol) goto 61
      icom = 2
      nsrch = 1
      return

   10 continue
      gd2 = zero
      do 11 i = 1,n
      gd2 = gd2 + g2(i) * del(i)
   11 continue
      if (iprint > 0) write(iprint,101) gd1,gd2

      if (abs(alpha) < 0.0001) goto 200
      if (abs(gd2/gd1) > qrter) goto 31

  200 if (iprint > 0) write(iprint,103)
  103 format(/6x,'change direction')
      alpha = alpha2 - alpha1
      dg = alpha * (gd2 - gd1)
      do 21 i = 1,n
      g1(i) = g2(i) - g1(i)
      del(i) = alpha * del(i)
   21 continue

      gwg = zero
      do 23 i = 1,n
      wrk(i) = zero
      do 22 j = 1,n
      wrk(i) = wrk(i) + w(i,j) * g1(j)
   22 continue
      gwg = gwg + g1(i) * wrk(i)
   23 continue

      do 24 i = 1,n
      g1(i) = g2(i)
      g2(i) = del(i)
      del(i) = zero
   24 continue

      q = one + gwg/dg
      gd1 = zero
      delmax = zero
      do 26 i = 1,n
      r = - g2(i)/dg
      s = (-wrk(i) + q * g2(i))/dg
      do 25 j = 1,n
      w(i,j) = w(i,j) + r * wrk(j) + s * g2(j)
      del(i) = del(i) - w(i,j) * g1(j)
   25 continue
      delmax = max(delmax,abs(del(i)))
      gd1 = gd1 + g1(i) * del(i)
   26 continue
      gd2 = zero
      if (iprint > 0) write(iprint,101) gd1,gd2

      if (delmax < tol) goto 61
      factor = min(one,dspmax/delmax)
      delmax = factor * delmax
      if (gd1 > zero) factor = - factor
      gd1 = factor * gd1
      do 27 i = 1,n
      del(i) = factor * del(i)
   27 continue
      if (iprint > 0) write(iprint,101) gd1,gd2
      alpha = one
      alpha1 = zero
      alpha2 = one
      if (iprint > 0) write(iprint,102) alpha,alpha1,alpha2
      icom = 2
      nsrch = 1
      return

   31 continue
      if (icom == 3) goto 51
      if (gd2 > gd1) goto 41

      if (iprint > 0) write(iprint,104)
  104 format(/6x,'negative curvature')
      nsrch = nsrch + 1
      if (nsrch > 5) goto 71

      gd1 = gd2
      do 32 i = 1,n
      g1(i) = g2(i)
   32 continue
      if (iprint > 0) write(iprint,101) gd1,gd2

      alpha = two * alpha
      if ((alpha * delmax) > dspmax) alpha = dspmax/delmax
      alpha1 = alpha2
      alpha2 = alpha1 + alpha
      if (iprint > 0) write(iprint,102) alpha,alpha1,alpha2
      return

   41 continue
      if (gd2 > zero) goto 45

      if (iprint > 0) write(iprint,105)
  105 format(/6x,'go further in this direction')
      alpha = (alpha1 - alpha2) * gd2/(gd2 - gd1)
      if ((alpha * delmax) > dspmax) alpha = dspmax/delmax
      alpha1 = alpha2
      alpha2 = alpha1 + alpha

      gd1 = gd2
      do 42 i = 1,n
      g1(i) = g2(i)
   42 continue
      if (iprint > 0) write(iprint,101) gd1,gd2
      if (iprint > 0) write(iprint,102) alpha,alpha1,alpha2
      return

   45 continue
      if (iprint > 0) write(iprint,106)
  106 format(/6x,'set bounds on minimum')
      alphal = alpha1
      alphau = alpha2
      if (iprint > 0) write(iprint,107) alphal,alphau
  107 format(6x,'lower bound = ',f10.5,10x,'upper bound = ',f10.5)
      icom = 3
      goto 55

   51 continue
      if (iprint > 0) write(iprint,108)
  108 format(/6x,'reset bounds')
      if (gd2 > zero) alphau = alpha2
      if (gd2 < zero) alphal = alpha2
      if (iprint > 0) write(iprint,107) alphal,alphau

   55 continue
      if (iprint > 0) write(iprint,109)
  109 format(/6x,'interpolate')
      alpha0 = ( - gd1 * alpha2 + gd2 * alpha1)/(gd2 - gd1)
      if (alpha0 > alphau) alpha0 = half * (alphau + alphal)
      if (alpha0 < alphal) alpha0 = half * (alphau + alphal)
      alpha = alpha0 - alpha2
      if (abs(gd2) > abs(gd1)) goto 57
      alpha1 = alpha2
      gd1 = gd2
      do 56 i = 1,n
      g1(i) = g2(i)
   56 continue
   57 continue
      alpha2 = alpha0
      if (iprint > 0) write(iprint,101) gd1,gd2
      if (iprint > 0) write(iprint,102) alpha,alpha1,alpha2
      return

   61 continue
      icom = 1
      if (iprint > 0) write(iprint,110)
  110 format(/6x,'valid minimisation.')
      return

   71 continue
      icom = 1
      write(*,111)
  111 format(/6x,'**** invalid minimisation - persistent positive curvature ****')
      end