 
   subroutine grdmat(grad,dr,z1,nstt1,z2,nstt2,scf,dsf,ham,ovl)

!
!    This is a routine to evaluate the gradients of the elements of the
!     tight-binding matrix.
!
!=============================================================================
!    Note that the angular momentum contributions are indexed as follows:
!
!		s		0
!
!		x		0
!		y		1
!		z		2
!
!		3zz-rr		0
!		xx-yy		1
!		xy		2
!		xz		3
!		yz		4
!
!=============================================================================
!
   use mod_precision
   use mod_const, only: mxnstat,ml,sqrt3
   use mod_conf, only: rtc,ham_conf_t,bond_conf_t
   use mod_atom_ar, only: zatype,asort,btype,bsort

   implicit none

   type(ham_conf_t),target :: ham
   type(bond_conf_t),pointer:: bnd
   
   integer, intent(in) :: z1,z2
   real(dp), intent(in) :: dr(3), scf(14), dsf(3,14)
   real(dp), intent(out) :: grad(3,mxnstat,mxnstat)
   real(dp) :: rfun(0:1,14)

   real(dp) :: vsss,vsps,vpss,vpps,vppp,vsds,vdss,vpds,vdps,vpdp,vdpp,vdds,vddp,vddd
   real(dp) :: dsss,dsps,dpss,dpps,dppp,dsds,ddss,dpds,ddps,dpdp,ddpp,ddds,dddp,dddd
   real(dp) :: magdr,l,m,n, lmn(3)
   real(dp) :: xl,xm,xn,yl,ym,yn,zl,zm,zn
   real(dp) :: xll,xlm,xln,xmm,xmn,xnn
   real(dp) :: yll,ylm,yln,ymm,ymn,ynn
   real(dp) :: zll,zlm,zln,zmm,zmn,znn
   real(dp) :: xlll,xllm,xlln,xlmm,xlmn,xlnn,xmmm,xmmn,xmnn,xnnn
   real(dp) :: ylll,yllm,ylln,ylmm,ylmn,ylnn,ymmm,ymmn,ymnn,ynnn
   real(dp) :: zlll,zllm,zlln,zlmm,zlmn,zlnn,zmmm,zmmn,zmnn,znnn
   real(dp) :: xllll,xlllm,xllln,xllmm,xllmn,xllnn,xlmmm,xlmmn
   real(dp) :: xlmnn,xlnnn,xmmmm,xmmmn,xmmnn,xmnnn,xnnnn
   real(dp) :: yllll,ylllm,yllln,yllmm,yllmn,yllnn,ylmmm,ylmmn
   real(dp) :: ylmnn,ylnnn,ymmmm,ymmmn,ymmnn,ymnnn,ynnnn
   real(dp) :: zllll,zlllm,zllln,zllmm,zllmn,zllnn,zlmmm,zlmmn
   real(dp) :: zlmnn,zlnnn,zmmmm,zmmmn,zmmnn,zmnnn,znnnn

   integer  :: nstt1,nstt2,nl1,nl2
   integer  :: i1,i2,j1,j2,i,j

   real(dp) :: dvsss(3),dvsps(3),dvpss(3),dvpps(3),dvppp(3)
   real(dp) :: dvsds(3),dvdss(3),dvpds(3),dvdps(3),dvpdp(3)
   real(dp) :: dvdpp(3),dvdds(3),dvddp(3),dvddd(3)

   integer :: llist1(ml)
   integer :: llist2(ml)
   
   logical :: ovl


!    print *, 'passed states', nstt1,nstt2
!    
!    call states(z1,nl1,nstt1,llist1)
!    call states(z2,nl2,nstt2,llist2)


      
      
!    print *, 'nl,nstt,llist'
!    print *, nl1,nstt1,llist1
!    print *, nl2,nstt2,llist2
   
   magdr = sqrt(dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3))



   if (magdr < 1.0e-6_dp) then

      grad(1:3, 1:nstt1, 1:nstt2) = 0.0_dp

   else
       
      nl1 = ham % a(asort(zatype(z1))) %b % norbs
      llist1 = ham % a(asort(zatype(z1))) %b % orblist
      
!       if (.not. dh) then
        nl2 = ham % a(asort(zatype(z2))) %b % norbs
        llist2 = ham % a(asort(zatype(z2))) %b % orblist
!       else
!         call states(z2,nl2,nstt2,llist2)
!       endif

      bnd => ham % b(bsort(btype(z1,z2)))
   
      call eval_rfun(z1, z2, magdr, 1, rfun,bnd,ovl)
      
      
      nullify(bnd)
!       
!       print *,'zero', rfun(0,:)
!       print *, 'one',rfun(1,:)


      lmn = dr/magdr
      
!             print *, lmn
      
      dvsss = rfun(1, 1) * lmn
      dvsps = rfun(1, 2) * lmn
      dvpss = rfun(1, 3) * lmn
      dvpps = rfun(1, 4) * lmn
      dvppp = rfun(1, 5) * lmn
      dvsds = rfun(1, 6) * lmn
      dvdss = rfun(1, 7) * lmn
      dvpds = rfun(1, 8) * lmn
      dvdps = rfun(1, 9) * lmn
      dvpdp = rfun(1,10) * lmn
      dvdpp = rfun(1,11) * lmn
      dvdds = rfun(1,12) * lmn
      dvddp = rfun(1,13) * lmn
      dvddd = rfun(1,14) * lmn
      
      vsss = rfun(0, 1)
      vsps = rfun(0, 2)
      vpss = rfun(0, 3)
      vpps = rfun(0, 4)
      vppp = rfun(0, 5)
      vsds = rfun(0, 6)
      vdss = rfun(0, 7)
      vpds = rfun(0, 8)
      vdps = rfun(0, 9)
      vpdp = rfun(0,10)
      vdpp = rfun(0,11)
      vdds = rfun(0,12)
      vddp = rfun(0,13)
      vddd = rfun(0,14)

      if (rtc%scf /= 0) then
         dvsss = dvsss * scf( 1) + vsss * dsf(:, 1)
         dvsps = dvsps * scf( 2) + vsps * dsf(:, 2)
         dvpss = dvpss * scf( 3) + vpss * dsf(:, 3)
         dvpps = dvpps * scf( 4) + vpps * dsf(:, 4)
         dvppp = dvppp * scf( 5) + vppp * dsf(:, 5)
         dvsds = dvsds * scf( 6) + vsds * dsf(:, 6)
         dvdss = dvdss * scf( 7) + vdss * dsf(:, 7)
         dvpds = dvpds * scf( 8) + vpds * dsf(:, 8)
         dvdps = dvdps * scf( 9) + vdps * dsf(:, 9)
         dvpdp = dvpdp * scf(10) + vpdp * dsf(:,10)
         dvdpp = dvdpp * scf(11) + vdpp * dsf(:,11)
         dvdds = dvdds * scf(12) + vdds * dsf(:,12)
         dvddp = dvddp * scf(13) + vddp * dsf(:,13)
         dvddd = dvddd * scf(14) + vddd * dsf(:,14)         
      
         vsss = vsss * scf( 1)
         vsps = vsps * scf( 2)
         vpss = vpss * scf( 3)
         vpps = vpps * scf( 4)
         vppp = vppp * scf( 5)
         vsds = vsds * scf( 6)
         vdss = vdss * scf( 7)
         vpds = vpds * scf( 8)
         vdps = vdps * scf( 9)
         vpdp = vpdp * scf(10)
         vdpp = vdpp * scf(11)
         vdds = vdds * scf(12)
         vddp = vddp * scf(13)
         vddd = vddd * scf(14)
         
!          print *, 'vscr',vdds,vddp,vddd
         
      end if

!         write(78,'("--------------------------------------------")')
!         write(78,'(4F8.4)') MAGDR,DR(1),DR(2),DR(3)     
!         write(78,'(''    DH      SCF    DSF1   DSF2   DSF3'')')
!         write(78,'(5F8.4)') DSSS, SCF(1), DVSSS(1), DVSSS(2), DVSSS(3)
!         write(78,'(5F8.4)') DSPS, SCF(2), DVSPS(1), DVSPS(2), DVSPS(3)
!         write(78,'(5F8.4)') DPSS, SCF(3), DVPSS(1), DVPSS(2), DVPSS(3)
!         write(78,'(5F8.4)') DPPS, SCF(4), DVPPS(1), DVPPS(2), DVPPS(3)
!         write(78,'(5F8.4)') DPPP, SCF(5), DVPPP(1), DVPPP(2), DVPPP(3)
!         write(78,'(5F8.4)') DSDS, SCF(6), DVSDS(1), DVSDS(2), DVSDS(3)
!         write(78,'(5F8.4)') DDSS, SCF(7), DVDSS(1), DVDSS(2), DVDSS(3)
!         write(78,'(5F8.4)') DPDS, SCF(8), DVPDS(1), DVPDS(2), DVPDS(3)
!         write(78,'(5F8.4)') DDPS, SCF(9), DVDPS(1), DVDPS(2), DVDPS(3)
!         write(78,'(5F8.4)') DPDP, SCF(10), DVPDP(1), DVPDP(2), DVPDP(3)
!         write(78,'(5F8.4)') DDPP, SCF(11), DVDPP(1), DVDPP(2), DVDPP(3)
!         write(78,'(5F8.4)') DDDS, SCF(12), DVDDS(1), DVDDS(2), DVDDS(3)
!         write(78,'(5F8.4)') DDDP, SCF(13), DVDDP(1), DVDDP(2), DVDDP(3)
!         write(78,'(5F8.4)') DDDD, SCF(14), DVDDD(1), DVDDD(2), DVDDD(3)


!
!       Set up the off-diagonal sub-matrix.
!

      l = lmn(1)
      m = lmn(2)
      n = lmn(3)

      if ((nstt1 > 4).or.(nstt2 > 4)) then

         xl = (1.0_dp - l*l)/magdr
         yl = -l*m/magdr
         zl = -l*n/magdr

         xm = -m*l/magdr
         ym = (1.0_dp - m*m)/magdr
         zm = -m*n/magdr

         xn = -n*l/magdr
         yn = -n*m/magdr
         zn = (1.0_dp - n*n)/magdr

         xll = 2*l*xl
         yll = 2*l*yl
         zll = 2*l*zl

         xmm = 2*m*xm
         ymm = 2*m*ym
         zmm = 2*m*zm

         xnn = 2*n*xn
         ynn = 2*n*yn
         znn = 2*n*zn

         xlm = xl*m + l*xm
         ylm = yl*m + l*ym
         zlm = zl*m + l*zm

         xln = xl*n + l*xn
         yln = yl*n + l*yn
         zln = zl*n + l*zn

         xmn = xm*n + m*xn
         ymn = ym*n + m*yn
         zmn = zm*n + m*zn

         xlll = 3*l*l*xl
         ylll = 3*l*l*yl
         zlll = 3*l*l*zl

         xllm = 2*l*m*xl + l*l*xm
         yllm = 2*l*m*yl + l*l*ym
         zllm = 2*l*m*zl + l*l*zm

         xlln = 2*l*n*xl + l*l*xn
         ylln = 2*l*n*yl + l*l*yn
         zlln = 2*l*n*zl + l*l*zn

         xlmm = 2*l*m*xm + m*m*xl
         ylmm = 2*l*m*ym + m*m*yl
         zlmm = 2*l*m*zm + m*m*zl

         xlmn = xl*m*n + l*xm*n + l*m*xn
         ylmn = yl*m*n + l*ym*n + l*m*yn
         zlmn = zl*m*n + l*zm*n + l*m*zn

         xlnn = 2*n*l*xn + n*n*xl
         ylnn = 2*n*l*yn + n*n*yl
         zlnn = 2*n*l*zn + n*n*zl

         xmmm = 3*m*m*xm
         ymmm = 3*m*m*ym
         zmmm = 3*m*m*zm

         xmmn = 2*m*n*xm + m*m*xn
         ymmn = 2*m*n*ym + m*m*yn
         zmmn = 2*m*n*zm + m*m*zn

         xmnn = 2*n*m*xn + n*n*xm
         ymnn = 2*n*m*yn + n*n*ym
         zmnn = 2*n*m*zn + n*n*zm

         xnnn = 3*n*n*xn
         ynnn = 3*n*n*yn
         znnn = 3*n*n*zn

         xllll = 4*l*l*l*xl
         yllll = 4*l*l*l*yl
         zllll = 4*l*l*l*zl

         xlllm = l*l*l*xm + m*xlll
         ylllm = l*l*l*ym + m*ylll
         zlllm = l*l*l*zm + m*zlll

         xllln = xlll*n + l*l*l*xn
         yllln = ylll*n + l*l*l*yn
         zllln = zlll*n + l*l*l*zn

         xllmm = l*l*xmm + m*m*xll
         yllmm = l*l*ymm + m*m*yll
         zllmm = l*l*zmm + m*m*zll

         xllmn = xll*m*n + l*l*xm*n + l*l*m*xn
         yllmn = yll*m*n + l*l*ym*n + l*l*m*yn
         zllmn = zll*m*n + l*l*zm*n + l*l*m*zn

         xllnn = xll*n*n + l*l*xnn
         yllnn = yll*n*n + l*l*ynn
         zllnn = zll*n*n + l*l*znn

         xlmmm = xl*m*m*m + l*xmmm
         ylmmm = yl*m*m*m + l*ymmm
         zlmmm = zl*m*m*m + l*zmmm

         xlmmn = xl*m*m*n + l*xmm*n + l*m*m*xn
         ylmmn = yl*m*m*n + l*ymm*n + l*m*m*yn
         zlmmn = zl*m*m*n + l*zmm*n + l*m*m*zn

         xlmnn = xl*m*n*n + l*xm*n*n + l*m*xnn
         ylmnn = yl*m*n*n + l*ym*n*n + l*m*ynn
         zlmnn = zl*m*n*n + l*zm*n*n + l*m*znn

         xlnnn = xl*n*n*n + l*xnnn
         ylnnn = yl*n*n*n + l*ynnn
         zlnnn = zl*n*n*n + l*znnn

         xmmmm = 4*m*m*m*xm
         ymmmm = 4*m*m*m*ym
         zmmmm = 4*m*m*m*zm

         xmmmn = xmmm*n + m*m*m*xn
         ymmmn = ymmm*n + m*m*m*yn
         zmmmn = zmmm*n + m*m*m*zn

         xmmnn = xmm*n*n + m*m*xnn
         ymmnn = ymm*n*n + m*m*ynn
         zmmnn = zmm*n*n + m*m*znn

         xmnnn = xm*n*n*n + m*xnnn
         ymnnn = ym*n*n*n + m*ynnn
         zmnnn = zm*n*n*n + m*znnn

         xnnnn = 4*n*n*n*xn
         ynnnn = 4*n*n*n*yn
         znnnn = 4*n*n*n*zn

      endif

      j1 = 1
      do i1 = 1,nl1

         if (llist1(i1) == 0) then

            j2 = 1
            do i2 = 1,nl2

               if (llist2(i2) == 0) then
! ss
                  grad(1,j1,j2) = dvsss(1)
                  grad(2,j1,j2) = dvsss(2)
                  grad(3,j1,j2) = dvsss(3)
                  j2 = j2 + 1

               elseif (llist2(i2) == 1) then
! sp
                  grad(1,j1,j2)   = l*(dvsps(1)-vsps*l/magdr) +vsps/magdr
                  grad(2,j1,j2)   = l*(dvsps(2)-vsps*m/magdr)
                  grad(3,j1,j2)   = l*(dvsps(3)-vsps*n/magdr)

                  grad(1,j1,j2+1) = m*(dvsps(1)-vsps*l/magdr)
                  grad(2,j1,j2+1) = m*(dvsps(2)-vsps*m/magdr) +vsps/magdr
                  grad(3,j1,j2+1) = m*(dvsps(3)-vsps*n/magdr)

                  grad(1,j1,j2+2) = n*(dvsps(1)-vsps*l/magdr)
                  grad(2,j1,j2+2) = n*(dvsps(2)-vsps*m/magdr)
                  grad(3,j1,j2+2) = n*(dvsps(3)-vsps*n/magdr) +vsps/magdr !!     BUG FIX 10/2/05 MJC ?
!                 grad(3,j1,j2+2) = n*(dvsps(3)-vsps*m/magdr) +vsps/magdr !!     BUG FIX 10/2/05 MJC ?
                  j2 = j2 + 3

               elseif (llist2(i2) == 2) then
! sd
                  grad(1,j1,j2)   = (n*n-0.5_dp*(l*l+m*m))*dvsds(1) + (xnn-0.5_dp*(xll+xmm))*vsds
                  grad(2,j1,j2)   = (n*n-0.5_dp*(l*l+m*m))*dvsds(2) + (ynn-0.5_dp*(yll+ymm))*vsds
                  grad(3,j1,j2)   = (n*n-0.5_dp*(l*l+m*m))*dvsds(3) + (znn-0.5_dp*(zll+zmm))*vsds

                  grad(1,j1,j2+1) = 0.5_dp*sqrt3*(l*l-m*m)*dvsds(1) + 0.5_dp*sqrt3*(xll-xmm)*vsds
                  grad(2,j1,j2+1) = 0.5_dp*sqrt3*(l*l-m*m)*dvsds(2) + 0.5_dp*sqrt3*(yll-ymm)*vsds
                  grad(3,j1,j2+1) = 0.5_dp*sqrt3*(l*l-m*m)*dvsds(3) + 0.5_dp*sqrt3*(zll-zmm)*vsds

                  grad(1,j1,j2+2) = sqrt3*l*m*dvsds(1) + sqrt3*xlm*vsds
                  grad(2,j1,j2+2) = sqrt3*l*m*dvsds(2) + sqrt3*ylm*vsds
                  grad(3,j1,j2+2) = sqrt3*l*m*dvsds(3) + sqrt3*zlm*vsds

                  grad(1,j1,j2+3) = sqrt3*l*n*dvsds(1) + sqrt3*xln*vsds
                  grad(2,j1,j2+3) = sqrt3*l*n*dvsds(2) + sqrt3*yln*vsds
                  grad(3,j1,j2+3) = sqrt3*l*n*dvsds(3) + sqrt3*zln*vsds

                  grad(1,j1,j2+4) = sqrt3*m*n*dvsds(1) + sqrt3*xmn*vsds
                  grad(2,j1,j2+4) = sqrt3*m*n*dvsds(2) + sqrt3*ymn*vsds
                  grad(3,j1,j2+4) = sqrt3*m*n*dvsds(3) + sqrt3*zmn*vsds

                  j2 = j2 + 5

               endif

            enddo

            j1 = j1 + 1

         elseif (llist1(i1) == 1) then

            j2 = 1
            do i2 = 1,nl2,1

               if (llist2(i2) == 0) then
! ps
                  grad(1,j1,j2)   = -l*(dvpss(1)-vpss*l/magdr) -vpss/magdr
                  grad(2,j1,j2)   = -l*(dvpss(2)-vpss*m/magdr)
                  grad(3,j1,j2)   = -l*(dvpss(3)-vpss*n/magdr)

                  grad(1,j1+1,j2) = -m*(dvpss(1)-vpss*l/magdr)
                  grad(2,j1+1,j2) = -m*(dvpss(2)-vpss*m/magdr) -vpss/magdr
                  grad(3,j1+1,j2) = -m*(dvpss(3)-vpss*n/magdr)

                  grad(1,j1+2,j2) = -n*(dvpss(1)-vpss*l/magdr)
                  grad(2,j1+2,j2) = -n*(dvpss(2)-vpss*m/magdr)
                  grad(3,j1+2,j2) = -n*(dvpss(3)-vpss*n/magdr) -vpss/magdr

                  j2 = j2 + 1

               elseif (llist2(i2) == 1) then
! pp
                  grad(1,j1,j2)     = l*l*(dvpps(1)-dvppp(1)-(2*l/magdr)*(vpps-vppp))+dvppp(1)+2.0_dp*l*(vpps-vppp)/magdr
                  grad(2,j1,j2)     = l*l*(dvpps(2)-dvppp(2)-(2*m/magdr)*(vpps-vppp))+dvppp(2)
                  grad(3,j1,j2)     = l*l*(dvpps(3)-dvppp(3)-(2*n/magdr)*(vpps-vppp))+dvppp(3)

                  grad(1,j1,j2+1)   = l*m*(dvpps(1)-dvppp(1)-(2*l/magdr)*(vpps-vppp))+m*(vpps-vppp)/magdr
                  grad(2,j1,j2+1)   = l*m*(dvpps(2)-dvppp(2)-(2*m/magdr)*(vpps-vppp))+l*(vpps-vppp)/magdr
                  grad(3,j1,j2+1)   = l*m*(dvpps(3)-dvppp(3)-(2*n/magdr)*(vpps-vppp))

                  grad(1,j1,j2+2)   = l*n*(dvpps(1)-dvppp(1)-(2*l/magdr)*(vpps-vppp))+n*(vpps-vppp)/magdr
                  grad(2,j1,j2+2)   = l*n*(dvpps(2)-dvppp(2)-(2*m/magdr)*(vpps-vppp))
                  grad(3,j1,j2+2)   = l*n*(dvpps(3)-dvppp(3)-(2*n/magdr)*(vpps-vppp))+l*(vpps-vppp)/magdr

                  grad(1,j1+1,j2)   = l*m*(dvpps(1)-dvppp(1)-(2*l/magdr)*(vpps-vppp))+m*(vpps-vppp)/magdr
                  grad(2,j1+1,j2)   = l*m*(dvpps(2)-dvppp(2)-(2*m/magdr)*(vpps-vppp))+l*(vpps-vppp)/magdr
                  grad(3,j1+1,j2)   = l*m*(dvpps(3)-dvppp(3)-(2*n/magdr)*(vpps-vppp))

                  grad(1,j1+1,j2+1) = m*m*(dvpps(1)-dvppp(1)-(2*l/magdr)*(vpps-vppp))+dvppp(1)
                  grad(2,j1+1,j2+1) = m*m*(dvpps(2)-dvppp(2)-(2*m/magdr)*(vpps-vppp))+dvppp(2)+2.0_dp*m*(vpps-vppp)/magdr
                  grad(3,j1+1,j2+1) = m*m*(dvpps(3)-dvppp(3)-(2*n/magdr)*(vpps-vppp))+dvppp(3)

                  grad(1,j1+1,j2+2) = m*n*(dvpps(1)-dvppp(1)-(2*l/magdr)*(vpps-vppp))
                  grad(2,j1+1,j2+2) = m*n*(dvpps(2)-dvppp(2)-(2*m/magdr)*(vpps-vppp))+n*(vpps-vppp)/magdr
                  grad(3,j1+1,j2+2) = m*n*(dvpps(3)-dvppp(3)-(2*n/magdr)*(vpps-vppp))+m*(vpps-vppp)/magdr

                  grad(1,j1+2,j2)   = l*n*(dvpps(1)-dvppp(1)-(2*l/magdr)*(vpps-vppp))+n*(vpps-vppp)/magdr
                  grad(2,j1+2,j2)   = l*n*(dvpps(2)-dvppp(2)-(2*m/magdr)*(vpps-vppp))
                  grad(3,j1+2,j2)   = l*n*(dvpps(3)-dvppp(3)-(2*n/magdr)*(vpps-vppp))+l*(vpps-vppp)/magdr

                  grad(1,j1+2,j2+1) = m*n*(dvpps(1)-dvppp(1)-(2*l/magdr)*(vpps-vppp))
                  grad(2,j1+2,j2+1) = m*n*(dvpps(2)-dvppp(2)-(2*m/magdr)*(vpps-vppp))+n*(vpps-vppp)/magdr
                  grad(3,j1+2,j2+1) = m*n*(dvpps(3)-dvppp(3)-(2*n/magdr)*(vpps-vppp))+m*(vpps-vppp)/magdr

                  grad(1,j1+2,j2+2) = n*n*(dvpps(1)-dvppp(1)-(2*l/magdr)*(vpps-vppp))+dvppp(1)
                  grad(2,j1+2,j2+2) = n*n*(dvpps(2)-dvppp(2)-(2*m/magdr)*(vpps-vppp))+dvppp(2)
                  grad(3,j1+2,j2+2) = n*n*(dvpps(3)-dvppp(3)-(2*n/magdr)*(vpps-vppp))+dvppp(3)+2.0_dp*n*(vpps-vppp)/magdr

                  j2 = j2 + 3

               elseif (llist2(i2) == 2) then
! pd
                  grad(1,j1,j2)   = (l*(n*n-0.5_dp*(l*l+m*m))*dvpds(1) & 
   &                                    - sqrt3*l*n*n*dvpdp(1)) & 
   &                                 + (xlnn-0.5_dp*(xlll+xlmm))*vpds & 
   &                                 - sqrt3*xlnn*vpdp
                  grad(2,j1,j2)   = (l*(n*n-0.5_dp*(l*l+m*m))*dvpds(2) & 
   &                                    - sqrt3*l*n*n*dvpdp(2)) & 
   &                                 + (ylnn-0.5_dp*(ylll+ylmm))*vpds & 
   &                                 - sqrt3*ylnn*vpdp
                  grad(3,j1,j2)   = (l*(n*n-0.5_dp*(l*l+m*m))*dvpds(3) & 
   &                                    - sqrt3*l*n*n*dvpdp(3)) & 
   &                                 + (zlnn-0.5_dp*(zlll+zlmm))*vpds & 
   &                                 - sqrt3*zlnn*vpdp

                  grad(1,j1,j2+1) = (0.5_dp*sqrt3*l*(l*l-m*m)*dvpds(1) & 
   &                                   + l*(1.0_dp-l*l+m*m)*dvpdp(1)) & 
   &                                 + 0.5_dp*sqrt3*(xlll-xlmm)*vpds & 
   &                                 + (xl-xlll+xlmm)*vpdp
                  grad(2,j1,j2+1) = (0.5_dp*sqrt3*l*(l*l-m*m)*dvpds(2) & 
   &                                   + l*(1.0_dp-l*l+m*m)*dvpdp(2)) & 
   &                                 + 0.5_dp*sqrt3*(ylll-ylmm)*vpds & 
   &                                 + (yl-ylll+ylmm)*vpdp
                  grad(3,j1,j2+1) = (0.5_dp*sqrt3*l*(l*l-m*m)*dvpds(3) & 
   &                                   + l*(1.0_dp-l*l+m*m)*dvpdp(3)) & 
   &                                 + 0.5_dp*sqrt3*(zlll-zlmm)*vpds & 
   &                                 + (zl-zlll+zlmm)*vpdp

                  grad(1,j1,j2+2)   = (sqrt3*l*l*m*dvpds(1) & 
   &                                 + m*(1.0_dp-2.0_dp*l*l)*dvpdp(1)) & 
   &                                 + sqrt3*xllm*vpds & 
   &                                 + (xm-2.0_dp*xllm)*vpdp
                  grad(2,j1,j2+2)   = (sqrt3*l*l*m*dvpds(2) & 
   &                                 + m*(1.0_dp-2.0_dp*l*l)*dvpdp(2)) & 
   &                                 + sqrt3*yllm*vpds & 
   &                                 + (ym-2.0_dp*yllm)*vpdp
                  grad(3,j1,j2+2)   = (sqrt3*l*l*m*dvpds(3) & 
   &                                 + m*(1.0_dp-2.0_dp*l*l)*dvpdp(3)) & 
   &                                 + sqrt3*zllm*vpds & 
   &                                 + (zm-2.0_dp*zllm)*vpdp

                  grad(1,j1,j2+3)   = (sqrt3*l*l*n*dvpds(1) & 
   &                                 + n*(1.0_dp-2.0_dp*l*l)*dvpdp(1)) & 
   &                                 + sqrt3*xlln*vpds & 
   &                                 + (xn-2.0_dp*xlln)*vpdp
                  grad(2,j1,j2+3)   = (sqrt3*l*l*n*dvpds(2) & 
   &                                 + n*(1.0_dp-2.0_dp*l*l)*dvpdp(2)) & 
   &                                 + sqrt3*ylln*vpds & 
   &                                 + (yn-2.0_dp*ylln)*vpdp
                  grad(3,j1,j2+3)   = (sqrt3*l*l*n*dvpds(3) & 
   &                                 + n*(1.0_dp-2.0_dp*l*l)*dvpdp(3)) & 
   &                                 + sqrt3*zlln*vpds & 
   &                                 + (zn-2.0_dp*zlln)*vpdp

                  grad(1,j1,j2+4)   = (sqrt3*l*m*n*dvpds(1) & 
   &                                    -2.0_dp*l*m*n*dvpdp(1)) & 
   &                                 + sqrt3*xlmn*vpds & 
   &                                 - 2.0_dp*xlmn*vpdp
                  grad(2,j1,j2+4)   = (sqrt3*l*m*n*dvpds(2) & 
   &                                    -2.0_dp*l*m*n*dvpdp(2)) & 
   &                                 + sqrt3*ylmn*vpds & 
   &                                 - 2.0_dp*ylmn*vpdp
                  grad(3,j1,j2+4)   = (sqrt3*l*m*n*dvpds(3) & 
   &                                    -2.0_dp*l*m*n*dvpdp(3)) & 
   &                                 + sqrt3*zlmn*vpds & 
   &                                 - 2.0_dp*zlmn*vpdp

                  grad(1,j1+1,j2) = (m*(n*n-0.5_dp*(l*l+m*m))*dvpds(1) & 
   &                                    - sqrt3*m*n*n*dvpdp(1)) & 
   &                                 + (xmnn-0.5_dp*(xllm+xmmm))*vpds & 
   &                                 - sqrt3*xmnn*vpdp
                  grad(2,j1+1,j2) = (m*(n*n-0.5_dp*(l*l+m*m))*dvpds(2) & 
   &                                    - sqrt3*m*n*n*dvpdp(2)) & 
   &                                 + (ymnn-0.5_dp*(yllm+ymmm))*vpds & 
   &                                 - sqrt3*ymnn*vpdp
                  grad(3,j1+1,j2) = (m*(n*n-0.5_dp*(l*l+m*m))*dvpds(3) & 
   &                                    - sqrt3*m*n*n*dvpdp(3)) & 
   &                                 + (zmnn-0.5_dp*(zllm+zmmm))*vpds & 
   &                                 - sqrt3*zmnn*vpdp

                  grad(1,j1+1,j2+1)=(0.5_dp*sqrt3*m*(l*l-m*m)*dvpds(1) & 
   &                                   - m*(1.0_dp+l*l-m*m)*dvpdp(1)) & 
   &                                 + 0.5_dp*sqrt3*(xllm-xmmm)*vpds & 
   &                                 - (xm+xllm-xmmm)*vpdp
                  grad(2,j1+1,j2+1)=(0.5_dp*sqrt3*m*(l*l-m*m)*dvpds(2) & 
   &                                   - m*(1.0_dp+l*l-m*m)*dvpdp(2)) & 
   &                                 + 0.5_dp*sqrt3*(yllm-ymmm)*vpds & 
   &                                 - (ym+yllm-ymmm)*vpdp
                  grad(3,j1+1,j2+1)=(0.5_dp*sqrt3*m*(l*l-m*m)*dvpds(3) & 
   &                                   - m*(1.0_dp+l*l-m*m)*dvpdp(3)) & 
   &                                 + 0.5_dp*sqrt3*(zllm-zmmm)*vpds & 
   &                                 - (zm+zllm-zmmm)*vpdp

                  grad(1,j1+1,j2+2) = (sqrt3*m*m*l*dvpds(1) & 
   &                                 + l*(1.0_dp-2.0_dp*m*m)*dvpdp(1)) & 
   &                                 + sqrt3*xlmm*vpds & 
   &                                 + (xl-2.0_dp*xlmm)*vpdp
                  grad(2,j1+1,j2+2) = (sqrt3*m*m*l*dvpds(2) & 
   &                                 + l*(1.0_dp-2.0_dp*m*m)*dvpdp(2)) & 
   &                                 + sqrt3*ylmm*vpds & 
   &                                 + (yl-2.0_dp*ylmm)*vpdp
                  grad(3,j1+1,j2+2) = (sqrt3*m*m*l*dvpds(3) & 
   &                                 + l*(1.0_dp-2.0_dp*m*m)*dvpdp(3)) & 
   &                                 + sqrt3*zlmm*vpds & 
   &                                 + (zl-2.0_dp*zlmm)*vpdp

                  grad(1,j1+1,j2+3) = (sqrt3*l*m*n*dvpds(1) & 
   &                                    - 2.0_dp*l*m*n*dvpdp(1)) & 
   &                                 + sqrt3*xlmn*vpds & 
   &                                 - 2.0_dp*xlmn*vpdp
                  grad(2,j1+1,j2+3) = (sqrt3*l*m*n*dvpds(2) & 
   &                                    - 2.0_dp*l*m*n*dvpdp(2)) & 
   &                                 + sqrt3*ylmn*vpds & 
   &                                 - 2.0_dp*ylmn*vpdp
                  grad(3,j1+1,j2+3) = (sqrt3*l*m*n*dvpds(3) & 
   &                                    - 2.0_dp*l*m*n*dvpdp(3)) & 
   &                                 + sqrt3*zlmn*vpds & 
   &                                 - 2.0_dp*zlmn*vpdp

                  grad(1,j1+1,j2+4) = (sqrt3*m*m*n*dvpds(1) & 
   &                                 + n*(1.0_dp-2.0_dp*m*m)*dvpdp(1)) & 
   &                                 + sqrt3*xmmn*vpds & 
   &                                 + (xn-2.0_dp*xmmn)*vpdp
                  grad(2,j1+1,j2+4) = (sqrt3*m*m*n*dvpds(2) & 
   &                                 + n*(1.0_dp-2.0_dp*m*m)*dvpdp(2)) & 
   &                                 + sqrt3*ymmn*vpds & 
   &                                 + (yn-2.0_dp*ymmn)*vpdp
                  grad(3,j1+1,j2+4) = (sqrt3*m*m*n*dvpds(3) & 
   &                                 + n*(1.0_dp-2.0_dp*m*m)*dvpdp(3)) & 
   &                                 + sqrt3*zmmn*vpds & 
   &                                 + (zn-2.0_dp*zmmn)*vpdp

                  grad(1,j1+2,j2) = (n*(n*n-0.5_dp*(l*l+m*m))*dvpds(1) & 
   &                                   + sqrt3*n*(l*l+m*m)*dvpdp(1)) & 
   &                                 + (xnnn-0.5_dp*(xlln+xmmn))*vpds & 
   &                                 + sqrt3*(xlln+xmmn)*vpdp
                  grad(2,j1+2,j2) = (n*(n*n-0.5_dp*(l*l+m*m))*dvpds(2) & 
   &                                   + sqrt3*n*(l*l+m*m)*dvpdp(2)) & 
   &                                 + (ynnn-0.5_dp*(ylln+ymmn))*vpds & 
   &                                 + sqrt3*(ylln+ymmn)*vpdp
                  grad(3,j1+2,j2) = (n*(n*n-0.5_dp*(l*l+m*m))*dvpds(3) & 
   &                                   + sqrt3*n*(l*l+m*m)*dvpdp(3)) & 
   &                                 + (znnn-0.5_dp*(zlln+zmmn))*vpds & 
   &                                 + sqrt3*(zlln+zmmn)*vpdp

                  grad(1,j1+2,j2+1)=(0.5_dp*sqrt3*n*(l*l-m*m)*dvpds(1) & 
   &                                    - n*(l*l-m*m)*dvpdp(1)) & 
   &                                 + 0.5_dp*sqrt3*(xlln-xmmn)*vpds & 
   &                                 - (xlln-xmmn)*vpdp
                  grad(2,j1+2,j2+1)=(0.5_dp*sqrt3*n*(l*l-m*m)*dvpds(2) & 
   &                                    - n*(l*l-m*m)*dvpdp(2)) & 
   &                                 + 0.5_dp*sqrt3*(ylln-ymmn)*vpds & 
   &                                 - (ylln-ymmn)*vpdp
                  grad(3,j1+2,j2+1)=(0.5_dp*sqrt3*n*(l*l-m*m)*dvpds(3) & 
   &                                    - n*(l*l-m*m)*dvpdp(3)) & 
   &                                 + 0.5_dp*sqrt3*(zlln-zmmn)*vpds & 
   &                                 - (zlln-zmmn)*vpdp

                  grad(1,j1+2,j2+2) = (sqrt3*l*m*n*dvpds(1) & 
   &                                    - 2.0_dp*l*m*n*dvpdp(1)) & 
   &                                 + sqrt3*xlmn*vpds & 
   &                                 - 2.0_dp*xlmn*vpdp
                  grad(2,j1+2,j2+2) = (sqrt3*l*m*n*dvpds(2) & 
   &                                    - 2.0_dp*l*m*n*dvpdp(2)) & 
   &                                 + sqrt3*ylmn*vpds & 
   &                                 - 2.0_dp*ylmn*vpdp
                  grad(3,j1+2,j2+2) = (sqrt3*l*m*n*dvpds(3) & 
   &                                    - 2.0_dp*l*m*n*dvpdp(3)) & 
   &                                 + sqrt3*zlmn*vpds & 
   &                                 - 2.0_dp*zlmn*vpdp

                  grad(1,j1+2,j2+3) = (sqrt3*n*n*l*dvpds(1) & 
   &                                + l*(1.0_dp-2.0_dp*n*n)*dvpdp(1)) & 
   &                                 + sqrt3*xlnn*vpds & 
   &                                 + (xl-2.0_dp*xlnn)*vpdp
                  grad(2,j1+2,j2+3) = (sqrt3*n*n*l*dvpds(2) & 
   &                                + l*(1.0_dp-2.0_dp*n*n)*dvpdp(2)) & 
   &                                 + sqrt3*ylnn*vpds & 
   &                                 + (yl-2.0_dp*ylnn)*vpdp
                  grad(3,j1+2,j2+3) = (sqrt3*n*n*l*dvpds(3) & 
   &                                + l*(1.0_dp-2.0_dp*n*n)*dvpdp(3)) & 
   &                                 + sqrt3*zlnn*vpds & 
   &                                 + (zl-2.0_dp*zlnn)*vpdp

                  grad(1,j1+2,j2+4) = (sqrt3*n*n*m*dvpds(1) & 
   &                                + m*(1.0_dp-2.0_dp*n*n)*dvpdp(1)) & 
   &                                 + sqrt3*xmnn*vpds & 
   &                                 + (xm-2.0_dp*xmnn)*vpdp
                  grad(2,j1+2,j2+4) = (sqrt3*n*n*m*dvpds(2) & 
   &                                + m*(1.0_dp-2.0_dp*n*n)*dvpdp(2)) & 
   &                                 + sqrt3*ymnn*vpds & 
   &                                 + (ym-2.0_dp*ymnn)*vpdp
                  grad(3,j1+2,j2+4) = (sqrt3*n*n*m*dvpds(3) & 
   &                                + m*(1.0_dp-2.0_dp*n*n)*dvpdp(3)) & 
   &                                 + sqrt3*zmnn*vpds & 
   &                                 + (zm-2.0_dp*zmnn)*vpdp

                  j2 = j2 + 5

               endif

            enddo

            j1 = j1 + 3

         elseif (llist1(i1) == 2) then

            j2 = 1
            do i2 = 1,nl2,1

               if (llist2(i2) == 0) then
! ds
                  grad(1,j1,j2)   = (n*n-0.5_dp*(l*l+m*m))*dvdss(1) & 
   &                               + (xnn-0.5_dp*(xll+xmm))*vdss
                  grad(2,j1,j2)   = (n*n-0.5_dp*(l*l+m*m))*dvdss(2) & 
   &                               + (ynn-0.5_dp*(yll+ymm))*vdss
                  grad(3,j1,j2)   = (n*n-0.5_dp*(l*l+m*m))*dvdss(3) & 
   &                               + (znn-0.5_dp*(zll+zmm))*vdss

                  grad(1,j1+1,j2) = 0.5_dp*sqrt3*(l*l-m*m)*dvdss(1) & 
   &                               + 0.5_dp*sqrt3*(xll-xmm)*vdss
                  grad(2,j1+1,j2) = 0.5_dp*sqrt3*(l*l-m*m)*dvdss(2) & 
   &                               + 0.5_dp*sqrt3*(yll-ymm)*vdss
                  grad(3,j1+1,j2) = 0.5_dp*sqrt3*(l*l-m*m)*dvdss(3) & 
   &                               + 0.5_dp*sqrt3*(zll-zmm)*vdss

                  grad(1,j1+2,j2) = sqrt3*l*m*dvdss(1) & 
   &                               + sqrt3*xlm*vdss
                  grad(2,j1+2,j2) = sqrt3*l*m*dvdss(2) & 
   &                               + sqrt3*ylm*vdss
                  grad(3,j1+2,j2) = sqrt3*l*m*dvdss(3) & 
   &                               + sqrt3*zlm*vdss

                  grad(1,j1+3,j2) = sqrt3*l*n*dvdss(1) & 
   &                               + sqrt3*xln*vdss
                  grad(2,j1+3,j2) = sqrt3*l*n*dvdss(2) & 
   &                               + sqrt3*yln*vdss
                  grad(3,j1+3,j2) = sqrt3*l*n*dvdss(3) & 
   &                               + sqrt3*zln*vdss

                  grad(1,j1+4,j2) = sqrt3*m*n*dvdss(1) & 
   &                               + sqrt3*xmn*vdss
                  grad(2,j1+4,j2) = sqrt3*m*n*dvdss(2) & 
   &                               + sqrt3*ymn*vdss
                  grad(3,j1+4,j2) = sqrt3*m*n*dvdss(3) & 
   &                               + sqrt3*zmn*vdss

                  j2 = j2 + 1

               elseif (llist2(i2) == 1) then
! dp
                  grad(1,j1,j2)  = -(l*(n*n-0.5_dp*(l*l+m*m))*dvdps(1) & 
   &                                    - sqrt3*l*n*n*dvdpp(1)) & 
   &                                 - (xlnn-0.5_dp*(xlll+xlmm))*vdps & 
   &                                 + sqrt3*xlnn*vdpp
                  grad(2,j1,j2)  = -(l*(n*n-0.5_dp*(l*l+m*m))*dvdps(2) & 
   &                                    - sqrt3*l*n*n*dvdpp(2)) & 
   &                                 - (ylnn-0.5_dp*(ylll+ylmm))*vdps & 
   &                                 + sqrt3*ylnn*vdpp
                  grad(3,j1,j2)  = -(l*(n*n-0.5_dp*(l*l+m*m))*dvdps(3) & 
   &                                    - sqrt3*l*n*n*dvdpp(3)) & 
   &                                 - (zlnn-0.5_dp*(zlll+zlmm))*vdps & 
   &                                 + sqrt3*zlnn*vdpp

                  grad(1,j1,j2+1)= -(m*(n*n-0.5_dp*(l*l+m*m))*dvdps(1) & 
   &                                    - sqrt3*m*n*n*dvdpp(1)) & 
   &                                 - (xmnn-0.5_dp*(xllm+xmmm))*vdps & 
   &                                 + sqrt3*xmnn*vdpp
                  grad(2,j1,j2+1)= -(m*(n*n-0.5_dp*(l*l+m*m))*dvdps(2) & 
   &                                    - sqrt3*m*n*n*dvdpp(2)) & 
   &                                 - (ymnn-0.5_dp*(yllm+ymmm))*vdps & 
   &                                 + sqrt3*ymnn*vdpp
                  grad(3,j1,j2+1)= -(m*(n*n-0.5_dp*(l*l+m*m))*dvdps(3) & 
   &                                    - sqrt3*m*n*n*dvdpp(3)) & 
   &                                 - (zmnn-0.5_dp*(zllm+zmmm))*vdps & 
   &                                 + sqrt3*zmnn*vdpp

                  grad(1,j1,j2+2)= -(n*(n*n-0.5_dp*(l*l+m*m))*dvdps(1) & 
   &                                  + sqrt3*n*(l*l+m*m)*dvdpp(1)) & 
   &                                 - (xnnn-0.5_dp*(xlln+xmmn))*vdps & 
   &                                 - sqrt3*(xlln+xmmn)*vdpp
                  grad(2,j1,j2+2)= -(n*(n*n-0.5_dp*(l*l+m*m))*dvdps(2) & 
   &                                  + sqrt3*n*(l*l+m*m)*dvdpp(2)) & 
   &                                 - (ynnn-0.5_dp*(ylln+ymmn))*vdps & 
   &                                 - sqrt3*(ylln+ymmn)*vdpp
                  grad(3,j1,j2+2)= -(n*(n*n-0.5_dp*(l*l+m*m))*dvdps(3) & 
   &                                  + sqrt3*n*(l*l+m*m)*dvdpp(3)) & 
   &                                 - (znnn-0.5_dp*(zlln+zmmn))*vdps & 
   &                                 - sqrt3*(zlln+zmmn)*vdpp

                  grad(1,j1+1,j2)= -(0.5_dp*sqrt3*l*(l*l-m*m)*dvdps(1) & 
   &                                  + l*(1.0_dp-l*l+m*m)*dvdpp(1)) & 
   &                                 - 0.5_dp*sqrt3*(xlll-xlmm)*vdps & 
   &                                 - (xl-xlll+xlmm)*vdpp
                  grad(2,j1+1,j2)= -(0.5_dp*sqrt3*l*(l*l-m*m)*dvdps(2) & 
   &                                  + l*(1.0_dp-l*l+m*m)*dvdpp(2)) & 
   &                                 - 0.5_dp*sqrt3*(ylll-ylmm)*vdps & 
   &                                 - (yl-ylll+ylmm)*vdpp
                  grad(3,j1+1,j2)= -(0.5_dp*sqrt3*l*(l*l-m*m)*dvdps(3) & 
   &                                  + l*(1.0_dp-l*l+m*m)*dvdpp(3)) & 
   &                                 - 0.5_dp*sqrt3*(zlll-zlmm)*vdps & 
   &                                 - (zl-zlll+zlmm)*vdpp

               grad(1,j1+1,j2+1) = -(0.5_dp*sqrt3*m*(l*l-m*m)*dvdps(1) & 
   &                                  - m*(1.0_dp+l*l-m*m)*dvdpp(1)) & 
   &                                 - 0.5_dp*sqrt3*(xllm-xmmm)*vdps & 
   &                                 + (xm+xllm-xmmm)*vdpp
               grad(2,j1+1,j2+1) = -(0.5_dp*sqrt3*m*(l*l-m*m)*dvdps(2) & 
   &                                  - m*(1.0_dp+l*l-m*m)*dvdpp(2)) & 
   &                                 - 0.5_dp*sqrt3*(yllm-ymmm)*vdps & 
   &                                 + (ym+yllm-ymmm)*vdpp
               grad(3,j1+1,j2+1) = -(0.5_dp*sqrt3*m*(l*l-m*m)*dvdps(3) & 
   &                                  - m*(1.0_dp+l*l-m*m)*dvdpp(3)) & 
   &                                 - 0.5_dp*sqrt3*(zllm-zmmm)*vdps & 
   &                                 + (zm+zllm-zmmm)*vdpp

               grad(1,j1+1,j2+2) = -(0.5_dp*sqrt3*n*(l*l-m*m)*dvdps(1) & 
   &                                    - n*(l*l-m*m)*dvdpp(1)) & 
   &                                 - 0.5_dp*sqrt3*(xlln-xmmn)*vdps & 
   &                                 + (xlln-xmmn)*vdpp
               grad(2,j1+1,j2+2) = -(0.5_dp*sqrt3*n*(l*l-m*m)*dvdps(2) & 
   &                                    - n*(l*l-m*m)*dvdpp(2)) & 
   &                                 - 0.5_dp*sqrt3*(ylln-ymmn)*vdps & 
   &                                 + (ylln-ymmn)*vdpp
               grad(3,j1+1,j2+2) = -(0.5_dp*sqrt3*n*(l*l-m*m)*dvdps(3) & 
   &                                    - n*(l*l-m*m)*dvdpp(3)) & 
   &                                 - 0.5_dp*sqrt3*(zlln-zmmn)*vdps & 
   &                                 + (zlln-zmmn)*vdpp

                  grad(1,j1+2,j2)   = -(sqrt3*l*l*m*dvdps(1) & 
   &                                + m*(1.0_dp-2.0_dp*l*l)*dvdpp(1)) & 
   &                                 - sqrt3*xllm*vdps & 
   &                                 - (xm-2.0_dp*xllm)*vdpp
                  grad(2,j1+2,j2)   = -(sqrt3*l*l*m*dvdps(2) & 
   &                                + m*(1.0_dp-2.0_dp*l*l)*dvdpp(2)) & 
   &                                 - sqrt3*yllm*vdps & 
   &                                 - (ym-2.0_dp*yllm)*vdpp
                  grad(3,j1+2,j2)   = -(sqrt3*l*l*m*dvdps(3) & 
   &                                + m*(1.0_dp-2.0_dp*l*l)*dvdpp(3)) & 
   &                                 - sqrt3*zllm*vdps & 
   &                                 - (zm-2.0_dp*zllm)*vdpp

                  grad(1,j1+2,j2+1) = -(sqrt3*m*m*l*dvdps(1) & 
   &                                + l*(1.0_dp-2.0_dp*m*m)*dvdpp(1)) & 
   &                                 - sqrt3*xlmm*vdps & 
   &                                 - (xl-2.0_dp*xlmm)*vdpp
                  grad(2,j1+2,j2+1) = -(sqrt3*m*m*l*dvdps(2) & 
   &                                + l*(1.0_dp-2.0_dp*m*m)*dvdpp(2)) & 
   &                                 - sqrt3*ylmm*vdps & 
   &                                 - (yl-2.0_dp*ylmm)*vdpp
                  grad(3,j1+2,j2+1) = -(sqrt3*m*m*l*dvdps(3) & 
   &                                + l*(1.0_dp-2.0_dp*m*m)*dvdpp(3)) & 
   &                                 - sqrt3*zlmm*vdps & 
   &                                 - (zl-2.0_dp*zlmm)*vdpp

                  grad(1,j1+2,j2+2) = -(sqrt3*l*m*n*dvdps(1) & 
   &                                    - 2.0_dp*l*m*n*dvdpp(1)) & 
   &                                 - sqrt3*xlmn*vdps & 
   &                                 + 2.0_dp*xlmn*vdpp
                  grad(2,j1+2,j2+2) = -(sqrt3*l*m*n*dvdps(2) & 
   &                                    - 2.0_dp*l*m*n*dvdpp(2)) & 
   &                                 - sqrt3*ylmn*vdps & 
   &                                 + 2.0_dp*ylmn*vdpp
                  grad(3,j1+2,j2+2) = -(sqrt3*l*m*n*dvdps(3) & 
   &                                    - 2.0_dp*l*m*n*dvdpp(3)) & 
   &                                 - sqrt3*zlmn*vdps & 
   &                                 + 2.0_dp*zlmn*vdpp

                  grad(1,j1+3,j2)   = -(sqrt3*l*l*n*dvdps(1) & 
   &                                + n*(1.0_dp-2.0_dp*l*l)*dvdpp(1)) & 
   &                                 - sqrt3*xlln*vdps & 
   &                                 - (xn-2.0_dp*xlln)*vdpp
                  grad(2,j1+3,j2)   = -(sqrt3*l*l*n*dvdps(2) & 
   &                                + n*(1.0_dp-2.0_dp*l*l)*dvdpp(2)) & 
   &                                 - sqrt3*ylln*vdps & 
   &                                 - (yn-2.0_dp*ylln)*vdpp
                  grad(3,j1+3,j2)   = -(sqrt3*l*l*n*dvdps(3) & 
   &                                + n*(1.0_dp-2.0_dp*l*l)*dvdpp(3)) & 
   &                                 - sqrt3*zlln*vdps & 
   &                                 - (zn-2.0_dp*zlln)*vdpp

                  grad(1,j1+3,j2+1) = -(sqrt3*l*m*n*dvdps(1) & 
   &                                    - 2.0_dp*l*m*n*dvdpp(1)) & 
   &                                 - sqrt3*xlmn*vdps & 
   &                                 + 2.0_dp*xlmn*vdpp
                  grad(2,j1+3,j2+1) = -(sqrt3*l*m*n*dvdps(2) & 
   &                                    - 2.0_dp*l*m*n*dvdpp(2)) & 
   &                                 - sqrt3*ylmn*vdps & 
   &                                 + 2.0_dp*ylmn*vdpp
                  grad(3,j1+3,j2+1) = -(sqrt3*l*m*n*dvdps(3) & 
   &                                    - 2.0_dp*l*m*n*dvdpp(3)) & 
   &                                 - sqrt3*zlmn*vdps & 
   &                                 + 2.0_dp*zlmn*vdpp

                  grad(1,j1+3,j2+2) = -(sqrt3*n*n*l*dvdps(1) & 
   &                                + l*(1.0_dp-2.0_dp*n*n)*dvdpp(1)) & 
   &                                 - sqrt3*xlnn*vdps & 
   &                                 - (xl-2.0_dp*xlnn)*vdpp
                  grad(2,j1+3,j2+2) = -(sqrt3*n*n*l*dvdps(2) & 
   &                                + l*(1.0_dp-2.0_dp*n*n)*dvdpp(2)) & 
   &                                 - sqrt3*ylnn*vdps & 
   &                                 - (yl-2.0_dp*ylnn)*vdpp
                  grad(3,j1+3,j2+2) = -(sqrt3*n*n*l*dvdps(3) & 
   &                                + l*(1.0_dp-2.0_dp*n*n)*dvdpp(3)) & 
   &                                 - sqrt3*zlnn*vdps & 
   &                                 - (zl-2.0_dp*zlnn)*vdpp

                  grad(1,j1+4,j2)   = -(sqrt3*l*m*n*dvdps(1) & 
   &                                    -2.0_dp*l*m*n*dvdpp(1)) & 
   &                                 - sqrt3*xlmn*vdps & 
   &                                 + 2.0_dp*xlmn*vdpp
                  grad(2,j1+4,j2)   = -(sqrt3*l*m*n*dvdps(2) & 
   &                                    -2.0_dp*l*m*n*dvdpp(2)) & 
   &                                 - sqrt3*ylmn*vdps & 
   &                                 + 2.0_dp*ylmn*vdpp
                  grad(3,j1+4,j2)   = -(sqrt3*l*m*n*dvdps(3) & 
   &                                    -2.0_dp*l*m*n*dvdpp(3)) & 
   &                                 - sqrt3*zlmn*vdps & 
   &                                 + 2.0_dp*zlmn*vdpp

                  grad(1,j1+4,j2+1) = -(sqrt3*m*m*n*dvdps(1) & 
   &                                + n*(1.0_dp-2.0_dp*m*m)*dvdpp(1)) & 
   &                                 - sqrt3*xmmn*vdps & 
   &                                 - (xn-2.0_dp*xmmn)*vdpp
                  grad(2,j1+4,j2+1) = -(sqrt3*m*m*n*dvdps(2) & 
   &                                + n*(1.0_dp-2.0_dp*m*m)*dvdpp(2)) & 
   &                                 - sqrt3*ymmn*vdps & 
   &                                 - (yn-2.0_dp*ymmn)*vdpp
                  grad(3,j1+4,j2+1) = -(sqrt3*m*m*n*dvdps(3) & 
   &                                + n*(1.0_dp-2.0_dp*m*m)*dvdpp(3)) & 
   &                                 - sqrt3*zmmn*vdps & 
   &                                 - (zn-2.0_dp*zmmn)*vdpp

                  grad(1,j1+4,j2+2) = -(sqrt3*n*n*m*dvdps(1) & 
   &                                + m*(1.0_dp-2.0_dp*n*n)*dvdpp(1)) & 
   &                                 - sqrt3*xmnn*vdps & 
   &                                 - (xm-2.0_dp*xmnn)*vdpp
                  grad(2,j1+4,j2+2) = -(sqrt3*n*n*m*dvdps(2) & 
   &                                + m*(1.0_dp-2.0_dp*n*n)*dvdpp(2)) & 
   &                                 - sqrt3*ymnn*vdps & 
   &                                 - (ym-2.0_dp*ymnn)*vdpp
                  grad(3,j1+4,j2+2) = -(sqrt3*n*n*m*dvdps(3) & 
   &                                + m*(1.0_dp-2.0_dp*n*n)*dvdpp(3)) & 
   &                                 - sqrt3*zmnn*vdps & 
   &                                 - (zm-2.0_dp*zmnn)*vdpp

                  j2 = j2 + 3

               elseif (llist2(i2) == 2) then
! dd
                  grad(1,j1,j2)     = (((n*n-0.5_dp*(l*l+m*m))**2)* & 
   &                            dvdds(1)+3.0_dp*n*n*(l*l+m*m)*dvddp(1) & 
   &                               +0.75_dp*((l*l+m*m)**2)*dvddd(1)) & 
   &                                 + 2.0_dp*(n*n-0.5_dp*(l*l+m*m))* & 
   &                                   (xnn-0.5_dp*(xll+xmm))*vdds & 
   &                                 + 3.0_dp*(xllnn+xmmnn)*vddp & 
   &                                 + 1.5_dp*(l*l+m*m)*(xll+xmm)*vddd
                  grad(2,j1,j2)     = (((n*n-0.5_dp*(l*l+m*m))**2)* & 
   &                            dvdds(2)+3.0_dp*n*n*(l*l+m*m)*dvddp(2) & 
   &                               +0.75_dp*((l*l+m*m)**2)*dvddd(2)) & 
   &                                 + 2.0_dp*(n*n-0.5_dp*(l*l+m*m))* & 
   &                                   (ynn-0.5_dp*(yll+ymm))*vdds & 
   &                                 + 3.0_dp*(yllnn+ymmnn)*vddp & 
   &                                 + 1.5_dp*(l*l+m*m)*(yll+ymm)*vddd
                  grad(3,j1,j2)     = (((n*n-0.5_dp*(l*l+m*m))**2)* & 
   &                            dvdds(3)+3.0_dp*n*n*(l*l+m*m)*dvddp(3) & 
   &                               +0.75_dp*((l*l+m*m)**2)*dvddd(3)) & 
   &                                 + 2.0_dp*(n*n-0.5_dp*(l*l+m*m))* & 
   &                                   (znn-0.5_dp*(zll+zmm))*vdds & 
   &                                 + 3.0_dp*(zllnn+zmmnn)*vddp & 
   &                                 + 1.5_dp*(l*l+m*m)*(zll+zmm)*vddd

                  grad(1,j1,j2+1)   = (0.5_dp*sqrt3*(l*l-m*m)* & 
   &                                    (n*n-0.5_dp*(l*l+m*m))*dvdds(1) & 
   &                                    + sqrt3*n*n*(m*m-l*l)*dvddp(1) & 
   &                                    + 0.25_dp*sqrt3*(1.0_dp+n*n)* & 
   &                                    (l*l-m*m)*dvddd(1)) & 
   &                                 + 0.5_dp*sqrt3*(xllnn-xmmnn- & 
   &                                   0.5_dp*(xllll-xmmmm))*vdds & 
   &                                 + sqrt3*(xmmnn-xllnn)*vddp & 
   &                                 + 0.25_dp*sqrt3*(xll-xmm+xllnn- & 
   &                                   xmmnn)*vddd
                  grad(2,j1,j2+1)   = (0.5_dp*sqrt3*(l*l-m*m)* & 
   &                                    (n*n-0.5_dp*(l*l+m*m))*dvdds(2) & 
   &                                    + sqrt3*n*n*(m*m-l*l)*dvddp(2) & 
   &                                    + 0.25_dp*sqrt3*(1.0_dp+n*n)* & 
   &                                    (l*l-m*m)*dvddd(2)) & 
   &                                 + 0.5_dp*sqrt3*(yllnn-ymmnn- & 
   &                                   0.5_dp*(yllll-ymmmm))*vdds & 
   &                                 + sqrt3*(ymmnn-yllnn)*vddp & 
   &                                 + 0.25_dp*sqrt3*(yll-ymm+yllnn- & 
   &                                   ymmnn)*vddd
                  grad(3,j1,j2+1)   = (0.5_dp*sqrt3*(l*l-m*m)* & 
   &                                    (n*n-0.5_dp*(l*l+m*m))*dvdds(3) & 
   &                                    + sqrt3*n*n*(m*m-l*l)*dvddp(3) & 
   &                                    + 0.25_dp*sqrt3*(1.0_dp+n*n)* & 
   &                                    (l*l-m*m)*dvddd(3)) & 
   &                                 + 0.5_dp*sqrt3*(zllnn-zmmnn- & 
   &                                   0.5_dp*(zllll-zmmmm))*vdds & 
   &                                 + sqrt3*(zmmnn-zllnn)*vddp & 
   &                                 + 0.25_dp*sqrt3*(zll-zmm+zllnn- & 
   &                                   zmmnn)*vddd

                  grad(1,j1,j2+2)   = (sqrt3*l*m*(n*n-0.5_dp* & 
   &                                 (l*l+m*m))*dvdds(1)-2.0_dp*sqrt3* & 
   &                                    l*m*n*n*dvddp(1)+0.5_dp*sqrt3* & 
   &                                    l*m*(1.0_dp+n*n)*dvddd(1)) & 
   &                                 + sqrt3*(xlmnn-0.5_dp*(xlllm+ & 
   &                                          xlmmm))*vdds & 
   &                                 - 2.0_dp*sqrt3*xlmnn*vddp & 
   &                                 + 0.5_dp*sqrt3*(xlm+xlmnn)*vddd
                  grad(2,j1,j2+2)   = (sqrt3*l*m*(n*n-0.5_dp* & 
   &                                 (l*l+m*m))*dvdds(2)-2.0_dp*sqrt3* & 
   &                                    l*m*n*n*dvddp(2)+0.5_dp*sqrt3* & 
   &                                    l*m*(1.0_dp+n*n)*dvddd(2)) & 
   &                                 + sqrt3*(ylmnn-0.5_dp*(ylllm+ & 
   &                                          ylmmm))*vdds & 
   &                                 - 2.0_dp*sqrt3*ylmnn*vddp & 
   &                                 + 0.5_dp*sqrt3*(ylm+ylmnn)*vddd
                  grad(3,j1,j2+2)   = (sqrt3*l*m*(n*n-0.5_dp* & 
   &                                 (l*l+m*m))*dvdds(3)-2.0_dp*sqrt3* & 
   &                                    l*m*n*n*dvddp(3)+0.5_dp*sqrt3* & 
   &                                    l*m*(1.0_dp+n*n)*dvddd(3)) & 
   &                                 + sqrt3*(zlmnn-0.5_dp*(zlllm+ & 
   &                                          zlmmm))*vdds & 
   &                                 - 2.0_dp*sqrt3*zlmnn*vddp & 
   &                                 + 0.5_dp*sqrt3*(zlm+zlmnn)*vddd

                  grad(1,j1,j2+3)   = (sqrt3*l*n*(n*n-0.5_dp* & 
   &                                   (l*l+m*m))*dvdds(1)+sqrt3* & 
   &                                    l*n*(l*l+m*m-n*n)*dvddp(1) & 
   &                                    -0.5_dp*sqrt3*l*n*(l*l+m*m)* & 
   &                                    dvddd(1)) & 
   &                                 + sqrt3*(xlnnn-0.5_dp*(xllln+ & 
   &                                          xlmmn))*vdds & 
   &                                 + sqrt3*(xllln+xlmmn-xlnnn)*vddp & 
   &                                 - 0.5_dp*sqrt3*(xllln+xlmmn)*vddd
                  grad(2,j1,j2+3)   = (sqrt3*l*n*(n*n-0.5_dp* & 
   &                                   (l*l+m*m))*dvdds(2)+sqrt3* & 
   &                                    l*n*(l*l+m*m-n*n)*dvddp(2) & 
   &                                    -0.5_dp*sqrt3*l*n*(l*l+m*m)* & 
   &                                    dvddd(2)) & 
   &                                 + sqrt3*(ylnnn-0.5_dp*(yllln+ & 
   &                                          ylmmn))*vdds & 
   &                                 + sqrt3*(yllln+ylmmn-ylnnn)*vddp & 
   &                                 - 0.5_dp*sqrt3*(yllln+ylmmn)*vddd
                  grad(3,j1,j2+3)   = (sqrt3*l*n*(n*n-0.5_dp* & 
   &                                   (l*l+m*m))*dvdds(3)+sqrt3* & 
   &                                    l*n*(l*l+m*m-n*n)*dvddp(3) & 
   &                                    -0.5_dp*sqrt3*l*n*(l*l+m*m)* & 
   &                                    dvddd(3)) & 
   &                                 + sqrt3*(zlnnn-0.5_dp*(zllln+ & 
   &                                          zlmmn))*vdds & 
   &                                 + sqrt3*(zllln+zlmmn-zlnnn)*vddp & 
   &                                 - 0.5_dp*sqrt3*(zllln+zlmmn)*vddd

                  grad(1,j1,j2+4)   = (sqrt3*m*n*(n*n-0.5_dp* & 
   &                                   (l*l+m*m))*dvdds(1)+sqrt3* & 
   &                                    m*n*(l*l+m*m-n*n)*dvddp(1) & 
   &                                    -0.5_dp*sqrt3*m*n*(l*l+m*m)* & 
   &                                    dvddd(1)) & 
   &                                 + sqrt3*(xmnnn-0.5_dp*(xllmn+ & 
   &                                          xmmmn))*vdds & 
   &                                 + sqrt3*(xllmn+xmmmn-xmnnn)*vddp & 
   &                                 - 0.5_dp*sqrt3*(xllmn+xmmmn)*vddd
                  grad(2,j1,j2+4)   = (sqrt3*m*n*(n*n-0.5_dp* & 
   &                                   (l*l+m*m))*dvdds(2)+sqrt3* & 
   &                                    m*n*(l*l+m*m-n*n)*dvddp(2) & 
   &                                    -0.5_dp*sqrt3*m*n*(l*l+m*m)* & 
   &                                    dvddd(2)) & 
   &                                 + sqrt3*(ymnnn-0.5_dp*(yllmn+ & 
   &                                          ymmmn))*vdds & 
   &                                 + sqrt3*(yllmn+ymmmn-ymnnn)*vddp & 
   &                                 - 0.5_dp*sqrt3*(yllmn+ymmmn)*vddd
                  grad(3,j1,j2+4)   = (sqrt3*m*n*(n*n-0.5_dp* & 
   &                                   (l*l+m*m))*dvdds(3)+sqrt3* & 
   &                                    m*n*(l*l+m*m-n*n)*dvddp(3) & 
   &                                    -0.5_dp*sqrt3*m*n*(l*l+m*m)* & 
   &                                    dvddd(3)) & 
   &                                 + sqrt3*(zmnnn-0.5_dp*(zllmn+ & 
   &                                          zmmmn))*vdds & 
   &                                 + sqrt3*(zllmn+zmmmn-zmnnn)*vddp & 
   &                                 - 0.5_dp*sqrt3*(zllmn+zmmmn)*vddd

                  grad(1,j1+1,j2)   = grad(1,j1,j2+1)
                  grad(2,j1+1,j2)   = grad(2,j1,j2+1)
                  grad(3,j1+1,j2)   = grad(3,j1,j2+1)

                  grad(1,j1+1,j2+1) = (0.75_dp*((l*l-m*m)**2)*dvdds(1) & 
   &                                 +(l*l+m*m-(l*l-m*m)**2)*dvddp(1) & 
   &                                   +(n*n+0.25_dp*(l*l-m*m)**2)* & 
   &                                   dvddd(1)) & 
   &                                 + 1.5_dp*(l*l-m*m)*(xll-xmm)*vdds & 
   &                                 + (xll+xmm-2.0_dp*(l*l-m*m)* & 
   &                                   (xll-xmm))*vddp & 
   &                                 + (xnn+0.5_dp*(l*l-m*m)* & 
   &                                   (xll-xmm))*vddd
                  grad(2,j1+1,j2+1) = (0.75_dp*((l*l-m*m)**2)*dvdds(2) & 
   &                                 +(l*l+m*m-(l*l-m*m)**2)*dvddp(2) & 
   &                                   +(n*n+0.25_dp*(l*l-m*m)**2)* & 
   &                                   dvddd(2)) & 
   &                                 + 1.5_dp*(l*l-m*m)*(yll-ymm)*vdds & 
   &                                 + (yll+ymm-2.0_dp*(l*l-m*m)* & 
   &                                   (yll-ymm))*vddp & 
   &                                 + (ynn+0.5_dp*(l*l-m*m)* & 
   &                                   (yll-ymm))*vddd
                  grad(3,j1+1,j2+1) = (0.75_dp*((l*l-m*m)**2)*dvdds(3) & 
   &                                 +(l*l+m*m-(l*l-m*m)**2)*dvddp(3) & 
   &                                   +(n*n+0.25_dp*(l*l-m*m)**2)* & 
   &                                   dvddd(3)) & 
   &                                 + 1.5_dp*(l*l-m*m)*(zll-zmm)*vdds & 
   &                                 + (zll+zmm-2.0_dp*(l*l-m*m)* & 
   &                                   (zll-zmm))*vddp & 
   &                                 + (znn+0.5_dp*(l*l-m*m)* & 
   &                                   (zll-zmm))*vddd

                  grad(1,j1+1,j2+2) = (1.5_dp*l*m*(l*l-m*m)*dvdds(1) & 
   &                                   +2.0_dp*l*m*(m*m-l*l)*dvddp(1) & 
   &                                 +0.5_dp*l*m*(l*l-m*m)*dvddd(1)) & 
   &                                 + 1.5_dp*(xlllm-xlmmm)*vdds & 
   &                                 + 2.0_dp*(xlmmm-xlllm)*vddp & 
   &                                 + 0.5_dp*(xlllm-xlmmm)*vddd
                  grad(2,j1+1,j2+2) = (1.5_dp*l*m*(l*l-m*m)*dvdds(2) & 
   &                                   +2.0_dp*l*m*(m*m-l*l)*dvddp(2) & 
   &                                 +0.5_dp*l*m*(l*l-m*m)*dvddd(2)) & 
   &                                 + 1.5_dp*(ylllm-ylmmm)*vdds & 
   &                                 + 2.0_dp*(ylmmm-ylllm)*vddp & 
   &                                 + 0.5_dp*(ylllm-ylmmm)*vddd
                  grad(3,j1+1,j2+2) = (1.5_dp*l*m*(l*l-m*m)*dvdds(3) & 
   &                                   +2.0_dp*l*m*(m*m-l*l)*dvddp(3) & 
   &                                 +0.5_dp*l*m*(l*l-m*m)*dvddd(3)) & 
   &                                 + 1.5_dp*(zlllm-zlmmm)*vdds & 
   &                                 + 2.0_dp*(zlmmm-zlllm)*vddp & 
   &                                 + 0.5_dp*(zlllm-zlmmm)*vddd

                  grad(1,j1+1,j2+3) = (1.5_dp*n*l*(l*l-m*m)*dvdds(1) & 
   &                                   +n*l*(1.0_dp-2.0_dp*(l*l-m*m))* & 
   &                                   dvddp(1) & 
   &                                   -n*l*(1.0_dp-0.5_dp*(l*l-m*m))* & 
   &                                   dvddd(1)) & 
   &                                 + 1.5_dp*(xllln-xlmmn)*vdds & 
   &                                 + (xln-2.0_dp*(xllln-xlmmn))*vddp & 
   &                                 - (xln-0.5_dp*(xllln-xlmmn))*vddd
                  grad(2,j1+1,j2+3) = (1.5_dp*n*l*(l*l-m*m)*dvdds(2) & 
   &                                   +n*l*(1.0_dp-2.0_dp*(l*l-m*m))* & 
   &                                   dvddp(2) & 
   &                                   -n*l*(1.0_dp-0.5_dp*(l*l-m*m))* & 
   &                                   dvddd(2)) & 
   &                                 + 1.5_dp*(yllln-ylmmn)*vdds & 
   &                                 + (yln-2.0_dp*(yllln-ylmmn))*vddp & 
   &                                 - (yln-0.5_dp*(yllln-ylmmn))*vddd
                  grad(3,j1+1,j2+3) = (1.5_dp*n*l*(l*l-m*m)*dvdds(3) & 
   &                                   +n*l*(1.0_dp-2.0_dp*(l*l-m*m))* & 
   &                                   dvddp(3) & 
   &                                   -n*l*(1.0_dp-0.5_dp*(l*l-m*m))* & 
   &                                   dvddd(3)) & 
   &                                 + 1.5_dp*(zllln-zlmmn)*vdds & 
   &                                 + (zln-2.0_dp*(zllln-zlmmn))*vddp & 
   &                                 - (zln-0.5_dp*(zllln-zlmmn))*vddd

                  grad(1,j1+1,j2+4) = (1.5_dp*m*n*(l*l-m*m)*dvdds(1) & 
   &                                   -m*n*(1.0_dp+2.0_dp*(l*l-m*m))* & 
   &                                   dvddp(1) & 
   &                                   +m*n*(1.0_dp+0.5_dp*(l*l-m*m))* & 
   &                                   dvddd(1)) & 
   &                                 + 1.5_dp*(xllmn-xmmmn)*vdds & 
   &                                 - (xmn+2.0_dp*(xllmn-xmmmn))*vddp & 
   &                                 + (xmn+0.5_dp*(xllmn-xmmmn))*vddd
                  grad(2,j1+1,j2+4) = (1.5_dp*m*n*(l*l-m*m)*dvdds(2) & 
   &                                   -m*n*(1.0_dp+2.0_dp*(l*l-m*m))* & 
   &                                   dvddp(2) & 
   &                                   +m*n*(1.0_dp+0.5_dp*(l*l-m*m))* & 
   &                                   dvddd(2)) & 
   &                                 + 1.5_dp*(yllmn-ymmmn)*vdds & 
   &                                 - (ymn+2.0_dp*(yllmn-ymmmn))*vddp & 
   &                                 + (ymn+0.5_dp*(yllmn-ymmmn))*vddd
                  grad(3,j1+1,j2+4) = (1.5_dp*m*n*(l*l-m*m)*dvdds(3) & 
   &                                   -m*n*(1.0_dp+2.0_dp*(l*l-m*m))* & 
   &                                   dvddp(3) & 
   &                                   +m*n*(1.0_dp+0.5_dp*(l*l-m*m))* & 
   &                                   dvddd(3)) & 
   &                                 + 1.5_dp*(zllmn-zmmmn)*vdds & 
   &                                 - (zmn+2.0_dp*(zllmn-zmmmn))*vddp & 
   &                                 + (zmn+0.5_dp*(zllmn-zmmmn))*vddd

                  grad(1,j1+2,j2)   = grad(1,j1,j2+2)
                  grad(2,j1+2,j2)   = grad(2,j1,j2+2)
                  grad(3,j1+2,j2)   = grad(3,j1,j2+2)

                  grad(1,j1+2,j2+1) = grad(1,j1+1,j2+2)
                  grad(2,j1+2,j2+1) = grad(2,j1+1,j2+2)
                  grad(3,j1+2,j2+1) = grad(3,j1+1,j2+2)

                  grad(1,j1+2,j2+2) = (3.0_dp*l*l*m*m*dvdds(1) & 
   &                                 +(l*l+m*m-4.0_dp*l*l*m*m)*dvddp(1) & 
   &                                   +(n*n+l*l*m*m)*dvddd(1)) & 
   &                                 + 3.0_dp*xllmm*vdds & 
   &                                 + (xll+xmm-4.0_dp*xllmm)*vddp & 
   &                                 + (xnn+xllmm)*vddd
                  grad(2,j1+2,j2+2) = (3.0_dp*l*l*m*m*dvdds(2) & 
   &                                 +(l*l+m*m-4.0_dp*l*l*m*m)*dvddp(2) & 
   &                                   +(n*n+l*l*m*m)*dvddd(2)) & 
   &                                 + 3.0_dp*yllmm*vdds & 
   &                                 + (yll+ymm-4.0_dp*yllmm)*vddp & 
   &                                 + (ynn+yllmm)*vddd
                  grad(3,j1+2,j2+2) = (3.0_dp*l*l*m*m*dvdds(3) & 
   &                                 +(l*l+m*m-4.0_dp*l*l*m*m)*dvddp(3) & 
   &                                   +(n*n+l*l*m*m)*dvddd(3)) & 
   &                                 + 3.0_dp*zllmm*vdds & 
   &                                 + (zll+zmm-4.0_dp*zllmm)*vddp & 
   &                                 + (znn+zllmm)*vddd

                  grad(1,j1+2,j2+3) = (3.0_dp*l*l*m*n*dvdds(1) & 
   &                                   +m*n*(1.0_dp-4.0_dp*l*l)*dvddp(1) & 
   &                                   +m*n*(l*l-1.0_dp)*dvddd(1)) & 
   &                                 + 3.0_dp*xllmn*vdds & 
   &                                 + (xmn-4.0_dp*xllmn)*vddp & 
   &                                 + (xllmn-xmn)*vddd
                  grad(2,j1+2,j2+3) = (3.0_dp*l*l*m*n*dvdds(2) & 
   &                                   +m*n*(1.0_dp-4.0_dp*l*l)*dvddp(2) & 
   &                                   +m*n*(l*l-1.0_dp)*dvddd(2)) & 
   &                                 + 3.0_dp*yllmn*vdds & 
   &                                 + (ymn-4.0_dp*yllmn)*vddp & 
   &                                 + (yllmn-ymn)*vddd
                  grad(3,j1+2,j2+3) = (3.0_dp*l*l*m*n*dvdds(3) & 
   &                                   +m*n*(1.0_dp-4.0_dp*l*l)*dvddp(3) & 
   &                                   +m*n*(l*l-1.0_dp)*dvddd(3)) & 
   &                                 + 3.0_dp*zllmn*vdds & 
   &                                 + (zmn-4.0_dp*zllmn)*vddp & 
   &                                 + (zllmn-zmn)*vddd

                  grad(1,j1+2,j2+4) = (3.0_dp*l*m*m*n*dvdds(1) & 
   &                                   +l*n*(1.0_dp-4.0_dp*m*m)*dvddp(1) & 
   &                                   +l*n*(m*m-1.0_dp)*dvddd(1)) & 
   &                                 + 3.0_dp*xlmmn*vdds & 
   &                                 + (xln-4.0_dp*xlmmn)*vddp & 
   &                                 + (xlmmn-xln)*vddd
                  grad(2,j1+2,j2+4) = (3.0_dp*l*m*m*n*dvdds(2) & 
   &                                   +l*n*(1.0_dp-4.0_dp*m*m)*dvddp(2) & 
   &                                   +l*n*(m*m-1.0_dp)*dvddd(2)) & 
   &                                 + 3.0_dp*ylmmn*vdds & 
   &                                 + (yln-4.0_dp*ylmmn)*vddp & 
   &                                 + (ylmmn-yln)*vddd
                  grad(3,j1+2,j2+4) = (3.0_dp*l*m*m*n*dvdds(3) & 
   &                                   +l*n*(1.0_dp-4.0_dp*m*m)*dvddp(3) & 
   &                                   +l*n*(m*m-1.0_dp)*dvddd(3)) & 
   &                                 + 3.0_dp*zlmmn*vdds & 
   &                                 + (zln-4.0_dp*zlmmn)*vddp & 
   &                                 + (zlmmn-zln)*vddd

                  grad(1,j1+3,j2)   = grad(1,j1,j2+3)
                  grad(2,j1+3,j2)   = grad(2,j1,j2+3)
                  grad(3,j1+3,j2)   = grad(3,j1,j2+3)

                  grad(1,j1+3,j2+1) = grad(1,j1+1,j2+3)
                  grad(2,j1+3,j2+1) = grad(2,j1+1,j2+3)
                  grad(3,j1+3,j2+1) = grad(3,j1+1,j2+3)

                  grad(1,j1+3,j2+2) = grad(1,j1+2,j2+3)
                  grad(2,j1+3,j2+2) = grad(2,j1+2,j2+3)
                  grad(3,j1+3,j2+2) = grad(3,j1+2,j2+3)

                  grad(1,j1+3,j2+3) = (3.0_dp*l*l*n*n*dvdds(1) & 
   &                                 +(l*l+n*n-4.0_dp*l*l*n*n)*dvddp(1) & 
   &                                   +(m*m+l*l*n*n)*dvddd(1)) & 
   &                                 + 3.0_dp*xllnn*vdds & 
   &                                 + (xll+xnn-4.0_dp*xllnn)*vddp & 
   &                                 + (xmm+xllnn)*vddd
                  grad(2,j1+3,j2+3) = (3.0_dp*l*l*n*n*dvdds(2) & 
   &                                 +(l*l+n*n-4.0_dp*l*l*n*n)*dvddp(2) & 
   &                                   +(m*m+l*l*n*n)*dvddd(2)) & 
   &                                 + 3.0_dp*yllnn*vdds & 
   &                                 + (yll+ynn-4.0_dp*yllnn)*vddp & 
   &                                 + (ymm+yllnn)*vddd
                  grad(3,j1+3,j2+3) = (3.0_dp*l*l*n*n*dvdds(3) & 
   &                                 +(l*l+n*n-4.0_dp*l*l*n*n)*dvddp(3) & 
   &                                   +(m*m+l*l*n*n)*dvddd(3)) & 
   &                                 + 3.0_dp*zllnn*vdds & 
   &                                 + (zll+znn-4.0_dp*zllnn)*vddp & 
   &                                 + (zmm+zllnn)*vddd

                  grad(1,j1+3,j2+4) = (3.0_dp*l*m*n*n*dvdds(1) & 
   &                                   +l*m*(1.0_dp-4.0_dp*n*n)*dvddp(1) & 
   &                                   +l*m*(n*n-1.0_dp)*dvddd(1)) & 
   &                                 + 3.0_dp*xlmnn*vdds & 
   &                                 + (xlm-4.0_dp*xlmnn)*vddp & 
   &                                 + (xlmnn-xlm)*vddd
                  grad(2,j1+3,j2+4) = (3.0_dp*l*m*n*n*dvdds(2) & 
   &                                   +l*m*(1.0_dp-4.0_dp*n*n)*dvddp(2) & 
   &                                   +l*m*(n*n-1.0_dp)*dvddd(2)) & 
   &                                 + 3.0_dp*ylmnn*vdds & 
   &                                 + (ylm-4.0_dp*ylmnn)*vddp & 
   &                                 + (ylmnn-ylm)*vddd
                  grad(3,j1+3,j2+4) = (3.0_dp*l*m*n*n*dvdds(3) & 
   &                                   +l*m*(1.0_dp-4.0_dp*n*n)*dvddp(3) & 
   &                                   +l*m*(n*n-1.0_dp)*dvddd(3)) & 
   &                                 + 3.0_dp*zlmnn*vdds & 
   &                                 + (zlm-4.0_dp*zlmnn)*vddp & 
   &                                 + (zlmnn-zlm)*vddd

                  grad(1,j1+4,j2)   = grad(1,j1,j2+4)
                  grad(2,j1+4,j2)   = grad(2,j1,j2+4)
                  grad(3,j1+4,j2)   = grad(3,j1,j2+4)

                  grad(1,j1+4,j2+1) = grad(1,j1+1,j2+4)
                  grad(2,j1+4,j2+1) = grad(2,j1+1,j2+4)
                  grad(3,j1+4,j2+1) = grad(3,j1+1,j2+4)

                  grad(1,j1+4,j2+2) = grad(1,j1+2,j2+4)
                  grad(2,j1+4,j2+2) = grad(2,j1+2,j2+4)
                  grad(3,j1+4,j2+2) = grad(3,j1+2,j2+4)

                  grad(1,j1+4,j2+3) = grad(1,j1+3,j2+4)
                  grad(2,j1+4,j2+3) = grad(2,j1+3,j2+4)
                  grad(3,j1+4,j2+3) = grad(3,j1+3,j2+4)

                  grad(1,j1+4,j2+4) = (3.0_dp*n*n*m*m*dvdds(1) & 
   &                                 +(n*n+m*m-4.0_dp*n*n*m*m)*dvddp(1) & 
   &                                   +(l*l+n*n*m*m)*dvddd(1)) & 
   &                                 + 3.0_dp*xmmnn*vdds & 
   &                                 + (xnn+xmm-4.0_dp*xmmnn)*vddp & 
   &                                 + (xll+xmmnn)*vddd
                  grad(2,j1+4,j2+4) = (3.0_dp*n*n*m*m*dvdds(2) & 
   &                                 +(n*n+m*m-4.0_dp*n*n*m*m)*dvddp(2) & 
   &                                   +(l*l+n*n*m*m)*dvddd(2)) & 
   &                                 + 3.0_dp*ymmnn*vdds & 
   &                                 + (ynn+ymm-4.0_dp*ymmnn)*vddp & 
   &                                 + (yll+ymmnn)*vddd
                  grad(3,j1+4,j2+4) = (3.0_dp*n*n*m*m*dvdds(3) & 
   &                                 +(n*n+m*m-4.0_dp*n*n*m*m)*dvddp(3) & 
   &                                   +(l*l+n*n*m*m)*dvddd(3)) & 
   &                                 + 3.0_dp*zmmnn*vdds & 
   &                                 + (znn+zmm-4.0_dp*zmmnn)*vddp & 
   &                                 + (zll+zmmnn)*vddd

                  j2 = j2 + 5
               endif
            enddo
            j1 = j1 + 3
         endif
      enddo
   endif

!    do i = 1,nstt1
!       do j = 1, nstt2
!           print *, 'grad', grad(:,i,j)
!       enddo
!    enddo
   
!      DO J1=1,5
!         WRITE(80,'(5("(",3F8.4,")",3X)')
!     +                    GRAD(1,J1,1),GRAD(2,J1,1),GRAD(3,J1,1),
!     +                    GRAD(1,J1,2),GRAD(2,J1,2),GRAD(3,J1,2),
!     +                    GRAD(1,J1,3),GRAD(2,J1,3),GRAD(3,J1,3),
!     +                    GRAD(1,J1,4),GRAD(2,J1,4),GRAD(3,J1,4),
!     +                    GRAD(1,J1,5),GRAD(2,J1,5),GRAD(3,J1,5)
!      ENDDO

   end subroutine grdmat

