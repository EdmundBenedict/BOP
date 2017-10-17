 
      subroutine bndinf(z1,z2,dr,mxnstat,ml,llist1,llist2, &
     &                  bnddat,bndscl,btype,nbtype,mbtype,bo, & 
     &                  vsss,vsps,vpss,vpps,vppp,vsds,vdss, & 
     &                  vpds,vdps,vpdp,vdpp,vdds,vddp,vddd, & 
     &                  bosss,bosps,bopss,bopps,boppp,bosds, & 
     &                  bodss,bopds,bodps,bopdp,bodpp,bodds, & 
     &                  boddp,boddd,besss,besps,bepss,bepps, & 
     &                  beppp,besds,bedss,bepds,bedps,bepdp, & 
     &                  bedpp,bedds,beddp,beddd)

          use mod_precision
!
!    This is a routine to evaluate the details of the bonding for one bond
!     in the frame of reference in which the z axis lies along the bond.
!

      implicit none

!
!    Declare the simple variables.
!

      real(dp) :: magdr
      real(dp) :: vsss,vsps,vpss,vpps,vppp
      real(dp) :: vsds,vdss,vpds,vdps,vpdp,vdpp
      real(dp) :: vdds,vddp,vddd
      real(dp) :: bosss,bosps,bopss,bopps,boppp
      real(dp) :: bosds,bodss,bopds,bodps,bopdp,bodpp
      real(dp) :: bodds,boddp,boddd
      real(dp) :: besss,besps,bepss,bepps,beppp
      real(dp) :: besds,bedss,bepds,bedps,bepdp,bedpp
      real(dp) :: bedds,beddp,beddd

      integer z1,z2
      integer mxnstat,mxz,ml
      integer nstt1,nstt2,nl1,nl2
      integer nbtype,mbtype
      integer il1,il2,l1,l2,lm1,lm2
      integer jl1,jl2,m1,m2
      integer i,j

!
!    Define the parameters.
!

      parameter (mxz = 103)

!
!    Declare the arrays.
!

      real(dp) :: dr(3)
      real(dp) :: bnddat(15,mbtype)
      real(dp) :: bndscl(14,15,mbtype)
      real(dp) :: bo(mxnstat,mxnstat)
      real(dp) :: fullbo(9,9)
      real(dp) :: rbo(9,9)

      integer llist1(ml)
      integer llist2(ml)
      integer btype(0:mxz,0:mxz)

!
!    Find some parameters.
!

      call states(z1,nl1,nstt1,llist1)
      call states(z2,nl2,nstt2,llist2)

!
!    Expand the bond order matrix out into the full 9x9 matrix.
!

      do lm1 = 1,9,1
         do lm2 = 1,9,1
            fullbo(lm1,lm2) = 0.0_dp
         enddo
      enddo

      jl1 = 0
      do il1 = 1,nl1,1
         l1 = llist1(il1)
         if (l1 == 0) then
            lm1 = 0
         elseif (l1 == 1) then
            lm1 = 1
         elseif (l1 == 2) then
            lm1 = 4
         endif
         do m1 = 1,2*l1+1,1
            jl1 = jl1 + 1
            lm1 = lm1 + 1
            jl2 = 0
            do il2 = 1,nl2,1
               l2 = llist2(il2)
               if (l2 == 0) then
                  lm2 = 0
               elseif (l2 == 1) then
                  lm2 = 1
               elseif (l2 == 2) then
                  lm2 = 4
               endif
               do m2 = 1,2*l2+1,1
                  jl2 = jl2 + 1
                  lm2 = lm2 + 1
                  fullbo(lm1,lm2) = bo(jl1,jl2)
               enddo
            enddo
         enddo
      enddo

!
!    Rotate the axes.
!

!     WRITE(6,'(9F7.3)') ((FULLBO(LM1,LM2),LM2=1,9,1),LM1=1,9,1)

      call rotbo(dr,fullbo,rbo)

!     WRITE(6,'(9F7.3)') ((RBO(LM1,LM2),LM2=1,9,1),LM1=1,9,1)

!
!    Extract the bond orders.
!

      bosss = rbo(1,1)

      bosps = rbo(1,4)
      bopss = rbo(4,1)

      bopps = rbo(4,4)
      boppp = 0.5_dp*(rbo(2,2)+rbo(3,3))

      bosds = rbo(1,5)
      bodss = rbo(5,1)

      bopds = rbo(4,5)
      bodps = rbo(5,4)
      bopdp = 0.5_dp*(rbo(2,8)+rbo(3,9))
      bodpp = 0.5_dp*(rbo(8,2)+rbo(9,3))

      bodds = rbo(5,5)
      boddp = 0.5_dp*(rbo(8,8)+rbo(9,9))
      boddd = 0.5_dp*(rbo(6,6)+rbo(7,7))

!
!    Evaluate the hopping integrals.
!

      magdr = sqrt(dr(1)**2 + dr(2)**2 + dr(3)**2)
      call radf(magdr,vsss,vsps,vpss,vpps,vppp, & 
     &          vsds,vdss,vpds,vdps,vpdp,vdpp, & 
     &          vdds,vddp,vddd,bnddat(1,btype(z1,z2)), & 
     &          bndscl(1,1,btype(z1,z2)))

!
!    Evaluate the bond energies.
!

      besss =       bosss*vsss

      besps =       bosps*vsps
      bepss =       bopss*vpss

      bepps =       bopps*vpps
      beppp = 2.0_dp*boppp*vppp

      besds =       bosds*vsds
      bedss =       bodss*vdss

      bepds =       bopds*vpds
      bedps =       bodps*vdps
      bepdp = 2.0_dp*bopdp*vpdp
      bedpp = 2.0_dp*bodpp*vdpp

      bedds =       bodds*vdds
      beddp = 2.0_dp*boddp*vddp
      beddd = 2.0_dp*boddd*vddd

      end









