 
      subroutine bonoav(ia)
          use mod_precision
          use mod_all_scalar
          use mod_const
          use mod_srt
          use ab_io
          use mod_ham
          use mod_chi
          use mod_funptr
!
!    This is a routine to write out the bondorder information.
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
!       include "Include/BondOrder.array"
!      include "Include/Force.array"
!       include "Include/Hamilt.array"
      include "Include/Moment.array"
      include "Include/NebList.array"
!      include "Include/NotSRT.array"
      include "Include/PosVel.array"
!      include "Include/RecurCoef.array"
!      include "Include/SRT.array"

!
!    Declare the simple variables.
!

      real(dp) :: bondeng,hopint
      real(dp) :: vsss,vsps,vpss,vpps,vppp
      real(dp) :: vsds,vdss,vpds,vdps,vpdp,vdpp
      real(dp) :: vdds,vddp,vddd
      real(dp) :: bosss,bosps,bopss,bopps,boppp
      real(dp) :: bosds,bodss,bopds,bodps,bopdp,bodpp
      real(dp) :: bodds,boddp,boddd
      real(dp) :: besss,besps,bepss,bepps,beppp
      real(dp) :: besds,bedss,bepds,bedps,bepdp,bedpp
      real(dp) :: bedds,beddp,beddd

      integer i,j,nmax,wrflag
      integer nla,nlb
      integer ia,ib,la,lb
      integer nstta,nsttb
      integer za,zb
      integer ja,ja0

      character*80 filename
      character*4 anum

!*** A displacement vector.
      real(dp) :: dr(3)
!*** Forces.
      real(dp) :: fab(3)
!*** Screening function
      real(dp) :: scf(14), scfcut(14),dsf(3,14)

!
!    Read in data.
!

      wrflag = 1
      do i = 1,ia,1
         call rdab(wrflag,i)
      enddo
      close(50)

!
!    Evaluate the susceptibilities.
!
      call assoc_ab(ia)
      call assoc_chi(ia)

      do la = 1,nchain,1
         nmax = lchain(la)
         if (term == 1) then
            call getchisrt(lef,nmax+1,la)
         elseif ((term == 2).or.(term == 3)) then
            call getchinull(lef,nmax,mrec,chia(0,la), & 
     &                      chib(0,la),diag(1,la), & 
     &                      eigvec(1,1,la),kt)
         endif
      enddo

!
!    Open the output files.
!

      write(anum,'(I4)') ia
      if (ia < 10) then
         anum(1:3) = '000'
      elseif (ia < 100) then
         anum(1:2) = '00'
      elseif (ia < 1000) then
         anum(1:1) = '0'
      endif

      filename = genfile(1:lengfn)//'_A'//anum//'.bo'
      open(unit = 10,file = filename,status = 'NEW')

      filename = genfile(1:lengfn)//'_A'//anum//'.rotbo'
      open(unit = 11,file = filename,status = 'NEW')

!
!    For each bond do ..
!

      write(10,'(''Atom # '',I3)') ia
      write(11,'(''Atom # '',I3)') ia

      za = z(ia)
      call states(za,nla,nstta,llista)

!
!    For each neighbor ....
!

      ja = aptr(ia)
      ja0 = ja
      do while (bptr(ja) /= eol)
         ib = bptr(ja)
         zb = z(ib)
         call states(zb,nlb,nsttb,llistb)
         if (ia /= ib) then

            write(10,'(''Neighboring atom # '',I3)') ib
            write(11,'(''Neighboring atom # '',I3)') ib

!
!          Evaluate the gradients of the hopping integrals.
!

            dr(1) = ad(1,ia) - ad(1,ib)
            dr(2) = ad(2,ia) - ad(2,ib)
            dr(3) = ad(3,ia) - ad(3,ib)


!
!          Evaluate the screening function and its derivative.
!

            if (scf_include == 1) then

               call screenf(ia,ib,scf)
               call scrcut(scf,scfcut)
               call dscreenf(ia,ib,scf,dsf)

            else
               do i=1,14,1
                  scfcut(i) = 1.0_dp
                  do j=1,3
                     dsf(j,i)=0.0_dp
                  enddo
               enddo
            endif


            call grdmat(grad,dr,za,nla,nstta,llista,zb,nlb,nsttb,llistb,mxnstat,ml,llista,llistb, & 
     &                  bnddat,bndscl,btype,nbtype,mbtype)

            if (momflg == 2) call utranv(transm(1,1,ia),grad, & 
     &                                   transm(1,1,map(ib)),trwork, & 
     &                                   mxnstat,nstta,nsttb)

!
!          Evaluate the bond orders.
!

            do la = 1,nstta,1
               nmax = lchain(la)
               do lb = 1,nsttb,1

                  bo(la,lb) = chia(0,la)* & 
     &                        darec(0,la,lb,ja-ja0+1)

                  do j = 1,nmax,1
                     bo(la,lb) = bo(la,lb) & 
     &                      + chia(j,la)* & 
     &                        darec(j,la,lb,ja-ja0+1) & 
     &                + 2.0_dp*chib(j,la)* & 
     &                        dbrec(j,la,lb,ja-ja0+1)
                  enddo

                  if (term == 1) bo(la,lb) = bo(la,lb) & 
     &                   + 2.0_dp*chib(nmax+1,la)* & 
     &                     dbrec(nmax+1,la,lb,ja-ja0+1)

                  bo(la,lb) = -2.0_dp*bo(la,lb)

                  write(10,'(''Bond Order ['',I2,'','',I2,''] = '', & 
     &                  G12.5,'' H = '',G12.5,'' DH = '',3G13.5)') & 
     &                  la,lb,bo(la,lb),darec(0,la,lb,ja-ja0+1), & 
     &                  grad(1,la,lb),grad(2,la,lb),grad(3,la,lb)

               enddo
            enddo

!
!          Evaluate the bond order etc for the rotated orbitals.
!

            call bndinf(za,zb,dr,mxnstat,ml,llista,llistb, & 
     &                  bnddat,bndscl,btype,nbtype,mbtype,bo, & 
     &                  vsss,vsps,vpss,vpps,vppp,vsds,vdss, & 
     &                  vpds,vdps,vpdp,vdpp,vdds,vddp,vddd, & 
     &                  bosss,bosps,bopss,bopps,boppp,bosds, & 
     &                  bodss,bopds,bodps,bopdp,bodpp,bodds, & 
     &                  boddp,boddd,besss,besps,bepss,bepps, & 
     &                  beppp,besds,bedss,bepds,bedps,bepdp, & 
     &                  bedpp,bedds,beddp,beddd)


!
!  Apply the screening function.
!

            vsss = vsss * scfcut(1)
            vsps = vsps * scfcut(2)
            vpss = vpss * scfcut(3)
            vpps = vpps * scfcut(4)
            vppp = vppp * scfcut(5)
            vsds = vsds * scfcut(6)
            vdss = vdss * scfcut(7)
            vpds = vpds * scfcut(8)
            vdps = vdps * scfcut(9)
            vpdp = vpdp * scfcut(10)
            vdpp = vdpp * scfcut(11)
            vdds = vdds * scfcut(12)
            vddp = vddp * scfcut(13)
            vddd = vddd * scfcut(14)
        

            write(11,'(''Overlap   Hopping Integral   Screening'', & 
     &                 '' Function   Bond Order    Bond Energy'')')
 1000       format(2x,a3,9x,f7.3,13x,f7.3,10x,f7.3,7x,f7.3)
            write(11,1000) 'SSS',vsss,scfcut(1),bosss,besss*scfcut(1)
            write(11,1000) 'SPS',vsps,scfcut(2),bosps,besps*scfcut(2)
            write(11,1000) 'PSS',vpss,scfcut(3),bopss,bepss*scfcut(3)
            write(11,1000) 'PPS',vpps,scfcut(4),bopps,bepps*scfcut(4)
            write(11,1000) 'PPP',vppp,scfcut(5),boppp,beppp*scfcut(5)
            write(11,1000) 'SDS',vsds,scfcut(6),bosds,besds*scfcut(6)
            write(11,1000) 'DSS',vdss,scfcut(7),bodss,bedss*scfcut(7)
            write(11,1000) 'PDS',vpds,scfcut(8),bopds,bepds*scfcut(8)
            write(11,1000) 'DPS',vdps,scfcut(9),bodps,bedps*scfcut(9)
            write(11,1000) 'PDP',vpdp,scfcut(10),bopdp,bepdp*scfcut(10)
            write(11,1000) 'DPP',vdpp,scfcut(11),bodpp,bedpp*scfcut(11)
            write(11,1000) 'DDS',vdds,scfcut(12),bodds,bedds*scfcut(12)
            write(11,1000) 'DDP',vddp,scfcut(13),boddp,beddp*scfcut(13)
            write(11,1000) 'DDD',vddd,scfcut(14),boddd,beddd*scfcut(14)


!
!          Calculate the bond energy.
!

            bondeng = 0.0_dp
            do la = 1,nstta,1
               do lb = 1,nsttb,1
                  hopint = darec(0,la,lb,ja-ja0+1)
                  bondeng = bondeng + hopint*bo(la,lb)
               enddo
            enddo

!
!          Calculate the bond forces.
!

            fab(1) = 0.0_dp
            fab(2) = 0.0_dp
            fab(3) = 0.0_dp
            do la = 1,nstta,1
               do lb = 1,nsttb,1
                  fab(1) = fab(1) - grad(1,la,lb)*bo(la,lb)
                  fab(2) = fab(2) - grad(2,la,lb)*bo(la,lb)
                  fab(3) = fab(3) - grad(3,la,lb)*bo(la,lb)
               enddo
            enddo

            write(10,'(''Totals:                          '', & 
     &            '' E = '',G12.5,''  F = '',3G13.5)') & 
     &            bondeng,fab

         endif

         ja = ja + 1

      enddo

      close(10)
      close(11)

      end

