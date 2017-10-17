 
      subroutine kbohuge(ia)
          use mod_precision


          use mod_all_scalar

          use mod_const

!
!    This is a rotuine to evaluate the band structure energy contribution
!     to the forces.
!

      implicit none

!
!    Include constants.
!

!      include "../Include/ALL.const"

!
!    Include scalars.
!

!      include "../Include/ALL.scalar"

!
!    Include arrays.
!

      include "../Include/Atom.array"
      include "../Include/BondOrder.array"
!      include "../Include/Force.array"
      include "../Include/KHamilt.array"
!      include "../Include/NebList.array"
      include "../Include/PosVel.array"

!
!    Declare the simple variables.
!

      real(dp) :: rcc,icc,rfac,ifac,dea
      real(dp) :: fac
      real(dp) :: phase
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

      integer nstta,nsttb,nla,nlb
      integer ik,in,iha,ihb
      integer ia,ib,la,lb
      integer za,zb

      character*80 filename
      character*4 anum

!
!    Declare the arrays.
!

!*** Atomic displacement.
      real(dp) :: dr(3)
!*** Forces.
      real(dp) :: fab(3)

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

      do ib = 1,nd,1
         zb = z(ib)
         call states(zb,nlb,nsttb,llistb)
         if (ia /= ib) then

            write(10,'(''Neighboring atom # '',I3)') ib
            write(11,'(''Neighboring atom # '',I3)') ib

!
!          Evaluate the hopping integrals and their gradients.
!

            dr(1) = ad(1,ia) - ad(1,ib)
            dr(2) = ad(2,ia) - ad(2,ib)
            dr(3) = ad(3,ia) - ad(3,ib)

            dea = de(ia)
            call matel(za,zb,dr,nstta,nsttb,ksubh,dea, & 
     &                 mxnstat,ml,llista,llistb)

!      Incorrect arglist
            call grdmat(kgrad,dr,za,zb,mxnstat,ml,llista,llistb)

!
!          Evaluate the bond orders.
!

            do la = 1,nstta,1
               do lb = 1,nsttb,1

                  bo(la,lb) = 0.0_dp

                  do ik = 1,nk,1

                     phase = kk(1,ik)*dr(1)+kk(2,ik)*dr(2) & 
     &                     + kk(3,ik)*dr(3)
                     rfac = cos(phase)
                     ifac = sin(phase)

                     in = 1
                     do while((in <= kmxh).and.(occ(in,ik) > 1.0e-6_dp))

                        fac = occ(in,ik)

                        iha = khpos(ia) + la
                        ihb = khpos(ib) + lb

                        rcc = (real(kpsi(iha,in,ik), dp)* & 
     &                         real(kpsi(ihb,in,ik), dp)+ & 
     &                         aimag(kpsi(iha,in,ik))* & 
     &                         aimag(kpsi(ihb,in,ik)))*fac
                        icc = (real(kpsi(iha,in,ik), dp)* & 
     &                         aimag(kpsi(ihb,in,ik))- & 
     &                         aimag(kpsi(iha,in,ik))* & 
     &                         real(kpsi(ihb,in,ik), dp))*fac

                        bo(la,lb) = bo(la,lb) + rcc*rfac - icc*ifac

                        in = in + 1

                     enddo

                  enddo

                  write(10,'(''Bond Order ['',I2,'','',I2,''] = '', & 
     &                  G12.5,'' H = '',G12.5,'' DH = '',3G13.5)') & 
     &                  la,lb,bo(la,lb),ksubh(la,lb), & 
     &                  kgrad(1,la,lb),kgrad(2,la,lb),kgrad(3,la,lb)

               enddo
            enddo

!
!          Evaluate the bond order etc for the rotated orbitals.
!

            call bndinf(za,zb,dr,mxnstat,ml,llista,llistb,bo, & 
     &                  vsss,vsps,vpss,vpps,vppp,vsds,vdss, & 
     &                  vpds,vdps,vpdp,vdpp,vdds,vddp,vddd, & 
     &                  bosss,bosps,bopss,bopps,boppp,bosds, & 
     &                  bodss,bopds,bodps,bopdp,bodpp,bodds, & 
     &                  boddp,boddd,besss,besps,bepss,bepps, & 
     &                  beppp,besds,bedss,bepds,bedps,bepdp, & 
     &                  bedpp,bedds,beddp,beddd)

            write(11,'(''Overlap   Hopping Integral   Bond Order'', & 
     &                 ''   Bond Energy'')')
 1000       format(2x,a3,9x,f7.3,9x,f7.3,6x,f7.3)
            write(11,1000) 'SSS',vsss,bosss,besss
            write(11,1000) 'SPS',vsps,bosps,besps
            write(11,1000) 'PSS',vpss,bopss,bepss
            write(11,1000) 'PPS',vpps,bopps,bepps
            write(11,1000) 'PPP',vppp,boppp,beppp
            write(11,1000) 'SDS',vsds,bosds,besds
            write(11,1000) 'DSS',vdss,bodss,bedss
            write(11,1000) 'PDS',vpds,bopds,bepds
            write(11,1000) 'DPS',vdps,bodps,bedps
            write(11,1000) 'PDP',vpdp,bopdp,bepdp
            write(11,1000) 'DPP',vdpp,bodpp,bedpp
            write(11,1000) 'DDS',vdds,bodds,bedds
            write(11,1000) 'DDP',vddp,boddp,beddp
            write(11,1000) 'DDD',vddd,boddd,beddd

!
!          Calculate the bond energy.
!

            bondeng = 0.0_dp
            do la = 1,nstta,1
               do lb = 1,nsttb,1
                  hopint = ksubh(la,lb)
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
                  fab(1) = fab(1) - kgrad(1,la,lb)*bo(la,lb)
                  fab(2) = fab(2) - kgrad(2,la,lb)*bo(la,lb)
                  fab(3) = fab(3) - kgrad(3,la,lb)*bo(la,lb)
               enddo
            enddo

            write(10,'(''Totals:                          '', & 
     &            '' E = '',G12.5,''  F = '',3G13.5)') & 
     &            bondeng,fab

         endif

      enddo

      close(10)
      close(11)

      end

