 
!***************************************************************************
!                                                                           *
!                                                                           *
!                               Matous Mrovec                               *
!              Department of Materials Science and Engineering              *
!                        University of Pennsylvania                         *
!                                                                           *
!                                                                           *
!****************************************************************************
!                                                                           *
!                                                                           *
! This subroutine creates plot file for scaling of screened bond integrals. *
!                                                                           *
!                                                                           *
!****************************************************************************


      subroutine plotscr( )
          use mod_precision


          use mod_all_scalar

          use mod_const
          use mod_atom_ar

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
!      include "Include/Hamilt.array"
!      include "Include/Moment.array"
!      include "Include/NebList.array"
      include "Include/NebListL.array"
      include "Include/PosVel.array"
      include "Include/BEC.array"


!
!    Common blocks introduced by A.Girshick.
!    
 
      include "Include/ag.conn"

!
!    Declare the simple variables.
!

      real(dp) :: rxab,ryab,rzab,rab,dalpha,vid
      real(dp) ::  alpha_lo, alpha_hi, alpha_step, alpha
      real(dp) ::  vsss,vsps,vpss,vpps,vppp
      real(dp) ::  vsds,vdss,vpds,vdps,vpdp,vdpp
      real(dp) ::  vdds,vddp,vddd
      real(dp) :: a1a2,a1a3,a2a3
      integer           i,j,k,it,bt
      integer           ia,ib
      integer           points
      real(dp) ::  rx, ry, rz, rpj, rh , mpw, getlamzero
      real(dp) ::  del_sum, pom_del, chi, e_env, venv
      real(dp) ::  dellam (2)
      integer  ja, la, ipick

!
!    Declare the arrays.
!
      real(dp) :: ta(3,3)
      real(dp) :: a1(3,3)
      real(dp) :: scf(14),scfcut(14)
      real(dp) :: v(10)

 10   write(*,*) "Screened bond integral scaling"
      write(6,'(''Enter first atom > '',$)')
      read(5,*) ia
      write(6,'(''Enter second atom > '',$)')
      read(5,*) ib

      if (ia == ib) goto 10


!  Define the mesh for scaling (usually 3 points).

! TMATRIX and number of points are defined in ucell.in
! in the same way as for elastic constants calculation

      alpha_hi = thi
      alpha_lo = tlo
      points = nt

      write(6,'("---"/)')
      write(6,'("*** Scaling of bond integrals ***")')
      write(6,'("Range from ",F6.3," to ",F6.3)') tlo,thi
      write(6,'("Number of points: ",I2/)') nt
!      WRITE(6,'(//)')

      if (points < 2) then
         alpha_hi = 0.0
         alpha_lo = 0.0
         alpha_step = 0.0
      else
         alpha_step = ( alpha_hi - alpha_lo ) / ( points - 1 )
      endif

      a1(1,1) = a(1,1)
      a1(1,2) = a(1,2)
      a1(1,3) = a(1,3)

      a1(2,1) = a(2,1)
      a1(2,2) = a(2,2)
      a1(2,3) = a(2,3)

      a1(3,1) = a(3,1)
      a1(3,2) = a(3,2)
      a1(3,3) = a(3,3)

      vid = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) & 
     &    + a(1,2)*(a(2,3)*a(3,1)-a(2,1)*a(3,3)) & 
     &    + a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))


      do i = 1,3,1
         do j = 1,3,1
            ta(i,j) = 0.0_dp
            do k = 1,3,1
               ta(i,j) = ta(i,j) + a1(i,k)*tmatrix(k,j)
            enddo
         enddo
      enddo

      dalpha = (thi - tlo)/real(nt-1, dp)
      talpha(1) = tlo

!      idebug=1
      open ( unit = 90, file = 'scrbond.dat')
      open ( unit = 91, file = 'env.dat')

      do it = 1,nt,1

!  Creating the strained block.

         do i = 1,3,1
            do j = 1,3,1
               a(j,i) = a1(j,i) + talpha(it)*ta(j,i)
            enddo
         enddo

         v(it) = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) & 
     &        + a(1,2)*(a(2,3)*a(3,1)-a(2,1)*a(3,3)) & 
     &        + a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))

         if (idebug == 1) write(9,'(''A = '',3F10.5)') a 
         do i = 1,nd,1
            ad(1,i) = a(1,1)*d(1,i)+a(2,1)*d(2,i)+a(3,1)*d(3,i)
            ad(2,i) = a(1,2)*d(1,i)+a(2,2)*d(2,i)+a(3,2)*d(3,i)
            ad(3,i) = a(1,3)*d(1,i)+a(2,3)*d(2,i)+a(3,3)*d(3,i)
            if (idebug == 1) write(9    ,'(''AD('',I4,'') = '', & 
     &           3G13.5)') i,ad(1,i),ad(2,i),ad(3,i)
         enddo
         do i = 1,ninert,1
            adinert(1,i) = a(1,1)*dinert(1,i)+a(2,1)*dinert(2,i)+ & 
     &                     a(3,1)*dinert(3,i)
            adinert(2,i) = a(1,2)*dinert(1,i)+a(2,2)*dinert(2,i)+ & 
     &                     a(3,2)*dinert(3,i)
            adinert(3,i) = a(1,3)*dinert(1,i)+a(2,3)*dinert(2,i)+ & 
     &                     a(3,3)*dinert(3,i)
         enddo

!         WRITE(*,*) "\\ *** --- *** \\"
!  Set up the neighbor lists for the bond part.

         a1a2 = abs(a(1,1)*a(2,1)+a(1,2)*a(2,2)+a(1,3)*a(2,3))
         a1a3 = abs(a(1,1)*a(3,1)+a(1,2)*a(3,2)+a(1,3)*a(3,3))
         a2a3 = abs(a(2,1)*a(3,1)+a(2,2)*a(3,2)+a(2,3)*a(3,3))
         if ((a1a2 < 1.0e-6_dp).and.(a1a3 < 1.0e-6_dp).and. & 
     &        (a2a3 < 1.0e-6_dp)) then
            nbtflg = 2
!            WRITE(*,*) "NBTFLG = 2  ->  unit cell is tetragonal"
         else
            nbtflg = 1
!            WRITE(*,*) "NBTFLG = 1  ->  unit cell is NOT tetragonal"
         endif

         nbtflg = 2
         call bldnebt()


!  Evaluate the screened bond integrals of the strained configuration.

         rxab = ad( 1, ib ) - ad( 1, ia )
         ryab = ad( 2, ib ) - ad( 2, ia )
         rzab = ad( 3, ib ) - ad( 3, ia )               
         rab = sqrt(rxab*rxab+ryab*ryab+rzab*rzab)


         if (scf_include == 1)then
            call screenf(ia,ib,scf)
            call scrcut(scf,scfcut)
         else
            do i=1,14,1
               scfcut(i) = 1.0_dp
            enddo
         endif

         write(74,'(4F12.8)') rab ,scfcut(12),scfcut(13),scfcut(14)

         bt=btype(z(ia),z(ib))

         call radf(rab,vsss,vsps,vpss,vpps,vppp, & 
     &             vsds,vdss,vpds,vdps,vpdp,vdpp, & 
     &             vdds,vddp,vddd,bnddat(1,bt), & 
     &             bndscl(1,1,bt))


!     Scaling for the repulsive environmental term

         call bldnebtfs()

         ja = aptrl( ia )

         do while ( bptrl( ja )  /=  eol )
            j = bptrl( ja )
            if ( j  /=  ia ) then
               ipick = ind( j )
               rx = ad( 1, ia ) - ad( 1, j )
               ry = ad( 2, ia ) - ad( 2, j )
               rz = ad( 3, ia ) - ad( 3, j )
               rpj = sqrt( rx*rx + ry*ry + rz*rz )
               
               del_sum = del_sum + chi( rpj, ipick)

            endif
            ja = ja + 1
         enddo

         dellam( 1 ) = exp( (1.0_dp / mpw(ipick)) * log( del_sum ) )

         ja = aptrl( ib )

         do while ( bptrl( ja )  /=  eol )
            j = bptrl( ja )
            if ( j  /=  ib ) then
               ipick = ind( j )
               rx = ad( 1, ib ) - ad( 1, j )
               ry = ad( 2, ib ) - ad( 2, j )
               rz = ad( 3, ib ) - ad( 3, j )
               rpj = sqrt( rx*rx + ry*ry + rz*rz )
               
               del_sum = del_sum + chi( rpj, ipick)

            endif
            ja = ja + 1
         enddo

         dellam( 2 ) = exp( (1.0_dp / mpw(ipick)) * log( del_sum ) )

         pom_del = 0.5_dp * ( dellam(1) + dellam(2) )
         e_env = venv(rab,ipick,pom_del)

         write(91,'(2F13.8)') rab, e_env

         write(*,*) dellam( ia ),dellam( ib )
         write(*,'(3F13.8)') rab, pom_del, e_env

         if (v(it)/vid > 0.99.and.v(it)/vid < 1.01) then

         write(6,'(" *** Results for volume ratio = 1 ***")')
         write(6,'("Atom A ",3F12.8)') ad(1,ia),ad(2,ia),ad(3,ia)
         write(6,'("Atom B ",3F12.8)') ad(1,ib),ad(2,ib),ad(3,ib)
         write(6,'("Interatomic distance = ",F12.6)') rab

         if (vsss /= 0.0) then
         write(6,'("SSS: SCRF =",F7.3," VSSS =",F7.3," -> SSS =",F7.3)')  &
     &        scfcut(1), vsss, vsss*scfcut(1)
         endif

         if (vsps /= 0.0) then
         write(6,'("SPS: SCRF =",F7.3," VSPS =",F7.3," -> SPS =",F7.3)')  &
     &        scfcut(2), vsps, vsps*scfcut(2)
         endif

         if (vpss /= 0.0) then
         write(6,'("PSS: SCRF =",F7.3," VPSS =",F7.3," -> SPS =",F7.3)')  &
     &        scfcut(3), vpss, vpss*scfcut(3)
         endif

         if (vpps /= 0.0) then
         write(6,'("PPS: SCRF =",F7.3," VPPS =",F7.3," -> PPS =",F7.3)')  &
     &        scfcut(4), vpps, vpps*scfcut(4)
         endif

         if (vppp /= 0.0) then
         write(6,'("PPP: SCRF =",F7.3," VPPP =",F7.3," -> PPP =",F7.3)')  &
     &        scfcut(5), vppp, vppp*scfcut(5)
         endif

         if (vsds /= 0.0) then
         write(6,'("SDS: SCRF =",F7.3," VSDS =",F7.3," -> SDS =",F7.3)')  &
     &        scfcut(6), vsds, vsds*scfcut(6)
         endif

         if (vdss /= 0.0) then
         write(6,'("DSS: SCRF =",F7.3," VDSS =",F7.3," -> DSS =",F7.3)')  &
     &        scfcut(7), vdss, vdss*scfcut(7)
         endif

         if (vpds /= 0.0) then
         write(6,'("PDS: SCRF =",F7.3," VPDS =",F7.3," -> PDS =",F7.3)') &
     &        scfcut(8),vpds, vpds*scfcut(8)
         endif

         if (vdps /= 0.0) then
         write(6,'("DPS: SCRF =",F7.3," VDPS =",F7.3," -> DPS =",F7.3)') &
     &        scfcut(9),vdps,vdps*scfcut(9)
         endif

         if (vpdp /= 0.0) then
         write(6,'("PDP: SCRF =",F7.3," VPDP =",F7.3," -> PDS =",F7.3)') &
     &        scfcut(10),vpdp,vpdp *scfcut(10)
         endif

         if (vdpp /= 0.0) then
         write(6,'("DPP: SCRF =",F7.3," VDPP =",F7.3," -> DPP =",F7.3)') &
     &        scfcut(11),vdpp,vdpp*scfcut(11)
         endif

         if (vdds /= 0.0) then
         write(6,'("DDS: SCRF =",F7.3," VDDS =",F7.3," -> DDS =",F7.3)') &
     &        scfcut(12),vdds,vdds*scfcut(12)
         endif

         if (vddp /= 0.0) then
         write(6,'("DDP: SCRF =",F7.3," VDDP =",F7.3," -> DDP =",F7.3)') &
     &        scfcut(13),vddp,vddp*scfcut(13)
         endif

         if (vddd /= 0.0) then
         write(6,'("DDD: SCRF =",F7.3," VDDD =",F7.3," -> DDD =",F7.3)') &
     &        scfcut(14),vddd,vddd*scfcut(14)
         endif

         write(6,'("--------------------")')

         endif

! Apply the screening function

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


!         VSSS = VSSS * ( 1 - SCF(1) )
!         VSPS = VSPS * ( 1 - SCF(2) )
!         VPSS = VPSS * ( 1 - SCF(3) )
!         VPPS = VPPS * ( 1 - SCF(4) )
!         VPPP = VPPP * ( 1 - SCF(5) )
!         VSDS = VSDS * ( 1 - SCF(6) )
!         VDSS = VDSS * ( 1 - SCF(7) )
!         VPDS = VPDS * ( 1 - SCF(8) )
!         VDPS = VDPS * ( 1 - SCF(9) )
!         VPDP = VPDP * ( 1 - SCF(10) )
!         VDPP = VDPP * ( 1 - SCF(11) )
!         VDDS = VDDS * ( 1 - SCF(12) )
!         VDDP = VDDP * ( 1 - SCF(13) )
!         VDDD = VDDD * ( 1 - SCF(14) )


!         WRITE(90,'(15F13.8)') RAB,VSSS,VSPS,VPSS,VPPS,VPPP,
!     &         VSDS,VDSS,VPDS,VDPS,VPDP,VDPP,VDDS,VDDP,VDDD

         write(90,'(5F13.8)') rab, vpps, vppp
!         WRITE(90,'(8F13.8)') RAB,VDDS,VDDP,VDDD
!         WRITE(90,'(3F13.8)') TALPHA(IT),VPPS,VPPP

         talpha(it+1) = talpha(it) + dalpha

      enddo

      close( 91 )
      close( 90 )

      a(1,1) = a1(1,1)
      a(1,2) = a1(1,2)
      a(1,3) = a1(1,3)

      a(2,1) = a1(2,1)
      a(2,2) = a1(2,2)
      a(2,3) = a1(2,3)

      a(3,1) = a1(3,1)
      a(3,2) = a1(3,2)
      a(3,3) = a1(3,3)


      return
      end

