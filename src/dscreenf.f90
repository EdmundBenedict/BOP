 
      subroutine dscreenf(ia,ib,scf,dsf)
          use mod_precision


          use mod_all_scalar

          use mod_const
      use mod_atom_ar, only : btype
!
!    This is a subroutine which calculates the derivative of
!    the screening function for a bond between atoms IA-IB.
!    Unfortunately it is quite a mess :-((.
!    (c) Matous Mrovec, Oxford/Penn, July 2000
!
!    Additionally messed up by Dimitar Pashov


      implicit none


      include "Include/Atom.array"
      include "Include/NebList.array"
      include "Include/PosVel.array"

!
!    Declare the simple variables.
!

      real(dp) :: rxab,ryab,rzab,rab,rik,rjk
      real(dp) :: rxik,ryik,rzik,rxjk,ryjk,rzjk
      real(dp) :: num,den,nod,theta,thetaj,c0,ca,cb
      real(dp) :: osss,opss,odss
      real(dp) :: gss,gps,gpp,gds,gdp,gdd
      real(dp) :: gssj,gpsj,gppj,gdsj,gdpj,gddj
      real(dp) :: vscale, dscale
      real(dp) :: hsss,hsps,hpss,hpps,hppp
      real(dp) :: hsds,hdss,hpds,hdps,hpdp
      real(dp) :: hdpp,hdds,hddp,hddd
      real(dp) :: hsssjk,hspsjk,hpssjk,hppsjk,hpppjk
      real(dp) :: hsdsjk,hdssjk,hpdsjk,hdpsjk,hpdpjk
      real(dp) :: hdppjk,hddsjk,hddpjk,hdddjk
      real(dp) :: hsssab,hspsab,hpssab,hppsab,hpppab
      real(dp) :: hsdsab,hdssab,hpdsab,hdpsab,hpdpab
      real(dp) :: hdppab,hddsab,hddpab,hdddab
      real(dp) :: dhsss,dhsps,dhpss,dhpps,dhppp
      real(dp) :: dhsds,dhdss,dhpds,dhdps,dhpdp
      real(dp) :: dhdpp,dhdds,dhddp,dhddd
      real(dp) :: dhsssab,dhspsab,dhpssab,dhppsab,dhpppab
      real(dp) :: dhsdsab,dhdssab,dhpdsab,dhdpsab,dhpdpab
      real(dp) :: dhdppab,dhddsab,dhddpab,dhdddab
      real(dp) :: dosss,dopss,dodss

      real(dp) :: vsss,vsps,vpss,vpps,vppp
      real(dp) :: vsds,vdss,vpds,vdps,vpdp,vdpp
      real(dp) :: vdds,vddp,vddd
      real(dp) :: dsss,dsps,dpss,dpps,dppp
      real(dp) :: dsds,ddss,dpds,ddps,dpdp,ddpp
      real(dp) :: ddds,dddp,dddd
      real(dp) :: dr(3)
      real(dp) :: dvsss(3),dvsps(3),dvpss(3),dvpps(3),dvppp(3)
      real(dp) :: dvsds(3),dvdss(3),dvpds(3),dvdps(3),dvpdp(3)
      real(dp) :: dvdpp(3),dvdds(3),dvddp(3),dvddd(3)


      real(dp) :: cos_i,sin_i,cos2_i,sin2_i
      real(dp) :: cos_j,sin_j,cos2_j,sin2_j
      real(dp) :: vecpr_i,vecpr_j

      real(dp) :: pom1,pom2,cden,cden2,ff,dff

      integer ia,ib
      integer ja,j,i,k
      integer bt,btj


!
!    Declare the local arrays.
!

      real(dp) :: scf(14),dsf(3,14)
      real(dp) :: mu2a(14),mu2b(14),oab(14),ojk(14)
      real(dp) :: mu2(14),mu3(14),c1(14)
      real(dp) :: doab(14),oab2(14),dtau(14)

      real(dp) :: dcos_i(3),dcos_j(3)
      real(dp) :: dsin_i(3),dsin_j(3)
      real(dp) :: dgss(3),dgps(3),dgpp(3)
      real(dp) :: dgds(3),dgdp(3),dgdd(3)
      real(dp) :: dgssj(3),dgpsj(3),dgppj(3)
      real(dp) :: dgdsj(3),dgdpj(3),dgddj(3)
      real(dp) :: dc1(14,3),dmu2(14,3),pdoab(14,3)
      real(dp) :: dmu2a(14,3),dmu2b(14,3)
      real(dp) :: dc1a(14,3),dc1bi(14,3),dc1bj(14,3)
!
!      EXTERNAL SCALE
!
! Initialization of variables


      do i=1,14
         mu2a(i) = 0.0_dp
         mu2b(i) = 0.0_dp
         do j=1,3
            dc1a(i,j) = 0.0_dp
            dc1bi(i,j) = 0.0_dp
            dc1bj(i,j) = 0.0_dp
            dmu2a(i,j) = 0.0_dp
            dmu2b(i,j) = 0.0_dp
            dmu2(i,j) = 0.0_dp
         enddo
      enddo

! Pair of atoms I-J

      rxab = ad( 1, ib ) - ad( 1, ia )
      ryab = ad( 2, ib ) - ad( 2, ia )
      rzab = ad( 3, ib ) - ad( 3, ia )               
      rab = sqrt(rxab*rxab + ryab*ryab + rzab*rzab)

      bt = btype(z(ia),z(ib))

      do i=1,14,1
         dtau(i) = bndscl(14,i,bt)
         if (bndscl(1,i,bt) == 0.0) then
            oab(i) = 0.0_dp 
            doab(i) = 0.0_dp 
         else
            oab(i) = vscale(rab,bndscl(1,i,bt),bndscl(2,i,bt), & 
     &                      bndscl(5,i,bt),bndscl(6,i,bt), & 
     &                      bndscl(3,i,bt),bndscl(4,i,bt), & 
     &                      bndscl(7,i,bt),bndscl(8,i,bt), & 
     &                      bndscl(9,i,bt),bndscl(10,i,bt), & 
     &                      bndscl(11,i,bt),bndscl(12,i,bt))

            oab(i) = oab(i)*bndscl(13,i,bt)

            doab(i) = dscale(rab,bndscl(1,i,bt),bndscl(2,i,bt), & 
     &                       bndscl(5,i,bt),bndscl(6,i,bt), & 
     &                       bndscl(3,i,bt),bndscl(4,i,bt), & 
     &                       bndscl(7,i,bt),bndscl(8,i,bt), & 
     &                       bndscl(9,i,bt),bndscl(10,i,bt), & 
     &                       bndscl(11,i,bt),bndscl(12,i,bt))

            doab(i) = doab(i)*bndscl(13,i,bt)

         endif
      enddo

      call radf(rab,hsssab,hspsab,hpssab,hppsab, & 
     &          hpppab,hsdsab,hdssab,hpdsab,hdpsab, & 
     &          hpdpab,hdppab,hddsab,hddpab,hdddab, & 
     &          bnddat(1,bt),bndscl(1,1,bt))


      call dradf(rab,dhsssab,dhspsab,dhpssab,dhppsab, & 
     &           dhpppab,dhsdsab,dhdssab,dhpdsab,dhdpsab, & 
     &           dhpdpab,dhdppab,dhddsab,dhddpab,dhdddab, & 
     &           bnddat(1,bt),bndscl(1,1,bt))


!
!    For each ion (A) do ...
!

      ja = aptr( ia )

      do while ( bptr( ja )  /=  eol )
         k = bptr( ja )
         if ((k /= ia).and.(k /= ib)) then

! Pair of atoms I-K
! Hop from i to k

            rxik = ad( 1, k ) - ad( 1, ia )
            ryik = ad( 2, k ) - ad( 2, ia )
            rzik = ad( 3, k ) - ad( 3, ia )               
            rik = sqrt(rxik*rxik + ryik*ryik + rzik*rzik)

            num = rxik*rxab + ryik*ryab + rzik*rzab
            den = rab*rik
            nod = num/den

            if (nod > 1.0_dp-1.0e-6_dp) then
               theta = 0.0_dp
               do i=1,3
                  dsin_i(i)=0.0
               enddo
            elseif (nod < -1.0_dp+1.0e-6_dp) then
               theta = pi
               do i=1,3
                  dsin_i(i)=0.0
               enddo
            else
               theta = acos(nod)
            endif

            cos_i = cos(theta)
            sin_i = sin(theta)
            cos2_i = cos(2.0_dp*theta)
            sin2_i = sin(2.0_dp*theta)

            gss = 1.0_dp
            gps = cos_i
            gpp = sin_i
            gds = 0.25_dp*(1.0_dp+3.0_dp*cos2_i)
            gdp = sqrt3/2.0_dp*sin2_i
            gdd = sqrt3/4.0_dp*(1.0_dp-cos2_i)


! Derivatives of Cos(Theta) 

            dcos_i(1) = ( rxab/(rab*rab) + rxik/(rik*rik) ) * cos_i - &
     &                  ( rxab + rxik ) / (rab*rik)
            dcos_i(2) = ( ryab/(rab*rab) + ryik/(rik*rik) ) * cos_i - &
     &                  ( ryab + ryik ) / (rab*rik)
            dcos_i(3) = ( rzab/(rab*rab) + rzik/(rik*rik) ) * cos_i - &
     &                  ( rzab + rzik ) / (rab*rik)

! Derivatives of Sin(Theta)

            vecpr_i = sqrt((ryab*rzik-ryik*rzab)*(ryab*rzik-ryik*rzab)  & 
     &           + (rzab*rxik-rzik*rxab)* (rzab*rxik-rzik*rxab)  & 
     &           + (rxab*ryik-rxik*ryab)*(rxab*ryik-rxik*ryab)) 

!            VECPR_I = DSQRT( (RYAB*RZIK-RYIK*RZAB)**2 +
!     &                      (RZAB*RXIK-RZIK*RXAB)**2 +
!     &                      (RXAB*RYIK-RXIK*RYAB)**2 )

          if ((nod > -1.0_dp+1.0e-6_dp).and.(nod < 1.0_dp-1.0e-6_dp)) then
               dsin_i(1) = -cos_i/sin_i*dcos_i(1)
               dsin_i(2) = -cos_i/sin_i*dcos_i(2)
               dsin_i(3) = -cos_i/sin_i*dcos_i(3)
          endif

! Derivatives of Gll'I with respect to x(=1), y(=2), z(=3)
!           DGSS(1-3) = 0.D0

            dgps(1) = dcos_i(1)
            dgps(2) = dcos_i(2)
            dgps(3) = dcos_i(3)

            dgpp(1) = dsin_i(1)
            dgpp(2) = dsin_i(2)
            dgpp(3) = dsin_i(3)

            dgds(1) = 3.0_dp*cos_i * dcos_i(1)
            dgds(2) = 3.0_dp*cos_i * dcos_i(2)
            dgds(3) = 3.0_dp*cos_i * dcos_i(3)

            dgdp(1) = sqrt3*(dsin_i(1)*cos_i+sin_i*dcos_i(1))
            dgdp(2) = sqrt3*(dsin_i(2)*cos_i+sin_i*dcos_i(2))
            dgdp(3) = sqrt3*(dsin_i(3)*cos_i+sin_i*dcos_i(3))

            dgdd(1) = -sqrt3*cos_i * dcos_i(1)
            dgdd(2) = -sqrt3*cos_i * dcos_i(2)
            dgdd(3) = -sqrt3*cos_i * dcos_i(3)


            bt = btype(z(ia),z(k))

            call radf(rik,hsss,hsps,hpss,hpps, & 
     &                hppp,hsds,hdss,hpds,hdps, & 
     &                hpdp,hdpp,hdds,hddp,hddd, & 
     &                bnddat(1,bt),bndscl(1,1,bt))


            call dradf(rik,dhsss,dhsps,dhpss,dhpps, & 
     &                 dhppp,dhsds,dhdss,dhpds,dhdps, & 
     &                 dhpdp,dhdpp,dhdds,dhddp,dhddd, & 
     &                 bnddat(1,bt),bndscl(1,1,bt))


! Pair of atoms J-K
! Hop from j to k

            rxjk = ad( 1, k ) - ad( 1, ib )
            ryjk = ad( 2, k ) - ad( 2, ib )
            rzjk = ad( 3, k ) - ad( 3, ib )               
            rjk = sqrt(rxjk*rxjk + ryjk*ryjk + rzjk*rzjk)

            num = rxjk*rxab + ryjk*ryab + rzjk*rzab
            den = rab*rjk
            nod = -num/den

            if (nod > 1.0_dp-1.0e-6_dp) then
               thetaj = 0.0_dp
               do i=1,3
                  dsin_j(i)=0.0
               enddo
            elseif (nod < -1.0_dp+1.0e-6_dp) then
               thetaj = pi
               do i=1,3
                  dsin_j(i)=0.0
               enddo
            else
               thetaj = acos(nod)
            endif

            cos_j = cos(thetaj)
            sin_j = sin(thetaj)
            cos2_j = cos(2.0_dp*thetaj)
            sin2_j = sin(2.0_dp*thetaj)

            gssj = 1.0_dp
            gpsj = cos_j
            gppj = sin_j
            gdsj = 0.25_dp*(1.0_dp+3.0_dp*cos2_j)
            gdpj = sqrt3/2.0_dp*sin2_j
            gddj = sqrt3/4.0_dp*(1.0_dp-cos2_j)

            vecpr_j = sqrt((ryab*rzjk-ryjk*rzab)*(ryab*rzjk-ryjk*rzab) & 
     &           + (rzab*rxjk-rzjk*rxab)*(rzab*rxjk-rzjk*rxab) & 
     &           + (rxab*ryjk-rxjk*ryab)*(rxab*ryjk-rxjk*ryab))

!            VECPR_J = SQRT( (RYAB*RZJK-RYJK*RZAB)**2 +
!     &                      (RZAB*RXJK-RZJK*RXAB)**2 +
!     &                      (RXAB*RYJK-RXJK*RYAB)**2 )

! Derivatives of Cos(ThetaJ) with respect to x(=1), y(=2), z(=3)

            dcos_j(1) = 1.0_dp/rab * (rxjk/rjk+rxab/rab*cos_j)
            dcos_j(2) = 1.0_dp/rab * (ryjk/rjk+ryab/rab*cos_j)
            dcos_j(3) = 1.0_dp/rab * (rzjk/rjk+rzab/rab*cos_j)

! Derivatives of Sin(ThetaJ)

          if ((nod > -1.0_dp+1.0e-6_dp).and.(nod < 1.0_dp-1.0e-6_dp)) then
               dsin_j(1) = -cos_j/sin_j*dcos_j(1)
               dsin_j(2) = -cos_j/sin_j*dcos_j(2)
               dsin_j(3) = -cos_j/sin_j*dcos_j(3)
          endif

! Derivatives of Gll'J with respect to x(=1), y(=2), z(=3)
!           DGSSJ(1-3) = 0.D0

            dgpsj(1) = dcos_j(1)
            dgpsj(2) = dcos_j(2)
            dgpsj(3) = dcos_j(3)

            dgppj(1) = dsin_j(1)
            dgppj(2) = dsin_j(2)
            dgppj(3) = dsin_j(3)

            dgdsj(1) = 3.0_dp*cos_j * dcos_j(1)
            dgdsj(2) = 3.0_dp*cos_j * dcos_j(2)
            dgdsj(3) = 3.0_dp*cos_j * dcos_j(3)

            dgdpj(1) = -sqrt3*cos2_j/sin_j*dcos_j(1)
            dgdpj(2) = -sqrt3*cos2_j/sin_j*dcos_j(2)
            dgdpj(3) = -sqrt3*cos2_j/sin_j*dcos_j(3)

            dgddj(1) = -sqrt3*cos_j * dcos_j(1)
            dgddj(2) = -sqrt3*cos_j * dcos_j(2)
            dgddj(3) = -sqrt3*cos_j * dcos_j(3)


            btj = btype(z(ib),z(k))

            do i=1,14,1
               if ((bndscl(1,i,btj) == 0.0).or.(rjk > rcut)) then
                  ojk(i) = 0.0_dp      
               else
                  ojk(i) = vscale(rjk,bndscl(1,i,btj),bndscl(2,i,btj), & 
     &                      bndscl(5,i,btj),bndscl(6,i,btj), & 
     &                      bndscl(3,i,btj),bndscl(4,i,btj), & 
     &                      bndscl(7,i,btj),bndscl(8,i,btj), & 
     &                      bndscl(9,i,btj),bndscl(10,i,btj), & 
     &                      bndscl(11,i,btj),bndscl(12,i,btj))

                  ojk(i) = ojk(i)*bndscl(13,i,btj)

               endif
            enddo
             
            call radf(rjk,hsssjk,hspsjk,hpssjk,hppsjk, & 
     &                hpppjk,hsdsjk,hdssjk,hpdsjk,hdpsjk, & 
     &                hpdpjk,hdppjk,hddsjk,hddpjk,hdddjk, & 
     &                bnddat(1,btj),bndscl(1,1,btj))


! Check whether we have s electrons

            if (bndscl(1,1,bt) == 0.0) then
               osss = 0.0_dp
               dosss = 0.0_dp
            else
               osss = vscale(rik,bndscl(1,1,bt),bndscl(2,1,bt), & 
     &                      bndscl(5,1,bt),bndscl(6,1,bt), & 
     &                      bndscl(3,1,bt),bndscl(4,1,bt), & 
     &                      bndscl(7,1,bt),bndscl(8,1,bt), & 
     &                      bndscl(9,1,bt),bndscl(10,1,bt), & 
     &                      bndscl(11,1,bt),bndscl(12,1,bt))

               osss = osss*bndscl(13,1,bt)


               dosss = dscale(rik,bndscl(1,1,bt),bndscl(2,1,bt), & 
     &                       bndscl(5,1,bt),bndscl(6,1,bt), & 
     &                       bndscl(3,1,bt),bndscl(4,1,bt), & 
     &                       bndscl(7,1,bt),bndscl(8,1,bt), & 
     &                       bndscl(9,1,bt),bndscl(10,1,bt), & 
     &                       bndscl(11,1,bt),bndscl(12,1,bt))

               dosss = dosss*bndscl(13,1,bt)

! Calculate second moment contributions

               mu2a(1) = mu2a(1) + osss*osss*gss*gss
               mu2a(2) = mu2a(1)
               mu2a(6) = mu2a(1)

! Calculate derivatives of MU2A(I)
               
! * SSS        MU2A(1) = MU2A(1) + OSSS*OSSS*GSS*GSS

               dmu2a(1,1) = dmu2a(1,1) + &
     &                      2.0_dp*osss*dosss*(-rxik/rik)
               dmu2a(1,2) = dmu2a(1,2) + &
     &                      2.0_dp*osss*dosss*(-ryik/rik)
               dmu2a(1,3) = dmu2a(1,3) + &
     &                      2.0_dp*osss*dosss*(-rzik/rik)


! * SPS        MU2A(2) = MU2A(2) + OSSS*OSSS*GSS*GSS

               dmu2a(2,1) = dmu2a(1,1)
               dmu2a(2,2) = dmu2a(1,2)
               dmu2a(2,3) = dmu2a(1,3)


! * SDS        MU2A(6) = MU2A(6) + OSSS*OSSS*GSS*GSS

               dmu2a(6,1) = dmu2a(1,1)
               dmu2a(6,2) = dmu2a(1,2)
               dmu2a(6,3) = dmu2a(1,3)


! Calculate derivatives of C1BI(I)

! *** SSS
               if (hsssab /= 0.0) then
                  cb = 0.5_dp*dtau(1)/hsssab* &
     &              (-hsss)*osss*oab(1)

                  dc1bi(1,1) = dc1bi(1,1) + 1.0_dp/hsssab * &
     &              ( dhsssab * cb * rxab/rab + 0.5_dp * dtau(1) * &
     &               ( -(dhsss*osss*oab(1)*(-rxik/rik)+ &
     &                   hsss*dosss*oab(1)*(-rxik/rik)+ &
     &                   hsss*osss*doab(1)*(-rxab/rab)) ))

                  dc1bi(1,2) = dc1bi(1,2) + 1.0_dp/hsssab * &
     &              ( dhsssab * cb * ryab/rab + 0.5_dp * dtau(1) * &
     &                ( -(dhsss*osss*oab(1)*(-ryik/rik)+ &
     &                   hsss*dosss*oab(1)*(-ryik/rik)+ &
     &                   hsss*osss*doab(1)*(-ryab/rab)) ))

                  dc1bi(1,3) = dc1bi(1,3) + 1.0_dp/hsssab * &
     &              ( dhsssab * cb * rzab/rab + 0.5_dp * dtau(1) * &
     &                ( -(dhsss*osss*oab(1)*(-rzik/rik)+ &
     &                   hsss*dosss*oab(1)*(-rzik/rik)+ &
     &                   hsss*osss*doab(1)*(-rzab/rab)) ))

               endif


! *** SPS
               if (hspsab /= 0.0) then
                  cb = 0.5_dp*dtau(2)/hspsab* &
     &              (-hsss)*osss*oab(2)

                  dc1bi(2,1) = dc1bi(2,1) + 1.0_dp/hspsab * &
     &              ( dhspsab * cb * rxab/rab + 0.5_dp * dtau(2) * &
     &                ( -(dhsss*osss*oab(2)*(-rxik/rik)+ &
     &                   hsss*dosss*oab(2)*(-rxik/rik)+ &
     &                   hsss*osss*doab(2)*(-rxab/rab)) ))

                  dc1bi(2,2) = dc1bi(2,2) + 1.0_dp/hspsab * &
     &              ( dhspsab * cb * ryab/rab + 0.5_dp * dtau(2) * &
     &                ( -(dhsss*osss*oab(2)*(-ryik/rik)+ &
     &                   hsss*dosss*oab(2)*(-ryik/rik)+ &
     &                   hsss*osss*doab(2)*(-ryab/rab)) ))

                  dc1bi(2,3) = dc1bi(2,3) + 1.0_dp/hspsab * &
     &              ( dhspsab * cb * rzab/rab + 0.5_dp * dtau(2) * &
     &                ( -(dhsss*osss*oab(2)*(-rzik/rik)+ &
     &                   hsss*dosss*oab(2)*(-rzik/rik)+ &
     &                   hsss*osss*doab(2)*(-rzab/rab)) ))

               endif

! *** SDS
               if (hsdsab /= 0.0) then
                  cb = 0.5_dp*dtau(6)/hsdsab* &
     &              (-hsss)*osss*oab(6)

                  dc1bi(6,1) = dc1bi(6,1) + 1.0_dp/hsdsab * &
     &              ( dhsdsab * cb * rxab/rab + 0.5_dp * dtau(6) * &
     &                ( -(dhsss*osss*oab(6)*(-rxik/rik)+ &
     &                   hsss*dosss*oab(6)*(-rxik/rik)+ &
     &                   hsss*osss*doab(6)*(-rxab/rab)) ))

                  dc1bi(6,2) = dc1bi(6,2) + 1.0_dp/hsdsab * &
     &              ( dhsdsab * cb * ryab/rab + 0.5_dp * dtau(6) * &
     &                ( -(dhsss*osss*oab(6)*(-ryik/rik)+ &
     &                   hsss*dosss*oab(6)*(-ryik/rik)+ &
     &                   hsss*osss*doab(6)*(-ryab/rab)) ))

                  dc1bi(6,3) = dc1bi(6,3) + 1.0_dp/hsdsab * &
     &              ( dhsdsab * cb * rzab/rab + 0.5_dp * dtau(6) * &
     &                ( -(dhsss*osss*oab(6)*(-rzik/rik)+ &
     &                   hsss*dosss*oab(6)*(-rzik/rik)+ &
     &                   hsss*osss*doab(6)*(-rzab/rab)) ))

               endif

!     ___  End of dC1BI(I)  ___


! Calculate derivatives of third moment contributions - C1A

               if (rjk < rcut) then

! *** SSS

               if (hsssab /= 0.0) then
                  ca = 0.5_dp*dtau(1)/hsssab* &
     &              (hsss*ojk(1)+osss*hsssjk)

!     &              ((HSSS*OJK(1)+OSSS*HSSSJK)
!     &              -HSSS*OSSS*OAB(1)
!     &              -OAB(1)*OJK(1)*HSSSJK)

                  dc1a(1,1) = dc1a(1,1) + 1.0_dp/hsssab * &
     &              ( dhsssab * ca * rxab/rab + 0.5_dp * dtau(1) * &
     &                 (dhsss*ojk(1)+dosss*hsssjk)*(-rxik/rik) )

!     &                ( (DHSSS*OJK(1)+DOSSS*HSSSJK)*(-RXIK/RIK)
!     &                 -(DHSSS*OSSS*OAB(1)*(-RXIK/RIK)+
!     &                   HSSS*DOSSS*OAB(1)*(-RXIK/RIK)+
!     &                   HSSS*OSSS*DOAB(1)*(-RXAB/RAB))
!     &                 -(DOAB(1)*(-RXAB/RAB)*OJK(1)*HSSSJK) ))

                  dc1a(1,2) = dc1a(1,2) + 1.0_dp/hsssab * &
     &              ( dhsssab * ca * ryab/rab + 0.5_dp * dtau(1) * &
     &                 (dhsss*ojk(1)+dosss*hsssjk)*(-ryik/rik) )

!     &                ( (DHSSS*OJK(1)+DOSSS*HSSSJK)*(-RYIK/RIK)
!     &                 -(DHSSS*OSSS*OAB(1)*(-RYIK/RIK)+
!     &                   HSSS*DOSSS*OAB(1)*(-RYIK/RIK)+
!     &                   HSSS*OSSS*DOAB(1)*(-RYAB/RAB))
!     &                 -(DOAB(1)*(-RYAB/RAB)*OJK(1)*HSSSJK) ))

                  dc1a(1,3) = dc1a(1,3) + 1.0_dp/hsssab * &
     &              ( dhsssab * ca * rzab/rab + 0.5_dp * dtau(1) * &
     &                 (dhsss*ojk(1)+dosss*hsssjk)*(-rzik/rik) )

!     &                ( (DHSSS*OJK(1)+DOSSS*HSSSJK)*(-RZIK/RIK)
!     &                 -(DHSSS*OSSS*OAB(1)*(-RZIK/RIK)+
!     &                   HSSS*DOSSS*OAB(1)*(-RZIK/RIK)+
!     &                   HSSS*OSSS*DOAB(1)*(-RZAB/RAB))
!     &                 -(DOAB(1)*(-RZAB/RAB)*OJK(1)*HSSSJK) ))

               endif

! *** SPS

               if (hspsab /= 0.0) then
                  ca = 0.5_dp*dtau(2)/hspsab* &
     &              (hsss*ojk(3)+osss*hpssjk)*gpsj

!     &              ((HSSS*OJK(3)+OSSS*HPSSJK)*GPSJ
!     &              -HSSS*OSSS*OAB(2)
!     &              -OAB(2)*OJK(3)*HPSSJK*GPSJ*GPSJ)

                  dc1a(2,1) = dc1a(2,1) + 1.0_dp/hspsab * &
     &              ( dhspsab * ca * rxab/rab + 0.5_dp * dtau(2) * &
     &                ( (dhsss*ojk(3)+dosss*hpssjk)*(-rxik/rik)*gpsj &
     &                 +(hsss*ojk(3)+osss*hpssjk)*dgpsj(1) ))

!     &                ( (DHSSS*OJK(3)+DOSSS*HPSSJK)*(-RXIK/RIK)*GPSJ
!     &                 +(HSSS*OJK(3)+OSSS*HPSSJK)*DGPSJ(1)
!     &                 -(DHSSS*OSSS*OAB(2)*(-RXIK/RIK)+
!     &                   HSSS*DOSSS*OAB(2)*(-RXIK/RIK)+
!     &                   HSSS*OSSS*DOAB(2)*(-RXAB/RAB))
!     &                 -(DOAB(2)*(-RXAB/RAB)*OJK(3)*HSSSJK*GPSJ*GPSJ)
!     &                 -(OAB(2)*OJK(3)*HPSSJK*2*GPSJ*DGPSJ(1)) ))

                  dc1a(2,2) = dc1a(2,2) + 1.0_dp/hspsab * &
     &              ( dhspsab * ca * ryab/rab + 0.5_dp * dtau(2) * &
     &                ( (dhsss*ojk(3)+dosss*hpssjk)*(-ryik/rik)*gpsj &
     &                 +(hsss*ojk(3)+osss*hpssjk)*dgpsj(2) ))

!     &                ( (DHSSS*OJK(3)+DOSSS*HPSSJK)*(-RYIK/RIK)*GPSJ
!     &                 +(HSSS*OJK(3)+OSSS*HPSSJK)*DGPSJ(2)
!     &                 -(DHSSS*OSSS*OAB(2)*(-RYIK/RIK)+
!     &                   HSSS*DOSSS*OAB(2)*(-RYIK/RIK)+
!     &                   HSSS*OSSS*DOAB(2)*(-RYAB/RAB))
!     &                 -(DOAB(2)*(-RYAB/RAB)*OJK(3)*HSSSJK*GPSJ*GPSJ)
!     &                 -(OAB(2)*OJK(3)*HPSSJK*2*GPSJ*DGPSJ(2)) ))

                  dc1a(2,3) = dc1a(2,3) + 1.0_dp/hspsab * &
     &              ( dhspsab * ca * rzab/rab + 0.5_dp * dtau(2) * &
     &                ( (dhsss*ojk(3)+dosss*hpssjk)*(-rzik/rik)*gpsj &
     &                 +(hsss*ojk(3)+osss*hpssjk)*dgpsj(3) ))

!     &                ( (DHSSS*OJK(3)+DOSSS*HPSSJK)*(-RZIK/RIK)*GPSJ
!     &                 +(HSSS*OJK(3)+OSSS*HPSSJK)*DGPSJ(3)
!     &                 -(DHSSS*OSSS*OAB(2)*(-RZIK/RIK)+
!     &                   HSSS*DOSSS*OAB(2)*(-RZIK/RIK)+
!     &                   HSSS*OSSS*DOAB(2)*(-RZAB/RAB))
!     &                 -(DOAB(2)*(-RZAB/RAB)*OJK(3)*HSSSJK*GPSJ*GPSJ)
!     &                 -(OAB(2)*OJK(3)*HPSSJK*2*GPSJ*DGPSJ(3)) ))

               endif

! *** SDS

               if (hsdsab /= 0.0) then
                  ca = 0.5_dp*dtau(6)/hsdsab* &
     &              (hsss*ojk(7)+osss*hdssjk)*gdsj

!     &              ((HSSS*OJK(7)+OSSS*HDSSJK)*GDSJ
!     &              -HSSS*OSSS*OAB(6)
!     &              -OAB(6)*OJK(7)*HDSSJK*GDSJ*GDSJ)

                  dc1a(6,1) = dc1a(6,1) + 1.0_dp/hsdsab * &
     &              ( dhsdsab * ca * rxab/rab + 0.5_dp * dtau(6) * &
     &                ( (dhsss*ojk(7)+dosss*hdssjk)*(-rxik/rik)*gdsj &
     &                 +(hsss*ojk(7)+osss*hdssjk)*dgdsj(1) ))

!     &                ( (DHSSS*OJK(7)+DOSSS*HDSSJK)*(-RXIK/RIK)*GDSJ
!     &                 +(HSSS*OJK(7)+OSSS*HDSSJK)*DGDSJ(1)
!     &                 -(DHSSS*OSSS*OAB(6)*(-RXIK/RIK)+
!     &                   HSSS*DOSSS*OAB(6)*(-RXIK/RIK)+
!     &                   HSSS*OSSS*DOAB(6)*(-RXAB/RAB))
!     &                 -(DOAB(6)*(-RXAB/RAB)*OJK(7)*HDSSJK*GDSJ*GDSJ)
!     &                 -(OAB(6)*OJK(7)*HDSSJK*2*GDSJ*DGDSJ(1)) ))

                  dc1a(6,2) = dc1a(6,2) + 1.0_dp/hsdsab * &
     &              ( dhsdsab * ca * ryab/rab + 0.5_dp * dtau(6) * &
     &                ( (dhsss*ojk(7)+dosss*hdssjk)*(-ryik/rik)*gdsj &
     &                 +(hsss*ojk(7)+osss*hdssjk)*dgdsj(2) ))

!     &                ( (DHSSS*OJK(7)+DOSSS*HDSSJK)*(-RYIK/RIK)*GDSJ
!     &                 +(HSSS*OJK(7)+OSSS*HDSSJK)*DGDSJ(2)
!     &                 -(DHSSS*OSSS*OAB(6)*(-RYIK/RIK)+
!     &                   HSSS*DOSSS*OAB(6)*(-RYIK/RIK)+
!     &                   HSSS*OSSS*DOAB(6)*(-RYAB/RAB))
!     &                 -(DOAB(6)*(-RYAB/RAB)*OJK(7)*HDSSJK*GDSJ*GDSJ)
!     &                 -(OAB(6)*OJK(7)*HDSSJK*2*GDSJ*DGDSJ(2)) ))

                  dc1a(6,3) = dc1a(6,3) + 1.0_dp/hsdsab * &
     &              ( dhsdsab * ca * rzab/rab + 0.5_dp * dtau(6) * &
     &                ( (dhsss*ojk(7)+dosss*hdssjk)*(-rzik/rik)*gdsj &
     &                 +(hsss*ojk(7)+osss*hdssjk)*dgdsj(3) ))

!     &                ( (DHSSS*OJK(7)+DOSSS*HDSSJK)*(-RZIK/RIK)*GDSJ
!     &                 +(HSSS*OJK(7)+OSSS*HDSSJK)*DGDSJ(3)
!     &                 -(DHSSS*OSSS*OAB(6)*(-RZIK/RIK)+
!     &                   HSSS*DOSSS*OAB(6)*(-RZIK/RIK)+
!     &                   HSSS*OSSS*DOAB(6)*(-RZAB/RAB))
!     &                 -(DOAB(6)*(-RZAB/RAB)*OJK(7)*HDSSJK*GDSJ*GDSJ)
!     &                 -(OAB(6)*OJK(7)*HDSSJK*2*GDSJ*DGDSJ(3)) ))

               endif

            endif

            endif


! Check whether we have p electrons.

            if (bndscl(1,3,bt) == 0.0) then
               opss = 0.0_dp
               dopss = 0.0_dp
            else
               opss = vscale(rik,bndscl(1,3,bt),bndscl(2,3,bt), & 
     &                      bndscl(5,3,bt),bndscl(6,3,bt), & 
     &                      bndscl(3,3,bt),bndscl(4,3,bt), & 
     &                      bndscl(7,3,bt),bndscl(8,3,bt), & 
     &                      bndscl(9,3,bt),bndscl(10,3,bt), & 
     &                      bndscl(11,3,bt),bndscl(12,3,bt))

               opss = opss*bndscl(13,3,bt)


               dopss = dscale(rik,bndscl(1,3,bt),bndscl(2,3,bt), & 
     &                       bndscl(5,3,bt),bndscl(6,3,bt), & 
     &                       bndscl(3,3,bt),bndscl(4,3,bt), & 
     &                       bndscl(7,3,bt),bndscl(8,3,bt), & 
     &                       bndscl(9,3,bt),bndscl(10,3,bt), & 
     &                       bndscl(11,3,bt),bndscl(12,3,bt))

               dopss = dopss*bndscl(13,3,bt)

! Calculate second moment contributions - MU2A(I)

               mu2a(3) = mu2a(3) + opss*opss*gps*gps
               mu2a(4) = mu2a(3)
               mu2a(5) = mu2a(5) + 0.5_dp*opss*opss*gpp*gpp
               mu2a(8) = mu2a(3)
               mu2a(10) = mu2a(5)

! Calculate derivatives of MU2A(I)

! * PSS        MU2A(3) = MU2A(3) + OPSS*OPSS*GPS*GPS

               dmu2a(3,1) = dmu2a(3,1) + 2.0_dp*opss*gps* &
     &                      (dopss*(-rxik/rik)*gps+opss*dgps(1))
               dmu2a(3,2) = dmu2a(3,2) + 2.0_dp*opss*gps* &
     &                      (dopss*(-ryik/rik)*gps+opss*dgps(2))
               dmu2a(3,3) = dmu2a(3,3) + 2.0_dp*opss*gps* &
     &                      (dopss*(-rzik/rik)*gps+opss*dgps(3))


! * PPS        MU2A(4) = MU2A(4) + OPSS*OPSS*GPS*GPS

               dmu2a(4,1) = dmu2a(3,1)
               dmu2a(4,2) = dmu2a(3,2)
               dmu2a(4,3) = dmu2a(3,3)


! * PPP        MU2A(5) = MU2A(5) + 0.5*OPSS*OPSS*GPP*GPP

               dmu2a(5,1) = dmu2a(5,1) + opss*gpp* &
     &                      (dopss*(-rxik/rik)*gpp+opss*dgpp(1))
               dmu2a(5,2) = dmu2a(5,2) + opss*gpp* &
     &                      (dopss*(-ryik/rik)*gpp+opss*dgpp(2))
               dmu2a(5,3) = dmu2a(5,3) + opss*gpp* &
     &                      (dopss*(-rzik/rik)*gpp+opss*dgpp(3))


! * PDS        MU2A(8) = MU2A(8) + OPSS*OPSS*GPS*GPS

               dmu2a(8,1) = dmu2a(3,1)
               dmu2a(8,2) = dmu2a(3,2)
               dmu2a(8,3) = dmu2a(3,3)


! * PDP        MU2A(10) = MU2A(10) + 0.5*OPSS*OPSS*GPP*GPP

               dmu2a(10,1) = dmu2a(5,1)
               dmu2a(10,2) = dmu2a(5,2)
               dmu2a(10,3) = dmu2a(5,3)

!    ___  End of dMU2A(I)  ___


! Calculate derivatives C1BI(I)

! *** PSS
               if (hpssab /= 0.0) then
                  cb = 0.5_dp*dtau(3)/hpssab* &
     &              (-hpss)*opss*oab(3)*gps*gps

                  dc1bi(3,1) = dc1bi(3,1) + 1.0_dp/hpssab * &
     &              ( dhpssab * cb * rxab/rab + 0.5_dp * dtau(3) * &
     &                ( -(dhpss*opss*oab(3)*(-rxik/rik)+ &
     &                   hpss*dopss*oab(3)*(-rxik/rik)+ &
     &                   hpss*opss*doab(3)*(-rxab/rab))*gps*gps ))

                  dc1bi(3,2) = dc1bi(3,2) + 1.0_dp/hpssab * &
     &              ( dhpssab * cb * ryab/rab + 0.5_dp * dtau(3) * &
     &                ( -(dhpss*opss*oab(3)*(-ryik/rik)+ &
     &                   hpss*dopss*oab(3)*(-ryik/rik)+ &
     &                   hpss*opss*doab(3)*(-ryab/rab))*gps*gps ))

                  dc1bi(3,3) = dc1bi(3,3) + 1.0_dp/hpssab * &
     &              ( dhpssab * cb * rzab/rab + 0.5_dp * dtau(3) * &
     &                ( -(dhpss*opss*oab(3)*(-rzik/rik)+ &
     &                   hpss*dopss*oab(3)*(-rzik/rik)+ &
     &                   hpss*opss*doab(3)*(-rzab/rab))*gps*gps ))

               endif

! *** PPS
               if (hppsab /= 0.0) then
                  cb = 0.5_dp*dtau(4)/hppsab* &
     &              (-hpss)*opss*oab(4)*gps*gps

                  dc1bi(4,1) = dc1bi(4,1) + 1.0_dp/hppsab * &
     &              ( dhppsab * cb * rxab/rab + 0.5_dp * dtau(4) * &
     &                ( -(dhpss*opss*oab(4)*(-rxik/rik)+ &
     &                   hpss*dopss*oab(4)*(-rxik/rik)+ &
     &                   hpss*opss*doab(4)*(-rxab/rab))*gps*gps &
     &                 -hpss*opss*oab(4)*2*gps*dgps(1) ))

                  dc1bi(4,2) = dc1bi(4,2) + 1.0_dp/hppsab * &
     &              ( dhppsab * cb * ryab/rab + 0.5_dp * dtau(4) * &
     &                ( -(dhpss*opss*oab(4)*(-ryik/rik)+ &
     &                   hpss*dopss*oab(4)*(-ryik/rik)+ &
     &                   hpss*opss*doab(4)*(-ryab/rab))*gps*gps &
     &                 -hpss*opss*oab(4)*2*gps*dgps(2) ))

                  dc1bi(4,3) = dc1bi(4,3) + 1.0_dp/hppsab * &
     &              ( dhppsab * cb * rzab/rab + 0.5_dp * dtau(4) * &
     &                ( -(dhpss*opss*oab(4)*(-rzik/rik)+ &
     &                   hpss*dopss*oab(4)*(-rzik/rik)+ &
     &                   hpss*opss*doab(4)*(-rzab/rab))*gps*gps &
     &                 -hpss*opss*oab(4)*2*gps*dgps(3) ))

               endif

! *** PPP
               if (hpppab /= 0.0) then
                  cb = -0.25_dp*dtau(5)/hpppab* &
     &              (-hpss)*opss*oab(5)*gpp*gpp

                  dc1bi(5,1) = dc1bi(5,1) + 1.0_dp/hpppab * &
     &              ( dhpppab * cb * rxab/rab - 0.25_dp * dtau(5) * &
     &                ( -(dhpss*opss*oab(5)*(-rxik/rik)+ &
     &                   hpss*dopss*oab(5)*(-rxik/rik)+ &
     &                   hpss*opss*doab(5)*(-rxab/rab))*gpp*gpp &
     &                 -hpss*opss*oab(5)*2.0_dp*gpp*dgpp(1) ))

                  dc1bi(5,2) = dc1bi(5,2) + 1.0_dp/hpppab * &
     &              ( dhpppab * cb * ryab/rab - 0.25_dp * dtau(5) * &
     &                ( -(dhpss*opss*oab(5)*(-ryik/rik)+ &
     &                   hpss*dopss*oab(5)*(-ryik/rik)+ &
     &                   hpss*opss*doab(5)*(-ryab/rab))*gpp*gpp &
     &                 -hpss*opss*oab(5)*2.0_dp*gpp*dgpp(2) ))

                  dc1bi(5,3) = dc1bi(5,3) + 1.0_dp/hpppab * &
     &              ( dhpppab * cb * rzab/rab - 0.25_dp * dtau(5) * &
     &                ( -(dhpss*opss*oab(5)*(-rzik/rik)+ &
     &                   hpss*dopss*oab(5)*(-rzik/rik)+ &
     &                   hpss*opss*doab(5)*(-rzab/rab))*gpp*gpp &
     &                 -hpss*opss*oab(5)*2.0_dp*gpp*dgpp(3) ))

               endif

! *** PDS
               if (hpdsab /= 0.0) then
                  cb = 0.5_dp*dtau(8)/hpdsab* &
     &              (-hpss)*opss*oab(8)*gps*gps

                  dc1bi(8,1) = dc1bi(8,1) + 1.0_dp/hpdsab * &
     &              ( dhpdsab * cb * rxab/rab + 0.5_dp * dtau(8) * &
     &                ( -(dhpss*opss*oab(8)*(-rxik/rik)+ &
     &                   hpss*dopss*oab(8)*(-rxik/rik)+ &
     &                   hpss*opss*doab(8)*(-rxab/rab))*gps*gps &
     &                 -hpss*opss*oab(8)*2.0_dp*gps*dgps(1) ))

                  dc1bi(8,2) = dc1bi(8,2) + 1.0_dp/hpdsab * &
     &              ( dhpdsab * cb * ryab/rab + 0.5_dp * dtau(8) * &
     &                ( -(dhpss*opss*oab(8)*(-ryik/rik)+ &
     &                   hpss*dopss*oab(8)*(-ryik/rik)+ &
     &                   hpss*opss*doab(8)*(-ryab/rab))*gps*gps &
     &                 -hpss*opss*oab(8)*2.0_dp*gps*dgps(2) ))

                  dc1bi(8,3) = dc1bi(8,3) + 1.0_dp/hpdsab * &
     &              ( dhpdsab * cb * rzab/rab + 0.5_dp * dtau(8) * &
     &                ( -(dhpss*opss*oab(8)*(-rzik/rik)+ &
     &                   hpss*dopss*oab(8)*(-rzik/rik)+ &
     &                   hpss*opss*doab(8)*(-rzab/rab))*gps*gps &
     &                 -hpss*opss*oab(8)*2.0_dp*gps*dgps(3) ))

               endif

! *** PDP
               if (hpdpab /= 0.0) then
                  cb = -0.25_dp*dtau(10)/hpdpab* &
     &              (-hpss)*opss*oab(10)*gpp*gpp

                  dc1bi(10,1) = dc1bi(10,1) + 1.0_dp/hpdpab * &
     &              ( dhpdpab * cb * rxab/rab - 0.25_dp * dtau(10) * &
     &                ( -(dhpss*opss*oab(10)*(-rxik/rik)+ &
     &                   hpss*dopss*oab(10)*(-rxik/rik)+ &
     &                   hpss*opss*doab(10)*(-rxab/rab))*gpp*gpp &
     &                 -hpss*opss*oab(10)*2.0_dp*gpp*dgpp(1) ))

                  dc1bi(10,2) = dc1bi(10,2) + 1.0_dp/hpdpab * &
     &              ( dhpdpab * cb * ryab/rab - 0.25_dp * dtau(10) * &
     &                ( -(dhpss*opss*oab(10)*(-ryik/rik)+ &
     &                   hpss*dopss*oab(10)*(-ryik/rik)+ &
     &                   hpss*opss*doab(10)*(-ryab/rab))*gpp*gpp &
     &                 -hpss*opss*oab(10)*2.0_dp*gpp*dgpp(2) ))

                  dc1bi(10,3) = dc1bi(10,3) + 1.0_dp/hpdpab * &
     &              ( dhpdpab * cb * rzab/rab - 0.25_dp * dtau(10) * &
     &                ( -(dhpss*opss*oab(10)*(-rzik/rik)+ &
     &                   hpss*dopss*oab(10)*(-rzik/rik)+ &
     &                   hpss*opss*doab(10)*(-rzab/rab))*gpp*gpp &
     &                 -hpss*opss*oab(10)*2.0_dp*gpp*dgpp(3) ))

               endif

!     ___  End of dC1BI(I)  ___


! Calculate C1A(I-J) derivatives.

               if (rjk < rcut) then

! *** PSS
               if (hpssab /= 0.0) then
                  ca = 0.5_dp*dtau(3)/hpssab* &
     &              (hpss*ojk(1)+opss*hsssjk)*gps

                  dc1a(3,1) = dc1a(3,1) + 1.0_dp/hpssab * &
     &              ( dhpssab * ca * rxab/rab + 0.5_dp * dtau(3) * &
     &                ( (dhpss*ojk(1)+dopss*hsssjk)*(-rxik/rik)*gps &
     &                 +(hpss*ojk(1)+opss*hsssjk)*dgps(1) ))

!     &                 -(DHPSS*OPSS*OAB(3)*(-RXIK/RIK)+
!     &                   HPSS*DOPSS*OAB(3)*(-RXIK/RIK)+
!     &                   HPSS*OPSS*DOAB(3)*(-RXAB/RAB))*GPS*GPS
!     &                 -HPSS*OPSS*OAB(3)*2*GPS*DGPS(1)
!     &                 -DOAB(3)*(-RXAB/RAB)*OJK(1)*HSSSJK ))

                  dc1a(3,2) = dc1a(3,2) + 1.0_dp/hpssab * &
     &              ( dhpssab * ca * ryab/rab + 0.5_dp * dtau(3) * &
     &                ( (dhpss*ojk(1)+dopss*hsssjk)*(-rxik/rik)*gps &
     &                 +(hpss*ojk(1)+opss*hsssjk)*dgps(2) ))

!     &                 -(DHPSS*OPSS*OAB(3)*(-RYIK/RIK)+
!     &                   HPSS*DOPSS*OAB(3)*(-RYIK/RIK)+
!     &                   HPSS*OPSS*DOAB(3)*(-RYAB/RAB))*GPS*GPS
!     &                 -HPSS*OPSS*OAB(3)*2*GPS*DGPS(2)
!     &                 -DOAB(3)*(-RYAB/RAB)*OJK(1)*HSSSJK ))

                  dc1a(3,3) = dc1a(3,3) + 1.0_dp/hpssab * &
     &              ( dhpssab * ca * rzab/rab + 0.5_dp * dtau(3) * &
     &                ( (dhpss*ojk(1)+dopss*hsssjk)*(-rzik/rik)*gps &
     &                 +(hpss*ojk(1)+opss*hsssjk)*dgps(3) ))

!     &                 -(DHPSS*OPSS*OAB(3)*(-RZIK/RIK)+
!     &                   HPSS*DOPSS*OAB(3)*(-RZIK/RIK)+
!     &                   HPSS*OPSS*DOAB(3)*(-RZAB/RAB))*GPS*GPS
!     &                 -HPSS*OPSS*OAB(3)*2*GPS*DGPS(3)
!     &                 -DOAB(3)*(-RZAB/RAB)*OJK(1)*HSSSJK ))

               endif

! *** PPS
               if (hppsab /= 0.0) then
                  ca = -0.5_dp*dtau(4)/hppsab* &
     &              (hpss*ojk(3)+opss*hpssjk)*gps*gpsj

                  dc1a(4,1) = dc1a(4,1) + 1.0_dp/hppsab * &
     &              ( dhppsab * ca * rxab/rab - 0.5_dp * dtau(4) * &
     &              ( (dhpss*ojk(3)+dopss*hpssjk)*(-rxik/rik)*gps*gpsj &
     &                 +(hpss*ojk(3)+opss*hpssjk)* &
     &                   (dgps(1)*gpsj+gps*dgpsj(1)) ))

!     &                 -(DHPSS*OPSS*OAB(4)*(-RXIK/RIK)+
!     &                   HPSS*DOPSS*OAB(4)*(-RXIK/RIK)+
!     &                   HPSS*OPSS*DOAB(4)*(-RXAB/RAB))*GPS*GPS
!     &                 -HPSS*OPSS*OAB(4)*2*GPS*DGPS(1)
!     &                 -DOAB(4)*(-RXAB/RAB)*OJK(3)*HPSSJK*GPSJ*GPSJ
!     &                 -OAB(4)*OJK(3)*HPSSJK*2*GPSJ*DGPSJ(1) ))

                  dc1a(4,2) = dc1a(4,2) + 1.0_dp/hppsab * &
     &              ( dhppsab * ca * ryab/rab - 0.5_dp * dtau(4) * &
     &                ( (dhpss*ojk(3)+dopss*hpssjk)*(-ryik/rik)*gps*gpsj &
     &                 +(hpss*ojk(3)+opss*hpssjk)* &
     &                   (dgps(2)*gpsj+gps*dgpsj(2)) ))

!     &                 -(DHPSS*OPSS*OAB(4)*(-RYIK/RIK)+
!     &                   HPSS*DOPSS*OAB(4)*(-RYIK/RIK)+
!     &                   HPSS*OPSS*DOAB(4)*(-RYAB/RAB))*GPS*GPS
!     &                 -HPSS*OPSS*OAB(4)*2*GPS*DGPS(2)
!     &                 -DOAB(4)*(-RYAB/RAB)*OJK(3)*HPSSJK*GPSJ*GPSJ
!     &                 -OAB(4)*OJK(3)*HPSSJK*2*GPSJ*DGPSJ(2) ))

                  dc1a(4,3) = dc1a(4,3) + 1.0_dp/hppsab * &
     &              ( dhppsab * ca * rzab/rab - 0.5_dp * dtau(4) * &
     &                ( (dhpss*ojk(3)+dopss*hpssjk)*(-rzik/rik)*gps*gpsj &
     &                 +(hpss*ojk(3)+opss*hpssjk)* &
     &                   (dgps(3)*gpsj+gps*dgpsj(3)) ))

!     &                 -(DHPSS*OPSS*OAB(4)*(-RZIK/RIK)+
!     &                   HPSS*DOPSS*OAB(4)*(-RZIK/RIK)+
!     &                   HPSS*OPSS*DOAB(4)*(-RZAB/RAB))*GPS*GPS
!     &                 -HPSS*OPSS*OAB(4)*2*GPS*DGPS(3)
!     &                 -DOAB(4)*(-RZAB/RAB)*OJK(3)*HPSSJK*GPSJ*GPSJ
!     &                 -OAB(4)*OJK(3)*HPSSJK*2*GPSJ*DGPSJ(3) ))

               endif

! *** PPP
               if (hpppab /= 0.0) then
                  ca = 0.25_dp*dtau(5)/hpppab* &
     &              (hpss*ojk(3)+opss*hpssjk)*gpp*gppj

                  dc1a(5,1) = dc1a(5,1) + 1.0_dp/hpppab * &
     &              ( dhpppab * ca * rxab/rab + 0.25_dp * dtau(5) * &
     &                ( (dhpss*ojk(3)+dopss*hpssjk)*(-rxik/rik)*gpp*gppj &
     &                 +(hpss*ojk(3)+opss*hpssjk)* &
     &                   (dgpp(1)*gppj+gpp*dgppj(1)) ))

!     &                 -(DHPSS*OPSS*OAB(5)*(-RXIK/RIK)+
!     &                   HPSS*DOPSS*OAB(5)*(-RXIK/RIK)+
!     &                   HPSS*OPSS*DOAB(5)*(-RXAB/RAB))*GPP*GPP
!     &                 -HPSS*OPSS*OAB(5)*2*GPP*DGPP(1)
!     &                 -DOAB(5)*(-RXAB/RAB)*OJK(3)*HPSSJK*GPPJ*GPPJ
!     &                 -OAB(5)*OJK(3)*HPSSJK*2*GPPJ*DGPPJ(1) ))

                  dc1a(5,2) = dc1a(5,2) + 1.0_dp/hpppab * &
     &              ( dhpppab * ca * ryab/rab + 0.25_dp * dtau(5) * &
     &                ( (dhpss*ojk(3)+dopss*hpssjk)*(-ryik/rik)*gpp*gppj &
     &                 +(hpss*ojk(3)+opss*hpssjk)* &
     &                   (dgpp(2)*gppj+gpp*dgppj(2)) ))

!     &                 -(DHPSS*OPSS*OAB(5)*(-RYIK/RIK)+
!     &                   HPSS*DOPSS*OAB(5)*(-RYIK/RIK)+
!     &                   HPSS*OPSS*DOAB(5)*(-RYAB/RAB))*GPP*GPP
!     &                 -HPSS*OPSS*OAB(5)*2*GPP*DGPP(2)
!     &                 -DOAB(5)*(-RYAB/RAB)*OJK(3)*HPSSJK*GPPJ*GPPJ
!     &                 -OAB(5)*OJK(3)*HPSSJK*2*GPPJ*DGPPJ(2) ))

                  dc1a(5,3) = dc1a(5,3) + 1.0_dp/hpppab * &
     &              ( dhpppab * ca * rzab/rab + 0.25_dp * dtau(5) * &
     &                ( (dhpss*ojk(3)+dopss*hpssjk)*(-rzik/rik)*gpp*gppj &
     &                 +(hpss*ojk(3)+opss*hpssjk)* &
     &                   (dgpp(3)*gppj+gpp*dgppj(3)) ))

!     &                 -(DHPSS*OPSS*OAB(5)*(-RZIK/RIK)+
!     &                   HPSS*DOPSS*OAB(5)*(-RZIK/RIK)+
!     &                   HPSS*OPSS*DOAB(5)*(-RZAB/RAB))*GPP*GPP
!     &                 -HPSS*OPSS*OAB(5)*2*GPP*DGPP(3)
!     &                 -DOAB(5)*(-RZAB/RAB)*OJK(3)*HPSSJK*GPPJ*GPPJ
!     &                 -OAB(5)*OJK(3)*HPSSJK*2*GPPJ*DGPPJ(3) ))

               endif

! *** PDS
               if (hpdsab /= 0.0) then
                  ca = 0.5_dp*dtau(8)/hpdsab* &
     &              (hpss*ojk(7)+opss*hdssjk)*gps*gdsj

                  dc1a(8,1) = dc1a(8,1) + 1.0_dp/hpdsab * &
     &              ( dhpdsab * ca * rxab/rab + 0.5_dp * dtau(8) * &
     &                ( (dhpss*ojk(7)+dopss*hdssjk)*(-rxik/rik)*gps*gdsj &
     &                 +(hpss*ojk(7)+opss*hdssjk)* &
     &                   (dgps(1)*gdsj+gps*dgdsj(1)) ))

!     &                 -(DHPSS*OPSS*OAB(8)*(-RXIK/RIK)+
!     &                   HPSS*DOPSS*OAB(8)*(-RXIK/RIK)+
!     &                   HPSS*OPSS*DOAB(8)*(-RXAB/RAB))*GPS*GPS
!     &                 -HPSS*OPSS*OAB(8)*2*GPS*DGPS(1)
!     &                 -DOAB(8)*(-RXAB/RAB)*OJK(7)*HDSSJK*GDSJ*GDSJ
!     &                 -OAB(8)*OJK(7)*HDSSJK*2*GDSJ*DGDSJ(1) ))

                  dc1a(8,2) = dc1a(8,2) + 1.0_dp/hpdsab * &
     &              ( dhpdsab * ca * ryab/rab + 0.5_dp * dtau(8) * &
     &                ( (dhpss*ojk(7)+dopss*hdssjk)*(-ryik/rik)*gps*gdsj &
     &                 +(hpss*ojk(7)+opss*hdssjk)* &
     &                   (dgps(2)*gdsj+gps*dgdsj(2)) ))

!     &                 -(DHPSS*OPSS*OAB(8)*(-RYIK/RIK)+
!     &                   HPSS*DOPSS*OAB(8)*(-RYIK/RIK)+
!     &                   HPSS*OPSS*DOAB(8)*(-RYAB/RAB))*GPS*GPS
!     &                 -HPSS*OPSS*OAB(8)*2*GPS*DGPS(2)
!     &                 -DOAB(8)*(-RYAB/RAB)*OJK(7)*HDSSJK*GDSJ*GDSJ
!     &                 -OAB(8)*OJK(7)*HDSSJK*2*GDSJ*DGDSJ(2) ))

                  dc1a(8,3) = dc1a(8,3) + 1.0_dp/hpdsab * &
     &              ( dhpdsab * ca * rzab/rab + 0.5_dp * dtau(8) * &
     &                ( (dhpss*ojk(7)+dopss*hdssjk)*(-rzik/rik)*gps*gdsj &
     &                 +(hpss*ojk(7)+opss*hdssjk)* &
     &                   (dgps(3)*gdsj+gps*dgdsj(3)) ))

!     &                 -(DHPSS*OPSS*OAB(8)*(-RZIK/RIK)+
!     &                   HPSS*DOPSS*OAB(8)*(-RZIK/RIK)+
!     &                   HPSS*OPSS*DOAB(8)*(-RZAB/RAB))*GPS*GPS
!     &                 -HPSS*OPSS*OAB(8)*2*GPS*DGPS(3)
!     &                 -DOAB(8)*(-RZAB/RAB)*OJK(7)*HDSSJK*GDSJ*GDSJ
!     &                 -OAB(8)*OJK(7)*HDSSJK*2*GDSJ*DGDSJ(3) ))

               endif

! *** PDP
               if (hpdpab /= 0.0) then
                  ca = -0.25_dp*dtau(10)/hpdpab* &
     &              (hpss*ojk(7)+opss*hdssjk)*gpp*gdpj

                  dc1a(10,1) = dc1a(10,1) + 1.0_dp/hpdpab * &
     &              ( dhpdpab * ca * rxab/rab - 0.25_dp * dtau(10) * &
     &                ( (dhpss*ojk(7)+dopss*hdssjk)*(-rxik/rik)*gpp*gdpj &
     &                 +(hpss*ojk(7)+opss*hdssjk)* &
     &                   (dgpp(1)*gdpj+gpp*dgdpj(1)) ))

!     &                 -(DHPSS*OPSS*OAB(10)*(-RXIK/RIK)+
!     &                   HPSS*DOPSS*OAB(10)*(-RXIK/RIK)+
!     &                   HPSS*OPSS*DOAB(10)*(-RXAB/RAB))*GPP*GPP
!     &                 -HPSS*OPSS*OAB(10)*2*GPP*DGPP(1)
!     &                 -DOAB(10)*(-RXAB/RAB)*OJK(7)*HDSSJK*GDPJ*GDPJ
!     &                 -OAB(10)*OJK(7)*HDSSJK*2*GDPJ*DGDPJ(1) ))

                  dc1a(10,2) = dc1a(10,2) + 1.0_dp/hpdpab * &
     &              ( dhpdpab * ca * ryab/rab - 0.25_dp * dtau(10) * &
     &                ( (dhpss*ojk(7)+dopss*hdssjk)*(-ryik/rik)*gpp*gdpj &
     &                 +(hpss*ojk(7)+opss*hdssjk)* &
     &                   (dgpp(2)*gdpj+gpp*dgdpj(2)) ))

!     &                 -(DHPSS*OPSS*OAB(10)*(-RYIK/RIK)+
!     &                   HPSS*DOPSS*OAB(10)*(-RYIK/RIK)+
!     &                   HPSS*OPSS*DOAB(10)*(-RYAB/RAB))*GPP*GPP
!     &                 -HPSS*OPSS*OAB(10)*2*GPP*DGPP(2)
!     &                 -DOAB(10)*(-RYAB/RAB)*OJK(7)*HDSSJK*GDPJ*GDPJ
!     &                 -OAB(10)*OJK(7)*HDSSJK*2*GDPJ*DGDPJ(2) ))

                  dc1a(10,3) = dc1a(10,3) + 1.0_dp/hpdpab * &
     &              ( dhpdpab * ca * rzab/rab - 0.25_dp * dtau(10) * &
     &                ( (dhpss*ojk(7)+dopss*hdssjk)*(-rzik/rik)*gpp*gdpj &
     &                 +(hpss*ojk(7)+opss*hdssjk)* &
     &                   (dgpp(3)*gdpj+gpp*dgdpj(3)) ))

!     &                 -(DHPSS*OPSS*OAB(10)*(-RZIK/RIK)+
!     &                   HPSS*DOPSS*OAB(10)*(-RZIK/RIK)+
!     &                   HPSS*OPSS*DOAB(10)*(-RZAB/RAB))*GPP*GPP
!     &                 -HPSS*OPSS*OAB(10)*2*GPP*DGPP(3)
!     &                 -DOAB(10)*(-RZAB/RAB)*OJK(7)*HDSSJK*GDPJ*GDPJ
!     &                 -OAB(10)*OJK(7)*HDSSJK*2*GDPJ*DGDPJ(3) ))

               endif

            endif

            endif


! Check whether we have d electrons.

            if (bndscl(1,7,bt) == 0.0) then
               odss = 0.0_dp
               dodss = 0.0_dp
            else
               odss = vscale(rik,bndscl(1,7,bt),bndscl(2,7,bt), & 
     &                      bndscl(5,7,bt),bndscl(6,7,bt), & 
     &                      bndscl(3,7,bt),bndscl(4,7,bt), & 
     &                      bndscl(7,7,bt),bndscl(8,7,bt), & 
     &                      bndscl(9,7,bt),bndscl(10,7,bt), & 
     &                      bndscl(11,7,bt),bndscl(12,7,bt))

               odss = odss*bndscl(13,7,bt)


               dodss = dscale(rik,bndscl(1,7,bt),bndscl(2,7,bt), & 
     &                       bndscl(5,7,bt),bndscl(6,7,bt), & 
     &                       bndscl(3,7,bt),bndscl(4,7,bt), & 
     &                       bndscl(7,7,bt),bndscl(8,7,bt), & 
     &                       bndscl(9,7,bt),bndscl(10,7,bt), & 
     &                       bndscl(11,7,bt),bndscl(12,7,bt))

               dodss = dodss*bndscl(13,7,bt)

! Calculate second moment contributions - MU2A(I)

               mu2a(7) = mu2a(7) + odss*odss*gds*gds
               mu2a(9) = mu2a(7)
               mu2a(11) = mu2a(11) + 0.5_dp*odss*odss*gdp*gdp
               mu2a(12) = mu2a(7)
               mu2a(13) = mu2a(11)
               mu2a(14) = mu2a(14) + 0.5_dp*odss*odss*gdd*gdd

! Calculate derivatives of MU2A(I)

! * DSS        MU2A(7) = MU2A(7) + ODSS*ODSS*GDS*GDS

               dmu2a(7,1) = dmu2a(7,1) + 2.0_dp*odss*gds* &
     &                      (dodss*(-rxik/rik)*gds+odss*dgds(1))
               dmu2a(7,2) = dmu2a(7,2) + 2.0_dp*odss*gds* &
     &                      (dodss*(-ryik/rik)*gds+odss*dgds(2))
               dmu2a(7,3) = dmu2a(7,3) + 2.0_dp*odss*gds* &
     &                      (dodss*(-rzik/rik)*gds+odss*dgds(3))


! * DPS        MU2A(9) = MU2A(9) + ODSS*ODSS*GDS*GDS

               dmu2a(9,1) = dmu2a(7,1)
               dmu2a(9,2) = dmu2a(7,2)
               dmu2a(9,3) = dmu2a(7,3)


! * DPP        MU2A(11) = MU2A(11) + 0.5*ODSS*ODSS*GDP*GDP

               dmu2a(11,1) = dmu2a(11,1) + odss*gdp* &
     &                      (dodss*(-rxik/rik)*gdp+odss*dgdp(1))
               dmu2a(11,2) = dmu2a(11,2) + odss*gdp* &
     &                      (dodss*(-ryik/rik)*gdp+odss*dgdp(2))
               dmu2a(11,3) = dmu2a(11,3) + odss*gdp* &
     &                      (dodss*(-rzik/rik)*gdp+odss*dgdp(3))


! * DDS        MU2A(12) = MU2A(12) + ODSS*ODSS*GDS*GDS

               dmu2a(12,1) = dmu2a(7,1)
               dmu2a(12,2) = dmu2a(7,2)
               dmu2a(12,3) = dmu2a(7,3)


! * DDP        MU2A(13) = MU2A(13) + 0.5*ODSS*ODSS*GDP*GDP

               dmu2a(13,1) = dmu2a(11,1)
               dmu2a(13,2) = dmu2a(11,2)
               dmu2a(13,3) = dmu2a(11,3)


! * DDD        MU2A(14) = MU2A(14) + 0.5*ODSS*ODSS*GDD*GDD

               dmu2a(14,1) = dmu2a(14,1) + odss*gdd* &
     &                      (dodss*(-rxik/rik)*gdd+odss*dgdd(1))
               dmu2a(14,2) = dmu2a(14,2) + odss*gdd* &
     &                      (dodss*(-ryik/rik)*gdd+odss*dgdd(2))
               dmu2a(14,3) = dmu2a(14,3) + odss*gdd* &
     &                      (dodss*(-rzik/rik)*gdd+odss*dgdd(3))


!    ___  End of dMU2A(I)  ___


! Calculate derivatives C1BI(I)

! *** DSS
               if (hdssab /= 0.0) then
                  cb = 0.5_dp*dtau(7)/hdssab* &
     &                 (-hdss)*odss*oab(7)*gds*gds

                  dc1bi(7,1) = dc1bi(7,1) + 1.0_dp/hdssab * &
     &              ( dhdssab * cb * rxab/rab + 0.5_dp * dtau(7) * &
     &                ( -(dhdss*odss*oab(7)*(-rxik/rik)+ &
     &                   hdss*dodss*oab(7)*(-rxik/rik)+ &
     &                   hdss*odss*doab(7)*(-rxab/rab))*gds*gds &
     &                 -hdss*odss*oab(7)*2.0_dp*gds*dgds(1) ))

                  dc1bi(7,2) = dc1bi(7,2) + 1.0_dp/hdssab * &
     &              ( dhdssab * cb * ryab/rab + 0.5_dp * dtau(7) * &
     &                ( -(dhdss*odss*oab(7)*(-ryik/rik)+ &
     &                   hdss*dodss*oab(7)*(-ryik/rik)+ &
     &                   hdss*odss*doab(7)*(-ryab/rab))*gds*gds &
     &                 -hdss*odss*oab(7)*2.0_dp*gds*dgds(2) ))

                  dc1bi(7,3) = dc1bi(7,3) + 1.0_dp/hdssab * &
     &              ( dhdssab * cb * rzab/rab + 0.5_dp * dtau(7) * &
     &                ( -(dhdss*odss*oab(7)*(-rzik/rik)+ &
     &                   hdss*dodss*oab(7)*(-rzik/rik)+ &
     &                   hdss*odss*doab(7)*(-rzab/rab))*gds*gds &
     &                 -hdss*odss*oab(7)*2.0_dp*gds*dgds(3) ))

               endif


! *** DPS
               if (hdpsab /= 0.0) then
                  cb = 0.5_dp*dtau(9)/hdpsab* &
     &                 (-hdss)*odss*oab(9)*gds*gds

                  dc1bi(9,1) = dc1bi(9,1) + 1.0_dp/hdpsab * &
     &              ( dhdpsab * cb * rxab/rab + 0.5_dp * dtau(9) * &
     &                ( -(dhdss*odss*oab(9)*(-rxik/rik)+ &
     &                   hdss*dodss*oab(9)*(-rxik/rik)+ &
     &                   hdss*odss*doab(9)*(-rxab/rab))*gds*gds &
     &                 -hdss*odss*oab(9)*2.0_dp*gds*dgds(1) ))

                  dc1bi(9,2) = dc1bi(9,2) + 1.0_dp/hdpsab * &
     &              ( dhdpsab * cb * ryab/rab + 0.5_dp * dtau(9) * &
     &                ( -(dhdss*odss*oab(9)*(-ryik/rik)+ &
     &                   hdss*dodss*oab(9)*(-ryik/rik)+ &
     &                   hdss*odss*doab(9)*(-ryab/rab))*gds*gds &
     &                 -hdss*odss*oab(9)*2.0_dp*gds*dgds(2) ))

                  dc1bi(9,3) = dc1bi(9,3) + 1.0_dp/hdpsab * &
     &              ( dhdpsab * cb * rzab/rab + 0.5_dp * dtau(9) * &
     &                ( -(dhdss*odss*oab(9)*(-rzik/rik)+ &
     &                   hdss*dodss*oab(9)*(-rzik/rik)+ &
     &                   hdss*odss*doab(9)*(-rzab/rab))*gds*gds &
     &                 -hdss*odss*oab(9)*2.0_dp*gds*dgds(3) ))

               endif


! *** DPP
               if (hdppab /= 0.0) then
                  cb = -0.25_dp*dtau(11)/hdppab* &
     &                 (-hdss)*odss*oab(11)*gdp*gdp

                  dc1bi(11,1) = dc1bi(11,1) + 1.0_dp/hdppab * &
     &              ( dhdppab * cb * rxab/rab - 0.25_dp * dtau(11) * &
     &                ( -(dhdss*odss*oab(11)*(-rxik/rik)+ &
     &                   hdss*dodss*oab(11)*(-rxik/rik)+ &
     &                   hdss*odss*doab(11)*(-rxab/rab))*gdp*gdp &
     &                 -hdss*odss*oab(11)*2.0_dp*gdp*dgdp(1) ))

                  dc1bi(11,2) = dc1bi(11,2) + 1.0_dp/hdppab * &
     &              ( dhdppab * cb * ryab/rab - 0.25_dp * dtau(11) * &
     &                ( -(dhdss*odss*oab(11)*(-ryik/rik)+ &
     &                   hdss*dodss*oab(11)*(-ryik/rik)+ &
     &                   hdss*odss*doab(11)*(-ryab/rab))*gdp*gdp &
     &                 -hdss*odss*oab(11)*2.0_dp*gdp*dgdp(2) ))

                  dc1bi(11,3) = dc1bi(11,3) + 1.0_dp/hdppab * &
     &              ( dhdppab * cb * rzab/rab - 0.25_dp * dtau(11) * &
     &                ( -(dhdss*odss*oab(11)*(-rzik/rik)+ &
     &                   hdss*dodss*oab(11)*(-rzik/rik)+ &
     &                   hdss*odss*doab(11)*(-rzab/rab))*gdp*gdp &
     &                 -hdss*odss*oab(11)*2.0_dp*gdp*dgdp(3) ))

               endif

! *** DDS
               if (hddsab /= 0.0) then
                  cb = 0.5_dp*dtau(12)/hddsab* &
     &                 (-hdss)*odss*oab(12)*gds*gds

                  dc1bi(12,1) = dc1bi(12,1) + 1.0_dp/hddsab * &
     &              ( dhddsab * cb * rxab/rab + 0.5_dp * dtau(12) * &
     &                ( -(dhdss*odss*oab(12)*(-rxik/rik)+ &
     &                   hdss*dodss*oab(12)*(-rxik/rik)+ &
     &                   hdss*odss*doab(12)*(-rxab/rab))*gds*gds &
     &                 -hdss*odss*oab(12)*2.0_dp*gds*dgds(1) ))

                  dc1bi(12,2) = dc1bi(12,2) + 1.0_dp/hddsab * &
     &              ( dhddsab * cb * ryab/rab + 0.5_dp * dtau(12) * &
     &                ( -(dhdss*odss*oab(12)*(-ryik/rik)+ &
     &                   hdss*dodss*oab(12)*(-ryik/rik)+ &
     &                   hdss*odss*doab(12)*(-ryab/rab))*gds*gds &
     &                 -hdss*odss*oab(12)*2.0_dp*gds*dgds(2) ))

                  dc1bi(12,3) = dc1bi(12,3) + 1.0_dp/hddsab * &
     &              ( dhddsab * cb * rzab/rab + 0.5_dp * dtau(12) * &
     &                ( -(dhdss*odss*oab(12)*(-rzik/rik)+ &
     &                   hdss*dodss*oab(12)*(-rzik/rik)+ &
     &                   hdss*odss*doab(12)*(-rzab/rab))*gds*gds &
     &                 -hdss*odss*oab(12)*2.0_dp*gds*dgds(3) ))

               endif


! *** DDP
               if (hddpab /= 0.0) then
                  cb = 0.25_dp*dtau(13)/hddpab* &
     &                 (-hdss)*odss*oab(13)*gdp*gdp

                  dc1bi(13,1) = dc1bi(13,1) + 1.0_dp/hddpab * &
     &              ( dhddpab * cb * rxab/rab + 0.25_dp * dtau(13) * &
     &                ( -(dhdss*odss*oab(13)*(-rxik/rik)+ &
     &                   hdss*dodss*oab(13)*(-rxik/rik)+ &
     &                   hdss*odss*doab(13)*(-rxab/rab))*gdp*gdp &
     &                 -hdss*odss*oab(13)*2.0_dp*gdp*dgdp(1) ))

                  dc1bi(13,2) = dc1bi(13,2) + 1.0_dp/hddpab * &
     &              ( dhddpab * cb * ryab/rab + 0.25_dp * dtau(13) * &
     &                ( -(dhdss*odss*oab(13)*(-ryik/rik)+ &
     &                   hdss*dodss*oab(13)*(-ryik/rik)+ &
     &                   hdss*odss*doab(13)*(-ryab/rab))*gdp*gdp &
     &                 -hdss*odss*oab(13)*2.0_dp*gdp*dgdp(2) ))

                  dc1bi(13,3) = dc1bi(13,3) + 1.0_dp/hddpab * &
     &              ( dhddpab * cb * rzab/rab + 0.25_dp * dtau(13) * &
     &                ( -(dhdss*odss*oab(13)*(-rzik/rik)+ &
     &                   hdss*dodss*oab(13)*(-rzik/rik)+ &
     &                   hdss*odss*doab(13)*(-rzab/rab))*gdp*gdp &
     &                 -hdss*odss*oab(13)*2.0_dp*gdp*dgdp(3) ))

               endif

! *** DDD
               if (hdddab /= 0.0) then
                  cb = 0.25_dp*dtau(14)/hdddab* &
     &                 (-hdss)*odss*oab(14)*gdd*gdd

                  dc1bi(14,1) = dc1bi(14,1) + 1.0_dp/hdddab * &
     &              ( dhdddab * cb * rxab/rab + 0.25_dp * dtau(14) * &
     &                ( -(dhdss*odss*oab(14)*(-rxik/rik)+ &
     &                   hdss*dodss*oab(14)*(-rxik/rik)+ &
     &                   hdss*odss*doab(14)*(-rxab/rab))*gdd*gdd &
     &                 -hdss*odss*oab(14)*2.0_dp*gdd*dgdd(1) ))

                  dc1bi(14,2) = dc1bi(14,2) + 1.0_dp/hdddab * &
     &              ( dhdddab * cb * ryab/rab + 0.25_dp * dtau(14) * &
     &                ( -(dhdss*odss*oab(14)*(-ryik/rik)+ &
     &                   hdss*dodss*oab(14)*(-ryik/rik)+ &
     &                   hdss*odss*doab(14)*(-ryab/rab))*gdd*gdd &
     &                 -hdss*odss*oab(14)*2.0_dp*gdd*dgdd(2) ))

                  dc1bi(14,3) = dc1bi(14,3) + 1.0_dp/hdddab * &
     &              ( dhdddab * cb * rzab/rab + 0.25_dp * dtau(14) * &
     &                ( -(dhdss*odss*oab(14)*(-rzik/rik)+ &
     &                   hdss*dodss*oab(14)*(-rzik/rik)+ &
     &                   hdss*odss*doab(14)*(-rzab/rab))*gdd*gdd &
     &                 -hdss*odss*oab(14)*2.0_dp*gdd*dgdd(3) ))

               endif

!     ___  End of dC1BI(I)  ___


! Calculate C1A(I-J) derivatives.

               if (rjk < rcut) then

! *** DSS
               if (hdssab /= 0.0) then
                  ca = 0.5_dp*dtau(7)/hdssab* &
     &                 (hdss*ojk(1)+odss*hsssjk)*gds

                  dc1a(7,1) = dc1a(7,1) + 1/hdssab * &
     &              ( dhdssab * ca * rxab/rab + 0.5_dp * dtau(7) * &
     &                ( (dhdss*ojk(1)+dodss*hsssjk)*(-rxik/rik)*gds &
     &                 +(hdss*ojk(1)+odss*hsssjk)*dgds(1) ))

!     &                 -(DHDSS*ODSS*OAB(7)*(-RXIK/RIK)+
!     &                   HDSS*DODSS*OAB(7)*(-RXIK/RIK)+
!     &                   HDSS*ODSS*DOAB(7)*(-RXAB/RAB))*GDS*GDS
!     &                 -HDSS*ODSS*OAB(7)*2*GDS*DGDS(1)
!     &                 -DOAB(7)*(-RXAB/RAB)*OJK(1)*HSSSJK ))

                  dc1a(7,2) = dc1a(7,2) + 1.0_dp/hdssab * &
     &              ( dhdssab * ca * ryab/rab + 0.5_dp * dtau(7) * &
     &                ( (dhdss*ojk(1)+dodss*hsssjk)*(-ryik/rik)*gds &
     &                 +(hdss*ojk(1)+odss*hsssjk)*dgds(2) ))

!     &                 -(DHDSS*ODSS*OAB(7)*(-RYIK/RIK)+
!     &                   HDSS*DODSS*OAB(7)*(-RYIK/RIK)+
!     &                   HDSS*ODSS*DOAB(7)*(-RYAB/RAB))*GDS*GDS
!     &                 -HDSS*ODSS*OAB(7)*2*GDS*DGDS(2)
!     &                 -DOAB(7)*(-RYAB/RAB)*OJK(1)*HSSSJK ))

                  dc1a(7,3) = dc1a(7,3) + 1.0_dp/hdssab * &
     &              ( dhdssab * ca * rzab/rab + 0.5_dp * dtau(7) * &
     &                ( (dhdss*ojk(1)+dodss*hsssjk)*(-rzik/rik)*gds &
     &                 +(hdss*ojk(1)+odss*hsssjk)*dgds(3) ))

!     &                 -(DHDSS*ODSS*OAB(7)*(-RZIK/RIK)+
!     &                   HDSS*DODSS*OAB(7)*(-RZIK/RIK)+
!     &                   HDSS*ODSS*DOAB(7)*(-RZAB/RAB))*GDS*GDS
!     &                 -HDSS*ODSS*OAB(7)*2*GDS*DGDS(3)
!     &                 -DOAB(7)*(-RZAB/RAB)*OJK(1)*HSSSJK ))

               endif


! *** DPS
               if (hdpsab /= 0.0) then
                  ca = 0.5_dp*dtau(9)/hdpsab* &
     &                 (hdss*ojk(3)+odss*hpssjk)*gds*gpsj

                  dc1a(9,1) = dc1a(9,1) + 1.0_dp/hdpsab * &
     &              ( dhdpsab * ca * rxab/rab + 0.5_dp * dtau(9) * &
     &                ( (dhdss*ojk(3)+dodss*hpssjk)*(-rxik/rik)*gds*gpsj &
     &                 +(hdss*ojk(3)+odss*hpssjk)* &
     &                   (dgds(1)*gpsj+gds*dgpsj(1)) ))

!     &                 -(DHDSS*ODSS*OAB(9)*(-RXIK/RIK)+
!     &                   HDSS*DODSS*OAB(9)*(-RXIK/RIK)+
!     &                   HDSS*ODSS*DOAB(9)*(-RXAB/RAB))*GDS*GDS
!     &                 -HDSS*ODSS*OAB(9)*2*GDS*DGDS(1)
!     &                 -DOAB(9)*(-RXAB/RAB)*OJK(3)*HPSSJK*GPSJ*GPSJ
!     &                 -OAB(9)*OJK(3)*HPSSJK*2*GPSJ*DGPSJ(1) ))

                  dc1a(9,2) = dc1a(9,2) + 1.0_dp/hdpsab * &
     &              ( dhdpsab * ca * ryab/rab + 0.5_dp * dtau(9) * &
     &                ( (dhdss*ojk(3)+dodss*hpssjk)*(-ryik/rik)*gds*gpsj &
     &                 +(hdss*ojk(3)+odss*hpssjk)* &
     &                   (dgds(2)*gpsj+gds*dgpsj(2)) ))

!     &                 -(DHDSS*ODSS*OAB(9)*(-RYIK/RIK)+
!     &                   HDSS*DODSS*OAB(9)*(-RYIK/RIK)+
!     &                   HDSS*ODSS*DOAB(9)*(-RYAB/RAB))*GDS*GDS
!     &                 -HDSS*ODSS*OAB(9)*2*GDS*DGDS(2)
!     &                 -DOAB(9)*(-RYAB/RAB)*OJK(3)*HPSSJK*GPSJ*GPSJ
!     &                 -OAB(9)*OJK(3)*HPSSJK*2*GPSJ*DGPSJ(2) ))

                  dc1a(9,3) = dc1a(9,3) + 1.0_dp/hdpsab * &
     &              ( dhdpsab * ca * rzab/rab + 0.5_dp * dtau(9) * &
     &                ( (dhdss*ojk(3)+dodss*hpssjk)*(-rzik/rik)*gds*gpsj &
     &                 +(hdss*ojk(3)+odss*hpssjk)* &
     &                   (dgds(3)*gpsj+gds*dgpsj(3)) ))

!     &                 -(DHDSS*ODSS*OAB(9)*(-RZIK/RIK)+
!     &                   HDSS*DODSS*OAB(9)*(-RZIK/RIK)+
!     &                   HDSS*ODSS*DOAB(9)*(-RZAB/RAB))*GDS*GDS
!     &                 -HDSS*ODSS*OAB(9)*2*GDS*DGDS(3)
!     &                 -DOAB(9)*(-RZAB/RAB)*OJK(3)*HPSSJK*GPSJ*GPSJ
!     &                 -OAB(9)*OJK(3)*HPSSJK*2*GPSJ*DGPSJ(3) ))

               endif


! *** DPP
               if (hdppab /= 0.0) then
                  ca = -0.25_dp*dtau(11)/hdppab* &
     &                 (hdss*ojk(3)+odss*hpssjk)*gdp*gppj

                  dc1a(11,1) = dc1a(11,1) + 1.0_dp/hdppab * &
     &              ( dhdppab * ca * rxab/rab - 0.25_dp * dtau(11) * &
     &                ( (dhdss*ojk(3)+dodss*hpssjk)*(-rxik/rik)*gdp*gppj &
     &                 +(hdss*ojk(3)+odss*hpssjk)* &
     &                   (dgdp(1)*gppj+gdp*dgppj(1)) ))

!     &                 -(DHDSS*ODSS*OAB(11)*(-RXIK/RIK)+
!     &                   HDSS*DODSS*OAB(11)*(-RXIK/RIK)+
!     &                   HDSS*ODSS*DOAB(11)*(-RXAB/RAB))*GDP*GDP
!     &                 -HDSS*ODSS*OAB(11)*2*GDP*DGDP(1)
!     &                 -DOAB(11)*(-RXAB/RAB)*OJK(3)*HPSSJK*GPPJ*GPPJ
!     &                 -OAB(11)*OJK(3)*HPSSJK*2*GPPJ*DGPPJ(1) ))

                  dc1a(11,2) = dc1a(11,2) + 1.0_dp/hdppab * &
     &              ( dhdppab * ca * ryab/rab - 0.25_dp * dtau(11) * &
     &                ( (dhdss*ojk(3)+dodss*hpssjk)*(-ryik/rik)*gdp*gppj &
     &                 +(hdss*ojk(3)+odss*hpssjk)* &
     &                   (dgdp(2)*gppj+gdp*dgppj(2)) ))

!     &                 -(DHDSS*ODSS*OAB(11)*(-RYIK/RIK)+
!     &                   HDSS*DODSS*OAB(11)*(-RYIK/RIK)+
!     &                   HDSS*ODSS*DOAB(11)*(-RYAB/RAB))*GDP*GDP
!     &                 -HDSS*ODSS*OAB(11)*2*GDP*DGDP(2)
!     &                 -DOAB(11)*(-RYAB/RAB)*OJK(3)*HPSSJK*GPPJ*GPPJ
!     &                 -OAB(11)*OJK(3)*HPSSJK*2*GPPJ*DGPPJ(2) ))

                  dc1a(11,3) = dc1a(11,3) + 1.0_dp/hdppab * &
     &              ( dhdppab * ca * rzab/rab - 0.25_dp * dtau(11) * &
     &                ( (dhdss*ojk(3)+dodss*hpssjk)*(-rzik/rik)*gdp*gppj &
     &                 +(hdss*ojk(3)+odss*hpssjk)* &
     &                   (dgdp(3)*gppj+gdp*dgppj(3)) ))

!     &                 -(DHDSS*ODSS*OAB(11)*(-RZIK/RIK)+
!     &                   HDSS*DODSS*OAB(11)*(-RZIK/RIK)+
!     &                   HDSS*ODSS*DOAB(11)*(-RZAB/RAB))*GDP*GDP
!     &                 -HDSS*ODSS*OAB(11)*2*GDP*DGDP(3)
!     &                 -DOAB(11)*(-RZAB/RAB)*OJK(3)*HPSSJK*GPPJ*GPPJ
!     &                 -OAB(11)*OJK(3)*HPSSJK*2*GPPJ*DGPPJ(3) ))

               endif


! *** DDS
               if (hddsab /= 0.0) then
                  ca = 0.5_dp*dtau(12)/hddsab* &
     &                 (hdss*ojk(7)+odss*hdssjk)*gds*gdsj

                  dc1a(12,1) = dc1a(12,1) + 1.0_dp/hddsab * &
     &              ( dhddsab * ca * rxab/rab + 0.5_dp * dtau(12) * &
     &                ( (dhdss*ojk(7)+dodss*hdssjk)*(-rxik/rik)*gds*gdsj &
     &                 +(hdss*ojk(7)+odss*hdssjk)* &
     &                   (dgds(1)*gdsj+gds*dgdsj(1)) ))

!     &                 -(DHDSS*ODSS*OAB(12)*(-RXIK/RIK)+
!     &                   HDSS*DODSS*OAB(12)*(-RXIK/RIK)+
!     &                   HDSS*ODSS*DOAB(12)*(-RXAB/RAB))*GDS*GDS
!     &                 -HDSS*ODSS*OAB(12)*2*GDS*DGDS(1)
!     &                 -DOAB(12)*(-RXAB/RAB)*OJK(7)*HDSSJK*GDSJ*GDSJ
!     &                 -OAB(12)*OJK(7)*HDSSJK*2*GDSJ*DGDSJ(1) ))

                  dc1a(12,2) = dc1a(12,2) + 1.0_dp/hddsab * &
     &              ( dhddsab * ca * ryab/rab + 0.5_dp * dtau(12) * &
     &                ( (dhdss*ojk(7)+dodss*hdssjk)*(-ryik/rik)*gds*gdsj &
     &                 +(hdss*ojk(7)+odss*hdssjk)* &
     &                   (dgds(2)*gdsj+gds*dgdsj(2)) ))

!     &                 -(DHDSS*ODSS*OAB(12)*(-RYIK/RIK)+
!     &                   HDSS*DODSS*OAB(12)*(-RYIK/RIK)+
!     &                   HDSS*ODSS*DOAB(12)*(-RYAB/RAB))*GDS*GDS
!     &                 -HDSS*ODSS*OAB(12)*2*GDS*DGDS(2)
!     &                 -DOAB(12)*(-RYAB/RAB)*OJK(7)*HDSSJK*GDSJ*GDSJ
!     &                 -OAB(12)*OJK(7)*HDSSJK*2*GDSJ*DGDSJ(2) ))

                  dc1a(12,3) = dc1a(12,3) + 1.0_dp/hddsab * &
     &              ( dhddsab * ca * rzab/rab + 0.5_dp * dtau(12) * &
     &                ( (dhdss*ojk(7)+dodss*hdssjk)*(-rzik/rik)*gds*gdsj &
     &                 +(hdss*ojk(7)+odss*hdssjk)* &
     &                   (dgds(3)*gdsj+gds*dgdsj(3)) ))

!     &                 -(DHDSS*ODSS*OAB(12)*(-RZIK/RIK)+
!     &                   HDSS*DODSS*OAB(12)*(-RZIK/RIK)+
!     &                   HDSS*ODSS*DOAB(12)*(-RZAB/RAB))*GDS*GDS
!     &                 -HDSS*ODSS*OAB(12)*2*GDS*DGDS(3)
!     &                 -DOAB(12)*(-RZAB/RAB)*OJK(7)*HDSSJK*GDSJ*GDSJ
!     &                 -OAB(12)*OJK(7)*HDSSJK*2*GDSJ*DGDSJ(3) ))

               endif


! *** DDP
               if (hddpab /= 0.0) then
                  ca = -0.25_dp*dtau(13)/hddpab* &
     &                 (hdss*ojk(7)+odss*hdssjk)*gdp*gdpj

                  dc1a(13,1) = dc1a(13,1) + 1.0_dp/hddpab * &
     &              ( dhddpab * ca * rxab/rab + 0.25_dp * dtau(13) * &
     &                ( (dhdss*ojk(7)+dodss*hdssjk)*(rxik/rik)*gdp*gdpj &
     &                 +(hdss*ojk(7)+odss*hdssjk)* &
     &                   (-dgdp(1)*gdpj-gdp*dgdpj(1)) ))

!     &                 -(DHDSS*ODSS*OAB(13)*(-RXIK/RIK)+
!     &                   HDSS*DODSS*OAB(13)*(-RXIK/RIK)+
!     &                   HDSS*ODSS*DOAB(13)*(-RXAB/RAB))*GDP*GDP
!     &                 -HDSS*ODSS*OAB(13)*2*GDP*DGDP(1)
!     &                 -DOAB(13)*(-RXAB/RAB)*OJK(7)*HDSSJK*GDPJ*GDPJ
!     &                 -OAB(13)*OJK(7)*HDSSJK*2*GDPJ*DGDPJ(1) ))

                  dc1a(13,2) = dc1a(13,2) + 1.0_dp/hddpab * &
     &              ( dhddpab * ca * ryab/rab + 0.25_dp * dtau(13) * &
     &                ( (dhdss*ojk(7)+dodss*hdssjk)*(ryik/rik)*gdp*gdpj &
     &                 +(hdss*ojk(7)+odss*hdssjk)* &
     &                   (-dgdp(2)*gdpj-gdp*dgdpj(2)) ))

!     &                 -(DHDSS*ODSS*OAB(13)*(-RYIK/RIK)+
!     &                   HDSS*DODSS*OAB(13)*(-RYIK/RIK)+
!     &                   HDSS*ODSS*DOAB(13)*(-RYAB/RAB))*GDP*GDP
!     &                 -HDSS*ODSS*OAB(13)*2*GDP*DGDP(2)
!     &                 -DOAB(13)*(-RYAB/RAB)*OJK(7)*HDSSJK*GDPJ*GDPJ
!     &                 -OAB(13)*OJK(7)*HDSSJK*2*GDPJ*DGDPJ(2) ))

                  dc1a(13,3) = dc1a(13,3) + 1.0_dp/hddpab * &
     &              ( dhddpab * ca * rzab/rab + 0.25_dp * dtau(13) * &
     &                ( (dhdss*ojk(7)+dodss*hdssjk)*(rzik/rik)*gdp*gdpj &
     &                 +(hdss*ojk(7)+odss*hdssjk)* &
     &                   (-dgdp(3)*gdpj-gdp*dgdpj(3)) ))

!     &                 -(DHDSS*ODSS*OAB(13)*(-RZIK/RIK)+
!     &                   HDSS*DODSS*OAB(13)*(-RZIK/RIK)+
!     &                   HDSS*ODSS*DOAB(13)*(-RZAB/RAB))*GDP*GDP
!     &                 -HDSS*ODSS*OAB(13)*2*GDP*DGDP(3)
!     &                 -DOAB(13)*(-RZAB/RAB)*OJK(7)*HDSSJK*GDPJ*GDPJ
!     &                 -OAB(13)*OJK(7)*HDSSJK*2*GDPJ*DGDPJ(3) ))

               endif

! *** DDD
               if (hdddab /= 0.0) then
                  ca = 0.25_dp*dtau(14)/hdddab* &
     &                    (hdss*ojk(7)+odss*hdssjk)*gdd*gddj

                  dc1a(14,1) = dc1a(14,1) + 1/hdddab * &
     &              ( dhdddab * ca * rxab/rab + 0.25_dp * dtau(14) * &
     &                ( (dhdss*ojk(7)+dodss*hdssjk)*(-rxik/rik)*gdd*gddj &
     &                 +(hdss*ojk(7)+odss*hdssjk)* &
     &                   (dgdd(1)*gddj+gdd*dgddj(1)) ))

!     &                 -(DHDSS*ODSS*OAB(14)*(-RXIK/RIK)+
!     &                   HDSS*DODSS*OAB(14)*(-RXIK/RIK)+
!     &                   HDSS*ODSS*DOAB(14)*(-RXAB/RAB))*GDD*GDD
!     &                 -HDSS*ODSS*OAB(14)*2*GDD*DGDD(1)
!     &                 -DOAB(14)*(-RXAB/RAB)*OJK(7)*HDSSJK*GDDJ*GDDJ
!     &                 -OAB(14)*OJK(7)*HDSSJK*2*GDDJ*DGDDJ(1) ))

                  dc1a(14,2) = dc1a(14,2) + 1.0_dp/hdddab * &
     &              ( dhdddab * ca * ryab/rab + 0.25_dp * dtau(14) * &
     &                ( (dhdss*ojk(7)+dodss*hdssjk)*(-ryik/rik)*gdd*gddj &
     &                 +(hdss*ojk(7)+odss*hdssjk)* &
     &                   (dgdd(2)*gddj+gdd*dgddj(2)) ))

!     &                 -(DHDSS*ODSS*OAB(14)*(-RYIK/RIK)+
!     &                   HDSS*DODSS*OAB(14)*(-RYIK/RIK)+
!     &                   HDSS*ODSS*DOAB(14)*(-RYAB/RAB))*GDD*GDD
!     &                 -HDSS*ODSS*OAB(14)*2*GDD*DGDD(2)
!     &                 -DOAB(14)*(-RYAB/RAB)*OJK(7)*HDSSJK*GDDJ*GDDJ
!     &                 -OAB(14)*OJK(7)*HDSSJK*2*GDDJ*DGDDJ(2) ))

                  dc1a(14,3) = dc1a(14,3) + 1.0_dp/hdddab * &
     &              ( dhdddab * ca * rzab/rab + 0.25_dp * dtau(14) * &
     &                ( (dhdss*ojk(7)+dodss*hdssjk)*(-rzik/rik)*gdd*gddj &
     &                 +(hdss*ojk(7)+odss*hdssjk)* &
     &                   (dgdd(3)*gddj+gdd*dgddj(3)) ))

!     &                 -(DHDSS*ODSS*OAB(14)*(-RZIK/RIK)+
!     &                   HDSS*DODSS*OAB(14)*(-RZIK/RIK)+
!     &                   HDSS*ODSS*DOAB(14)*(-RZAB/RAB))*GDD*GDD
!     &                 -HDSS*ODSS*OAB(14)*2*GDD*DGDD(3)
!     &                 -DOAB(14)*(-RZAB/RAB)*OJK(7)*HDSSJK*GDDJ*GDDJ
!     &                 -OAB(14)*OJK(7)*HDSSJK*2*GDDJ*DGDDJ(3) ))

               endif

            endif

            endif

         endif
         ja = ja + 1
      enddo     



!
!       Calculate the second moment derivatives on atom IB (also J)
!

      ja = aptr( ib )

      do while ( bptr( ja )  /=  eol )
         k = bptr( ja )

         if ((k /= ia).and.(k /= ib)) then
            rxjk = ad( 1, k ) - ad( 1, ib )
            ryjk = ad( 2, k ) - ad( 2, ib )
            rzjk = ad( 3, k ) - ad( 3, ib )               
            rjk = sqrt(rxjk*rxjk + ryjk*ryjk + rzjk*rzjk)
            num = rxjk*rxab + ryjk*ryab + rzjk*rzab
            den = rab*rjk
            nod = -num/den

            if (nod > 1.0_dp-1.0e-6_dp) then
               thetaj = 0.0_dp
            elseif (nod < -1.0_dp+1.0e-6_dp) then
               thetaj = pi
            else
               thetaj = acos(nod)
            endif

            bt = btype(z(ib),z(k))

            call radf(rjk,hsss,hsps,hpss,hpps, & 
     &                hppp,hsds,hdss,hpds,hdps, & 
     &                hpdp,hdpp,hdds,hddp,hddd, & 
     &                bnddat(1,bt),bndscl(1,1,bt))

            call dradf(rjk,dhsss,dhsps,dhpss,dhpps, & 
     &                 dhppp,dhsds,dhdss,dhpds,dhdps, & 
     &                 dhpdp,dhdpp,dhdds,dhddp,dhddd, & 
     &                 bnddat(1,bt),bndscl(1,1,bt))


            cos_j = cos(thetaj)
            sin_j = sin(thetaj)
            cos2_j = cos(2.0*thetaj)
            sin2_j = sin(2.0*thetaj)

            gssj = 1.0_dp
            gpsj = cos_j
            gppj = sin_j
            gdsj = 0.25_dp*(1.0_dp+3.0_dp*cos2_j)
            gdpj = sqrt3/2.0_dp*sin2_j
            gddj = sqrt3/4.0_dp*(1.0_dp-cos2_j)

            vecpr_j = sqrt((ryab*rzjk-ryjk*rzab)*(ryab*rzjk-ryjk*rzab) & 
     &           + (rzab*rxjk-rzjk*rxab)*(rzab*rxjk-rzjk*rxab) & 
     &           + (rxab*ryjk-rxjk*ryab)*(rxab*ryjk-rxjk*ryab))

!            VECPR_J = SQRT( (RYAB*RZJK-RYJK*RZAB)**2 +
!     &                      (RZAB*RXJK-RZJK*RXAB)**2 +
!     &                      (RXAB*RYJK-RXJK*RYAB)**2 )

! Derivatives of Cos(ThetaJ) with respect to x(=1), y(=2), z(=3)

            dcos_j(1) = 1.0_dp/rab * (rxjk/rjk+rxab/rab*cos_j)
            dcos_j(2) = 1.0_dp/rab * (ryjk/rjk+ryab/rab*cos_j)
            dcos_j(3) = 1.0_dp/rab * (rzjk/rjk+rzab/rab*cos_j)

          if ((nod > -1.0_dp+1.0e-6_dp).and.(nod < 1.0_dp-1.0e-6_dp)) then
               dsin_j(1) = -cos_j/sin_j*dcos_j(1)
               dsin_j(2) = -cos_j/sin_j*dcos_j(2)
               dsin_j(3) = -cos_j/sin_j*dcos_j(3)
          endif

! Derivatives of Gll'J with respect to x(=1), y(=2), z(=3)
!           DGSSJ(1-3) = 0.D0

            dgssj(1) = 0.0_dp
            dgssj(2) = 0.0_dp
            dgssj(3) = 0.0_dp

            dgpsj(1) = dcos_j(1)
            dgpsj(2) = dcos_j(2)
            dgpsj(3) = dcos_j(3)

            dgppj(1) = dsin_j(1)
            dgppj(2) = dsin_j(2)
            dgppj(3) = dsin_j(3)

            dgdsj(1) = 3.0_dp*cos_j * dcos_j(1)
            dgdsj(2) = 3.0_dp*cos_j * dcos_j(2)
            dgdsj(3) = 3.0_dp*cos_j * dcos_j(3)

            dgdpj(1) = sqrt3*(dsin_j(1)*cos_j+sin_j*dcos_j(1))
            dgdpj(2) = sqrt3*(dsin_j(2)*cos_j+sin_j*dcos_j(2))
            dgdpj(3) = sqrt3*(dsin_j(3)*cos_j+sin_j*dcos_j(3))

            dgddj(1) = -sqrt3*cos_j * dcos_j(1)
            dgddj(2) = -sqrt3*cos_j * dcos_j(2)
            dgddj(3) = -sqrt3*cos_j * dcos_j(3)


! Check whether we have s electrons

            if (bndscl(1,1,bt) == 0.0) then
               osss = 0.0_dp
               dosss = 0.0_dp
            else
               osss = vscale(rjk,bndscl(1,1,bt),bndscl(2,1,bt), & 
     &                      bndscl(5,1,bt),bndscl(6,1,bt), & 
     &                      bndscl(3,1,bt),bndscl(4,1,bt), & 
     &                      bndscl(7,1,bt),bndscl(8,1,bt), & 
     &                      bndscl(9,1,bt),bndscl(10,1,bt), & 
     &                      bndscl(11,1,bt),bndscl(12,1,bt))

               osss = osss*bndscl(13,1,bt)

               dosss = dscale(rjk,bndscl(1,1,bt),bndscl(2,1,bt), & 
     &                       bndscl(5,1,bt),bndscl(6,1,bt), & 
     &                       bndscl(3,1,bt),bndscl(4,1,bt), & 
     &                       bndscl(7,1,bt),bndscl(8,1,bt), & 
     &                       bndscl(9,1,bt),bndscl(10,1,bt), & 
     &                       bndscl(11,1,bt),bndscl(12,1,bt))

               dosss = dosss*bndscl(13,1,bt)

               mu2b(1) = mu2b(1) + osss*osss*gssj*gssj
               mu2b(3) = mu2b(1)
               mu2b(7) = mu2b(1)

! * SSS        MU2B(1) = MU2B(1) + OSSS*OSSS*GSSJ*GSSJ

               dmu2b(1,1) = dmu2b(1,1) + &
     &                      2.0_dp*osss*osss*gssj*dgssj(1)
               dmu2b(1,2) = dmu2b(1,2) + &
     &                      2.0_dp*osss*osss*gssj*dgssj(2)
               dmu2b(1,3) = dmu2b(1,3) + &
     &                      2.0_dp*osss*osss*gssj*dgssj(3)


! * SPS(PSS)   MU2B(3) = MU2B(3) + OSSS*OSSS*GSSJ*GSSJ

               dmu2b(3,1) = dmu2b(1,1)
               dmu2b(3,2) = dmu2b(1,2)
               dmu2b(3,3) = dmu2b(1,3)


! * SDD(DSS)   MU2B(7) = MU2B(7) + OSSS*OSSS*GSSJ*GSSJ

               dmu2b(7,1) = dmu2b(1,1)
               dmu2b(7,2) = dmu2b(1,2)
               dmu2b(7,3) = dmu2b(1,3)

!    ___  End of derivatives of MU2B  ___


! Derivatives of C1BJ(J)

! *** SSS

               if (hsssab /= 0.0) then
                  cb = -0.5_dp*dtau(1)/hsssab* &
     &              hsss*osss*oab(1)

                  dc1bj(1,1) = dc1bj(1,1) + 1.0_dp/hsssab * &
     &              ( dhsssab * cb * rxab/rab + 0.5_dp * dtau(1) * &
     &                 (doab(1)*(rxab/rab)*osss*hsss) )

                  dc1bj(1,2) = dc1bj(1,2) + 1.0_dp/hsssab * &
     &              ( dhsssab * cb * ryab/rab + 0.5_dp * dtau(1) * &
     &                 (doab(1)*(ryab/rab)*osss*hsss) )

                  dc1bj(1,3) = dc1bj(1,3) + 1.0_dp/hsssab * &
     &              ( dhsssab * cb * rzab/rab + 0.5_dp * dtau(1) * &
     &                 (doab(1)*(rzab/rab)*osss*hsss) )

               endif


! *** SPS (=PSS)

               if (hpssab /= 0.0) then
                  cb = -0.5_dp*dtau(3)/hpssab* &
     &              hsss*osss*oab(3)

                  dc1bj(3,1) = dc1bj(3,1) + 1.0_dp/hpssab * &
     &              ( dhpssab * cb * rxab/rab + 0.5_dp * dtau(3) * &
     &                 doab(3)*(rxab/rab)*osss*hsss )

                  dc1bj(3,2) = dc1bj(3,2) + 1.0_dp/hpssab * &
     &              ( dhpssab * cb * ryab/rab + 0.5_dp * dtau(3) * &
     &                 doab(3)*(ryab/rab)*osss*hsss )

                  dc1bj(3,3) = dc1bj(3,3) + 1.0_dp/hpssab * &
     &              ( dhpssab * cb * rzab/rab + 0.5_dp * dtau(3) * &
     &                 doab(3)*(rzab/rab)*osss*hsss )

               endif


! *** SDS (=DSS)

               if (hdssab /= 0.0) then
                  cb = -0.5_dp*dtau(7)/hdssab* &
     &              hsss*osss*oab(7)

                  dc1bj(7,1) = dc1bj(7,1) + 1.0_dp/hdssab * &
     &              ( dhdssab * cb * rxab/rab + 0.5_dp * dtau(7) * &
     &                 doab(7)*(rxab/rab)*osss*hsss )

                  dc1bj(7,2) = dc1bj(7,2) + 1.0_dp/hdssab * &
     &              ( dhdssab * cb * ryab/rab + 0.5_dp * dtau(7) * &
     &                 doab(7)*(ryab/rab)*osss*hsss )

                  dc1bj(7,3) = dc1bj(7,3) + 1.0_dp/hdssab * &
     &              ( dhdssab * cb * rzab/rab + 0.5_dp * dtau(7) * &
     &                 doab(7)*(rzab/rab)*osss*hsss )

               endif

            endif


! Check whether we have p electrons

            if (bndscl(1,3,bt) == 0.0) then
               opss = 0.0_dp
               dopss = 0.0_dp
            else
               opss = vscale(rjk,bndscl(1,3,bt),bndscl(2,3,bt), & 
     &                      bndscl(5,3,bt),bndscl(6,3,bt), & 
     &                      bndscl(3,3,bt),bndscl(4,3,bt), & 
     &                      bndscl(7,3,bt),bndscl(8,3,bt), & 
     &                      bndscl(9,3,bt),bndscl(10,3,bt), & 
     &                      bndscl(11,3,bt),bndscl(12,3,bt))

               opss = opss*bndscl(13,3,bt)

               dopss = dscale(rjk,bndscl(1,3,bt),bndscl(2,3,bt), & 
     &                       bndscl(5,3,bt),bndscl(6,3,bt), & 
     &                       bndscl(3,3,bt),bndscl(4,3,bt), & 
     &                       bndscl(7,3,bt),bndscl(8,3,bt), & 
     &                       bndscl(9,3,bt),bndscl(10,3,bt), & 
     &                       bndscl(11,3,bt),bndscl(12,3,bt))

               dopss = dopss*bndscl(13,3,bt)


               mu2b(2) = mu2b(2) + opss*opss*gpsj*gpsj
               mu2b(4) = mu2b(2)
               mu2b(5) = mu2b(5) + 0.5_dp*opss*opss*gppj*gppj
               mu2b(9) = mu2b(2)
               mu2b(11) = mu2b(5)


! * PSS(SPS)   MU2B(2) = MU2B(2) + OPSS*OPSS*GPSJ*GPSJ

               dmu2b(2,1) = dmu2b(2,1) + &
     &                      2.0_dp*opss*opss*gpsj*dgpsj(1)
               dmu2b(2,2) = dmu2b(2,2) + &
     &                      2.0_dp*opss*opss*gpsj*dgpsj(2)
               dmu2b(2,3) = dmu2b(2,3) + &
     &                      2.0_dp*opss*opss*gpsj*dgpsj(3)


! * PPS        MU2B(4) = MU2B(4) + OPSS*OPSS*GPSJ*GPSJ

               dmu2b(4,1) = dmu2b(2,1)
               dmu2b(4,2) = dmu2b(2,2)
               dmu2b(4,3) = dmu2b(2,3)


! * PPP        MU2B(5) = MU2B(5) + 0.5*OPSS*OPSS*GPPJ*GPPJ

               dmu2b(5,1) = dmu2b(5,1) + &
     &                      opss*opss*gppj*dgppj(1)
               dmu2b(5,2) = dmu2b(5,2) + &
     &                      opss*opss*gppj*dgppj(2)
               dmu2b(5,3) = dmu2b(5,3) + &
     &                      opss*opss*gppj*dgppj(3)


! * PDS(DPS)   MU2B(9) = MU2B(9) + OPSS*OPSS*GPSJ*GPSJ

               dmu2b(9,1) = dmu2b(2,1)
               dmu2b(9,2) = dmu2b(2,2)
               dmu2b(9,3) = dmu2b(2,3)


! * PDP(DPP)   MU2B(11) = MU2B(11) + 0.5*ODSS*ODSS*GPPJ*GPPJ

               dmu2b(11,1) = dmu2b(5,1)
               dmu2b(11,2) = dmu2b(5,2)
               dmu2b(11,3) = dmu2b(5,3)

!    ___  End of derivatives of MU2B  ___


! Derivatives of C1BJ(J)

! *** PSS (=SPS)
               if (hspsab /= 0.0) then
                  cb = -0.5_dp*dtau(2)/hspsab* &
     &              hpss*opss*oab(2)*gpsj*gpsj

                  dc1bj(2,1) = dc1bj(2,1) + 1.0_dp/hspsab * &
     &              ( dhspsab * cb * rxab/rab + 0.5_dp * dtau(2) * &
     &               ( -doab(2)*(-rxab/rab)*opss*hpss*gpsj*gpsj &
     &                  -oab(2)*opss*hpss*2.0_dp*gpsj*dgpsj(1) ))

                  dc1bj(2,2) = dc1bj(2,2) + 1.0_dp/hspsab * &
     &              ( dhspsab * cb * ryab/rab + 0.5_dp * dtau(2) * &
     &               ( -doab(2)*(-ryab/rab)*opss*hpss*gpsj*gpsj &
     &                 -oab(2)*opss*hpss*2.0_dp*gpsj*dgpsj(2) ))

                  dc1bj(2,3) = dc1bj(2,3) + 1.0_dp/hspsab * &
     &              ( dhspsab * cb * rzab/rab + 0.5_dp * dtau(2) * &
     &               ( -doab(2)*(-rzab/rab)*opss*hpss*gpsj*gpsj &
     &                 -oab(2)*opss*hpss*2.0_dp*gpsj*dgpsj(3) ))

               endif


! *** PPS (=PPS)
               if (hppsab /= 0.0) then
                  cb = -0.5_dp*dtau(4)/hppsab* &
     &              hpss*opss*oab(4)*gpsj*gpsj

                  dc1bj(4,1) = dc1bj(4,1) + 1.0_dp/hppsab * &
     &              ( dhppsab * cb * rxab/rab + 0.5_dp * dtau(4) * &
     &               ( -doab(4)*(-rxab/rab)*opss*hpss*gpsj*gpsj &
     &                 -oab(4)*opss*hpss*2.0_dp*gpsj*dgpsj(1) ))

                  dc1bj(4,2) = dc1bj(4,2) + 1.0_dp/hppsab * &
     &              ( dhppsab * cb * ryab/rab + 0.5_dp * dtau(4) * &
     &               ( -doab(4)*(-ryab/rab)*opss*hpss*gpsj*gpsj &
     &                 -oab(4)*opss*hpss*2.0_dp*gpsj*dgpsj(2) ))

                  dc1bj(4,3) = dc1bj(4,3) + 1.0_dp/hppsab * &
     &              ( dhppsab * cb * rzab/rab + 0.5_dp * dtau(4) * &
     &               ( -doab(4)*(-rzab/rab)*opss*hpss*gpsj*gpsj &
     &                 -oab(4)*opss*hpss*2.0_dp*gpsj*dgpsj(3) ))

               endif


! *** PPP (=PPP)
               if (hpppab /= 0.0) then
                  cb = 0.25_dp*dtau(5)/hpppab* &
     &              hpss*opss*oab(5)*gppj*gppj

                  dc1bj(5,1) = dc1bj(5,1) + 1.0_dp/hpppab * &
     &              ( dhpppab * cb * rxab/rab - 0.25_dp * dtau(5) * &
     &               ( -doab(5)*(-rxab/rab)*opss*hpss*gppj*gppj &
     &                 -oab(5)*opss*hpss*2.0_dp*gppj*dgppj(1) ))

                  dc1bj(5,2) = dc1bj(5,2) + 1.0_dp/hpppab * &
     &              ( dhpppab * cb * ryab/rab - 0.25_dp * dtau(5) * &
     &               ( -doab(5)*(-ryab/rab)*opss*hpss*gppj*gppj &
     &                 -oab(5)*opss*hpss*2.0_dp*gppj*dgppj(2) ))

                  dc1bj(5,3) = dc1bj(5,3) + 1.0_dp/hpppab * &
     &              ( dhpppab * cb * rzab/rab - 0.25_dp * dtau(5) * &
     &               ( -doab(5)*(-rzab/rab)*opss*hpss*gppj*gppj &
     &                 -oab(5)*opss*hpss*2.0_dp*gppj*dgppj(3) ))

               endif


! *** PDS (=DPS)
               if (hdpsab /= 0.0) then
                  cb = -0.5_dp*dtau(9)/hdpsab* &
     &              hpss*opss*oab(9)*gpsj*gpsj

                  dc1bj(9,1) = dc1bj(9,1) + 1.0_dp/hdpsab * &
     &              ( dhdpsab * cb * rxab/rab + 0.5_dp * dtau(9) * &
     &               ( -doab(9)*(-rxab/rab)*opss*hpss*gpsj*gpsj &
     &                 -oab(9)*opss*hpss*2.0_dp*gpsj*dgpsj(1) ))

                  dc1bj(9,2) = dc1bj(9,2) + 1.0_dp/hdpsab * &
     &              ( dhdpsab * cb * ryab/rab + 0.5_dp * dtau(9) * &
     &               ( -doab(9)*(-ryab/rab)*opss*hpss*gpsj*gpsj &
     &                 -oab(9)*opss*hpss*2.0_dp*gpsj*dgpsj(2) ))

                  dc1bj(9,3) = dc1bj(9,3) + 1.0_dp/hdpsab * &
     &              ( dhdpsab * cb * rzab/rab + 0.5_dp * dtau(9) * &
     &               ( -doab(9)*(-rzab/rab)*opss*hpss*gpsj*gpsj &
     &                 -oab(9)*opss*hpss*2.0_dp*gpsj*dgpsj(3) ))

               endif


! *** PDP (=DPP)
               if (hdppab /= 0.0) then
                  cb = 0.25_dp*dtau(11)/hdppab* &
     &              hpss*opss*oab(11)*gppj*gppj

                  dc1bj(11,1) = dc1bj(11,1) + 1.0_dp/hdppab * &
     &              ( dhdppab * cb * rxab/rab - 0.25_dp * dtau(11) * &
     &               ( -doab(11)*(-rxab/rab)*opss*hpss*gppj*gppj &
     &                 -oab(11)*opss*hpss*2.0_dp*gppj*dgppj(1) ))

                  dc1bj(11,2) = dc1bj(11,2) + 1.0_dp/hdppab * &
     &              ( dhdppab * cb * ryab/rab - 0.25_dp * dtau(11) * &
     &               ( -doab(11)*(-ryab/rab)*opss*hpss*gppj*gppj &
     &                 -oab(11)*opss*hpss*2.0_dp*gppj*dgppj(2) ))

                  dc1bj(11,3) = dc1bj(11,3) + 1.0_dp/hdppab * &
     &              ( dhdppab * cb * rzab/rab - 0.25_dp * dtau(11) * &
     &               ( -doab(11)*(-rzab/rab)*opss*hpss*gppj*gppj &
     &                 -oab(11)*opss*hpss*2.0_dp*gppj*dgppj(3) ))

               endif


            endif


! Check whether we have d electrons


            if (bndscl(1,7,bt) == 0.0) then
               odss = 0.0_dp
               dodss = 0.0_dp
            else
               odss = vscale(rjk,bndscl(1,7,bt),bndscl(2,7,bt), & 
     &                      bndscl(5,7,bt),bndscl(6,7,bt), & 
     &                      bndscl(3,7,bt),bndscl(4,7,bt), & 
     &                      bndscl(7,7,bt),bndscl(8,7,bt), & 
     &                      bndscl(9,7,bt),bndscl(10,7,bt), & 
     &                      bndscl(11,7,bt),bndscl(12,7,bt))

               odss = odss*bndscl(13,7,bt)

               dodss = dscale(rjk,bndscl(1,7,bt),bndscl(2,7,bt), & 
     &                       bndscl(5,7,bt),bndscl(6,7,bt), & 
     &                       bndscl(3,7,bt),bndscl(4,7,bt), & 
     &                       bndscl(7,7,bt),bndscl(8,7,bt), & 
     &                       bndscl(9,7,bt),bndscl(10,7,bt), & 
     &                       bndscl(11,7,bt),bndscl(12,7,bt))

               dodss = dodss*bndscl(13,7,bt)


               mu2b(6) = mu2b(6) + odss*odss*gdsj*gdsj
               mu2b(8) = mu2b(6)
               mu2b(10) = mu2b(10) + 0.5_dp*odss*odss*gdpj*gdpj
               mu2b(12) = mu2b(6)
               mu2b(13) = mu2b(10)
               mu2b(14) = mu2b(14) + 0.5_dp*odss*odss*gddj*gddj


! * DSS        MU2B(6) = MU2B(6) + ODSS*ODSS*GDSJ*GDSJ

               dmu2b(6,1) = dmu2b(6,1) + &
     &                      2.0_dp*odss*odss*gdsj*dgdsj(1)
               dmu2b(6,2) = dmu2b(6,2) + &
     &                      2.0_dp*odss*odss*gdsj*dgdsj(2)
               dmu2b(6,3) = dmu2b(6,3) + &
     &                      2.0_dp*odss*odss*gdsj*dgdsj(3)


! * DPS(PDS)   MU2B(8) = MU2B(8) + ODSS*ODSS*GDSJ*GDSJ

               dmu2b(8,1) = dmu2b(6,1)
               dmu2b(8,2) = dmu2b(6,2)
               dmu2b(8,3) = dmu2b(6,3)


! * DPP(PDP)   MU2B(10) = MU2B(10) + 0.5*ODSS*ODSS*GDPJ*GDPJ

               dmu2b(10,1) = dmu2b(10,1) + &
     &                      odss*odss*gdpj*dgdpj(1)
               dmu2b(10,2) = dmu2b(10,2) + &
     &                      odss*odss*gdpj*dgdpj(2)
               dmu2b(10,3) = dmu2b(10,3) + &
     &                      odss*odss*gdpj*dgdpj(3)


! * DDS        MU2B(12) = MU2B(12) + ODSS*ODSS*GDSJ*GDSJ

               dmu2b(12,1) = dmu2b(6,1)
               dmu2b(12,2) = dmu2b(6,2)
               dmu2b(12,3) = dmu2b(6,3)


! * DDP        MU2B(13) = MU2B(13) + 0.5*ODSS*ODSS*GDPJ*GDPJ

               dmu2b(13,1) = dmu2b(10,1)
               dmu2b(13,2) = dmu2b(10,2)
               dmu2b(13,3) = dmu2b(10,3)


! * DDD        MU2B(14) = MU2B(14) + 0.5*ODSS*ODSS*GDDJ*GDDJ

               dmu2b(14,1) = dmu2b(14,1) + &
     &                      odss*odss*gddj*dgddj(1)
               dmu2b(14,2) = dmu2b(14,2) + &
     &                      odss*odss*gddj*dgddj(2)
               dmu2b(14,3) = dmu2b(14,3) + &
     &                      odss*odss*gddj*dgddj(3)

!    ___  End of derivatives of MU2B  ___


! Derivatives of C1BJ(J)

! *** DSS (=SDS)
               if (hsdsab /= 0.0) then
                  cb = -0.5_dp*dtau(6)/hsdsab* &
     &              hdss*odss*oab(6)*gdsj*gdsj

                  dc1bj(6,1) = dc1bj(6,1) + 1.0_dp/hsdsab * &
     &              ( dhsdsab * cb * rxab/rab + 0.5_dp * dtau(6) * &
     &               ( -doab(6)*(-rxab/rab)*odss*hdss*gdsj*gdsj &
     &                 -oab(6)*odss*hdss*2.0_dp*gdsj*dgdsj(1) ))

                  dc1bj(6,2) = dc1bj(6,2) + 1.0_dp/hsdsab * &
     &              ( dhsdsab * cb * ryab/rab + 0.5_dp * dtau(6) * &
     &               ( -doab(6)*(-ryab/rab)*odss*hdss*gdsj*gdsj &
     &                 -oab(6)*odss*hdss*2.0_dp*gdsj*dgdsj(2) ))

                  dc1bj(6,3) = dc1bj(6,3) + 1.0_dp/hsdsab * &
     &              ( dhsdsab * cb * rzab/rab + 0.5_dp * dtau(6) * &
     &               ( -doab(6)*(-rzab/rab)*odss*hdss*gdsj*gdsj &
     &                 -oab(6)*odss*hdss*2.0_dp*gdsj*dgdsj(3) ))

               endif


! *** DPS (=PDS)
               if (hpdsab /= 0.0) then
                  cb = -0.5_dp*dtau(8)/hpdsab* &
     &                 hdss*odss*oab(8)*gdsj*gdsj

                  dc1bj(8,1) = dc1bj(8,1) + 1.0_dp/hpdsab * &
     &              ( dhpdsab * cb * rxab/rab + 0.5_dp * dtau(8) * &
     &               ( -doab(8)*(-rxab/rab)*odss*hdss*gdsj*gdsj &
     &                 -oab(8)*odss*hdss*2.0_dp*gdsj*dgdsj(1) ))

                  dc1bj(8,2) = dc1bj(8,2) + 1.0_dp/hpdsab * &
     &              ( dhpdsab * cb * ryab/rab + 0.5_dp * dtau(8) * &
     &               ( -doab(8)*(-ryab/rab)*odss*hdss*gdsj*gdsj &
     &                 -oab(8)*odss*hdss*2.0_dp*gdsj*dgdsj(2) ))

                  dc1bj(8,3) = dc1bj(8,3) + 1.0_dp/hpdsab * &
     &              ( dhpdsab * cb * rzab/rab + 0.5_dp * dtau(8) * &
     &               ( -doab(8)*(-rzab/rab)*odss*hdss*gdsj*gdsj &
     &                 -oab(8)*odss*hdss*2.0_dp*gdsj*dgdsj(3) ))

               endif


! *** DPP (=PDP)
               if (hpdpab /= 0.0) then
                  cb = 0.25_dp*dtau(10)/hpdpab* &
     &                 hdss*odss*oab(10)*gdpj*gdpj

                  dc1bj(10,1) = dc1bj(10,1) + 1.0_dp/hpdpab * &
     &              ( dhpdpab * cb * rxab/rab - 0.25_dp * dtau(10) * &
     &               ( -doab(10)*(-rxab/rab)*odss*hdss*gdpj*gdpj &
     &                 -oab(10)*odss*hdss*2.0_dp*gdpj*dgdpj(1) ))

                  dc1bj(10,2) = dc1bj(10,2) + 1.0_dp/hpdpab * &
     &              ( dhpdpab * cb * ryab/rab - 0.25_dp * dtau(10) * &
     &               ( -doab(10)*(-ryab/rab)*odss*hdss*gdpj*gdpj &
     &                 -oab(10)*odss*hdssjk*2.0_dp*gdpj*dgdpj(2) ))

                  dc1bj(10,3) = dc1bj(10,3) + 1.0_dp/hpdpab * &
     &              ( dhpdpab * cb * rzab/rab - 0.25_dp * dtau(10) * &
     &               ( -doab(10)*(-rzab/rab)*odss*hdss*gdpj*gdpj &
     &                 -oab(10)*odss*hdss*2.0_dp*gdpj*dgdpj(3) ))

               endif


! *** DDS
               if (hddsab /= 0.0) then
                  cb = -0.5_dp*dtau(12)/hddsab* &
     &                 hdss*odss*oab(12)*gdsj*gdsj

                  dc1bj(12,1) = dc1bj(12,1) + 1.0_dp/hddsab * &
     &              ( dhddsab * cb * rxab/rab + 0.5_dp * dtau(12) * &
     &               ( -doab(12)*(-rxab/rab)*odss*hdss*gdsj*gdsj &
     &                 -oab(12)*odss*hdss*2.0_dp*gdsj*dgdsj(1) ))

                  dc1bj(12,2) = dc1bj(12,2) + 1.0_dp/hddsab * &
     &              ( dhddsab * cb * ryab/rab + 0.5_dp * dtau(12) * &
     &               ( -doab(12)*(-ryab/rab)*odss*hdss*gdsj*gdsj &
     &                 -oab(12)*odss*hdss*2.0_dp*gdsj*dgdsj(2) ))

                  dc1bj(12,3) = dc1bj(12,3) + 1.0_dp/hddsab * &
     &              ( dhddsab * cb * rzab/rab + 0.5_dp * dtau(12) * &
     &               ( -doab(12)*(-rzab/rab)*odss*hdss*gdsj*gdsj &
     &                 -oab(12)*odss*hdss*2.0_dp*gdsj*dgdsj(3) ))

               endif


! *** DDP
               if (hddpab /= 0.0) then
                  cb = -0.25_dp*dtau(13)/hddpab* &
     &                 hdss*odss*oab(13)*gdpj*gdpj

                  dc1bj(13,1) = dc1bj(13,1) + 1.0_dp/hddpab * &
     &              ( dhddpab * cb * rxab/rab + 0.25_dp * dtau(13) * &
     &               ( -doab(13)*(-rxab/rab)*odss*hdss*gdpj*gdpj &
     &                 -oab(13)*odss*hdss*2.0_dp*gdpj*dgdpj(1) ))

                  dc1bj(13,2) = dc1bj(13,2) + 1.0_dp/hddpab * &
     &              ( dhddpab * cb * ryab/rab + 0.25_dp * dtau(13) * &
     &               ( -doab(13)*(-ryab/rab)*odss*hdss*gdpj*gdpj &
     &                 -oab(13)*odss*hdss*2.0_dp*gdpj*dgdpj(2) ))

                  dc1bj(13,3) = dc1bj(13,3) + 1.0_dp/hddpab * &
     &              ( dhddpab * cb * rzab/rab + 0.25_dp * dtau(13) * &
     &               ( -doab(13)*(-rzab/rab)*odss*hdss*gdpj*gdpj &
     &                 -oab(13)*odss*hdss*2.0_dp*gdpj*dgdpj(3) ))

               endif

! *** DDD
               if (hdddab /= 0.0) then
                  cb = -0.25_dp*dtau(14)/hdddab* &
     &                 hdss*odss*oab(14)*gddj*gddj

                  dc1bj(14,1) = dc1bj(14,1) + 1.0_dp/hdddab * &
     &              ( dhdddab * cb * rxab/rab + 0.25_dp * dtau(14) * &
     &               ( -doab(14)*(-rxab/rab)*odss*hdss*gddj*gddj &
     &                 -oab(14)*odss*hdss*2.0_dp*gddj*dgddj(1) ))

                  dc1bj(14,2) = dc1bj(14,2) + 1.0_dp/hdddab * &
     &              ( dhdddab * cb * ryab/rab + 0.25_dp * dtau(14) * &
     &               ( -doab(14)*(-ryab/rab)*odss*hdss*gddj*gddj &
     &                 -oab(14)*odss*hdss*2.0_dp*gddj*dgddj(2) ))

                  dc1bj(14,3) = dc1bj(14,3) + 1.0_dp/hdddab * &
     &              ( dhdddab * cb * rzab/rab + 0.25_dp * dtau(14) * &
     &               ( -doab(14)*(-rzab/rab)*odss*hdss*gddj*gddj &
     &                 -oab(14)*odss*hdss*2.0_dp*gddj*dgddj(3) ))

               endif

            endif

         endif
         ja = ja + 1

      enddo



! Let's finally calculate the derivative of the screening function.

      do i=1,14

         oab2(i) = oab(i)*oab(i)
         mu2(i) = mu2a(i)+mu2b(i)

         dmu2(i,1) = dmu2a(i,1) + dmu2b(i,1)
         dmu2(i,2) = dmu2a(i,2) + dmu2b(i,2)
         dmu2(i,3) = dmu2a(i,3) + dmu2b(i,3)

         cden = 1.0_dp - oab2(i) - mu2(i)
         cden2 = cden * cden

         pdoab(i,1) = 2.0_dp*oab(i)*doab(i)*(-rxab/rab)
         pdoab(i,2) = 2.0_dp*oab(i)*doab(i)*(-ryab/rab)
         pdoab(i,3) = 2.0_dp*oab(i)*doab(i)*(-rzab/rab)

         dc1(i,1) = dc1a(i,1) + dc1bi(i,1) + dc1bj(i,1)
         dc1(i,2) = dc1a(i,2) + dc1bi(i,2) + dc1bj(i,2)
         dc1(i,3) = dc1a(i,3) + dc1bi(i,3) + dc1bj(i,3)


! The final formula varies according to 
! various cut-offs for the screening function.
! If the SCF_CUT=0 there is no cut-off.

         if ( scf_cut == 0 ) then

            do j=1,3

               dsf(j,i) = -(dc1(i,j)-pdoab(i,j)-0.50_dp*dmu2(i,j))/cden &
     &                    + scf(i)*(-pdoab(i,j)-dmu2(i,j))/cden 

            enddo

! Step cut-off: if SCF>1 => SCF=1

         elseif ( scf_cut == 1 ) then

            do j=1,3
               if ( scf(i)  >  1.0 ) then
                  dsf(j,i) = 0.0_dp
               else

               dsf(j,i) = -(dc1(i,j)-pdoab(i,j)-0.50_dp*dmu2(i,j))/cden &
     &                    + scf(i)*(-pdoab(i,j)-dmu2(i,j))/cden 
               
               endif
            enddo

         elseif ( scf_cut == 2 ) then

            if (scf(i) > 20.0) then
               ff = 0.0_dp
               dff = 0.0_dp
            elseif (scf(i) < -20.0) then
               ff = 1.0_dp
               dff = 0.0_dp
               write(6,'(/" *** WARNING ***")')
               write(6,'("Screening function has a large negative value")')
               write(6,'(" SF = ",F7.3/)') scf(i)
            else
               ff = 1.0_dp/(1.0_dp+exp(10.0_dp*(scf(i)-1.5_dp))) 
               dff = -10.0*exp(10.0_dp*(scf(i)-1.5_dp))
            endif


            do j=1,3

               dsf(j,i) =(-(dc1(i,j)-pdoab(i,j)-0.50_dp*dmu2(i,j))/cden &
     &                    + scf(i)*(-pdoab(i,j)-dmu2(i,j))/cden) &
     &                    *ff*(1.0_dp+(1.0_dp-scf(i))*ff*dff)
            enddo

! Hyperbolic tangent cut-off.

         elseif ( scf_cut == 3 ) then

            do j=1,3

               dsf(j,i) =(-(dc1(i,j)-pdoab(i,j)-0.50_dp*dmu2(i,j))/cden &
     &                    + scf(i)*(-pdoab(i,j)-dmu2(i,j))/cden) &
     &                    /cosh(scf(i))**2.0_dp
            enddo

         endif



! The following is for testing purposes
!            DO J=1,3
!               DSF(j,i) = -DMU2A(I,J)
!              DSF(j,i) = -DC1BI(I,J)
!            ENDDO

      enddo


      end

