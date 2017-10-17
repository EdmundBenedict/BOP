 
      subroutine screenf(ia,ib,sf)
          use mod_precision


          use mod_all_scalar

          use mod_const
          use mod_atom_ar
!
!    This is a subroutine which calculates the screening function
!    for a bond between atoms A-B.
!    Unfortunately it is quite a mess :-((.
!    (c) Matous Mrovec, Oxford, June 2000
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
!      include "Include/Hamilt.array"
!      include "Include/Moment.array"
      include "Include/NebList.array"
      include "Include/PosVel.array"

!
!    Declare the simple variables.
!

      real(dp) :: rxab,ryab,rzab,rab,rik,rjk
      real(dp) :: rxik,ryik,rzik,rxjk,ryjk,rzjk
      real(dp) :: num,den,nod,theta,thetaj,c0,ca,cb
      real(dp) :: osss,opss,osps,odss
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

      real(dp) :: cos_i,sin_i,cos2_i,sin2_i
      real(dp) :: cos_j,sin_j,cos2_j,sin2_j

      real(dp) :: pom1,pom2
      real(dp) :: ff

      real(dp) :: am_c3,am_c4,am_c5,am_theta,am_s

      integer ia,ib
      integer ja,j,i,k,idbg
      integer bt,btj


!
!    Declare the local arrays.
!

      real(dp) :: sf(14)
      real(dp) :: mu2a(14),mu2b(14),oab(14),ojk(14),oab2(14)
      real(dp) :: mu2(14),c1(14),c1a(14),c1bi(14),c1bj(14),dtau(14)
!
!      EXTERNAL SCALE
!
! Initialization of variables

      idbg = 0

      pom1=0.0_dp
      pom2=0.0_dp

      do i=1,14
         mu2a(i) = 0.0_dp
         mu2b(i) = 0.0_dp
         c1(i) = 0.0_dp
         c1a(i) = 0.0_dp
         c1bi(i) = 0.0_dp
         c1bj(i) = 0.0_dp
        oab2(i) = 0.0_dp
      enddo

! Pair of atoms I-J

      rxab = ad( 1, ib ) - ad( 1, ia )
      ryab = ad( 2, ib ) - ad( 2, ia )
      rzab = ad( 3, ib ) - ad( 3, ia )               
      rab = sqrt(rxab*rxab+ryab*ryab+rzab*rzab)

      bt = btype(z(ia),z(ib))

      do i=1,14,1
         dtau(i) = bndscl(14,i,bt)
         if (bndscl(1,i,bt) == 0.0_dp) then
            oab(i) = 0.0_dp 
         else
            oab(i) = vscale(rab,bndscl(1,i,bt),bndscl(2,i,bt), & 
     &                      bndscl(5,i,bt),bndscl(6,i,bt), & 
     &                      bndscl(3,i,bt),bndscl(4,i,bt), & 
     &                      bndscl(7,i,bt),bndscl(8,i,bt), & 
     &                      bndscl(9,i,bt),bndscl(10,i,bt), & 
     &                      bndscl(11,i,bt),bndscl(12,i,bt))

            oab(i) = oab(i)*bndscl(13,i,bt)

         endif
      enddo

      call radf(rab,hsssab,hspsab,hpssab,hppsab, & 
     &          hpppab,hsdsab,hdssab,hpdsab,hdpsab, & 
     &          hpdpab,hdppab,hddsab,hddpab,hdddab, & 
     &          bnddat(1,bt),bndscl(1,1,bt))


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
            rik = sqrt(rxik*rxik+ryik*ryik+rzik*rzik)

            num = rxik*rxab+ryik*ryab+rzik*rzab
            den = rab*rik
            nod = num/den

            if (nod > 1.0_dp-1.0e-6_dp) then
               theta = 0.0_dp
            elseif (nod < -1.0_dp+1.0e-6_dp) then
               theta = pi
            else
               theta = acos(nod)
            endif

            cos_i = cos(theta)
            sin_i = sin(theta)
            cos2_i = cos(2.0_dp*theta)
            sin2_i = sin(2.0_dp*theta)


            bt = btype(z(ia),z(k))

            call radf(rik,hsss,hsps,hpss,hpps, & 
     &                hppp,hsds,hdss,hpds,hdps, & 
     &                hpdp,hdpp,hdds,hddp,hddd, & 
     &                bnddat(1,bt),bndscl(1,1,bt))



! Pair of atoms J-K
! Hop from j to k

            rxjk = ad( 1, k ) - ad( 1, ib )
            ryjk = ad( 2, k ) - ad( 2, ib )
            rzjk = ad( 3, k ) - ad( 3, ib )               
            rjk = sqrt(rxjk*rxjk+ryjk*ryjk+rzjk*rzjk)

            num = rxjk*rxab+ryjk*ryab+rzjk*rzab
            den = rab*rjk
            nod = -num/den

            if (nod > 1.0_dp-1.0e-6_dp) then
               thetaj = 0.0_dp
            elseif (nod < -1.0_dp+1.0e-6_dp) then
               thetaj = pi
            else
               thetaj = acos(nod)
            endif

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


            btj = btype(z(ib),z(k))
!            BTJ = BTYPE(Z(K),Z(IB))

            do i=1,14,1
               if ((bndscl(1,i,btj) == 0.0_dp).or.(rjk > rcut)) then
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


            if (idbg == 1) then
            if ((rik < rcut).or.(rjk < rcut)) then
               write(72,'("-------------")')
               write(72,'(3I6)') ia,ib,k
               write(72,'("Rij=",F6.4," Rik=",F6.4," Rjk=",F6.4)')  &
     &                      rab,rik,rjk
               write(72,'("Bond type ik=",I1," Bond type jk=",I1)')  &
     &                      bt,btj
            endif
            endif


! Check whether we have s electrons

            if (bndscl(1,1,bt) == 0.0_dp) then
               osss = 0.0_dp
            else
               osss = vscale(rik,bndscl(1,1,bt),bndscl(2,1,bt), & 
     &                      bndscl(5,1,bt),bndscl(6,1,bt), & 
     &                      bndscl(3,1,bt),bndscl(4,1,bt), & 
     &                      bndscl(7,1,bt),bndscl(8,1,bt), & 
     &                      bndscl(9,1,bt),bndscl(10,1,bt), & 
     &                      bndscl(11,1,bt),bndscl(12,1,bt))

               osss = osss*bndscl(13,1,bt)

               gss = 1.0_dp

! Calculate second moment contributions around atom I - MU2A(I)

               mu2a(1) = mu2a(1) + osss*osss*gss*gss
               mu2a(2) = mu2a(1)
               mu2a(6) = mu2a(1)

! Calculate second moment contributions around atom I - C1BI(I)

! *** SSS
               if (hsssab /= 0.0_dp) then
                  cb = 0.5_dp*dtau(1)/hsssab* &
     &                 (-hsss)*osss*oab(1)
                  c1bi(1) = c1bi(1)+cb
               endif

! *** SPS
               if (hspsab /= 0.0_dp) then
                  cb = 0.5_dp*dtau(2)/hspsab* &
     &                 (-hsss)*osss*oab(2)
                  c1bi(2) = c1bi(2) + cb
               endif

! *** SDS
               if (hsdsab /= 0.0_dp) then
                  cb = 0.5_dp*dtau(6)/hsdsab* &
     &                 (-hsss)*osss*oab(6)
                  c1bi(6) = c1bi(6) + cb
               endif


! Calculate third moment contributions - C1A

               if (rjk < rcut) then

! *** SSS
               if (hsssab /= 0.0_dp) then
                  ca = 0.5_dp*dtau(1)/hsssab* &
     &                 (hsss*ojk(1)+osss*hsssjk)
                  c1a(1) = c1a(1)+ca
               endif

! *** SPS
               if (hspsab /= 0.0_dp) then
                  ca = 0.5_dp*dtau(2)/hspsab* &
     &                 (hsss*ojk(3)+osss*hpssjk)*gpsj
                  c1a(2) = c1a(2) + ca
               endif

! *** SDS
               if (hsdsab /= 0.0) then
                  ca = 0.5_dp*dtau(6)/hsdsab* &
     &                 (hsss*ojk(7)+osss*hdssjk)*gdsj
                  c1a(6) = c1a(6) + ca
               endif

            endif

            endif


! Check whether we have p electrons

            if (bndscl(1,3,bt) == 0.0_dp) then
               opss = 0.0_dp
            else
               opss = vscale(rik,bndscl(1,3,bt),bndscl(2,3,bt), & 
     &                      bndscl(5,3,bt),bndscl(6,3,bt), & 
     &                      bndscl(3,3,bt),bndscl(4,3,bt), & 
     &                      bndscl(7,3,bt),bndscl(8,3,bt), & 
     &                      bndscl(9,3,bt),bndscl(10,3,bt), & 
     &                      bndscl(11,3,bt),bndscl(12,3,bt))

               opss = opss*bndscl(13,3,bt)

            if (bndscl(1,2,bt) == 0.0_dp) then
               osps = 0.0_dp
            else
               osps = vscale(rik,bndscl(1,2,bt),bndscl(2,2,bt), & 
     &                      bndscl(5,2,bt),bndscl(6,2,bt), & 
     &                      bndscl(3,2,bt),bndscl(4,2,bt), & 
     &                      bndscl(7,2,bt),bndscl(8,2,bt), & 
     &                      bndscl(9,2,bt),bndscl(10,2,bt), & 
     &                      bndscl(11,2,bt),bndscl(12,2,bt))

               osps = osps*bndscl(13,2,bt)

            endif

               gps = cos_i
               gpp = sin_i

! Calculate second moment contributions around atom I - MU2A(I)

               mu2a(3) = mu2a(3) + opss*opss*gps*gps
               mu2a(4) = mu2a(3)
               mu2a(5) = mu2a(5) + 0.5*opss*opss*gpp*gpp
               mu2a(8) = mu2a(3)
               mu2a(10) = mu2a(5)

! Calculate second moment contributions around atom I - C1BI(I)

! *** PSS
               if (hpssab /= 0.0_dp) then
                  cb = 0.5_dp*dtau(3)/hpssab* &
     &                 (-hpss)*opss*oab(3)*gps*gps
                  c1bi(3) = c1bi(3) + cb
               endif

! *** PPS
               if (hppsab /= 0.0_dp) then
                  cb = 0.5_dp*dtau(4)/hppsab* &
     &                 (-hpss)*opss*oab(4)*gps*gps
                  c1bi(4) = c1bi(4) + cb
               endif

! *** PPP
               if (hpppab /= 0.0_dp) then
                  cb = -0.25_dp*dtau(5)/hpppab* &
     &                 (-hpss)*opss*oab(5)*gpp*gpp
                  c1bi(5) = c1bi(5) + cb
               endif

! *** PDS
               if (hpdsab /= 0.0_dp) then
                  cb = 0.5_dp*dtau(8)/hpdsab* &
     &                 (-hpss)*opss*oab(8)*gps*gps
                  c1bi(8) = c1bi(8) + cb
               endif

! *** PDP
               if (hpdpab /= 0.0_dp) then
                  cb = -0.25_dp*dtau(10)/hpdpab* &
     &                 (-hpss)*opss*oab(10)*gpp*gpp
                  c1bi(10) = c1bi(10) + cb
               endif


! Calculate third moment contributions - C1A

               if (rjk < rcut) then

! *** PSS
               if (hpssab /= 0.0_dp) then
                  ca = 0.5_dp*dtau(3)/hpssab* &
     &                 (hpss*ojk(1)+opss*hsssjk)*gps
                  c1a(3) = c1a(3) + ca
               endif

! *** PPS
               if (hppsab /= 0.0_dp) then
                  ca = -0.5_dp*dtau(4)/hppsab* &
     &                 (hpss*ojk(3)+opss*hpssjk)*gps*gpsj
                  c1a(4) = c1a(4) + ca
               endif

! *** PPP
               if (hpppab /= 0.0_dp) then
                  ca = 0.25_dp*dtau(5)/hpppab* &
     &                 (hpss*ojk(3)+opss*hpssjk)*gpp*gppj
                  c1a(5) = c1a(5) + ca
               endif

! *** PDS
               if (hpdsab /= 0.0_dp) then
                  ca = 0.5_dp*dtau(8)/hpdsab* &
     &                 (hpss*ojk(7)+opss*hdssjk)*gps*gdsj
                  c1a(8) = c1a(8) + ca
               endif

! *** PDP
               if (hpdpab /= 0.0_dp) then
                  ca = -0.25_dp*dtau(10)/hpdpab* &
     &                 (hpss*ojk(7)+opss*hdssjk)*gpp*gdpj
                  c1a(10) = c1a(10) + ca
               endif

            endif

            endif


! Check whether we have d electrons

            if (bndscl(1,7,bt) == 0.0_dp) then
               odss = 0.0_dp
            else
               odss = vscale(rik,bndscl(1,7,bt),bndscl(2,7,bt), & 
     &                      bndscl(5,7,bt),bndscl(6,7,bt), & 
     &                      bndscl(3,7,bt),bndscl(4,7,bt), & 
     &                      bndscl(7,7,bt),bndscl(8,7,bt), & 
     &                      bndscl(9,7,bt),bndscl(10,7,bt), & 
     &                      bndscl(11,7,bt),bndscl(12,7,bt))

               odss = odss*bndscl(13,7,bt)


               gds = 0.25_dp*(1.0_dp+3.0_dp*cos2_i)
               gdp = sqrt3/2.0_dp*sin2_i
               gdd = sqrt3/4.0_dp*(1.0_dp-cos2_i)

! Calculate second moment contributions around atom I - MU2A(I)

               mu2a(7) = mu2a(7) + odss*odss*gds*gds
               mu2a(9) = mu2a(7)
               mu2a(11) = mu2a(11) + 0.5_dp*odss*odss*gdp*gdp
               mu2a(12) = mu2a(7)
               mu2a(13) = mu2a(11)
               mu2a(14) = mu2a(14) + 0.5_dp*odss*odss*gdd*gdd


! Calculate second moment contributions around atom I - C1BI(I)

! *** DSS
               if (hdssab /= 0.0_dp) then
                  cb = 0.5_dp*dtau(7)/hdssab* &
     &                 (-hdss)*odss*oab(7)*gds*gds
                  c1bi(7) = c1bi(7) + cb
               endif

! *** DPS
               if (hdpsab /= 0.0_dp) then
                  cb = 0.5_dp*dtau(9)/hdpsab* &
     &                 (-hdss)*odss*oab(9)*gds*gds
                  c1bi(9) = c1bi(9) + cb
!                  write(72,*) "DPS C1BI: ",cb,c1b(9)
               endif

! *** DPP
               if (hdppab /= 0.0_dp) then
                  cb = -0.25_dp*dtau(11)/hdppab* &
     &                 (-hdss)*odss*oab(11)*gdp*gdp
                  c1bi(11) = c1bi(11) + cb
!                  write(72,*) "DPP C1BI: ",cb,c1b(11)

               endif

! *** DDS
               if (hddsab /= 0.0_dp) then
                  cb = 0.5_dp*dtau(12)/hddsab* &
     &                 (-hdss)*odss*oab(12)*gds*gds
                  c1bi(12) = c1bi(12) + cb
!                  write(72,*) "DDS C1BI(I): ",cb,c1b(12)
               endif

! *** DDP
               if (hddpab /= 0.0_dp) then
                  cb = 0.25_dp*dtau(13)/hddpab* &
     &                 (-hdss)*odss*oab(13)*gdp*gdp
                  c1bi(13) = c1bi(13) + cb
!                  write(72,*) "DDP C1BI(I): ",cb,c1b(13)
               endif

! *** DDD
               if (hdddab /= 0.0_dp) then
                  cb = 0.25_dp*dtau(14)/hdddab* &
     &                 (-hdss)*odss*oab(14)*gdd*gdd
                  c1bi(14) = c1bi(14) + cb
!                  write(72,*) "DDD C1BI(I): ",cb,c1b(14)
               endif



! Calculate third moment contributions - C1A

               if (rjk < rcut) then

! *** DSS
                  if (hdssab /= 0.0_dp) then
                     ca = 0.5_dp*dtau(7)/hdssab* &
     &                 (hdss*ojk(1)+odss*hsssjk)*gds
                     c1a(7) = c1a(7) + ca
                  endif

! *** DPS
                  if (hdpsab /= 0.0_dp) then
                     ca = 0.5_dp*dtau(9)/hdpsab* &
     &                 (hdss*ojk(3)+odss*hpssjk)*gds*gpsj
                     c1a(9) = c1a(9) + ca
!                     write(72,*) "DPS ",ca,c1a(9)
                  endif

! *** DPP
                  if (hdppab /= 0.0_dp) then
                     ca = -0.25_dp*dtau(11)/hdppab* &
     &                    (hdss*ojk(3)+odss*hpssjk)*gdp*gppj
                     c1a(11) = c1a(11) + ca
!                     write(72,*) "DPP ",ca,c1a(11)
                  endif

! *** DDS
                  if (hddsab /= 0.0_dp) then
                     ca = 0.5_dp*dtau(12)/hddsab* &
     &                 (hdss*ojk(7)+odss*hdssjk)*gds*gdsj
                     c1a(12) = c1a(12) + ca
!                     write(72,*) "DDS C1A: ",ca,c1a(12)
                  endif

! *** DDP
                  if (hddpab /= 0.0_dp) then
                     ca = 0.25_dp*dtau(13)/hddpab* &
     &                 (-(hdss*ojk(7)+odss*hdssjk))*gdp*gdpj
                     c1a(13) = c1a(13) + ca
!                     write(72,*) "DDP C1A: ",ca,c1a(13)
                  endif
               
! *** DDD
                  if (hdddab /= 0.0_dp) then
                     ca = 0.25_dp*dtau(14)/hdddab* &
     &                    (hdss*ojk(7)+odss*hdssjk)*gdd*gddj
                     c1a(14) = c1a(14) + ca
!                     write(72,*) "DDD C1A: ",ca,c1a(14)
                  endif

               endif

            endif

         endif
         ja = ja + 1
      enddo     




!
!       For each ion (B) do almost the same as with (A)
!

      ja = aptr( ib )

      do while ( bptr( ja )  /=  eol )
         k = bptr( ja )

         if ((k /= ia).and.(k /= ib)) then
            rxik = ad( 1, k ) - ad( 1, ib )
            ryik = ad( 2, k ) - ad( 2, ib )
            rzik = ad( 3, k ) - ad( 3, ib )               
            rik = sqrt(rxik*rxik+ryik*ryik+rzik*rzik)
            num = rxik*rxab+ryik*ryab+rzik*rzab
            den = rab*rik
            nod = -num/den

            if (nod > 1.0_dp-1.0e-6_dp) then
               theta = 0.0_dp
            elseif (nod < -1.0_dp+1.0e-6_dp) then
               theta = pi
            else
               theta = acos(nod)
            endif

            bt = btype(z(ib),z(k))

            call radf(rik,hsss,hsps,hpss,hpps, & 
     &                hppp,hsds,hdss,hpds,hdps, & 
     &                hpdp,hdpp,hdds,hddp,hddd, & 
     &                bnddat(1,bt),bndscl(1,1,bt))


! Check whether we have s electrons

            if (bndscl(1,1,bt) == 0.0_dp) then
               osss = 0.0_dp
            else
               osss = vscale(rik,bndscl(1,1,bt),bndscl(2,1,bt), & 
     &                      bndscl(5,1,bt),bndscl(6,1,bt), & 
     &                      bndscl(3,1,bt),bndscl(4,1,bt), & 
     &                      bndscl(7,1,bt),bndscl(8,1,bt), & 
     &                      bndscl(9,1,bt),bndscl(10,1,bt), & 
     &                      bndscl(11,1,bt),bndscl(12,1,bt))

               osss = osss*bndscl(13,1,bt)

               gss = 1.0_dp

! Calculate second moment contributions around atom J - MU2B(J)
               
               mu2b(1) = mu2b(1) + osss*osss*gss*gss
               mu2b(3) = mu2b(1)
               mu2b(7) = mu2b(1)

! Calculate second moment contributions around atom J - C1BJ(J)

! *** SSS
               if (hsssab /= 0.0_dp) then
                  cb = 0.5_dp*dtau(1)/hsssab* &
     &                 (-hsss)*osss*oab(1)
                  c1bj(1) = c1bj(1)+cb
               endif

! *** SPS (=PSS)
               if (hpssab /= 0.0_dp) then
                  cb = 0.5_dp*dtau(3)/hpssab* &
     &                 (-hsss)*osss*oab(3)
                  c1bj(3) = c1bj(3) + cb
               endif

! *** SDS (=DSS)
               if (hdssab /= 0.0_dp) then
                  cb = 0.5_dp*dtau(7)/hdssab* &
     &                 (-hsss)*osss*oab(7)
                  c1bj(7) = c1bj(7) + cb
               endif

            endif


! Check whether we have p electrons

            if (bndscl(1,3,bt) == 0.0_dp) then
               opss = 0.0_dp
            else
               opss = vscale(rik,bndscl(1,3,bt),bndscl(2,3,bt), & 
     &                      bndscl(5,3,bt),bndscl(6,3,bt), & 
     &                      bndscl(3,3,bt),bndscl(4,3,bt), & 
     &                      bndscl(7,3,bt),bndscl(8,3,bt), & 
     &                      bndscl(9,3,bt),bndscl(10,3,bt), & 
     &                      bndscl(11,3,bt),bndscl(12,3,bt))

               opss = opss*bndscl(13,3,bt)

               gps = cos(theta)
               gpp = sin(theta)

! Calculate second moment contributions around atom J - MU2A(J)

               mu2b(2) = mu2b(2) + opss*opss*gps*gps
               mu2b(4) = mu2b(2)
               mu2b(5) = mu2b(5) + 0.5_dp*opss*opss*gpp*gpp
               mu2b(9) = mu2b(2)
               mu2b(11) = mu2b(5)

! Calculate second moment contributions around atom J - C1BJ(J)

! *** PSS (=SPS)
               if (hspsab /= 0.0_dp) then
                  cb = 0.5_dp*dtau(2)/hspsab* &
     &                 (-hpss)*opss*oab(2)*gps*gps
                  c1bj(2) = c1bj(2) + cb
               endif

! *** PPS (=PPS)
               if (hppsab /= 0.0_dp) then
                  cb = 0.5_dp*dtau(4)/hppsab* &
     &                 (-hpss)*opss*oab(4)*gps*gps
                  c1bj(4) = c1bj(4) + cb
               endif

! *** PPP (=PPP)
               if (hpppab /= 0.0_dp) then
                  cb = -0.25_dp*dtau(5)/hpppab* &
     &                 (-hpss)*opss*oab(5)*gpp*gpp
                  c1bj(5) = c1bj(5) + cb
               endif

! *** PDS (=DPS)
               if (hdpsab /= 0.0_dp) then
                  cb = 0.5_dp*dtau(9)/hdpsab* &
     &                 (-hpss)*opss*oab(9)*gps*gps
                  c1bj(9) = c1bj(9) + cb
               endif

! *** PDP (=DPP)
               if (hdppab /= 0.0_dp) then
                  cb = -0.25_dp*dtau(11)/hdppab* &
     &                 (-hpss)*opss*oab(11)*gpp*gpp
                  c1bj(11) = c1bj(11) + cb
               endif

            endif


! Check whether we have d electrons

            if (bndscl(1,7,bt) == 0.0_dp) then
               odss = 0.0_dp
            else
               odss = vscale(rik,bndscl(1,7,bt),bndscl(2,7,bt), & 
     &                      bndscl(5,7,bt),bndscl(6,7,bt), & 
     &                      bndscl(3,7,bt),bndscl(4,7,bt), & 
     &                      bndscl(7,7,bt),bndscl(8,7,bt), & 
     &                      bndscl(9,7,bt),bndscl(10,7,bt), & 
     &                      bndscl(11,7,bt),bndscl(12,7,bt))

               odss = odss*bndscl(13,7,bt)

               gds = 0.25_dp*(1.0_dp+3.0_dp*cos(2.0_dp*theta))
               gdp = sqrt3/2.0_dp*sin(2.0_dp*theta)
               gdd = sqrt3/4.0_dp*(1.0_dp-cos(2.0_dp*theta))


! Calculate second moment contributions around atom J - MU2A(J)

               mu2b(6) = mu2b(6) + odss*odss*gds*gds
               mu2b(8) = mu2b(6)
               mu2b(10) = mu2b(10) + 0.5_dp*odss*odss*gdp*gdp
               mu2b(12) = mu2b(6)
               mu2b(13) = mu2b(10)
               mu2b(14) = mu2b(14) + 0.5_dp*odss*odss*gdd*gdd

! Calculate second moment contributions around atom J - C1BJ(J)

! *** DSS (=SDS)
               if (hsdsab /= 0.0_dp) then
                  cb = 0.5_dp*dtau(6)/hsdsab* &
     &                 (-hdss)*odss*oab(6)*gds*gds
                  c1bj(6) = c1bj(6) + cb
               endif

! *** DPS (=PDS)
               if (hpdsab /= 0.0_dp) then
                  cb = 0.5_dp*dtau(8)/hpdsab* &
     &                 (-hdss)*odss*oab(8)*gds*gds
                  c1bj(8) = c1bj(8) + cb
!                  write(72,*) "DPS C1BJ: ",cb,c1b(9)
               endif

! *** DPP (=PDP)
               if (hpdpab /= 0.0_dp) then
                  cb = -0.25_dp*dtau(10)/hpdpab* &
     &                 (-hdss)*odss*oab(10)*gdp*gdp
                  c1bj(10) = c1bj(10) + cb
!                  write(72,*) "DPP C1BJ: ",cb,c1b(11)

               endif

! *** DDS
               if (hddsab /= 0.0_dp) then
                  cb = 0.5_dp*dtau(12)/hddsab* &
     &                 (-hdss)*odss*oab(12)*gds*gds
                  c1bj(12) = c1bj(12) + cb
!                  write(72,*) "DDS C1BJ: ",cb,c1b(12)
               endif

! *** DDP
               if (hddpab /= 0.0_dp) then
                  cb = 0.25_dp*dtau(13)/hddpab* &
     &                 (-hdss)*odss*oab(13)*gdp*gdp
                  c1bj(13) = c1bj(13) + cb
!                  write(72,*) "DDP C1BJ: ",cb,c1b(13)
               endif

! *** DDD
               if (hdddab /= 0.0_dp) then
                  cb = 0.25_dp*dtau(14)/hdddab* &
     &                 (-hdss)*odss*oab(14)*gdd*gdd
                  c1bj(14) = c1bj(14) + cb
!                  write(72,*) "DDD C1BJ: ",cb,c1b(14)
               endif


            endif

         endif
         ja = ja + 1

      enddo


! Let's finally calculate the screening function

      do i=1,14

         oab2(i) = oab(i)*oab(i)
         mu2(i) = mu2a(i)+mu2b(i)
         c1(i) = c1a(i)+c1bi(i)+c1bj(i)

         sf(i) = ( c1(i) - oab2(i) - 0.50_dp*mu2(i) ) / &
     &           ( 1 - oab2(i) - mu2(i) )

! The following is for testing purposes
!         SF(I) = C1(I) - MU2(I)
!         SF(I) = C1BI(I)
!         SF(I) = MU2A(I)

      enddo
!      PRINT *, "OAB4", OAB(4)
!      PRINT *, "MU2A", MU2A(4)
!      PRINT *, "MU2B", MU2B(4)
!      PRINT *, "C1A", C1A(4)
!      PRINT *, "C1BI", C1BI(4)
!      PRINT *, "C1BJ", C1BJ(4)
!      PRINT *, "SF(4)", SF(4)
      end


