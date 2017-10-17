 
      subroutine moldyn()
          use mod_precision
          use mod_all_scalar
          use mod_const
          use mod_atom_ar, only : mass

!
!    This is a subroutine to carry out molecular dynamics.
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
      include "Include/BEC.array"
!      include "Include/BondOrder.array"
      include "Include/Force.array"
!       include "Include/Hamilt.array"
      include "Include/Misc.array"
!      include "Include/Moment.array"
      include "Include/NebList.array"
!      include "Include/NotSRT.array"
      include "Include/PosVel.array"
!      include "Include/RecurCoef.array"
      include "Include/Relax.array"
!      include "Include/SRT.array"

!
!    Declare the simple variables.
!

      real(dp) :: dmu
      real(dp) :: x
      real(dp) :: eout, u0
      real(dp) :: t,p,vm,vmm
      real(dp) :: instem,insprs,avennb
      real(dp) :: op, orderp
      real(dp) :: msd_disp, msd

!     DOUBLE PRECISION TRET

      integer usrexit
      integer i,ia,flag

      character*80 filename

!
!    Carry out molecular dynamics.
!
!    Open the configuration file.
!

      if ((writeper > 0).and.(nequil < mxiter)) then
         filename = genfile(1:lengfn)//'.dat'
         open(unit = 2,file = filename,status = 'NEW')
!        WRITE(2,'(''ND''/I6)') ND
!        WRITE(2,'(''DT''/G12.5)') DT
!        WRITE(2,'(''RLIM''/3I5)') RLIM
         if (datform == 1) then
            write(2,'(''no_of_atoms ='',I6)') nd
            write(2,'(''no_of_frames ='',I4)') (mxiter-nequil)/writeper
            write(2,'(''data ='')')
         endif
      endif

!
!    Open the trace file.
!

      if (trace == 1) then
         filename = genfile(1:lengfn)//'.trace'
         open(unit = 3,file = filename,status = 'NEW')
      endif

!
!   Open the order parameter file.
!
      if (mda_op == 1) then
         filename = genfile(1:lengfn)//'.op'
         open(unit = 4,file = filename,status = 'NEW')
      endif


! 
!   Open the mean squared displacement file.
!
      if (mda_msd > 0) then
         filename = genfile(1:lengfn)//'.msd'
         open(unit = 7,file = filename,status = 'NEW')
         if (mda_msd == 1) then
            write(7,'(A,I5,A,3f12.5)') '# Original position of atom ', & 
     &           msd_atom,' was ',(admsda(i),i=1,3)
         else if (mda_msd == 2) then
            write(7,'(A)') '# MSD averaged over all atoms.'
            write(7,'(A)') '# Units of Angstroms^2.'
         endif
      endif

!
!    Evaluate the energy and forces of the initial configuration.
!

      if (mvflag == 4) then
         call setpress()
      else
         flag = 1
         call getetot(flag)
      endif

      call erasab()

      if ((mvflag == 2).or.(mvflag == 3)) then
         do i = 1,nd,1
            newftot(1,i) = ftot(1,i)
            newftot(2,i) = ftot(2,i)
            newftot(3,i) = ftot(3,i)
         enddo
      endif

      t = instem(vel,nd,mass)
      p = insprs(nd,t,vol,rdfbs+rdfpp)
      vm = vol
      vmm = vol

      do while ((iter <= mxiter).and.(usrexit() == 0))

!        IF (ITER.EQ.2)  THEN
!           CALL PCLOCK(TRET)
!           WRITE(9,'(/''Time for first step = '',G12.5)') TRET
!        ENDIF

         if (mvflag == 2) then
            call mdnve(newftot,ftot,dad,vel,nd,mass,dt,1)
         elseif (mvflag == 3) then
            call mdnve(newftot,ftot,dad,vel,nd,mass,dt,1)
!           CALL MDNVT(FTOT,DAD,VEL,ND,MASS,DT,TEMP)
         elseif (mvflag == 4) then
            call mdnpt(ftot,ad,dad,vel,nd,mass,dt,temp, & 
     &                 vol,vm,vmm,press,p,mpiston,a)
         elseif (mvflag == 5) then
            call constr(ftot,nd,cnst_a,cnst_v,cnst_n,mcnst_n)
            call mdnvt(ftot,dad,vel,nd,mass,dt,temp)
            if (iter > nequil) then
               temp = temp*dtemp
               t = instem(vel,nd,mass)
               x = sqrt(temp/t)
               call rescale(vel,nd,x)
            endif
         endif

         if (dndm > 1.0e-6_dp) then
            dmu = (locc-totnia)/dndm
            do ia = 1,nd,1
               dmu = dmu + dmdl(1,ia)*dad(1,ia) & 
     &                   + dmdl(2,ia)*dad(2,ia) & 
     &                   + dmdl(3,ia)*dad(3,ia)
            enddo
            lef = lef + dmu
         endif

         do ia = 1,nd,1
            do i = 1,3
               ad(i,ia)    = ad(i,ia)    + dad(i,ia)
               if ((mda_adav > 0).or.(mda_msd > 0)) then
                  adnopb(i,ia)= adnopb(i,ia)+ dad(i,ia)
                  adsum(i,ia) = adsum(i,ia) + adnopb(i,ia)
               endif
            enddo
         enddo

         if (abs(locc-totnia)/real(nd, dp) > 0.05_dp) then
            flag = 1
         else
            flag = 2
         endif

         call getetot(flag)

!        IF (MVFLAG.EQ.2) CALL MDNVE(NEWFTOT,FTOT,DAD,VEL,ND,MASS,DT,2)
         if ((mvflag == 2).or.(mvflag == 3)) & 
     &                    call mdnve(newftot,ftot,dad,vel,nd,mass,dt,2)

         t = instem(vel,nd,mass)
         if (((mvflag == 3).or.(mvflag == 4)).and. & 
     &       (abs(t-temp)/temp > ttol)) then
            x = sqrt(temp/t)
            call rescale(vel,nd,x)
            t = temp
         endif
         p = insprs(nd,t,vol,rdfbs+rdfpp)
         ke = 1.5_dp*kb*t*real(nd-1, dp)

         eout = eprom + ebond + epair - eatom + ke - uent
         u0   = eout  + 0.5_dp*uent

         if (trace == 1) write(3,'(I5,'' '',8G13.5)') & 
     &                     iter,u0,eout,eout-ke,ke,t,p,vol,totnia

         if (mda_op == 1) then
            op = orderp()
            write(4,'(I5,F12.5)') iter, op
         endif

         if (mda_msd > 0) then
            if (mda_msd == 1) then
               msd_disp = (adnopb(1,msd_atom)-admsda(1))**2 + & 
     &                    (adnopb(2,msd_atom)-admsda(2))**2 +  & 
     &                    (adnopb(3,msd_atom)-admsda(3))**2
               write(7,'(I7,4F12.5)') iter, & 
     &              (adnopb(i,msd_atom),i=1,3),msd_disp
            else if (mda_msd == 2) then
               msd_disp = msd()
               write(7,'(I7,F12.5)') iter,msd_disp
            endif
         endif

         if ((mvflag >= 2).and.(mvflag <= 4).and. & 
     &       (iter > nequil)) then
            sume = sume + eout
            sumee = sumee + eout**2
            sumt = sumt + t
            sumtt = sumtt + t**2
            sump = sump + p
            sumpp = sumpp + p**2
            sumv = sumv + vol
            sumvv = sumvv + vol**2
            if (mda_gr > 0) call addgr(nd,totnd,ad,gr_bins, & 
     &                             mda_gr_n,mda_gr_r,mda_gr_dr)
            if (mda_g3 > 0) call addg3(ad,totnd,g3_bins,mda_g3_n, & 
     &                                  aptr,nd,bptr,mxnb,mda_g3_r)
            if (mda_sw > 0) call addsw(ad,totnd,nd,sw_bins,mda_sw_n, & 
     &                                  aptr,bptr,mxnb,mda_sw_de)
         endif

         if ((writeper > 0).and.(iter > nequil)) then
            if (mod(iter,writeper) == 0) then
!              WRITE(2,'(/''Iteration # '',I6)') ITER
!              WRITE(2,'(''A'')')
!              WRITE(2,'(3G14.5)') A(1,1),A(1,2),A(1,3),
!    +                             A(2,1),A(2,2),A(2,3),
!    +                             A(3,1),A(3,2),A(3,3)
!              WRITE(2,'(''RV'')')
!              WRITE(2,'(6G14.5)') (AD(1,I),AD(2,I),AD(3,I),
!    +                              VEL(1,I),VEL(2,I),VEL(3,I),
!    +                              I=1,ND,1)
!              WRITE(2,'(''Energy = '',G22.15)') EOUT
!              WRITE(2,'(''Temperature = '',G22.15)') T
!              WRITE(2,'(''Pressure = '',G22.15)') P
!              WRITE(2,'(''Volume = '',G22.15)') VOL
               if (datform == 1) then
                  do i = 1,nd,1
                     write(2,'(3(E10.5,1X),4(E7.2,1X))') & 
     &                    ad(1,i),ad(2,i),ad(3,i), & 
     &                    avs(1,z(i)),avs(2,z(i)), & 
     &                    avs(3,z(i)),avs(4,z(i))
                  enddo 
               else
                  write(2,'(I6)') nd
                  write(2,'(A)') ' '
                  do i = 1, nd, 1
                     write(2,'(A,3F12.5)') symb(i), & 
     &                    ad(1,i),ad(2,i),ad(3,i)
                  enddo
               endif
            endif
         endif

         if (monitorper > 0) then
            if (mod(iter,monitorper) == 0) then
               write(6,'(/''Iteration # '',I6)') iter
               write(6,'(''Average number of nearest'', & 
     &                   '' neighbours = '',G12.5)') & 
     &                   avennb(aptr,nd)
               write(6,'(''Average cluster size = '',G12.5)') & 
     &                   aveclusiz
               write(6,'(''Free energy, U(T=0) = '',2G22.15)')eout,u0
               write(6,'(''Number of electrons = '',G22.15)') totnia
               write(6,'(''Temperature = '',G22.15)') t
               write(6,'(''Pressure = '',G22.15)') p
               write(6,'(''Volume = '',G22.15)') vol
            endif
         endif

         if (autosave > 0) then
            if (mod(iter,autosave) == 0) call dump()
         endif
         iter = iter + 1

      enddo

      end

