 
      subroutine convrg()
          use mod_precision
          use mod_all_scalar
          use mod_const
          use mod_srt
          use ab_io
!           use mod_conf, only : run_conf_t
          use mod_ham
      use mod_atom_ar, only : btype
!        
!    This is a routine to carry out a convergence test.
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
!      include "Include/Misc.array"
!      include "Include/Moment.array"
      include "Include/NebList.array"
!      include "Include/NotSRT.array"
      include "Include/PosVel.array"
!      include "Include/RecurCoef.array"
!      include "Include/Relax.array"
!      include "Include/SRT.array"
!
      include "Include/ag.conn"
!
!    Declare the simple variables.
!
!       type(run_conf_t), intent(in) :: rtconf
      real(dp) :: eplus,eminus,dz,fhf,fnum,ebulk,numel
      real(dp) :: ktlo,kthi,dkt,efhi,def,ef,nel
!      DOUBLE PRECISION KTLO,KTHI,DKT,EFLO,EFHI,DEF,EF,NEL
      real(dp) :: ezero

      integer cflag,flag,za,zb,la,lb,nstta,nsttb
      integer iatom,ibatom,ix
      integer nbaseold
      integer i,j
      integer ia,ib,ja,ja0,jb,n,nkt,nef
      integer getz, zin
!
      character*80 filename
      character*2 chemsymb

!
!    Declare local arrays.
!

      real(dp) :: dr(3)
      real(dp) :: dv(3)
      real(dp) :: scf(14),scfcut(14),dsf(3,14)
      real(dp) :: dsf0(14),sf0(14),sfplus(14)
      real(dp) :: sfminus(14),dsfnum(14)
      real(dp) :: nmfile
!
      external getz

!
!    Prompt the user for the required convergence test.
!

 1     continue
      if (.not. quiet) then
        write(6,'(/''You have the following convergence tests:'')')
        write(6,'('' 0) Exit'')')
        write(6,'('' 1) Total energy'')')
        write(6,'('' 3) Force'')')
    !      WRITE(6,'('' 4) Derivative of chemical potential w.r.t.'',
    !     +          '' number of electrons'')')
    !      WRITE(6,'('' 5) Derivative of chemical potential w.r.t.'',
    !     +          '' atomic position'')')
    !      WRITE(6,'('' 6) Gradient of Hamiltonian matrix'')')
    !      WRITE(6,'('' 7) Read/write to /tmp'')')
    !      WRITE(6,'('' 8) Stability of finite temperature integrals'')')
    !      WRITE(6,'('' 9) Create output file for XBS plotting program'')')
    !      WRITE(6,'(''Enter your choice > '',$)')

    !Dimitar Pashov edit
         READ(5,*) CFLAG
!         write(6,'(a)') ' CFLAG = 1 automatically chosen :)'
      end if
       cflag = 1
!end edit

!      CFLAG=1


!      IF ((CFLAG.LT.0).OR.(CFLAG.GT.9)) THEN
!         WRITE(6,'(''Unknown option. Please try again.'')')
!         GOTO 1
!      ENDIF

      if (cflag == 0) return

      goto (10,20,30,40,50,60,70,80) cflag-1

!=================================================================================
!    Energy convergence test.
!    CFLAG = 0
!=================================================================================

      filename = genfile(1:lengfn)//'.eng'
!      CHEMSYMB = "Ir"
!      ZIN = GETZ(CHEMSYMB)
!      IF (ZC(ZIN).LT.10.0) THEN 
!         WRITE(FILENAME,'(A3,".",F3.1,".energy")') STRTYP, ZC(ZIN)
!      ELSEIF (ZC(ZIN).GE.10.0) THEN
!         WRITE(FILENAME,'(A3,".",F4.1,".energy")') STRTYP, ZC(ZIN)
!      ENDIF
!

! Dimitar: commented the writing of tial.eng
!       open(unit = 3,file = filename,status = 'NEW',err=9998)
      
      flag = 1
      call getetot(flag)
!      DO I = 1, ND
!         PRINT*, FTOT(1,I), FTOT(2,I), FTOT(3,I)
!      ENDDO

      etot = eprom+ebond+epair-eatom-uent
      if (.not. quiet) then
    !      PRINT *, "EPROM=", EPROM
    !      PRINT *, "EPROMS=", EPROMS
    !      PRINT *, "EPROMP=", EPROMP
    !      PRINT *, "EPROMD=", EPROMD

    !      WRITE (3,'(2F16.9)') LENA(1)*LENA(2)*LENA(3), ETOT
    !      WRITE (3,'(2F16.9)') LENA(3)/LENA(1), ETOT

    ! Dimitar: commented the writing of tial.eng and the closing of 3 after a paragraph
    !       write(3,'(7F16.9)') ebond/nd, eproms/nd, epromp/nd, epromd/nd,  & 
    !      &     eprom/nd, (ebond+eprom-eatom)/nd, etot/nd



    !      WRITE(3,'(8G16.9)') ZC(ZIN), EBOND/ND, EPROMS/ND, EPROMP/ND,
    !     +     EPROMD/ND, EPROM/ND,(EBOND+EPROM-EATOM)/ND, ETOT/ND
        write(*,'(/10X,"CALCULATION OF TOTAL ENERGY"//)')
        write(*,'("   Etot  = ",G17.10/)') etot
        write(*,'("   Etot/atom  = ",G17.10//)') etot/nd

        write(*,'("   Ebond = ",G17.10)') ebond
        write(*,'("   Eprom = ",G17.10)') eprom
        write(*,'("   Epair = ",G17.10)') epair
        write(*,'("   Eatom = ",G17.10)') eatom
        write(*,'("   Eent  = ",G17.10)') uent 

        
    !       close(3)
    !      CALL BOAVG(2)

        call outfil(1)
!         if (.not. rtconf % paral) then
            call outfil(2)
            call outfil(3)
!         end if
      end if

      return

!=================================================================================
!    Evaluate defect energies.
!    Note that a bulk energy convergence test must be carried out first.
!    CFLAG = 1
!=================================================================================

 10   filename = genfile(1:lengfn)//'.eng'
      open(unit = 2,file = filename,status = 'OLD',err=9999)

      filename = genfile(1:lengfn)
      nbaseold = 0

 12   read(2,*,end=11) nbase,nrec,ebulk

      write(6,'(''NBASE = '',I2,'' NREC = '',I2)') nbase,nrec

      if (nbase /= nbaseold) then
         if (nbaseold > 0) close(3)
         write(filename(lengfn+1:lengfn+2),'(''.'',I1)') nbase
         open(unit = 3,file = filename,status = 'NEW',err=9998)
      endif

      nbaseold = nbase

      flag = 1
      call getetot(flag)

      etot = eprom+ebond+epair-eatom-uent

      write(3,'(I4,G16.7)') 2*nrec+1,etot-(real(nd, dp)/64.0_dp)*ebulk

      goto 12

 11   close(2)
      close(3)

      return

!=============================================================================
!     Test of the derivative of the screening function.
!    CFlag = 2
!     Important: pbc must be switched or the block must be sufficiently
!     large so that the tested atoms do not interact with their images
!     for this test to be correct !!!
!=============================================================================

 20   filename = genfile(1:lengfn)//'.cnv'
      open(unit = 2,file = filename,status = 'OLD',err=9999)

      filename = genfile(1:lengfn)//'.F'
      open(unit = 3,file = filename,status = 'NEW',err=9998)

      write(6,'(''Enter atom number > '',$)')
      read(5,*) iatom

      write(6,'(''Enter neighbor atom number > '',$)')
      read(5,*) ibatom

      write(6,'(''Enter force component (1,2 or 3) > '',$)')
      read(5,*) ix

      if     ( ix  ==  1 )  then
          write(3,'("Force evaluation: for atom ",I2, &
     &              " in X direction")' )  iatom 
      elseif ( ix  ==  2 )  then
          write(3,'("Force evaluation: for atom ",I2, &
     &              " in Y direction")' )  iatom
      elseif ( ix  ==  3 )  then
          write(3,'("Force evaluation: for atom ",I2, &
     &              " in Z direction")' )  iatom
      else
          write(3,'(//"Do not act stupid!"//)' )
      endif


      dz = 1.0e-3_dp

      call bldnebt()


      call screenf(iatom,ibatom,scf)
      call scrcut(scf,scfcut)
      call dscreenf(iatom,ibatom,scf,dsf)

      do i=1,14
         sf0(i)=scfcut(i)
         dsf0(i)=dsf(ix,i)
      enddo

      write(79,*) " *** Calculating  +++DZ *** "
      ad(ix,iatom) = ad(ix,iatom) + dz

      call bldnebt()

      call screenf(iatom,ibatom,scf)
      call scrcut(scf,scfcut)
      do i=1,14
         sfplus(i)=scfcut(i)
      enddo

      ad(ix,iatom) = ad(ix,iatom) - dz

      write(79,*) " *** Calculating  ---DZ *** "
      ad(ix,iatom) = ad(ix,iatom) - dz

      call bldnebt()

      call screenf(iatom,ibatom,scf)
      call scrcut(scf,scfcut)
      do i=1,14
         sfminus(i)=scfcut(i)
      enddo

      ad(ix,iatom) = ad(ix,iatom) + dz


      do i=1,14
         dsfnum(i) = (sfplus(i)-sfminus(i))/(dz+dz)

         write(3,'(3F8.5)') sf0(i),sfplus(i),sfminus(i)
      enddo


      write(6,'(/"        DSF        NUMERIC       ABS.DIFF" &
     &           "    REL.DIFF. %"/)')
      do i=1,14
         write(6,'(''    '',2E13.5,2(3X,F10.4))') dsf0(i), &
     &      dsfnum(i),dsf0(i)-dsfnum(i),  &
     &      abs((dsf0(i)-dsfnum(i))/dsfnum(i))*100.0
      enddo

!      GOTO 22

 21   write(3,'(//)')
      close(2)
      close(3)

      return

!=================================================================================
!    Derivative of chemical potential with respect to number of electrons.
!    CFlag = 3
!=================================================================================

 30   filename = genfile(1:lengfn)//'.cnv'
      open(unit = 2,file = filename,status = 'OLD',err=9999)

      dz = 1.0e-3_dp

      filename = genfile(1:lengfn)//'.dmu'
      open(unit = 3,file = filename,status = 'NEW',err=9998)

      write(3,'('' M    Analytic      Numeric'')')

 32   read(2,*,end=31) nbase,nrec

      write(6,'(''NBASE = '',I2,'' NREC = '',I2)') nbase,nrec

      flag = 1
      call getetot(flag)
      fhf = 1.0_dp/dndm

      locc = locc + dz
      flag = 0
      call getetot(flag)
      eplus = lef
      locc = locc - dz

      locc = locc - dz
      flag = 0
      call getetot(flag)
      eminus = lef
      locc = locc + dz

      fnum = (eplus-eminus)/(dz+dz)

      write(3,'(I2,''  '',2G13.5)') 2*nrec+1,fhf,fnum

      goto 32

 31   close(2)
      close(3)

      return

!=================================================================================
!    Derivative of chemical potential with respect to atomic position.
!    CFlag = 4
!=================================================================================

 40   filename = genfile(1:lengfn)//'.cnv'
      open(unit = 2,file = filename,status = 'OLD',err=9999)

      filename = genfile(1:lengfn)//'.dmudl'
      open(unit = 3,file = filename,status = 'NEW',err=9998)

      write(6,'(''Enter atom number > '',$)')
      read(5,*) iatom

      write(6,'(''Enter displacement component (1,2 or 3) > '',$)')
      read(5,*) ix

      dz = 1.0e-3_dp

      write(3,'('' M    Analytic      Numeric'')')

 42   read(2,*,end=41) nbase,nrec

      write(6,'(''NBASE = '',I2,'' NREC = '',I2)') nbase,nrec

      flag = 1
      call getetot(flag)
      fhf = dmdl(ix,iatom)

      ad(ix,iatom) = ad(ix,iatom) + dz
      flag = 1
      call getetot(flag)
      eplus = lef
      ad(ix,iatom) = ad(ix,iatom) - dz

      ad(ix,iatom) = ad(ix,iatom) - dz
      flag = 1
      call getetot(flag)
      eminus = lef
      ad(ix,iatom) = ad(ix,iatom) + dz

      fnum = (eplus-eminus)/(dz+dz)

      write(3,'(I2,''  '',2G13.5)') 2*nrec+1,fhf,fnum

      goto 42

 41   close(2)
      close(3)

      return
      
! The calls to grdmat were long out of date before I changed its interface. To be commented out only after update
! Dimitar
!=================================================================================
!    Test GRDMAT
!    CFlag = 5
!=================================================================================
! 
!  50   write(6,'(''Enter atomic numbers of atoms a and b > '',$)')
!       read(5,*) za,zb
! 
!       write(6,'(''Enter atomic displacement vector > '',$)')
!       read(5,*) dr
! 
!       dz = 1.0e-3_dp
! 
!       filename = genfile(1:lengfn)//'.grd'
!       open(unit = 3,file = filename,status = 'NEW',err=9999)
! 
!       call grdmat(grad,dr,za,zb,mxnstat,ml,llista,llistb, & 
!      &            bnddat,btype,nbtype,mbtype)
! 
!       do la = 1,nstt(za),1
!          do lb = 1,nstt(zb),1
! 
!             dr(1) = dr(1) + dz
!             call matel(za,zb,dr,nstta,nsttb,subh,0.0_dp, & 
!      &                 mxnstat,ml,llista,llistb,bnddat, & 
!      &                 btype,nbtype,mbtype)
!             eplus = subh(la,lb)
!             dr(1) = dr(1) - dz
!             dr(1) = dr(1) - dz
!             call matel(za,zb,dr,nstta,nsttb,subh,0.0_dp, & 
!      &                 mxnstat,ml,llista,llistb,bnddat, & 
!      &                 btype,nbtype,mbtype)
!             eminus = subh(la,lb)
!             dr(1) = dr(1) + dz
!             dv(1) = (eplus-eminus)/(dz+dz)
! 
!             dr(2) = dr(2) + dz
!             call matel(za,zb,dr,nstta,nsttb,subh,0.0_dp, & 
!      &                 mxnstat,ml,llista,llistb,bnddat, & 
!      &                 btype,nbtype,mbtype)
!             eplus = subh(la,lb)
!             dr(2) = dr(2) - dz
!             dr(2) = dr(2) - dz
!             call matel(za,zb,dr,nstta,nsttb,subh,0.0_dp, & 
!      &                 mxnstat,ml,llista,llistb,bnddat, & 
!      &                 btype,nbtype,mbtype)
!             eminus = subh(la,lb)
!             dr(2) = dr(2) + dz
!             dv(2) = (eplus-eminus)/(dz+dz)
! 
!             dr(3) = dr(3) + dz
!             call matel(za,zb,dr,nstta,nsttb,subh,0.0_dp, & 
!      &                 mxnstat,ml,llista,llistb,bnddat, & 
!      &                 btype,nbtype,mbtype)
!             eplus = subh(la,lb)
!             dr(3) = dr(3) - dz
!             dr(3) = dr(3) - dz
!             call matel(za,zb,dr,nstta,nsttb,subh,0.0_dp, & 
!      &                 mxnstat,ml,llista,llistb,bnddat, & 
!      &                 btype,nbtype,mbtype)
!             eminus = subh(la,lb)
!             dr(3) = dr(3) + dz
!             dv(3) = (eplus-eminus)/(dz+dz)
! 
!             write(3,'(''['',I1,'','',I1,'']  '',3G13.5)') & 
!      &            la,lb,grad(1,la,lb),grad(2,la,lb),grad(3,la,lb)
!             write(3,'(''       '',3G13.5)') dv
! 
!          enddo
!       enddo

      close(3)

      return


!=================================================================================
!    Check disk read/write
!    CFlag = 6
!=================================================================================

 60   flag = 2
      call getetot(flag)

      ia = nd
      call assoc_ab(ia)
      call assoc_ham(ia)
      filename = genfile(1:lengfn)//'.v1'
      open(unit = 3,file = filename,status = 'NEW',err=9999)

      write(3,*) nchain
      do la = 1,nchain,1
         write(3,*) lchain(la)
         write(3,*) (arec(n,la),brec(n,la),n=0,lchain(la),1)
         write(3,*) lainf(la),lbinf(la)
         ja = aptr(ia)
         ja0 = ja
         do while(bptr(ja) /= eol)
            jb = ja-ja0+1
            ib = bptr(ja)
            do lb = 1,nstt(z(ib)),1
               write(3,*) (darec(n,la,lb,jb),dbrec(n,la,lb,jb), & 
     &                    n=0,lchain(la),1)
            enddo
            ja = ja + 1
         enddo
         if (term == 1) then
            write(3,*) forder(la)
            write(3,*) (root(i,la),f1(i,la),f2(i,la), & 
     &                 i=1,forder(la),1)
         elseif ((term == 2).or.(term == 3)) then
            write(3,*) (diag(i,la),i=1,lchain(la)+1,1)
            write(3,*) ((eigvec(i,j,la),i=1,lchain(la)+1,1), & 
     &                  j=1,lchain(la)+1,1)
         endif
      enddo

      close(3)

      call wrtab(1,ia)
      close(50)
      call rdab(1,ia)
      close(50)

      filename = genfile(1:lengfn)//'.v2'
      open(unit = 3,file = filename,status = 'NEW',err=9999)

      write(3,*) nchain
      do la = 1,nchain,1
         write(3,*) lchain(la)
         write(3,*) (arec(n,la),brec(n,la),n=0,lchain(la),1)
         write(3,*) lainf(la),lbinf(la)
         ja = aptr(ia)
         ja0 = ja
         do while(bptr(ja) /= eol)
            jb = ja-ja0+1
            ib = bptr(ja)
            do lb = 1,nstt(z(ib)),1
               write(3,*) (darec(n,la,lb,jb),dbrec(n,la,lb,jb), & 
     &                    n=0,lchain(la),1)
            enddo
            ja = ja + 1
         enddo
         if (term == 1) then
            write(3,*) forder(la)
            write(3,*) (root(i,la),f1(i,la),f2(i,la), & 
     &                 i=1,forder(la),1)
         elseif ((term == 2).or.(term == 3)) then
            write(3,*) (diag(i,la),i=1,lchain(la)+1,1)
            write(3,*) ((eigvec(i,j,la),i=1,lchain(la)+1,1), & 
     &                  j=1,lchain(la)+1,1)
         endif
      enddo

      close(3)

      return

!=============================================================================
!    Check integrals
!    CFlag = 7
!=============================================================================

 70   write(6,'(''Enter lowest and highest electron'', & 
     &          '' temperatures > '',$)')
      read(5,*) ktlo,kthi
       
      write(6,'(''Enter the number of temperatures > '',$)')
      read(5,*) nkt

      write(6,'(''Enter lowest and highest chemical potentials > '',$)')
      read(5,*) eflo,efhi
       
      write(6,'(''Enter the number of chemical potentials > '',$)')
      read(5,*) nef

      filename = genfile(1:lengfn)//'.int'
      open(unit = 3,file = filename,status = 'NEW',err=9999)

      ef = 0.0_dp
      kt = 0.5_dp
      flag = 1
      call bldnebt()
      nel = numel(flag,ef)

      dkt = (kthi-ktlo)/real(nkt-1, dp)
      def = (efhi-eflo)/real(nef-1, dp)

      kt = ktlo
      do j = 1,nkt,1
         write(6,'(''KT = '',G12.5)') kt
         write(3,'(/''KT = '',G12.5)') kt
         ef = eflo
         do i = 1,nef,1
            nel = numel(flag,ef)
            write(3,'(''EF = '',F8.2,'' N = '',G12.5)') ef,nel
            ef = ef + def
         enddo
         kt = kt + dkt
      enddo

      close(3)

      return

!=============================================================================
!    Create input file for XBS program.
!=============================================================================

!80    CALL DISTANCE
!80      CALL FORCECHECK

 80      filename = genfile(1:lengfn)//'.cnv'
      open(unit = 2,file = filename,status = 'OLD',err=9999)

      filename = genfile(1:lengfn)//'.F'
      open(unit = 3,file = filename,status = 'NEW',err=9998)

      write(6,'(''Enter atom number > '',$)')
      read(5,*) iatom

      write(6,'(''Enter force component (1,2 or 3) > '',$)')
      read(5,*) ix

      if     ( ix  ==  1 )  then
          write(3,'("Force evaluation: for atom ",I2, &
     &              " in X direction")' )  iatom 
      elseif ( ix  ==  2 )  then
          write(3,'("Force evaluation: for atom ",I2, &
     &              " in Y direction")' )  iatom
      elseif ( ix  ==  3 )  then
          write(3,'("Force evaluation: for atom ",I2, &
     &              " in Z direction")' )  iatom
      else
          write(3,'(//"Do not act stupid!"//)' )
      endif


      dz = 1.0e-3_dp

      write(3,'(/" M        HF         NUMERIC       ABS.DIFF" &
     &           "    REL.DIFF. %"/)')

 82      read(2,*,end=81) nbase,nrec

      write(6,'(''NBASE = '',I2,'' NREC = '',I2)') nbase,nrec

      flag = 1

      call getetot(flag)

      etot = eprom+ebond+epair-eatom-uent
      ezero = etot

      fhf = ftot(ix,iatom)

      ad(ix,iatom) = ad(ix,iatom) + dz
      flag = 1
      call getetot(flag)
      eplus = eprom+ebond+epair-eatom-uent
      ad(ix,iatom) = ad(ix,iatom) - dz

      ad(ix,iatom) = ad(ix,iatom) - dz
      flag = 1
      call getetot(flag)
      eminus = eprom+ebond+epair-eatom-uent
      ad(ix,iatom) = ad(ix,iatom) + dz

      fnum = -(eplus-eminus)/(dz+dz)

      write(3,'(I2,''  '',2E13.5,2(3X,F10.6)/)') 2*nrec+1,fhf,fnum, &
     &      fhf-fnum, abs((fhf-fnum)/fnum)*100.0

      write(6,'(/" M        HF         NUMERIC       ABS.DIFF" &
     &           "    REL.DIFF. %"/)')
      write(6,'(I2,''  '',2E13.5,2(3X,F10.6))') 2*nrec+1,fhf,fnum, &
     &      fhf-fnum, abs((fhf-fnum)/fnum)*100.0

      write(3,*) ezero
      write(3,*) eplus
      write(3,*) eminus
      write(3,*) eplus-eminus


      goto 82

 81      write(3,'(//)')
      close(2)
      close(3)

     
      return


!=================================================================================
!    Error opening a file.
!=================================================================================

 9998 close(2)
 9999 write(6,'(''Unable to open file '',A)') filename
      goto 1

      end


