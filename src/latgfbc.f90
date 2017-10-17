 
      subroutine latgfbc()
          use mod_precision


          use mod_all_scalar

          use mod_const
!
!     Before you start to use this subroutine to determine the lattice
!     Greens function boundary conditions, PLEASE read the paper by
!     Woodward, Rao et al (Phil Mag A, 77 (1998), p 231-256)!!!!
!
!     Work started on this code June 2003 at Northwestern University
!     with the kind help of Chris Woodward. 
!
!     This is a subroutine to calculate the lattice Greens
!     functions using a variable metric relaxation. A specified
!     force is applied to the atom at the origin and the displacements
!     of all active atoms is determined after relaxation. This is then
!     used to determine the elements Gij(r) for u_i(r) and f_j
!
!     The input block must contain the displacements from the elastic
!     Greens function solution corresponding to a specified line force
!     in the appropriate direction. You will need the package of 
!     programs in GFutils for this, both for determining the elastic
!     solution and then applying it your block of atoms. You'll need
!     elast.f and latticeloop.f for this - GOOD LUCK!!!
!
!     M. Cawkwell 16/6/2003

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

!      include "Include/Atom.array"
      include "Include/BEC.array"
!      include "Include/BondOrder.array"
      include "Include/Force.array"
!      include "Include/Hamilt.array"
!      include "Include/Misc.array"
!      include "Include/Moment.array"
!      include "Include/NebList.array"
!      include "Include/NotSRT.array"
      include "Include/PosVel.array"
!      include "Include/RecurCoef.array"
      include "Include/Relax.array"
!      include "Include/SRT.array"

!
!    Declare the simple variables.
!

      real(dp) :: ein,eout,ered,u0
      real(dp) :: alpha,maxf,magf,dmu
      real(dp) :: gd1,alpha1,alpha2,alphau,alphal
      real(dp) :: gg,ggold,gamma
      real(dp) :: munrl(3,nd), latgf(3,3,nd), testf, linef
      real(dp) :: mdisp(3,nd), dispor(nd), tmpdis, meadj
      real(dp) :: mcad(3,nd)
!
      integer usrexit,status, stat
      integer i,j,flag
      integer icom, nsrch, morig, mfcomp, mczr
!
      character*80 filename, gfoutput
      character*2 mcele
!
!     We need the coordinates of atoms in a block with no displacement
!     field applied in order to properly evaluate displacements. This
!     is read from unit 94.
!
      call fndrec(94,'D',stat)
      read(94,*) (mcad(1,i), mcad(2,i), mcad(3,i), mcele, mczr, i=1,nd)
!
      do i = 1,nd
         do j = 1,3
            mcad(j,i)=mcad(j,i)*lena(j)
         enddo
      enddo
!
!     Direction in which the line force is to be applied
!
      call fndrec(8,'MFCOMP',stat)
      read(8,*) mfcomp
!
!     Read in the line force. This is in eV / A^2 so we need to covert
!     it into eV / A for use in the code. TESTF = LENA(3)*LINEF
!
      call fndrec(8,'LINEF',stat)
      read(8,*) linef
!
      testf = lena(3) * linef
!
!     First read in unrelaxed coordinates of active atoms and store in
!     order to calculate the displacements. From these we also need to
!     find the atom at the origin.
!
      if (iter == 1) then
         do i=1,nd,1
            do j=1,3
               munrl(j,i)=ad(j,i)
            enddo
            if(abs(mcad(1,i)) < 1.0e-5_dp.and.abs(mcad(2,i)) < 1.0e-5_dp & 
     &      .and. (abs(mcad(3,i)) < 1.0e-5_dp)) morig = i
         enddo
      endif
      print *, "MORIG=", morig
      print *, "TESTF=", testf
      print *, "MFCOMP=", mfcomp

!
!    Initialize the Hessian matrix.
!

      if (rlxflg == 1) then
         if (hessmx < 3*nd) then
            write(6,'(''Increase the size of HESSMX to at least '', & 
     &                I5)') 3*nd
            call panic()
         endif
         do i = 1,3*nd,1
            do j = 1,3*nd,1
               hess(j,i) = 0.0_dp
            enddo
            hess(i,i) = 1.0_dp
         enddo
      endif
      if (rlxflg == 3) then
         ggold = 0.0_dp
      endif

!
!    Relax structure.
!

      icom = 0
      status = 0

      flag = 1
      call getetot(flag)
      call erasab()

      ein = eprom + ebond + epair - eatom - uent
      eout = ein
      u0   = eout + 0.5_dp*uent
      ered = 1.0e30_dp

!
!     Applying a force of TESTF eV/A to atom at origin
!
      ftot(mfcomp,morig) = ftot(mfcomp,morig)+testf
!
      maxf = 0.0_dp
      do i = 1,nd,1
         magf = sqrt(ftot(1,i)**2+ftot(2,i)**2+ftot(3,i)**2)
         if (magf > maxf) maxf = magf
      enddo

      if (rlxflg == 1) then
         do i = 1,nd,1
            newftot(1,i) = -ftot(1,i)
            newftot(2,i) = -ftot(2,i)
            newftot(3,i) = -ftot(3,i)
         enddo
      endif

!
!     This is the main loop where the variable metric relaxation is
!     performed.
!
      do while ((iter <= mxiter).and.(maxf > ftol).and. & 
     &          (icom /= 1).and.(usrexit() == 0).and. & 
     &          (step > stpmin).and.(status == 0))
!
         write(6,'(/''Iteration #             '',I5)') iter
         write(6,'(''Free energy, U(T=0) (start) = '', & 
     &             2G22.12)') ein,u0
         write(6,'(''Fermi energy (start)        = '',G12.5)') lef
         write(6,'(''Number of electrons (start) = '',G12.5)') totnia
         write(6,'(''Maximum force (start)       = '',G12.5)') maxf
!
!     Variable metric relaxation
!     
         if (rlxflg == 1) then
            call mstmin(3*nd,3*mxnd,hess,newftot,ftot,cg,newad, & 
     &                  step,0.0_dp,alpha,icom, & 
     &                  gd1,alpha1,alpha2,alphau,alphal,nsrch,6)
            if (dndm > 1.0e-6_dp) then
               dmu = (locc-totnia)/dndm
               call move(ad,nd,cg,alpha,dmdl,dmu)
               lef = lef + dmu
            else
               call move(ad,nd,cg,alpha,dmdl,dmu)
            endif
         endif
!
!     Conjugate gradient relaxation
!
         if (rlxflg == 3) then
            if (dndm > 1.0e-6_dp) lef = lef + (locc-totnia)/dndm
            gg = 0.0_dp
            do i = 1,nd,1
               gg = gg + ftot(1,i)**2+ftot(2,i)**2+ftot(3,i)**2
            enddo
            if ((abs(ggold) < 1.0e-6_dp).or. & 
     &          (mod(iter,5) == 0)) then
               gamma = 0.0_dp
            else
               gamma = gg/ggold
            endif
            write(6,'(''GAMMA = '',G12.5)') gamma
            ggold = gg
 1          do i = 1,nd,1
               cg(1,i) = gamma*cg(1,i) + ftot(1,i)
               cg(2,i) = gamma*cg(2,i) + ftot(2,i)
               cg(3,i) = gamma*cg(3,i) + ftot(3,i)
               newad(1,i) = ad(1,i)
               newad(2,i) = ad(2,i)
               newad(3,i) = ad(3,i)
            enddo
            call lgfsafemin(newad,cg,ein,eout,status,testf,morig,mfcomp, & 
     &           mcad)
            if ((status == 1).and.(abs(ggold) > 1.0e-6_dp)) then
               ggold = 0.0_dp
               gamma = 0.0_dp
               goto 1
            endif
         endif        

         if (rlxflg == 1) then
            flag = 2
            call getetot(flag)
            eout = eprom + ebond + epair - eatom - uent
!
!     Adjusting total energy due to application of the test force
!
            tmpdis = ad(mfcomp,morig)-mcad(mfcomp,morig)
            meadj = -1.0_dp*testf*tmpdis
            eout = eout + meadj
            ered = ein-eout
            ftot(mfcomp,morig) = ftot(mfcomp,morig)+testf            
!
            if (ered < 0.0_dp) then
               write(6,'(/''Energy increasing => Reducing step.'')')
               write(6,'(''Old step = '',G12.5)') step
!               STEP = STEP/2.0D0
               step = 0.75_dp * step
               write(6,'(''New step = '',G12.5)') step
            else
               write(6,'(/''Energy decreasing => Increasing step.'')')
               write(6,'(''Old step = '',G12.5)') step
               step = 1.1_dp*step
               write(6,'(''New step = '',G12.5)') step
            endif
            if (abs(locc-totnia)/real(nd, dp) > 0.001_dp) then
               flag = 1
               call getetot(flag)
               ftot(mfcomp,morig) = ftot(mfcomp,morig)+testf
            endif
            do i = 1,nd
               ftot(1,i) = -ftot(1,i)
               ftot(2,i) = -ftot(2,i)
               ftot(3,i) = -ftot(3,i)
            enddo
            eout = eprom + ebond + epair - eatom - uent + meadj
         endif

         u0   = eout  + 0.5_dp*uent

         maxf = 0.0_dp
         do i = 1,nd,1
            magf = sqrt(ftot(1,i)**2+ftot(2,i)**2+ftot(3,i)**2)
            if (magf > maxf) maxf = magf
         enddo

         write(6,'(''Free energy, U(T=0) (finish) = '',2G22.12)')eout,u0
         write(6,'(''Fermi energy (finish)        = '',G12.5)') lef
         write(6,'(''Number of electrons (finish) = '',G12.5)') totnia
         write(6,'(''Maximum force (finish)       = '',G12.5)') maxf
         ein = eout
         iter = iter + 1

         if (autosave > 0) then
            if (mod(iter,autosave) == 0) call dump()
         endif

      enddo
      
      if (rlxflg == 1) then
         do i = 1,nd,1
            ftot(1,i) = -ftot(1,i)
            ftot(2,i) = -ftot(2,i)
            ftot(3,i) = -ftot(3,i)
         enddo
      endif

!
!     Determining Gij(r) when the relaxation has finished
!
      do i=1,nd
         do j=1,3
            mdisp(j,i)=ad(j,i)-mcad(j,i)
         enddo
!
!     Adding condition for whether atoms jump from top of block to bottom...
!
         if (mdisp(3,i) > 0.5_dp*lena(3)) mdisp(3,i)=mdisp(3,i) & 
     &        -lena(3)
         if (mdisp(3,i) < -0.5_dp*lena(3)) mdisp(3,i)=mdisp(3,i) & 
     &        +lena(3)
      enddo
!
      do i = 1,nd
         do j = 1,3
            if (abs(mcad(j,i)) < 1.0e-5_dp) mcad(j,i) = 0.0_dp
         enddo
      enddo
!
      do i=1,nd
         dispor(i)=sqrt(((mcad(1,i)-mcad(1,morig))**2.0_dp)+ & 
     &        ((mcad(2,i)-mcad(2,morig))**2.0_dp)+  & 
     &        ((mcad(3,i)-mcad(3,morig))**2.0_dp))
      enddo
!     
      do i=1,nd
         do j=1,3
!     
            latgf(j,mfcomp,i) = mdisp(j,i)/linef
!     
         enddo
      enddo
!     
      write(gfoutput,'("latGFcomp.",I1)') mfcomp
      open (unit=92, status='UNKNOWN', file=gfoutput)
      open (unit=93, status='UNKNOWN', file='displGF.out')
      do i=1,nd
!
!     Printing lattice greens functions for n.n. up to the distance
!     given in the condition below
!
         if (dispor(i)  <  5.70_dp) then
            write(92,*) "1"
            write(92,20) (mcad(j,i),j=1,3)
            write(92,20) (latgf(j,mfcomp,i),j=1,3)
         endif
!
         write(93,*) i, dispor(i)
         write(93,*) "DISP.", (mdisp(j,i), j=1,3)
         write(93,*) "FTOT.", (ftot(j,i), j=1,3)
         write(93,*) "*********************************************"
      enddo
 20   format(3f15.7)
!
      call diffdisp(ad,munrl,mcad,dispor,nd,mfcomp,lena,ftot)
!
      end
