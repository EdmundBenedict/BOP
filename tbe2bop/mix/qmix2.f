      subroutine qmix2(lc,nbas,nsp,nlmq,ltb,dclabl,ipc,it,itmax,cnvg,
     .                 cnvm,mmix,nkill,nelts,beta,betav,tm,tj1,tj2,qmpol,
     .                 qits,mmom,mmom0,mits,a1,a2,rms,wc,broy,wt,magmix,rms2)
! C- Mixing multpole moments for TB-L (l-independent U)
! C ----------------------------------------------------------------------
! Ci Inputs:
! Ci   lc: switch, if T mix only l=0 monopole moments
! Ci   nbas,nsp,nlmq,ltb,it,itmax,cnvg,cnvm,mmix,nkill,beta,betav,tm
! Ci   beta and betav mix multipoles and magnetic moments respectively
! Ci   qmpol: multipole moments from tbmpol
! Ci   mmom: magnetic moments by site
! Ci   qits, mits, tj1, tj2, a1, a2: work arrays
! Co Outputs:
! Co   qmpol (and mmom) mixed, rms
! Cr Remarks
! Cr   If U is independent of l (but site dependent), i.e., Uav=T then
! Cr   the most efficient is to mix the Q_RL and magnetic moments.
! Cr
! Cr   Mixing works like this:
! Cr   At some iteration, i, we have a density (represented here simply
! Cr   by multipoles, rather like the ASA) rho_i. This makes a hamiltonian
! Cr   and the Schrodinger equation (bndtb, tbfrce) acts like a black box
! Cr   to produce a function f(rho_i). The Anderson mixing consults all
! Cr   previous rho_i and f(rho_i) and proposes a new density rho_i+1
! Cr   which is fed back into the black box until rho_i+i is equal to
! Cr   within an rms difference to rho_i.
! Cr
! Cr   qits is a work array which holds all previous rho_i, f(rho_i) and
! Cr   rho_i+1, keeping just the previous mmix iterations. At each entry
! Cr   to qmix these are rolled back in the array and iteration mmix+1 is
! Cr   discarded so that after rolling the array is structured like this:
! Cr
! Cr                  i=0     i=1    2    3   ...
! Cr   qits(..,i,1)  empty*   rho_i from previous iterations**
! Cr   amix notation          x_0   x_1  x_2 ...
! Cr
! Cr    * mixed Q_RL will be put in here after amix
! Cr   ** note i=1 is most recent, i=2 next most recent, etc.
! Cr      if it=1 then for i=1 these are zero or Q_RL from disc (restart)
! Cr
! Cr   qits(..,i,2)  i=0     i=1    2    3   ...
! Cr                 empty*   f(rho_i) from previous iterations**
! Cr   amix notation f(x_0) f(x_1) f(x_2) ...
! Cr
! Cr    * the most recent rho_i from tbmpol is copied into here before
! Cr      mixing
! Cr   ** note i=1 is most recent, i=2 next most recent, etc.
! Cr
! Cr   amix requires as input f(x_i) and d(x_i) = f(x_i) - x_i
! Cr   hence the call to daxpy in constructing the work array, a.
! Cr
! Cr   If nsp=2, then qits(1,.) holds the charge, the magnetic moments
! Cr   are maintained in mits. These are mixed in separate calls to amix.
! Cr   For notes on Broyden mixing see qmix
! Cr   Included a new form of "clock mixing", when wa=1 the procedure for
! Cr   both broy and anders is to mix the charge until it's self consistent
! Cr   then mix the magnetic moment for one iteration, and repeat this procedure
! Cr   until the magnetic moment becomes self consistent.
! Cu Updates
! Cu   qmix2 mixes the charges and the magnetic moments separately, which
! Cu   allows different beta for each.
! Cu   Broyden mixing uses wt(1) and (2) to scale the charge and magnetic 
! Cu   moment respectively, wt(3) is unused and setting either to zero 
! Cu   freezes that moment.
! Cu   Added "clock mixing", where the charge mixes to sc before each iteration
! cu   of mag moment mixing, this is enabled by setting wa=1, a seperate cnvg
! Cu   parameter may be specified for the magnetic moments with 'elind=', 
! Cu   if this left at zero it will use cnvg.
! Cb Bugs
! Cb   No doubt keeping both qits and a is redundant, maybe someone clever
! Cb   can save memory by getting rid of qits and building "a" directly
! Cb   from qmpol.
! C ----------------------------------------------------------------------
      implicit none
! C Passed Parameters
      logical lc
      integer nbas,nsp,nlmq,ltb,it,itmax,mmix,nkill,ipc(nbas),nelts,magmix,broy
      double precision qmpol(nlmq,nbas),qits(nlmq,nbas,0:mmix+1,2),
     .                 dclabl(1),
     .                 mmom(nbas),mmom0(nbas),mits(nbas,0:mmix+1,2),
     .                 a1(nelts,0:mmix+1,2),a2(nbas,0:mmix+1,2),wt(3)
      double precision cnvg,cnvm,beta,betav,tm,tj1(*),tj2(*),rms,wc
! C Local Variables
      integer neltsm,nmix1,nmix2,amix,npmix,i,ipr,iprint,
     .        ido,ib,ic,ierr,jump,i1mach,mgmix,crgmix
! !       integer onorm,okpvt
      real(8), allocatable :: norm(:)
      integer, allocatable :: kpvt(:)
      double precision rms1,rms2,b,qmp(nlmq),dqtot,d1mach,dabs,dsum,wctrue
      character clabl*8, outs*20
      logical T, F, IO, kill, bittst, cmdopt, ChgA, ChgB, SpinA, SpinB, ClockB
! C For MPI
      integer procid,master,mpipid
      logical mlog
! C Local iteration count
      integer LOCIT,jmix,j,magit,crgit,chrg
      save LOCIT,magit,crgit,rms1,ChgA,ChgB,SpinA,SpinB
      data  T, F / .true., .false. /

!       procid = mpipid(1)
!       master = 0
!       mlog = cmdopt('--mlog',6,0,outs)

!       if (iprint() .gt. 10) then
!         call awrit1(' QMIX2: rms=%d',' ',120,i1mach(2),rms)
!       endif

!       IO = bittst(ltb,2**16)
      jump = 1
!       if (lc) then
!         jump = nlmq
!       endif
!       neltsm = 0
!       if (nsp .eq. 2) then
!         neltsm = nbas
!       endif
     
! C --- "zero'th iteration" ---
!       if (it .eq. 0) then
!         ierr = -1
!         if (IO) then
!           if (procid == master) then
!             call ioqm(nsp,nlmq*nbas,neltsm,qmpol,mmom,ierr)
!           endif
!           call mpibc1(ierr,1,2,mlog,'qmix2','ierr')
!           if (ierr .ne. -1) then
!             call mpibc1(qmpol,nlmq*nbas,4,mlog,'qmix2','qmpol')
!             call mpibc1(mmom,neltsm,4,mlog,'qmix2','mmom')
!           endif
!         endif
!         if (ierr .eq. -1) then
!           call dcopy(nlmq*nbas,0d0,0,qmpol,1)
!           if (nsp .eq. 2) then
!             call dcopy(nbas,mmom0,1,mmom,1)
!           endif
!         endif
!         if (nsp .ne. 2 .and. wt(3) == 1) then
!           call rx('qmix2: cannot clock mix without magnetic moments')
!         endif
!         rms1 = 0
!         ChgA = T
!         ChgB = T
!         SpinA = T
!         SpinB = T
!         return
!       endif
      
! C --- kill "mix files" according to nkill or rms (see parmxp) ---
!       kill = .true.
!       if (cmdopt('--nomixkill',11,0,outs)) then
!         kill = .false.
!       endif

      if (it .eq. 1) then
        LOCIT = 0
        magit = 1
        crgit = 0
      endif
      if ( ( nkill .lt. 0 .or. rms .lt. 0d0 .or.
     .     ( nkill .gt. 0 .and. mod(it,nkill) .eq. 0 ) )
     .    .and. kill ) then
        LOCIT = 1
      else
        LOCIT = LOCIT + 1
      endif
      npmix = min(LOCIT-1,mmix)

! C --- check for charge neutrality ---
      dqtot = dsum(nbas,qmpol,nlmq)
      if (dabs(dqtot) .gt. d1mach(3)) then
        if (iprint() .gt. 40 .or. (iprint() .ge. 30 .and.
     .      dabs(dqtot) .gt. 1d-4)) then
          call awrit1(' QMIX2: input qtot=%;2e',' ',120,i1mach(2),dqtot)
        endif
      endif
      
!       if ((wt(3) == 1) .and. (it <= 2)) then
!           chrg = 1
!           rms = 1
!       endif
        
! C --- Roll back previous iterations ---
!       if ((wt(3) == 0) .or. it == 1 ) then
!         if (beta > 0 .or. it ==1 ) then
!           do  i = mmix, 0, -1
!             call dcopy(nelts,qits(1,1,i,1),jump,qits(1,1,i+1,1),jump)
!             call dcopy(nelts,qits(1,1,i,2),jump,qits(1,1,i+1,2),jump)
!           enddo
!         endif
!         if (betav > 0 .or. it == 1)  then
!           do  i = mmix, 0, -1
!             call dcopy(neltsm,mits(1,i,1),1,mits(1,i+1,1),1)
!             call dcopy(neltsm,mits(1,i,2),1,mits(1,i+1,2),1)
!           enddo
!         endif
!       elseif ((broy == 1) .and. (wt(3) == 1) .and. (rms > cnvg) ) then


        do  i = mmix, 0, -1
          call dcopy(nelts,qits(1,1,i,1),jump,qits(1,1,i+1,1),jump)
          call dcopy(nelts,qits(1,1,i,2),jump,qits(1,1,i+1,2),jump)
        enddo
        
        
        
!       elseif ((broy == 1) .and. (wt(3) == 1) .and. (rms < cnvg) ) then
!         do  i = mmix, 0, -1
!           call dcopy(neltsm,mits(1,i,1),1,mits(1,i+1,1),1)
!           call dcopy(neltsm,mits(1,i,2),1,mits(1,i+1,2),1)
!         enddo
!       endif
       
! C --- reset flags ---
!       ChgA = F
!       ChgB = F
!       SpinA = F
!       SpinB = F
!       ClockB = F

! C --- Copy new Q_RL and mmom for this iteration from input qmpol, mmom
!       if (it .eq. 1) then
!         ierr = -1
!         if (IO) then
!           call pshprt(0)
!           if (procid. eq. master) then
!             call ioqm(nsp,nlmq*nbas,neltsm,qits(1,1,1,1),mits(1,1,1),
!      .                ierr)
!           endif
!           call mpibc1(ierr,1,2,mlog,'qmix2','ierr')
!           if (ierr .ne. -1) then
!             call mpibc1(qits(1,1,1,1),nlmq*nbas,4,mlog,'qmix2','qits')
!             call mpibc1(mits(1,1,1),neltsm,4,mlog,'qmix2','mits')
!           endif
!           call popprt
!         endif
!         if (ierr .eq. -1) then
!           call dcopy(nlmq*nbas,0d0,0,qits(1,1,1,1),1)
!           call dcopy(neltsm,mmom0,1,mits(1,1,1),1)
!         endif
!       endif
      
      
      if ((wt(3) == 0) .or. (it == 1)) then
        if (beta > 0 .or. it == 1) then
          call dcopy(nlmq*nbas,qmpol,1,qits(1,1,0,2),1)
        endif
!         if (betav > 0 .or. it == 1) then
!           call dcopy(neltsm,mmom,1,mits(1,0,2),1)
!         endif
      elseif ((broy == 1) .and.(wt(3)==1) .and. (rms > cnvg)) then 
        call dcopy(nlmq*nbas,qmpol,1,qits(1,1,0,2),1)
!       elseif ((broy == 1) .and.(wt(3)==1) .and. (rms < cnvg)) then 
!         call dcopy(neltsm,mmom,1,mits(1,0,2),1)
      endif
          
! C --- Build work arrays for amix ---
      if ( it <= 1 ) then
        if (beta > 0 .and. betav > 0 .or. it == 1) then
          do  i = 0, npmix
            call dcopy(nelts,qits(1,1,i+1,1),jump,a1(1,i,2),1)
            call dcopy(nelts,qits(1,1,i,2),jump,a1(1,i,1),1)
          enddo
        endif      
      elseif ((broy == 1) .and. (wt(3)==1) .and. (rms > cnvg)) then 
        chrg=1
        do  i = 0, npmix
          call dcopy(nelts,qits(1,1,i+1,1),jump,a1(1,i,2),1)
          call dcopy(nelts,qits(1,1,i,2),jump,a1(1,i,1),1)
!           call dcopy(neltsm,mits(1,i,1),1,a2(1,i,2),1)
!           call dcopy(neltsm,mits(1,i,2),1,a2(1,i,1),1)
        enddo
      endif
            
! C --- Mix ---
      if (broy == 1) then
        if (wt(3) == 1) then
          if (it <= 1) then
            rms2 = 0
          endif
! C --- Clock mix charge to self consistency and one mix of magnetic moments
!           if (rms > cnvg .or. it <= 2) then
            crgmix = min(crgit,mmix)
            b = beta
            wc = wt(1)
            call qmixb(nelts,crgmix,mmix,b,wc,rms1,a1,wctrue,nmix1)
            nmix2 = 0
            rms = rms1
            crgit = crgit +1
!             ChgB = T
!           elseif (rms < cnvg) then
!             mgmix = min(magit,mmix)
!             magmix = 1  
!             wc = wt(2)
!             b = betav
!             call qmixb(neltsm,mgmix,mmix,b,wc,rms2,a2,wctrue,nmix2)
!             magit = magit + 1
!             nmix1 = 0 
!             rms = rms2
!             SpinB = T
!           endif
!           ClockB = T
        endif
      endif

! C --- Get new Q_RL (and mmom) from work array ---
! C      if (wt(3) == 0) then
! C        if (ChgA .or. ChgB) then
! C          call dcopy(nelts,a1(1,0,2),1,qits(1,1,0,1),jump)
! C          call dcopy(nlmq*nbas,0d0,0,qmpol,1)
! C          call dcopy(nelts,a1(1,0,2),1,qmpol,jump)
! C        endif
! C        if (SpinA .or. SpinB) then
! C          call dcopy(neltsm,a2(1,0,2),1,mits(1,0,1),1)
! C          call dcopy(neltsm,a2(1,0,2),1,mmom,1)
! C        endif
! C      elseif (wt(3) == 1) then
       call dcopy(nelts,a1(1,0,2),1,qits(1,1,0,1),jump)
       call dcopy(nlmq*nbas,0d0,0,qmpol,1)
       call dcopy(nelts,a1(1,0,2),1,qmpol,jump)
!        call dcopy(neltsm,a2(1,0,2),1,mits(1,0,1),1)
!        call dcopy(neltsm,a2(1,0,2),1,mmom,1)
! C      endif 

! C --- check for charge neutrality ---
      dqtot = dsum(nbas,qmpol,nlmq)
      if (dabs(dqtot) .gt. d1mach(3)) then
        if (iprint() .gt. 40 .or. (iprint() .ge. 30 .and.
     .      dabs(dqtot) .gt. 1d-4)) then
          call awrit1(' QMIX2: adding q=%;2e to conserve charge',' ',
     .      120,i1mach(2),-dqtot)
        endif
        dqtot = dqtot / nbas
        call daxpy(nbas,dqtot,-1d0,0,qmpol,nlmq)
      endif
        
! C --- write moments to disc ---
!       if (procid .eq. master) then
!         if (iprint() .ge. 30) then
!           print *, ' '
!           print *,'QMIX2: writing moments to disc..'
!         endif
!         ierr = 1
!         call ioqm(nsp,nlmq*nbas,neltsm,qmpol,mmom,ierr)
!       endif

! C --- Printout ---
!       if (iprint() .lt. 10) return
!       print 100
! 
!       if (ChgA .and. SpinA) then
!         call awrit7(
!      .  ' Anderson iteration %,4i: %,4i charges; mixed %,2i of %,2i, beta=%,4d '
!      .  //'rms diff: %g (tol: %g)',' ',128,i1mach(2),it,nelts,nmix1,npmix,beta,
!      .  rms1,cnvg)
!         if (nmix1 .gt. 0) write (*,110) (tj1(i),i=1,nmix1)
!         call awrit6('                          %,4i spins;'//
!      .                  '   mixed %,2i of %,2i, beta=%,4d '
!      .    //'rms diff: %g (tol: %g)',' ',128,i1mach(2),neltsm,nmix2,npmix,betav,rms2,cnvm)
!       endif
!       if (ChgA .and. .not. SpinA) then
!         call awrit7(
!      .  ' Anderson iteration %,4i: %,4i charges; mixed %,2i of %,2i, beta=%,4d '
!      .  //'rms diff: %g (tol: %g)',' ',128,i1mach(2),it,nelts,nmix1,npmix,beta,
!      .  rms1,cnvg)
!         if (nmix1 .gt. 0) write (*,110) (tj1(i),i=1,nmix1)
!       endif
!       if (.not. ChgA .and. SpinA) then
!         call awrit7(
!      .  ' Anderson iteration %,4i: %,4i spins;   mixed %,2i of %,2i, beta=%,4d '
!      .  //'rms diff: %g (tol: %g)',' ',128,i1mach(2),it,nelts,nmix2,npmix,betav,rms2,cnvm)
!         if (nmix2 .gt. 0) write (*,110) (tj2(i),i=1,nmix2)
!       endif        
!       if (ChgB .and. SpinB) then
!         call awrit6(
!      .  ' Broyden iteration %,4i: %,4i charges; mixed %,2i of %,2i, '
!      .  //'rms diff: %g (tol: %g)',' ',128,i1mach(2),it,nelts,nmix1,npmix,
!      .  rms1,cnvg)
!         call awrit5('                         %,4i spins;'//
!      .    '   mixed %,3i of %,3i, '
!      .          //'rms diff: %g (tol: %g)',' ',128,i1mach(2),neltsm,nmix2,npmix,rms2,cnvm)
!       endif
!       if (ChgB .and. .not. SpinB) then
!         call awrit6(
!      .  ' Broyden iteration %,4i: %,4i charges; mixed %,2i of %,2i, '
!      .  //'rms diff: %g (tol: %g)',' ',128,i1mach(2),it,nelts,nmix1,npmix,
!      .  rms1,cnvg)
!       endif
!       if (.not. ChgB .and. SpinB) then
!         call awrit6(
!      .        ' Broyden iteration %,4i: %,4i spins;   mixed %,2i of %,2i, '
!      .       //'rms diff: %g (tol: %g)',' ',128,i1mach(2),it,nelts,nmix2,npmix,rms2,cnvm)
!       endif        
!       if (iprint() .lt. 40) return
!       do  ib = 1, nbas
!         ic = ipc(ib)
!         call r8tos8(dclabl(ic),clabl)
!         call awrit1(' Atom %i '//clabl//'%cmultipole moments:',
!      .        ' ',180,i1mach(2),ib)
!         call dcopy(nlmq,qits(1,ib,1,1),1,qmp,1)
!         call awrit3('        Q(in) %d, %3:1d, %5:1d',' ',180,
!      .               i1mach(2),qmp,qmp(2),qmp(5))
!         call dcopy(nlmq,qits(1,ib,0,2),1,qmp,1)
!         call awrit3('       Q(out) %d, %3:1d, %5:1d',' ',180,
!      .               i1mach(2),qmp,qmp(2),qmp(5))
!         call dcopy(nlmq,qits(1,ib,0,1),1,qmp,1)
!         call awrit3('     Q(mixed) %d, %3:1d, %5:1d',' ',180,
!      .               i1mach(2),qmp,qmp(2),qmp(5))
!         if (nsp .eq. 2) then
!           call awrit3('    magnetic moment in %d, out %d, mixed %d',' ',
!      .                120,i1mach(2),mits(ib,1,1),mits(ib,0,2),mmom(ib))
!         endif
!       enddo
!   100 format(' QMIX2 mixing multipole moments:')
!   110 format(28x,' t_j :',10f8.4)
      end

 subroutine qmixb(neltst,npmix,mmix,beta,wc,rms,a,wctrue,jmix)
! C- Broyden mixing of a vector, Duane Johnson's approach
! C ------------------------------------------------------------------
! Ci  npmix: number of iterates available to mix
! Ci  a:    (*,i,1)  output values for prev. iteration i
! Ci        (*,i,2)  input  values for prev. iteration i
! Cio mmix: mmix > 0: number of iter to try and mix
! Ci        mmix < 0: use npmix instead of mmix.
! Co  npmix: (abs)  number of iter actually mixed.
! Co        (sign) <0, intended that caller update npmix for next call.
! Cr  Notations:
! Cr  x^(m): input vector for iteration m
! Cr  F^(m): difference between output and input vector in iteration m
! Cu Updates
! C ------------------------------------------------------------------
      implicit none
! C ... Passed parameters
      integer neltst,npmix,mmix
      double precision beta,rms,wctrue,a(neltst,0:mmix+1,2)
! C ... Dynamically allocated local arrays
      real(8), allocatable :: xmp1(:)
      real(8), allocatable :: dx(:)
      real(8), allocatable :: wk(:)
! C ... Local parameters
      double precision ddot,dval,wc
      integer im,km,i,iprint,imix,stdo,nglob,broyj
      integer, intent(out) :: jmix
      real(8) :: dxnorm, neltstsq

! C --- Setup ---

      allocate(xmp1(neltst))
      allocate(dx(neltst))
      neltstsq = sqrt(real(neltst,8))
! C ... imix is a local copy of npmix
      imix = npmix
! C --- Starting from iteration jmix, build the Jacobian matrix ---
   10 continue
      jmix = min(mmix,iabs(imix))
! C     call defdr(owk,neltst*2*(jmix+2))
      allocate(wk(neltst*2*(jmix+2)))
      do  km = 1, jmix
! C   ... this loops from most-distant to most-recent
        im = jmix-km+1
        call dcopy(neltst,a(1,im-1,1),1,dx,1)
        call daxpy(neltst,-1d0,a(1,im-1,2),1,dx,1)
        dxnorm = sqrt(ddot(neltst,dx,1,dx,1))
        rms = dxnorm/neltstsq
! C ---   Determine wc_true if wc < 0 ---
!         if (wc .lt. 0) then
          wctrue = -0.01_8*wc/dxnorm
          wctrue = min(max(wctrue,1d0),1d4)
!         else
!           wctrue = wc
!         endif
        if (km .eq. 1) wctrue = .01d0

        i = iprint()
        if (km .ne. jmix) i = i-20
        i = broyj(neltst,a(1,im-1,2),dx,km,0,i,beta,0d0,0d0,0d0,wctrue,wk,neltst,xmp1)
      enddo

      deallocate(dx)

! ! C --- Check for interactive change of npmix ---
! !       im = imix
! ! C      if (iprint() .gt. 30) call query('redo, npmix=',2,imix)
! !       if (iabs(imix) .gt. mmix .and. imix .ne. im .and. iprint() .gt. 30)
! !      .  call awrit1(' (warning) only %i iter available',' ',80,mmix)
! !       if (im .ne. imix) then
! !         deallocate(wk)
! !         goto 10
! !       endif
      npmix = imix
! C ... If no prior iter allowed, give up on npmix
! C --- Save x^(m+2) into a(*,0,2) and exit ---
      if (npmix /= 0) call dcopy(neltst,xmp1,1,a(1,0,2),1)


      deallocate(xmp1)

      end
