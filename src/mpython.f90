 
      subroutine mpython()
          use mod_precision


          use mod_all_scalar

          use mod_const
!
!     This is a subroutine for Monte Carlo simulation (Monte Python, 
!     geddit??!). The Metropolis scheme has been implemented with few
!     embellishments. I am interested in using this code to attempt to
!     look at surface reconstructions, so it is really for finding 
!     optimised structures and not for getting thermodynamical data, 
!     although such things could easily be added later.
!
!     As we all know, Monte Carlo schemes are a bit fiddly and there 
!     are lots of adjustable parameters to help things proceed properly
!     so it will be written so that all of these things can be read in
!     from an external file.
!
!     M. Cawkwell, 11th July 2003.
!
      implicit none

!
!     Include constants.
!

!      include "Include/ALL.const"
      
!     
!     Include scalars.
!

!      include "Include/ALL.scalar"

!
!     Include arrays.
!

!      include "Include/Atom.array"
      include "Include/BEC.array"
!      include "Include/BondOrder.array"
!      include "Include/Force.array"
!      include "Include/Hamilt.array"
!      include "Include/Misc.array"
!      include "Include/Moment.array"
!      include "Include/NebList.array"
!      include "Include/NotSRT.array"
      include "Include/PosVel.array"
!      include "Include/RecurCoef.array"
!      include "Include/Relax.array"
!      include "Include/SRT.array"

!
!     Declare the simple variables.
!
      real(dp) :: howhot, kboltz, malpha, acprat
      real(dp) :: oldad(3,nd), pretot, nwetot, mcedif, mcprob
      real(dp) :: surad(3,100)
      parameter (kboltz = 8.617387467e-5_dp)
!
      integer selatm, mcstep, stat, rlxstp, mctot, flag, rlxcnt
      integer noaccp, ccount
      integer rndseed
      integer nsurad, sureqv(100)
      integer i, j
!      
      real rvec(1)
!
!     First we need to read in parameters from fort.8. These are
!     the temperaure, the number of MC steps for the simulation, the
!     acceptance ratio wanted (ALPHA will be adjusted during the
!     simulation to achieve this).
!
      call fndrec(8,'MCTOT',stat)
      read(8,*) mctot
      print *, "Total number of MC steps in simulation = ", mctot
!
      call fndrec(8,'ACPRAT',stat)
      read(8,*) acprat
      print *, "Acceptance ratio for simulation will be = ", acprat
!
      call fndrec(8,'HOWHOT',stat)
      read(8,*) howhot
      print *, "Temperature set to ", howhot,"K"
!
      call fndrec(8,'ALPHA',stat)
      read(8,*) malpha
      print *, "The starting ALPHA is = ", malpha
!
      call fndrec(8,'RLXSTP',stat)
      read(8,*) rlxstp
      print *, "Structure will be relaxed every ",rlxstp," interations"
!
!     Read in the seed for the random number generator
!
      call fndrec(8,'RNDSEED',stat)
      read(8,*) rndseed
      print *, "The seed for the random number generator = ", rndseed
!
!     Now we have these parameters, we can get started
!
      open(unit=97, status='UNKNOWN', file='mcprogress.dat')
!
      mcstep = 0
      rlxcnt = 0
      noaccp = 0
      ccount = 0
!
!     Initialize the random number generator
!
      call rluxgo(4,rndseed,0,0)
!
      do while (mcstep  <=  mctot)
!
!     Determining energy of starting structure
!
         flag = 1
         call getetot(flag)
         call erasab()
! 
         pretot = eprom + ebond + epair - eatom - uent
!
!     Counting the number of iterations of various things
!
 20      mcstep = mcstep + 1
         print *, "MCSTEP = ", mcstep
         rlxcnt = rlxcnt + 1
         ccount = ccount + 1
!
!     We only want to move the atoms in the surface layer using the
!     Monte Carlo method but we want several layers to relax and 
!     determine any movement of charge due to the lower coordination
!     at the surface. Therefore we will select the surface atoms, move
!     these and then put the new configuration back into the AD(3,ND)
!     array for calculating energies and forces. OK?
!
         nsurad = 0
         do i = 1,nd
            if (abs(ad(3,i)) < 1.0e-1_dp) then
               nsurad = nsurad + 1
               sureqv(nsurad) = i
               do j = 1,3
                  surad(j,nsurad) = ad(j,i)
               enddo
            endif
         enddo
!         PRINT *, "NSURAD =", NSURAD
!         DO I = 1,NSURAD
!            PRINT *, I, SUREQV(I)
!         ENDDO
!
!     Select atom to displace
!     
         call ranlux(rvec,1)
         selatm = int(rvec(1)*nsurad)
!
!     Saving the original configuration
!
         do i = 1,nd
            do j = 1,3
               oldad(j,i) = ad(j,i)
            enddo
         enddo
!
!     We only want to move the surface atoms in 2 dimensions (the
!     relaxation every step with take care of displacements normal
!     to the surface
         call ranlux(rvec,1)
         surad(1,selatm) = surad(1,selatm)  & 
     &        + ((2.0_dp*rvec(1))-1.0_dp)*malpha
         call ranlux(rvec,1)
         surad(2,selatm) = surad(2,selatm)  & 
     &        + ((2.0_dp*rvec(1))-1.0_dp)*malpha
!         CALL RANLUX(RVEC,1)
!         SURAD(3,SELATM) = SURAD(3,SELATM) 
!     +        + ((2.0D0*RVEC(1))-1.0D0)*MALPHA
!
!     Now we have moved one of the surface atoms, we can put the 
!     SURAD(3,NSURAD) array bach into the AD(3,ND) array
!
         do i = 1, nsurad
            do j = 1, 3
               ad(j,sureqv(i)) = surad(j,i)
            enddo
         enddo
!
!     Relaxing this new structure before calculating its energy
!
         if (rlxcnt  ==  rlxstp) then
            rlxcnt = 0
            iter = 0
            call fndrec(8,'STEP',stat)
            read(8,*) step
            call relax()
         endif
!     
!     Now getting energy of this new stucture
!
         flag = 1
         call getetot(flag)
         call erasab()
!     
         nwetot = eprom + ebond + epair - eatom - uent
!     
!     We now have to decide if the new configuration will be accepted
!     
         mcedif = nwetot - pretot
         mcprob = exp(-1.0_dp*mcedif / (kboltz*howhot))
!     
!     Accept the move if the energy is lower
!
         if ( mcedif  <=  0.0_dp ) then
            noaccp = noaccp + 1
            print *, "MOVE ACCEPTED"
            goto 10
         endif
!     
!     Accept the move if the energy increases but the probability
!     of the move is greater than a random number 0 -> 1
!
         call ranlux(rvec,1)
         if ( mcprob  >  rvec(1) ) then 
            noaccp = noaccp + 1
            print *, "MOVE ACCEPTED"
            goto 10
         endif
!
!     Otherwise we go back to the original configuration
!
         do i = 1,nd
            do j = 1,3
               ad(j,i) = oldad(j,i)
!
!     If the step wasn't accepted, we can cut out a few of the
!     prelimnaries - we already have the energy of this structure
!     (PRETOT) and we can just select another atom to move, so...
!
            enddo
         enddo
         print *, "MOVE REJECTED"
         goto 20
!
!     Now evaluating the acceptance ratio and adjusting ALPHA
!     We will adjust ALPHA to achieve the desired acceptance ratio
!     every 1000 MC steps
!
 10      if (ccount  ==  10) then
            if ( (noaccp/mcstep)  >  acprat ) malpha = malpha + 0.05_dp
            if ( (noaccp/mcstep)  <  acprat ) malpha = malpha - 0.05_dp
!
!     We will write out the coordinates of the atoms every 1000
!     interations too
!
            call outblock(0)
            ccount = 0
         endif
!
!     We will need a file to write the total energy vs. the number of
!     MC steps (fort.97 is available)
!
         write(97,100) mcstep, pretot
 100     format(1x,i7,1x,f12.7)
!
      enddo
!
      end
