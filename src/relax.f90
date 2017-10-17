   subroutine relax()
 
! This is MM's version with small changes. The older version is commented out by the end.

      
! C
! C    This is a subroutine to relax atomic configurations.
! C    RLXFLG = 1  ==> Variable metric relaxation
! C    RLXFLG = 2  ==> Steepest descent relaxation
! C    RLXFLG = 3  ==> Conjugate gradient relaxation
! C    RLXFLG = 4  ==> FIRE
! C    RLXFLG = 5  ==> Steepest descent relaxation (without SAFEMIN)
! C



      use mod_precision
      use mod_all_scalar
      use mod_const

      implicit none


      include "Include/Atom.array"
      include "Include/BEC.array"
      include "Include/Force.array"
      include "Include/Misc.array"
      include "Include/PosVel.array"
      include "Include/Relax.array"
      include "Include/ag.conn"
      

      real(dp) :: ein,eout,ered,u0
      real(dp) :: alpha,maxf,magf,dmu
      real(dp) :: gd1,alpha1,alpha2,alphau,alphal
      real(dp) :: gg,ggold,gamma,displ
      real(dp), parameter :: maxdad = 0.20_dp
      real(dp) :: dis( 3 )

      integer usrexit,status
      integer i,j,flag
      integer icom,nsrch,keep(1)

      character(len=80) :: filename

! C     FIRE
      real(dp) :: max_a,max_v,cur_a,cur_v
      real(dp) :: vf,vg_dot_vg,fg_dot_fg
      real(dp) :: mix,help,incfac,decfac
      real(dp) :: maxdr,maxdt,mix_in,mixdec,dot
      integer v_mix,minsteps,cut,cuts


      
      if ( iter  >  mxiter )  then
         if (.not. quiet) write( 9, '("Maximum number of iterations exceeded.")')
         return
      endif

      if (.not. quiet) then
          filename = 'rel.dat'
         if ( mvflag == 17 ) then            
            open(unit = 72,file = filename,position='append')
            write ( 72, '(/"# stress: ",f8.5)') totstr
         else
            open(unit = 72,file = filename,status = 'replace')
         endif

      end if
      
! C
! C    Initialize FIRE parameters.
! C
      max_a=0.0_dp
      max_v=0.0_dp
      maxdr=0.1_dp
      v_mix=1
      mix_in=0.1_dp
      minsteps=5
      incfac=1.1_dp
      decfac=0.5_dp
      cut=0
      cuts=0
      mix=mix_in
      mixdec=0.99_dp
! C      DT=0.05_dp

! C
! C    Initialize the Hessian matrix.
! C

      if (rlxflg == 1) then
         hess = 0.0_dp
      elseif (rlxflg == 3) then
         ggold = 0.0_dp
      endif


! C
! C    Relax structure.
! C

      icom = 0
      status = 0
      
      if (.not. quiet) then
         write(9,'(/"============================================")')
         if (rlxflg == 1) then
            write(9,'("Starting variable metric relaxation with:")')
         elseif (rlxflg == 2) then
            write(9,'("Starting steepest descent relaxation with:")')
         elseif (rlxflg == 3) then
            write(9,'("Starting conjugate gradient relaxation with:")')
         elseif (rlxflg == 4) then
            write(9,'("Starting FIRE relaxation with:")')
         elseif (rlxflg == 5) then
            write(9,'("starting steepest descent relaxation (no SAFEMIN) with:")')
         endif
         
         write(9,'("ftol   = ",g12.5)') ftol
         write(9,'("mxiter = ",i5)')    mxiter
         write(9,'("stpmin = ",g12.5)') stpmin

      end if
      
      
      flag = 1
      call getetot(flag)

      ein = eprom + ebond + epair + eenv - eatom - uent
      eout = ein
      u0   = eout + 0.5_dp*uent
      ered = 1.0e+99_dp

      call constr(ftot,nd,cnst_a,cnst_v,cnst_n,mcnst_n)

      maxf = sqrt(maxval(sum(ftot(1:3,:nd)*ftot(1:3,:nd), dim = 1)))      

      if (.not. quiet) then
         write(filename, '("rel_",i6.6)') iter
         call writecfg('cfg/'//trim(filename)//'.cfg')
         call writecell('cell/'//trim(filename)//'.cell.in')
      
         if ((mvflag == 16 ).or.(mvflag == 17)) then
            write(filename,'("st.interm_",f5.3,".",i5.5)') totstr,iter
            call writeddp('ddplot/'//trim(filename)//'.dat')
! c         call plot( 0 )
         endif
      end if    

      if (rlxflg == 1) then
         newftot(:,:nd) = -ftot(:,:nd)         
      elseif (rlxflg == 3) then
         cg = 0.0_dp
      elseif (rlxflg == 4) then         
! c        vel(:,:nd) = ftot(:,:nd)*dt
         vel(:,:nd) = 0.0_dp
      endif

! C
! C     Main loop
! C
      do while ((iter < mxiter) .and. (maxf > ftol) .and. (icom /= 1) &
                  & .and. (usrexit() == 0) .and. (step > stpmin).and.(status == 0))

         if (.not. quiet) then         
            write(9,'(/"Iteration #             ",i5)') iter
            write(9,'("Free energy, U(T=0) (start) = ", 2(x,g22.12))') ein,u0
            write(9,'("Maximum force (start)       = ",g12.5)') maxf
            write(9,'(f12.6)') ftot(3,11)
            write(72,'(i6,2f15.5)') iter,ein,maxf
            
            write(6,'(/"RELAX Iteration #             ",i5)') iter
            write(6,'("RELAX Free energy, U(T=0) (start) = ", 2(x,g22.12))') ein,u0
            write(6,'("RELAX Maximum force (start)       = ",g12.5)') maxf
            write(6,'("RELAX FTOT",f12.6)') ftot(3,11)
         end if 

         if (rlxflg == 1) then
! C     Variable metric (BFGS) relaxation 

            call mstmin(3*nd,3*nd,hess,newftot,ftot,cg,newad, step, 0.0_dp, &
                        & alpha,icom, gd1,alpha1,alpha2,alphau,alphal,nsrch,9)
            call move(ad,nd,cg,alpha,dmdl,dmu)

! C            WRITE(77,'("Iteration : ",I5)') ITER
! C            WRITE(77,'(6F12.6)') HESS(1:3*ND,1:3*ND)

            flag = 1
            call getetot(flag)
            eout = eprom + ebond + epair + eenv - eatom - uent

            ered = ein-eout
            if (ered < 0.0_dp) then
               if (.not. quiet) write(9,'(/"Energy increasing => Reducing step.")')
               if (.not. quiet) write(9,'("Old step = ",G12.5)') step
               step = step * 0.5_dp
               if (.not. quiet) write(9,'("New step = ",G12.5)') step
            else
               if (.not. quiet) write(9,'(/"Energy decreasing => Increasing step.")')
               if (.not. quiet) write(9,'("Old step = ",G12.5)') step
               step = 1.1_dp * step
               if (.not. quiet) write(9,'("New step = ",G12.5)') step
            endif
!             
! C            IF (ABS(LOCC-TOTNIA)/DBLE(ND) > 0.001_dp) THEN
! C               FLAG = 1
! c               CALL GETETOT(FLAG)
! C            ENDIF

            ftot(:,:nd) = -ftot(:,:nd)
            eout = eprom + ebond + epair + eenv - eatom - uent

         elseif (rlxflg == 2) then
! C     Steepest descent relaxation 

            cg(:,:nd) = ftot(:,:nd)
            newad(:,:nd) = ad(:,:nd)
            
            call safemin(newad,cg,ein,eout,status)


         elseif (rlxflg == 3) then
! C     Conjugate gradient relaxation 
            
            gg = sum(ftot(:,:nd)*ftot(:,:nd))

            if ((abs(ggold) < 1.0e-6_dp) .or. (mod(iter,5) == 0)) then
               gamma = 0.0_dp
            else
               gamma = gg/ggold
            endif
            
            if (.not. quiet) write(9,'("gamma = ",g12.5)') gamma
            ggold = gg
            
            do
               cg(:,:nd) = gamma * cg(:,:nd) + ftot(:,:nd) 
               newad(:,:nd) = ad(:,:nd)
               
               call safemin(newad,cg,ein,eout,status)

               if ((status /= 1) .or. abs(ggold) <= 1.0e-6_dp) exit
               
               ggold = 0.0_dp               
               gamma = 0.0_dp
               
            end do

         elseif (rlxflg == 4) then

! C     FIRE relaxation

            vel(:,:nd) = vel(:,:nd) + 0.5_dp * dt * ftot(:,:nd)

            vf = 0.0_dp
            vg_dot_vg = 0.0_dp
            fg_dot_fg = 0.0_dp
            max_v = 0.0_dp
            max_a = 0.0_dp

            do i=1,nd
               vf = vf + dot_product(vel(1:3,i), ftot(1:3,i))
               
               cur_v      = dot_product(vel(1:3,i), vel(1:3,i))
               vg_dot_vg  = vg_dot_vg + cur_v
               if (cur_v > max_v)  max_v = cur_v

               cur_a      = dot_product(ftot(1:3,i), ftot(1:3,i))
               fg_dot_fg  = fg_dot_fg + cur_a
               if (cur_a > max_a)  max_a = cur_a
            enddo
            
            max_v = sqrt(max_v)
            max_a = sqrt(max_a)
            

! C     Cut the velocities if the total power done by forces is negative.
            
            if( vf <= 0.0_dp ) then
               if (.not. quiet) write (9, *)  "VF: Setting all velocities to zero "
               
               vel = 0.0_dp
               
               cut   = iter
               dt    = dt*decfac
               mix   = mix_in
               cuts  = cuts + 1
            else
! C     Otherwise mix the velocity and force.
               if (v_mix /= 0) then
                  help = mix * sqrt(vg_dot_vg/fg_dot_fg)
                  vel(:,:nd) =  (1-mix)*vel(:,:nd) + help*ftot(:,:nd)
               endif

! C     Adjust the maximum achievable timestep
     
               maxdt = max_v/(2*max_a) +  sqrt((max_v*max_v)/(2*max_a*max_a) + maxdr/max_a)
! C              WRITE (9, *)  "VF = ", VF
! C               WRITE (9, *)  "Maximum time step set to: ", MAXDT

               if (maxdt > 0.2)  maxdt=0.2

               if( (iter - cut)  >  minsteps ) then
                  dt = min(dt * incfac, maxdt)
                  mix = mix * mixdec
               endif
            endif

            if (.not. quiet) write(9,'("vf=",f12.8,"  dt=",f12.8,"  maxdt=",f12.8)') vf,dt,maxdt

            vel(:,:nd) = vel(:,:nd)+ 0.5_dp*dt*ftot(:,:nd)
            
! C     Move the atoms

            newftot(:,:nd) = ftot(:,:nd)

            if (mvflag == 17) then               
               unshad(:,:nd) = unshad(:,:nd) + dt*vel(:,:nd)
            else
! c             ad(:,:nd) = ad(:,:nd) + dt*vel(:,:nd) + 0.5_dp*dt*dt*ftot(:,:nd)
               ad(:,:nd) = ad(:,:nd) + dt*vel(:,:nd)             
            endif

            flag = 1
            call getetot(flag)
            eout = eprom + ebond + epair + eenv - eatom - uent

!             vel(:,:nd) = vel(:,:nd) + 0.5_dp * dt * (ftot(:,:nd) + newftot(:,:nd))

         elseif (rlxflg == 5) then
! C     Steepest descent relaxation (without SAFEMIN)

            do i = 1, nd
               dis(1:3) = ftot(1:3,i)*step

               displ = sqrt(sum(dis*dis))
               
               if ( displ > maxdad ) dis = maxdad * dis / displ

               if (mvflag == 17) then
                  unshad(1:3, i ) = unshad(1:3, i) + dis
               else
                  ad(1:3,i) = ad(1:3,i) + dis 
               endif
            enddo

            flag = 1
            call getetot(flag)
            eout = eprom + ebond + epair + eenv - eatom - uent

         endif


         u0 = eout + 0.5_dp*uent

         call constr(ftot,nd,cnst_a,cnst_v,cnst_n,mcnst_n)

         
         keep = maxloc(sum(ftot(:,:nd)*ftot(:,:nd), dim=1))
         maxf = sqrt(sum(ftot(:,keep(1))*ftot(:,keep(1))))
         
         if (.not. quiet) write(9,'("Free energy, U(T=0) (finish) = ",2G22.12)') eout,u0
         if (.not. quiet) write(9,'("Maximum force (finish)       = ",G12.5)') maxf
         
         if (.not. quiet) write(6,'("RELAX Free energy, U(T=0) (finish) = ",2G22.12)') eout,u0
         if (.not. quiet) write(6,'("RELAX Maximum force (finish)       = ",G12.5)') maxf

! C         WRITE(3,'(I6,3F15.5)') ITER,EIN,EOUT,MAXF
! C         WRITE(3,'(I6,2F15.5,G12.5)') AD(1,1),AD(2,1),AD(3,1)
! C         WRITE(3,'(I6,2F15.5,G12.5)') AD(1,2),AD(2,2),AD(3,2)

         ein = eout
         iter = iter + 1

         if ((writeper > 0) .and. (mod(iter,writeper) == 0).and.(.not. quiet)) then

            call outblock( 0 )

            write(filename, '("rel_",i6.6)') iter
            call writecfg('cfg/'//trim(filename)//'.cfg')
            
            call writecell('cell/'//trim(filename)//'.cell.in')         

            if ((mvflag == 16 ).or.(mvflag == 17)) then
               write(filename,'("st.interm_",f5.3,".",i5.5)') totstr,iter
               call writeddp('ddplot/'//trim(filename)//'.dat')
! c               call plot( 0 )
            endif
         endif
         
         if (autosave > 0 .and. mod(iter,autosave) == 0) call dump()

      enddo

      if (rlxflg == 1) ftot(1:3,1:nd) = -ftot(1:3,1:nd)

      if (.not. quiet) then
         write(9,'("Stopping relaxation with:")')
         if ((icom == 1) .or. (maxf < ftol)) then
            write(9,'("Relaxation successfully completed.")')
            write(9,'("ICOM = 1")')
         elseif (iter > mxiter) then
            write(9,'("Maximum number of iterations exceeded.")')
         elseif (step < stpmin) then
            write(9,'("Minimum step size reached.")')
         elseif (status == 1) then
            write(9,'("Forces not consistent with energy.")')
         endif

         write(9,'("Maximum force = ",G12.5,"eV/A")') maxf
         write(9,'("============================================"/)')
         
         write(72,'(i6,2f15.5)') iter,ein,maxf

         if (writeper > 0) close(2)
         close(72)
      end if

! C  If this is stress application, put back the unstressed coordinates.

      if (mvflag == 17) then
         ad(1:3,1:nd) = unshad(1:3,1:nd)
         adinert(1:3,1:ninert) = unshind(1:3,1:ninert)
! c         write (*,'("relaxds:",f15.9)') ad(3,1)
      endif

      
      if (.not. quiet) then
         write(*,'(/10X,"CALCULATION OF TOTAL ENERGY"/)')
         write(*,'("   Etot  = ",g16.10 ,"(",f10.6,")"/)') eout , eout/nd
         write(*,'("   Ebond = ",g14.7,2x,"(",f10.6,")")') ebond, ebond/nd
         write(*,'("   Eprom = ",g14.7,2x,"(",f10.6,")")') eprom, eprom/nd
         write(*,'("   Epair = ",g14.7,2x,"(",f10.6,")")') epair, epair/nd
         write(*,'("   Eenv  = ",g14.7,2x,"(",f10.6,")")') eenv , eenv/nd
         write(*,'("   Eatom = ",g14.7,2x,"(",f10.6,")")') eatom, eatom/nd
         write(*,'("   Eent  = ",g14.7,2x,"(",f10.6,")")') uent , uent/nd

         write(*,'(/" Maximum force is ",E14.7," acting on atom ",I5/)') maxf, keep

         filename = genfile(1:lengfn)//'.eng'
         open(unit = 3,file = filename)

         write(3,'(/10X,"CALCULATION OF TOTAL ENERGY"/)')
         write(3,'("   Etot  = ",g16.10  ,"(",f10.6,")"/)') eout,eout/nd
         write(3,'("   Ebond = ",g14.7,2x,"(",f10.6,")")') ebond,ebond/nd
         write(3,'("   Eprom = ",g14.7,2x,"(",f10.6,")")') eprom,eprom/nd
         write(3,'("   Epair = ",g14.7,2x,"(",f10.6,")")') epair,epair/nd
         write(3,'("   Eenv  = ",g14.7,2x,"(",f10.6,")")') eenv,eenv/nd
         write(3,'("   Eatom = ",g14.7,2x,"(",f10.6,")")') eatom,eatom/nd
         write(3,'("   Eent  = ",g14.7,2x,"(",f10.6,")")') uent, uent/nd
         write(3,'(/" Maximum force is ",E14.7," acting on atom ",I5/)') maxf, keep

         close(3)
         
         if ( mvflag == 16 .or. mvflag == 17 ) then
            write(filename,'("plot_",f5.3)') totstr
            call writeddp('ddplot/'//trim(filename)//'.dat')
            call writecfg('cfg/'//trim(filename)//'.cfg')
! c         call plot( 1 )
         endif
         

!    C      IF ( MVFLAG == 17 ) THEN
!    C         WRITE(FILENAME,'("Cfg/stress_",F5.3,".dat")') TOTSTR
!    C         CALL WRITECFG(FILENAME)
!    C         CALL PLOT( 1 )
!    C      ENDIF

      end if

   end subroutine relax

      
      
      
      
      
      
      
      
!       
!       
!       
!       
!       
!       
!       
!       
!       
!       
!       
!       
!       
!       
!       
!       
!       
!       
! 
!       
!       
!       
! !       
! !       
! !       
! !       
! !       
! !       
! !       
! !       
! !             subroutine relax()
! !       
! !       
! !           use mod_precision
! ! 
! ! 
! !           use mod_all_scalar
! ! 
! !           use mod_const
! ! 
! ! !
! ! !    This is a subroutine to relax atomic configurations.
! ! !    RLXFLG = 1  ==> Variable metric relaxation
! ! !    RLXFLG = 2  ==> Steepest descent relaxation
! ! !    RLXFLG = 3  ==> Conjugate gradient relaxation
! ! !
! ! 
! !       implicit none
! ! 
! ! !
! ! !    Include constants.
! ! !
! ! 
! ! !      include "Include/ALL.const"
! ! 
! ! !
! ! !    Include scalars.
! ! !
! ! 
! ! !      include "Include/ALL.scalar"
! ! 
! ! !
! ! !    Include arrays.
! ! !
! ! 
! !       include "Include/Atom.array"
! !       include "Include/BEC.array"
! ! !      include "Include/BondOrder.array"
! !       include "Include/Force.array"
! ! !      include "Include/Hamilt.array"
! !       include "Include/Misc.array"
! ! !      include "Include/Moment.array"
! ! !      include "Include/NebList.array"
! ! !      include "Include/NotSRT.array"
! !       include "Include/PosVel.array"
! ! !      include "Include/RecurCoef.array"
! !       include "Include/Relax.array"
! ! !      include "Include/SRT.array"
! ! 
! ! !
! ! !    Declare the simple variables.
! ! !
! ! 
! !       real(dp) :: ein,eout,ered,u0
! !       real(dp) :: alpha,maxf,magf,dmu
! !       real(dp) :: gd1,alpha1,alpha2,alphau,alphal
! !       real(dp) :: gg,ggold,gamma
! ! 
! !       integer usrexit,status
! !       integer i,j,flag
! !       integer icom,nsrch
! ! 
! !       character*80 filename, filename2
! ! 
! !       if ( elplusrel  ==  1 ) then
! !          iter = 1
! !          step = 0.05_dp
! !       endif
! !          
! ! !
! ! !    Open the configuration file.
! ! !
! ! 
! !       if (writeper > 0) then
! !          filename = genfile(1:lengfn)//'.dat'
! !          open(unit = 2,file = filename,status = 'NEW')
! !          if (datform == 1) then
! !             write(2,'(''no_of_atoms ='',I6)') nd
! !             write(2,'(''no_of_frames ='',I4)') (mxiter-nequil)/writeper
! !             write(2,'(''data ='')')
! !          endif
! !       endif
! ! 
! ! !
! ! !    Initialize the Hessian matrix.
! ! !
! ! 
! !       if (rlxflg == 1) then
! !          if (hessmx < 3*nd) then
! !             write(6,'(''Increase the size of HESSMX to at least '', & 
! !      &                I5)') 3*nd
! !             call panic()
! !          endif
! !          do i = 1,3*nd,1
! !             do j = 1,3*nd,1
! !                hess(j,i) = 0.0_dp
! !             enddo
! !             hess(i,i) = 1.0_dp
! !          enddo
! !       elseif (rlxflg == 3) then
! !          ggold = 0.0_dp
! !       endif
! ! 
! ! !
! ! !    Relax structure.
! ! !
! ! 
! !       icom = 0
! !       status = 0
! ! 
! !       write(9,'(/''============================================'')')
! !       if (rlxflg == 1) then
! !          write(9,'(''Starting variable metric relaxation with:'')')
! !       elseif (rlxflg == 2) then
! !          write(9,'(''Starting steepest descent relaxation with:'')')
! !       elseif (rlxflg == 3) then
! !          write(9,'(''Starting conjugate gradient relaxation with:'')')
! !       endif
! !       write(9,'(''FTOL   = '',G12.5)') ftol
! !       write(9,'(''MXITER = '',I5)')    mxiter
! !       write(9,'(''STPMIN = '',G12.5)') stpmin
! ! 
! !       flag = 1
! !       call getetot(flag)
! !       call erasab()
! ! 
! !       ein = eprom + ebond + epair - eatom - uent
! !       eout = ein
! !       u0   = eout + 0.5_dp*uent
! !       ered = 1.0e30_dp
! ! 
! !       call constr(ftot,nd,cnst_a,cnst_v,cnst_n,mcnst_n)
! !       
! !       maxf = 0.0_dp
! !       do i = 1,nd,1
! !          magf = sqrt(ftot(1,i)**2+ftot(2,i)**2+ftot(3,i)**2)
! !          if (magf > maxf) maxf = magf
! !       enddo
! !       
! !       if (rlxflg == 1) then
! !          newftot(:,:nd) = -ftot(:,:nd)
! !       endif
! ! !       
! ! !      write (6,*) 'ITER', iter
! ! !        
! ! !      write (6,*)'MXITER',mxiter
! ! !         
! ! !      write (6,*)'MAXF',maxf
! ! !         
! ! !      write (6,*)'FTOL' ,ftol 
! ! !      write (6,*)'ICOM',icom
! ! !      write (6,*)'USREXIT()',usrexit()
! ! !      write (6,*)'STEP',step
! ! !      write (6,*)'STPMIN',stpmin
! ! !      write (6,*)'STATUS',status
! ! !       
! ! !       stop
! !       do while ((iter <= mxiter).and.(maxf > ftol).and. & 
! !      &          (icom /= 1).and.(usrexit() == 0).and. & 
! !      &          (step > stpmin).and.(status == 0))
! !          
! !          write(6,'(/''Iteration #             '',I5)') iter
! !          write(6,'(''Free energy, U(T=0) (start) = '', & 
! !      &             2G22.12)') ein,u0
! !          write(6,'(''Fermi energy (start)        = '',G12.5)') lef
! !          write(6,'(''Number of electrons (start) = '',G12.5)') totnia
! !          write(6,'(''Maximum force (start)       = '',G12.5)') maxf
! !         
! !          if (rlxflg == 1) then
! !             call mstmin(3*nd,3*mxnd,hess,newftot,ftot,cg,newad, & 
! !      &                  step,0.0_dp,alpha,icom, & 
! !      &                  gd1,alpha1,alpha2,alphau,alphal,nsrch,6)
! !             if (dndm > 1.0e-6_dp) then
! !                dmu = (locc-totnia)/dndm
! !                call move(ad,nd,cg,alpha,dmdl,dmu)
! !                lef = lef + dmu
! !             else
! !                call move(ad,nd,cg,alpha,dmdl,dmu)
! !             endif
! !          elseif (rlxflg == 2) then
! !             if (dndm > 1.0e-6_dp) lef = lef + (locc-totnia)/dndm
! !             do i = 1,nd,1
! !                cg(1,i) = ftot(1,i)
! !                cg(2,i) = ftot(2,i)
! !                cg(3,i) = ftot(3,i)
! !                newad(1,i) = ad(1,i)
! !                newad(2,i) = ad(2,i)
! !                newad(3,i) = ad(3,i)
! !             enddo
! !             call safemin(newad,cg,ein,eout,status)
! !          elseif (rlxflg == 3) then
! !             if (dndm > 1.0e-6_dp) lef = lef + (locc-totnia)/dndm
! !             gg = 0.0_dp
! !             do i = 1,nd,1
! !                gg = gg + ftot(1,i)**2+ftot(2,i)**2+ftot(3,i)**2
! !             enddo
! !             if ((abs(ggold) < 1.0e-6_dp).or. & 
! !      &          (mod(iter,5) == 0)) then
! !                gamma = 0.0_dp
! !             else
! !                gamma = gg/ggold
! !             endif
! !             write(6,'(''GAMMA = '',G12.5)') gamma
! !             ggold = gg
! !  1          do i = 1,nd,1
! !                cg(1,i) = gamma*cg(1,i) + ftot(1,i)
! !                cg(2,i) = gamma*cg(2,i) + ftot(2,i)
! !                cg(3,i) = gamma*cg(3,i) + ftot(3,i)
! !                newad(1,i) = ad(1,i)
! !                newad(2,i) = ad(2,i)
! !                newad(3,i) = ad(3,i)
! !             enddo
! !             call safemin(newad,cg,ein,eout,status)
! !             if ((status == 1).and.(abs(ggold) > 1.0e-6_dp)) then
! !                ggold = 0.0_dp
! !                gamma = 0.0_dp
! !                goto 1
! !             endif
! !          endif
! ! 
! !          if (rlxflg == 1) then
! !             flag = 2
! !             call getetot(flag)
! !             eout = eprom + ebond + epair - eatom - uent
! !             ered = ein-eout
! !             if (ered < 0.0_dp) then
! !                write(6,'(/''Energy increasing => Reducing step.'')')
! !                write(6,'(''Old step = '',G12.5)') step
! !                step = step/2.0_dp
! !                write(6,'(''New step = '',G12.5)') step
! !             else
! !                write(6,'(/''Energy decreasing => Increasing step.'')')
! !                write(6,'(''Old step = '',G12.5)') step
! !                step = 1.1_dp*step
! !                write(6,'(''New step = '',G12.5)') step
! !             endif
! !             if (abs(locc-totnia)/real(nd, dp) > 0.001_dp) then
! !                flag = 1
! !                call getetot(flag)
! !             endif
! !             ftot(:,:nd) = -ftot(:,:nd)
! !             eout = eprom + ebond + epair - eatom - uent
! !          endif
! !          
! !          u0   = eout  + 0.5_dp*uent
! ! 
! !          call constr(ftot,nd,cnst_a,cnst_v,cnst_n,mcnst_n)
! ! 
! !          maxf = 0.0_dp
! !          do i = 1,nd,1
! !             magf = sqrt(ftot(1,i)**2+ftot(2,i)**2+ftot(3,i)**2)
! !             if (magf > maxf) maxf = magf
! !          enddo
! ! 
! !          write(6,'(''Free energy, U(T=0) (finish) = '',2G22.12)')eout,u0
! !          write(6,'(''Fermi energy (finish)        = '',G12.5)') lef
! !          write(6,'(''Number of electrons (finish) = '',G12.5)') totnia
! !          write(6,'(''Maximum force (finish)       = '',G12.5)') maxf
! ! 
! !          ein = eout
! !          iter = iter + 1
! ! 
! !          if ((writeper > 0).and.(iter > nequil)) then
! !             if (mod(iter,writeper) == 0) then
! !                if (datform == 1) then
! !                   do i = 1,nd,1
! !                      write(2,'(3(e12.5,1X),4(e9.2,1X))') & 
! !      &                    ad(1,i),ad(2,i),ad(3,i), & 
! !      &                    avs(1,z(i)),avs(2,z(i)), & 
! !      &                    avs(3,z(i)),avs(4,z(i))
! !                   enddo
! !                else
! !                   write(2,'(I6)') nd
! !                   write(2,'(A)') ' '
! !                   do i = 1, nd, 1
! !                      write(2,'(A,3F12.5)') symb(i), & 
! !      &                    ad(1,i),ad(2,i),ad(3,i)
! !                   enddo
! !                endif
! !             endif
! !          endif
! ! 
! !          if (autosave > 0) then
! !             if (mod(iter,autosave) == 0) call dump()
! !          endif
! ! 
! !       enddo
! ! 
! !       if (rlxflg == 1) then
! !          do i = 1,nd,1
! !             ftot(1,i) = -ftot(1,i)
! !             ftot(2,i) = -ftot(2,i)
! !             ftot(3,i) = -ftot(3,i)
! !          enddo
! !       endif
! ! 
! !       write(9,'(''Stopping relaxation with:'')')
! !       if ((icom == 1).or.(maxf < ftol)) then
! !          write(9,'(''Relaxation successfully completed.'')')
! !       elseif (iter > mxiter) then
! !          write(9,'(''Maximum number of iterations exceeded.'')')
! !       elseif (step < stpmin) then
! !          write(9,'(''Minimum step size reached.'')')
! !       elseif (status == 1) then
! !          write(9,'(''Forces not consistent with energy.'')')
! !       endif
! !       write(9,'(''Maximum force = '',G12.5,''eV/A'')') maxf
! !       write(9,'(''============================================''/)')
! ! !
! !       filename2 = genfile(1:lengfn)//'.eng'
! !       open(unit = 3,file = filename2,status = 'UNKNOWN')
! ! !      WRITE (3, '(2F16.9)') LENA(1)*LENA(2)*LENA(3), EOUT
! !       write (3,'(2F16.9)') lena(3)/lena(1), eout
! !       close(3)
! ! !
! !       if (writeper > 0) close(2)
! ! !
! !       return
! ! !
! !       end
! ! 
