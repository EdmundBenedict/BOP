 
!****************************************************************************
!                                                                           *
!                                                                           *
!                             Alexey Girshick                               *
!              Department of Materials Science and Engineering              *
!                        University of Pennsylvania                         *
!                                                                           *
!                                                                           *
!****************************************************************************
!                                                                           *
!                                                                           *
!      This subroutine relaxes atomic configuration of the dislocation.     *
!                                                                           *
!                                                                           *
!****************************************************************************


   subroutine relax_ds( )
      use mod_precision
      use mod_all_scalar
      use mod_const
      use topologia
      
      implicit none

 
      include "Include/Atom.array"
      include "Include/PosVel.array"
      include "Include/Force.array"
      include "Include/Relax.array"
      include "Include/ag.conn"


      real(dp) :: eout, u0
      real(dp) :: maxf
      real(dp) :: ezero, surf

      integer usrexit
      integer i, j, k, flag
      integer stat
      character(len=80) :: filename


!       if (.not. quiet) call plot(0)
      
      
      flag = 1
      call getetot(flag)
      call erasab()
      eout = eprom + ebond + epair - eatom - uent
!      PRINT *, "EOUT=", EOUT
!      DO I=1,ND
!         PRINT *, J, I, FTOT(1,I), FTOT(2,I), FTOT(3,I)
!         PRINT *, "AD", AD(1,I)/3.839, AD(2,I)/3.839, AD(3,I)/3.839
!      ENDDO
      u0   = eout + 0.5_dp*uent
!       if (.not. quiet) write( 72, '("Iteration ",I4,".      Energy: ",F16.12)' )  iter-1, eout / nd

      call constr(ftot,nd,cnst_a,cnst_v,cnst_n,mcnst_n)

      call maxf_ds( maxf )
      
      if (.not. quiet) then
         write(72,*) iter-1, eout, maxf

         write(filename,'("st.interm_",f5.3,".",i5.5)') totstr,iter
         
         call writeddp('ddplot/'//trim(filename)//'.dat')
         call writexyz('xyz/'//trim(filename)//'.xyz', filename)         
         call writecfg('cfg/'//trim(filename)//'.cfg')
         
      end if

!  Iteration loop.
!       if (.not. quiet) print *, 'mxiter,ftol:' , mxiter, ftol

      do while ((iter <= mxiter) .and. (maxf > ftol) .and. (usrexit() == 0) )
       
!        Move the atoms.

         call move_ds( ) 

!        plot the intermediate configuration, if necessary.

         if ( mod( iter, writeper )  ==  0  .and.  (.not. quiet) ) then
            call outblock( 0 )
            
            write(filename,'("st.interm_",f5.3,".",i5.5)') totstr,iter
         
            call writeddp('ddplot/'//trim(filename)//'.dat')
            call writexyz('xyz/'//trim(filename)//'.xyz', filename)         
            call writecfg('cfg/'//trim(filename)//'.cfg')
            call writexbs('bs/'//trim(filename)//'.bs')            
         end if


!        Recalculate the energy.

!          if (dndm > 1.0e-6_dp) lef = lef + (locc-totnia)/dndm      
!          flag = 2

         call getetot(flag)

         eout = eprom + ebond + epair - eatom - uent
         u0   = eout  + 0.5_dp*uent
         if (.not. quiet) write( 72, '("Iteration ",I4,".      Energy: ",F16.12)' ) iter, eout / nd
         if (.not. quiet) call flush( 72 )

         call constr(ftot,nd,cnst_a,cnst_v,cnst_n,mcnst_n)

         call maxf_ds( maxf )
         if (.not. quiet) write(72,*) iter, eout, maxf

         iter = iter + 1

      enddo

      if (.not. quiet) then
         if     ( maxf  <  ftol )  then
            write( 72, '(''Relaxation successfully completed.'')')
         elseif ( iter  >  mxiter )  then
            write( 72, '(''Maximum number of iterations exceeded.'')')
         endif
      
         
         call writeddp('ddplot/relstrain.dat')
         call writexyz('xyz/relstrain.xyz', 'relstrain')
         call writecfg('cfg/relstrain.cfg')
         call writexbs('bs/relstrain.bs')
         call writecell('cell/relstrain.cell.in')
         
!      WRITE ( *, '("APPLIED STRESS = ",F8.4)') TOTSTR
!
!     Added by M. Cawkwell 11th September 2003. If the GFBC's are turned
!     on, we call subroutine GFBC_OUTP to write the file of coordinates
!     and forces needed by greenlat.c.
!
!      CALL FNDREC(8,'GFBCON',STAT)
!      READ(8,*) GFBCON
!
! C      CALL OUTBLOCK( 1 )
! C      WRITE (*,'(" ******* relaxds:",f15.9)') ad(3,1)

         if (mvflag  ==  16 .and. gfbcon == 1) then
            call gfbc_outp()
         endif
      
      end if

      end subroutine relax_ds

