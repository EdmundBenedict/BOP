
   subroutine bop(rtconf)
      use mod_precision
      use mod_all_scalar
      use mod_const
      use mod_conf
      use mod_ham, only : init_ham, free_ham
      use ab_io, only : init_ab
      use mod_clock
      use mod_chi, only : init_chi, free_chi
      use mod_kspace, only : kspace_destroy
      use mod_atom_ar
      use topologia, only : iproc, mpmap, master
      use mod_neb, only : nebf

!           use ssa
      implicit none

      type(run_conf_t), intent(inout), target :: rtconf


! This thing should be integrated better.
      include "Include/ag.conn"


      real(dp) :: curvature, slope
      integer(8) :: c1,c2
      integer :: i

      rtc => rtconf

      
      

!
!    Set up arrays, and read in parameters.
!
!    call system_clock(c1)
      call resetup( rtconf )
!    call system_clock(c2)
!    print *,'reset',real(c2-c1,dp)/real(cr,dp)
      if (rtconf % fs_only == 0 .and. rtconf % pot_flg == 2) then
         call init_ham(mpmap(iproc)+1, mpmap(iproc+1))
         call init_ab (mpmap(iproc)+1, mpmap(iproc+1), nsp)
         call init_chi(mpmap(iproc)+1, mpmap(iproc+1), nsp)
      end if


!
!    Find total number of electrons.
!

      call getocc(1,locc)


!
!    Carry out a simulation.
!

      if (mvflag == -1) then

!
!       Convergence test.
!
!       call system_clock(c1)
!          call convrg(rtconf)
         call getetot(1)
!       call system_clock(c2)
!       print *,'conv',real(c2-c1,dp)/real(cr,dp)
!          call outblock(1)
         call writecell('en-cell.out')

      elseif (mvflag == -5) then
! Dos plotting
         call getetot(1)
         select case (rtconf%pot_flg)
         case (2)
            call dosplot()
         case (4)
            call kdosplot()
         end select

      elseif (mvflag == -2) then

         call forcecheck(rtconf)

      elseif (mvflag == 0) then

!
!       Evaluate the energy and forces of a configuration.
!
      stop 'analyz temporarily removed until updated to magnetic. use mvflag -1'
!
!          call analyz()
!          call outfil(1)
!          call outfil(2)
!          call outfil(3)
!          call tag( 9 )


      elseif (mvflag == 1) then

!
!       Relax the atoms.
!

!          call outfil(1)
         call relax()
         iter = iter - 1
!          call outfil(2)
!          call outfil(3)
!          call dump()
!

      elseif ((mvflag >= 2).and.(mvflag <= 5)) then

!
!       Carry out molecular dynamics.
!

         call moldyn()
         call outfil(1)
         iter = iter - 1
         call outfil(3)
         call wrthis()
         call dump()

      elseif (mvflag == 6) then

!
!       Evaluate the elastic constants.
!

         call getelast( 0, curvature, slope )
         call outfil(1)
         call outfil(3)


      elseif (mvflag == 7) then

!
!       Perform a diffusion barrier calculation.
!

         call diffuse()


      elseif (mvflag == 9) then

!
!     Calculate the force constant matrix
!

         call getvib()


      elseif (mvflag == 10) then

!
!     Calculate the gamma surface (Alex Bratkovski's code).
!

         call gammas()
         call dump()
         call outfil(1)
         call outfil(3)


      elseif (mvflag == 11) then

!
!     Calculate the unit cell dimensions A and C (Alex Girshick's code).
!
         if      ( strtyp  ==  'hcp' )  then

             call findminhcp()

         else if ( strtyp  ==  'fcc' )  then

             call findminfcc()

         else if ( strtyp  ==  'bcc' )  then

             call findminfcc()

         else if ( strtyp  ==  'l10' )  then

             call findminhcp()

         endif

         call outfil(1)
         call outfil(3)


      elseif (mvflag == 12) then



!          do i = 0, 2000
! !            rtconf%ham_conf%b(0)%nu = 0.0 + i*0.01
!           rtconf%ham_conf%b(2)%c  = 1.0 + i*0.1

          if      ( trim(strtyp)  ==  'hcp' ) then
              call elcon_hcp()
          else if ( trim(strtyp)  ==  'hsh' ) then
              call elcon_hcp_short()
          else if ( trim(strtyp)  ==  'fcc' ) then
              call elcon_fcc()
          else if ( trim(strtyp)  ==  'bcc' ) then
              call elcon_bcc()
          else if ( trim(strtyp)  ==  'l10' ) then
              call elcon_l10()
          else if ( trim(strtyp)  ==  'lsh' ) then
              call elcon_hcp_short()
          endif

!          print *, 'elst',0, rtconf%ham_conf%b(0)%nu
! !           print *, 'elst',2, rtconf%ham_conf%b(2)%c
!           call print_elastic_consts(rtc%elastic_const_conf, " ", 6)
!          end do


      elseif (mvflag == 36) then

          if   (      trim(strtyp) == 'd19' &
               & .or. trim(strtyp) == 'l10' &
               & .or. trim(strtyp) == 'd22' &
               & .or. trim(strtyp) == 'hcp' &
               & .or. trim(strtyp) == 'b19' &
               & .or. trim(strtyp) == 'd24' ) then
              call strs_tetr_short(rtconf)
          else if (   trim(strtyp) == 'l12' &
               & .or. trim(strtyp) == 'fcc' ) then
              call strs_cube_short(rtconf)
          end if


      elseif (mvflag == 13) then

!
!     Calculate the energy surface vs. unit cell dimensions A and C.
!                                    (Alex Girshick's code)

         call getensurf()
         call outfil(1)
         call outfil(3)


      elseif (mvflag == 27) then

!   energy volume curves


         call getenvol()
         call outfil(1)
         call outfil(3)



      elseif (mvflag == 14) then

!
!     Calculate the gamma surface (Alex Girshick's code).
!

         call gammasurf()
         call outfil(1)
         call outfil(3)


      elseif (mvflag == 15) then

!
!     Calculate the grain boundary structure (Alex Girshick's code).
!

         call grainbound()
         call outfil(1)
         call outfil(3)


      elseif (mvflag == 16) then

!
!     Calculate the dislocation structure (Alex Girshick's code).
!

         call strain()
         call relax_ds()

         if (iproc == master) then
            call writexyz('xyz/last.xyz','last')
            call writexbs('bs/last.bs')
            call writecfg('cfg/last.cfg')
            call writecell('strain-cell.out')
            call writeddp('ddplot/last.dat')

            call outfil(1)
            call outfil(3)
         end if

! Relaxation like with "MVFLAG=1"
!         CALL OUTFIL(1)
!         CALL RELAX()
!         ITER = ITER - 1
!         CALL OUTFIL(3)
!         CALL DUMP()


      elseif (mvflag == 17) then

!
!     Calculate the dislocation structure under applied stress
!            (Alex Girshick's code).
!

         call stress()
         call writecell('stress-cell.out')
!         CALL PLOT(1)
!         CALL RELAX_DS()
!         CALL OUTFIL(1)
!         CALL OUTFIL(3)


      elseif (mvflag == 18) then

!
!     Build the round blocks needed to calculate the dislocation
!     self-energy (Alex Girshick's code).
!

         call makeblocks()


      elseif (mvflag == 19) then

!
!     Calculate the logariphmic plot for the dislocation self-energy
!     (Alex Girshick's code).
!

         call logplot()


      elseif (mvflag == 20) then

!
!     Calculate the Rose universal equation of state curve
!     (Alex Girshick's code).
!

         call rose_curve()


      elseif (mvflag == 21) then

!
!     Calculate the shortest distance between the atoms in the relaxed
!     configuration (Alex Girshick's code).
!

!          call distance()


      elseif (mvflag == 22) then

!
!     Screeened bond integral scaling.
!

!          call plotscr()
      continue


      elseif (mvflag == 23) then

!
!     Screeened bond integral scaling.
!

         call volrelax()

      elseif (mvflag == 24) then
!
!     Determining lattice Greens function boundary conditions (M.J.C.)
!
         call outfil(1)
         call latgfbc()
         iter = iter - 1
         call outfil(3)
         call dump()
!
      elseif (mvflag == 25) then
!
!     Performing Monte Carlo structure optimisation (M.J.C.)
!
         call mpython()
         call outfil(1)
         iter = iter - 1
         call outfil(3)
         call dump()
!
      elseif (mvflag == 26) then
!
        stop 'rfit removed!!!'
!         call env_fit()
!          call ssiman(rtconf)


      else if (mvflag == 28) then

        if (rtconf%pot_flg == 4) then
            stop 'not yet'
!             call bandfit() ! not yet
            continue
        else
            print *, 'thou shalt obey mvflag == 28 and pot_flg == 4 to use bandfit !!!'
            call panic()
        end if

      else if (mvflag == 29) then

         call print_bond_scalings()

      else if (mvflag == 30) then
         call nebf()

      endif

!
!    Close the output files.
!
      if ((mvflag > 1).and.(writeper > 0).and.(nequil < mxiter)) close(2)
      if ((mvflag > 1).and.(trace > 0)) close(3)
      if (mda_op == 1) close(4)
      if (mda_msd > 0) close(7)


      if (rtconf % fs_only == 0) then
        select case (rtconf % pot_flg)
            case (2)
                call free_ham()
                call free_chi()
                call erasab()
            case(4)
                call kspace_destroy()
        end select

        if (allocated(dq     )) deallocate(dq     )
        if (allocated(mass   )) deallocate(mass   )
        if (mag) then
         if (allocated(mg     )) deallocate(mg     )
!          if (allocated(mginert)) deallocate(mginert)
        end if
      end if


      end


