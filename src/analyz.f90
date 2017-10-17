 
      subroutine analyz()
          use mod_precision


          use mod_all_scalar

          use mod_const

          use ab_io
          use mod_ham, only : assoc_ham
          use mod_kspace

!
!    This is a routine for analyzing a structure.
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
!       include "Include/KHamilt.array"

!
!    Declare the simple variables.
!

      real(dp) :: elo,ehi

      integer la,ne,ldos,flag,rptflg
      integer ia,choice
      integer i,j,k

!
!    Evaluate the energy and forces of a configuration.
!

      flag = 1
      call getetot(flag)




!      WRITE(6,'(/''Extra analysis option is turned off temporarily'')')
!
!      RETURN





!
!    Perform any final extra analysis called for.
!

 1    write(6,'(/''You have the following analysis options:'')')
      write(6,'('' 0) Exit'')')
      write(6,'('' 1) Calculate the total density of states'')')
      write(6,'('' 2) Calculate the partial density of states'')')
      write(6,'('' 3) Calculate the bond orders about one site'')')
      write(6,'('' 4) Calculate the weight associated with a'', & 
     &          '' range of states'')')
      write(6,'('' 5) Calculate the energy versus occupancy curves'')')
      write(6,'('' 6) Dump values'')')
      write(6,'('' 7) Write out bandstructure data'')')
      write(6,'(/''Enter your choice > '',$)')
      read(5,*) choice
!      CHOICE = 3


      if (choice == 0) then

         return

      elseif (choice == 1) then

         write(6,'(/''Enter number of intervals for the density'', & 
     &              '' of states > '',$)')
         read(5,*) ne

!Dimitar
!          if (potflg == 2) then
!             call getdos(ne) ! see dosplot
!          elseif (potflg == 4) then
!             call kdos(kpsi,ne) !see kdosplot
!          endif

      elseif (choice == 2) then

         write(6,'(/''You have the following choice of LDOS'', & 
     &              '' calculations:'')')
         write(6,'('' 1) Calculate the LDOS for one orbital.'')')
         write(6,'('' 2) Calculate the LDOS for one site.'')')
         write(6,'(/''Enter your choice > '',$)')
         read(5,*) ldos
         if (ldos == 1) then
            write(6,'(/''Enter atom number > '',$)')
            read(5,*) ia
            if ((ia <= 0).or.(ia > nd)) then
               write(6,'(/''Atom number out of range.'')')
            else
               call assoc_ab(ia,1)
               call assoc_ham(ia)
               write(6,'(/''Enter orbital number > '',$)')
               read(5,*) la
               call getnch(ia)
               if ((la <= 0).or.(la > nchain)) then
                  write(6,'(/''Orbital number out of range.'')')
               else
                  write(6,'(/''Enter number of intervals for'', & 
     &                       '' the density of states > '',$)')
                  read(5,*) ne
                  call getldos(ia,la,ne,ldos)
               endif
            endif
         elseif (ldos == 2) then
            write(6,'(/''Enter atom number > '',$)')
            read(5,*) ia
            if ((ia <= 0).or.(ia > nd)) then
               write(6,'(/''Atom number out of range.'')')
            else
               call assoc_ab(ia,1)
               call assoc_ham(ia)
               write(6,'(/''Enter number of intervals for'', & 
     &                    '' the density of states > '',$)')
               read(5,*) ne
               call getldos(ia,la,ne,ldos)
            endif
         else
            write(6,'(/''Unknown choice. Please try again.'')')
         endif

      elseif (choice == 3) then

         write(6,'(/''Enter atom number > '',$)')
         read(5,*) ia
!         IA = 1
         if (potflg  ==  2) then
            if ((ia <= 0).or.(ia > nd)) then
               write(6,'(/''Atom number out of range.'')')
            else
               call assoc_ab(ia,1)
               call assoc_ham(ia)
               if (momflg == 1) then
                  call boavg(ia)
!                  CALL BOAVG_MJC(IA)
               elseif ((momflg == 2).or.(momflg == 3)) then
                  call bonoav(ia)
               endif
            endif
         elseif (potflg  ==  4) then
            if ((ia <= 0).or.(ia > nd)) then
               write(6,'(/''Atom number out of range.'')')
            else
               call assoc_ab(ia,1)
               call assoc_ham(ia)
               stop 'kbo excluded until done well'
!                call kbo(ia)
            endif
         endif
!
      elseif (choice == 4) then

         write(6,'(''Enter lower energy bound > '',$)')
         read(5,*) elo
         write(6,'(''Enter upper energy bound > '',$)')
         read(5,*) ehi
         call getrho(elo,ehi)
      elseif (choice == 5) then
         write(6,'(/''Enter the number of occupancies to be'', & 
     &              '' considered > '',$)')
         read(5,*) nfermi
         if (nfermi < 2) then
            write(6,'(/''Must have at least two occupancies.'')')
         elseif (nfermi > mfermi) then
            write(6,'(/''Too many occupancies. Maximum is '', & 
     &                I4,''.'')') mfermi
         else
            stop 'getevsocc was never really well integrated. needs to be brought up to date'
!             call getevsocc()
         endif

      elseif (choice == 6) then

         write(6,'(/'' 1) - Check that the onsite and intersite'')')
         write(6,'(''      rec. coefficients obey the sum rule.'')')
         write(6,'('' 2) - Write out recursion coefficients.'')')
         write(6,'('' 3) - Write out eigenvalues and vectors '')')
         write(6,'(''      for no terminator.'')')
         write(6,'('' 4) - Write out susceptibilities.'')')
         write(6,'('' 5) - Write out some recursion coefficients.'')')
         write(6,'(/''Enter report flag value > '',$)')
         read(5,*) rptflg
         call report(rptflg)

      elseif (choice == 7) then
         
         if (potflg /= 4) then
            write(6,'(/''Bandstructure data not available for this'', & 
     &           '' potential.'')')            
         else
            open(71,file='BS.dat')            
            do j = 1, nk, 1
               do i = 1, kmxh
                  write(71,'(I6,F12.5)') j,enk(i,j)
               enddo            
            enddo
            close(71)
            write(6,'(/''Bandstructure data written to BS.dat.'')')
         endif

      else
         write(6,'(/''Unknown choice. Please try again.'')')
      endif

      goto 1

      end

