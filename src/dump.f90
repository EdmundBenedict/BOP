 
      subroutine dump()
          use mod_precision


          use mod_all_scalar

          use mod_const
          

!
!    This is a subroutine to write essential information out to a backup file.
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
!      include "Include/Force.array"
!      include "Include/Hamilt.array"
      include "Include/Misc.array"
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

      integer i,loop,status



     call dump_current()



!
!    Open the backup file.
!

      write(6,'(/''Dumping essential data to BOP.dump'')')

      status = 1
      loop = 1

      do while ((status /= 0).and.(loop <= 100))
         open(unit=99,file='BOP.dump',status='UNKNOWN',err=1)
         status = 0
         goto 2
 1       write(6,'(''Problem opening dump file. Trying again ...'')')
 2       loop = loop + 1
      enddo

      if (status /= 0) then
         write(6,'(''Unable to dump data. Continuing program.'')')
         return
      endif

!
!    Write out the primitive translation vectors.
!

      write(99,'(/''A'')')
      write(99,'(3G24.15)') a(1,1),a(1,2),a(1,3), & 
     &                      a(2,1),a(2,2),a(2,3), & 
     &                      a(3,1),a(3,2),a(3,3)

!
!    Write out the final atomic positions and velocities.
!

      write(99,'(/''RV'')')
      write(99,'(6G24.15)') (ad(1,i),ad(2,i),ad(3,i), & 
     &                       vel(1,i),vel(2,i),vel(3,i),i=1,nd,1)

!
!    Write out onsite energy shifts.
!

      write(99,'(''DE'')')
      write(99,'(6G14.5)') (de(i),i=1,nd,1)

!
!    Write out the atomic numbers.
!

      write(99,'(''Z'')')
      write(99,'(16I5)') (z(i),i=1,nd,1)

!
!    Write out constraint information.
!

      write(99,'(''CNST_A'')')
      write(99,'(10I5)') (cnst_a(i),i=1,cnst_n,1)

!
!    Write out running totals.
!

      write(99,'(/''ITER'')') 
      write(99,*) iter
      write(99,'(''SUME'')')
      write(99,*) sume
      write(99,'(''SUMEE'')')
      write(99,*) sumee
      write(99,'(''SUMT'')')
      write(99,*) sumt
      write(99,'(''SUMTT'')')
      write(99,*) sumtt
      write(99,'(''SUMP'')')
      write(99,*) sump
      write(99,'(''SUMPP'')')
      write(99,*) sumpp
      write(99,'(''SUMV'')')
      write(99,*) sumv
      write(99,'(''SUMVV'')')
      write(99,*) sumvv

!
!    Write out histograms.
!

      if (mda_gr > 0) then
         write(99,'(/''GR'')')
         do i = 1,mda_gr_n,1
            write(99,*) gr_bins(i)
         enddo
      endif

      if (mda_g3 > 0) then
         write(99,'(/''G3'')')
         do i = 1,mda_g3_n,1
            write(99,*) g3_bins(i)
         enddo
      endif

      if (mda_sw > 0) then
         write(99,'(/''SW'')')
         do i = 1,mda_sw_n,1
            write(99,*) sw_bins(i)
         enddo
      endif

      if ((mda_adav > 0).or.(mda_msd > 0)) then
         write(99,'(/''NOPBSUM'')')
         write(99,'(6G24.15)') (adnopb(1,i),adnopb(2,i),adnopb(3,i), & 
     &        adsum(1,i),adsum(2,i),adsum(3,i),i=1,nd,1)
      endif

      if (mda_msd == 1) then
         write(99,'(/''ADMSDA'')')
         write(99,'(3G24.15)') (admsda(i),i=1,3,1)
      endif

      if (mda_msd == 2) then
         write(99,'(/''DISPO'')')
         write(99,'(3G24.15)') (dispo(i),i=1,nd,1)
      endif

!
!    Close the output file.
!

      status = 1
      loop = 1

      do while ((status /= 0).and.(loop <= 100))
         close(99,iostat=status)
         if (status /= 0) write(6,'(''Problem closing dump file.'', & 
     &                              '' Trying again ...'')')
         loop = loop + 1
      enddo

      if (status /= 0) then
         write(6,'(''Unable to close dump file. Continuing program.'')')
         return
      endif

      end



        subroutine dump_current()
            use mod_conf
            use mod_io, only : get_new_unit, close_unit
            integer :: conf_io, cell_io
        
            conf_io = get_new_unit()
            open(unit=conf_io, file='dumpconf.dump', action='write')
            call print_run_conf(rtc, '', conf_io)
            call close_unit(conf_io)
            
            
            cell_io = get_new_unit()
            open(unit=cell_io, file='dumpcell.dump', action='write')
            call print_cell(rtc % cell, '', cell_io)
            call close_unit(cell_io)
        
            
        
        end subroutine dump_current
        



