 
      subroutine wrthis()
          use mod_precision


          use mod_all_scalar

          use mod_const

!
!    This is a subroutine to write out the histograms.
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
      include "Include/Misc.array"
      include "Include/PosVel.array"

!
!    Declare the simple variables.
!

      real(dp) :: norm
      real(dp) :: niter

      integer i

      character*80 filename

!
!    Write out histograms.
!

      if (iter > nequil) then

         niter = real(iter-nequil, dp)

!
!       Write out pair correlation function.
!

         if (mda_gr > 0) then
            filename = genfile(1:lengfn)//'.gr'
            open(unit = 3,file = filename,status = 'NEW')
            do i = 1,mda_gr_n,1
               norm = 4.0_dp*pi*((real(i, dp)-0.5_dp)**2)*(mda_gr_dr**3)* & 
     &                niter*(real(nd, dp)**2)/avgv
               write(3,*) (real(i, dp)-0.5_dp)*mda_gr_dr, & 
     &                     real(gr_bins(i), dp)/norm
            enddo
            close(3)
         endif

!
!       Write out angular correlation function.
!

         if (mda_g3 > 0) then
            filename = genfile(1:lengfn)//'.g3'
            open(unit = 3,file = filename,status = 'NEW')
            do i = 1,mda_g3_n,1
               write(3,*) (real(i, dp)-0.5_dp)*180.0_dp/real(mda_g3_n, dp), & 
     &                     real(g3_bins(i), dp)
            enddo
            close(3)
         endif

!
!       Write out the distribution of Stillinger-Weber three-body
!        energies.
!

         if (mda_sw > 0) then
            filename = genfile(1:lengfn)//'.sw'
            open(unit = 3,file = filename,status = 'NEW')
            do i = 1,mda_sw_n,1
               norm = mda_sw_de*niter*real(nd, dp)
               write(3,*) (real(i, dp)-0.5_dp)*mda_sw_de, & 
     &                     real(sw_bins(i), dp)/norm
            enddo
            close(3)
         endif

      endif


!
!       Write out time averaged positons.
!
 
      if (mda_adav > 0) then
         filename = genfile(1:lengfn)//'.adav'
         open(unit = 3,file = filename,status = 'NEW')
         do i = 1, nd, 1
            write(3,'(3F12.5)') adsum(1,i)/iter,adsum(2,i)/iter, &
     &           adsum(3,i)/iter
         enddo
         close(3)
      endif

      end

