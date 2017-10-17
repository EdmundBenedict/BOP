
   function kfndef(locc,ef0) result(ef)

   !
   !    This is a routine to find the band occupancies.
   !
      use mod_precision
      use mod_const
      use mod_all_scalar, only : totnia, quiet

      implicit none

      real(dp), intent(in) :: locc, ef0
      real(dp) :: ef
      real(dp) :: eflo,efhi,nello,nelhi
      real(dp) :: nel, knumel


      real(dp), parameter :: def = 1.0_dp, maxdn = 1.0e-12_dp



!       do it = 1, 100
!           if (it > 1) nel = knumel(ef)
!
!           if (.not. quiet) then
!             write(6,'(/,"it:",i3)') it
!             write(6,'("ef,nf:",4(x,f22.14))') ef_prev, ef, dnf_prev, dnf
!           end if
!
!           if (abs(dnf) <= 1.0e-10_dp) exit
!           def = dnf * (ef - ef_prev) / ( dnf - dnf_prev )
!
!           ef_prev = ef
!           dnf_prev = dnf
!           ef = ef - def
!
!       end do

!
!    Evaluate the occupancy of the bands as a function of k point.
!
!
      ef = ef0
      nel = knumel(ef)
      if (.not. quiet) print *, 'ef,nel:',ef,nel
      if (nel < locc) then
         eflo = ef
         nello = nel
         efhi = eflo + def
   1       nelhi = knumel(efhi)
         if (.not. quiet) print *, 'efhi,nelhi:', efhi,nelhi
         if (nelhi < locc) then
            eflo = efhi
            nello = nelhi
            efhi = efhi+def
            goto 1
         endif
      else
         efhi = ef
         nelhi = nel
         eflo = efhi - def
   4       nello = knumel(eflo)
         if (.not. quiet) print *, 'eflo,nello:',eflo,nello
         if (nello > locc) then
            efhi = eflo
            nelhi = nello
            eflo = eflo-def
            goto 4
         endif
      endif



!
!    Use  binary sections to find the Fermi energy.
!

   2    ef = (eflo + efhi)*0.5_dp
      nel = knumel(ef)
      if (.not. quiet) print *, 'ef,nel:',ef,nel
      if (nel > locc) then
         nelhi = nel
         efhi = ef
      else
         nello = nel
         eflo = ef
      endif
      if (abs(nel-locc) > maxdn) goto 2

      totnia = nel

   end

!--------------------------------------------------------------------------------

   function knumel(ef) result (nel)

      use mod_precision
      use mod_all_scalar, only : kt, nsp, mag
      use mod_const
      use mod_kspace
      use topologia, only : iproc, mpmap, nproc
      use mpi, only : mpi_in_place, mpi_real8, mpi_sum, mpi_comm_world

      implicit none


      real(dp), intent(in) :: ef

      real(dp) :: nel,knel
      real(dp) :: theta


      integer :: ik,in,ierr, isp

      nel = 0.0_dp

      do isp = 1, nsp
         do ik = mpmap(iproc)+1, mpmap(iproc+1)
            knel = 0.0_dp
            do in=1,knh
               knel = knel + theta(enk(in,ik,isp)-ef,kt)
            end do
            nel = nel + wtk(ik)*knel
         end do
      end do

      if (.not. mag) nel = 2*nel


      if (nproc > 1) call mpi_allreduce(mpi_in_place, nel, 1, mpi_real8, mpi_sum, mpi_comm_world, ierr)

!     print *, nel

   end function knumel


   subroutine occupy(ef, knh, nk, enk, wtk,  occ)

      use mod_precision
      use mod_all_scalar, only : kt, quiet
      use mod_const
!       use topologia, only : iproc, mpmap

      implicit none

      integer, intent(in) :: knh, nk
      real(dp), intent(in) :: ef, enk(knh, nk), wtk(nk)
      real(dp), intent(out) :: occ(knh,nk)

      integer :: ik, in, isp

      real(dp) :: theta

      do ik = 1, nk
         do in=1, knh
            occ(in,ik) = wtk(ik)*theta(enk(in,ik)-ef,kt)
         end do
      end do



   end subroutine occupy
