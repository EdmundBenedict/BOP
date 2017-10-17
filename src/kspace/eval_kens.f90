   subroutine eval_kens(ef)

   use mod_precision
   use mod_const
   use mod_all_scalar
   use mod_kspace
   use topologia
   use mpi
   implicit none

   real(dp), intent(in) :: ef
   real(dp) :: bseng, keprom, kentropy, kemag, ddot

   integer, parameter :: nseats = 3
   real(dp) :: train(nseats)
   integer :: ierr




   !    eband = bseng() !    Evaluate the band energy.
   eband = ddot((mpmap(iproc+1)-mpmap(iproc))*knh*nsp, enk, 1, occ, 1)

   if (.not. mag) eband = 2*eband

   eprom = keprom() !    Evaluate the onsite energy.
   uent = kentropy(ef) !    Evaluate the electron entropy
   if (mag) emag = kemag()



   if (nproc > 1) then
       train = [eband, eprom, uent] !eprom is already made for all kpts due to rho
       call mpi_allreduce(mpi_in_place, train, nseats, mpi_real8, mpi_sum, mpi_comm_world, ierr)
       eband = train(1)
       eprom = train(2)
       uent  = train(3)
   end if

   ebond = eband - eprom !    Evaluate the bond energy.


   print *, 'eband: ', eband
   print *, 'ebond: ', ebond
   print *, 'eprom: ', eprom
   print *, 'uent:  ', uent
   if (mag) print *, 'emag:  ', emag

   end subroutine eval_kens
