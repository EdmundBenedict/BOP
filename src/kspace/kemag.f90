

   function kemag()
      use mod_precision
      use mod_const
      use mod_all_scalar, only : nd
      use mod_atom_ar, only : mg
      use topologia, only : nproc
      use mpi

      implicit none

      real(dp) :: kemag

      include '../Include/Atom.array'

      integer :: ia, ierr

      if (nproc > 1) call mpi_allreduce(mpi_in_place, mg, nd, mpi_real8, mpi_sum, mpi_comm_world, ierr)

      kemag = 0.0_dp

      do ia = 1, nd
         kemag = kemag + istn(z(ia))*(mg(ia)*mg(ia))
      end do

      kemag = 0.25_dp*kemag

   end function kemag
