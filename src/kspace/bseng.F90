 
    function bseng()
    use mod_precision
    use mod_kspace
    use topologia, only : iproc, mpmap
#ifdef MPI         
    use mpi, only : mpi_in_place, mpi_real8, mpi_sum, mpi_comm_world
#endif
    implicit none
    real(dp) :: bseng, ddot
    integer :: ierr

    bseng = ddot((mpmap(iproc+1)-mpmap(iproc))*knh, enk, 1, occ, 1)

#ifdef MPI    
    call mpi_allreduce(mpi_in_place, bseng, 1, mpi_real8, mpi_sum, mpi_comm_world, ierr)
#endif      
    end

