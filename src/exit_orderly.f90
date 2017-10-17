
    subroutine exit_orderly()
        use mod_io
        use ab_io  , only : free_ab
        use mod_ham, only : free_ham
        use mod_chi, only : free_chi
        use topologia
        use mpi
        implicit none

        integer :: ierror

        call free_ab()
        call free_ham()
        call free_chi()
        call close_io()
        call close_all_files()
        call mpi_finalize(ierror)
        print *, 'proc ', iproc, 'exitting...'
        stop

    end subroutine exit_orderly






