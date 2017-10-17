!reset ot mene drugoto ot Alin
module distrib
  use mpi
  implicit none
   
  private
  public :: CreateCounter, ResetCounter, DestroyCounter, GetCounter

 contains
  

  subroutine ResetCounter(winRank,init,win)
    integer, intent(in) :: winRank,win
    integer :: init, ierror
!     print *, 'predi reseta ' 
    call MPI_Win_lock(MPI_LOCK_EXCLUSIVE,winRank,0,win,ierror)
    call MPI_Put(init, 1, MPI_INTEGER, winRank, 0, 1, MPI_INTEGER, win,ierror)
    call MPI_Win_unlock(winRank,win,ierror);
!     print *, 'sled reseta '
  end subroutine ResetCounter
  
  subroutine CreateCounter(rank, winRank, counter, init, win, universe)
    integer, intent(in) :: rank, winRank, init
    integer, intent(inout),allocatable:: counter(:)
    integer, intent(inout):: win
    integer, intent(in) :: universe

    integer :: winSize,ierror


    winSize=0
    if (rank==winRank) then
      winSize=1
      allocate(counter(winSize))
      counter=init
!       write(*,*)"Process",winRank," has the following value for counter:",counter
    endif
    call MPI_Win_create(counter, winSize, 1, MPI_INFO_NULL,universe, win,ierror);

  end subroutine CreateCounter

  subroutine DestroyCounter(myrank,winRank,counter,win)
    integer, intent(in) :: myrank,winRank
    integer, intent(inout),allocatable :: counter(:)
    integer, intent(inout) :: win

    integer :: ierror
    call MPI_Win_free(win,ierror)
    if (winRank==myrank) then
      deallocate(counter)
    endif
  end subroutine DestroyCounter
  
  integer function GetCounter(winRank,increment,win)
    integer, intent(in) :: winRank,increment,win
    integer :: myCounter,ierror

    call MPI_Win_lock(MPI_LOCK_EXCLUSIVE,winRank,0,win,ierror)
    call MPI_Get(myCounter, 1, MPI_INTEGER, winRank, 0, 1, MPI_INTEGER, win,ierror)
    call MPI_Accumulate(increment, 1, MPI_INTEGER, winRank, 0, 1, MPI_INTEGER, MPI_SUM, win,ierror)
    call MPI_Win_unlock(winRank,win,ierror);
    GetCounter=myCounter;
  end function GetCounter
end module distrib


