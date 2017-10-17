   module mpi

! MPI replacement for when getting a real MPI library is too much hassle.
! This simulates the availability of MPI environment running with 1 process.
! The mpif.h file is also available for those objecting to the use of modules.
!
! Dimitar Pashov <d.pashov@gmail.com>

      use iso_c_binding, only : c_ptr, c_size_t

      implicit none

      logical :: fmpi_on = .false., fmpi_off = .false.

      integer :: mpi_comm_world = 0, &
                 mpi_real8      = 0, &
                 mpi_complex16  = 0, &
                 mpi_real       = 0, &
                 mpi_logical    = 0, &
                 mpi_integer    = 0, &
                 mpi_2integer   = 0, &
                 mpi_character  = 0, &
                 mpi_max        = 0, &
                 mpi_maxloc     = 0, &
                 mpi_sum        = 0, &
                 mpi_datatype_null = 0, &
                 mpi_order_fortran = 0, &
                 mpi_request_null  = 0, &
                 mpi_max_processor_name = 0, &
                 mpi_distribute_cyclic  = 0

      integer, parameter :: mpi_status_size = 0

!       integer, parameter :: mpi_address_kind = c_size_t
      integer, parameter :: mpi_address_kind = 8 ! it shalt be c_size_t but the mpif.h has no sane way of getting hold of it.
                                                 ! this will now only work for 64bit machines I guess... or it may not matter
                                                 ! what it is as long as it is the same in mpif.h


! These are special values defined by the standard, only tentatively represented by 'integer' type in the real implementations.
! More appropriate representation is 'type(c_ptr), pointer'
! Update: Since this is noo longer a wrapper but rather a replacement we may as well not bother with c_ptr
!       type(c_ptr), pointer ::  mpi_bottom          => null(), &
!                            & mpi_status_ignore   => null(), &
!                            & mpi_errcodes_ignore => null(), &
!                            & mpi_in_place        => null(), &
!                            & mpi_argv_null       => null(), &
!                            & mpi_argvs_null      => null(), &
!                            & mpi_statuses_ignore => null(), &
!                            & mpi_unweighted      => null()

      integer :: mpi_bottom          = 0, &
                 mpi_status_ignore   = 0, &
                 mpi_errcodes_ignore = 0, &
                 mpi_in_place        = 0, &
                 mpi_argv_null       = 0, &
                 mpi_argvs_null      = 0, &
                 mpi_statuses_ignore = 0, &
                 mpi_unweighted      = 0


      integer(8) :: clock_init, clock_rate, clock_max
      real(8) :: inv_clock_rate

      contains

      function mpi_wtime() result(r)
!          use mpi
!          implicit none
         integer(8) :: ticks
         real(8) :: r
         call system_clock(ticks)
         r = real(ticks - clock_init, 8)*inv_clock_rate
      end function mpi_wtime

   end module mpi


   subroutine mpi_initialized(flag, ierror)
      use mpi
      implicit none
      logical :: flag
      integer :: ierror
      flag = fmpi_on
      ierror = 0
   end subroutine mpi_initialized


   subroutine mpi_init(ierror)
      use mpi
      use iso_c_binding, only : c_loc, c_f_pointer
      implicit none
      integer :: ierror
      ierror = 0
      fmpi_on = .true.

      call system_clock(count_rate = clock_rate, count_max = clock_max)
      call system_clock(clock_init)
      inv_clock_rate = 1.0_8/clock_rate
   end subroutine mpi_init


   subroutine mpi_barrier(comm, ierror)
      use mpi
      implicit none
      integer :: comm, ierror
      ierror = 0
   end subroutine mpi_barrier


   subroutine mpi_finalized(flag, ierror)
      use mpi
      implicit none
      logical :: flag
      integer :: ierror
      flag = fmpi_off
      ierror = 0
   end subroutine mpi_finalized


   subroutine mpi_finalize(ierror)
      use mpi
      implicit none
      integer :: ierror
      ierror = 0
      fmpi_off = .true.
   end subroutine mpi_finalize


   subroutine mpi_comm_rank(comm, rank, ierror)
      use mpi
      implicit none
      integer :: comm, rank, ierror
      rank = 0
      ierror = 0
   end subroutine mpi_comm_rank


   subroutine mpi_comm_size(comm, size, ierror)
      use mpi
      implicit none
      integer :: comm, size, ierror
      size = 1
      ierror = 0
   end subroutine mpi_comm_size


   subroutine mpi_bcast(buffer, count, datatype, root, comm, ierror)
      use mpi
      implicit none
      type(c_ptr) :: buffer(*)
      integer ::  count, datatype, root, comm, ierror
      ierror = 0
   end subroutine mpi_bcast


   subroutine  mpi_comm_split(comm, color, key, newcomm, ierror)
      use mpi
      implicit none

      integer :: comm, color, key, newcomm, ierror
      ierror = 0
      newcomm = comm
   end subroutine mpi_comm_split


   subroutine mpi_get_processor_name(name, resultlen, ierror)
      use mpi
      implicit none
      character*(*) ::  name
      integer :: resultlen, ierror
      ierror = 0
      resultlen = 0
      name = ''
   end subroutine mpi_get_processor_name

   subroutine mpi_allreduce(sendbuf, recvbuf, count, datatype, op, comm, ierror)
      use mpi
      implicit none
      integer ::  count, datatype, op, comm, ierror
      type(c_ptr) :: sendbuf(*), recvbuf(*)
      ierror = 0
   end subroutine mpi_allreduce

   subroutine mpi_allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, ierror)
      use mpi
      implicit none
      type(c_ptr) ::  sendbuf (*), recvbuf (*)
      integer ::  sendcount, sendtype, recvcount, recvtype, comm, ierror
      ierror = 0
   end subroutine mpi_allgather

   subroutine mpi_allgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcount, displs, recvtype, comm, ierror)
      use mpi
      implicit none
      type(c_ptr) ::  sendbuf (*), recvbuf (*)
      integer ::  sendcount, sendtype, recvcount(*), displs(*), recvtype, comm, ierror
      ierror = 0
   end subroutine mpi_allgatherv

   subroutine mpi_comm_free(comm, ierror)
      use mpi
      implicit none
      integer ::  comm,  ierror
      ierror = 0
   end subroutine mpi_comm_free

   subroutine mpi_type_contiguous(count, oldtype, newtype, ierror)
      use mpi
      implicit none
      integer :: count, oldtype, newtype, ierror
      ierror = 0
   end subroutine mpi_type_contiguous

   subroutine mpi_type_commit(datatype, ierror)
      use mpi
      implicit none
      integer :: datatype, ierror
      ierror = 0
   end subroutine mpi_type_commit


   subroutine mpi_type_free(datatype, ierror)
      use mpi
      implicit none
      integer :: datatype, ierror
      ierror = 0
   end subroutine mpi_type_free



   subroutine mpi_cartdim_get(comm, ndims, ierror)
!    DO NOT set ndims to anything since it could be anything positive even with 1 process.
      use mpi
      implicit none
      integer ::  comm, ndims, ierror
      ierror = 0
   end subroutine mpi_cartdim_get


   subroutine mpi_cart_get(comm, maxdims, dims, periods, coords, ierror)
      use mpi
      implicit none
      integer :: comm, maxdims, dims(*), coords(*), ierror
      logical :: periods(*)
      ierror = 0
      dims(1:maxdims) = 1
      periods(1:maxdims) = .false.
      coords(1:maxdims) = 0
   end subroutine mpi_cart_get



   subroutine mpi_dims_create(nnodes, ndims, dims, ierror)
      use mpi
      implicit none
      integer :: nnodes, ndims, dims(*), ierror
      ierror = 0
   end subroutine mpi_dims_create

   subroutine mpi_cart_create(comm_old, ndims, dims, periods, reorder, comm_cart, ierror)
      use mpi
      implicit none
      integer :: comm_old, ndims, dims(*), comm_cart, ierror
      logical :: periods(*), reorder
      ierror = 0
      comm_cart = comm_old
   end subroutine mpi_cart_create


   subroutine mpi_cart_sub(comm, remain_dims, comm_new, ierror)
      use mpi
      implicit none
      integer :: comm, comm_new, ierror
      logical :: remain_dims(*)
      ierror = 0
      comm_new = comm
   end subroutine mpi_cart_sub


   subroutine mpi_cart_rank(comm, coords, rank, ierror)
      use mpi
      implicit none
      integer :: comm, coords(*), rank, ierror
      ierror = 0
      rank = 0
   end subroutine mpi_cart_rank


   subroutine mpi_cart_coords(comm, rank, maxdims, coords, ierror)
      use mpi
      implicit none
      integer :: comm, rank, maxdims, coords(*), ierror
      ierror = 0
      coords(1:maxdims) = 0
   end subroutine mpi_cart_coords

   subroutine mpi_gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, ierror)
      use mpi
      implicit none
      type(c_ptr) ::    sendbuf(*), recvbuf(*)
      integer :: sendcount, sendtype, recvcount, recvtype, root
      integer :: comm, ierror
      ierror = 0
   end subroutine mpi_gather

   subroutine mpi_gatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, root, comm, ierror)
      use mpi
      implicit none
      type(c_ptr) :: sendbuf(*), recvbuf(*)
      integer :: sendcount, sendtype, recvcounts(*), displs(*)
      integer :: recvtype, root, comm, ierror
      ierror = 0
   end subroutine mpi_gatherv

   subroutine mpi_scatterv(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount, recvtype, root, comm, ierror)
      use mpi
      implicit none
      type(c_ptr) :: sendbuf(*), recvbuf(*)
      integer :: sendcounts(*), displs(*), sendtype
      integer :: recvcount, recvtype, root, comm, ierror
      ierror = 0
   end subroutine mpi_scatterv


   subroutine mpi_alltoallw(sendbuf, sendcounts, sdispls, sendtypes, recvbuf, recvcounts, rdispls, recvtypes, comm, ierror)
      use mpi
      implicit none
      type(c_ptr) :: sendbuf(*), recvbuf(*)
      integer :: sendcounts(*), sdispls(*), sendtypes(*)
      integer :: recvcounts(*), rdispls(*), recvtypes(*)
      integer :: comm, ierror
      ierror = 0
   end subroutine mpi_alltoallw



   subroutine mpi_type_create_darray(sz, rank, ndims, array_of_gsizes, array_of_distribs, &
                                          & array_of_dargs, array_of_psizes, order,oldtype, newtype, ierror)
      use mpi
      implicit none
      integer ::  sz, rank, ndims, array_of_gsizes(*), array_of_distribs(*), &
            & array_of_dargs(*), array_of_psizes(*), order, oldtype, newtype, ierror
      ierror = 0
   end subroutine mpi_type_create_darray


   subroutine mpi_sendrecv_replace(buf, count, datatype, dest, sendtag, source, recvtag, comm, status, ierror)
      use mpi
      implicit none
      type(c_ptr) :: buf(*)
      integer :: count, datatype, dest, sendtag
      integer :: source, recvtag, comm
      integer :: status(mpi_status_size), ierror
      ierror = 0
   end subroutine mpi_sendrecv_replace


   subroutine mpi_irecv(buf, count, datatype, source, tag, comm, request, ierror)
      use mpi
      implicit none
      type(c_ptr) ::   buf(*)
      integer :: count, datatype, source, tag, comm, request, ierror
      ierror = 0
   end subroutine


   subroutine mpi_irsend(buf, count, datatype, dest, tag, comm, request, ierror)
      use mpi
      implicit none
      type(c_ptr) ::  buf(*)
      integer :: count, datatype, dest, tag, comm, request, ierror
      ierror = 0
   end subroutine mpi_irsend

   subroutine mpi_wait(request, status, ierror)
      use mpi
      implicit none
      integer :: request, status(mpi_status_size), ierror
      ierror = 0
   end subroutine mpi_wait

   subroutine mpi_waitall(count, array_of_requests, array_of_statuses, ierror)
      use mpi
      implicit none
      integer :: count, array_of_requests(*)
      integer :: array_of_statuses(mpi_status_size,*), ierror
      ierror = 0
   end subroutine mpi_waitall


!    subroutine mpi_abort(comm, errorcode, ierror)
!       use mpi
!       implicit none
!       integer   comm, errorcode, ierror
!       ierror = 0
!    end subroutine mpi_abort


   subroutine mpi_get_version(version, subversion, ierror)
      use mpi
      implicit none
      integer :: version, subversion, ierror
      version = 2
      subversion = 2
      ierror = 0
   end subroutine mpi_get_version



