
   module mod_par
      use mod_precision
      use mpi
         
      implicit none
      
    
      type trans_pack_t
         integer :: nloca, req, stat(mpi_status_size)
         real(dp), allocatable :: vals(:)
      end type
   
   contains
   
      subroutine prep_vlvec_ir_allgather(tpack,iproc,aproc,root,mpmap,ierr)

         type(trans_pack_t), intent(inout), allocatable :: tpack(:)
         integer, intent(in) :: iproc, aproc, root, mpmap(0:) 
         integer, intent(out) :: ierr

         integer :: i 

         if (iproc /= root) then
            if (mpmap(iproc+1) /= mpmap(iproc)) then
               allocate(tpack(1))
               tpack(1) % nloca = mpmap(iproc+1) - mpmap(iproc)
               allocate(tpack(1) % vals(tpack(1) % nloca ) )
            end if
         else if (aproc>1) then
      !          tm1 = mpi_wtime()
            
            allocate(tpack(aproc-1))
            do i=1,aproc-1
               tpack(i) % nloca = mpmap(i+1) - mpmap(i)
               allocate(tpack(i) % vals(tpack(i) % nloca) )
               call mpi_irecv(tpack(i) % vals, tpack(i) % nloca, mpi_real8, i, i, mpi_comm_world, tpack(i) % req, ierr )
            end do         
         
      !          tm2 = mpi_wtime()
      !          print '("m alloc",f22.14)',tm2 - tm1
         end if



      end subroutine prep_vlvec_ir_allgather
      
      
      subroutine exec_vlvec_ir_allgather(tpack,iproc,aproc,root,mpmap,n,vlvec,ierr)

         type(trans_pack_t), intent(inout), allocatable :: tpack(:)
         integer, intent(in) :: iproc, aproc, root, mpmap(0:), n
         real(dp), intent(inout) :: vlvec(n)
         integer, intent(out) :: ierr

         integer :: i 
         
         
         if (iproc /= root) then
            if (mpmap(iproc) /= mpmap(iproc+1)) then
               tpack(1)%vals =  vlvec( mpmap(iproc) + 1 : mpmap(iproc + 1)) 
               call mpi_irsend(tpack(1)%vals,tpack(1)%nloca, mpi_real8, root, iproc, mpi_comm_world, tpack(1)%req, ierr ) 
               call mpi_wait(tpack(1)%req, tpack(1)%stat, ierr)
            end if
            call mpi_bcast(vlvec, n, mpi_real8, root, mpi_comm_world, ierr)
            if (mpmap(iproc) /= mpmap(iproc+1)) then
               deallocate(tpack(1) % vals)         
               deallocate(tpack)
            end if

         else if (aproc>1) then
            do i=1,aproc-1
               call mpi_wait(tpack(i)%req, tpack(i)%stat, ierr)
               vlvec( mpmap(i) + 1 : mpmap(i + 1)) = tpack(i)%vals
            end do
            call mpi_bcast(vlvec, n, mpi_real8, root, mpi_comm_world, ierr)

            do i=1,aproc-1            
               deallocate(tpack(i) % vals)
            end do
            deallocate(tpack)
         end if

         end subroutine exec_vlvec_ir_allgather



   end module mod_par
