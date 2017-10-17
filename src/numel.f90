
   function numel(efermi)
      use mod_precision
      use mod_all_scalar
      use mod_const
      use topologia, only : iproc, mpmap, aproc, master, nproc
      use ab_io
      use mod_funptr, only : numelsrt
      use mod_atom_ar, only : dq, mg
!       use mod_clock

      use mpi


      implicit none

      real(dp), intent(in) :: efermi
      real(dp) :: numel

      integer :: ia,la,isp
      real(dp) :: nia
      integer :: ierr

      procedure(real(dp)) :: numelnull



      include "Include/Atom.array"




! ! ! #ifdef MPI
! ! !       integer :: i,ierr
! ! !
! ! !
! ! !       real(dp), allocatable :: buf(:,:)
! ! !       integer, allocatable :: req(:)
! ! !       integer, allocatable :: stat(:,:)
! ! ! !       integer, save :: cnt = 0
! ! !
! ! ! !       cnt = cnt +1
! ! !       integer :: n
! ! !       n = 2
! ! !
! ! !       if (iproc /= master) then
! ! !          allocate(buf(n,1),req(1),stat(mpi_status_size,1))
! ! !       else
! ! !         if (aproc>1) then
! ! !
! ! !          allocate(buf(n,(aproc-1)),req(aproc-1),stat(mpi_status_size,aproc-1))
! ! !          do i=1,aproc-1
! ! !             call mpi_irecv(buf(1,i), n, mpi_real8, i, i, mpi_comm_world, req(i), ierr )
! ! !          end do
! ! !         end if
! ! !       end if
! ! ! #endif

      numel = 0.0_dp
      dqmax = 0.0_dp

      do isp = 1, nsp
         do ia = mpmap(iproc)+1, mpmap(iproc+1)

            call assoc_ab(ia,isp)

      !       Evaluate the number of electrons.
            nia = 0.0_dp
            if (term == 1) then
                  do la = 1,nchain
!                      print *, "la,wt(la), numelsrt(efermi,la)",la,wt(la), numelsrt(efermi,la)
                     nia = nia + wt(la) * numelsrt(efermi,la)
                  enddo
            elseif ((term == 2).or.(term == 3)) then
                  do la = 1,nchain
                     nia = nia + wt(la) * numelnull(efermi,lchain(la), mrec,diag(1,la), eigvec(1,1,la),kt)
                  enddo
            endif


            if (isp == 1) then
!                   print *, "atom", ia,dq(ia)
                  dq(ia) = nia
            else ! means isp == nsp and mag == true
                  dq(ia) = dq(ia) + nia
                  dq(ia) = 0.5_dp*dq(ia)
!                   print *, "atom", ia,dq(ia)
                  mg(ia) = dq(ia) - nia  ! equiv to "up+down - 2*down = up-down" bar some numerical err on all these +-
            endif
!
!             if (isp==nsp) then ! this has been moved to spnbop just after the ef is found
!                   dq(ia) = dq(ia) - zc(z(ia))
!                   if (abs(dq(ia)) > dqmax) dqmax = abs(dq(ia))
!             end if

            numel = numel + nia
         enddo
      end do

      if (mag) numel = 0.5_dp*numel

      if (nproc > 1) call mpi_allreduce(mpi_in_place, numel, 1, mpi_real8, mpi_sum, mpi_comm_world, ierr)

! ! ! #ifdef MPI
! ! !
! ! ! ! write(6,*) 'efermi', cnt, efermi
! ! !
! ! !       if (iproc /= master) then
! ! !          if (mpmap(iproc) /= mpmap(iproc+1)) then
! ! !             buf(:,1) = (/numel, dqmax/)
! ! !             call mpi_irsend(buf,n, mpi_real8, master, iproc, mpi_comm_world, req, ierr )
! ! !             call mpi_wait(req(1),stat(:,1),ierr)
! ! !
! ! !          end if
! ! !          call mpi_bcast(buf,n, mpi_real8, master, mpi_comm_world, ierr )
! ! !
! ! !          numel = buf(1,1)
! ! !          dqmax = buf(2,1)
! ! !
! ! !          deallocate(buf,req,stat)
! ! !
! ! !       else
! ! !         if (aproc>1) then
! ! !
! ! ! !          call mpi_waitall(aproc-1,req,stat,ierr)
! ! !          do i=1,aproc-1
! ! !             call mpi_wait(req(i),stat(:,i),ierr)
! ! !
! ! !             numel = numel + buf(1,i)
! ! !             if (buf(2,i)>dqmax) dqmax = buf(2,i)
! ! !
! ! !          end do
! ! !
! ! !          buf(:,1) = (/numel, dqmax /)
! ! !
! ! !          call mpi_bcast(buf(1,1),n, mpi_real8, master, mpi_comm_world, ierr )
! ! !
! ! !
! ! !          deallocate(buf,req,stat)
! ! !         end if
! ! !       end if
! ! ! !    write(6,*)'numel',cnt,numel
! ! !
! ! ! #endif

   end function numel

