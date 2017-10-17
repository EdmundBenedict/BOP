   subroutine eval_bsens(efermi)
      use mod_precision
      use mod_const
      use mod_all_scalar
      use mod_chi
      use ab_io
      use mod_ham
      use mod_funptr
      use topologia, only : iproc, mpmap, aproc, master, nproc

      use mpi

      implicit none

      include "Include/Force.array"

      real(dp), intent(in) :: efermi

      integer :: ia, la, i, ierr, isp
      real(dp) :: bsons, error
      procedure(real(dp)) :: entropy,febond,  dndmfn, femag
!       real(dp) :: dndmfn

      real(dp), allocatable :: buf(:,:)
      integer, allocatable :: req(:)
      integer, allocatable :: stat(:,:)

      integer :: n

      if (nproc > 1) then
          n = 6

          if (iproc /= master) then
             allocate(buf(n,1),req(1),stat(mpi_status_size,1))
          else if (aproc>1) then
    !          tm1 = mpi_wtime()

             allocate(buf(n,aproc-1),req(aproc-1),stat(mpi_status_size,aproc-1))
             do i=1,aproc-1
                call mpi_irecv(buf(1,i), n, mpi_real8, i, i, mpi_comm_world, req(i), ierr )
             end do
          end if
      endif



      dndm  = 0.0_dp
      eprom = 0.0_dp
      ebond = 0.0_dp
      uent  = 0.0_dp
      bsons = 0.0_dp
      emag  = 0.0_dp

      do isp = 1, nsp

         do ia = mpmap(iproc)+1, mpmap(iproc+1)

            call assoc_ab(ia, isp)
            call assoc_chi(ia, isp)
            call assoc_ham(ia)


      ! !       Evaluate the rate of change of number of electrons with chemical potential.
      !          do la = 1,nchain
      !             dndm = dndm + dndmfn(efermi,la)
      !          enddo

            eprom = eprom + feprom(ia,efermi,isp)
            ebond = ebond + febond(ia,efermi)
            uent  = uent  + entropy(efermi)
            bsons = bsons + ebsos(ia,efermi)

         end do

      end do

      if (mag) then
         eprom = 0.5_dp*eprom
         ebond = 0.5_dp*ebond
         uent  = 0.5_dp*uent
         bsons = 0.5_dp*bsons
         emag  = femag()
      end if

!       write(6,'("predi bond, prom, onsite:", 3(x,f12.5))') ebond, eprom, bsons

      if (nproc > 1) then

          if (iproc /= master) then
             if (mpmap(iproc) /= mpmap(iproc+1)) then
                buf(:,1) = (/ dndm, eprom, ebond, uent, bsons, emag /)
                call mpi_irsend(buf, n, mpi_real8, master, iproc, mpi_comm_world, req(1), ierr )
                call mpi_wait(req(1),stat(:,1),ierr)
             end if

             call mpi_bcast(buf,n, mpi_real8, master, mpi_comm_world, ierr )

             dndm  = buf(1,1)
             eprom = buf(2,1)
             ebond = buf(3,1)
             uent  = buf(4,1)
             bsons = buf(5,1)
             emag  = buf(6,1)

             deallocate(buf,req,stat)

          else if (aproc>1) then
    !          tm1 = mpi_wtime()

             call mpi_waitall(aproc-1,req,stat,ierr)
             do i=1,aproc-1
                dndm  = dndm  + buf(1,i)
                eprom = eprom + buf(2,i)
                ebond = ebond + buf(3,i)
                uent  = uent  + buf(4,i)
                bsons = bsons + buf(5,i)
                emag  = emag  + buf(6,i)
             end do

             buf(:,1) = [dndm, eprom, ebond, uent, bsons, emag]

             call mpi_bcast(buf(1,1), n, mpi_real8, master, mpi_comm_world, ierr )


             deallocate(buf,req,stat)
          end if


      endif

! write(6,'("sled bond, prom, onsite:", 3(x,f12.5))')  ebond, eprom, bsons

!
!    Test the internal consistency.
!
!          print '("ia,bond, prom, onsite:", i0,3(x,f12.5))', ia, ebond, eprom, bsons
      eband = eprom + ebond
      error = abs((bsons-eband)/(bsons+eband))



      if ((.not. quiet) .and. (iproc == master)) then
        print *, 'bsons: ', bsons
        print *, 'eband: ', eband
        print *, 'ebond: ', ebond
        print *, 'eprom: ', eprom
        if (mag) print *, 'emag : ', emag

        print *, 'bserr:', error
!
         if (error > 10.0_dp*maxerr) then
            write(6,'("Unacceptably large uncertainty in the band structure energy.")')
            write(6,'("Uncertainty = ",f12.2,"%,"," on/inter",2(x,f20.12))') error*100.0_dp, bsons, eband
            write(6,'(''==> Terminating calculation.'')')
!             call panic()
         elseif (error > maxerr) then
            write(6,'("WARNING: Large uncertainty in the band structure energy.")')
            write(6,'("Uncertainty = ",f12.2,"%,"," on/inter",2(x,f20.12))') error*100.0_dp, bsons, eband
         endif

      end if



   end subroutine eval_bsens
