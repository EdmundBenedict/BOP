
   subroutine forcecheck(rtconf)
      use mod_precision
      use mod_all_scalar
      use mod_const
      use mod_conf
      use mod_io
      use mod_atom_ar
      use topologia
      use mpi

      implicit none

!      Prints out a table with energies and forces in forcecheck.out for each displacement of the specified atom along the specified coordinate.
!     The table may be useful for checking whether the energies and forces are consistent. There is a python script available to plot it all.


      type(run_conf_t), intent(inout) :: rtconf


      include "Include/Atom.array"
      include "Include/Force.array"
      include "Include/PosVel.array"
      include "Include/NebList.array"
      include "Include/ag.conn"



      integer :: ia, crd, np, i, u
      real(dp) :: mn, mx, dlt, shift, p0

      integer :: itrain(3), err
      real(dp) :: rtrain(2)

      if (iproc == master) then
         write(6, '("atom, coord,  min, max, npoints?:")') !, advance='no')
         read(5, *) ia, crd, mn, mx, np
      end if

      if (nproc > 1) then
          if (iproc == master) then
             itrain = [ia, crd, np]
             rtrain = [mn, mx]
          endif

          call mpi_bcast(itrain, 3, mpi_integer, master, mpi_comm_world, err)
          call mpi_bcast(rtrain, 2, mpi_real8, master, mpi_comm_world, err)

          if (iproc /= master) then
             ia  = itrain(1)
             crd = itrain(2)
             np  = itrain(3)
             mn  = rtrain(1)
             mx  = rtrain(2)
          end if
      end if

      print *, 'iproc,ia,crd,mn,mx,np:', iproc,ia,crd,mn,mx,np
      if (ia > nd .or. crd > 3 .or. crd < 1) stop 'ia > nd .or. crd > 3 .or. crd < 1'

      dlt = (mx - mn)/real(np - 1, dp)
      ad(:,:nd) = matmul(transpose(a), d(:,:nd))
      p0 = ad(crd,ia)

      if (.not. quiet) then
         u = get_new_unit()
         open(u, file='forcecheck.out', action='write')
         write(u,'(a)') 'forcecheck'
         write(u,'(3x,a,1x,i4,1x,a,1x,a,/,3x,a,8x,a,17x,a,16x,a,17x,a,16x,a)') &
            & 'atom:', ia, 'coord:', crdnam(crd), 'n', 'shift', 'ebond', 'fbond', 'epair', 'fpair'
      end if

      do i = 0, np - 1
         shift = mn + i*dlt
         ad(crd,ia) = p0 + shift
!          de = 0.0_dp
!          mg = 2.85_dp
         call getetot(0)
         if (.not. quiet) write(u, '(i4,x,es20.12,2(x,2(x,es20.12)))') &
         & i, shift, eelct, fbs(crd,ia), epair, fpp(crd,ia)
      end do

      if (.not. quiet) call close_unit(u)


   end subroutine forcecheck





