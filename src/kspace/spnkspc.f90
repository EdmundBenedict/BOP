

   subroutine spnkspc()

      use mod_precision
      use mod_const
      use mod_all_scalar, only : lef, quiet, dqmax, nd, locc, totnia, nsp, mag
      use mod_kspace
      use topologia
      use mod_io
      use mod_atom_ar, only : dq, mg
      use mpi


      implicit none

      include '../Include/Atom.array'

      character(len=25) :: ffmt
      real(dp) :: ef
      procedure(real(dp)) :: kfndef
      logical :: gpt
      integer :: isp, ik, ia, ierr, lnk, mxidx(1), ih
      real(dp) :: v1, v2

      ef = lef

      if (.not. quiet) write(ffmt,"(a,i0,a)") '(a,":",', nd, '(x,f22.14))'

      do isp = 1, nsp
         do ik = mpmap(iproc)+1, mpmap(iproc+1)
            call shiftons(isp, knh, kham(1,1,ik), kpsi(1,1,ik,isp))
            gpt = all(kk(:,ik) == 0.0_dp)
            call kdiag(gpt, knh, enk(1,ik,isp), kpsi(1,1,ik,isp),kovl(1,1,ik))
         end do
      end do

      ef = kfndef(locc,lef)
      lef = ef

      do isp = 1, nsp
         ik = mpmap(iproc)+1
         lnk = mpmap(iproc+1) - mpmap(iproc)
         call occupy(ef, knh, lnk, enk(1,ik,isp), wtk,  occ(1,ik,isp))
         call krhodiag(knh, lnk, occ(1,ik,isp), kpsi(1,1,ik,isp),  krhod(1,isp),kovl(1,1,ik))
      end do

!       print out gamma point enk, occ
      if (.not. quiet) then
         do ik = mpmap(iproc)+1, mpmap(iproc+1)
            if (all(kk(:,ik) == 0.0_dp)) then
               write(6, '("occ  en")')
               do ih = 1, knh
                  do isp = 1, nsp
                     write(6, '(2(x,es20.12))', advance='no') enk(ih,ik,isp), occ(ih,ik,isp)
                  end do
                  write(6,'()')
               end do
            end if
         end do
      end if

!       call bsinter()

      call atqsym(nd, knh, krhod(1,1),  dq)
      if (mag) then
         call atqsym(nd, knh, krhod(1,2),  mg)
         do ia = 1, nd
           
            v1 = dq(ia) + mg(ia)
            v2 = dq(ia) - mg(ia)
            dq(ia) = v1
            mg(ia) = v2
         end do
      else
         dq = 2*dq
      end if

      if (nproc > 1) call mpi_allreduce(mpi_in_place, dq, nd, mpi_real8, mpi_sum, mpi_comm_world, ierr)

      do ia = 1, nd
          print *, ia,dq(ia)
         dq(ia) = dq(ia) - zc(z(ia))
      end do

      mxidx = maxloc(abs(dq))
      dqmax = dq(mxidx(1))


      if (.not. quiet) then
         write(6, ffmt) 'dq', dq(:nd)
         write(6,'("dqmax:",x,f22.14,2x,"ef:",x,f22.14,2x,"numel:",x,f22.14,x,f22.14)') dqmax, ef, totnia, sum(dq(:nd))
      end if


      call eval_kens(ef)



   end subroutine spnkspc
