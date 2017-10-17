
! mostly taken from Matous

   module mod_neb


   use mod_precision
   use mod_const
   use mod_io
   use mod_atom_ar, only : demi, mg
   use mod_all_scalar, only : nd, ninert, ebind, mag, quiet, mxiter, rlxflg, step, ftol
   use topologia, only : iproc

!    use ieee_exceptions

   use mpi

   implicit none

   private

   procedure(real(dp)) :: ddot
   procedure(integer) :: usrexit

   integer :: nimg, climb, nx, u, ierr
   real(dp) :: spring_k

   real(dp), allocatable :: neb_ad(:,:), neb_adinert(:,:), neb_f(:,:), &
               & neb_e(:), neb_de(:,:), neb_mg(:,:), neb_demi(:,:), ks(:), alph(:), tau(:,:), neb_dist(:)

   integer :: current_img

   public :: nebf

   contains

   subroutine nebf()

      integer :: it, i, j
      real(dp) :: fmaxtot, diffe(2), emax, emin, ezero, invl

      integer :: iemax(1)

      real(dp), allocatable :: velneb(:,:)


!       these shall be put somewhere else:

      integer :: cut, cuts, minsteps, v_mix, fmax_idx
      real(dp) :: kmin, kmax, cur_a, cur_v, max_a, max_v, vf, vg_dot_vg, &
                  & fg_dot_fg, kd, ft, disp, help, dt, incfac, decfac, dispmax, &
                  & maxdr, maxdt, mix, mix_in, mixdec, fmax, loc, df

      character(len=18), parameter :: imgfmt = '(i3, 8(x, f22.16))'

      kmin = 0.0_dp
      kmax = 45.0_dp
      dispmax = 0.02_dp
      
!kmax=0.02_dp
!       FIRE init
      max_a = 0.0_dp
      max_v = 0.0_dp
      maxdr = 0.1_dp
      v_mix = 1
      mix_in = 0.1_dp
      minsteps = 5
      incfac = 1.1_dp
      decfac = 0.5_dp
      cut = 0
      cuts = 0
      mix = mix_in
      mixdec = 0.99_dp
      dt = 0.1_dp

      nx = 3*nd
      print *, 'KMAX EQUALS', kmax
      call init()

      if (mxiter < 1) return

      u = iproc + 5640

      it = 0

      if (rlxflg == 4) then
         allocate(velneb(nx, 2:nimg-1))
         velneb = 0.0_dp
      end if

      call exec_img(1, it)
      call exec_img(nimg, it)

      ezero = min(neb_e(1), neb_e(nimg))

      do
         do i = 2, nimg - 1
            call exec_img(i, it)
         end do

         iemax = maxloc(neb_e)
         emax = neb_e(iemax(1))
         emin = minval(neb_e)


         do i = 2, nimg-1
            diffe = [neb_e(i) - neb_e(i-1), neb_e(i+1) - neb_e(i)]
            if ( all(diffe > 0.0_dp) ) then
               tau(:, i) = neb_ad(:, i+1) - neb_ad(:, i)
            elseif ( all(diffe < 0.0_dp) ) then
               tau(:, i) = neb_ad(:, i) - neb_ad(:, i-1)
            else
               diffe = [minval(abs(diffe)), maxval(abs(diffe))]
               if (neb_e(i+1) > neb_e(i-1)) diffe = cshift(diffe, 1)
               tau(:, i) = diffe(1)*(neb_ad(:,i+1) - neb_ad(:, i)) + diffe(2)*(neb_ad(:, i) - neb_ad(:, i-1))
            endif
            invl = 1.0_dp/sqrt(ddot(nx, tau(1, i), 1, tau(1, i), 1))
            tau(:, i) = invl * tau(:, i)
         enddo



         do i = 2, nimg

            if (neb_e(i) > ezero) then
               ks(i) = kmax
!               ks(i) = kmin + (kmax - kmin)*(emax - neb_e(i))/(emax - ezero)
            else
               ks(i) = kmax
            endif
!             df = neb_ad(:,i) - neb_ad(:,i-1)
!             neb_dist(i) = sqrt(ddot(nx, df, 1, df, 1))
            print *, ks(i)
            loc = 0.0_dp
            do j = 1, nx
               df = neb_ad(j,i) - neb_ad(j,i-1)
               loc = loc + df*df
            end do
            neb_dist(i) = sqrt(loc)

!             scalxy = 0.0_dp
!             do k = 1, nx, 3
!                scalxy = scalxy + (neb_ad(k,i)-neb_ad(k,i-1))**2 + (neb_ad(k+1,i)-neb_ad(k+1,i-1))**2
!             enddo
!             distxy(i) = sqrt(scalxy)

!             write(6000+iproc, '(3(x,g0))') i, loc, neb_dist(i)
         enddo







         call mpi_barrier(mpi_comm_world, ierr)
!
!          call panic()

         write(u,"('it: ', i5, /, 'i', 7x, 'e', 12x, 'e - e_1', 6x, 'e - e_min', 4x, 'k', 12x, &
                                             & 'd', 12x, 'max(Fn)', 6x, 'Fs', 11x, 'max(Fsp)', 5x, 'max(F)')") it

         i = 1
         write(u, imgfmt, advance='no') i, neb_e(i), neb_e(i) - neb_e(1), &
                                 & neb_e(i) - emin, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp

         fmaxtot = -1.0_dp

         do i = 2, nimg-1
            ft = ddot(nx, neb_f(1, i), 1, tau(1, i), 1)
            kd = ks(i+1)*neb_dist(i+1) - ks(i)*neb_dist(i)

            write(u, imgfmt, advance='no') i, neb_e(i), neb_e(i) - neb_e(1), &
               & neb_e(i) - emin, ks(i), neb_dist(i), maxval(abs(neb_f(:,i) - ft*tau(:, i))), kd, maxval(abs(kd*tau(:, i)))

            neb_f(:, i) = neb_f(:, i) + (kd - ft)*tau(:, i)

!             if (.not. quiet) then
            fmax = maxval(abs(neb_f(:,i)))
            if (fmax > fmaxtot) then
               fmaxtot = fmax
               fmax_idx = i
            end if
            write(u, '(x,f22.16)') fmax
!             end if
         end do

         i = nimg
         kd = ks(i+1)*neb_dist(i+1) - ks(i)*neb_dist(i)
         write(u, imgfmt, advance='no') i, neb_e(i), neb_e(i) - neb_e(1), &
            & neb_e(i) - emin, ks(i), neb_dist(i), 0.0_dp, kd, 0.0_dp

         write(u,"('fmaxtot: ', f22.16, ' on image: ',i3,/)") fmaxtot, fmax_idx


         if ((fmaxtot < ftol) .or. (it > mxiter) .or. (usrexit() == 1)) exit



         if (rlxflg == 1) then

            do i = 2, nimg-1
               do j = 1, nx
                  disp = alph(j) * neb_f(j, i)
                  if (abs(disp) > dispmax) disp = sign(dispmax, disp)
                  neb_ad(j, i) = neb_ad(j, i) + disp
               end do
            end do

         elseif (rlxflg == 4) then
!          Fire for NEB

            velneb(:, 2:nimg-1) = velneb(:, 2:nimg-1) + 0.5_dp * dt * neb_f(:, 2:nimg-1)

!             vf = ddot(nx*(nimg-2), velneb(1,2), 1, neb_f(1,2), 1)
            vf = 0.0_dp
            vg_dot_vg = 0.0_dp
            fg_dot_fg = 0.0_dp
            max_v = 0.0_dp
            max_a = 0.0_dp

            do i=2, nimg-1
               do j = 1, nx
                  vf = vf + velneb(j,i)*neb_f(j,i)
                  cur_v      = velneb(j,i)*velneb(j,i)
                  vg_dot_vg  = vg_dot_vg + cur_v
                  cur_a      = neb_f(j,i)*neb_f(j,i)
                  fg_dot_fg  = fg_dot_fg + cur_a
                  if (cur_v > max_v) max_v = cur_v
                  if (cur_a > max_a) max_a = cur_a
               enddo
            enddo

            max_v = sqrt(max_v)
            max_a = sqrt(max_a)

!         Cut the velocities if the total power done by forces is negative.

            if( vf <= 0.0_dp ) then
               write (9, *)  "vf: setting all velocities to zero "
               velneb(:, 2:nimg-1) = 0.0_dp
               cut   = it
               dt    = dt*decfac
               mix   = mix_in
               cuts  = cuts + 1
            else
!         Otherwise mix the velocity and force.

               if (v_mix /= 0) then
                  help = mix*sqrt(vg_dot_vg/fg_dot_fg)
                  do i = 2, nimg-1
                     do j = 1, nx
                        velneb(j,i) = (1-mix)*velneb(j,i) + help*neb_f(j,i)
                     enddo
                  enddo
               endif

!         Then adjust the maximum achievable timestep

               loc = max_v/(2*max_a)
               maxdt = loc + sqrt(loc*loc + maxdr/max_a)

               if (maxdt > 0.2_dp) maxdt = 0.2_dp

               if ( (it-cut) > minsteps ) then
                  dt = min(dt*incfac, maxdt)
                  mix = mix*mixdec
               end if

            end if

            write(9,'("vf=",f12.8,"  dt=",f12.8,"  maxdt=",f12.8)') vf,dt,maxdt

            velneb(:, 2:nimg-1) = velneb(:, 2:nimg-1) + 0.5_dp * dt * neb_f(:, 2:nimg-1)


!         Move the atoms

            do i = 2, nimg-1
               do j = 1, nx
                  disp = dt*velneb(j,i)
                  if (abs(disp) > dispmax) disp = sign(dispmax,disp)
                  neb_ad(j,i) = neb_ad(j,i) + disp
               enddo
            enddo

         else
            write(*,'("relax flag",i2," not implemented for NEB")') rlxflg
            call panic()
         endif

         it = it + 1
      end do


      deallocate(neb_ad, neb_adinert, neb_f, neb_e, neb_de, neb_mg, neb_demi, tau, alph, ks, neb_dist)
      if (allocated(velneb)) deallocate(velneb)

      call close_unit(u)

   end subroutine nebf


   subroutine exec_img(img, it)

      integer, intent(in) :: img, it

      call set_img(img)

      call getetot(1)

      call get_img(img)

      if (.not. quiet) call print_img(img, it)

   end subroutine exec_img


   subroutine set_img(img)

      integer, intent(in) :: img

      include "Include/Atom.array"
      include "Include/PosVel.array"


      ad(1:3, 1:nd) = reshape(neb_ad(1:nx, img), [3, nd])

      de(1:nd) = neb_de(1:nd, img)
      adinert(1:3, 1:ninert) = reshape(neb_adinert(1:3*ninert, img), [3, ninert])

      if (mag) then
         mg(1:nd) = neb_mg(1:nd, img)
         demi(1:ninert) = neb_demi(1:ninert, img)
      end if

      current_img = img

   end subroutine set_img


   subroutine get_img(img)

      integer, intent(in) :: img

      include "Include/Atom.array"
      include "Include/PosVel.array"
      include "Include/Force.array"


      neb_de(1:nd, img) = de(1:nd)
      if (mag) neb_mg(1:nd, img) = mg(1:nd)

      neb_f(1:nx, img) = reshape(ftot(1:3, 1:nd), shape(neb_f(1:nx, img)))
      neb_e(img) = ebind

   end subroutine get_img



   subroutine init()

      use mod_conf, only : cell_t, read_cell, rtc


      integer :: ia, i, j, c
      real(dp) :: frac, dx

      type(ifl_t) :: cell_fin_fl
      type(cell_t) :: ucell_fin

      real(dp) :: a_fin(3,3), lena_fin(3)!, loc, df,df2
      real(dp), allocatable :: ndf(:)

      include "Include/Atom.array"
      include "Include/PosVel.array"


      nimg = rtc % neb_conf % nimg
      climb = rtc % neb_conf % climb
      spring_k = rtc % neb_conf % spring_k

!       u = get_new_unit()
!       open(u, file='neb.log')

      allocate(neb_ad(nx, nimg), neb_f(nx, nimg), neb_de(nd, nimg), neb_adinert(3*ninert, nimg), &
                              &  neb_e(nimg), alph(nx), ks(nimg), tau(nx, 2:nimg-1), neb_dist(nimg))

      alph = step

      if (mag) allocate(neb_mg(nd, nimg), neb_demi(ninert, nimg))

      neb_ad(1:nx, 1) = reshape(ad(1:3, 1:nd), shape(neb_ad(1:nx, 1)))
      neb_de(1:nd, 1) = de(1:nd)
      neb_adinert(1:3*ninert, 1) = reshape( adinert(1:3, 1:ninert), shape(neb_adinert(1:3*ninert, 1)))

      call fl_to_ifl('cell_neb_fin.in', cell_fin_fl)
      
      if (rtc % pot_flg == 9 .or. rtc % pot_flg == 10) then
          call read_cell(ucell_fin, cell_fin_fl,.true.)
      else
          call read_cell(ucell_fin, cell_fin_fl,.false.)
      endif
          
          
      a_fin = ucell_fin % a
      lena_fin =  ucell_fin % lena

      do i = 1,3
         a_fin(i,1:3) = a_fin(i,1:3)*lena_fin(i)/sqrt(sum(a_fin(i,1:3)*a_fin(i,1:3)))
      enddo

      do ia = 1, nd
         d(1:3, ia) = ucell_fin % d(ia) % crds(1:3)
         neb_de(ia, nimg) = ucell_fin % d(ia) % de
      end do

      do ia = 1, ninert
         dinert(1:3, ia) = ucell_fin % dinert(ia) % crds(1:3)
      end do

      do ia = 1, nd
         ad(1,ia) = a_fin(1,1)*d(1,ia)+a_fin(2,1)*d(2,ia)+a_fin(3,1)*d(3,ia)
         ad(2,ia) = a_fin(1,2)*d(1,ia)+a_fin(2,2)*d(2,ia)+a_fin(3,2)*d(3,ia)
         ad(3,ia) = a_fin(1,3)*d(1,ia)+a_fin(2,3)*d(2,ia)+a_fin(3,3)*d(3,ia)
      enddo

      neb_ad(1:nx, nimg) = reshape(ad(1:3, 1:nd), shape(neb_ad(1:nx, nimg)))

      do ia = 1, ninert
         adinert(1,ia) = a_fin(1,1)*dinert(1,ia)+a_fin(2,1)*dinert(2,ia)+a_fin(3,1)*dinert(3,ia)
         adinert(2,ia) = a_fin(1,2)*dinert(1,ia)+a_fin(2,2)*dinert(2,ia)+a_fin(3,2)*dinert(3,ia)
         adinert(3,ia) = a_fin(1,3)*dinert(1,ia)+a_fin(2,3)*dinert(2,ia)+a_fin(3,3)*dinert(3,ia)
      end do

      neb_adinert(1:3*ninert, nimg) = reshape(adinert(1:3, 1:ninert), shape(neb_adinert(1:3*ninert, nimg)))


      do c = 0, 2
         if ( rtc%rlim(c+1) > 0) then
            do j = 1, nx, 3
               dx = neb_ad(j+c, nimg) - neb_ad(j+c, 1)
               if (abs(dx) > 0.5_dp*lena(c+1)) neb_ad(j+c, nimg) = neb_ad(j+c, nimg) - sign(lena(c+1), dx)
            end do
            do j = 1, 3*ninert, 3
               dx = neb_adinert(j+c, nimg) - neb_adinert(j+c, 1)
               if (abs(dx) > 0.5_dp*lena(c+1)) neb_adinert(j+c, nimg) = neb_adinert(j+c, nimg) - sign(lena(c+1), dx)
            end do
         end if
      end do


      if (mag) then
         neb_mg(1:nd, 1) = mg(1:nd)
         neb_demi(1:ninert, 1) = demi(1:ninert)

         do ia = 1, nd
            neb_mg(ia, nimg) = ucell_fin % d(ia) % mg
         end do

         do ia = 1, ninert
            neb_demi(ia, nimg) = 0.5_dp * ucell_fin % dinert(ia) % mg * istn(zinert(ia))
         end do
      end if


      do i = 2, nimg-1
         frac = real(i-1, dp) / real(nimg-1, dp)

         neb_ad(1:nx, i) = neb_ad(1:nx, 1)*(1.0_dp-frac) + neb_ad(1:nx, nimg)*frac
         neb_adinert(1:3*ninert, i) = neb_adinert(1:3*ninert,1)*(1.0-frac) + neb_adinert(1:3*ninert,nimg)*frac

         neb_de(1:nd, i) = neb_de(1:nd, 1)*(1.0_dp-frac) + neb_de(1:nd, nimg)*frac
      end do




!       call ieee_set_flag(ieee_overflow, .true.)
!       call ieee_set_halting_mode(ieee_overflow, .true.)
!
!       do i = 2, nimg
!          do j = 1, nx
!             write(8000+iproc, '(3(x,g0))') i, j, neb_ad(j,i) - neb_ad(j,i-1)
!          end do
!       end do
!
!       allocate(ndf(nx))
!
!       do i = 2, nimg
! !             write(7000+iproc, '(g0)') i
! !             loc = 0.0_dp
! !             do j = 1, nx
! !                df = neb_ad(j,i) - neb_ad(j,i-1)
! ! !                df2 = df*df
! !                loc = loc + df*df
! ! !                write(7000+iproc, '(2x,4(x,g0))') i, j, df, df2, loc
! !             end do
! !             neb_dist(i) = sqrt(loc)
!             ndf = neb_ad(:,i) - neb_ad(:,i-1)
!
!             neb_dist(i) = sqrt(ddot(nx, ndf, 1, ndf, 1))
!
! !             scalxy = 0.0_dp
! !             do k = 1, nx, 3
! !                scalxy = scalxy + (neb_ad(k,i)-neb_ad(k,i-1))**2 + (neb_ad(k+1,i)-neb_ad(k+1,i-1))**2
! !             enddo
! !             distxy(i) = sqrt(scalxy)
!
!             write(6000+iproc, '(3(x,g0))') i, loc, neb_dist(i)
!          enddo
!
!
!
!
!
!
!
!
! #ifdef MPI
!          call mpi_barrier(mpi_comm_world, ierr)
!          call mpi_finalize(i)
! #endif
!          call panic()
! !
!
!
!





      if (mag) then
         do i = 2, nimg - 1
            frac = real(i-1, dp) / real(nimg-1, dp)
            neb_mg(1:nd, i) = neb_mg(1:nd, 1)*(1.0_dp-frac) + neb_mg(1:nd, nimg)*frac
            neb_demi(1:ninert, i) = neb_demi(1:ninert,1)*(1.0_dp-frac) + neb_demi(1:ninert,nimg)*frac
         end do
      end if

      ks = spring_k

      if (mxiter < 1) then
         do i = 1, nimg
            call print_img(i, 0)
         end do
      end if
!
   end subroutine init

   subroutine print_img(img, it)

      integer, intent(in) :: img, it
      character(len=100) :: filename

      if (img /= current_img) call set_img(img)

      if (quiet) return

      write(filename,'("neb.it",i6.6,".img",i2.2)') it, img

      call writeddp('ddplot/'//trim(filename)//'.dat')
      call writexyz('xyz/'//trim(filename)//'.xyz', filename)
      call writecfg('cfg/'//trim(filename)//'.cfg')
      call writexbs('bs/'//trim(filename)//'.bs')
      call writecell('cell/'//trim(filename)//'.cell.in')

   end subroutine print_img

   function inclmod(a,b) result (r)

      real(dp), intent(in) :: a, b
      real(dp) :: t, r, s
! plot it if it is not clear

      t = a/b
      s = sign(1.0_dp, t)
      t = t + 0.5_dp*s
      r = aint(t)

      if (t == r) r = r - s

      r = a - r*b

! Alternative implementation. Should give the same result
!       s = sign(1.0_dp, a)
!       r = -s*(modulo(-s*a+0.5_dp*b, b) - 0.5_dp*b)

   end function inclmod


   end module mod_neb
