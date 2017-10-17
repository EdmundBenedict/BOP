

   module mod_kspace
      use mod_io
      use mod_precision
      use mod_const
      use mod_all_scalar
      use mod_conf
      use topologia, only : mpmap, iproc, aproc, sdistrib, nproc_kspace => nproc
      use mod_atom_ar
      implicit none

      private

      integer :: nk, knh


   !*** The derivative of the Hamiltonian.
      real(dp) :: kgrad(3,mxnstat,mxnstat)

      integer, allocatable :: khpos(:)

   !*** The K points, and their weights.
      real(dp), allocatable :: kk(:,:)
      real(dp), allocatable :: wtk(:)

   !*** The occupancies
      real(dp), allocatable :: occ(:,:,:)

   !*** Eigenvalues
      real(dp), allocatable :: enk(:,:,:)

   !*** Eigenvectors
      complex(dp), allocatable :: kpsi(:,:,:,:), kham(:,:,:), ccw(:,:)
      complex(dp), allocatable :: kovl(:,:,:)

      real(dp), allocatable :: krhod(:,:)

      integer, allocatable :: ateqv(:)
      logical :: symon


      public :: kspace_setup, kspace_destroy, nk, knh, kgrad, khpos, kk, wtk, occ, enk, kpsi, kham, krhod, ccw, ateqv, symon, kovl

   contains
   subroutine kspace_setup()

      integer :: ia

      include "../Include/Atom.array"


      call load_kpoints()

      allocate(khpos(nd+1))

      khpos(1) = 0
      do ia = 2,nd+1
         khpos(ia) = khpos(ia-1) + nstt(z(ia-1))
      enddo
      knh = khpos(nd+1)

      if (.not. quiet) print *, 'khpos:', khpos
      if (.not. quiet) print *, 'khn:', knh

      allocate(occ(knh, mpmap(iproc)+1:mpmap(iproc+1), nsp), &
             & enk(knh, mpmap(iproc)+1:mpmap(iproc+1), nsp))

      allocate( kpsi(knh,knh,mpmap(iproc)+1:mpmap(iproc+1), nsp), &
              & kham(knh,knh,mpmap(iproc)+1:mpmap(iproc+1)), &
              & krhod(knh,nsp)) !,ccw(knh,knh))
              
      allocate(kovl(knh,knh,mpmap(iproc)+1:mpmap(iproc+1)))

!          These shoul be really temporarily here

      if (.not. allocated(dq) ) then
         allocate(dq(nd))
      else
         if (nd /= size(dq)) then
            deallocate(dq)
            allocate(dq(nd))
         end if
      end if


      if (.not. allocated(mass) ) then
         allocate(mass(nd))
      else
         if (nd /= size(mass)) then
            deallocate(mass)
            allocate(mass(nd))
         end if
      end if


      if (mag) then
         if (.not. allocated(mg) ) then
            allocate(mg(nd))
         else
            if (nd /= size(mg)) then
               deallocate(mg)
               allocate(mg(nd))
            end if
         end if
         do ia = 1, nd
            mg(ia) = rtc % cell % d(ia) % mg
         enddo

!          if (.not. allocated(dem) ) then
!             allocate(dem(nd))
!          else
!             if (nd /= size(dem)) then
!                deallocate(dem)
!                allocate(dem(nd))
!             end if
!          end if

      end if

   end subroutine kspace_setup


   subroutine kspace_destroy()
      deallocate(khpos,kk,wtk,occ,enk,kpsi,kham,krhod,ateqv)
   end subroutine kspace_destroy


   subroutine load_kpoints()
#ifndef NOSYM
      use spglib_f08
#endif
      use mod_clock

      integer :: i, j, k, ikx, iky, ikz, ik, irk

      real(dp) :: pvects(3,3), invmesh(3), invnk, symprec
      real(dp), allocatable :: lkk(:,:), crds(:,:)
      integer, allocatable :: atypes(:), kki(:,:), kkmap(:), kmltp(:)
      integer :: ku, intshift(3), kmode, time_rev, nne, nkt

      character(len=20) :: international_sym, schoenflies_sym
      integer :: space_group

      type(cell_t), pointer :: cell
      type(kspace_conf_t), pointer :: conf
#ifndef NOSYM
      type(SpglibDataset) :: symset
#endif
      integer(8) :: c1, c2

      call system_clock(c1)

      if (.not. associated(rtc % cell % kspace_conf)) stop 'kspace config not associated! exitting..'

      cell => rtc % cell
      conf => cell % kspace_conf

      kmode = conf % kmode
      kbase = conf % kbase
      symon = conf % symon


      if (kmode == 0 .or. symon) then
         do i=1,3
            pvects(i,:) = cell % a(i,:) * cell % lena(i)/sqrt(sum(cell%a(i,:)*cell%a(i,:)))
         end do
      end if

      space_group = 0

#ifndef NOSYM
      if (symon) then

         symprec = 0.00001_dp

         nne = cell % nd + cell % ninert

         allocate(ateqv(nne),crds(3,nne),atypes(nne))

         do i=1,cell%nd
            crds(1:3,i) = modulo(cell%d(i)%crds(1:3),1.0_dp)
         end do

         do i=1,cell%ninert
            crds(1:3,cell%nd+i) = modulo(cell%dinert(i)%crds(1:3),1.0_dp)
         end do


         do j=1,cell%nd
            i = 0
            do while (i /= rtc % ham_conf % na .and. rtc % ham_conf % a (i) % symb /= cell%d(j)%symb)
               i = i + 1
            end do
            atypes(j) = i
         end do

         do j=1,cell%ninert
            i = 0
            do while (i /= rtc % ham_conf % na .and. rtc % ham_conf % a (i) % symb /= cell%dinert(j)%symb)
               i = i + 1
            end do
            atypes(cell%nd+j) = i+1
         end do

         symset = spg_get_dataset(pvects, crds, atypes, nne, symprec)
         space_group = symset % spacegroup_number
         if (space_group /= 0) ateqv = symset % equivalent_atoms + 1
         if (space_group /= 0 .and. .not. quiet) print *, 'ateqv:', ateqv
      end if
#endif



      if (kmode == 0) then

         nkt = product(conf % vnkpts)
         nk = nkt
         invmesh = 1.0_dp/real(2*conf % vnkpts, dp)
         invnk = 1.0_dp/real(nk, dp)

         if (space_group /= 0) then
#ifndef NOSYM
   !   kkmap as obtained from spg_get_ir_reciprocal_mesh maps indices starting at 0,
   !    to avoid having to add 1 the arrays are allocated from 0 here
            allocate(kki(3,0:nk-1), kkmap(0:nk-1))

            intshift = 0
            where ( conf % offcentre ) intshift = 1
            time_rev = 1


            conf % nk = spg_get_ir_reciprocal_mesh(kki, kkmap, conf % vnkpts, &
                     & intshift, time_rev, pvects, crds, atypes, nne, symprec)

            deallocate(crds, atypes)

            allocate( conf % kpts(3, conf % nk),  conf % wtk(conf % nk), kmltp(0:nk-1))

! find the multiplicity then the weights
            kmltp = 0
            do ik = 0, nk-1
               kmltp(kkmap(ik)) = kmltp(kkmap(ik)) + 1
            end do


            irk = 0
            do ik = 0, nk-1
               if ( ik == kkmap(ik) ) then
                  irk = irk + 1
                  conf % kpts(1:3, irk) = (intshift + 2*kki(1:3, ik))*invmesh
                  conf % wtk(irk) =  kmltp(ik) * invnk          !count(kkmap == ik)
               end if
            end do

            if (irk /= conf % nk) stop 'Inconsistency in the number of irreducible points! probably a bug in the symmetry library'

!             print *, 'mul:', kmltp
            deallocate(kkmap, kki, kmltp)

            nk = conf % nk
#else
            stop 'The symmetry finder is disabled at compile time! Exitting..'
#endif
         else !either symmetry not found or not enabled
            conf % nk = nk

            allocate( conf % kpts(3, nk),  conf % wtk(nk))

            conf % wtk = invnk

            do k = 0, conf % vnkpts(3) - 1
               ikz = k
               if (ikz >= conf % vnkpts(3)/2) ikz = ikz - conf % vnkpts(3)
               do j = 0, conf % vnkpts(2) -1
                  iky = j
                  if (iky >= conf % vnkpts(2)/2) iky = iky - conf % vnkpts(2)
                  do i = 0, conf % vnkpts(1) -1
                     ik = 1 + i + j * conf % vnkpts(1) + k * conf % vnkpts(2)
                     ikx = i
                     if (ikx >= conf % vnkpts(1)/2) ikx = ikx - conf % vnkpts(1)
                     conf % kpts(1:3, ik) = (intshift + 2*[ikx,iky,ikz]) * invmesh
                  end do
               end do
            end do
         end if
      end if





      call sdistrib(nk,nproc_kspace,mpmap,aproc)


      allocate(kk(3,mpmap(iproc)+1:mpmap(iproc+1)), wtk(mpmap(iproc)+1:mpmap(iproc+1)))


      kk  = conf % kpts(1:3, mpmap(iproc)+1:mpmap(iproc+1))
      wtk = conf % wtk(mpmap(iproc)+1:mpmap(iproc+1))

      if (.not. quiet .and. kmode == 0) then
         print *, 'tb'
         print *, 'aproc:', aproc
         print *, 'mpmap:', mpmap

         write(6,'("nk:",x,i0)') nk
         do i=mpmap(iproc)+1,mpmap(iproc+1)
               write(6,'(i5,x,3(x,f9.6),2x,i2,2x,f9.6)') i, kk(:,i), nint(nkt*wtk(i)), wtk(i)
         end do
!          ku = get_new_unit()
!          open(ku,file='qpts.out',action='write')
!          write(ku,'("   nkp=",i0,"   nkabc=",3(x,i0),"  lshft=",3(x,i0))') nk, conf % vnkpts, intshift
!          do i = 1, nk
!                write(ku,'(i5,4(x,e19.12e2))') i, conf % kpts(:,i),  conf % wtk(i)
!          end do
!          call close_unit(ku)
      end if



      if ( kbase == 0 .or. kmode == 0) then
         if (.not. quiet) write(6,'("Transforming k points to cartesian coords")')

         call inv3x3(pvects)
         do i = 1,3
         print *, pvects(i,:)
         end do
         allocate(lkk(3,mpmap(iproc)+1:mpmap(iproc+1)))
         lkk = kk
         call dgemm('n','n',3,size(kk,2),3,2*pi,pvects,3,lkk,3,0.0_dp,kk,3)
!             kk(:,:) = matmul(2*pi*pvects(:,:), lkk(:,:))
         deallocate(lkk)
         conf%kpts(:,mpmap(iproc)+1:mpmap(iproc+1)) = kk(:,:)
      else
         if (.not. quiet) write(6,'("k point already in cartesian coords")')
      end if


      if (.not. quiet) then
         write(6,'("nk:",x,i0)') nk
         do i=mpmap(iproc)+1,mpmap(iproc+1)
               write(6,'(i5,x,3(x,f9.6),2x,i2,2x,f9.6)') i, kk(:,i), nint(nkt*wtk(i)), wtk(i)
         end do
!          ku = get_new_unit()
!          open(ku,file='qptsc.out',action='write')
!          write(ku,'("   nkp=",i0,"   nkabc=",3(x,i0),"  lshft=",3(x,i0))') nk, conf % vnkpts, intshift
!          do i = 1, nk
!                write(ku,'(i5,4(x,e19.12e2))') i, conf % kpts(:,i),  2*conf % wtk(i)
!          end do
!          call close_unit(ku)
      end if

      call system_clock(c2)

      if (.not. quiet) print *, 't(kpts+symms):', real(c2-c1,dp)/real(cr,dp)


   end subroutine load_kpoints



   end module mod_kspace























