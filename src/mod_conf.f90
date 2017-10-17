   module mod_conf
   
   
!  
   use mod_precision
   use mod_io, only : ifl_t, fout, ferr, goto_next_uncmt, read_m
   use mod_tail, only : tail_t
   use mod_pft, only : pft_t

   
   
   implicit none
   
   
   type :: atom_b_t
      logical :: lorbs(0:2) = .false.
      integer :: norbs = 0, orblist(3) = 0
      real(dp) :: cc = 0.0_dp, es = 0.0_dp, ep = 0.0_dp, ed = 0.0_dp, st = 0.0_dp, xtra_coefs(4) = 0.0_dp
   end type atom_b_t
   
         
   type :: atom_e_t
      real(dp) :: lam_0 = 0.0_dp, m = 0.0_dp, r_core = 0.0_dp
   end type atom_e_t
   
   type :: atom_conf_t
      integer :: id
      character(len=2) :: symb
      integer :: z = 0
      real(dp) :: ravs = 0.0_dp, col(3) = 0.0_dp, qpol(10), uh = 0.0_dp
      type(atom_b_t) :: b 
      type(atom_e_t) :: e 
   end type atom_conf_t
   
   type :: atom_t
      real(dp) :: crds(3) = 0.0_dp, de = 0.0_dp, mg = 0.0_dp
      character(len=2) :: symb
!         integer :: z
!         type(atom_b_t), pointer :: b_conf => null()
!         type(atom_e_t), pointer :: e_conf => null()
   end type atom_t

   type :: hop_t
      real(dp) :: v = 0.0_dp, dtau = 0.0_dp, Oscr = 0.0_dp ! r0, d, n, nc
!         logical :: fit = .false.
      type(pft_t)  :: sc
      type(tail_t) :: tl
   end type hop_t
   
   type :: pwp_t
      type(pft_t)  :: fp        
      type(tail_t) :: tl
   end type pwp_t

   type :: env_t
      real(dp) :: a = 0.0_dp, c = 0.0_dp, nu = 0.0_dp, rnu = 0.0_dp
      type(tail_t) :: tl1, tl2
   end type env_t
   
   type :: bond_conf_t
      integer :: id
      character(len=4) :: name
      type(atom_conf_t), pointer :: atom0 => null(), &
                                 & atom1 => null()
      type(hop_t), pointer :: h(:) => null(),o(:) => null()
      type(hop_t), pointer :: sss => null(), sps => null(), pps => null(), ppp => null(), sds => null(), &
                           & pds => null(), pdp => null(), dds => null(), ddp => null(), ddd => null()
      type(bond_conf_t), pointer :: rev => null()
      type(env_t) :: env
      type(pwp_t) :: pwp
   end type bond_conf_t

   type :: ham_conf_t
      integer :: na = 0, nb = 0
      logical :: ovl = .false., sctb
      type(atom_conf_t), allocatable :: a(:)
      type(bond_conf_t), allocatable :: b(:)
   end type ham_conf_t
   
   type :: energy_volume_t
      real(dp) :: volrange(3) = 0.0_dp, ratiorange(3) = 0.0_dp
   end type energy_volume_t
   
   type :: energy_surface_t
      real(dp) :: esurfrange(3) = 0.0_dp
   end type energy_surface_t
   
   type :: relaxation_t
      integer :: writeper = 0, datform = 0, reltype = 0, rlxflg = 0, cnst_n = 0, mxiter = 0
      integer :: autosave = 0        !logical :: autosave
      real(dp) :: ok = 0.0_dp, ftol = 0.0_dp, step = 0.0_dp
      real(dp), allocatable :: cnst_v(:,:)
      integer, allocatable :: cnst_a(:)
   end type relaxation_t
   
   type :: dislocation_t
      integer :: gfbcon
      real(dp) :: radiusi, plotlims(4)
   end type dislocation_t

   type :: md_t
      real(dp) :: trace = 0.0_dp, dt = 0.0_dp, ttol = 0.0_dp, nequil = 0.0_dp, temp = 0.0_dp, &
               & mda_gr_r = 0.0_dp, mda_g3_r = 0.0_dp, mda_sw_e = 0.0_dp
      integer :: monitorper = 0, &
               & mda_gr = 0, mda_gr_n = 0, &
               & mda_g3 = 0, mda_g3_n = 0, &
               & mda_sw = 0, mda_sw_n = 0, &
               & mda_adav = 0, mda_op = 0
   end type md_t
   
   type :: neb_t
      integer :: nimg = 0, climb = 0
      real(dp) :: spring_k = 0.0_dp      
   end type neb_t
   
   type :: gamma_surface_t
      integer :: indrel = 0, gamxnp = 0, gamynp = 0, sym = 0
      real(dp) :: gamxmin = 0.0_dp, gamxmax = 0.0_dp, gamymin = 0.0_dp, gamymax = 0.0_dp
   end type gamma_surface_t

   type :: grain_boundary_t
      integer :: gbmoveopt = 0
      real(dp) :: xshiftgb = 0.0_dp, yshiftgb = 0.0_dp, zshiftgb = 0.0_dp
   end type grain_boundary_t

   type :: elastic_const_contribs_t
      real(dp) :: tot = 0.0_dp, elc = 0.0_dp, env = 0.0_dp, pwp = 0.0_dp
   end type elastic_const_contribs_t

   type :: elastic_const_t
      integer :: nt = 0, polyord = 0
      real(dp) :: tmatrix(3,3) = 0.0_dp, tlo = 0.0_dp, thi = 0.0_dp
      type(elastic_const_contribs_t) :: s11, s33, c11, c33, c44, c66, c12, c13, c1266, c1344

   end type elastic_const_t

   type :: bop_conf_t
      integer :: momflg = 0, nrec = 0, nbase = 0, term=0,v1flag=0,chi_meth=0,mfac=0,mxit=0,mgit=0,nit=-1
      real(dp) :: kt = 0.0_dp, etail = 0.0_dp, qerr = 0.0_dp, merr = 0.0_dp, eftol = 0.0_dp,b=0.0_dp,wc=0.0_dp
      logical :: mag = .false.
   end type bop_conf_t
   
   type :: kspace_conf_t
      integer :: kmode = 0, vnkpts(3) = 0, kbase = 0, nk = 0
      logical :: offcentre(3) = .false., symon = .true.
      real(dp), allocatable :: kpts(:,:), wtk(:)
   end type kspace_conf_t
      

   type :: ens_t
   real(dp) :: bind = 0.0_dp, elct = 0.0_dp, band = 0.0_dp, prom = 0.0_dp, bond = 0.0_dp, & 
               & atom = 0.0_dp, ent  = 0.0_dp, clas = 0.0_dp, pair = 0.0_dp, env  = 0.0_dp, &
               & mag  = 0.0_dp
   end type ens_t

   type :: cell_t
      character(len=10) :: name, str
      real(dp) :: a(3,3), lena(3), latpar, nullheight, x_center, y_center, r_min, r_max, e_coh, log_rad, eflo
      integer :: nd = 0, nbp =-1, ninert = 0, n_mesh = 0, iter = 0
      integer :: if_relaxed = 0
!         logical :: if_relaxed
      type(atom_t), allocatable :: d(:), dinert(:), unrld(:)
      type(kspace_conf_t), pointer :: kspace_conf => null()
      type(dislocation_t), pointer :: dislocation_conf => null()
   end type cell_t
   
   type :: run_conf_t
      character(len=10) :: genfile
      type(cell_t), pointer :: cell => null()
      integer :: pot_flg = 0, mvflag = 0
      integer :: scf = 0, writescale = 0, elplusrel = 0, rlim(3) = 0, &
               & fs = 0, env = 0, vpair = 0, fs_only = 0, emb = 0, &
               & pot_scl = 0, idebug = 0, restart = 0
      !logical :: scf, writescale, elplusrel, rlim(3), fs, env, vpair, fs_only, emb, pot_scl, idebug, restart
      real(dp) :: scf_cut = 0.0_dp, rcut(2) = 0.0_dp, rprune = 0.0_dp
      
      type(bop_conf_t), pointer :: bop_conf => null()
      type(ham_conf_t) :: ham_conf,tbham_conf
      
      type(energy_volume_t) , pointer :: energy_volume_conf  => null()
      type(energy_surface_t), pointer :: energy_surface_conf => null()
      type(gamma_surface_t) , pointer :: gamma_surface_conf  => null() 
      type(grain_boundary_t), pointer :: grain_boundary_conf => null()
      type(elastic_const_t) , pointer :: elastic_const_conf  => null()
      
      type(relaxation_t), pointer :: relaxation_conf => null()
      type(md_t) , pointer :: md_conf  => null()
      type(neb_t), pointer :: neb_conf => null()
      type(ens_t) :: tote, avre
      
      logical :: quiet, paral, forces
      
   end type run_conf_t
   
   interface print_ens
      module procedure print_ens_1, print_ens_a
   end interface
   
   type(run_conf_t), pointer :: rtc

   contains
   
   subroutine print_atom_b(atom_b, l, u)
      type(atom_b_t), intent(in) :: atom_b
      character(len=*), intent(in) :: l
      integer, intent(in) :: u
      
      write(u, '(a,x,f8.4,2x,i0,x, 3(x,i0),x,3(x,f12.6),x,x,f8.4,x,4(x,f12.6),a)') &
      & l, atom_b%cc, atom_b%norbs, atom_b%orblist, atom_b%es, atom_b%ep, atom_b%ed, atom_b%st, atom_b%xtra_coefs, &
      & '  !cc norbs orblist(3i) es ep ed st poly_coeffs(4r)'
      
   end subroutine print_atom_b
   
   
   subroutine print_atom_e(atom_e, l, u) 
      type(atom_e_t), intent(in) :: atom_e
      character(len=*), intent(in) :: l
      integer, intent(in) :: u
      write(u, '(a,3(x,f8.4),a)')l, atom_e%lam_0, atom_e%m, atom_e%r_core, '  !lam0 m r_core'
   end subroutine print_atom_e
   
   
   
   subroutine print_atom_conf(atom_conf, l, u,sctb)
      type(atom_conf_t), intent(in) :: atom_conf
      character(len=*), intent(in) :: l
      integer, intent(in) :: u
      logical :: sctb
      
      write(u, '(a,x,a,2x,i0,x,f12.6,x,3(x,f12.6),a)') &
      & l, trim(atom_conf%symb), atom_conf%z, atom_conf%ravs, atom_conf%col, '  !symb z ravs col(3r)'
      call print_atom_b(atom_conf%b, l//'    ', u)
      call print_atom_e(atom_conf%e, l//'    ', u)
      if (sctb) then
          write(u, '(a,10(x,f12.6),a)') l, atom_conf%qpol, '   !qpol'
          write(u, '(a,x,f12.6,a)') l, atom_conf%uh, '   !Uh'
      endif
      
   end subroutine print_atom_conf


   subroutine print_atom(atom,l,u)
      type(atom_t), intent(in) :: atom
      character(len=*), intent(in) :: l
      integer, intent(in) :: u
      
      write(u, '(a,x,a,x,3(x,f20.12),x,2(x,f20.12),a)') l, atom%symb, atom%crds, atom%de, atom%mg, '  !symb crds(3r) de mg'
      
   end subroutine print_atom

   subroutine print_tail(t, l, u, adv)
      type(tail_t), intent(in) :: t
      character(len=*), intent(in) :: l
      integer, intent(in) :: u
      logical, intent(in), optional :: adv
      logical :: la
      
      la = .true.
      if (present(adv)) la = adv
      
      if (la) then
         write(u, '(a,x,i0,x,4(x,f20.12),x,a)') l, t%t, t%r1, t%rc, t%n, t%m, '! tt r1 rcut n m'
      else
         write(u, '(a,x,i0,x,4(x,f20.12))',advance='no') l, t%t, t%r1, t%rc, t%n, t%m
      end if
   end subroutine print_tail

   subroutine print_pft(p,l,u,adv,shr)
      type(pft_t), intent(in) :: p
      character(len=*), intent(in) :: l
      integer, intent(in) :: u
      logical, intent(in), optional :: adv, shr
      logical :: la, sh
      
      sh = .false.
      la = .true.
      if (present(shr)) sh = shr
      if (present(adv)) la = adv
      
      write(u, '(a)',advance='no') l
      if (.not. sh) then
         write(u, '(x,i0)', advance='no') p%t
         write(u, '(x,i0)', advance='no') p%n
      end if      
      write(u, '(20(:,x,f20.12))', advance='no') p % a(:p%n)
      
      if (la) write(u, '("")')
      
   end subroutine print_pft
   
   
   subroutine print_hop(hop, l, u)
      type(hop_t), intent(in) :: hop
      character(len=*), intent(in) :: l
      integer, intent(in) :: u
      
      write(u, '(a,7(x,f20.12))',advance='no') l,hop%v, hop%dtau, hop%Oscr
      call print_pft (hop%sc,' ',u, adv=.false., shr=.true. )
      call print_tail(hop%tl,' ',u)        
      
   end subroutine print_hop
   
   
   subroutine print_pwp(pwp, l, u)
      type(pwp_t), intent(in) :: pwp
      character(len=*), intent(in) :: l
      integer, intent(in) :: u
      
      call print_pft (pwp%fp, l//'pwp', u)
      call print_tail(pwp%tl, l//'    ptl', u)
      
   end subroutine print_pwp
   
   subroutine print_env(e, l, u)
      type(env_t), intent(in) :: e
      character(len=*), intent(in) :: l
      integer, intent(in) :: u
   
      write(u, '(a,4(x,f20.12))') l//'env', e%a, e%c, e%nu, e%rnu
      call print_tail(e%tl1, l//'    et1', u) 
      call print_tail(e%tl2, l//'    et2', u) 
      
   end subroutine print_env
   
   subroutine print_bond_conf(bond_conf, l, u,ovl)
      type(bond_conf_t), intent(in) :: bond_conf
      character(len=*), intent(in) :: l
      integer, intent(in) :: u
      logical,  intent(in) :: ovl
      
      write(u,'(a,x,a)') l, bond_conf%name
      
      write(u,'(a,13x,a)') l, '! Vscr                 DTau                 Oscr                   '&
                           & //'r0                   rc                   n                    nc             '&
                           & //'tt        r1                   rcut                 n                    m'
      call print_hop(bond_conf%sss, l//'    sss', u)
      call print_hop(bond_conf%sps, l//'    sps', u)
      call print_hop(bond_conf%pps, l//'    pps', u)
      call print_hop(bond_conf%ppp, l//'    ppp', u)
      call print_hop(bond_conf%sds, l//'    sds', u)
      call print_hop(bond_conf%pds, l//'    pds', u)
      call print_hop(bond_conf%pdp, l//'    pdp', u)
      call print_hop(bond_conf%dds, l//'    dds', u)
      call print_hop(bond_conf%ddp, l//'    ddp', u)
      call print_hop(bond_conf%ddd, l//'    ddd', u)
      
      write(u,'(a,a,x,l,a)') l,'ovl', ovl 
      if (ovl) then 
          call print_hop(bond_conf%o(0), l//'    sss', u)
          call print_hop(bond_conf%o(1), l//'    sps', u)
          call print_hop(bond_conf%o(2), l//'    pps', u)
          call print_hop(bond_conf%o(3), l//'    ppp', u)
          call print_hop(bond_conf%o(4), l//'    sds', u)
          call print_hop(bond_conf%o(5), l//'    pds', u)
          call print_hop(bond_conf%o(6), l//'    pdp', u)
          call print_hop(bond_conf%o(7), l//'    dds', u)
          call print_hop(bond_conf%o(8), l//'    ddp', u)
          call print_hop(bond_conf%o(9), l//'    ddd', u)
      endif
      
      write(u,'(a,13x,a)') l, '! A                    C                    nu                   rnu' 
      call print_env(bond_conf%env , l//'    ', u)
      write(u,'(a,6x,a)') l, '!tp sz       args...' 
      call print_pwp(bond_conf%pwp, l//'    ', u)
      
      
   end subroutine print_bond_conf

   subroutine print_ham_conf(ham_conf, l, ham,u)
      type(ham_conf_t), intent(in) :: ham_conf
      character(len=*), intent(in) :: l
      integer, intent(in) :: u
      character(len=*), intent(in) :: ham
      
      integer :: i
      write(u, '(a,x,a)') l, ham
      write(u,'(a,l,a)') l , ham_conf % sctb, ' ! sctb'
      write(u, '(a,2(x,i0),a)') l, ham_conf%na, ham_conf%nb, ' ! natypes nbtypes'
      do i=0,ham_conf%na-1
         call print_atom_conf(ham_conf%a(i), l//'    ', u,ham_conf % sctb)
      end do
      
      do i = 0, ham_conf % nb - 1
         call print_bond_conf(ham_conf%b(i), l//'    ', u,ham_conf % ovl)
      end do        
   end subroutine print_ham_conf

   subroutine print_cell(cell, l, u)
      type(cell_t), intent(in) :: cell
      character(len=*), intent(in) :: l
      integer, intent(in) :: u
      
      integer :: i
      
      write(u, '(a,x,a,x,a,a)') l, cell%name, cell%str, '  !name str'
      
      do i=1,3
         write(u, '(a,3(x,f20.12))') l, cell%a(i,1:3)
      end do
      
      write(u, '(a)') l
      
      write(u, '(a,x,a,3(x,f20.12))') l, 'len', cell%lena 
      write(u, '(a)') l
      write(u, '(a,x,a,x,f20.12)') l, 'latpar', cell%latpar
      write(u, '(a)') l
      
      write(u, '(a,x,a,x,i0)') l, 'nd', cell%nd
      if (cell%nbp > -1) write(u, '(a,x,a,x,i0)') l, 'nbp', cell%nbp
      do i=1,cell%nd
         call print_atom(cell%d(i),l//'    ',u)
      end do
      
      write(u, '(a)') l
      
      if (cell%ninert>0) then
         write(u, '(a,x,a,x,i0)') l, 'ninert', cell%ninert
         do i=1,cell%ninert
            call print_atom(cell%dinert(i),l//'    ',u)
         end do
         write(u, '(a)') l
      end if
      
      if (allocated(cell%unrld)) then
         write(u, '(a,x,a)') l, 'unrld'
         do i=1,cell%nd
            call print_atom(cell%unrld(i),l//'    ',u)
         end do
         write(u, '(a)') l
      end if
      
      
      write(u,'(a,x,a,x,f20.12)') l, 'nullheight', cell%nullheight 
      write(u, '(a)') l
      
      if (associated(cell%kspace_conf)) call print_kspace_conf_t(cell%kspace_conf, l, u)
      write(u, '(a)') l
      
      if (associated(cell%dislocation_conf)) call print_dislocation_conf(cell%dislocation_conf, l, u)
      write(u, '(a)') l
      
   end subroutine print_cell

   subroutine print_kspace_conf_t(kspace_conf, l, u)

      type(kspace_conf_t), intent(in) :: kspace_conf
      character(len=*), intent(in) :: l
      integer, intent(in) :: u

      integer :: i

      write(u, '(a,x,a)') l, 'kspace_conf'
      write(u, '(a,x,i0,x,l,x,a)') l, kspace_conf % kmode, kspace_conf % symon, '!symmmode, symetry_enabled'

      select case(kspace_conf % kmode)
      case(0)
            write(u,'(a,x,3(x,i0),x,a)') l, kspace_conf % vnkpts, '!vnkpts'
            write(u,'(a,x,3(x,l),x,a)') l, kspace_conf % offcentre , '!offcentre'              
      case(1)

            write(u,'(a,x,i0,x,a)') l, kspace_conf % kbase, '!kbase'
            write(u,'(a,x,i0,x,a)') l, kspace_conf % nk, '!nk : kpts, wtk'
            do i=1,kspace_conf % nk
               write(u,'(a,x,3(x,f20.12),x,x,f20.12)') l, kspace_conf % kpts(:,i), kspace_conf % wtk(i)
            end do      
      end select

   end subroutine print_kspace_conf_t


   subroutine print_elastic_const_conf(ecc, l, u)
      type(elastic_const_t), intent(in) :: ecc
      character(len=*), intent(in) :: l
      integer, intent(in) :: u
      
      integer :: i,j
      write(u, '(a,a)') l,'elastic_const_conf'        
      do i=1,3
         write(u, '(a,3(x,f20.12))') l, (ecc%tmatrix(i,j), j=1,3)
      end do
      write(u, '(a,x,i0,a)') l,ecc % polyord, '  !polynomial order'
      write(u, '(a,2(x,f20.12),x,i0)') l,  ecc % tlo, ecc % thi, ecc % nt

   end subroutine print_elastic_const_conf


   subroutine print_elastic_consts(ecc, l, u)
      type(elastic_const_t), intent(in) :: ecc
      character(len=*), intent(in) :: l
      integer, intent(in) :: u
      
      write (u, '(a,x,a,/,10(a,x,a,4(x,f20.12),/))') &
      & l, '             total                electronstr          environmental        pairwise', &
      & l, 's11  ', ecc%s11%tot  , ecc%s11%elc  , ecc%s11%env  , ecc%s11%pwp  , & 
      & l, 's33  ', ecc%s33%tot  , ecc%s33%elc  , ecc%s33%env  , ecc%s33%pwp  , &
      & l, 'c11  ', ecc%c11%tot  , ecc%c11%elc  , ecc%c11%env  , ecc%c11%pwp  , &
      & l, 'c33  ', ecc%c33%tot  , ecc%c33%elc  , ecc%c33%env  , ecc%c33%pwp  , &
      & l, 'c12  ', ecc%c12%tot  , ecc%c12%elc  , ecc%c12%env  , ecc%c12%pwp  , &
      & l, 'c13  ', ecc%c13%tot  , ecc%c13%elc  , ecc%c13%env  , ecc%c13%pwp  , &
      & l, 'c44  ', ecc%c44%tot  , ecc%c44%elc  , ecc%c44%env  , ecc%c44%pwp  , &
      & l, 'c66  ', ecc%c66%tot  , ecc%c66%elc  , ecc%c66%env  , ecc%c66%pwp  , &
      & l, 'c1344', ecc%c1344%tot, ecc%c1344%elc, ecc%c1344%env, ecc%c1344%pwp, &
      & l, 'c1266', ecc%c1266%tot, ecc%c1266%elc, ecc%c1266%env, ecc%c1266%pwp


      write(u, '(a,x,a,/,(a,x,a,4(x,f20.12),/))') l, 'diffs', l, 'elc  ', &
      & ecc%c1344%tot-ecc%c1266%tot, ecc%c1344%elc-ecc%c1266%elc, &
      & ecc%c1344%env-ecc%c1266%env, ecc%c1344%pwp-ecc%c1266%pwp

   end subroutine print_elastic_consts

   
   

   subroutine print_bop_conf(bop_conf, l, u)
      type(bop_conf_t), intent(in) :: bop_conf
      character(len=*), intent(in) :: l
      integer, intent(in) :: u
      
      write(u, '(a,x,a)') l, 'bop_conf'
      write(u, '(a,7(x,i0),2(x,f20.12),a)') l, &
      & bop_conf%momflg, bop_conf%nrec, bop_conf%nbase, bop_conf%term, &
      & bop_conf%v1flag, bop_conf%chi_meth, bop_conf%mfac, bop_conf%kt, &
      & bop_conf%etail, ' ! momflg nrec nbase term v1flag chi_meth mfac kt etail'
      
      write(u, '(a,x,a)') l, 'sc_conf'
      if (bop_conf%nit == -1) then
        write(u, '(a,x,f20.12,x,i0,x,f20.12,x,l,x,f20.12,x,i0,a)') l, &
        & bop_conf%qerr, bop_conf%mxit, bop_conf%eftol, &
        & bop_conf%mag , bop_conf%merr, bop_conf%mgit, &
        & ' ! qerr mxit eftol mag merr mgit'
      else
        write(u, '(a,x,f20.12,x,i0,x,f20.12,x,l,x,f20.12,x,i0,x,i0,x,&
        &    f20.12,x,f20.12,a)') l, &
        & bop_conf%qerr, bop_conf%mxit, bop_conf%eftol, &
        & bop_conf%mag , bop_conf%merr, bop_conf%mgit,  &
        & bop_conf%nit , bop_conf%b, bop_conf%wc,&
        & ' ! qerr mxit eftol mag merr mgit nit b wc'
      endif
   
   
   end subroutine print_bop_conf


   subroutine print_run_conf(run_conf, l, u)
      type(run_conf_t), intent(in) :: run_conf
      character(len=*), intent(in) :: l
      integer, intent(in) :: u
      
      write(u, '(a,x,a)') l, run_conf%genfile
      
      write(u,'(a,2(x,i0),2(x,l),a)') l, run_conf%pot_flg, run_conf%mvflag, run_conf%quiet, run_conf%paral, &
         & '  !pot_flg mvflag quiet paral'
      write(u,'(a,3(x,i0),a)') l, run_conf%scf, run_conf%writescale, run_conf%elplusrel, '  !scf writescale elplusrel'
      write(u,'(a,3(x,i0),a)') l, run_conf%rlim, '  !rlim(3i)'
      write(u,'(a,6(x,i0),x,l,a)') l, run_conf%fs, run_conf%env, run_conf%vpair, run_conf%fs_only, run_conf%emb, &
         & run_conf%pot_scl, run_conf%forces, '  !fs env vpair fs_only emb pot_scl forces'
      write(u,'(a,2(x,f20.12),a)') l, run_conf%scf_cut, run_conf%rprune, '  !scf_cut rprune'
      write(u,'(a,2(x,f20.12),a)') l, run_conf%rcut, '  !rcut(2r)'
      write(u,'(a,2(x,i0),a)') l, run_conf%idebug, run_conf%restart, '  !idebug restart'
      
      if (associated(run_conf%bop_conf)) call print_bop_conf(run_conf%bop_conf, l//'    ',u)
      call print_ham_conf(run_conf%ham_conf, l//'    ','ham_conf', u)
      
      if  (allocated(run_conf%tbham_conf%a)) call print_ham_conf(run_conf%tbham_conf, l//'    ','tbham_conf', u)
      
      if (associated(run_conf%elastic_const_conf)) &
         & call print_elastic_const_conf(run_conf%elastic_const_conf, l//'    ',u)
      
      
      if (associated(run_conf%relaxation_conf)) call print_relaxation_conf(run_conf%relaxation_conf, l//'    ', u)
      
      if (associated(run_conf%neb_conf)) call print_neb_conf(run_conf%neb_conf, l//'    ', u)
      
      if (associated(run_conf%md_conf)) call print_md_t(run_conf%md_conf, l//'    ', u)
      
   end subroutine print_run_conf
      
   
   subroutine print_ens_1(ens, l, u)
      type(ens_t), intent(in) :: ens
      character(len=*), intent(in) :: l
      integer, intent(in) :: u
      
      write(u, '(11(a,x,a,x,f20.12,/))') &
      & l, 'ebond', ens%bond, &
      & l, 'eprom', ens%prom, &
      & l, 'eband', ens%band, &
      & l, 'eatom', ens%atom, &
      & l, 'emag ', ens%mag,  &
      & l, 'uent ', ens%ent,  &
      & l, 'eelct', ens%elct, &
      & l, 'eenv ', ens%env,  &
      & l, 'epair', ens%pair, &
      & l, 'eclas', ens%clas, &
      & l, 'ebind', ens%bind
                                    
         
      
   end subroutine print_ens_1
   
   subroutine print_ens_a(enst, ensa, l, u)
      type(ens_t), intent(in) :: enst, ensa
      character(len=*), intent(in) :: l
      integer, intent(in) :: u
      
      write(u, "(11(a,x,a,2(x,es20.12),/))") &
      & l, 'ebond', enst % bond, ensa % bond, &
      & l, 'eprom', enst % prom, ensa % prom, &
      & l, 'eband', enst % band, ensa % band, &
      & l, 'eatom', enst % atom, ensa % atom, &
      & l, 'emag ', enst % mag,  ensa % mag,  &
      & l, 'uent ', enst % ent,  ensa % ent,  &
      & l, 'eelct', enst % elct, ensa % elct, &
      & l, 'eenv ', enst % env,  ensa % env,  &
      & l, 'epair', enst % pair, ensa % pair, &
      & l, 'eclas', enst % clas, ensa % clas, &
      & l, 'ebind', enst % bind, ensa % bind
      
   end subroutine print_ens_a
   
   subroutine print_neb_conf(neb, l, u)
      type(neb_t), intent(in) :: neb
      character(len=*), intent(in) :: l
      integer, intent(in) :: u
      write(u, "(a,x,a)") l, 'neb_conf'
      write(u, "(a,x,2(x,i0),x,es20.12,x,a)") l, neb % nimg, neb % climb, neb % spring_k, '! nimg climb k'
   end subroutine print_neb_conf
   
   subroutine print_relaxation_conf(rlx, l, u)
      type(relaxation_t), intent(in) :: rlx
      character(len=*), intent(in) :: l
      integer, intent(in) :: u
      integer :: i
      write(u, "(a,x,a)") l, 'relaxation_conf'
      write(u, "(a,x,3(x,i0),2(x,es20.12),x,a)") l, rlx%reltype, rlx%rlxflg, rlx%mxiter, rlx%ftol, rlx%step, &
         & '! reltype rlxflg mxiter ftol step'
      write(u, '(a,x,3(x,i0),x,es20.12,x,a)') l, rlx%writeper, rlx%datform, rlx%autosave, rlx%ok, &
         & '! writeper datform autosave ok'
      write(u, '(a,x,i0,x,a)') l, rlx%cnst_n, '! cnst_n'
      do i = 1, rlx%cnst_n
         write(u, '(a,x,x,i0,x,3(x,es20.12))') rlx%cnst_a(i), rlx%cnst_v(:,i)
      end do 
   end subroutine print_relaxation_conf

   subroutine print_md_t(md, l, u)
      type(md_t), intent(in) :: md
      character(len=*), intent(in) :: l
      integer, intent(in) :: u
      
      write(u, '(a,x,a)') l, 'md_conf' 
      write(u, '(a,x,5(x,es20.12),3(x,i0),x,a)') l , md%trace, md%dt, md%ttol, md%nequil, md%temp, md%monitorper, &
         & md%mda_adav, md%mda_op, '! trace dt ttol nequil temp monitorper mda_adav mda_op'
      write(u, '(a,x,2(x,i0),x,es20.12,x,a)') l, md%mda_gr, md%mda_gr_n, md%mda_gr_r, '! mda_gr mda_gr_n mda_gr_r'
      write(u, '(a,x,2(x,i0),x,es20.12,x,a)') l, md%mda_g3, md%mda_g3_n, md%mda_g3_r, '! mda_g3 mda_g3_n mda_g3_r'
      write(u, '(a,x,2(x,i0),x,es20.12,x,a)') l, md%mda_sw, md%mda_sw_n, md%mda_sw_e, '! mda_sw mda_sw_n mda_sw_e'
   end subroutine print_md_t
   
   
   subroutine print_dislocation_conf(dsl, l, u)
      type(dislocation_t), intent(in) :: dsl
      character(len=*), intent(in) :: l
      integer, intent(in) :: u
      
      write(u, '(a,x,a)') l, 'dislocation_conf'
      write(u, '(a,x,x,i0,x,es20.12,x,a)') l, dsl%gfbcon, dsl%radiusi, ' ! gfbcon, radiusi'
      write(u, '(a,x,4(x,es20.12),x,a)') l, dsl%plotlims, '! plxmin plxmax plymin plymax'
   end subroutine print_dislocation_conf
   
   
   
!     reading




   subroutine read_md_t(md, u)
      type(md_t), intent(inout) :: md
      type(ifl_t), intent(inout) :: u
      character(len=15) :: nam
      
      call read_m(u, md%trace)
      call read_m(u, md%dt)
      call read_m(u, md%ttol)
      call read_m(u, md%nequil)
      call read_m(u, md%temp)
      call read_m(u, md%monitorper)
      call read_m(u, md%mda_adav)
      call read_m(u, md%mda_op)
      call read_m(u, md%mda_gr)
      call read_m(u, md%mda_gr_n)
      call read_m(u, md%mda_gr_r)
      call read_m(u, md%mda_g3)
      call read_m(u, md%mda_g3_n)
      call read_m(u, md%mda_g3_r)
      call read_m(u, md%mda_sw)
      call read_m(u, md%mda_sw_n)
      call read_m(u, md%mda_sw_e)
   end subroutine read_md_t

   subroutine read_neb_conf(neb, u)
      type(neb_t), intent(inout) :: neb
      type(ifl_t), intent(inout) :: u
      character(len=15) :: nam
      
      call read_m(u, neb % nimg)
      call read_m(u, neb % climb)
      call read_m(u, neb % spring_k)
   end subroutine read_neb_conf   
   
   subroutine read_relaxation_conf(rlx, u)
      type(relaxation_t), intent(inout) :: rlx
      type(ifl_t), intent(inout) :: u
      character(len=15) :: nam
      integer :: i
      
      call read_m(u, rlx%reltype)
      call read_m(u, rlx%rlxflg)
      call read_m(u, rlx%mxiter)
      call read_m(u, rlx%ftol)
      call read_m(u, rlx%step)
      call read_m(u, rlx%writeper)
      call read_m(u, rlx%datform)
      call read_m(u, rlx%autosave)
      call read_m(u, rlx%ok)
      call read_m(u, rlx%cnst_n)
      
      allocate(rlx%cnst_a(rlx%cnst_n), rlx%cnst_v(3,rlx%cnst_n))
      
      do i = 1, rlx%cnst_n
         call read_m(u, rlx%cnst_a(i))
         call read_m(u, rlx%cnst_v(:,i))
      end do
   end subroutine read_relaxation_conf
   
   
   subroutine read_atom_b(atom_b, u)
      type(atom_b_t), intent(inout) :: atom_b
      type(ifl_t), intent(inout) :: u
      
      call read_m(u, atom_b%cc)
      
      call read_m(u, atom_b%norbs)
      call read_m(u, atom_b%orblist)
      call read_m(u, atom_b%es)
      call read_m(u, atom_b%ep)
      call read_m(u, atom_b%ed)
      call read_m(u, atom_b%st)
      call read_m(u, atom_b%xtra_coefs)

      atom_b % lorbs(atom_b%orblist(:atom_b%norbs)) = .true.
   end subroutine read_atom_b
   
   
   subroutine read_atom_e(atom_e, u) 
      type(atom_e_t), intent(inout) :: atom_e
      type(ifl_t), intent(inout) :: u
      call read_m(u, atom_e%lam_0)
      call read_m(u, atom_e%m)
      call read_m(u, atom_e%r_core)
   end subroutine read_atom_e
   
   
   
   subroutine read_atom_conf(atom_conf, u,sctb)
      type(atom_conf_t), intent(inout) :: atom_conf
      type(ifl_t), intent(inout) :: u
      logical :: sctb
      
      call read_m(u, atom_conf%symb)
      call read_m(u, atom_conf%z)
      call read_m(u, atom_conf%ravs)
      call read_m(u, atom_conf%col)
      call read_atom_b(atom_conf%b, u)
      call read_atom_e(atom_conf%e, u)
      if (sctb) then
          call read_m(u,atom_conf%qpol)
          call read_m(u, atom_conf%uh)
      endif
   end subroutine read_atom_conf


   subroutine read_atom(atom, u)
      type(atom_t), intent(inout) :: atom
      type(ifl_t), intent(inout) :: u
      
      call read_m(u, atom%symb)
      call read_m(u, atom%crds)
      call read_m(u, atom%de)
      call read_m(u, atom%mg)
      
   end subroutine read_atom

   subroutine read_tail(t, u)
      type(tail_t), intent(inout) :: t
      type(ifl_t), intent(inout) :: u
      character(len=3) :: l
      
      call read_m(u, t % t)
      
!       select case(t%t)
!          case (1);  
!             t%f   => polytail
!          case (2);
!             t%mlt = .true.
!             t%f   => binmtail
!          case (3);   
!             t%mlt = .true.
!             t%f   => dstptail
!          case default
!          stop 'invalid tail choice'
!       end select
      
      call read_m(u, t % r1)
      call read_m(u, t % rc)
      call read_m(u, t % n)
      call read_m(u, t % m)
      
   end subroutine read_tail

   subroutine read_pft(p,u,t,n)
      type(pft_t), intent(inout) :: p
      type(ifl_t), intent(inout) :: u
      integer, intent(in), optional :: t, n
      
      if (present(t)) then
         p % t = t
      else
         call read_m(u, p%t)
      end if
      
      if (present(n)) then
         p % n = n
      else
         call read_m(u, p%n)
      end if
      
!       select case (p%t)
!          case (1)
!             p%f => gsph
!          case (2)
!             p%f => spln
!          case (3)
!             p%f => rexs
!          case (4)
!             p%f => plnm
!       end select
      
      allocate(p%a(p%n))
      call read_m(u, p%a)
      
   end subroutine read_pft
         
   subroutine read_hop(hop, u)
      type(hop_t), intent(inout) :: hop
      type(ifl_t), intent(inout) :: u
      character(len=3) :: l
      
      
      call read_m(u, l)
      call read_m(u, hop % v)
      call read_m(u, hop % dtau)
      call read_m(u, hop % Oscr)
!         call read_m(u, hop % fit)
      call read_pft(hop % sc, u, 1, 4)
      
      call read_tail(hop % tl, u)
   end subroutine read_hop

   subroutine read_pwp(pwp, u)
      type(pwp_t), intent(inout) :: pwp
      type(ifl_t), intent(inout) :: u
      character(len=3) :: l
      
      call read_m(u, l)
      call read_pft (pwp%fp, u)
      call read_m(u, l)
      call read_tail(pwp%tl, u)
      
   end subroutine read_pwp
   
   subroutine read_env(e, u)
      type(env_t), intent(inout) :: e
      type(ifl_t), intent(inout) :: u
      character(len=3) :: l
   
      call read_m(u, l)
      call read_m(u, e%a)
      call read_m(u, e%c)
      call read_m(u, e%nu)
      call read_m(u, e%rnu)
      
      call read_m(u, l)
      call read_tail(e%tl1, u)
      call read_m(u, l)
      call read_tail(e%tl2, u)
      
   end subroutine read_env
   
   subroutine read_bond_conf(bond_conf, u)
      type(bond_conf_t), intent(inout) :: bond_conf
      type(ifl_t), intent(inout) :: u
      integer :: i
      logical :: ovl
      character(len=10) :: l(10)
      
      call read_m(u, bond_conf%name)
      
      allocate(bond_conf % h(0:9))

      do i = 0, 9
         call read_hop(bond_conf % h(i), u)
      end do
      
      bond_conf % sss => bond_conf % h(0)
      bond_conf % sps => bond_conf % h(1)
      bond_conf % pps => bond_conf % h(2)
      bond_conf % ppp => bond_conf % h(3)
      bond_conf % sds => bond_conf % h(4)
      bond_conf % pds => bond_conf % h(5)
      bond_conf % pdp => bond_conf % h(6)
      bond_conf % dds => bond_conf % h(7)
      bond_conf % ddp => bond_conf % h(8)
      bond_conf % ddd => bond_conf % h(9)
      
      if (goto_next_uncmt(u, 'ovl')) then
         call read_m(u,ovl)
         if (ovl) then
            allocate(bond_conf % o(0:9))
            do i = 0, 9
              call read_hop(bond_conf % o(i), u)
            end do
          endif    
      endif
      
      call read_env(bond_conf % env, u)
      call read_pwp(bond_conf % pwp, u)
   end subroutine read_bond_conf

   subroutine read_ham_conf(ham_conf, u)
      use mod_const, only : hnam
      
      type(ham_conf_t), intent(inout), target :: ham_conf
      type(ifl_t), intent(inout) :: u
      integer :: i, j, k, l
      integer, parameter :: symm_orb_offsets(0:2) = [0,2,7]
      
      type(bond_conf_t), pointer :: b(:) => null()
      
      call read_m(u, ham_conf%sctb)
      
      call read_m(u, ham_conf%na)

      call read_m(u, ham_conf%nb)
      if (.not. allocated(ham_conf % a)) then
         allocate( ham_conf % a(0: ham_conf%na-1))
      else
         if (size(ham_conf % a) /=  ham_conf%na) then
               deallocate (ham_conf % a )
               allocate( ham_conf % a (0: ham_conf%na-1))
         end if
      end if
      do i=0,ham_conf%na-1
         call read_atom_conf(ham_conf%a(i), u,ham_conf%sctb)
         ham_conf%a(i) % id = i
      end do
      

      if (.not. allocated(ham_conf % b)) then
         allocate( ham_conf % b(0: ham_conf%nb-1))
      else
         if (size(ham_conf % b) /= ham_conf%nb) then
               deallocate (ham_conf % b )
               allocate( ham_conf % b (0: ham_conf%nb-1))
         end if
      end if
      
      b => ham_conf % b
      
      
      do i=0,ham_conf%nb-1
          
         call read_bond_conf(b(i), u)
         
         if (associated(ham_conf%b(i)%o)) then
    !        write(6, '("OVERLAP = TRUE")')
            ham_conf%ovl = .true.
         endif
         b(i) % id = i 
         do j=0,ham_conf%na-1 
               if (trim(ham_conf % a(j) % symb) == trim(b(i) % name(1:2))) then                 
                  b(i) % atom0 => ham_conf % a(j)
               end if              
               if (trim(ham_conf % a(j) % symb) == trim(b(i) % name(3:4))) then
                  b(i) % atom1 => ham_conf % a(j)
               end if
         end do

      end do
      
      

           
      do i=0,ham_conf%nb-1
!          print *, b(i) % atom0 % symb
!          print *, b(i) % atom1 % symb
         if (associated(b(i) % rev)) cycle
!          print *, b(i) % atom0 % symb
!          print *, b(i) % atom1 % symb
         if (b(i) % atom0 % z == b(i) % atom1 % z) then
               b(i) % rev => b(i)
         else
               do j = i+1, ham_conf % nb - 1
!                 print *, j, size(ham_conf % b)
!                   print *, b(j) % atom0 % z
!                   print *, b(j) % atom1 % z
                  
                  
                  if (associated(b(j) % rev)) cycle
                  if         (b(j) % atom0 % z == b(i) % atom1 % z &
                     & .and. b(j) % atom1 % z == b(i) % atom0 % z) then
                     b(i) % rev => b(j)
                     b(j) % rev => b(i)
                     exit
                  end if
               end do
         end if
      end do

      do i=0,ham_conf%nb-1
         if (b(i) % atom0 % z < b(i) % atom1 % z) then     
               do j = 0,2
                  if (b(i) % atom0 % b % lorbs(j) .and. b(i) % atom1 % b % lorbs(j)) then
                     do k = 0,j
                           l = symm_orb_offsets(j)+k
                           write(fout,'(a)') b(i) % rev % name//'.'//hnam(l) &
                                       & //' => '//b(i) % name//'.'//hnam(l)
                           if (b(i)%rev%h(l)%v /= b(i)%h(l)%v) call hop_mismatch_warning(b(i),hnam(l))
                           b(i)%rev%h(l) =  b(i)%h(l)
                     end do
                  end if
               end do
         end if
      end do
      
      
   contains
      subroutine hop_mismatch_warning(b,hop_name)
!             use mod_io
         implicit none
         
         type(bond_conf_t), intent(in) :: b
         character(len=*), intent(in) :: hop_name
         
         
         write(ferr,"(a)") 'WARNING!! '//b % rev % name//'.'//hop_name//' != ' &
                           & //b % name//'.'//hop_name//' using '//b % name//"'s"
      end subroutine hop_mismatch_warning
      
   end subroutine read_ham_conf

   subroutine read_cell(cell, u,emb)
      type(cell_t), intent(inout) :: cell
      type(ifl_t), intent(inout) :: u
      
      integer :: i
      logical :: l,emb
      
      call read_m(u, cell%name)
      call read_m(u, cell%str)
      
      call read_m(u, cell%a)
      cell % a = transpose(cell % a)
      
      l = goto_next_uncmt(u, 'len')
      call read_m(u, cell%lena)
      l = goto_next_uncmt(u, 'latpar')
      call read_m(u, cell%latpar)
      
      l = goto_next_uncmt(u, 'nd')
      call read_m(u, cell%nd)

      if (.not. goto_next_uncmt(u, 'nbp') .and. emb) stop 'embedding requires nbp (N. bop atoms) to be set after nd'
      if (goto_next_uncmt(u, 'nbp'))  call read_m(u, cell%nbp)

      
      if (.not. allocated(cell % d)) then
         allocate( cell % d(cell % nd))
      else
         if (size(cell % d) /= cell % nd) then
               deallocate (cell % d )
               allocate( cell % d(cell % nd))
         end if
      end if
      
      do i=1,cell%nd
         call read_atom(cell%d(i), u)
      end do
      
      if (goto_next_uncmt(u, 'ninert')) then
         call read_m(u, cell%ninert)
         if (cell%ninert > 0) then
            if (.not. allocated(cell % dinert)) then
               allocate( cell % dinert(cell % ninert))
            else
               if (size(cell % dinert) /= cell % ninert) then
                     deallocate (cell % dinert )
                     allocate( cell % dinert(cell % ninert))
               end if
            end if
            
            do i=1,cell%ninert
               call read_atom(cell%dinert(i), u)
            end do
         end if
      end if
      
      
      if (goto_next_uncmt(u, 'unrld')) then
         if (.not. allocated(cell % unrld)) then
            allocate( cell % unrld(cell % nd))
         else
            if (size(cell % unrld) /= cell % nd) then
                  deallocate (cell % unrld )
                  allocate( cell % unrld(cell % nd))
            end if
         end if
         
         do i=1,cell%nd
            call read_atom(cell%unrld(i), u)
         end do
      end if
      
      cell%nullheight = 0.0_dp
      if (goto_next_uncmt(u, 'nullheight')) call read_m(u, cell%nullheight)
      
      
      if (goto_next_uncmt(u, 'kspace_conf')) then 
         if (.not. associated(cell % kspace_conf)) allocate(cell % kspace_conf)
         call read_kspace_conf_t(cell % kspace_conf, u)
      end if
      
      if (goto_next_uncmt(u, 'dislocation_conf')) then
         if (.not. associated(cell % dislocation_conf)) allocate(cell % dislocation_conf)
         call read_dislocation_conf(cell % dislocation_conf, u)
      end if
      
   end subroutine read_cell


   subroutine read_kspace_conf_t(kspace_conf, u)

      type(kspace_conf_t), intent(inout) :: kspace_conf
      type(ifl_t), intent(inout) :: u
      integer :: i

      call read_m(u, kspace_conf % kmode)
      call read_m(u, kspace_conf % symon)
      
      select case(kspace_conf % kmode)
      case(0)

            call read_m(u, kspace_conf % vnkpts)            
            call read_m(u, kspace_conf % offcentre)                
      case(1)

            call read_m(u, kspace_conf % kbase)
            call read_m(u, kspace_conf % nk)
            if (allocated(kspace_conf % kpts)) deallocate(kspace_conf % kpts, kspace_conf % wtk)
            allocate(kspace_conf % kpts(3,kspace_conf % nk),kspace_conf % wtk(kspace_conf % nk))
            do i=1,kspace_conf % nk
               call read_m(u, kspace_conf % kpts(:,i))
               call read_m(u, kspace_conf % wtk(i))
            end do
      case default
            write(6,'("Illegal value passed to KMODE: ")') kspace_conf % kmode
            stop
      end select

   end subroutine read_kspace_conf_t



   subroutine read_elastic_const_conf(elastic_const_conf, u)
      type(elastic_const_t), intent(inout) :: elastic_const_conf
      type(ifl_t), intent(inout) :: u
      
      integer :: i
      
      do i=1,3
         call read_m(u, elastic_const_conf%tmatrix(i,:))
      end do
      call read_m(u, elastic_const_conf % polyord)
      call read_m(u, elastic_const_conf % tlo)
      call read_m(u, elastic_const_conf % thi)
      call read_m(u, elastic_const_conf % nt)
      
      
   end subroutine read_elastic_const_conf

   subroutine read_bop_conf(bop_conf,mix, u)
      type(bop_conf_t), intent(inout) :: bop_conf
      type(ifl_t), intent(inout) :: u
      logical :: mix
      
      if (.not. goto_next_uncmt(u, 'bop_conf')) stop 'bop_conf section not found! exitting ..'
      call read_m(u, bop_conf%momflg)
      call read_m(u, bop_conf%nrec)
      call read_m(u, bop_conf%nbase)
      call read_m(u, bop_conf%term)
      call read_m(u, bop_conf%v1flag)
      call read_m(u, bop_conf%chi_meth)
      call read_m(u, bop_conf%mfac)
      call read_m(u, bop_conf%kt)
      call read_m(u, bop_conf%etail)
      
      if (.not. goto_next_uncmt(u, 'sc_conf')) stop 'sc_conf section not found! exitting ..'
      call read_m(u, bop_conf%qerr)
      call read_m(u, bop_conf%mxit)
      call read_m(u, bop_conf%eftol)
      call read_m(u, bop_conf%mag)
      call read_m(u, bop_conf%merr)
      call read_m(u, bop_conf%mgit)
      
      if (mix) then
          call read_m(u, bop_conf%nit)
          call read_m(u, bop_conf%b)
          call read_m(u, bop_conf%wc)
      endif 
   
   end subroutine read_bop_conf


   subroutine read_run_conf(run_conf, u)
      type(run_conf_t), intent(inout) :: run_conf
      type(ifl_t), intent(inout) :: u
      logical :: avail
      integer :: i
      call read_m(u, run_conf%genfile)
      
      call read_m(u, run_conf%pot_flg)
      call read_m(u, run_conf%mvflag)
      call read_m(u, run_conf%quiet)
      call read_m(u, run_conf%paral)
      call read_m(u, run_conf%scf)
      call read_m(u, run_conf%writescale)
      call read_m(u, run_conf%elplusrel)
      call read_m(u, run_conf%rlim)
      call read_m(u, run_conf%fs)
      call read_m(u, run_conf%env)
      call read_m(u, run_conf%vpair)
      call read_m(u, run_conf%fs_only)
      call read_m(u, run_conf%emb)
      call read_m(u, run_conf%pot_scl)
      call read_m(u, run_conf%forces)
      call read_m(u, run_conf%scf_cut)
      call read_m(u, run_conf%rprune)
      call read_m(u, run_conf%rcut)
      call read_m(u, run_conf%idebug)
      call read_m(u, run_conf%restart)
      
      if (.not. associated(run_conf%bop_conf)) allocate(run_conf%bop_conf)
      
      if (run_conf%pot_flg == 9 .or. run_conf%pot_flg == 10) then
          call read_bop_conf(run_conf%bop_conf,.true., u)
      else
          call read_bop_conf(run_conf%bop_conf,.false., u)
      endif
      
      if (goto_next_uncmt(u, 'ham_conf')) then
 !         write(6, '("BOP ham_conf found...")')
          call read_ham_conf(run_conf%ham_conf, u)
      endif
      
      if (goto_next_uncmt(u, 'tbham_conf')) then
 !         write(6, '("TB ham_conf found...")')
          call read_ham_conf(run_conf%tbham_conf, u)
      endif      
      
      if (goto_next_uncmt(u, 'elastic_const_conf')) then
         if (.not. associated(run_conf%elastic_const_conf)) allocate(run_conf%elastic_const_conf)
         call read_elastic_const_conf(run_conf%elastic_const_conf, u)
      end if

      if (goto_next_uncmt(u, 'relaxation_conf')) then
         if (.not. associated(run_conf % relaxation_conf)) allocate(run_conf % relaxation_conf)
         call read_relaxation_conf(run_conf % relaxation_conf, u)
      end if
      
      if (goto_next_uncmt(u, 'md_conf')) then
         if (.not. associated(run_conf%md_conf)) allocate(run_conf%md_conf)
         call read_md_t(run_conf%md_conf, u)
      end if
      
      if (goto_next_uncmt(u, 'neb_conf')) then
         if (.not. associated(run_conf%neb_conf)) allocate(run_conf%neb_conf)
         call read_neb_conf(run_conf%neb_conf, u)
      end if
     
      
   end subroutine read_run_conf
   
   
   
   
   subroutine read_dislocation_conf(dsl, u)
      type(dislocation_t), intent(inout) :: dsl
      type(ifl_t), intent(inout) :: u
      
      call read_m(u, dsl % gfbcon) 
      call read_m(u, dsl % radiusi)
      call read_m(u, dsl % plotlims)
   end subroutine read_dislocation_conf
   
      
   end module mod_conf
    
    
    
!     
!     program test_mod_conf
!     
!     
!         use mod_conf
!         
!         type(cell_t), target :: cell, cellin
!         type(bop_conf_t), target :: bop
!         type(elastic_const_t), target :: elcon
!         type(ham_conf_t), target :: ham
!         type(run_conf_t) :: rt, rtin
!         
!         cell % name = 'cle'
!         cell % str = 'hcp'
!         
!         cell % a(1,:) = (/ 1,2,3 /)
!         cell % a(2,:) = (/ 4,5,6 /)
!         cell % a(3,:) = (/ 7,8,9 /)
!         
!         cell % lena = (/ 2.3, 4.5, 6.7/)
!         cell % latpar = 1.0_dp
!         
!         cell % nd = 2
!         cell % ninert = 0
!         allocate(cell % d(cell % nd), cell % dinert(cell % ninert))
!         
!         cell % d(1) % crds = (/ 1.0_dp , 3.4_dp, 2.1_dp/)
!         cell % d(2) % crds = (/ .0_dp , .4_dp, .1_dp/)
!         cell % d(1) % symb = 'Al'
!         cell % d(2) % symb = 'Ti'
!         
!         cell % nullheight = 0.0
!         
!         
!         ham % na = 2
!         ham % nb = 3
!         allocate( ham % a( 0:ham % na-1 ), ham % b( 0: ham % nb-1 ) )
!         
!         ham % a(0) % symb = 'Al'
!         ham % a(0) % ravs = .1
!         ham % a(0) % col = (/ 0.1 , 0.3 ,0.2/)
!         ham % a(0) % b % norbs = 1
!         ham % a(0) % b % orblist = (/1,0,0/)
!         ham % a(0) % b % cc = 2.5
!         ham % a(0) % b % es = 0.
!         ham % a(0) % b % ep = 0.3421
!         ham % a(0) % b % ed = 0.
!         ham % a(0) % b % xtra_coefs = (/ 0.5, 0., 0., 0./)
!         
!         ham % a(0) % e % lam_0 = 1.4
!         ham % a(0) % e % m = 2.0
!         ham % a(0) % e % r_core = 0.84
!     
!         ham % a(1) % symb = 'Ti'
!         ham % a(1) % ravs = .2
!         ham % a(1) % col = (/ 0.5 , 0.1 ,0.7/)
!         ham % a(1) % b % norbs = 1
!         ham % a(1) % b % orblist = (/2,0,0/)
!         ham % a(1) % b % cc = 2.1
!         ham % a(1) % b % es = 0.
!         ham % a(1) % b % ep = 0.
!         ham % a(1) % b % ed = 0.6725
!         ham % a(1) % b % xtra_coefs = (/ 0.5, 0., 0., 0./)
!        
!         ham % a(1) % e % lam_0 = 1.8
!         ham % a(1) % e % m = 2.0
!         ham % a(1) % e % r_core = 0.845
!         
!         ham % b(0) % name = 'AlAl'
!         ham % b(0) % pps % v = 1.3
!         ham % b(0) % ppp % v = -1.0
!         
!         ham % b(0) % a = 1231.
!         ham % b(0) % c = 11.
!         ham % b(0) % nu = 3.2
!         
!         ham % b(1) % name = 'AlTi'
!         ham % b(1) % pds % v = 0.9
!         ham % b(1) % pdp % v = 0.6
!         
!         ham % b(1) % a = 310.
!         ham % b(1) % c = 14.
!         ham % b(1) % nu = 3.3
!         
!         ham % b(2) % name = 'TiTi'
!         ham % b(2) % dds % v = 1.1
!         ham % b(2) % ddp % v = -0.2
!         ham % b(2) % ddd % v = .01
!         
!         ham % b(2) % a = 561.
!         ham % b(2) % c = 16.
!         ham % b(2) % nu = 3.1
!         
!         ham % b(0) % atom0 => ham % a(0)
!         ham % b(0) % atom1 => ham % a(0)
!         ham % b(1) % atom0 => ham % a(0)
!         ham % b(1) % atom1 => ham % a(1)
!         ham % b(2) % atom0 => ham % a(1)
!         ham % b(2) % atom1 => ham % a(1)
!         
!         
!         
!     
! 
!         
!         elcon % nt = 5
!         elcon % polyord = 3
!         
!         elcon % tmatrix(1,:) = (/1,0,0/)
!         elcon % tmatrix(2,:) = (/0,1,0/)
!         elcon % tmatrix(3,:) = (/0,0,1/)
!         
!         elcon % tlo = -0.005
!         elcon % thi = 0.005
!         
!         
!         bop % momflg = 1
!         bop % nrec = 4
!         bop % nbase = 4
!         bop % term = 1
!         bop % v1flag = 0
!         bop % chi_meth = 2
!         bop % mfac = 100
!         bop % kt = 0.3_dp
!         bop % etail = 1.0_dp
!         bop % qerr = 0.01_dp
!     
!     
!     
!         
!         rt % genfile = 'tial'
!         rt % cell => cell
!         rt % pot_flg = 2
!         rt % mvflag = -1
!         rt % scf = 0
!         rt % writescale = 0
!         rt % elplusrel = 0
!         rt % rlim = (/1,1,1/)
!         rt % fs = 1
!         rt % env = 1
!         rt % vpair = 0
!         rt % fs_only = 0
!         rt % emb = 0
!         rt % pot_scl = 0
!         rt % idebug = 0
!         rt % restart = 0
!         rt % scf_cut = 0.
!         rt % rcut = (/4.80_dp, 5.50_dp/)
!         rt % rprune = 20.
!         rt % bop_conf => bop
!         rt % ham_conf = ham
!         rt % elastic_const_conf => elcon
!         
!         
!         
!         
!         
! !         call print_run_conf( rt, '', 9)
! !         call print_cell( rt % cell ,'', 10)
!         
!         
!     
!         deallocate( ham % a, ham % b)
!         deallocate(cell % d ,cell % dinert)
!         
! !         flush(9)
! !         flush(10)
! !         close(9)
! !         close(10)
!         
!         call read_cell( cellin, 10)
!         call read_run_conf(rtin, 9)
!         rtin % cell => cellin
!         
!         call print_run_conf( rtin, '', 99)
!         call print_cell( rtin % cell ,'', 100)
!         
!         flush(99)
!         flush(100)
!         close(99)
!         close(100)
!         
!         
!         
!         
!         
!     end program test_mod_conf


