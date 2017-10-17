
   ! 
   !     
   !     subroutine loadhop( hop, k, i  )
   !         use mod_precision
   !         use mod_all_scalar
   !         use mod_const
   !         use mod_gsp
   !         use mod_conf
   !         
   !         
   !         implicit none
   !         
   !         type(hop_t), intent(in) :: hop
   !         integer, intent(in) :: k, i
   ! 
   ! 
   !         include "Include/Atom.array"
   ! 
   ! !         real(dp) :: phi0,oscr,vscr,deltau
   ! !         real(dp) :: r0,rc,r1,rcutin,n,nc
   ! !         real(dp) :: d0,dc,d1,dcutin,m,mc
   ! !         real(dp) :: delta
   ! !         real(dp) :: cv0,cv1,cv2,cv3,cv4,cv5
   !         real(dp) :: y0
   !         
   ! 
   ! 
   ! !         vscr   = hop % v 
   ! !         deltau = hop % dtau
   ! !         oscr   = hop % oscr
   ! !         r0     = hop % r0
   ! !         rc     = hop % rc 
   ! !         r1     = hop % r1
   ! !         rcutin = hop % rcut
   ! !         n      = hop % n 
   ! !         nc     = hop % nc 
   ! 
   ! ! v * (r/r0)**(-n) * exp(-n*((r/rc)**nc - (r0/rc)**nc ))
   ! ! v/(r0**(-n)*exp(-n*(r0/rc)**nc)) * (r**(-n)*exp(-n*(r/rc)**nc))
   ! 
   ! !         if (hop % n==0.0_dp) then
   ! !             if (hop % nc == 1.0_dp) then
   ! !                bnddat(i,k)    = hop % v  * exp(hop % r0/hop % rc)
   ! !             else if (nc == 0.0_dp) then
   ! !                bnddat(i,k)    = hop % v  * exp(1.0_dp)
   ! !             else
   ! !                bnddat(i,k)    = hop % v  * exp((hop % r0/hop % rc)**hop % nc)
   ! !             end if
   ! !         else if (hop % n == 1.0_dp) then
   ! !             if (hop % nc == 1.0_dp) then
   ! !                bnddat(i,k)    = hop % v * hop % r0 * exp(hop % r0/hop % rc)
   ! !             else if (nc == 0.0_dp) then
   ! !                bnddat(i,k)    = hop % v * hop % r0 * exp(1.0_dp)
   ! !             else
   ! !                bnddat(i,k)    = hop % v * hop % r0 * exp((hop % r0/hop % rc)**hop % nc)
   ! !             end if
   ! !         else
   ! !             if (hop % nc == 1.0_dp) then
   ! !                bnddat(i,k)    = hop % v * hop % r0**hop % n * exp(hop % n * hop % r0/hop % rc)
   ! !             else if (hop % nc == 0.0_dp) then
   ! !                bnddat(i,k)    = hop % v * hop % r0**hop % n * exp(hop % n)
   ! !             else
   ! !                bnddat(i,k)    = hop % v * hop % r0**hop % n * exp(hop % n*(hop % r0/hop % rc)**hop % nc)
   ! !             end if               
   ! !         end if
   ! 
   ! !         print *, hop % fit, k, i
   ! 
   !         bnddat(i,k)    = hop % v 
   ! 
   !         if (hop%r0 > hop%r1) then
   !             y0 = (hop%r0-hop%r1)/(hop%rcut-hop%r1)
   !             bnddat(i,k) = bnddat(i,k)/(((6.0_dp*y0 + 3.0_dp)*y0 + 1.0_dp) * (1.0_dp - y0)**3)
   ! !             print *, 'y0', bnddat(i,k), hop % v *  (((6.0_dp*y0 + 3.0_dp)*y0 + 1.0_dp) * (1.0_dp - y0)**3)
   !             
   !         end if
   ! 
   !         bndscl(13,i,k) = hop % oscr
   !         bndscl(14,i,k) = hop % dtau
   !         bndscl(1,i,k)  = hop % r0
   !         bndscl(2,i,k)  = hop % rc
   !         bndscl(3,i,k)  = hop % r1
   !         bndscl(4,i,k)  = hop % rcut
   !         bndscl(5,i,k)  = hop % n
   !         bndscl(6,i,k)  = hop % nc
   ! 
   ! ! These are the polynomial coefs.
   ! ! Set to zero because they are largely useless with the multiplicative cutoff 
   ! ! and only cause formatting problems in outfill.f90. The labels there are wrong anyway.
   ! 
   !                bndscl(7,i,k) = 0.0_dp   
   !                bndscl(8,i,k) = 0.0_dp
   !                bndscl(9,i,k) = 0.0_dp
   !                bndscl(10,i,k) = 0.0_dp
   !                bndscl(11,i,k) = 0.0_dp
   !                bndscl(12,i,k) = 0.0_dp
   ! 
   ! ! Commented out see the note above
   ! !         ! Calculation of 5th order spline polynomial
   ! !         delta = rcutin - r1
   ! !         if (delta == 0.0) then
   ! !             cv0=0.0
   ! !             cv1=0.0
   ! !             cv2=0.0
   ! !             cv3=0.0
   ! !             cv4=0.0
   ! !             cv5=0.0
   ! !             
   ! !         else
   ! !             cv0 = gsp(r1,r0,rc,n,nc)
   ! !             cv1 = dgsp(r1,r0,rc,n,nc)
   ! !             cv2 = d2gsp(r1,r0,rc,n,nc)/2
   ! ! 
   ! !             cv3 = -10.0_dp*(cv0/(delta**3))-6.0_dp*(cv1/(delta**2)) - 3.0_dp*(cv2/delta)
   ! ! 
   ! !             cv4 = 15.0_dp*(cv0/(delta**4))+8.0_dp*(cv1/(delta**3)) + 3.0_dp*(cv2/(delta**2))
   ! ! 
   ! !             cv5 = -6.0_dp*(cv0/(delta**5))-3.0_dp*(cv1/(delta**4)) - (cv2/delta**3)
   ! ! 
   ! !         !     write(6,*) bndscl(1,1,1)
   ! !         !     write(6,*) bndscl(1,i,k)
   ! ! 
   ! !         end if
   ! !         
   ! !         bndscl(7,i,k) = cv0   
   ! !         bndscl(8,i,k) = cv1
   ! !         bndscl(9,i,k) = cv2
   ! !         bndscl(10,i,k) = cv3
   ! !         bndscl(11,i,k) = cv4
   ! !         bndscl(12,i,k) = cv5
   ! 
   !         !        write(*,*) r0,rc,r1,rcutin,n,nc
   ! 
   !         !        write(*,*) cv0
   !         !        write(*,*) cv1
   !         !        write(*,*) cv2
   !         !        write(*,*) cv3
   !         !        write(*,*) cv4
   !         !        write(*,*) cv5
   !         !
   ! 
   ! 
   ! 
   ! 
   !     end subroutine loadhop

      



      
   subroutine reloadtb( ham_conf )
         
      use mod_precision
      use mod_all_scalar
      use mod_const
      use mod_gsp
      use mod_conf
      use mod_tail
      use mod_pft
      use topologia, only : iproc, master

      implicit none

      include "Include/Atom.array"

      type(ham_conf_t), intent(in), target :: ham_conf
      type(bond_conf_t), pointer :: b
      type(hop_t), pointer :: h,o
      integer :: zin, i, j, k, rk, li, lj, oi, oj, bi, u
      real(dp) :: fv(0:2), r
      logical :: ovl
      
      ovl = ham_conf % ovl
!       do j = 0,natype
        do j = 0,rtc % ham_conf % na-1
         zin = ham_conf % a(j) % z
         
         zc(zin) = ham_conf % a(j) % b % cc 
         nl(zin) = ham_conf % a(j) % b % norbs
         llist(:,zin) = ham_conf % a(j) % b % orblist
         nstt(zin) = 0
         
         if (mag) istn(zin) = ham_conf % a(j) % b % st
         
         do i = 1,nl(zin)
            nstt(zin) = nstt(zin) + 2*llist(i,zin)+1
         enddo
         
         es(zin) = ham_conf % a(j) % b % es
         ep(zin) = ham_conf % a(j) % b % ep
         ed(zin) = ham_conf % a(j) % b % ed
   !          embed(:,zin) = ham_conf % a(j) % b % useless_coefs(:)
         
      end do
      
      
   !       
      do i = 0, ham_conf % nb - 1
         b => ham_conf % b(i)
         do oj = 1, b % atom1 % b % norbs
            lj = b % atom1 % b % orblist(oj)
            do oi = 1, b % atom0 % b % norbs
               li = b % atom0 % b % orblist(oi)
               if (li > lj) cycle
               do bi = 0, li 
   !                     rescale here to speedup gsp later
                  h => b % h(hsidx(li,lj,bi))
                  call init_pfts(h % sc, h % tl)
                  
                  
                  if (ovl) then
                     o => b % o(hsidx(li,lj,bi))
                     call init_pfts(o % sc, o % tl)
                  endif
                  
               end do
            end do
         end do 
         if (rtc%vpair /= 0) call init_pfts( b % pwp % fp, b % pwp % tl)
      end do

      
      nullify(b,h)
      if (ovl) nullify(o)
   !       
   !       
   !       bndscl = 0.0_dp
   !       
   !       do i = 0, ham_conf % nb - 1
   !             b => ham_conf % b(i)
   !             k  = btype(b % atom0 % z, b % atom1 % z)
   !             rk = btype(b % atom1 % z, b % atom0 % z)
   !             
   !             do oj = 1, b % atom1 % b % norbs
   !                 lj = b % atom1 % b % orblist(oj)
   !                 do oi = 1, b % atom0 % b % norbs
   !                     li = b % atom0 % b % orblist(oi)
   !                     if (li > lj) cycle
   !                     do bi = 0, li  !relies on correct repeating info in b%h for ss, pp, dd, ...
   !                         call loadhop(b % h(hsidx(li,lj,bi)),  k, hlidx(li,lj,bi)) 
   !                         if (li /= lj) call loadhop(b % h(hsidx(li,lj,bi)),  rk, hlidx(lj,li,bi)) 
   !                     end do
   !                 end do
   !             end do 
   !       end do
   !       
   ! 
   !       if ((.not. quiet) .and. (iproc == master) ) then
   !         do i=1, ham_conf % nb
   !             k = btype(atype(mod(i-1,2)),atype((i-1)/2))
   !             print *, atype(mod(i-1,2)),atype((i-1)/2),k
   !             do j=1,14
   !                 print *, hlnam(j), bnddat(j,k)
   !             end do
   !         end do
   !       end if

      end subroutine reloadtb
      
      

   subroutine reinit_block( cell )

      use mod_precision
      use mod_all_scalar
      use mod_const
      use mod_conf
      use topologia, only : mpmap, iproc, aproc, sdistrib, nproc_bop => nproc ! this and the following call sdistrib(...) are of dubious quality
      use mod_atom_ar
      
      implicit none

      interface
         function getzquick(symb)
            implicit none
            character(len=2), intent(in) :: symb
            integer :: getzquick
         end function getzquick
      end interface 


      type(cell_t), intent(in) :: cell

      include "Include/Atom.array"
      include "Include/Misc.array"
      include "Include/PosVel.array"
      include "Include/ag.conn"


      !
      !    Declare the simple variables.
      !


      integer if_relaxed
      integer i,j, num
      integer stat,restart
      integer reclen,tokbgn,tokend,error
      
      real(dp) :: getmass

      character(len=2) :: chemsymb


      !
      !    Read in the cell parameters.
      !

      latpar = cell % latpar

      a = cell % a

      lena = cell % lena

      nd = cell % nd
      if (nd > mxnd) then
         write(6,'(''Too many atoms. Increase MXND.'')')
         stop
      endif

      call sdistrib(nd,nproc_bop,mpmap,aproc)
      
      if (.not.quiet) then
         print *, 'aproc:', aproc
         print *, 'mpmap:', mpmap
      end if
      
      
   !         d = cell % d % crds


   !         There should be separation between atom parallelism and kpoints parallelism in the TB branch
   !        may be introduce kpmap and replace mpmap with apmap.Later the atom parallelism for BOP should
   !        be made linear by carefully distributing to the cluster only without repeating the overlapping atoms.
   !        This seems quite tricky and I haven't thought it out yet. Back to kspace: dq, mass and mg will be reallocated in kspace_setup

      if (.not. allocated(dq) ) then
            allocate(dq(mpmap(iproc)+1:mpmap(iproc+1)))
      else
            if ((mpmap(iproc+1) - mpmap(iproc)) /= size(dq)) then
               deallocate(dq)
               allocate(dq(mpmap(iproc)+1:mpmap(iproc+1)))
            end if
      end if
      
   !       if (.not. allocated(mass) ) then
   !             allocate(mass(mpmap(iproc)+1:mpmap(iproc+1)))
   !       else
   !             if ((mpmap(iproc+1) - mpmap(iproc)) /= size(mass)) then
   !                deallocate(mass)
   !                allocate(mass(mpmap(iproc)+1:mpmap(iproc+1)))
   !             end if
   !       end if
      
      if (.not. allocated(mass) ) allocate(mass(nd))
      
      do i = 1,nd
         symb(i) = cell % d(i) % symb
         d(1:3,i) = cell % d(i) % crds(1:3)
   !             z(i) = getz(cell % d(i) % symb)
         z(i) = getzquick(cell % d(i) % symb)
      end do    
      
      
      do i = 1, nd            
         de(i) = cell % d(i) % de            
      enddo
      
      
      if (mag) then
         
   !          if (.not. allocated(mg) ) then
   !                allocate(mg(mpmap(iproc)+1:mpmap(iproc+1)))
   !          else
   !                if ((mpmap(iproc+1) - mpmap(iproc)) /= size(mg)) then
   !                   deallocate(mg)
   !                   allocate(mg(mpmap(iproc)+1:mpmap(iproc+1)))
   !                end if
   !          end if
   !          do i = mpmap(iproc)+1, mpmap(iproc+1)            
   !             mg(i) = cell % d(i) % mg            
   !          enddo
         
   !          the above shall be looked at later
         if (.not. allocated(mg) ) allocate(mg(nd))
         do i = 1, nd            
            mg(i) = cell % d(i) % mg            
         enddo
         
      end if
      
      ninert = cell % ninert
      
      if (ninert > minert) then
         write(6,'(''Too many inert atoms. Increase MINERT.'')')
         stop
      endif

      if (ninert > 0) then

         do i = 1,ninert
               symi(i) = cell % dinert(i) % symb
               dinert(:,i) = cell % dinert(i) % crds(:)
   !                 zinert(i) = getz(cell % dinert(i) % symb)
               zinert(i) = getzquick(cell % dinert(i) % symb)
         enddo
         
         if (mag) then
            if (.not. allocated(demi) ) then
                  allocate(demi(ninert))
            else
                  if (ninert /= size(demi)) then
                     deallocate(demi)
                     allocate(demi(ninert))
                  end if
            end if
            
!             This is not the full demi, it is updated after istn is loaded by reloadtb.
            do i = 1, ninert
               demi(i) = 0.5_dp * cell % dinert(i) % mg
            end do
         end if
                  
      endif


      if   ( mvflag == 14  .or.  mvflag == 15 ) zeroheight = cell % nullheight


      if ( mvflag == 16 .or. mvflag == 17 .or. mvflag == 18 &
                        & .or. mvflag == 19 .or. mvflag == 30 ) then
         
         do i = 1,nd
            unrld(:,i) = cell % unrld(i) % crds(:)
         enddo

         if_relaxed = cell % if_relaxed
         
         do i = 1,nd
            mass(i) = getmass(z(i))
         enddo

      endif
      !
      !     Modified by M. Cawkwell. Dislocation line now paralled to z.
      !     29th July 2003
      !
      if ( mvflag == 18 )  then
         x_center = cell % x_center
         y_center = cell % y_center
         r_min    = cell % r_min
         r_max    = cell % r_max
         n_mesh   = cell % n_mesh
      endif

      if ( mvflag == 18 .or. mvflag == 19 ) e_coh = cell % e_coh
      if ( mvflag == 19 ) log_rad = cell % log_rad

      !
      !    Initialise some variables.
      !

      iter = 1
      eflo = 0.d0
      
      if (mvflag == 18) eflo = cell % eflo
      

      ! zmena - puvodne coment na prvni podmince

      if ( ((mvflag == 16) .and. (if_relaxed == 1)) ) then

         !      IF ( ((MVFLAG .EQ. 16) .AND. (IF_RELAXED .EQ. 1)) .OR.
         !     &      (MVFLAG .EQ. 17) .OR.  (MVFLAG .EQ. 18)   ) THEN
         iter = cell % iter
         eflo = cell % eflo
      endif

   end subroutine reinit_block

      
      
      
      subroutine reset_rlx(rlxc)
         use mod_conf, only : relaxation_t
         use mod_all_scalar
         implicit none
         
         type(relaxation_t), intent(in) :: rlxc
      
         include 'Include/Relax.array'
      
         integer :: i
        

       
         autosave = rlxc % autosave
         mxiter   = rlxc % mxiter  
         rlxflg   = rlxc % rlxflg  
         step     = rlxc % step    
         ftol     = rlxc % ftol    
         cnst_n   = rlxc % cnst_n  

         if (cnst_n > mcnst_n) then
            write(6,'(''Too many constraints.'')')
            write(6,'(''CNST_N = '',I5,'' MCNST_N = '',I5)') cnst_n,mcnst_n
            stop
         endif

         if (cnst_n > 0) then
            
            do i = 1,cnst_n
                  cnst_a(i) = rlxc % cnst_a(i)
                  cnst_v(:,i) = rlxc % cnst_v(i,:)
                  cnst_v(:,i) = cnst_v(:,i)/sqrt(sum(cnst_v(1:3,i)*cnst_v(1:3,i)))
            enddo
         endif

         
         writeper = rlxc % writeper
         datform = rlxc % datform
      
      end subroutine reset_rlx

      
      subroutine resetup(run_conf)
         
         use mod_precision
         use mod_all_scalar
         use mod_const
         use mod_srt
         use pair_coefficients
   !         use env_coefficients
         use ab_io
         use mod_gsp
         use mod_conf
         use mod_ham
         use mod_funptr
         use mod_kspace, only : kspace_setup
         use mod_atom_ar
         
         implicit none

         
         interface 
            subroutine reloadtb(ham_conf)
                  import ham_conf_t
                  type(ham_conf_t), intent(in), target :: ham_conf
            end subroutine reloadtb
         end interface
         
         
         type(run_conf_t), intent(in),target :: run_conf


   !         include "Include/Parallel.const"

         include "Include/Atom.array"
         include "Include/BEC.array"

   !         include "Include/Hamilt.array"
         include "Include/Misc.array"
         include "Include/NebList.array"
         include "Include/PosVel.array"
         include "Include/Relax.array"
         include "Include/Vib.array"
         
   !         include "Include/Parallel.array"
         include "Include/ag.conn"
         include 'Include/Envir.array'



         !  Declare functions.

         !       real(dp) ::  gsp, dgsp, d2gsp, 
         real(dp) :: vee, dvee, venv, chi, vcore
         type(ham_conf_t),pointer :: hamc


         !
         !    Declare the simple variables.
         !

         real(dp) :: rcheck, vscale, dscale, scalecore, scaledbond
         real(dp) :: scalebond, scalerep, scalechi, scaleenv
         real(dp) :: fcc_chi, nnb1, nnb2, nnb3, mpw
         real(dp) :: x,instem,trs
         real(dp) :: getmass,geteatom
         real(dp) :: xstart,ystart,reltype
         real(dp) :: diff, dum1, dum2, dum3
         real(dp) :: bs(300,15)
         real(dp) :: testf

         integer :: writescale
         integer :: ptr,bt
         integer :: ng1,ng2
         integer :: i,j,k
         integer :: stat,restart,nla
         integer :: reclen,tokbgn,tokend,error
         integer :: mstat,nstata,cnstat
         integer :: mfcomp

         character(len=80) :: filename
         character(len=1) :: answer

         logical :: ex

         real(dp) :: v3(3)

         !
         !    Declare other common blocks.
         !

         common /gamcom/ xstart,ystart,ng1,ng2
         !
         !      EXTERNAL SCALE
         !
         !    Open the input files.
         !
         quiet = run_conf % quiet
         genfile = run_conf % genfile
         lengfn = index(genfile,' ') - 1
         !

         !  Read in whether we want to print out the scaling dependence.

         et = 0.0_dp
         
         writescale = run_conf % writescale

         !  Read in the structural type.

         
         strtyp = trim(run_conf % cell % str)

         rlim = run_conf % rlim

         rcut = run_conf % rcut(1)
         rcutl = run_conf % rcut(2)

         ! new k-space flags 

         potflg = run_conf % pot_flg


         forces = run_conf % forces

         fs_include = run_conf % fs
         vpair_include = run_conf % vpair
         env_include = run_conf % env
                        
         scf_include = run_conf % scf
         if (scf_include == 0 ) then
            if (.not. quiet) write(*,*)  "NO SCREENING !!!"
         endif

         scf_cut = run_conf % scf_cut

         fs_only = run_conf % fs_only

         emb_include = run_conf % emb

         if ( fs_include == 0 )  then
            fs_only = 0
            emb_include = 0
         endif 

         pot_scl = run_conf % pot_scl

         momflg = run_conf % bop_conf % momflg
         if (momflg == 2) then
            write(6,'(/''Locally oriented orbitals not implemented.'')')
            write(6,'(''==> Setting MOMFLG = 3.'')')
            momflg = 3
         endif
         
         call select_feprom(momflg)
         
         nrec = run_conf % bop_conf % nrec
         
         if (nrec > mrec) then
            write(6,'(''Too many recursion levels asked for.'')')
            write(6,'(''Setting NREC = '',I2)') mrec
            nrec = mrec
         endif

         nbase = run_conf % bop_conf % nbase
         
         if (nbase > nrec) then
            write(6,'(''NBASE > NREC. This is wasteful of memory.'')')
            write(6,'(''=> Setting NBASE = NREC.'')')
            nbase = nrec
         endif

         qerr = run_conf % bop_conf % qerr

         eftol = run_conf % bop_conf % eftol
         
         mag = run_conf % bop_conf % mag
         if (mag) nsp = 2
         
         
         
         term = run_conf % bop_conf % term

         kt  = run_conf % bop_conf % kt
         
         if (kt < 1.0e-3_dp) then
            write(6,'(''KT is too small for numerical stability.'')')
            write(6,'(''Setting KT = 0.001eV.'')')
            kt = 1.0e-3_dp
         endif

         restart = run_conf % restart

         rprune = run_conf % rprune

         v1flag = run_conf % bop_conf % v1flag


         chi_meth = run_conf % bop_conf % chi_meth
         if ((chi_meth <= 0).or.(chi_meth >= 3)) then
            write(6,'(''Invalid CHI_METH flag found.'')')
            write(6,'(''Setting CHI_METH=1 (analytic integration).'')')
            chi_meth = 1
         endif
         if ((chi_meth == 2).and.(term /= 1)) then
            write(6,'(''Numerical integration for response functions  '')')
            write(6,'(''only valid for square root termintor (TERM=1).'')')
            write(6,'(''Setting CHI_METH=1 (analytic integration).'')')         
            chi_meth = 1
         endif
         if (chi_meth == 2) then
            if (.not. quiet) print *,'**** NUMERICAL INTEGRATION FOR CHIs ***'
         endif

         call select_sri(chi_meth)
         
         mfac = run_conf % bop_conf % mfac 

         idebug = run_conf  % idebug


         !
         !    Read in atom mover parameters.
         !

         mvflag = run_conf % mvflag


         !  Read in the block.

   !         call inblock( )
         call reinit_block( run_conf % cell )
         
         


         !  Read the rest of the parameters from the input file (MVFLAG specific).

         do_subblock = 0
         if ( mvflag == 14  .or.  mvflag == 15 )  then
            do_subblock = 1
            call sort_out()
         endif


         if (mvflag == 0) then

            etail = run_conf % bop_conf % etail

   ! Relaxation
   !
   !     Calculating lattice greens functions - direction of test force
   !
         elseif (mvflag == 1 .or. mvflag == 24 .or. mvflag == 25) then

         !         CALL FNDREC(8,'MFCOMP',STAT)
         !         READ(8,*) MFCOMP

         !         CALL FNDREC(8,'TESTF',STAT)
         !         READ(8,*) TESTF

         !
            call reset_rlx(run_conf % relaxation_conf)
            
         elseif ((mvflag == 2).or.(mvflag == 3)) then

            trace = run_conf % md_conf % trace

            autosave = run_conf % relaxation_conf % autosave

            mxiter   = run_conf % relaxation_conf % mxiter

            
            dt = run_conf % md_conf % dt
            
            temp = run_conf % md_conf % temp

            if (mvflag == 3) then
                  ttol = run_conf % md_conf % ttol
            endif

            nequil = run_conf % md_conf % nequil

            datform  = run_conf % relaxation_conf % datform
            autosave = run_conf % relaxation_conf % autosave

            monitorper = run_conf % md_conf % monitorper

            do i = 1,nd
               mass(i) = getmass(z(i))
            enddo

            mda_gr = run_conf % md_conf % mda_gr

            if (mda_gr > 0) then
                  mda_gr_r = run_conf % md_conf % mda_gr_r
                  mda_gr_n = run_conf % md_conf % mda_gr_n
                  if ((mda_gr_n < 1).or.(mda_gr_n > mbins)) mda_gr_n = mbins
                  mda_gr_dr = mda_gr_r/real(mda_gr_n, dp)
                  do i = 1,mda_gr_n,1
                     gr_bins(i) = 0
                  enddo
            endif

            mda_g3 = run_conf % md_conf % mda_g3

            if (mda_g3 > 0) then
                  mda_g3_r = run_conf % md_conf % mda_g3_r
                  mda_g3_n = run_conf % md_conf % mda_g3_n
                  if ((mda_g3_r < 0.0_dp).or.(mda_g3_r > rcut)) mda_g3_r = rcut
                  if ((mda_g3_n < 1).or.(mda_g3_n > mbins)) mda_g3_n = mbins
                  do i = 1,mda_g3_n,1
                     g3_bins(i) = 0
                  enddo
            endif

            mda_sw = run_conf % md_conf % mda_sw

            if (mda_sw > 0) then
                  mda_sw_e = run_conf % md_conf % mda_sw_e
                  mda_sw_n = run_conf % md_conf % mda_sw_n
                  if ((mda_sw_n < 1).or.(mda_sw_n > mbins)) mda_sw_n = mbins
                  mda_sw_de = mda_sw_e/real(mda_sw_n, dp)
                  do i = 1,mda_sw_n,1
                     sw_bins(i) = 0
                  enddo
            endif

            mda_adav = run_conf % md_conf % mda_adav

            mda_op = run_conf % md_conf % mda_op

   !!!!!!!!!NO MORE PATIENCE !!!!!!!!   SHALT BE FINISHED WHEN NEEDED         
            
            read(8,*) latparam,(kcpt(i),i=1,3,1)

            call fndrec(8,'MDA_MSD',stat)
            read(8,*) mda_msd
            if (mda_msd == 1) then
                  read(8,*) msd_atom
            endif

         elseif (mvflag == 4) then

            trace = run_conf % md_conf % trace

            autosave = run_conf % relaxation_conf % autosave

            mxiter   = run_conf % relaxation_conf % mxiter

            dt = run_conf % md_conf % dt

            temp = run_conf % md_conf % temp

            ttol = run_conf % md_conf % ttol

            call fndrec(8,'PRESS',stat)
            read(8,*) press

            call fndrec(8,'MPISTON',stat)
            read(8,*) mpiston
            mpiston = mpiston*1.0e8_dp

            nequil= run_conf % md_conf % nequil

            datform  = run_conf % relaxation_conf % datform
            autosave = run_conf % relaxation_conf % autosave

            call fndrec(8,'MONITORPER',stat)
            read(8,*) monitorper

            do i = 1,nd,1
               mass(i) = getmass(z(i))
            enddo

            call fndrec(8,'MDA_GR',stat)
            read(8,*) mda_gr

            if (mda_gr > 0) then
            read(8,*) mda_gr_r,mda_gr_n
            if ((mda_gr_n < 1).or.(mda_gr_n > mbins)) & 
            &         mda_gr_n = mbins
            mda_gr_dr = mda_gr_r/real(mda_gr_n, dp)
            do i = 1,mda_gr_n,1
            gr_bins(i) = 0
            enddo
            endif

            call fndrec(8,'MDA_G3',stat)
            read(8,*) mda_g3

            if (mda_g3 > 0) then
            read(8,*) mda_g3_r,mda_g3_n
            if ((mda_g3_r < 0.0_dp).or.(mda_g3_r > rcut)) & 
            &         mda_g3_r = rcut
            if ((mda_g3_n < 1).or.(mda_g3_n > mbins)) & 
            &         mda_g3_n = mbins
            do i = 1,mda_g3_n,1
            g3_bins(i) = 0
            enddo
            endif

            call fndrec(8,'MDA_SW',stat)
            read(8,*) mda_sw

            if (mda_sw > 0) then
            read(8,*) mda_sw_e,mda_sw_n
            if ((mda_sw_n < 1).or.(mda_sw_n > mbins)) & 
            &         mda_sw_n = mbins
            mda_sw_de = mda_sw_e/real(mda_sw_n, dp)
            do i = 1,mda_sw_n,1
            sw_bins(i) = 0
            enddo
            endif

         elseif (mvflag == 5) then

            trace = run_conf % md_conf % trace

            autosave = run_conf % relaxation_conf % autosave

            mxiter   = run_conf % relaxation_conf % mxiter

            dt = run_conf % md_conf % dt

            nequil = run_conf % md_conf % nequil

            datform  = run_conf % relaxation_conf % datform
            autosave = run_conf % relaxation_conf % autosave

            call fndrec(8,'MONITORPER',stat)
            read(8,*) monitorper

            call fndrec(8,'TMIN',stat)
            read(8,*) tmin

            call fndrec(8,'TMAX',stat)
            read(8,*) tmax

            temp = tmax
            if (mxiter <= nequil) then
            dtemp = 1.0_dp
            else
            dtemp = (tmin/tmax)**(1.0_dp/real(mxiter-nequil-1, dp))
            endif

            do i = 1,nd,1
               mass(i) = getmass(z(i))
            enddo

         elseif ((mvflag == 6).or.(mvflag == 12).or.(mvflag == 22) & 
         &        .or. (mvflag == 20) .or. (mvflag == 26) .or. (mvflag == 36)) then

            !  Elastic constants run.

            tmatrix = run_conf % elastic_const_conf % tmatrix

            tlo = run_conf % elastic_const_conf % tlo
            thi = run_conf % elastic_const_conf % thi
            nt = run_conf % elastic_const_conf % nt

            polyord = run_conf % elastic_const_conf % polyord

            !
            !     Added by M. Cawkwell 2nd November 2004 
            !     Needed if you want to relax internal degrees of freedom
            !     when calculating the elastic constants - important for C11b!!!
            !
            elplusrel = run_conf % elplusrel
            if ( elplusrel == 1 ) call reset_rlx(run_conf % relaxation_conf)

         elseif (mvflag == 7) then

         !  Diffusion barrier height run.

            call reset_rlx(run_conf % relaxation_conf)


            call fndrec(8,'D_ATOM',stat)
            read(8,*) d_atom

            call fndrec(8,'RC_DIFF',stat)
            read(8,*) rc_diff

            call fndrec(8,'DISP_VEC',stat)
            read(8,*) disp_vec

            call fndrec(8,'N_DISP',stat)
            read(8,*) n_disp

         elseif (mvflag == 9) then

         !     calculate vibrational frequencies

            call fndrec(8,'NPER',stat)
            read(8,*) nper

            call fndrec(8,'DELTA',stat)
            read(8,*) delta

         elseif (mvflag == 10) then

         !     calculate gamma surface.

            call reset_rlx(run_conf % relaxation_conf)

            call fndrec(8,'GAMMA',stat)
            read(8,*)ng1,ng2
            read(8,*)xstart,ystart

            !
            !     Cartesian offset for the gamma-S
            !

            xstart = xstart*lena(1)
            ystart = ystart*lena(2)



         elseif (mvflag == 13) then

            !  Calculate the energy surface.  (A.G.)

            call fndrec(8,'ESURFRANGE',stat)
            read(8,*) zmina, zmaxa, numpointsa
            read(8,*) zminc, zmaxc, numpointsc
            elplusrel = run_conf % elplusrel
            
            if (elplusrel == 1) call reset_rlx(run_conf % relaxation_conf)

         elseif (mvflag == 14) then   

            !  Calculate the gamma-surface.   (A.G.)

            call fndrec(8,'INDREL',stat)
            read(8,*) indrelax

            call fndrec(8,'GAMSURFRANGE',stat)
            read(8,*) gamxmin, gamxmax, gamxnp
            read(8,*) gamymin, gamymax, gamynp

            call fndrec(8,'SYM',stat) 
            read(8,*) req 

            ok = run_conf % relaxation_conf % ok

            call fndrec(8,'RELTYPE',stat)
            read(8,*) reltype
            if      ( reltype == 1 )  then
            const_volume = 0
            else if ( reltype == 2 )  then
            const_volume = 1 
            endif

            call reset_rlx(run_conf % relaxation_conf)
         
         elseif (mvflag == 15) then

            !  Calculate the grain boundary structure.  (A.G.)

            call fndrec(8,'XSHIFTGB',stat)
            read(8,*)      xshiftgb
            write( 6,'("xshift",F8.5//)' ) xshiftgb

            call fndrec(8,'YSHIFTGB',stat)  
            read(8,*)      yshiftgb           
            write( 6,'("yshift",F8.5//)' ) yshiftgb

            call fndrec(8,'ZSHIFTGB',stat)    
            read(8,*)      zshiftgb
            write( 6,'("zshift",F8.5//)' ) zshiftgb

            call fndrec(8,'GBMOVEOPT',stat)  
            read(8,*)      gbmoveopt          

            ok = run_conf % relaxation_conf % ok

            call fndrec(8,'RELTYPE',stat)
            read(8,*) reltype
            if      ( reltype == 1 )  then
               const_volume = 0
            else if ( reltype == 2 )  then  
               const_volume = 1
            endif


            call reset_rlx(run_conf % relaxation_conf)

         elseif (mvflag == 16 .or. mvflag == 17 .or. mvflag == 18 ) then

         !     Relaxing dislocation with(out) GFBCs

            gfbcon = run_conf % cell % dislocation_conf % gfbcon
            
            if (gfbcon  /=  0) then
               print *, "Relaxing dislocation with GFBCs"
               radiusi = run_conf % cell % dislocation_conf % radiusi
               print *, "Radius of region I =", radiusi
            endif

         !  Calculate the dislocation structure.  (A.G.)

            ok = run_conf % relaxation_conf % ok

            plxmin = run_conf % cell % dislocation_conf % plotlims(1)
            plxmax = run_conf % cell % dislocation_conf % plotlims(2)
            plymin = run_conf % cell % dislocation_conf % plotlims(3)
            plymax = run_conf % cell % dislocation_conf % plotlims(4)


            call reset_rlx(run_conf % relaxation_conf)
            
         elseif (mvflag == 23) then

         !  Volume relaxation run.

            autosave = run_conf % relaxation_conf % autosave
            mxiter   = run_conf % relaxation_conf % mxiter
            rlxflg   = run_conf % relaxation_conf % rlxflg
            step     = run_conf % relaxation_conf % step    
            ftol     = run_conf % relaxation_conf % ftol

            tmatrix = run_conf % elastic_const_conf % tmatrix

            tlo = run_conf % elastic_const_conf % tlo
            thi = run_conf % elastic_const_conf % thi
            nt = run_conf % elastic_const_conf % nt

            polyord = run_conf % elastic_const_conf % polyord


         elseif (mvflag == 27) then

            !  Calculate the energy/volume surf  (D.P.)

            call fndrec(8,'VOLRANGE',stat)
            read(8,*) zminv, zmaxv, numpointsv

            call fndrec(8,'RATIORANGE',stat)
            read(8,*) zminrt, zmaxrt, numpointsrt

            elplusrel = run_conf % elplusrel
            if (elplusrel == 1) call reset_rlx(run_conf % relaxation_conf)
            
         else if (mvflag == 30) then
         
            call reset_rlx(run_conf % relaxation_conf)
!             do not call init of neb here because ad etc are not set up yet

         endif
         
         !
         !    Read from restart file if requested.
         !

         if (restart == 1) then
            inquire(file = 'BOP.dump', exist = ex)
            if (ex) then
                  write(6,'(/''BOP.dump found. Do you wish to use it ? '')')
                  read(5,'(A)') answer
                  if ((answer == 'y').or.(answer == 'Y')) then
                     open(unit = 1,file = 'BOP.dump',status = 'OLD')
                  else
                     restart = 0
                  endif
            else
                  restart = 0
            endif
         elseif (restart == 2) then
            inquire(file = 'BOP.dump', exist = ex)
            if (ex) then
                  open(unit = 1,file = 'BOP.dump',status = 'OLD')
                  write(6,'(/''BOP.dump found =>'', & 
                  &                 '' Reading in restart information.''/)')
            else
                  restart = 0
            endif
         elseif (restart == 3) then
   35          write(6,'('' Enter the name of the restart file > '',$)')
            read(5,'(A)') filename
            open(unit = 1,file = filename,status = 'OLD',err = 37)
            goto 36
   37          write(6,'('' Unable to open the restart file : '',A)') filename
            write(6,'('' Please try again.'')')
            goto 35
         else
            restart = 0
         endif

   36    if (restart /= 0) then

            call fndrec(8,'A',stat)
            read(8,*) a(1,1),a(1,2),a(1,3), & 
            &             a(2,1),a(2,2),a(2,3), & 
            &             a(3,1),a(3,2),a(3,3)

            call fndrec(8,'RV',stat)
            read(8,*) (ad(1,i),ad(2,i),ad(3,i), & 
            &              vel(1,i),vel(2,i),vel(3,i),i=1,nd,1)

            call fndrec(8,'DE',stat)
            read(8,*) (de(i),i=1,nd,1)

            call fndrec(8,'Z',stat)
            read(8,*) (z(i),i=1,nd,1)

            if (mvflag > 0) then

               do i = 1,nd,1
                  mass(i) = getmass(z(i))
               enddo

               if (cnst_n > 0) then
                  call fndrec(8,'CNST_A',stat)
                  read(8,'(10I5)') (cnst_a(i),i=1,cnst_n,1)
               endif

               call fndrec(8,'ITER',stat)
               read(8,*) iter
               iter = iter + 1

               call fndrec(8,'SUME',stat)
               read(8,*) sume

               call fndrec(8,'SUMEE',stat)
               read(8,*) sumee

               call fndrec(8,'SUMT',stat)
               read(8,*) sumt

               call fndrec(8,'SUMTT',stat)
               read(8,*) sumtt

               call fndrec(8,'SUMP',stat)
               read(8,*) sump

               call fndrec(8,'SUMPP',stat)
               read(8,*) sumpp

               call fndrec(8,'SUMV',stat)
               read(8,*) sumv

               call fndrec(8,'SUMVV',stat)
               read(8,*) sumvv

               if ((mda_gr > 0).and.(iter > nequil)) then
                  call fndrec(8,'GR',stat)
                  read(8,*) (gr_bins(i),i=1,mda_gr_n,1)
               endif

               if ((mda_g3 > 0).and.(iter > nequil)) then
                  call fndrec(8,'G3',stat)
                  read(8,*) (g3_bins(i),i=1,mda_g3_n,1)
               endif

               if ((mda_sw > 0).and.(iter > nequil)) then
                  call fndrec(8,'SW',stat)
                  read(8,*) (sw_bins(i),i=1,mda_sw_n,1)
               endif

               if ((mda_adav > 0).or.(mda_msd > 0)) then
                  call fndrec(8,'NOPBSUM',stat)
                  if (stat == 0) then
                     read(8,*) (adnopb(1,i),adnopb(2,i),adnopb(3,i), & 
                     &                 adsum(1,i),adsum(2,i),adsum(3,i),i=1,nd,1)
                  else
                     print *,'INITIALISING NOPB COORDS'
                     do i = 1, nd
                        do j = 1, 3
                           adnopb(j,i) = ad(j,i)
                           adsum(j,i) = 0.0_dp
                        enddo
                     enddo                  
                  endif                  
               endif

               if (mda_msd == 1) then
                  if (stat == 0) then
                     call fndrec(8,'ADMSDA',stat)
                     read(8,*) (admsda(i),i=1,3,1)
                  else
                     print *,'INITIALISING MSD ATOM'
                     do j = 1, 3, 1
                        admsda(j) = ad(j,msd_atom)
                     enddo                  
                  endif
               endif

               if (mda_msd == 2) then
                  call fndrec(8,'DISPO',stat)
                  if (stat == 0) then
                     read(8,*) (dispo(i),i=1,nd,1)
                  else
                     print *,'INITIALISING DISPO ARRAY'
                     do i = 1, nd, 1
                           dispo(i)=sqrt(ad(1,i)**2+ad(2,i)**2+ad(3,i)**2)
                     enddo
                  endif
               endif

            endif

            close(1)

            if ((mvflag >= 2).and.(mvflag <= 4)) then
                  trs = instem(vel,nd,mass)
                  if (trs == 0.0_dp) then
                     call initvel(vel,nd,temp,mass)
                     trs = instem(vel,nd,mass)
                  endif
                  x = sqrt(temp/trs)
                  call rescale(vel,nd,x)
            endif

         else

            !
            !       Initialise some variables.
            !

            sume = 0.0_dp
            sumee = 0.0_dp
            sumt = 0.0_dp
            sumtt = 0.0_dp
            sump = 0.0_dp
            sumpp = 0.0_dp
            sumv = 0.0_dp
            sumvv = 0.0_dp


            !
            !       Scale the primitive translation vectors.
            !
            

            do i = 1,3
               a(i,1:3) = a(i,1:3)*lena(i)/sqrt(sum(a(i,1:3)*a(i,1:3)))
            enddo

            !
            !       Set up the atomic positions in cartesian coordinates.
            !

            do i = 1,nd
                
               ad(1,i) = a(1,1)*d(1,i)+a(2,1)*d(2,i)+a(3,1)*d(3,i)
               ad(2,i) = a(1,2)*d(1,i)+a(2,2)*d(2,i)+a(3,2)*d(3,i)
               ad(3,i) = a(1,3)*d(1,i)+a(2,3)*d(2,i)+a(3,3)*d(3,i)
               
            enddo
            if (mvflag == 16 .or. mvflag == 17 .or. mvflag == 18 &
                           &  .or. mvflag == 19 .or. mvflag == 30) then
               do i = 1, nd
                  call mul3x3(a,unrld(1,i),v3)
                  unrld(:,i) = v3
               end do   
!                unrld(1:3,:nd) = matmul(transpose(a), unrld(1:3,:nd))
            end if   
!                
               
            if ((mda_adav == 1).or.(mda_msd > 0)) then
               do i = 1, nd, 1
                  do j = 1, 3, 1
                     adnopb(j,i) = ad(j,i)
                     adsum(j,i) = 0.0_dp
                  enddo
               enddo
            endif

            if (mda_msd == 1) then
               do j = 1, 3, 1
                  admsda(j) = ad(j,msd_atom)
               enddo
            endif

            if (mda_msd == 2) then
               do i = 1, nd, 1
                  dispo(i) = sqrt(ad(1,i)**2 + ad(2,i)**2 + ad(3,i)**2)
               enddo
            endif

         !
         !       Initialise velocities, if needed.
         !

            if ((mvflag >= 2).and.(mvflag <= 5)) then
               call initvel(vel,nd,temp,mass)
               x = sqrt(temp/instem(vel,nd,mass))
               call rescale(vel,nd,x)
            endif

         endif

         !
         !    Rescale RLIM according to RPRUNE.
         !

         do i = 1,3,1
            if (rlim(i) /= 0) then
               rlim(i) = int(rprune/lena(i))+1
            endif
         enddo

         !
         !    Evaluate the volume of the cell.
         !

         vol = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) & 
         &   + a(1,2)*(a(2,3)*a(3,1)-a(2,1)*a(3,3)) & 
         &   + a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))

         !
         !    Set up the atomic positions in cartesian coordinates for the inert atoms.
         !

         do i = 1,ninert
            adinert(1,i) = a(1,1)*dinert(1,i)+a(2,1)*dinert(2,i)+a(3,1)*dinert(3,i)
            adinert(2,i) = a(1,2)*dinert(1,i)+a(2,2)*dinert(2,i)+a(3,2)*dinert(3,i)
            adinert(3,i) = a(1,3)*dinert(1,i)+a(2,3)*dinert(2,i)+a(3,3)*dinert(3,i)
         enddo

         !
         !    Find out which atom types have been used.
         !

   ! Note        
   ! zatype(z) -> type
   ! atype(type) -> z

         atype = 0

         do i = 1,nd
            atype(z(i)) = 1
         enddo

         do i = 1,ninert
            atype(zinert(i)) = 1
         enddo

         
         natype = -1
         do i = 0,mxz
            if (atype(i) == 1) then
               natype = natype + 1
               atype(natype) = i
            endif
         enddo

         do i = natype+1,matype-1
            atype(i) = 0
         enddo
         
         do i = 1,natype
         zatype(atype(i)) = i
         end do

         !    Build list of bond type pointers.
         
         btype = 0
         nbtype = 0
         do i = 0,natype
            do j = 0,natype
                  nbtype = nbtype + 1
                  if (nbtype > mbtype) then
                     write(6,'(''Too many types of bond. Increase MBTYPE.'')')
                     call panic()
                  endif
   !                 print *, i,j, atype(j),atype(i), nbtype
                  btype(atype(j),atype(i)) = nbtype
            enddo
         enddo

! This will make BOP use the tbham to produce bond types etc if tbham exists, so it must have the same
! atom and bond ordering as bop ham, although bop ham need not have all the atoms or bonds
         if (potflg == 9 .or. potflg == 10) then
            hamc => run_conf % tbham_conf
         else
            hamc => run_conf % ham_conf
         endif
   ! asort(z_ord_type) -> conf_atom_idx
   ! bsort(z_ord_btype) -> conf_bond_idx
         asort = -1
         do i = 0, hamc % na - 1
            do j = 0,natype
                  if (atype(j) == hamc % a(i) % z) then
                     asort(j) = i
                     exit
                  endif
            end do
         end do
         
         bsort = -1
         
         do i = 0, hamc % nb - 1
            bsort(btype(hamc%b(i)%atom0%z, hamc%b(i)%atom1%z)) = i
         enddo
         
         
         nullify(hamc)
   !         print *,'asort',asort
   !         print *,'bsort',bsort
         
   !         stop
         !
         !    Load in Tight Binding data.
         !
! This is just needed for the BOP part to set its zc etc in the common block, so doesn't use hamc
         call reloadtb(run_conf % ham_conf )

!   This is not in reinit_block yet because of a circular dependence: the demi values depend on reloadtb to load istn. 
!   reloadtb on the other hand depends on zinert through atype but zinert is loaded in reinit_block
         if (mag) then
            do i = 1, ninert
               demi(i) = istn(zinert(i)) * demi(i)
            end do
         end if   
         
         !
         !    Make sure not too many states appear on one atom.
         !

         mstat = 0

         cnstat = 0 

         do i = 1,nd
            call states(z(i),nla,nstata,llista)
            cnstat = cnstat + nstata
            mstat = max(mstat,nstata)
         enddo

         do i = 1,ninert
            call states(zinert(i),nla,nstata,llista)
            mstat = max(mstat,nstata)
         enddo

        if (potflg == 9 .or. potflg == 10) then
            do i = 0,natype
                nstata = 0
                do j = 1, run_conf % tbham_conf % a(i) % b % norbs
                    nstata = nstata +  2*(run_conf % tbham_conf % a(i) % b % orblist(j))+1
                enddo
                mstat = max(mstat,nstata)
            enddo 
        endif

         if (mstat > mxnstat) then
            write(6,'(/''Too many orbitals on one atom.'')')
            write(6,'(''Increase MXNSTAT to '',I2)') mstat
            call panic()
         elseif (mstat < mxnstat) then
            write(6,'(/''MXNSTAT is bigger than is needed.'')')
            write(6,'(''You can decrease MXNSTAT to '',I2)') mstat
         endif

   !        kmxh = cnstat

         !
         !    Evaluate the total band energy for the isolated atoms.
         !


         if (potflg == 4) then 
               call kspace_setup()
         end if
         
   !          if (iproc >= aproc) call exit_orderly()
         
         
         eatom = geteatom(z,nd)

         !
         !    Set roots for finite temperature with square root terminator.
         !

         !     CALL ATNRTS()

         !
         !    Set up polynomial coefficients for Fermi-Dirac distribution.
         !

         call fdpoly()
         call pascal()


         !  Set the starting value of EFLO to zero.

         if(mvflag /= 18) then
            eflo = 0.0_dp
         endif



   !         !
   !         !  This piece writes out the scaling dependences of hopping integrals
   !         !  and the repulsive potentials to the channels 33 and 34.
   !         !

   !         if ( writescale  /=  0 )    then
   !             open( unit = 33, file='bond.dat')
   !             open( unit = 34, file='pair.sc')
   !             open( unit = 35, file='chi.sc')
   !             open( unit = 36, file='env.sc')
   !             open( unit = 37, file='core.sc')
   ! 
   !             do i = 1, 300, 1
   !                 k=1
   !                 rcheck = 2.0_dp + 0.01*real( i , dp)
   !                 bs(i,k) = rcheck
   !                 do j=1,14
   !                     bt = writescale
   !                     if (bndscl(1,j,bt) /= 0.0) then
   !                         k=k+1
   !                         scalebond = vscale(rcheck,bndscl(1,j,bt),bndscl(2,j,bt), & 
   !                         &                      bndscl(5,j,bt),bndscl(6,j,bt), & 
   !                         &                      bndscl(3,j,bt),bndscl(4,j,bt), & 
   !                         &                      bndscl(7,j,bt),bndscl(8,j,bt), & 
   !                         &                      bndscl(9,j,bt),bndscl(10,j,bt), & 
   !                         &                      bndscl(11,j,bt),bndscl(12,j,bt))
   ! 
   !                         bs(i,k) = scalebond*bnddat(j,bt)
   !                     endif
   !                 enddo
   ! 
   !                 write(33,'(14F12.6)') bs(i,1),bs(i,2),bs(i,3),bs(i,4),bs(i,5), &
   !                 &    bs(i,6),bs(i,7),bs(i,8),bs(i,9),bs(i,10),bs(i,11),bs(i,13), & 
   !                 &            bs(i,14)
   ! 
   ! 
   ! !                 scalerep = vee( rcheck, 0 )
   ! !                 scalecore = vcore(rcheck,1)
   !                 scalechi = chi(rcheck,0)
   !                 !
   !                 !     NNB's for fcc Ir
   !                 !
   !                 nnb1 = 2.714582933_dp
   !                 nnb2 = 3.839_dp
   !                 nnb3 = 4.701795561_dp
   !                 !
   !                 fcc_chi = 12.0_dp*chi(nnb1,0) + 6.0_dp*chi(nnb2,0) & 
   !                 &              + 24.0_dp*chi(nnb3,0)
   !                 fcc_chi = exp((1.0_dp/mpw(0))*log(fcc_chi))
   !                 scaleenv = venv(rcheck,0,fcc_chi)
   ! 
   !                 if ( rcheck  >=  rcut )   then
   !                     scalebond = 0._dp
   !                     !                  SCALEREP  = 0._dp
   !                     !                  SCALECHI = 0._dp
   !                     !                  SCALEENV = 0._dp
   !                     !                  SCALECORE = 0._dp
   !                 endif
   ! 
   ! !                 write ( 34, '(8X,F16.10,8X,F16.10)') rcheck, scalerep
   !                 write ( 35, '(8X,F16.10,8X,F16.10)') rcheck, scalechi
   !                 write ( 36, '(8X,F16.10,8X,F16.10)') rcheck, scaleenv
   ! !                 write ( 37, '(8X,F16.10,8X,F16.10)') rcheck, scalecore
   !             enddo
   !             close ( 33 )
   !             close ( 34 )
   !             close ( 35 )
   !             close ( 36 )
   !             close ( 37 )
   !         endif



      end subroutine resetup

   function getzquick(symb)
      use mod_conf, only : rtc
      implicit none
      character(len=2), intent(in) :: symb
      integer :: getzquick
      integer :: i

      
      if (rtc % pot_flg /=9 .and. rtc % pot_flg /= 10) then
            i = 0
            do while (i /= rtc % ham_conf % na .and. rtc % ham_conf % a (i) % symb /= symb)   
            i = i + 1
            end do
            getzquick = rtc % ham_conf % a (i) % z
      else
         i = 0 
         do while (i /= rtc % tbham_conf % na .and. rtc % tbham_conf % a (i) % symb /= symb)
            i = i + 1
         end do
         getzquick = rtc % tbham_conf % a (i) % z
         
      endif
      


   end function getzquick


