


   module tbbop_emb

   use mod_precision
   use mod_const
   use mod_conf
   use mod_clock
   use topologia

#ifdef MPI
   use mpi
#endif

   implicit none

   integer, parameter :: nr = 2  ! number of regions (only TB and BOP so far)

   
   type :: tbconf
       real(dp) :: istn(0:mxz),zc(0:mxz),es(0:mxz),ep(0:mxz),ed(0:mxz)
       integer :: nstt(0:mxz),nl(0:mxz),llist(ml,mxz)
   end type tbconf

   type mbtb_t
      logical :: ovl,sctb,mxgf,allbop        ! overlap, SC charge transfer and GF type switches
      integer :: rlst(nr+1) ! atom offsets for regions
      integer :: n(nr+1)    ! number of atoms according to region
      integer :: nd         ! total number of atoms
      integer :: thsz, bhsz, ahsz ! sizes for TB, BOP and global Hamiltonians respectively
!                         posham,   aptr_tb   bptr_tb   aptr_dh  bptr_dh 
      integer, allocatable :: pth(:), pbh(:), pah(:) ! atom offsets for Hamiltonians
      real(dp), allocatable :: ad(:,:)   ! coordinates of dynamic atoms
      integer, allocatable :: apt(:), bpt(:) ! neighbour lists for TB -> TB
      integer, allocatable :: apd(:), bpd(:) ! neighbour lists for TB -> BOP
      integer, allocatable :: apb(:), bpb(:) ! neighbour lists for BOP -> BOP

      real(dp), allocatable :: h0(:,:),h(:,:,:), s(:,:) ! TB input and updated hamiltonians, and overlap matrices in dense form
      real(dp), allocatable :: dh(:,:,:)      ! TB->BOP matrix elements in sparse form
      real(dp), allocatable :: q(:,:),qmpol(:,:),rho(:,:,:,:),dqchi(:) ! global atomic charges, tb multipole elements, density matrix in sparse form and Eb elements for force in N.O. basis
      real(dp), allocatable ::en(:,:) ! eigenvectors and eigenvalues for GF calc in Dyson embedding
      !, v(:), ! rho(mxnstat,mxnstat,jb,ia)

      complex(dp), allocatable :: a(:,:), ac(:,:), gt(:,:)!,hp(:,:)!In short:
!      a = G_BOP (-dh')
!      gt = (z - (h - dh a))^{-1}
!      ac = a gt
! see the notes in 'formuli.pdf' pages: 4-6 for details.

      complex(dp), allocatable :: cn(:,:,:), ggd(:) ! global Green matrix diagonal elements
      !, w(:,:), occ(:), gg(:,:)
      type(tbconf) :: tbc
      
      real(dp), allocatable :: vm(:,:),vu(:)
      integer,allocatable :: orth(:)
   end type
   
   type :: ctsx_t
      real(dp) :: gaunt(9,9,25)
      real(dp), allocatable :: strx(:,:,:,:), dstrx(:,:,:,:),U(:)
      integer, allocatable :: struxidx(:,:), dstrxidx(:,:)
      real(dp), allocatable :: struxd(:), dstrxd(:)
   end type ctsx_t
   
   
   contains

   subroutine emb_sc()
!    self consistency loop.
!    This is the entry point for the embedding program.
   
      use mod_all_scalar, only : mag, quiet, dqmax, forces, mgmax, nd, qerr, nsp ,locc
      use mod_atom_ar, only : dq, mg
      use mod_chi
      use topologia, only : iproc, master
      
      include "../Include/Atom.array"
      include "../Include/Force.array" 
      include "../Include/NebList.array"    
      
!       REMOVE THIS
      include "../Include/PosVel.array"


      type(mbtb_t) :: mbc ! embedding conf
      type(ctsx_t) :: ctsx ! charge transfer constants

      integer :: it, mxit, ia, mgit, mit
      real(dp) :: etott,bpons
      real(qp) :: avdde, sumz, merr
      real(qp), allocatable :: dec(:)
      real(dp), allocatable :: mg_o(:),qmpol_prev(:)

      character(len=25) :: ffmt

      integer,parameter :: mixing_type = 2

      integer(8) :: c1,c2

      integer :: i,j
      real(dp) :: ef,efb

      integer :: ib, ierr,gib,gia
      real(dp) :: dqmax_old, aa, aamin, aamax,det,deb,dea
      

      integer :: lb,mmix,git,nat
      logical :: calF
      
      integer :: nlmq1,nlmq
      
      real(dp), allocatable :: qmx(:,:,:)
      
      real(dp) :: beta,wc,rms,ecorr

      ef = 0.0_dp
      efb = 0.0_dp
      de =0.0_dp
      dem=0.0_dp
      bpons = 0.0_dp
      
      calF = .false.
      
      !       mxit = 500
      mxit = rtc % bop_conf % mxit 
      if (mag) then
        mgit = rtc % bop_conf % mgit 
        merr = rtc % bop_conf % merr
        mit = 1
        if (merr < 0.0_dp) mgit = 1
      end if

      
      if (.not. allocated(dq) ) then
         allocate(dq(nd))
      else
         if (nd /= size(dq)) then
            deallocate(dq)
            allocate(dq(nd))
         end if
      end if
      
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
       mbc % rlst(1) = 1
       mbc % rlst(2) = rtc % cell % nbp + 1
       mbc % rlst(3) = rtc % cell % nd + 1

       mbc % ovl = rtc % tbham_conf % ovl
       mbc % sctb = rtc % tbham_conf % sctb
       
!  This switch decides whether the GF for each part should be made separately and then correlating blocks 
!  introduced in the construction of the entire GF (potflg 9) or whether the system GF should be calculated 
!  with BOP and then the peturbation of switching the TB region from BOP to TB be calculated with Dyson's eqn
!  (potflg 10).
       if (rtc % pot_flg == 9) then
          if (.not. quiet) write(6,'("CARRYING OUT INVERSION EMBEDDING")')
          mbc % mxgf = .true.
       else
          if (.not. quiet) write(6,'("CARRYING OUT DYSON EMBEDDING")')
          mbc % mxgf = .false.
       endif
       mbc % allbop = .false.
       if (.not. quiet .and. mbc%ovl)  write(6,'(/,''OVERLAP = TRUE'')')      

! The cell file will specify the numbers nt, nb, and nd then list all atoms ordered as follows: 
! first the bop atoms then the tb atoms and afterwards the pure bop atoms. Pure bop atoms are not consdered yet 
! but should be kept in mind
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      call loadtb(rtc%tbham_conf,mbc%tbc,mbc%ovl)

      allocate( dec(nd))
      dec = real(de(:nd),qp)

      if (forces) then
          fbs(:,:nd) = 0.0_dp
      endif
      
      
      call mbldh(mbc,ctsx)
      
      
      if ((qerr < 0.0_dp).or.(mixing_type == 0)) mxit = 1
      if (.not. quiet ) then
         if (.not. mbc%sctb) write(9,"('Imposing LCN ...')")
         write(ffmt,"(a,i0,a)") '(a,":",', mbc%nd, '(x,f22.14))' 
      end if
      

      sumz = sum( (/ (real(zc(z(ia)), qp), ia = 1,mbc%rlst(2)-1) /) )
      

      if (mbc%sctb) then
          mmix = rtc%bop_conf%nit
          beta = rtc%bop_conf%b
          wc = rtc%bop_conf%wc
          nat = mbc%n(2) - sum(mbc%orth)
          if (.not. allocated(qmx)) allocate(qmx(nat,0:mmix+1,2))
          if (.not. allocated(qmpol_prev)) allocate(qmpol_prev(nat))
          qmx(:,:,:) = 0.0_dp
          nlmq = 1
          if (.not. forces) then
              nlmq1 = 1
          else
              nlmq1 = 4
          endif
          mbc % qmpol = 0.0_dp
          qmpol_prev = 0.0_dp
          mbc%dqchi = 0.0_dp
      endif
      
      
      
      
      do ia = 1,mbc%n(2)
          if (mbc%orth(ia) == 0) then
              sumz = sumz + mbc%tbc%zc(z(ia+mbc%n(1)))
              
          else
              sumz = sumz + zc(z(ia+mbc%n(1)))
          endif
      enddo
      
      locc = sumz
      
      
      allocate(mg_o(mbc%nd))
      if (mag) mg_o = mg      

      if (.not. mbc%mxgf) then
          call bop_gf(mbc,efb,mxit,mgit,merr)
          ef= efb
          mg_o(:mbc%n(1))=mg(:mbc%n(1))
          mg(mbc%n(1)+1:mbc%nd)=mg_o(mbc%n(1)+1:mbc%nd)
          
!           if (mbc%sctb) then
!               de(mbc%n(1)+1:mbc%nd) =0.0_dp
!               mbc%dqchi(mbc%n(1)+1:mbc%nd) = 0.0_dp
!           endif
          
      endif

      if (.not.(quiet .or. mbc%mxgf)) write(6,'(/,''Calculating TB portion of GF.'')')
      
      do 
         if (mag) then
            if (.not. quiet) write(6,'(/,"mit:",i3)') mit   
            

            mg_o = mg
            
            if (mbc%mxgf) then
                do ia = 1, mbc%n(1)
                    dem(ia) = 0.5_dp*istn(z(ia))*mg(ia)               
                enddo
            endif
            
            do ia = mbc%n(1)+1, mbc%nd
                if (mbc%orth(ia-mbc%n(1)) == 0) then
                  dem(ia) = 0.5_dp*mbc%tbc%istn(z(ia))*mg_o(ia)
                else
                  dem(ia) = 0.5_dp*istn(z(ia))*mg_o(ia) 
                endif 
            enddo
         end if

         dqmax_old = 0.0_dp
         aamax=1.0_dp
         aamin=1.0e-4_dp
         aa=10.0_dp
         

         it = 1
         do               
            if (.not. mbc % sctb) then
            
                deb =  sum( (/ (dec(ia)*zc(z(ia)), ia = 1,mbc%rlst(2)-1) /) )
                det =  sum( (/ (dec(ia)*mbc%tbc%zc(z(ia)), ia = mbc%rlst(2),mbc%rlst(3)-1) /) )
                do i = 1,nsp
                    
                    mbc%h(:,:,i) = mbc % h0(:,:)
                enddo
                
                if (it > 1) de(:mbc%nd) = real(dec - (det+deb)/sumz, dp)           
                do i = 1, nsp
                    do ib = 1, mbc % n(2)
                        gib = (mbc % rlst(2)-1) + ib
                        dea = de(gib)                      
                        if (mag) dea = dea + oppm(i)*dem(gib)
                        do lb = mbc % pth(ib)+1, mbc % pth(ib+1)
                            mbc % h(lb,lb,i) = mbc % h(lb,lb,i) + dea
                        end do
                    enddo
                end do
            else if ( it == 1) then
                do i = 1, nsp
                    mbc%h(:,:,i) = mbc % h0(:,:)
                    do ib = 1, mbc % n(2)
                        dea = 0.0_dp
                        gib = (mbc % rlst(2)-1) + ib
                        if (mbc%orth(ib) == 1) dea = dea + de(gib)   
!                         dea = de(gib)                      
                        if (mag) then
                            if (mbc%orth(ib) == 1) then
                                dea = dea + oppm(i)*dem(gib)
                            else
                                dea = dea + oppm(i)* 0.5_dp*mbc%tbc%istn(z(gib))*mg_o(gib)
                            endif
                        endif
                        do lb = mbc % pth(ib)+1, mbc % pth(ib+1)
                            mbc % h(lb,lb,i) = mbc % h(lb,lb,i) + dea      
                        end do
                    enddo
                end do
            endif
            
            if (mbc%mxgf .and. mbc%sctb) de(:mbc%n(1)) = de(:mbc%n(1))+mbc % dqchi(:mbc%n(1)) 
            
!             sumzde = 0.0_dp
!             sumzde = sum((/ (de(ia)*zc(z(ia)), ia = 1,mbc%n(1)) /))
!             
!             if (.not. quiet) print *, 'sumzde bop:',sumzde
!             do ia = 1,mbc%n(2)
!                 if (mbc%orth(ia) == 1) sumzde = sumzde + de(ia+mbc%n(1))*zc(z(ia+mbc%n(1)))
!             enddo
!             if (.not. quiet) print *, 'sumzde:',sumzde,sumzde/sumzo
!             do ia = 1,mbc%n(2)
!                 if (mbc%orth(ia) == 1) de(ia+mbc%n(1)) = de(ia+mbc%n(1)) - sumzde/sumzo
!             enddo
            
            if (.not. quiet) then
               write(6,'(/,"it:",i3)') it
!               write(6, ffmt) 'de ', de(:mbc%nd)
!               if (mag) write(6, ffmt) 'dem', dem(:mbc%nd)
            end if
            
            
            
 
 !            if (mbc % mpol) call tbmpol(mbc,ef)
           
            call mbeq(mbc,ef,sumz)
            
            
            if (.not. mbc % sctb) then
                if (     (dqmax_old*dqmax > 0.0_dp &
                  &  .and. abs(dqmax_old) > abs(dqmax)) &
                &  .or. (dqmax_old*dqmax < 0.0_dp)) then
                    aa=aa*abs(dqmax_old)/abs(dqmax_old-dqmax)
                end if

                aa=max(aamin,aa)
                aa=min(aamax,aa)

                if (abs(dqmax) <= qerr) aa=aamin


!                 print *, redb, 'aa:', aa, endc
            
                dec(:mbc%nd) = de(:mbc%nd) + aa*dq(:mbc%nd)

                dqmax_old = dqmax
            else
!                 if (mbc%mxgf) de(:mbc%n(1)) = de(:mbc%n(1)) + mbc % dqchi(:mbc%n(1)) 
                de(mbc%n(1)+1:mbc%nd) = de(mbc%n(1)+1:mbc%nd) + mbc % dqchi(mbc%n(1)+1:mbc%nd) 
!                 print *,iproc, 'TB DE',de(mbc%n(1)+1:mbc%nd)

!             endif
            
            
!        
            if (mbc%n(2)/= 0) then
                        
                do  i = mmix, 0, -1
                    call dcopy(nat,qmx(1,i,1),1,qmx(1,i+1,1),1)
                enddo
                
                
                
!                 if (.not. quiet) print *, 'qmpol before mix',mbc%qmpol(1,:)
!  This is done to isolate the charges of the LCN orth atoms from the SCTB others
                i = 1
                do ia = 1,mbc%n(2)
                    if (mbc%orth(ia) == 0) then
                        qmx(i,0,1) = mbc % qmpol(1,ia)
                        i = i + 1
                    endif
                enddo
!                 print *, mbc%qmpol(1,:)
!                 
! !                 qmx(:,0,1) = mbc % qmpol(1,:)
!                 print *, 'qm1'
!                 
!                 
!                 
!                 do i = 1, nat
!                     print *, i,qmx(i,:,1)
!                 enddo
!                 print *, 'qm2'
!                 do i = 1, nat
!                     print *,i,qmx(i,:,2)         
!                 enddo

                
                call qmix(qmx,nat,nsp,nlmq,rms,it,mmix,beta,wc,.true.)
                
!                 print *, 'rms,mmix,beta,wc'
!                 print *,  rms,mmix,beta,wc
                
                do  i = mmix, 1, -1
                    call dcopy(nat,qmx(1,i,2),1,qmx(1,i+1,2),1)
                enddo
                
                  
                
                qmx(:,1,2) = qmpol_prev(:)
                qmpol_prev(:) = qmx(:,0,2)
                
!                 mbc % qmpol(1,:) = qmx(:,0,2)

                mbc % qmpol(1,:) = 0.0_dp
                i = 1
                do ia = 1,mbc%n(2)
                    if (mbc%orth(ia) == 0) then
                        mbc % qmpol(1,ia) = qmx(i,0,2)
                        i = i + 1
                    endif
                enddo

                
!                 print *, 'qm2'
!                 do i = 1, nat
!                     print *,i,qmx(i,:,2)         
!                 enddo
!                 
!                 print *, 'qmpol'
!                 print *, mbc%qmpol(1,:)
                
    !                     do i = 1, mbc%n(2)
    !                         mbc % dqchi(i+mbc%rlst(2)-1) = mbc % qmpol(1,i)/mbc % dqchi(i+mbc%rlst(2)-1)
    !                     enddo
            
                
                do i = 1,nsp
!                     print *,'isp',i
                    mbc%h(:,:,i) = mbc % h0(:,:)
                enddo
                
                if (forces .and. (abs(dqmax) <= qerr .or. it == mxit)) then
                    if (mag) then
                        mgmax = maxval(abs(mg - mg_o))
                        if (mgmax <= merr .or. mit == mgit) calF = .true.
                    else
                        calF = .false.
                    endif
                endif
                
                call inctbh(ctsx,mbc,ef,calF,mg_o,nlmq,ecorr,it)
                
               endif
            endif

            
      !                      dec(:nd) = de(:nd)
  
            if (iproc == master) print *, 'dqmax:',dqmax
!             print *, redb, 'dqmax:', dqmax, endc
            if (abs(dqmax) <= qerr .or. it == mxit) exit


            

            it = it + 1
            
         end do

         if (.not. mag) exit
         mgmax = maxval(abs(mg - mg_o))

!         if (.not. quiet) then
!            write(6,ffmt) 'mgo', mg_o
!            write(6,ffmt) 'mg ', mg
!         end if
         if (iproc == master) print *,  'mgmax:', mgmax
         
         if (mgmax <= merr .or. mit == mgit) exit
         mit = mit + 1
      end do

      if (mag) deallocate(mg_o)
      deallocate(dec)

      
     if (mbc%mxgf) then
        efb = ef
     else
        mbc % apb = 0
        mbc % bpb = 0
        if (sum(rlim) == 0) then
            nat = mbc% nd
        else
            nat = mxtotnd-1
        endif
        call extract_boplist(aptr, bptr,mbc % rlst, mbc % apb, mbc % bpb,nat,.true.)
     endif
     
!      print *, 'final apb:',mbc%apb
!      print *, 'final bpb:',mbc%bpb(:mbc%apb(nat+1)+30)

      call mbenergy(mbc, ef,efb,ecorr)

      if (forces) call mbf(mbc,ef,efb)
      
      if (rtc % vpair /=0) call mbpair(mbc) 
      
      
!     if (iproc == master) then
! !      print *, 'main iproc=',iproc,'master=',master
!       print *, 'FORCES'
!       Print *, 'BOP ATOMS'
!       do ia = 1,mbc%n(1)
!           print *, 'ATOM:',ia
!           print *, 'COORD:',ad(:,ia)
!           print *, 'FBS:',fbs(:,ia)
!           print *, 'FPP',fpp(:,ia)
!           print *, 'FTOT:',fbs(:,ia)+fpp(:,ia)
!           print *,' '
!       enddo
!       print *,' '
!       print *,' '
!       PRINT *, 'TB ATOMS'
!       do ia = 1,mbc%n(2)
!           git = ia + mbc%n(1)
!           print *, 'ATOM:',git
!           print *, 'COORD:',ad(:,git)
!           print *, 'FBS:',fbs(:,git)
!           print *, 'FPP',fpp(:,git)
! ! !           IF(MBC%SCTB) then
! ! !               print *, 'FOVL', fnou(:,ia)+fnom(:,ia)
! ! !               print *, 'FTOT:',fbs(:,git)+fpp(:,git)+fnou(:,ia)+fnom(:,ia)
! ! !           else 
!               print *, 'FTOT:',fbs(:,git)+fpp(:,git)
! !           endif
!           print *,' '
!       enddo
!       endif
      
   end subroutine emb_sc



   subroutine mbf(mbc,ef,efb)
      use mod_atom_ar
      use mod_precision
      use mod_all_scalar
      use mod_const
      use topologia, only : iproc, mpmap, aproc, master, nproc
      use mpi
     
     
      include "../Include/Atom.array"
      include "../Include/Force.array"
      include "../Include/PosVel.array"
      
      include "../Include/NebList.array"    
     
      type(mbtb_t), intent(inout) :: mbc
      real(dp), intent(in) :: ef,efb
      real(dp),allocatable :: eb(:,:,:,:),f(:,:)
      integer :: isp, la, ia,i

      complex(dp) :: zp,zep,zepb,zpn
      real(dp) :: phase,w0
!       ,fnou(3,mbc%n(2)),fnom(3,mbc%n(2))
      integer :: m, p,ierr
      
      integer :: gib,ib,it,git,ja,ja0
!       real(dp) ::  vm(9,mbc%n(2)),vu(mbc%n(2))


!      M     = NINT(MFAC*0.25D0*(EF-LAINF(LA)+2.0D0*LBINF(LA))/KT)
!        m     = nint(mfac*0.25_dp*(ef-lainf(la)+2*lbinf(la))/kt)
!       print *,ef, lainf(la), lbinf(la), (ef-lainf(la)+2.0_dp*lbinf(la))
!       m     = nint(mfac*0.25_dp*1.8_dp/kt)
      m = mfac
!  When ovl==true the tb part requires more neighbours due to greater range of s-orb
      if (.not. mbc % ovl) then
        if (.not. allocated(mbc%rho)) allocate(mbc % rho(mxnstat, mxnstat,(mxnnb+1) , mbc % nd))
      else
        if (.not. allocated(mbc%rho)) allocate(mbc % rho(mxnstat, mxnstat,80 , mbc % nd))
                  
        if (.not. allocated(eb))  allocate(eb(mxnstat, mxnstat,80 , mbc % n(2)))
        
      endif
      
      if (.not. allocated(f))  allocate(f(3, mbc % n(2)))
      
      f(:,:) = fbs(:,mbc%n(1)+1:mbc%nd)
      fbs = 0.0_dp
      
      rdfbs = 0.0_dp

      do isp = 1, nsp         
          mbc % rho(:,:,:,:) = 0.0_dp
          if(mbc%ovl) eb(:,:,:,:) = 0.0_dp
          w0    = kt*(2*m)
          phase = pi/real(2*m, dp)
          zp    = cmplx(cos(phase),sin(phase), kind=dp)

          do p = mpmap(iproc), mpmap(iproc+1)-1
              zpn = zp**(2*p+1)
              zep  = ef+w0*(zpn-1)
              zepb = efb+w0*(zpn-1)
              
              call gg_rho(mbc,eb, zep,zepb, zpn, isp)
          enddo

          call mbforce(mbc,eb)
      
      end do
      
      if (nproc > 1) then
          call mpi_allreduce(mpi_in_place,  fbs, 3*mbc%nd, mpi_real8, mpi_sum, mpi_comm_world, ierr)
          call mpi_allreduce(mpi_in_place, rdfbs,   1, mpi_real8, mpi_sum, mpi_comm_world, ierr)

      end if
      
      fbs(:,mbc%n(1)+1:mbc%nd) = fbs(:,mbc%n(1)+1:mbc%nd) + f(:,:)
      
      if (allocated(eb))  deallocate(eb)
      
!       if (iproc == master) then
! !         print *, 'iproc=',iproc,'master=',master
! !         print *, 'FORCES'
! !         do ia = 1, mbc%nd
! !             print *,'ATOM:',ia
! !             ja = aptr(ia)+1
! !             ib = bptr(ja)
! !             ja0=ja-2
! !             do while (ib /= eol)
! !                 print *,'    NEIGH:',ib,' F: ',force(:,ia,ib)
! !         
! !                 ja = ja + 1
! !                 ib = bptr(ja)
! !             end do
! !         enddo
! ! !         do ia = 1, mbc%n(1)
! ! !             print *,'bpATOM:',ia
! ! !             ja = mbc % apb(ia)+1
! ! !             ib = mbc % bpb(ja)
! ! !             ja0=ja-2
! ! !             do while (ib /= eol)
! ! !                 print *,'    NEIGH:',ib,' F: ',force(:,ia,ib)
! ! !         
! ! !                 ja = ja + 1
! ! !                 ib = mbc % bpb(ja)
! ! !             end do
! ! !         enddo
! ! !         
! ! !         do ia = 1, mbc%n(2)
! ! !             git = mbc % rlst(2) -1 + ia
! ! !             print *,'tbATOM:',git
! ! !             ja = mbc % apt(ia)+1
! ! !             ib = mbc % bpt(ja)
! ! !             ja0=ja-2
! ! !             do while (ib /= eol)
! ! !                 gib = mbc % rlst(2) -1 + ib
! ! !                 print *,'    NEIGH:',gib,' F: ',force(:,git,gib)
! ! !         
! ! !                 ja = ja + 1
! ! !                 ib = mbc % bpt(ja)
! ! !             end do
! ! !         enddo
! ! !         
! ! !         do ia = 1, mbc%n(2)
! ! !             git = mbc % rlst(2) -1 + ia
! ! !             print *,'dATOM:',git
! ! !             ja = mbc % apd(ia)+1
! ! !             ib = mbc % bpd(ja)
! ! !             ja0=ja-2
! ! !             do while (ib /= eol)
! ! !                 print *,'    NEIGH:',ib,' F: ',force(:,git,ib)
! ! !         
! ! !                 ja = ja + 1
! ! !                 ib = mbc % bpd(ja)
! ! !             end do
! ! !         enddo
!         print *, 'band forces'
!         do ia = 1, mbc%n(2)
!            print *, ia
!            print *, 'ham:',fham(:,ia)
!            print *, 'ovl:',fovl(:,ia)
!            print*, 'tot:',(fbs(:,ia+mbc%n(1))-(fsc(:,ia,1)+fsc(:,ia,2)+f(:,ia)))
! 
!         enddo
!         do ia = 1, mbc%n(2)
!           print *,'FORCES ON ATOM:',ia
!           print *, 'from bands   :  ', (fbs(:,ia+mbc%n(1))-(fsc(:,ia,1)+fsc(:,ia,2)+f(:,ia)))
!           print *, 'from estrx   :  ', (f(:,ia))
!           print *, 'from overlap (U) :  ',fsc(:,ia,1)
!          print *, 'from overlap (Mad) :  ',fsc(:,ia,2)
!          print *, ' '
!         enddo
! !         
! ! !         print *, 'SUM'
! ! !         do i =1,mbc%nd
! ! !             print *, i,fbs(:,i)
! ! !         enddo
!       endif
      
      
   end subroutine mbf


!    subroutine mbrho(mbc, ef)
!       use mod_all_scalar, only : kt, mfac, nsp
!       use mod_atom_ar
! 
!       type(mbtb_t), intent(inout) :: mbc   
!       real(dp), intent(in) :: ef
!       integer :: isp, la, ia
! 
!       real(dp) :: cnel
! 
!       complex(dp) :: zp,zfac,zep
!       real(dp) :: phase,w0
!       integer :: m, p
! 
! 
! !      M     = NINT(MFAC*0.25D0*(EF-LAINF(LA)+2.0D0*LBINF(LA))/KT)
! !        m     = nint(mfac*0.25_dp*(ef-lainf(la)+2*lbinf(la))/kt)
! !       print *,ef, lainf(la), lbinf(la), (ef-lainf(la)+2.0_dp*lbinf(la))
! !       m     = nint(mfac*0.25_dp*1.8_dp/kt)
!       m = mfac
!         if (.not. allocated(mbc%rho)) allocate(mbc % rho(mxnstat, mxnstat,(mxnnb+1) , mbc % nd))
!       mbc % rho = 0.0_dp
! 
!       do isp = 1, nsp         
! 
!          w0    = kt*(2*m)
!          phase = pi/real(2*m, dp)
!          zp    = cmplx(cos(phase),sin(phase), kind=dp)
!          zfac  = zp*zp
! 
!          do p = 0, m-1
! !             print *, 'p', p
!             zep  = ef+w0*(zp-1)
!             call gg_rho(mbc, zep, zp, isp)
!             zp   = zfac*zp
!          enddo
! 
!       end do
! 
!    end subroutine mbrho
   

   subroutine gg_diag(mbc, ez,isp)
      use ab_io
      use mod_ham

      type(mbtb_t), intent(inout) :: mbc
      complex(dp), intent(in) :: ez
      integer, intent(in) :: isp
      integer :: i, j,ia,it,gia
      complex(dp) :: gab(mxnstat, mxnstat)
      
      call a_and_ac(mbc, ez, isp)
      
      if (mbc%mxgf) then
!       calc the diagonal here; do over atoms and orbs  gb(ia,ia) +-  sum_ib(ac(ia,ib) * a(ia,ib))
          do i = 1, mbc % bhsz  
              do j = 1, mbc % thsz
                  mbc % ggd(i) = mbc % ggd(i) + mbc % ac(i,j)*(mbc % a(i,j))
              end do
          end do
      endif

      do i = 1, mbc % thsz
          mbc % ggd(mbc % bhsz + i) = mbc % gt(i,i)
      end do
      
   end subroutine gg_diag   



   subroutine gg_rho(mbc,eb, ez,ezb, zp, isp)
      use ab_io
      use mod_ham
      
      type(mbtb_t), intent(inout) :: mbc
      real(dp),allocatable :: eb(:,:,:,:)
      complex(dp), intent(in) :: ez,ezb, zp
      integer, intent(in) :: isp

      integer :: ia, gia, ib, gib, it, git, ja,ja0, jb, jt, i, j, jt0,jtl

      complex(dp) :: gab(mxnstat, mxnstat)

!      complex(dp) :: gab1(mxnstat, mxnstat),gab2(mxnstat, mxnstat)

      include "../Include/NebList.array"

    call a_and_ac(mbc, ez, isp)
     

    do ia = 1, mbc % n(1)
!          print *, ia
        call assoc_ab(ia,isp)
        call assoc_ham(ia) ! for cluster, pcluster and decipher

        gia = mbc % rlst(1) -1 + ia
                                
!          BOP onsites
        gab = 0.0_dp
        
!     STILL NOT SURE WHY IM CALCULATING THIS (the onsite density matrix)
        call get_gab(ezb, ia, ia, 1, gab,.false.)

        if (mbc%mxgf) then
            do i = 1, mbc % pbh(ia+1) - mbc % pbh(ia)
                do j = 1, mbc % pbh(ia+1) - mbc % pbh(ia)
                    gab(i,j) = gab(i,j) + sum(mbc % ac(mbc % pbh(ia) + i, :) *(mbc % a(mbc % pbh(ia) + j, :)))
                end do
            end do
        endif
        mbc % rho(:,:,1,gia) = mbc % rho(:,:,1,gia) + real(zp*gab)

!          BOP intersites

        ja = mbc % apb(ia)+1
        ib = mbc % bpb(ja)
        ja0 = ja - 2

        do while (ib /= eol)
!          print *, "BOP",ia,ib
            jb = decipher(ib)
            if (jb /= 1) then !jb == 1 => ib == ia
                gib = mbc%rlst(1)-1 + ib !

                gab = 0.0_dp
!            gab1 = 0.0_dp
!            gab2 = 0.0_dp
                call get_gab(ezb, ia, ib, jb, gab,.false.)

                if (map(ib) == ib .and. mbc%mxgf) then
                    do i = 1, mbc % pbh(ia+1) - mbc % pbh(ia)
                        do j = 1, mbc % pbh(ib+1) - mbc % pbh(ib)
                            gab(i,j) = gab(i,j) + sum(mbc % ac(mbc % pbh(ia) + i, :) * (mbc % a(mbc % pbh(ib) + j, :)))
                        end do
                    end do
                endif 
                
                mbc % rho(:,:,ja-ja0,gia) = mbc % rho(:,:,ja-ja0,gia) + (real(zp*gab))
                
            end if
            ja = ja + 1
            ib = mbc % bpb(ja)
        end do

    end do

    do it = 1, mbc % n(2)

        git = mbc % rlst(2) -1 + it

!          TB all
!          jt = mbc % apt(it)+1 IM NOT SURE IF I SHOULD BE DOING ONSITES
    
    
        jt = mbc % apt(it)
        jt0 = jt - 1
        ib = mbc % bpt(jt)

        do while (ib /= eol)
            gib = mbc % rlst(2) -1 + ib
!             print *, "tb",it,jt-jt0,gib
            gab = 0.0_dp
            gab(1 : mbc % pth (it + 1) - mbc % pth(it), 1 : mbc % pth (ib + 1) - mbc % pth(ib)) =  &
            & mbc % gt(mbc % pth(it) + 1 : mbc % pth(it + 1), mbc % pth(ib) + 1 : mbc % pth(ib + 1))
            
            mbc % rho(:,:,jt-jt0,git) = mbc % rho(:,:,jt-jt0,git) + real(zp*gab)
            
            if (mbc % ovl) eb(:,:,jt-jt0,it) = eb(:,:,jt-jt0,it)  + real(zp*ez*gab)
            

            jt = jt + 1
            ib = mbc % bpt(jt)
            jtl = jt-(jt0+1)
        end do
!          TB->BOP intersites

            if (.not. mbc%mxgf) then    
                call assoc_ab(git,isp)
                call assoc_ham(git) ! for cluster, pcluster and decipher
            endif
            
            jt = mbc % apd(it)
            jt0 = jt - 1
    !        ib = mbc % bpd(ja)
            ib = mbc % bpd(jt)
    ! 	 print *, "TBBOP jt:",jt
            do while (ib /= eol)
            
    ! should this be changed so that the global verison is rlst(1)-1+ib as they are BOP atoms?
                gib = mbc%rlst(1)-1 + ib
    !             print *, "delta",git,(jt-jt0)+jtl,ib
                gab = 0.0_dp
                
                
                if (mbc%mxgf) then 
                    gab(1 : mbc % pth (it + 1) - mbc % pth(it), 1 : mbc % pbh (ib + 1) - mbc % pbh(ib)) =  &
                    & - (transpose(mbc % ac(mbc % pbh(ib) + 1 : mbc % pbh(ib + 1), mbc % pth(it) + 1 : mbc % pth(it + 1))))
                else
                    jb = decipher(ib)
                    call get_gab(ez, git, ib, jb, gab,.false.)
                endif
                
    !           if (it == 1 .and. ib == 22) then
    !               do i=1,mbc % pbh (ib + 1) - mbc % pbh(ib)
    !                   print *, gab(:,i)
    !               enddo
    !           endif
                mbc % rho(:,:,(jt-jt0)+jtl,git) = mbc % rho(:,:,(jt-jt0)+jtl,git) + real(zp*gab)
                
                jt = jt + 1
                ib = mbc % bpd(jt)


            end do

    end do

   end subroutine gg_rho    
   
   
   subroutine mbforce(mbc,eb)
      use mod_precision
      use mod_all_scalar
      use mod_const
      use ab_io
      use mod_ham, only : cluster, pcluster, assoc_ham ! assoc_ham for clustr and pcluster
      use mod_chi
      use mod_funptr
!       use topologia, only : iproc, mpmap, aproc, master, nproc

!       use mod_atom_ar, only : btype

      use mpi
   !
   !    This is a subroutine to evaluate the bond contribution to the
   !     band structure energy, and the band structure forces.
   !

      implicit none


      include "../Include/Atom.array"
      include "../Include/Force.array"
      include "../Include/Moment.array"
      include "../Include/NebList.array"
      include "../Include/PosVel.array"

      real(dp) :: bondeng,hopint
      type(mbtb_t), intent(inout) :: mbc
      real(dp),allocatable :: eb(:,:,:,:)
      integer nmax,isp
      integer nla,nlb,n
      integer la,lb, mib
      integer nstta,nsttb
      integer za,zb
      integer lla0,nma,ma
      integer :: nmb, llb0, mb, nmaxb, mlistb(mxnstat)
      integer :: ia, gia, ib, gib, it, git
      integer :: ja,ja0, jb, jb2, jt, i, j, jt0,jtl
      type(ham_conf_t),pointer :: hamtb,hambp

   !
   !    Declare the local arrays.
   !

   !*** A displacement vector.
      real(dp) :: dr(3)
   !*** Forces.
      real(dp) :: fab(3)

      real(dp) :: grad(3,mxnstat,mxnstat), sgrad(3,mxnstat,mxnstat)

      integer :: ierr,off1,off2

      real(dp) :: scfcut(14),dsf(14,3)
      real(dp) ::  drhos(3)
! vm(9,mbc%n(2)),vu(mbc%n(2)),fnou(3,mbc%n(2)),fnom(3,mbc%n(2))
      scfcut = 1.0_dp
      dsf = 0.0_dp

      
      hambp => rtc % ham_conf 
!       print *, 'COORD:'
      do ia = 1, mbc % n(1)

!       should this assoc_ab be here? there is no value for isp
!         call assoc_ab(ia,isp)
!             call assoc_chi(ia,isp)
!             call assoc_ham(ia)

            za = z(ia)
            nstta = nstt(za)
!             call states(za,nla,nstta,llista)

            
!             I don't think any of this stuff is actually used
!             if (momflg == 1) then
!                lla0 = 0
!                do la = 1,nla
!                   nma = 2*llista(la) + 1
!                   mlista(lla0+1:lla0+nma) = la
!                   lla0 = lla0 + nma
!                enddo
!             else
!                do la = 1, nstta
!                   mlista(la) = la
!                end do
!             endif      
!       
         gia = mbc % rlst(1) -1 + ia
                                    
!          BOP onsITES
! DO I NEED ONSITE FORCES?

!          BOP intersites
!         print *, ia,ad(:,ia)

      ja = mbc % apb(ia)+1
      ib = mbc % bpb(ja)
!        NOT SURE IF THIS JA IS CORRECT
      ja0=ja-2
      do while (ib /= eol)
            gib = mbc % rlst(1) -1 + ib
            
  !                 mib = map(ib)
            zb = z(ib)
!              call states(zb,nlb,nsttb,llistb)
            nsttb = nstt(zb)

            dr = ad(1:3,ia) - ad(1:3,ib)
            call grdmat(grad,dr,za,nstta,zb,nsttb,scfcut,dsf,hambp,.false.)
         

!          do n =1,3
!             print *, 'GRAD',n
!             do j = 1, mxnstat
!                 print *, (grad(n,i,j),i = 1, mxnstat)
!             enddo
!          enddo
               
!                 if (momflg == 2) call utranv(transm(1,1,ia),grad,transm(1,1,mib),trwork,mxnstat,nstta,nsttb *CHECK WHAT THIS WAS FOR*
            fab  = 0.0_dp
               

            do i = 1, mbc % pbh(ia+1) - mbc % pbh(ia)
                do j = 1, nsttb
                  fab=fab+mbc % rho(i,j,(ja-ja0),gia)*grad(1:3,i,j)
!                   print *, mbc % rho(i,j,(ja-ja0),gia), grad(1:3,i,j)
                end do
            end do
            
            if (.not. mag) then
              fab = 2*fab
            end if

            fab=-2*kt*fab
            
            if (map(ib) == 0) then
	            fbs(:,gia) = fbs(:,gia) + 2*fab
            else
            	fbs(:,gia) = fbs(:,gia) + fab
                fbs(:,map(ib)) = fbs(:,map(ib)) - fab
                
                
            endif
            
!             print *, "BOP,ia,ib", ia,ib
!             print *, 'dr',dr
!             print *, fab
!             print *, ' '
            
            rdfbs = rdfbs - sum(dr*fab)
            
            ja = ja + 1
            ib = mbc % bpb(ja)
         end do

      end do
      
      
      hamtb => rtc % tbham_conf 
!       print *, 'TB atoms'
      do it = 1, mbc % n(2)
         git = mbc % rlst(2) -1 + it
         za = z(git)
         nstta = mbc % tbc % nstt(za)
!          llista = mbc % tbc % llist(:,za)
!          nla = mbc % tbc % nl(zb)
!          call states(za,nla,nstta,llista)
!          TB all
!          jt = mbc % apt(it)+1 IM NOT SURE IF I SHOULD BE DOING ONSITES
        
!             print *, git,ad(:,git)
            
            jt = mbc % apt(it)+1
            jt0 = jt - 2
            ib = mbc % bpt(jt)

            do while (ib /= eol)
                gib = mbc % rlst(2) -1 + ib
    !             print *, 'neigh =',ib,'jt=',jt,'jt-jt0=',jt-jt0
    !             mapb = map(gib)
!     
!                  if (mbc%mxgf .or. .not. mbc%allbop) then
                zb = z(gib)
                nsttb = mbc % tbc % nstt(zb)
    !             llistb = mbc % tbc % llist(:,zb)
    !             nlb = mbc % tbc % nl(zb)
    !             call states(zb,nlb,nsttb,llistb)
                dr = ad(1:3,git) - ad(1:3,gib)
                
                call grdmat(grad,dr,za,nstta,zb,nsttb,scfcut,dsf,hamtb,.false.)
                if (mbc % ovl) call grdmat(sgrad,dr,za,nstta,zb,nsttb,scfcut,dsf,hamtb,.true.)
        
    !             do i = 1,3
    !                 gradt(i,:mbc%pth(it+1)-mbc%pth(it), :mbc%pth(ib+1)-mbc%pth(ib)) =    &
    !                 &             grad(i,1+mbc%orth(it):nstta,1+mbc%orth(ib):nsttb)
    !                 if (mbc % ovl) sgradt(i,:mbc%pth(it+1)-mbc%pth(it), :mbc%pth(ib+1)-mbc%pth(ib)) =    &
    !                 &             sgrad(i,1+mbc%orth(it):nstta,1+mbc%orth(ib):nsttb)
    !             enddo
                
                ! TO CUT OUT THE sd PART OF GRAD
!                 n = mbc%orth(it) + mbc%orth(ib)
!                 if (n /= 0 .and. n /= 2) then
!                     if (mbc%orth(it) == 1) then
!                         grad(1:3,:,1) = 0.0_dp
!                         if (mbc % ovl) then
!                             sgrad(1:3,:,1) = 0.0_dp
!                         endif
!                     else
!                         grad(1:3,1,:) = 0.0_dp
!                         if (mbc % ovl) then
!                             sgrad(1:3,1,:) = 0.0_dp
!                         endif
!                     endif
!                 endif

                off1 = 0
                off2 = 0
                if ( mbc%orth(it) + mbc%orth(ib) /= 0) then
                    off1 = nstta - (mbc % pth(it+1) - mbc % pth(it))
                    off2 = nsttb - (mbc % pth(ib+1) - mbc % pth(ib))
                    grad(1:3,:nstta-nstt(za),:) = 0.0_dp
                    grad(1:3,:,:nsttb-nstt(zb)) = 0.0_dp
                    if (mbc%ovl) then
                        sgrad(1:3,:nstta-nstt(za),:) = 0.0_dp
                        sgrad(1:3,:,:nsttb-nstt(zb)) = 0.0_dp
                    endif
                endif

                fab=0.0_dp
                
                
                drhos = 0.0_dp
                
!                 print *,'frce, it,ito,ib,ibo',it,mbc%orth(it),ib,mbc%orth(ib)
!                 print *, 'rho'
!                 do i = 1, mbc % pth(it+1) - mbc % pth(it)
!                   print *,mbc%rho(i,:,jt-jt0,git)
!                 enddo

!                 if (.not. quiet) then
!                      print *, it,ib
!                     print *,'rho'
!                     do i = 1, mbc % pth(it+1) - mbc % pth(it)
!                       print *,mbc%rho(i,:,jt-jt0,git)*4*kt
!                     enddo
!                     
!                     if (mbc% ovl) then
!                         print *,'eb',off1,off2
!                         do i = 1, mbc % pth(it+1) - mbc % pth(it)
!                           print *,eb(i,:,jt-jt0,it)*4*kt
!                         enddo
!                         
!                  !       print *,'hgrad'
!                  !       do j = 1,3
!                  !           do i = 1, mbc % pth(it+1) - mbc % pth(it)
!                  !             print *,grad(j,i,:)
!                  !           enddo
!                    !     enddo
!                         
!                    !     print *,'sgrad'
!                    !     do j = 1,3
!                    !         do i = 1, mbc % pth(it+1) - mbc % pth(it)
!                   !            print *,sgrad(j,i,:)
!                   !          enddo
!                   !      enddo
!                     endif
!                     
!                     
!                 endif
                
                do i = 1, mbc % pth(it+1) - mbc % pth(it)
                  do j = 1, mbc % pth(ib+1) - mbc % pth(ib)
!                      PRINT *, i,j,'rho', mbc%rho(i,j,jt-jt0,git),'grad',(grad(1,i+off1,j+off2))
                    fab = fab + mbc%rho(i,j,jt-jt0,git) * (grad(1:3,i+off1,j+off2))
                    if (mbc % ovl) then
                        fab = fab - (eb(i,j,jt-jt0,it) * sgrad(1:3,i+off1,j+off2))
                        if (mbc % sctb) drhos(1:3)=drhos(1:3) + mbc%rho(i,j,jt-jt0,git) * sgrad(1:3,i+off1,j+off2)
                    endif
    !                  print *, 'rho', mbc%rho(i,j,jt-jt0,git) * grad(1:3,i,j)
                  end do
                end do            
                
                if (.not. mag) then
                  fab = 2*fab
                end if

                fab=-2*kt*fab
!                 if (.not. quiet) then
!                     print *, isp,"TB,git,gib", git,gib
! !             print *, 'dr',dr
!                     print *, fab
!                     print *, ' '
!                 endif 
                
                if (mbc%sctb .and. mbc % ovl) then
    !       Producing Hubbard and Madelung contibutions to forces from overlap.
                    drhos =  drhos*2*kt
                    
                    if (.not. mag) drhos =  drhos*2
                    
!                     fnou(1:3,it) = fnou(1:3,it) - (mbc%vu(it) + mbc%vu(ib)) * drhos(1:3)
!                     fnom(1:3,it) = fnom(1:3,it) - (mbc%vm(1,it) + mbc%vm(1,ib)) * drhos(1:3)

                      fbs(:,git) = fbs(:,git) - (mbc%vu(it) + mbc%vu(ib)) * drhos(1:3)
                      fbs(:,git) = fbs(:,git) - (mbc%vm(1,it) + mbc%vm(1,ib)) * drhos(1:3)
                    
!                     if (.not. quiet) then
!                        
! !                         print *,'u',mbc%vu(it),mbc%vu(ib)
! !                         print *, 'm', mbc%vm(1,it) , mbc%vm(1,ib)
!                         print *, drhos(1:3)
!                     endif
                    
                endif
                
                fbs(:,git) = fbs(:,git) + fab
                fbs(:,gib) = fbs(:,gib) - fab
                

                rdfbs = rdfbs - sum(dr*fab)
                jt = jt + 1
                ib = mbc % bpt(jt)
                jtl = jt-(jt0+1)
                

            end do
            
              
!             if (mbc%sctb) fbs(:,git) = fbs(:,git) + fnou(:,it) + fnom(:,it)

          
!            if (mbc%sctb) then
!               fsc(:,it,1) = fsc(:,it,1)  + fnou(:,it)
!               fsc(:,it,2) = fsc(:,it,2) + fnom(:,it)
!           endif
    !          TB->BOP intersites
            jt = mbc % apd(it)
            jt0 = jt - 1
    !          LOOK INTO WHETHER THIS IS RIGHT TOO
            
            call states(za,nla,nstta,llista)
            ib = mbc % bpd(jt)
              
            
            do while (ib /= eol)
            
                gib = mbc%rlst(1)-1 + ib
                
    !             mapb = map(gib)
                zb = z(gib)
                call states(zb,nlb,nsttb,llistb)
                dr = ad(1:3,git) - ad(1:3,gib)
                grad = 0.0_dp
    !             print *, 'ia,ib,nstta,nsttb,za,zb',ia,ib,nstta,nsttb,za,zb
                call grdmat(grad,dr,za,nstta,zb,nsttb,scfcut,dsf,hambp,.false.)
                
                
                fab=0.0_dp
                
                do i = 1, nstta
                  do j = 1, nsttb
                    fab = fab + mbc%rho(i,j,(jt-jt0)+jtl,git) * grad(1:3,j,i)
                  end do
                end do
                
                if (.not. mag) then
                  fab = 2*fab
                end if

                fab=-2*kt*fab
                
!                 if (.not. quiet) then
!                     print *, "DH,git,ib", git,ib
! !             print *, 'dr',dr
!                     print *, fab
!                     print *, ' '
!                 endif
                
                fbs(:,git) = fbs(:,git) + 2*fab
                if (map(gib) /= 0) fbs(:,map(gib)) = fbs(:,map(gib)) - 2*fab
                
!                 force(:,git,map(gib)) = force(:,git,map(gib)) + 2*fab
!                 force(:,map(gib),git) = force(:,map(gib),git) - 2*fab
                
                rdfbs = rdfbs - sum(dr*fab)
                
    !             print *, 'rho', mbc%rho(:,:,(jt-jt0)+jtl,git)
                jt = jt + 1
                ib = mbc % bpd(jt)
                
    !             print *, fab
            end do
      end do
      
      

      
!       if (iproc == master) then
!           print *, 'FORCES'
!           do ia = 1, mbc%nd
!               print *, ia, fbs(:,ia)
!           enddo
!       endif
   
   end subroutine mbforce
   

   subroutine mbenergy(mbc, efermi,efermb,ecorr)
      use ab_io
      use mod_ham
      use mod_precision
      use mod_const
      use mod_funptr
      use mod_all_scalar
      use topologia, only : iproc, mpmap, aproc, master, nproc
      use mpi
      use mod_clock 
      
      integer(8) :: c1, c2, c3, c4, c5
      real(dp), intent(in) :: efermi,efermb
      type(mbtb_t), intent(inout) :: mbc
      integer :: isp,nsttb
      procedure(real(dp)) :: entropy
!       , femag
      integer :: ia, gia, ib, gib, it, git, ja, jb,jt, i, j, jt0,ierr,m,p

      complex(dp) :: gab(mxnstat, mxnstat),zp,zep,zepb,zpn
      complex(dp), allocatable :: si(:,:),gti(:,:)
!,a2(mbc%thsz,mbc%bhsz),gab1(mxnstat, mxnstat),gab2(mxnstat, mxnstat)
       real(dp) :: bsons, error,phase,w0,dea,ecorr

      
      integer :: info,liwork,lwork
      integer, allocatable :: iwork(:)
      real(dp), allocatable :: cnw(:,:),ovlm(:,:),w(:),work(:),q(:),ql(:,:)
      
!      real(dp), allocatable :: buf(:,:)
      real(dp) :: buf(4)
!      integer, allocatable :: req(:)
!      integer, allocatable :: stat(:,:)
!      integer :: n
      
      complex(dp) :: z_one,z_zero
       
       integer :: llistt(3),nlt,no,lla,zt
       
       
       real(dp) :: eprom0
       
       
      include '../Include/Atom.array'
      include "../Include/NebList.array"


      call system_clock(c1)
      
      if (mbc%ovl) then
          if (.not. allocated(si)) allocate(si(mbc%thsz,mbc%thsz))
          if (.not. allocated(gti)) allocate(gti(mbc%thsz,mbc%thsz))
      endif
      
      if (.not. allocated(q)) allocate(q(mbc%n(1)))
      if (.not. allocated(ql)) allocate(ql(3,mbc%n(2)))
      
!      if (nproc > 1) then
!          n = 4
!
!          if (iproc /= master) then
!             allocate(buf(n,1),req(1),stat(mpi_status_size,1))
!          else if (aproc>1) then
!             allocate(buf(n,aproc-1),req(aproc-1),stat(mpi_status_size,aproc-1))
!             do i=1,aproc-1
!                call mpi_irecv(buf(1,i), n, mpi_real8, i, i, mpi_comm_world, req(i), ierr )
!             end do
!          end if
!      endif
      
      eprom = 0.0_dp
      eprom0 = 0.0_dp
      ebond = 0.0_dp
      uent  = 0.0_dp
      bsons = 0.0_dp
      if (.not. mbc%sctb) emag  = 0.0_dp
      
      do isp = 1, nsp
          q = 0.0_dp
          ql = 0.0_dp
      
!  Calculate TB GF from expansion coefficients in Dyson case          
          if (.not. mbc% mxgf) then
              call system_clock(c3)
              allocate(cnw(mbc%thsz,mbc%thsz),w(mbc%thsz))
              
              if (mbc%ovl) allocate(ovlm(mbc%thsz,mbc%thsz))
              allocate(work(4))
              work = 1.0_dp
              lwork = -1
              
              cnw(:,:) = mbc%h(:,:,isp)
              
              if (mbc% ovl) then
                  allocate(iwork(4))
                  iwork = 1.0_dp
                  liwork = -1
                  
                  ovlm(:,:) = mbc%s(:,:)
                  call dsygvd(1,'V','U',mbc%thsz,cnw,mbc%thsz,ovlm,mbc%thsz,w,work,lwork,iwork,liwork,info) 
                  if (info /= 0) stop 'dsygvd i info /= 0'
              else
                  call dsyev('V','U',mbc%thsz,cnw,mbc%thsz,w,work,lwork,info)
                  if (info /= 0) stop 'dsyev i info /= 0'    
              endif
              
              
              lwork = int(work(1))
              deallocate(work)
              allocate(work(lwork))
              
              cnw(:,:) = mbc%h(:,:,isp)
              
              if (mbc%ovl) then
                  liwork = int(iwork(1))
                  deallocate(iwork)
                  allocate(iwork(liwork))
                  ovlm(:,:) = mbc%s(:,:)
                  call dsygvd(1,'V','U',mbc%thsz,cnw,mbc%thsz,ovlm,mbc%thsz,w,work,lwork,iwork,liwork,info)  
                  if (info /= 0) stop 'dsygvd r info /= 0'
                  deallocate(iwork)
              else
                  call dsyev('V','U',mbc%thsz,cnw,mbc%thsz,w,work,lwork,info)
                  if (info /= 0) stop 'dsyev i info /= 0' 
              endif
              
              mbc%cn(:,:,isp) = cmplx(cnw(:,:),0.0_dp,kind=dp)
              mbc%en(:,isp) = w(:)
              deallocate(work)
              deallocate(cnw,w)
              if (mbc%ovl) deallocate (ovlm)
              
!               if (.not. quiet) then
!                 print *, 'eigvals',isp
!                 print *, mbc%en(:,isp)
!               endif
              call system_clock(c4)
              if (.not. quiet) print *, 't_tben', real(c4-c3,8)/cr, isp
          endif

          m = mfac

          w0    = kt*(2*m)
          phase = pi/real(2*m, dp)
          zp    = cmplx(cos(phase),sin(phase), kind=dp)
!          zfac  = zp*zp

          call system_clock(c3)

          do p = mpmap(iproc), mpmap(iproc+1)-1
              zpn = zp**(2*p+1)
              zep  = efermi+w0*(zpn-1)
              zepb  = efermb+w0*(zpn-1)

!               call a_and_ac(mbc, zep, isp,.false.)
              call gg_diag(mbc,zep,isp)
            
              do ia = 1, mbc % n(1)
                  call assoc_ab(ia,isp)
                  call assoc_ham(ia) ! for cluster, pcluster and decipher

                  gia = mbc % rlst(1) -1 + ia
                  
! Calculate onsite energy and charges for promotion energy
                  if (.not. mbc%mxgf) then
                      call get_gab(zepb, ia, ia, 1, gab,.true.)
                      do i = 1, mbc%pah(ia+1)-mbc%pah(ia) 
                          bsons = bsons +real(zepb * zpn * gab(i,i))
                          q(ia) = q(ia) + real(zpn * gab(i,i))
                      end do
                  else
                      do i = mbc%pah(ia)+1, mbc%pah(ia+1)
                          bsons = bsons +real(zep*zpn * mbc % ggd(i))
                          q(ia) = q(ia) + real(zpn * mbc % ggd(i))
                      end do
                  endif
!          BOP intersites

                  ja = mbc % apb(ia)
                  ib = mbc % bpb(ja)

                  do while (ib /= eol)
                      jb = decipher(ib)

                      if (jb /= 1) then !jb == 1 => ib == ia THIS HAS BEEN CHANGED TO JB /=1
                          gib = mbc%rlst(1)-1 + ib !

                          gab = 0.0_dp

                          call get_gab(zepb, ia, ib, jb, gab,.false.)
                          
                          if (map(ib) == ib .and. mbc%mxgf) then
                              do i = 1, mbc % pbh(ia+1) - mbc % pbh(ia)
                                  do j = 1, mbc % pbh(ib+1) - mbc % pbh(ib)
                                      gab(i,j) = gab(i,j) + &
                                      &    sum(mbc % ac(mbc % pbh(ia) + i, :) * (mbc % a(mbc % pbh(ib) + j, :)))
                                      gab(i,j) = 4*kt*real(zpn*gab(i,j))
                                      ebond = ebond + gab(i,j)*darec(0,i,j,jb)
                                      
                                      
                                  end do
                              end do
                          else
                              nsttb = nstt(z(ib))
                              do i = 1, mbc % pbh(ia+1) - mbc % pbh(ia)
                                  do j = 1, nsttb
                                      gab(i,j) = 4*kt*real(zpn*gab(i,j))
                                      ebond = ebond + gab(i,j)*darec(0,i,j,jb)
                                      
!                                       ebondb = ebondb  + gab(i,j)*darec(0,i,j,jb)
                                  end do
                              end do
                          endif
                      end if
                      ja = ja + 1
                      ib = mbc % bpb(ja)
                  end do

              end do
          

  !          TB all
              do it = 1, mbc % n(2)
                  git = mbc % rlst(2) -1 + it
                  jt = mbc % apt(it)+1

                  ib = mbc % bpt(jt)
                
                  do while (ib /= eol)
                      gib = mbc % rlst(2) -1 + ib

                      gab = 0.0_dp
                      gab(1 : mbc % pth (it + 1) - mbc % pth(it), 1 : mbc % pth (ib + 1) - mbc % pth(ib)) =  &
                          & mbc % gt(mbc % pth(it) + 1 : mbc % pth(it + 1), mbc % pth(ib) + 1 : mbc % pth(ib + 1))
                          
                          
                      if (.not. mbc % ovl .or. (mbc%orth(it) + mbc%orth(ib) /= 0)) then
                          do i = 1, mbc % pth(it+1) - mbc % pth(it)
                              do j = 1, mbc % pth(ib+1) - mbc % pth(ib)
                                    gab(i,j) = 4*kt*real(zpn*gab(i,j))   
                                    ebond = ebond + gab(i,j)*(mbc%h(mbc%pth(it)+i,mbc%pth(ib)+j,isp))
!                                     ebondt = ebondt  + gab(i,j)*darec(0,i,j,jb)
                              end do
                          end do
                      else
            !       In the N.O. case this becomes the covalent energy instead of the bond energy, see Finnis textbook pg 204      
                          do i = 1, mbc % pth(it+1) - mbc % pth(it)
                              do j = 1, mbc % pth(ib+1) - mbc % pth(ib)                                        
                                  gab(i,j) = 4*kt*real(zpn*gab(i,j))                             
                                  ebond = ebond + gab(i,j)*(mbc%h(mbc%pth(it)+i,mbc%pth(ib)+j,isp) - &
                                          & mbc %s(mbc%pth(it)+i,mbc%pth(ib)+j)*(mbc%h(mbc%pth(it)+i,mbc%pth(it)+i,isp)))             
                                          
        
 
                                          
                              end do
                          end do
                      endif
                        
                      jt = jt + 1
                      ib = mbc % bpt(jt)
                  end do
            !             endif

            !          TB->BOP intersites
            !             if (mbc%mxgf .or. .not. mbc%allbop) then
                  if (.not. mbc%mxgf) then    
                      call assoc_ab(git,isp)
                      call assoc_ham(git) ! for cluster, pcluster and decipher
                  endif

                  jt = mbc % apd(it)
                  ib = mbc % bpd(jt)

                  if (mbc% bhsz /=0) then
                      do while (ib /= eol)
                          gib = mbc%rlst(1)-1 + ib
                          gab = 0.0_dp
                          if (mbc%mxgf) then
                              gab(1 : mbc % pth (it + 1) - mbc % pth(it), 1 : mbc % pbh (ib + 1) - mbc % pbh(ib)) = &
                              &  -transpose((mbc%ac(mbc%pbh(ib)+1:mbc%pbh(ib+1),mbc%pth(it)+1 :mbc%pth(it+1))))
                              
                              nsttb = mbc % pbh(ib+1) - mbc % pbh(ib)
                          else 
                              jb = decipher(ib)
                              call get_gab(zep, git, ib, jb, gab,.false.)
                              
                              nsttb = nstt(z(ib))
                          endif
!                 & - conjg(transpose(mbc % ac(mbc % pbh(ib) + 1 : mbc % pbh(ib + 1), mbc % pth(it) + 1 : mbc % pth(it + 1))))

                        
                          if (map(ib) == ib) gab = gab * 2
                          
                          do i = 1, mbc % pth(it+1) - mbc % pth(it)
                              do j = 1, nsttb
                                  gab(i,j) = 4*kt*real(zpn*gab(i,j))
                                  ebond = ebond + gab(i,j)*mbc%dh(i,j, jt)
!                                   ebondd = ebondd  + gab(i,j)*darec(0,i,j,jb)
                              end do
                          end do
                        
                          jt = jt + 1
                          ib = mbc % bpd(jt)
                      end do
                  endif
              end do

!    Multiply density matrix by overlap for onsite energy and charges
              if (mbc%ovl .and. mbc % thsz > 0) then
                  z_one = 1
                  z_zero = 0
                  gti(:,:) = mbc % gt(:,:)
                  si(:,:) = mbc % s(:,:)
                  call zgemm ('N', 'N', mbc % thsz, mbc % thsz, mbc % thsz, z_one, &
                              & si, mbc % thsz, gti, mbc % thsz, &
                              & z_zero, mbc%gt, mbc % thsz)
              endif    
!   Calculate onsite energy and charges for promotion energy
              do it = 1, mbc%n(2)
                  git = mbc % rlst(2) -1 + it
                  i = mbc%pth(it)
                  zt = z(git)
                  if (mbc%orth(it) == 1) then
                      llistt(:) = llist(:,zt)
                      nlt = nl(zt)
                  else
                      llistt(:) = mbc % tbc % llist(:,zt)
                      nlt = mbc % tbc % nl(zt)
                  endif
                  do no = 1,nlt
                      lla = llistt(no)
                      do i = i+1,i+2*lla+1

                        bsons = bsons +real(zep*zpn * mbc % gt(i,i))
                        ql(lla+1,it) = ql(lla+1,it) + real(zpn * mbc % gt(i,i))
                      enddo 
                      i=i-1
                  enddo
              
              enddo
          end do

          call system_clock(c4)
          if (.not. quiet) print *, 't_enpar', real(c4-c3,8)/cr, isp

          call system_clock(c3)
          q = 4 * kt * q
          ql(:,:) = 4 * kt * ql(:,:) 
          
          
          call mbprom(mbc,q,ql,isp,eprom,eprom0)

          call system_clock(c4)
          if (.not. quiet) print *, 't_mbprom', real(c4-c3,8)/cr, isp
          
      enddo

      bsons = 4*kt*bsons
!       call mbosprom(mbc,efermi,efermb,bsons,eprom)
      call system_clock(c3)
      uent  = mbent(mbc,efermi,efermb)
      call system_clock(c4)
      if (.not. quiet) print *, 't_mbent', real(c4-c3,8)/cr
 
      call system_clock(c3)
      buf = [ebond, eprom, bsons,eprom0]
      call mpi_allreduce(mpi_in_place, buf, 4, mpi_real8, mpi_sum, mpi_comm_world, ierr)
      
      ebond = buf(1)
      eprom = buf(2)
      bsons = buf(3)
      eprom0 = buf(4)

      call system_clock(c4)
      if (.not. quiet) print *, 't_allred', real(c4-c3,8)/cr
      
!      if (nproc > 1) then
!!           call mpi_allreduce(mpi_in_place,  eprom0, 1, mpi_real8,mpi_sum, mpi_comm_world, ierr)
!
!          if (iproc /= master) then
!             if (mpmap(iproc) /= mpmap(iproc+1)) then
!!                 buf(:,1) = (/ dndm, eprom, ebond, uent, bsons, emag /)
!                 buf(:,1) = (/ ebond, eprom, bsons,eprom0 /)
!                call mpi_irsend(buf, n, mpi_real8, master, iproc, mpi_comm_world, req(1), ierr )
!                call mpi_wait(req(1),stat(:,1),ierr)
!             end if
!
!             call mpi_bcast(buf,n, mpi_real8, master, mpi_comm_world, ierr )
!
!             
!             ebond = buf(1,1)
!             eprom = buf(2,1)
!             bsons = buf(3,1)
!             eprom0  = buf(4,1)
!             
!!              emag  = buf(5,1)
!             
!             
!             
!             deallocate(buf,req,stat)
!
!          else if (aproc>1) then
!!               ebsos = 0.0_dp
!    !          tm1 = mpi_wtime()
!
!             call mpi_waitall(aproc-1,req,stat,ierr)
!             do i=1,aproc-1
!                
!                ebond = ebond + buf(1,i)
!                eprom = eprom + buf(2,i)
!                bsons = bsons + buf(3,i)
!                eprom0  = eprom0  + buf(4,i)
!                
!!                 emag  = emag  + buf(5,i)
!             end do
!
!!              buf(:,1) = [ebond, eprom, uent, bsons, emag]
!            buf(:,1) = [ebond, eprom, bsons,eprom0]
!!             buf(1,1) = ebond
!             call mpi_bcast(buf(1,1), n, mpi_real8, master, mpi_comm_world, ierr )
!
!
!             deallocate(buf,req,stat)
!          end if
!
!
!      endif
      
      if (mbc% ovl) deallocate(gti,si)
      deallocate(q,ql)
      
      if (mag) then
          eprom = 0.5_dp*eprom
          ebond = 0.5_dp*ebond
          uent  = 0.5_dp*uent
          bsons = 0.5_dp*bsons
          emag  = emag + mbmag(mbc)
          eprom0 = 0.5_dp*eprom0
      end if

      eband = eprom + ebond
      error = abs((bsons-eband)/(bsons+eband))
      
      if ((.not. quiet)) then
      
          print *, 'bsons: ', bsons
          print *, 'eband: ', eband
          print *, 'ebond: ', ebond
          print *, 'eprom: ', eprom
          print *, 'uent: ', uent
          if (mag) print *, 'emag: ', emag
          print *, 'eprom0: ', eprom0
           
          print *, 'ecorr:',ecorr
          print *, 'bserr:', error
          
      endif
      
      ebond = ebond + ecorr

      eprom = eprom0

      if (iproc == master) then
          if (error > 10.0_dp*maxerr) then
              write(6,'("Unacceptably large uncertainty in the band structure energy.")')
              write(6,'("Uncertainty = ",f12.2,"%,"," on/inter",2(x,f20.12))') error*100.0_dp, bsons, eband
              write(6,'(''==> Terminating calculation.'')')
              call panic()
          elseif (error > maxerr) then
              write(6,'("WARNING: Large uncertainty in the band structure energy.")')
              write(6,'("Uncertainty = ",f12.2,"%,"," on/inter",2(x,f20.12))') error*100.0_dp, bsons, eband
          endif
      endif
   call system_clock(c2)

   if (.not. quiet) print *, 't_mbenergy', real(c2-c1,8)/real(cr,8)
      
   end subroutine mbenergy   
   
   subroutine mbprom(mbc,q,ql,isp,mbeprom,mbeprom0) 
      use mod_all_scalar, only: mag
!       use mod_atom_ar
      use mod_precision          
!       use mod_const
!       use mod_srt
!       use ab_io
!       use mod_ham

      include "../Include/Atom.array"
      
      type(mbtb_t), intent(in) :: mbc   
      integer :: isp, la, ia,tia
      real(dp) :: ela, deia
      real(dp), intent(inout) :: mbeprom

      integer lla
      integer nstta,nla
      integer za,n
      real(dp) :: ql(3,mbc%n(2)),q(mbc%n(1))
       
      integer :: llistt(3),nlt
       
      real(dp), intent(inout) :: mbeprom0
      real(dp) :: ela0

!       mbeprom = 0.0_dp

      do ia = 1, mbc % n(1)
          za = z(ia)
          deia = de(ia)
          
          if (mag) deia = deia + oppm(isp)*dem(ia)
          
          call states(za,nla,nstta,llista)
          do la = 1,nla
            lla = llista(la)
            if (lla == 0) then
                ela = es(za) + deia
                ela0 = es(za)
            elseif (lla == 1) then
                ela = ep(za) + deia
                ela0 = ep(za)
            elseif (lla == 2) then
                ela = ed(za) + deia
                ela0 = ed(za)
            endif
            mbeprom = mbeprom + q(ia)*ela
            mbeprom0 = mbeprom0 + q(ia)*ela0
          enddo    
      enddo
                
      do ia = mbc%n(1)+1, mbc % nd
          za = z(ia)
          
          tia = ia - mbc%n(1)
          if (mbc%orth(tia) == 1) then
              llistt(:) = llist(:,za)
              nlt = nl(za)
!               print *, iproc,'orth',ql(:,Tia)
          else
              llistt(:) = mbc % tbc % llist(:,za)
              nlt = mbc % tbc % nl(za)
!               print *,iproc, 'nonorth',ql(:,Tia)
          endif            
          n = 1
          do la = 1,nlt
              lla = llistt(la)
              
              if (lla == 0) then
                  ela = (mbc%h(mbc%pth(tia)+1,mbc%pth(tia)+1,isp))*ql(1,tia)
                  ela0 = (mbc%h0(mbc%pth(tia)+1,mbc%pth(tia)+1))*ql(1,tia)

                  n = 2
              elseif (lla == 1) then
                  ela = (mbc%h(mbc%pth(tia)+n,mbc%pth(tia)+n,isp))*ql(2,tia)
                  ela0 = (mbc%h0(mbc%pth(tia)+n,mbc%pth(tia)+n))*ql(2,tia)
                  n = n + 3
              elseif (lla == 2) then
                  ela =  (mbc%h(mbc%pth(tia)+n,mbc%pth(tia)+n,isp))*ql(3,tia)
                  ela0 =  (mbc%h0(mbc%pth(tia)+n,mbc%pth(tia)+n))*ql(3,tia)
             
              endif
              mbeprom = mbeprom + ela
              mbeprom0 = mbeprom0 + ela0

          enddo             
      enddo
          
   end subroutine mbprom

   
!    function mbeprom(mbc,ef) result(mbprom)
!           use mod_precision
!           use mod_all_scalar    
!           use mod_const
!           use mod_funptr
!           use ab_io, only : diag, eigvec, lchain
! 
! !
! !    This is a subroutine to evaluate the promotion energy.
! 
!       include "../Include/Atom.array"
!       
!       type(mbtb_t), intent(inout) :: mbc  
!       integer :: ia,isp
!       real(dp), intent(in) :: ef
!       
!       real(dp) :: mbprom
!       real(dp) :: qla(mbc%thsz)
!       real(dp) :: ela, deia
!       complex(dp) :: zp,zfac,zep
!       real(dp) :: phase,w0      
!       integer :: m,p
!       integer lla
!       integer nstta,nla,la
!       integer za,st,n
! 
!          mbc%q = 0
!          m = mfac
! 
!       print *, 'm:',m
!       mbprom = 0.0_dp
!       do isp = 1, nsp         
!          mbc % q = 0.0_dp
!          mbc % ggd = 0.0_dp
! 
!          w0    = kt*(2*m)
!          phase = pi/real(2*m, dp)
!          zp    = cmplx(cos(phase),sin(phase), kind=dp)  
!          zfac  = zp*zp
! 
!          do p = 0, m-1
!             zep  = ef+w0*(zp-1)
! !             print *, p
!             call gg_diag(mbc, zep, isp)
! 
!             do ia = 1, mbc % nd
!                do la = mbc%pah(ia)+1, mbc%pah(ia+1)
!                   mbc % q(ia) = mbc % q(ia) + real(zp * mbc % ggd(la))
! !                   if (ia > mbc%n(1)) qla((mbc%pah(ia)-mbc%bhsz)+la) = qla((mbc%pah(ia)-mbc%bhsz)+la) + real(zp * mbc % ggd(la))
!                end do
!             end do      
! 
!             zp   = zfac*zp
!          enddo
! 
!          mbc % q = 4 * kt * mbc % q
!          
! 
!         
!         
!         do ia = 1, mbc % nd
!         
!         
! !       do ia = 1,mbc%n(1)
! !          if (ia <= mbc%n(1)) then
! 
! 	    za = z(ia)
! 	    
! 	    deia = de(ia)
! 	    if (mag) deia = deia + oppm(isp)*dem(ia)
! 	    
! 	    call states(za,nla,nstta,llista)
! 	    
! 	    do la = 1,nla
! 	      lla = llista(la)
! 	      if (lla == 0) then
! 		  ela = es(za) + deia
! 	      elseif (lla == 1) then
! 		  ela = ep(za) + deia
! 	      elseif (lla == 2) then
! 		  ela = ed(za) + deia
! 	      endif
!               mbprom = mbprom + mbc%q(ia)*ela
!             enddo    
!            enddo
!          enddo
!       end function mbeprom

      function mbent(mbc,ef,efb) result(entropy)
        use mod_all_scalar, only : kt, mfac, nsp,nrec, quiet
        use mod_atom_ar
        use mod_precision          
        use mod_const
        use mod_srt
        use ab_io
        use mod_ham
        use mod_clock
        use topologia, only : iproc,nproc, master, sdistrib

      integer(8) :: c1, c2
      integer :: ierr, aproc
      integer, allocatable :: szs(:), displs(:)
      type(mbtb_t), intent(inout) :: mbc   
      real(dp), intent(in) :: ef,efb
      integer :: isp, la, ia,m,i,it
      integer, parameter :: mx = 1111
      real(dp) :: entropy,simp,x0,x1,x,dx,entden
!      real(dp) :: y(mbc%nd,mx),yn(mx)
      real(dp), allocatable :: y(:,:),yn(:)
      complex(dp) :: zi,zib,gab(mxnstat,mxnstat)
      
      complex(dp) :: z_one,z_zero
      complex(dp),allocatable :: si(:,:),gti(:,:)
      

      call system_clock(c1)
      
      if (mbc%ovl) then
          if (.not. allocated(si)) allocate(si(mbc%thsz,mbc%thsz))
          if (.not. allocated(gti)) allocate(gti(mbc%thsz,mbc%thsz))
      endif

      m = min(mx,max(111,nrec*nint(500*kt)))
      m = 2*(m/2)+1

      allocate(displs(0:nproc))
      call sdistrib(m, nproc, displs, aproc)

      if (iproc == master) then
          allocate(y(mbc%nd,0:m-1))
      else
          allocate(y(mbc%nd,displs(iproc):displs(iproc+1)-1))
      end if
      y=0.0_dp

      x0=-10_dp
      x1=10_dp
      dx = (x1-x0)/real(m-1, dp)

      do isp = 1, nsp         
         

!          w0    = kt*(2*m)
!          phase = pi/real(2*m, dp)
!          zp    = cmplx(cos(phase),sin(phase), kind=dp)
!          zfac  = zp*zp
          
         do i = displs(iproc), displs(iproc+1)-1
            x = x0 + i*dx

            zi = cmplx(ef +kt*x,kt, kind=dp)           
            zib = cmplx(efb +kt*x,kt, kind=dp)

            mbc % ggd = 0.0_dp
            call gg_diag(mbc, zi, isp)
            
            if (mbc%ovl .and. mbc % thsz > 0) then
                  z_one = 1
                  z_zero = 0
                  gti(:,:) = mbc % gt(:,:)
                  si(:,:) = mbc % s(:,:)
                  call zgemm ('N', 'N', mbc % thsz, mbc % thsz, mbc % thsz, z_one, &
                            & si, mbc % thsz, gti, mbc % thsz, &
                              & z_zero, mbc%gt, mbc % thsz)
            endif    
            
            

            do ia = 1, mbc % nd
               if (.not. mbc%mxgf .and. ia <= mbc%n(1)) then
                  call assoc_ab(ia,isp)
                  call assoc_ham(ia) ! for cluster, pcluster and decipher
                  call get_gab(zib, ia, ia, 1, gab,.false.)
                  do la = 1, mbc%pah(ia+1) - mbc%pah(ia)
                    y(ia,i) = y(ia,i) + aimag(gab(la,la))*entden(kt*x,kt)
                  end do
               elseif (mbc%mxgf .and. ia <= mbc%n(1)) then
                  do la = mbc%pah(ia)+1, mbc%pah(ia+1)
                      y(ia,i) = y(ia,i) + aimag(mbc % ggd(la))*entden(kt*x,kt)
                  end do
               else
                  it  = ia - mbc % n(1)
                  do la = mbc%pth(it)+1, mbc%pth(it+1)
                      y(ia,i) = y(ia,i) + aimag(mbc % gt(la,la))*entden(kt*x,kt)
                  end do 
               endif
            end do
         enddo
      end do

      call system_clock(c2)       
      if (.not. quiet) print *, 't_entr1', real(c2-c1,8)/real(cr,8)

      call system_clock(c1)
      allocate(szs(0:nproc-1))
      displs = displs * mbc % nd
      szs = displs(1:nproc) - displs(0:nproc-1)

      if (iproc == master) then
          call mpi_gatherv(mpi_in_place, 0, mpi_datatype_null, y, szs, &
                 displs, mpi_real8, master, mpi_comm_world, ierr)

          allocate(yn(mx))

          entropy = 0.0_dp
          do ia = 1, mbc % nd
              yn=y(ia,:)  
              entropy = entropy + simp(dx,yn,m)
          end do

          entropy = 2*kt*kt*entropy/pi

      else
          call mpi_gatherv(y, szs(iproc), mpi_real8, 0, 0, &
                 0, mpi_datatype_null, master, mpi_comm_world, ierr)
          deallocate(y,szs,displs)
      end if
      
      call system_clock(c2)
      if (.not. quiet) print *, 't_entrgthr', real(c2-c1,8)/real(cr,8)

      call system_clock(c1)
      call mpi_bcast(entropy, 1, mpi_real8, master, mpi_comm_world, ierr)
      call system_clock(c2)
      if (.not. quiet)  print *, 't_entbcast', real(c2-c1,8)/real(cr,8)

      if (iproc == master) deallocate(yn,y,szs,displs) 
!        


!       entropy = 0.0_dp
!       y=0.0_dp
!       do isp = 1, nsp   
! 
!                m     = mfac
!               m     = nint(mfac*0.25_dp*1.8_dp/kt)
!                w0    = kt*(2*m)
!                phase = pi/real(2*m, dp)
!                zp    = cmplx(cos(phase),sin(phase), kind=dp)
!                zfac  = zp*zp
!                
!                do p = 0, m-1
!                   zep = ef + w0*(zp-1)
!                   
!                   call gg_diag(mbc, zep, isp)
!                   
!                   do la = 1, mbc % ahsz
!                     entropy = entropy &
!                                 & + real(zp*zentden(zep)*mbc%ggd(la))
!                   enddo
!                   zp  = zp*zfac                                    
!                end do
!                
!                
! !             
!                enddo
!               entropy = -4*kt*entropy
!                 
!              contains
!       
!          function zentden(x)
!             use mod_precision
! 
! !             evluate complex entropy density function.
! 
!             implicit none
!             complex(dp), intent(in) :: x
!             complex(dp) :: zentden
!             complex(dp) :: n,m
! 
!             n = 1.0_dp/(1.0_dp+exp(x))
!             m = 1.0_dp - n
! 
!             zentden = n*log(n) + m*log(m)
!          end function zentden
      end function mbent
      
      
      subroutine mbpair(mbc)
         use mod_conf
         use mod_all_scalar 
         use mod_atom_ar, only : asort, bsort, btype
         use mod_pft, only : tailed_fun
      
      type(mbtb_t), intent(inout) :: mbc
      type(ham_conf_t),pointer :: hamtb,hambp
      real(dp) :: rcutt,ett(3,3)
      logical :: stressed
      integer :: rlimt(3),maxneigh,df,bt,nat
      integer :: ia,ib,it,gia,git,gib,ja,jt,i
      type(pwp_t), pointer :: pwp 
      real(dp) :: dr(3),r,fv(0:1)
      integer, allocatable :: aptl(:),bptl(:)

      real(dp) :: epaira
      
      include "../Include/PosVel.array"
      include "../Include/NebList.array"    
      include "../Include/Atom.array"
      include "../Include/ag.conn"
      include "../Include/NebListL.array"
      include "../Include/Force.array"
      

      allocate(aptl(mxtotnd),bptl(mxnbl))
      aptl = 0
      bptl = 0
      
      
      ett = 0.0_dp
      stressed = totstr /= 0.0_dp .and. any(et /= 0.0_dp)      
      if (stressed) ett = et*totstr
      
      if (stressed) then
          ad(:,:nd) = unshad(:,:nd)
          adinert(:,:ninert) = unshind(:,:ninert)
      end if
         
         
      if (rtc % rcut(1) /= rtc % rcut(2)) then
          mbc % apb = 0
          mbc % bpb = 0
          mbc % apd = 0
          mbc % bpd = 0
          call bldnebt(.true., rlim, lena, a, rprune, rcutl, stressed, ett, &
                      &  nd, ad, z, ninert, adinert, zinert, &
                      &  map, aptl, bptl,  maxneigh, totnd, nne, mapi,mxnnb)
               nat = mbc % n(1)
               if (sum(rlim) /= 0) nat = mxnbl
          call extract_boplist(aptl, bptl,mbc % rlst, mbc % apb, mbc % bpb,nat,.true.)
          call extract_tblist(aptl, bptl, 2, 1, mbc % rlst, mbc % apd, mbc % bpd,mbc%mxgf)
      endif
        
      
        
      rcutt = rtc % rcut(2)    
      do i = 0, rtc % tbham_conf % nb -1
         if (rtc % tbham_conf % b(i) % pwp % tl % rc > rcutt) then
            rcutt = rtc % tbham_conf % b(i) % pwp % tl % rc
         endif
      enddo
      
      if ((rtc % rcut(1) == rtc % rcut(2)) .or. (rcutt /= rtc % rcut(2))) then      
          aptl = 0
          bptl = 0
          rlimt = 0
          call bldnebt(.false., rlimt, lena, a, rprune, rcutt, stressed, ett, &
                      &  nd, ad, z, ninert, adinert, zinert, &
                      &  map, aptl, bptl,  maxneigh, totnd, nne, mapi,mxnnb) 
      endif
       
      mbc % apt = 0
      mbc % bpt = 0
      call extract_tblist(aptl, bptl, 2, 2, mbc % rlst, mbc % apt, mbc % bpt,mbc%mxgf)
      

      deallocate(aptl,bptl)

      hambp => rtc % ham_conf 
      hamtb => rtc % tbham_conf 
      
       df = 0
      if (forces) df = 1


      do ia = 1, mbc % n(1)                                  
          ja = mbc % apb(ia)+1
          ib = mbc % bpb(ja)
          do while (ib /= eol)
              bt = bsort(btype(z(ia),z(ib)))
              if (bt == -1) bt = bsort(btype(z(ib),z(ia)))
              pwp => hambp % b(bt) % pwp 
              
              dr = ad(1:3,ib) - ad(1:3,ia)
              r = sqrt(sum(dr*dr))

              call tailed_fun(r, pwp % fp, pwp % tl, df, fv)
            
              epair = epair + fv(0)
                  
              if (forces) fpp(:,ia) = fpp(:,ia) + fv(1)*dr/r
              
!             if (ia > 16 .and. ib < 17) then
!                 print *, ia,ib,'force', fv(1)*dr/r
!             endif
            
              ja = ja + 1
              ib = mbc % bpb( ja )         
            end do
            

            
        end do



      do it = 1, mbc % n(2)
        
         epaira = 0.0_dp
 
         git = mbc % rlst(2) -1 + it

         jt = mbc % apt(it)+1
         ib = mbc % bpt(jt)

         do while (ib /= eol)
            gib = mbc % rlst(2) -1 + ib
            bt = bsort(btype(z(git),z(gib)))
            
            if (bt == -1) bt = bsort(btype(z(gib),z(git)))
!             pwp => hamtb % b(bt) % pwp 

            if (mbc%orth(it) + mbc%orth(ib) == 0) then
                pwp => hamtb % b(bt) % pwp 
            else
                pwp => hambp % b(bt) % pwp 
            endif 
            
            dr = ad(1:3,gib) - ad(1:3,git)
            
            r = sqrt(sum(dr*dr))

            call tailed_fun(r, pwp % fp, pwp % tl, df, fv)
           
            epair = epair + fv(0)
            epaira =  epaira + fv(0)                  

            if (forces) fpp(:,git) = fpp(:,git) + fv(1)*dr/r
            
!             if (git > 16 .and. gib < 17) then
!                 print *, git,gib,'force', fv(1)*dr/r
!             endif
            
            jt = jt + 1
            ib = mbc % bpt(jt)
            

         end do

  !       if (.not. quiet) print *, 'PP E:', it, epaira
      
         jt = mbc % apd(it)

         ib = mbc % bpd(jt)
       
         do while (ib /= eol)
         
            bt = bsort(btype(z(git),z(ib)))
            if (bt == -1) bt = bsort(btype(z(ib),z(git)))
            
!    dH atoms have BOP pwp as they'd have BOP interactions in dH, but
!     having TB gives results closer to pure BOP/TB, needs to be checked

!             pwp => hamtb % b(bt) % pwp
            pwp => hambp % b(bt) % pwp
            
            gib = map(ib)

            dr = ad(1:3,ib) - ad(1:3,git)

            r = sqrt(sum(dr*dr))

            call tailed_fun(r, pwp % fp, pwp % tl, df, fv)
           
            epair = epair + fv(0)
            if (gib == ib) epair = epair + fv(0)

            if (forces) then
                fpp(:,git) = fpp(:,git) + fv(1)*dr/r
!                 print *, git,gib,'force', fv(1)*dr/r
!                 pwp => hambp % b(bt) % pwp 
!                 call tailed_fun(r, pwp % fp, pwp % tl, df, fv)

                if (gib /= 0) fpp(:,gib) = fpp(:,gib) - fv(1)*dr/r
            endif
            
            jt = jt + 1
            ib = mbc % bpd(jt)
            
         end do


      end do
      
      epair = 0.5_dp * epair
      eclas = epair
      
!       do ia = 1, nd
!          print *, 'PAIR',ia,fpp(:,ia)
!       enddo
      
      end subroutine mbpair
   
  
      function mbmag(mbc) result(mbemag)
        
        use mod_precision
        use mod_const        
        use topologia
        use mod_atom_ar, only : mg
        
        implicit none
        
        type(mbtb_t), intent(inout) :: mbc   
        real(dp) :: mbemag
        integer :: ia
        
        include '../Include/Atom.array'
        
        
        mbemag = 0.0_dp
        
        do ia = 1,mbc % nd
            if (ia <= mbc%n(1)) then
                mbemag = mbemag + istn(z(ia))*(mg(ia)*mg(ia))
            else if (.not. mbc%sctb) then
                if (mbc%orth(ia - mbc%n(1)) == 0) then
                    mbemag = mbemag + mbc%tbc%istn(z(ia))*(mg(ia)*mg(ia))
                else 
                    mbemag = mbemag + istn(z(ia))*(mg(ia)*mg(ia))
                endif
            endif
        end do
        
        mbemag = 0.25_dp*mbemag

    end function mbmag
   
   
   
   
   
!    subroutine gtb(mbc, ez, isp)
!       use mod_all_scalar, only : mag
!
!
!       type(mbtb_t), intent(inout) :: mbc
!       complex(dp), intent(in) :: ez
!       integer, intent(in) :: isp
!
!       integer :: nstt1, nstt2, nsttb, l1, l2, lb, jt, it1, it2, ib, pt1, pt2, pb, it, git, lt
!       integer :: i, off
!
!       real(dp) :: dea
!       complex, parameter :: z_one = 1.0_dp, z_zero = 0.0_dp
!
!
!       integer :: info, lwork
!       real(dp), allocatable :: rwork(:)
!       complex(dp), allocatable :: work(:)
!
!       include '../Include/Atom.array'
!
!       ! (z*I*S - (H + dH'*Gb*dH))^-1 = W (z*I - V)^-1 W'
!
!       mbc % w = mbc % h
! !       mbc % w = 0.0_dp
!
!       off = mbc % rlst(2)-1
!
!       do it = 1, mbc % n(2)
!          git = off + it
!
!          dea = de(git)
!          if (mag) dea = dea + oppm(isp)*dem(git)
!
!          do lt = mbc % pth(it)+1, mbc % pth(it+1)
!             mbc % w(lt,lt) = mbc % w(lt,lt) + dea
!          end do
!       end do
!
!
!       do it1 = 1, mbc % n(2)
!          pt1 = mbc % pth(it1)
!          nstt1 = mbc % pth(it1+1) - pt1
!
!          do it2 = 1, mbc % n(2)
!             pt2 = mbc % pth(it2)
!             nstt2 = mbc % pth(it2+1) - pt2
!
!             jt = mbc % apd(it1)
!             ib = mbc % bpd(jt)
!             do while ( ib /= eol)
!                pb = mbc % pbh(ib)
!                nsttb = mbc % pbh(ib+1) - pb
!
!                do l1 = 1, nstt1
!                   do l2 = 1, nstt2
!                      do lb = 1, nsttb
!                         mbc % w(pt1+l1,pt2+l2) =  mbc % w(pt1+l1,pt2+l2) &
!                            & + mbc % dh(l1, lb, jt) * mbc % a (pb+lb,pt2+l2)
!                      end do
!                   end do
!                end do
!
!                jt = jt + 1
!                ib = mbc % bpd(jt)
!             end do
!
!          end do
!       end do
!
!       call print_c(mbc%w,500,'sgm'//achar(10))
!
! !       stopgg_diag
!       allocate(work(1), rwork(3*mbc % thsz-2))
!       lwork = -1
!
!       if (.not. mbc%ovl) then
!          call zheev('v', 'l', mbc % thsz, mbc % w, mbc % thsz, mbc % v, work, lwork, rwork, info)
!       else
!          call zhegv(1, 'v', 'l', mbc % thsz, mbc % w, mbc % thsz, mbc % s, mbc % thsz, &
!                                                    & mbc % v, work, lwork, rwork, info)
!       end if
!
!       lwork = work(1)
!       deallocate(work)
!       allocate(work(lwork))
!
!       if (.not. mbc%ovl) then
!          call zheev('v', 'l', mbc % thsz, mbc % w, mbc % thsz, mbc % v, work, lwork, rwork, info)
!       else
!          call zhegv(1, 'v', 'l', mbc % thsz, mbc % w, mbc % thsz, mbc % s, mbc % thsz, &
!                                                    & mbc % v, work, lwork, rwork, info)
!       end if
!
!       deallocate(work,rwork)
!
! !       print *, info 
!
!       call print_c(mbc%v,400,'v'//achar(10))
!
!       do i = 1, mbc % thsz
!          mbc % w (:,i) = mbc % w(:,i)/sqrt(ez - mbc % v(i))
!       end do
!
!       call print_c(mbc%w,400,'w'//achar(10))
!
!
!
!       call zgemm ('n', 'c', mbc % thsz, mbc % thsz, mbc % thsz, z_one, &
!                            & mbc % w, mbc % thsz, mbc % w, mbc % thsz, &
!                            & z_zero, mbc % gt, mbc % thsz)
!
!       stop 'remove me when you find out whats wrong with zgemm'
!
!       call print_c(mbc%gt,400,'gt'//achar(10))
!
!
! !       AC, C: gt. may be try zhemm?
!       call zgemm ('n', 'n', mbc % bhsz, mbc % thsz, mbc % thsz, z_one, &
!                            & mbc % a, mbc % bhsz, mbc % gt, mbc % thsz, &
!                            & z_zero, mbc % ac, mbc % bhsz)
!
!    end subroutine gtb
!gg_diag

   subroutine gtdhgb(mbc, ez, isp)
      use mod_all_scalar, only : mag
      use ab_io
      use mod_ham

      type(mbtb_t), intent(inout) :: mbc
      complex(dp), intent(in) :: ez
      integer, intent(in) :: isp

      integer :: nstt1, nstt2, nsttb, l1, l2, lb, jt, it1, it2, ib, pt1, pt2
      integer :: i, off, j,git2,jb

      real(dp) :: dea
      complex(dp) :: z_one  , z_zero  
!       ,si(mbc%thsz,mbc%thsz)
      complex(dp) :: gab(mxnstat, mxnstat)
      complex(dp),allocatable :: cm(:,:),cz(:,:), gti(:,:)

!       integer :: info, lwork
!       real(dp), allocatable :: rwork(:)
!       complex(dp), allocatable :: work(:)
!       integer, allocatable :: ipiv(:)
      

      include '../Include/Atom.array' ! for nstt & z only
      

      ! Gt = gt + gt * dh12 * G21
      if (.not. allocated(cm)) allocate(cm(mbc%thsz,mbc%thsz))
      if (.not. allocated(cz)) allocate(cz(mbc%thsz,mbc%thsz))
      if (.not. allocated(gti)) allocate(gti(mbc%thsz,mbc%thsz))

      gti = 0.0_dp

!       CALC gt AS (S*ez - H)^-1
!       mbc % gt = 0.0_dp
!       mbc % gt = mbc % h(:,:,isp)
! 
!        mbc % gt = - mbc % gt
! 
!       if (.not. mbc % ovl) then
!          do i = 1, mbc % thsz
!             mbc % gt(i,i) = mbc % gt(i,i) + ez
!          end do
!       else
!          mbc % gt = mbc % gt + (mbc % s * ez)
!       end if
! 
! 
!             
!       if (mbc % thsz /= 0) then
!         allocate(ipiv(mbc % thsz))
! 
!         call zgetrf (mbc % thsz, mbc % thsz, mbc % gt, mbc % thsz, ipiv, info)
!         if (info /= 0) stop 'zgetrf info /= 0'
! 
!         allocate(work(1))
!         lwork = -1
! 
!         call zgetri(mbc % thsz, mbc % gt, mbc % thsz, ipiv, work, lwork, info)
!         if (info /= 0) stop 'zgetri i info /= 0'
! 
!         lwork = work(1)
!         deallocate(work)
!         allocate(work(lwork))
! 
!         call zgetri(mbc % thsz, mbc % gt, mbc % thsz, ipiv, work, lwork, info)
!         if (info /= 0) stop 'zgetri r info /= 0'
! 
!         deallocate(work,ipiv)
!       endif

! CALC gt AS Cn((ez-en)^-1)Cn*

      do i = 1, mbc%thsz
          cm(:,i) = mbc%cn(:,i,isp) * (1.0_dp/(ez-mbc%en(i,isp)))
      enddo
      
      cz(:,:) = mbc%cn(:,:,isp)
      
       z_one= 1.0_dp
       z_zero =0.0_dp

       call zgemm ('N', 'T', mbc % thsz, mbc % thsz, mbc % thsz, z_one, &
                      & cm, mbc % thsz, cz, mbc % thsz, &
                      & z_zero, mbc % gt, mbc % thsz)
      
      
      do it1 = 1, mbc % n(2)
         pt1 = mbc % pth(it1)
         nstt1 = mbc % pth(it1+1) - pt1

         do it2 = 1, mbc % n(2)
            git2 = it2 + mbc % n(1)
            pt2 = mbc % pth(it2)
            nstt2 = mbc % pth(it2+1) - pt2
            
            call assoc_ab(git2,isp)
            call assoc_ham(git2)
            
            jt = mbc % apd(it1)
            ib = mbc % bpd(jt)
            do while ( ib /= eol)
            
                
               
               jb = decipher(ib)
               
!                print *, 'it1,it2,ib,decipher'
!                print *, it1,it2,ib,decipher(ib)
               if (jb /= 0) then
                  nsttb = nstt(z(ib))
                  call get_gab(ez, git2,ib, jb, gab,.false.)
                  

                  
                  do l1 = 1, nstt1
                      do l2 = 1, nstt2
                          do lb = 1, nsttb
                              gti(pt1+l1,pt2+l2) = gti(pt1+l1,pt2+l2)+mbc % dh(l1, lb, jt) * gab(l2,lb)
                          end do
                      end do
                  end do
               endif

               jt = jt + 1
               ib = mbc % bpd(jt)
            end do

         end do
      end do
      

      cm(:,:) = mbc%gt(:,:)
      call zgemm ('N', 'N', mbc % thsz, mbc % thsz, mbc % thsz, z_one, &
                      & cm, mbc % thsz, gti, mbc % thsz, &
                      & z_one, mbc % gt, mbc % thsz)
      
      deallocate(cz,cm,gti)

!       if (mbc%sctb) then
!         do it1 = 1, mbc % n(2)          
!             do i = mbc%pth(it1)+1,mbc%pth(it1+1)
!                 mbc % dqchi(it1+mbc%n(1)) = mbc % dqchi(it1+mbc%n(1)) + &
!                         & real(zp * mbc % gt(i,i) * mbc % gt(i,i))
!             enddo
!         enddo
!       endif
     

!   This is to adjust the density matrix with the overlap for onsite elements required for charges, eprom, etc
!    see Finnis (2007) eqn. 21
!       if (mbc%ovl .and. rho .and. mbc % thsz > 0) then
!           gti(:,:) = mbc % gt(:,:)
!           si(:,:) = mbc % s(:,:)
!           call zgemm ('N', 'N', mbc % thsz, mbc % thsz, mbc % thsz, z_one, &
!                       & si, mbc % thsz, gti, mbc % thsz, &
!                       & z_zero, mbc%gt, mbc % thsz)
!       endif



   end subroutine gtdhgb
   
   
   subroutine bop_gf(mbc,ef,mxit,mgit,merr)
    use mod_all_scalar, only : mag, quiet, dqmax, forces, mgmax, nd, qerr, nsp,kt,mfac,lef
    use mod_atom_ar, only : dq, mg
    use mod_chi
    use mod_conf
    use mod_const
!     use ab_io
    use mod_ham
    use topologia, only : mpmap, iproc,nproc, master

    include "../Include/NebList.array"
    include "../Include/PosVel.array"
    include "../Include/Atom.array"
    include "../Include/Force.array" 



    type(mbtb_t), intent(inout) :: mbc
    integer :: it, mxit, ia,nat, mgit, mit
    real(dp) :: ef
    real(qp) ::  sumz ,merr,qerrb
    complex(dp) :: gab(mxnstat, mxnstat)
    integer :: ib, jb, mp,jt
    integer :: m, p,isp
    real(dp), allocatable :: mg_p(:),dec(:)
    character(len=25) :: ffmt
     real(dp) :: dqmax_old, aa, aamin, aamax

      m=mfac
      
      write(ffmt,"(a,i0,a)") '(a,":",', mbc%nd, '(x,f22.14))'

      sumz = sum( (/ (real(zc(z(ia)), qp), ia = 1,mbc%nd) /) )

      mbc%allbop = .true.
      
      allocate(mg_p(mbc%nd))
      if (.not. quiet) write(6,'(/,''Calculating full BOP GF'')')

      if (.not. mbc%sctb) then
          dqmax_old = 0.0_dp
          aamax=1.0_dp
          aamin=1.0e-4_dp
          aa=10.0_dp
          allocate(dec(nd))
          dec = real(de(:nd),qp)
      endif

      qerrb = qerr/100
      
      if (mag) mit = 1
      do 
          if (mag) then
              if (.not. quiet) write(6,'(/,"mit:",i3)') mit   
              
              
              mg_p = mg
              do ia = 1, mbc%nd
                  if (zc(z(ia)) /= 0) dem(ia) = 0.5_dp*istn(z(ia))*mg(ia)   
              enddo
          end if
          it = 1
          do           
              if (.not. quiet) write(6,'(/,"it:",i3)') it
              if (mbc%sctb) then
! ISSUE: how do i modify de() without every proc changing it?
                  de(:mbc%nd) = de(:mbc%nd) + mbc%dqchi(:mbc%nd) 
!                   print *, 'DQCHI',mbc%dqchi(:mbc%nd) 
!                   print *, iproc,'de:',de(:mbc%nd)
              else
                  de(:mbc%nd) = real(dec - sum( (/ (dec(ia)*zc(z(ia)), ia = 1,mbc%nd) /) ) / sumz, dp)
              endif
              
!              if (.not. quiet .and. mbc%nd < 500) then
!                  if (.not. quiet) write(6, ffmt) 'de', de(:mbc%nd)
!                  if (.not. quiet) write(6, ffmt) 'dem', dem(:mbc%nd)
!              endif 
!               call sdistrib(mbc%nd, nproc, mpmap, aproc)
              call recurse(.true.,mbc%nd) 
                  
              do ia =  1, mbc % nd
                  call assoc_ham(ia)
                  do jb = 1, pcluster(2)-1
                      decipher(cluster(jb)) = jb
                  end do
              end do

              
              
              ef = embeff(mbc,sumz,lef) 
              
              lef = ef
              
              dqmax = 0.0_dp
              
              do ia = 1, mbc % nd
                  dq(ia) = dq(ia) - zc(z(ia))
                  if (abs(dq(ia)) > abs(dqmax)) dqmax = dq(ia)
                  if (mbc%sctb) mbc % dqchi(ia) = dq(ia)/mbc % dqchi(ia)
!                   print *, 'bop',ia, dq(ia) + zc(z(ia)), dq(ia)
              enddo 
              if (.not. mbc%sctb) then
                  if (     (dqmax_old*dqmax > 0.0_dp &
                    &  .and. abs(dqmax_old) > abs(dqmax)) &
                  &  .or. (dqmax_old*dqmax < 0.0_dp)) then
                      aa=aa*abs(dqmax_old)/abs(dqmax_old-dqmax)
                  end if
                  aa=max(aamin,aa)
                  aa=min(aamax,aa)
                  if (abs(dqmax) <= qerrb) aa=aamin
                  dec(:mbc%nd) = de(:mbc%nd) + aa*dq(:mbc%nd)
                  dqmax_old = dqmax
              endif

              
              

!               if (iproc == master) print *, redb, 'bop dqmax:', dqmax, endc

!              if (.not. quiet) write(6, ffmt) 'dq', dq(:mbc%nd)
              if (iproc == master) print *, 'BOP dqmax:', dqmax
              
              if (abs(dqmax) <= qerrb .or. it == mxit) exit

              it = it + 1
          end do

        if (.not. mag) exit
        mgmax = maxval(abs(mg - mg_p))
        
!        if (.not. quiet) then
!            write(6, ffmt) 'mg_o', mg_p(:mbc%nd)
!            write(6, ffmt) 'mg', mg(:mbc%nd)
!        endif
        
        if (iproc == master) print *, 'BOP mgmax:',mgmax
!         
!         print *, iproc, 'mg',mg(:mbc%nd)
!         print *, iproc,'mg_o',mg_p(:mbc%nd)
!         do ia = 1, mbc%nd
!             print *, 'MG',mg(ia),'MGO',mg_p(ia)
!         enddo
    !          if (.not. quiet) then
    !             write(6,ffmt) 'mgo', mg_p
    !             write(6,ffmt) 'mg ', mg
    !             write(6,*) 'mgmax:', mgmax
    !          end if
        if (mgmax <= merr .or. mit == mgit) exit
        mit = mit + 1
        

      end do          

      deallocate(mg_p)
      
      mbc%allbop = .false.
   
   
   end subroutine bop_gf

   subroutine gtbi(mbc, ez, isp)
      use mod_all_scalar, only : mag

      type(mbtb_t), intent(inout) :: mbc
      complex(dp), intent(in) :: ez
      integer, intent(in) :: isp

      integer :: nstt1, nstt2, nsttb, l1, l2, lb, jt, it1, it2, ib, pt1, pt2, pb, it, git, lt
      integer :: i, off, j

      real(dp) :: dea
      complex(dp) :: z_one  , z_zero  
!        ,si(mbc%thsz,mbc%thsz)

      integer :: info, lwork
      complex(dp), allocatable :: work(:)
!       , gti(:,:)
      integer, allocatable :: ipiv(:)
      

      include '../Include/Atom.array' ! for dem only
      
!       if (.not. allocated(gti)) allocate(gti(mbc%thsz,mbc%thsz))
! !       
! !       if (mbc % thsz /= 0 .and. mbc%ovl) then
! !         si(:,:) = mbc%s(:,:)
! !         allocate(ipiv(mbc % thsz))
! ! 
! !         call zgetrf (mbc % thsz, mbc % thsz,  si, mbc % thsz, ipiv, info)
! !         if (info /= 0) stop 'zgetrf info /= 0'
! ! 
! !         allocate(work(1))
! !         lwork = -1
! ! 
! !         call zgetri(mbc % thsz, si, mbc % thsz, ipiv, work, lwork, info)
! !         if (info /= 0) stop 'zgetri i info /= 0'
! ! 
! !         lwork = work(1)
! !         deallocate(work)
! !         allocate(work(lwork))
! ! 
! !         call zgetri(mbc % thsz, si, mbc % thsz, ipiv, work, lwork, info)
! !         if (info /= 0) stop 'zgetri r info /= 0'
! ! 
! !         deallocate(work,ipiv)
! !       endif
      
      

      ! (z*I*S - (H + dH'*Gb*dH))^-1 = W (z*I - V)^-1 W'
       mbc % ac = 0.0_dp
      mbc % gt = 0.0_dp
      mbc % gt = mbc % h(:,:,isp)
!         do i = 1, mbc % thsz
!           do j = 1,mbc % thsz
!             mbc%gt(i,j) = sum( si(i,:)* mbc%h(:,j))
!           enddo
!         enddo
!       mbc%hp(:,:)=mbc%gt(:,:)

!       off = mbc % rlst(2)-1

!       do it = 1, mbc % n(2)
!          git = off + it
! 
!          dea = de(git)
!           if (mag) dea = dea + oppm(isp)*dem(git)
! !           print *, 'dea, oppm(isp)*dem(git)',dea, oppm(isp)*dem(git)
!          do lt = mbc % pth(it)+1, mbc % pth(it+1)
!             mbc % gt(lt,lt) = mbc % gt(lt,lt) + dea
!          end do
!       end do


      do it1 = 1, mbc % n(2)
         pt1 = mbc % pth(it1)
         nstt1 = mbc % pth(it1+1) - pt1

         do it2 = 1, mbc % n(2)
            pt2 = mbc % pth(it2)
            nstt2 = mbc % pth(it2+1) - pt2
            
            jt = mbc % apd(it1)
            ib = mbc % bpd(jt)
!             print *, it1,ib
!             print *, mbc % dh(:, :, jt)
            do while ( ib /= eol)
               pb = mbc % pbh(ib)
               nsttb = mbc % pbh(ib+1) - pb
!                print *, "nsttb",nsttb
!                   print *, "it1,ib,jt      :",it1,ib,jt
               do l1 = 1, nstt1
                  do l2 = 1, nstt2
                     do lb = 1, nsttb
                        mbc % gt(pt1+l1,pt2+l2) =  mbc % gt(pt1+l1,pt2+l2) &
                                    & - mbc % dh(l1, lb, jt) * mbc % a (pb+lb,pt2+l2)
                     end do
                  end do
               end do

               jt = jt + 1
               ib = mbc % bpd(jt)
            end do

         end do
      end do

       mbc % gt = - mbc % gt

      if (.not. mbc % ovl) then
         do i = 1, mbc % thsz
            mbc % gt(i,i) = mbc % gt(i,i) + ez
         end do
      else
         mbc % gt = mbc % gt+ (mbc % s * ez)
      end if


            
      if (mbc % thsz /= 0) then
          allocate(ipiv(mbc % thsz))

          call zgetrf (mbc % thsz, mbc % thsz, mbc % gt, mbc % thsz, ipiv, info)
          if (info /= 0) stop 'zgetrf info /= 0'

          allocate(work(1))
          lwork = -1

          call zgetri(mbc % thsz, mbc % gt, mbc % thsz, ipiv, work, lwork, info)
          if (info /= 0) stop 'zgetri i info /= 0'

          lwork = work(1)
          deallocate(work)
          allocate(work(lwork))

          call zgetri(mbc % thsz, mbc % gt, mbc % thsz, ipiv, work, lwork, info)
          if (info /= 0) stop 'zgetri r info /= 0'

          deallocate(work,ipiv)
      endif
      
!       if (mbc%sctb) then
!           do it1 = 1, mbc % n(2)          
!               do i = mbc%pth(it1)+1,mbc%pth(it1+1)
!                   mbc % dqchi(it1+mbc%n(1)) = mbc % dqchi(it1+mbc%n(1)) + &
!                         & real(zp * mbc % gt(i,i) * mbc % gt(i,i))
!               enddo
!           enddo
!       endif
           
       z_one= 1
       z_zero =0


          

! !         
!         gti(:,:) = mbc % gt(:,:)
!         mbc % gt = 0.0_dp
!         do i = 1, mbc % thsz
!           do j = 1,mbc % thsz
!             mbc%gt(i,j) = sum( gti(i,:)*mbc%s(:,j))
!           enddo
!         enddo
! !         
! !         gti(:,:) = mbc % gt(:,:)
! !         do i = 1, mbc % thsz
! !           do j = 1,mbc % thsz
! !             mbc%gt(i,j) = sum( gti(i,:)*mbc%sr(:,j))
! !           enddo
! !         enddo
! 
     
!       
      if (mbc % bhsz > 0 .and. mbc % thsz > 0) then
        call zgemm ('n', 'n', mbc % bhsz, mbc % thsz, mbc % thsz, z_one, &
                      & mbc % a, mbc % bhsz, mbc% gt, mbc % thsz, &
                      & z_zero, mbc % ac, mbc % bhsz)
      endif
      
      
!       !   This is to adjust the density matrix with the overlap for onsite elements required for charges, eprom, etc
! ! see Finnis (2007) eqn. 21
!       if (mbc%ovl .and. rho .and. mbc % thsz > 0) then
!           gti(:,:) = mbc % gt(:,:)
!           si(:,:) = mbc % s(:,:)
! 
!           call zgemm ('N', 'N', mbc % thsz, mbc % thsz, mbc % thsz, z_one, &
!                       & si, mbc % thsz, gti, mbc % thsz, &
!                       & z_zero, mbc %gt, mbc % thsz)
!       endif
! CARRIED OUT IN MBQ NOW

   end subroutine gtbi


   subroutine gb_x_dh(mbc,ez,isp)
      use ab_io
      use mod_ham

      type(mbtb_t), intent(inout) :: mbc
      complex(dp), intent(in) :: ez
      integer, intent(in) :: isp

      integer :: ib, ia, it, jb, jt, la, lb, lt, gia, gib, git,i,j

      complex(dp) :: gab(mxnstat, mxnstat)
!, gab2(mxnstat, mxnstat), gab1(mxnstat, mxnstat)

      mbc%a = 0.0_dp
      mbc%ggd = 0.0_dp
      do ia = 1, mbc % n(1)
         call assoc_ab(ia,isp)
         call assoc_ham(ia) ! for cluster, pcluster and decipher

         gia = mbc%rlst(1)-1 + ia ! gia == ia, not sure why am i doing the extra addition. consistency?

         call get_gab(ez, ia, ia, 1, gab,.true.)
!          IF (ISP ==2) THEN
!          print *, gab
!          ENDIF
         do la = 1, mbc%pah(gia+1)-mbc%pah(gia)      ! it is a little dirty to sed ggd and gg here
            mbc % ggd(mbc%pah(gia) + la) = gab(la,la)  ! but in this way gab is reused and another call get_gab is avoided
         end do


         do it = 1, mbc % n(2)
            jt = mbc % apd(it)
            ib = mbc % bpd(jt)

            git = mbc % rlst(2) -1 + it
!                       print *, 'ia,git:', ia,git
            do while ( ib /= eol)
               jb = decipher(ib) ! ib is local for region 1 as is decipher
!                print *, "ia,ib,jb", ia,ib,jb
!                 print *, 'jb:        ', jb
               if (jb /= 0) then

                  gib = mbc%rlst(1)-1 + ib !
                  gab = 0.0_dp
           !       gab1 = 0.0_dp
                                              
                  call get_gab(ez, ia, ib, jb, gab,.false.)
                                    
!          !          print *, "jbs",ia,ib,jb

           !        if (ia /= ib) then
           !          gab2 = 0.0_dp
           !          call assoc_ab(ib,isp)
           !          call assoc_ham(ib)
           !          jb2 = decipher(ia)
           !          call get_gab(ez, ib, ia, jb2, gab2)
           !          call assoc_ab(ia,isp)
           !          call assoc_ham(ia)
           !          
           !          
           !          do i = 1,mxnstat
           !             do j = 1,mxnstat
           !                 gab(i,j)=(gab1(i,j)+(gab2(j,i)))/2
           !             enddo
           !       enddo
                  
  !                 else
 !                    gab=gab1
!                   endif

                  
                  do lt = mbc % pth(it) + 1, mbc % pth(it + 1)
                     do la = mbc % pbh(ia)+1, mbc % pbh(ia + 1)
                        do lb = 1, mbc % pbh(ib + 1) - mbc % pbh(ib)
! CHANGE TO NEG
                           mbc % a(la, lt) = mbc % a(la, lt) &
                              & - (gab(la - mbc % pbh(ia), lb) )* mbc%dh(lt - mbc % pth(it),lb,jt)
!                               print *, ia,git,ib,(gab(la - mbc % pbh(ia), lb) ),mbc%dh(lt - mbc % pth(it),lb,jt)
                        end do
                     end do
                  end do

               end if
               jt = jt + 1
               ib = mbc % bpd(jt)
            end do
         end do
         
         
         
      end do
      
      
!               print *, "A matrix"
!          do i= mbc % pbh(22)+1,mbc % pbh (22 + 1)
!             print *, (mbc% a(i,j),j=mbc % pth(1)+1,mbc % pth(2))
!          enddo

!               print *, mbc%a
!        mbc%a=conjg(mbc%a)
      
   end subroutine gb_x_dh



   subroutine get_gab(ez, ia, ib, jb, gab,gdiag)
      use ab_io      
      use mod_g0n
      use mod_all_scalar, only : momflg, term

      complex(dp), intent(in) :: ez
      integer, intent(in) :: ia, ib, jb
      complex(dp), intent(out) :: gab(mxnstat,mxnstat)

      complex(dp) :: g0n(0:mrec+2), hia(0:mrec+1), hib(1:mrec+1)

      integer :: la, lb, nla, nstta, nsttb, lla, nma, ma, nmax, nrec, za, zb
      logical :: gdiag

      include "../Include/Atom.array"

      za = z(ia) ! offset = 0 region 1  ! These are inside because the number of cross-neighbours is expected to be low. 
      call states(za, nla, nstta, llista) !

      if (momflg == 1) then
         lla = 0
         do la = 1,nla
            nma = 2*llista(la) + 1
            mlista(lla+1:lla+nma) = la
            lla = lla + nma
         enddo
      else
         do la = 1,nstta
            mlista(la) = la
         enddo 
      endif

      gab = 0.0_dp

      nsttb = nstt(z(ib))

      do la = 1, nstta

         ma = mlista(la)
         nrec = lchain(ma)
         nmax = nrec + 1
         
!            CHANGE BREC FROM 0:NREC TO NREC+1
         call getg0n(ez, arec(0:nrec,ma), brec(0:nrec+1,ma), g0n, nmax, nrec, lainf(ma), lbinf(ma))

         hia(0:nmax) = g0n(0:nmax  ) * g0n(0:nmax)
         hib(1:nmax) = g0n(0:nmax-1) * g0n(1:nmax)

         if (.not. gdiag) then
            do lb = 1, nsttb
                if (la /= lb .or. jb /= 1) then 
                  gab(la,lb) =   sum(hia(0:nrec) * darec(0:nrec,la,lb,jb)) &
                            & + 2*sum(hib(1:nrec) * dbrec(1:nrec,la,lb,jb))
                  if (term == 1) gab(la,lb) = gab(la,lb) + 2 * hib(nmax) * dbrec(nmax,la,lb,jb)
                else
                
    !             CHANGE SAME
                  gab(la,la) = g00(ez, arec(0:nrec,ma), brec(0:nrec+1,ma), nrec, lainf(ma), lbinf(ma))
                end if
            end do
        else
            gab(la,la) = g00(ez, arec(0:nrec,ma), brec(0:nrec+1,ma), nrec, lainf(ma), lbinf(ma))
        endif

      end do

!       write(390,*) 'ez, ia, ib, jb',ez, ia, ib, jb
!       call print_c(gab, 390, 'gab'//achar(10))

   end subroutine get_gab


!    subroutine gbop(mbc, ez, isp)
!       use ab_io
!       use mod_ham
!       use mod_all_scalar, only : nbase
! 
!       type(mbtb_t), intent(inout) :: mbc
!       complex(dp), intent(in) :: ez
!       integer, intent(in) :: isp
! 
!       complex(dp) :: gab(mxnstat,mxnstat)
!       complex(dp), allocatable :: g(:,:)
! 
!       integer :: ia, ib, la, lb, jb, pa, pb
! 
!       allocate(g(mbc % bhsz, mbc % bhsz))
! 
!       g = 0.0_dp
! 
!       do ia = 1, mbc % n(1)
!          call assoc_ab(ia,isp)
!          call assoc_ham(ia)
! 
!          pa = mbc % pbh(ia)
! 
!          do jb = 1, pcluster(nbase+1)-1
!             ib = cluster(jb)
!             pb = mbc % pbh(ib)
! 
!             call get_gab(ez, ia, ib, jb, gab,.false.)
! 
!             do la = 1, mbc % pbh(ia+1) - pa
!                do lb = 1, mbc % pbh(ib+1) - pb
!                   g(pa+la,pb+lb) = gab(la,lb)
!                end do
!             end do
! 
! !             write(700,*),'ia,ib,jb,ez:', ia,ib,jb, ez
! !             call print_c(gab,700,'')
! !
!          end do
!       end do
! 
! !       write(600,*) 'ez', ez
! !       call print_c(g,600,'g'//achar(10),f='f8.3')
! 
!       deallocate(g)
! 
!    end subroutine gbop

   subroutine a_and_ac(mbc, ez, isp)
   

      type(mbtb_t), intent(inout) :: mbc
      complex(dp), intent(in) :: ez
      integer, intent(in) :: isp

!       call gbop(mbc,ez,isp)

      if (mbc%mxgf) then
          call gb_x_dh(mbc,ez,isp)
          call gtbi(mbc, ez, isp)
      else
!       if (.not. mbc%allbop) then 
          call gtdhgb(mbc,ez,isp)
      endif
      
      
   end subroutine a_and_ac


   subroutine mbeq(mbc,ef,sumz)
      use mod_all_scalar, only : mag, quiet, lef, dqmax, nbase
      use mod_atom_ar
      use ab_io
      use mod_ham

      type(mbtb_t), intent(inout) :: mbc

      real(dp) :: ef,zct
      
!       ,bpq(mbc%n(1))
      real(qp) :: sumz
      character(len=25) :: ffmt
      integer :: ia, jb
      include '../Include/Atom.array' ! for dq, zc, z

      
    if (mbc%mxgf .and. mbc%bhsz /= 0) then
        call recurse(.false.,mbc%n(1))  ! require full cluster rather than first shells used in BOP
        do ia = 1, mbc % n(1)
            call assoc_ham(ia)
            do jb = 1, pcluster(nbase+1)-1
                decipher(cluster(jb)) = jb
            end do
        end do
!     elseif (.not. mbc%allbop) then 
!         call recurse(.false.)
!         do ia =  1, mbc % nd
!             call assoc_ham(ia)
!             do jb = 1, pcluster(2)-1
!                 decipher(cluster(jb)) = jb
!             end do
!         end do
    endif 

!    if (.not. quiet) Print *, 'old dq',dq(mbc%n(1)+1:mbc%nd)
      
      ef = embeff(mbc,sumz,lef) 
      lef = ef

      
      dqmax = 0.0_dp
      if (.not. mbc%sctb) then
          do ia = 1, mbc % n(1)
!               print *, ia, dq(ia)
              dq(ia) = dq(ia) - zc(z(ia))
              if (abs(dq(ia)) > abs(dqmax)) dqmax = dq(ia)
          enddo 
          do ia = mbc % n(1)+1, mbc % nd
              if (mbc%orth(ia-mbc%n(1)) == 0) then
                  dq(ia) = dq(ia) - mbc%tbc%zc(z(ia))
!                   print *, 'tb',ia, dq(ia),dq(ia) + mbc%tbc%zc(z(ia))
              else
                  dq(ia) = dq(ia) - zc(z(ia))
!                   print *, 'tb',ia, dq(ia),dq(ia) + zc(z(ia))
              endif
               if (abs(dq(ia)) > abs(dqmax)) dqmax = dq(ia)
          enddo
      else 
         do ia = 1, mbc % n(1)
    !          if (ia <= mbc%n(1)) then
              
              dq(ia) = dq(ia) - zc(z(ia))
              if (abs(dq(ia)) > abs(dqmax)) dqmax = dq(ia)
!               if (abs(dq(ia)) > abs(dqmaxb)) dqmaxb = dq(ia)
              if (mbc%mxgf) mbc % dqchi(ia) = dq(ia)/mbc % dqchi(ia)
!               print *, 'bop',ia, dq(ia) + zc(z(ia)), dq(ia)
          enddo 
      
!           do ia = 1, mbc % n(1)
!               dq(ia) = dq(ia) - bpq(ia)
!               bpq(ia) = dq(ia) + bpq(ia)
!               print *, 'bop',ia,bpq(ia),dq(ia)
!               if (abs(dq(ia)) > abs(dqmax)) dqmax = dq(ia)
!               if (abs(bpq(ia)-zc(z(ia))) > abs(dqmaxb)) dqmaxb = bpq(ia)-zc(z(ia))
!           enddo 
          
          
!           dqmaxb = dqmax
          do ia = mbc % n(1)+1, mbc % nd 
!               qin(ia-mbc%rlst(2)+1) = mbc % qmpol(1,ia-mbc%rlst(2)+1)+mbc%tbc%zc(z(ia))
!               mbc % qmpol(1,ia-mbc%rlst(2)+1) = dq(ia) - mbc%tbc%zc(z(ia))
!               dq(ia) =  (dq(ia) - qin(ia-mbc%rlst(2)+1))
              if (mbc % orth(ia - mbc%n(1)) == 1) then
                  zct = zc(z(ia))
              else
                  zct = mbc % tbc % zc(z(ia))
              endif
              dq(ia) =  dq(ia) - (mbc % qmpol(1,ia-mbc%rlst(2)+1) + zct)              
              mbc % qmpol(1,ia-mbc%rlst(2)+1) = dq(ia) + mbc % qmpol(1,ia-mbc%rlst(2)+1)
              
              
!              if (abs(dq(ia)) > 1.0_dp ) then
!                  if (.not. quiet) then
!                      print *, 'atom',ia,'has large dq:',dq(ia)
!                      print *, 'scale by',(2*abs(dq(ia))+zc(z(ia)))/zc(z(ia))
!                      print *, 'dqchi',mbc % dqchi(ia),'=>',((10*abs(dq(ia))+zc(z(ia)))/zc(z(ia)))*mbc % dqchi(ia)
!                  endif
                  
!                  mbc % dqchi(ia) = ((10*abs(dq(ia))+zc(z(ia)))/zc(z(ia)))*mbc % dqchi(ia)
!              endif
     
              mbc % dqchi(ia) = mbc%qmpol(1,ia - mbc%n(1))/mbc % dqchi(ia)
!               mbc % dqchi(ia) = mbc % qmpol(1,ia-mbc%rlst(2)+1)/mbc % dqchi(ia)
!              if (.not. quiet) then
                    
!                    print *, 'ia',ia,ia - mbc%n(1),ia-mbc%rlst(2)+1
!                    print *,'orth,chi',mbc % orth(ia - mbc%n(1)),mbc % dqchi(ia)/mbc%qmpol(1,ia - mbc%n(1))
!                    print *, 'dq,qmpol,q',dq(ia) , mbc %qmpol(1,ia-mbc%n(1)) , (mbc % qmpol(1,ia-mbc%rlst(2)+1)+zct)
!                    print *, zct
!                    print *, ' '
!                    if (abs(dq(ia)) > 6) then
!                         print *,'HUGE DQ, DUMPING COORD '
!                         do jb = 1, mbc % nd
!                           print *,jb,mbc%ad(:,jb),de(jb),mg(jb)
!                         enddo
!              !           stop
!                    endif
!              endif
!               print *, 'tb',ia,mbc % qmpol(1,ia-mbc%rlst(2)+1),zct+mbc % qmpol(1,ia-mbc%rlst(2)+1),dq(ia)
              if (abs(dq(ia)) > abs(dqmax)) dqmax = dq(ia)
!               if (abs(mbc % qmpol(1,ia-mbc%rlst(2)+1)) > abs(dqmaxb)) &
!                               & dqmaxb = mbc % qmpol(1,ia-mbc%rlst(2)+1)
          enddo        
      endif

!       print *, 'mg',mg

!      if (.not. quiet) then
!          write(ffmt,"(a,i0,a)") '(a,":",', mbc % nd, '(x,f22.14))'
!          write(6,'(a)',advance='no') grnb
!          write(6, ffmt,advance='no') 'dq', dq(: mbc % nd)
!          write(6,'(a)') endc
!          if (mag) write(6, ffmt) 'dm', dem(:mbc % nd)
   !                if (mag) write(6, ffmt) 'mg', mg(:mbc%nd)
!          write(ffmt,"(a,i0,a)") '(a,":",', mbc%nd, '(x,f22.14))'
!          if(mbc%n(1) < 300 .and. mbc%n(1) > 0) write(6, ffmt) 'BOP dq', dq(:mbc%n(1))
!          if(mbc%n(2) < 300.and. mbc%n(2) > 0) write(6, ffmt) 'TB dq', dq(mbc%n(1)+1:mbc%nd)
          
!          print *,'mg =',mg(:mbc%n(2))
        
!      end if
!       call mbq

   end subroutine mbeq

   function embeff(mbc,locc,ef0) result(ef)
      use mod_all_scalar, only : totnia,nsp,mfac,quiet
      use topologia, only : mpmap, iproc,sdistrib, nproc, master

      type(mbtb_t), intent(inout) :: mbc
      real(qp), intent(in) :: locc
      real(dp) :: ef, ef0
      real(dp) :: eflo,efhi,nello,nelhi
      real(dp) :: nel,diff
      real(dp), parameter :: def = 1.0_dp, maxdn = 1.0e-8_dp
      
      integer :: info,liwork,lwork,isp
      integer, allocatable :: iwork(:)
      real(dp), allocatable :: cnw(:,:),ovlm(:,:),w(:),work(:)

 !     ncall = 0
!       do it = 1, 100
!           if (it > 1) nel = mbq(mbc, ef)
!
!           if (.not. quiet) then
!             write(6,'(/,"it:",i3)') it
!             write(6,'("ef,nf:",4(x,f22.14))') ef_prev, ef, dnf_prev, dnf
!           end if
!
!           if (abs(dnf) <= 1.0e-10_dp) exit
!           def = dnf * (ef - ef_prev) / ( dnf - dnf_prev )
!
!           ef_prev = ef
!           dnf_prev = dnf
!           ef = ef - def
!
!       end do

!       call sdistrib(mfac, nproc, mpmap, aproc)
      
      
      if (.not. mbc% mxgf .and. .not. mbc%allbop) then
          allocate(cnw(mbc%thsz,mbc%thsz),w(mbc%thsz))
          
          if (mbc%ovl) allocate(ovlm(mbc%thsz,mbc%thsz))
          do isp = 1,nsp
!               print *, 'ISP =',isp
              allocate(work(4))
              work = 1.0_dp
              lwork = -1
              
              cnw(:,:) = mbc%h(:,:,isp)
              
              if (mbc% ovl) then
                  allocate(iwork(4))
                  iwork = 1.0_dp
                  liwork = -1
                  
                  ovlm(:,:) = mbc%s(:,:)
                  call dsygvd(1,'V','U',mbc%thsz,cnw,mbc%thsz,ovlm,mbc%thsz,w,work,lwork,iwork,liwork,info) 
                  if (info /= 0) stop 'dsygvd i info /= 0'
              else
                  call dsyev('V','U',mbc%thsz,cnw,mbc%thsz,w,work,lwork,info)
                  if (info /= 0) stop 'dsyev i info /= 0'    
              endif
              
              
              lwork = int(work(1))
              deallocate(work)
              allocate(work(lwork))
              
              cnw(:,:) = mbc%h(:,:,isp)
              
              if (mbc%ovl) then
                  liwork = int(iwork(1))
                  deallocate(iwork)
                  allocate(iwork(liwork))
                  ovlm(:,:) = mbc%s(:,:)
                  call dsygvd(1,'V','U',mbc%thsz,cnw,mbc%thsz,ovlm,mbc%thsz,w,work,lwork,iwork,liwork,info)  
                  if (info /= 0) stop 'dsygvd r info /= 0'
                  deallocate(iwork)
              else
                  call dsyev('V','U',mbc%thsz,cnw,mbc%thsz,w,work,lwork,info)
                  if (info /= 0) stop 'dsyev i info /= 0' 
              endif
              
              mbc%cn(:,:,isp) = cmplx(cnw(:,:),0.0_dp,kind=dp)
              mbc%en(:,isp) = w(:)
              deallocate(work)
              
!             if (.not. quiet)  print *,'eigvals',isp, mbc%en(:,isp)
          enddo
          deallocate(cnw,w)
          if (mbc%ovl) deallocate(ovlm)
      endif
       

      if (.not. quiet) write(6,'(''LOCC = '',G12.5)') locc
      ef = ef0
      nel = mbq(mbc, 0, ef)
      if (.not. quiet) print *, 'ef,nel:',ef,nel
      if (nel < locc) then
         eflo = ef
         nello = nel
         efhi = eflo + def
   1       nelhi = mbq(mbc, 0, efhi)
         if (.not. quiet) print *, 'efhi,nelhi:', efhi,nelhi
         if (nelhi < locc) then
            eflo = efhi
            nello = nelhi
            efhi = efhi+def
            goto 1
         endif
      else
         efhi = ef
         nelhi = nel
         eflo = efhi - def
   4       nello = mbq(mbc, 0, eflo)
         if (.not. quiet) print *, 'eflo,nello:',eflo,nello
         if (nello > locc) then
            efhi = eflo
            nelhi = nello
            eflo = eflo-def
            goto 4
         endif
      endif


!    Use  binary sections to find the Fermi energy.

 2    ef = (eflo + efhi)*0.5_dp
 
!         if (mbc % sctb .and. diff < 1d6*maxdn) chi = .true.
 
      nel = mbq(mbc, 0, ef)
      if (.not. quiet) print *, 'ef,nel:',ef,nel
      if (nel > locc) then
         nelhi = nel
         efhi = ef
      else
         nello = nel
         eflo = ef
      endif
      diff = abs(nel-locc)
!      if (diff > maxdn .and. ncall < 8) goto 2
      if (diff > maxdn) goto 2

      nel = mbq(mbc, 1, ef)


      totnia = nel

     

!      if (.not. quiet) print *, 'ef,nel:',ef,nel
      
!       if (mbc % sctb .and. .not. chi) then
!           if (.not. quiet) write(6,'("ERROR: chis not calculated, increase diff tolerance in emb_sc.F90.")')
!           call panic()
!       endif

   end function embeff


   function mbq(mbc, mode, ef) result(cnel)
      use mod_all_scalar, only : kt, mfac, nsp,mag,quiet
      use mod_atom_ar
      use mod_ham
      use ab_io
      use topologia
      use mpi 
      use mod_clock

    include "../Include/Atom.array"

      type(mbtb_t), intent(inout) :: mbc   
      integer, intent(in) :: mode
      real(dp), intent(in) :: ef
      integer :: isp, la, ia

      real(dp) :: cnel,qi

      complex(dp) :: zp,zep,zpn
      real(dp) :: phase,w0,q(mbc%nd)
      integer :: m, p,nstta,ierr
      complex (dp) :: gab(mxnstat,mxnstat)
      complex (dp),allocatable :: rho(:,:),si(:,:)
      complex(dp) :: z_one  , z_zero 

      integer(8) :: c1, c2

      call system_clock(c1)

!      M     = NINT(MFAC*0.25D0*(EF-LAINF(LA)+2.0D0*LBINF(LA))/KT)
!        m     = nint(mfac*0.25_dp*(ef-lainf(la)+2*lbinf(la))/kt)
!       print *,ef, lainf(la), lbinf(la), (ef-lainf(la)+2.0_dp*lbinf(la))
!       m     = nint(mfac*0.25_dp*1.8_dp/kt)
      m = mfac
      if (.not. allocated(rho)) allocate(rho(mbc%thsz,mbc%thsz))
      if (.not. allocated(si)) allocate(si(mbc%thsz,mbc%thsz))

!       print *, 'm:',m
!       write(100+iproc,*) ef, lainf(la), lbinf(la)
!       write(100+iproc,*) '      ',ef, arec(0:lchain(la),la), brec(0:lchain(la)+1,la)
      if (mbc%sctb .and. mode == 1) mbc % dqchi = 0.0_dp
      
      dq = 0.0_dp
      do isp = 1, nsp        
         rho = 0.0_dp
         q = 0.0_dp
         mbc % ggd = 0.0_dp

         w0    = kt*(2*m)
         phase = pi/real(2*m, dp)
         zp    = cmplx(cos(phase),sin(phase), kind=dp)
!          zfac  = zp*zp

         do p = mpmap(iproc), mpmap(iproc+1)-1
         
            zpn = zp**(2*p + 1)
            zep  = ef+w0*(zpn-1)
            if (.not. mbc%allbop) then
                call gg_diag(mbc, zep, isp)
                
                rho(:,:) = rho(:,:) + zpn * mbc%gt(:,:)
            
                do ia = 1, mbc % nd
                  do la = mbc%pah(ia)+1, mbc%pah(ia+1)
                      if (mbc%mxgf.and. ia < mbc%rlst(2)) q(ia) = q(ia) + real(zpn * mbc % ggd(la))
                      if (mbc%sctb .and. mode == 1) mbc % dqchi(ia) = mbc%dqchi(ia)+real(zpn*mbc%ggd(la)*mbc%ggd(la))
                  end do
                end do
            else
                do ia = 1, mbc % nd
                    if (zc(z(ia)) /= 0) then
                        call assoc_ab(ia,isp)
                        call assoc_ham(ia) ! for cluster, pcluster and decipher

                        nstta = nstt(z(ia))
                        
                        call get_gab(zep, ia, ia, 1, gab,.true.)
                        
                        do la = 1,nstta      
                            q(ia) = q(ia) + real(zpn * gab(la,la))
                            if (mbc%sctb .and. mode == 1) mbc%dqchi(ia) = mbc%dqchi(ia) + real(zpn * gab(la,la) * gab(la,la))
                              
                        end do
                    endif
                    
                end do
            
            endif
            
!             zp   = zfac*zp
         enddo
          
        if (.not. mbc%allbop .and. mbc%thsz/=0) then
!   This is to adjust the density matrix with the overlap for "bond charge" 
!   required for charges, eprom, etc. see Finnis (2007) eqn. 21
            if (mbc%ovl) then
                z_one= 1.0_dp
               z_zero =0.0_dp
                si(:,:) = mbc % s(:,:)
!                 call zgemm ('N', 'N', mbc % thsz, mbc % thsz, mbc % thsz, z_one, &
!                     & si, mbc % thsz, rho, mbc % thsz, &
!                     & z_zero, mbc%gt, mbc % thsz)
                call zgemm ('N', 'N', mbc % thsz, mbc % thsz, mbc % thsz, z_one, &
                    & si, mbc % thsz, rho, mbc % thsz, &
                    & z_zero, mbc%gt, mbc % thsz)
            else
                mbc%gt(:,:) = rho(:,:)
            endif
            do ia = 1, mbc%n(2)
                do la = mbc%pth(ia)+1, mbc%pth(ia+1)
                    q(ia+mbc%n(1)) = q(ia+mbc%n(1)) + real(mbc%gt(la,la))
                enddo
            enddo
        endif
       
        
       if (mode == 0) then

         dq(:mbc%nd) = dq(:mbc%nd) + q

       else 
!    potential ISSUE: may need to reduce rho before i calc. q
        if (nproc > 1) then
            call mpi_allreduce(mpi_in_place,  q, mbc%nd, mpi_real8, mpi_sum, mpi_comm_world, ierr)
        endif
          
         q = 4 * kt * q
         
!          if(iproc == master) print *, 'spn',isp, 'charges',q(:mbc%nd)
         
         if (.not. mbc%mxgf) then
            if (mbc % allbop) then
                mbc%q(isp,:mbc%n(1)) = q(:mbc%n(1))
            else
                q(:mbc%n(1)) = mbc%q(isp,:mbc%n(1))
            endif
         endif
         
         dq(:mbc%nd) = dq(:mbc%nd) + q
         if (isp /= 1) then ! mag == true and isp == 2
            dq(:mbc%nd) = 0.5_dp * dq(:mbc%nd)
            mg(:mbc%nd) = dq(:mbc%nd) - q(:mbc%nd)
         end if
       end if
      end do
      
      if (mode == 1) then      
      if (mbc%sctb) then
          if (nproc > 1) then
              call mpi_allreduce(mpi_in_place,  mbc%dqchi, mbc%nd, mpi_real8, mpi_sum, mpi_comm_world, ierr)
          endif
          mbc % dqchi = -2 * kt * mbc % dqchi
          if (nsp == 1) mbc % dqchi = 2 * mbc % dqchi
      
      endif
     end if  
     deallocate(rho,si)
      
!       print *, iproc,'after dqchi',mbc%dqchi(:mbc%nd)
     if (.not. mbc%mxgf) then 
       cnel = sum(dq(mbc%n(1)+1:mbc%nd))
     else
       cnel = sum(dq(:mbc%nd))
     endif

     cnel = (4/nsp) * kt * cnel

     call mpi_allreduce(mpi_in_place,  cnel, 1, mpi_real8, mpi_sum, mpi_comm_world, ierr)

     if (.not. mbc%mxgf) cnel = cnel + sum(mbc%q)/nsp



      call system_clock(c2)

      if (.not. quiet) print *, 't_mbq:', real(c2-c1,8)/real(cr,8)

   end function mbq

   subroutine inctbh(ctsx,mbc,ef,force,mgo,nlmq,ecorr,itt)
      use ab_io
      use mod_ham
      use mod_atom_ar
      use mod_precision
      use mod_all_scalar
      use mod_const
      use topologia, only: iproc, master
!   Produces correction to ham. from self consistent charge transfer, higher multipoles have
!   not been included, but the possibility has been left fairly open.
      type(mbtb_t), intent(inout) :: mbc
      type(ctsx_t), intent(in) :: ctsx
      integer :: isp
      real(dp), intent(in) :: ef,mgo(mbc%nd)
      integer :: it, git,i,ilm,ilmp,ilmpp,lp,lpp,l,p,pp,nlm,nlmi1,nlmi,nlmj
      real(dp) :: U,Dr
      real(dp) ::  M,sumU,sumJ,sumM,ecorr,e2
      integer ::  zt,r,x,y,ib,jt,j,nlmq,nlmq1
      logical :: force
!       Br(mbc%n(2)),vm(9,mbc%n(2)),
      integer          :: kz(nlmq),kx1(nlmq),kx2(nlmq),ky1(nlmq),ky2(nlmq)
      integer   :: il,ilm1,ilm2,nlm1
      
      real(dp) :: cz(nlmq),cx1(nlmq),cx2(nlmq),cy1(nlmq),cy2(nlmq),cc
      
      real(dp), allocatable :: vmg(:,:,:), vj(:,:),pot0(:),vhub(:,:)
      real(dp), allocatable :: Ar(:),f(:,:) 
      
      integer :: llistt(3),nlt
      real(dp) :: istnt
      
      integer :: ierr,itt
      

      include "../Include/NebList.array"
      include "../Include/Atom.array"
      include "../Include/Force.array"
      include "../Include/PosVel.array"
      
      
      
      if (.not. allocated(vj)) allocate(vj(2,mbc%n(2)))
      if (.not. allocated(pot0)) allocate(pot0(mbc%n(2)))
      if (.not. allocated(vhub)) allocate(vhub(0:2,mbc%n(2)))
      if (mbc%ovl) then
          if (.not. allocated(Ar)) allocate(Ar(mbc%n(2)))
          Ar = 0.0_dp
      endif

      mbc%vm = 0.0_dp
      vj = 0.0_dp
      mbc%vu = 0.0_dp
      vhub = 0.0_dp
      pot0 = 0.0_dp
      sumJ = 0.0_dp
      sumU = 0.0_dp
      summ = 0.0_dp
      
!   value of e**2 in BOP units ( = Eryd*2a0)   
      e2=14.399645351469879_dp

!       mbc % hrs(:,:,1) = mbc % h0(:,:)
!       mbc % hrs(:,:,2) = mbc % h0(:,:)
!     if (.not. quiet) then
!         print *, 'qmpol:'
!         print *, mbc % qmpol(1,:)
!         print *, 'mg:'
!         print *, mgo(:)
!     endif
      
    if (.not. quiet) write(6,'(/,"inctbh: calculating new shifts to Ham.")')
      do it = 1, mbc % n(2)
          git = mbc % rlst(2) -1 + it
          zt = z(git)
          
          if (mbc%orth(it) == 1) then
              llistt(:) = llist(:,zt)
              nlt = nl(zt)
              istnt = istn(zt)
          else
              llistt(:) = mbc % tbc % llist(:,zt)
              nlt = mbc % tbc % nl(zt)
              istnt = mbc%tbc%istn(zt)
          endif
          
          U = ctsx %U(it)
          nlm = 1
          
          do l = 1, nlt
              vhub(llistt(l),it) = U * mbc % qmpol(1,it)
          enddo
          sumU = sumU + U * (mbc % qmpol(1,it))**2
          if (nsp == 2) then
              vj(1,it) = -0.5*istnt*(mbc % qmpol(1,it)+mgo(git))
              vj(2,it) = -0.5*istnt*(mbc % qmpol(1,it)-mgo(git))
              sumJ = sumJ - istnt * (0.25d0 * (mbc % qmpol(1,it) + mgo(git))**2 &
              &                +  0.25d0 * (mbc % qmpol(1,it) - mgo(git))**2)
             
          endif
!            print *, 'vj', it, vj(:,it)
!     Only calculates monopoles (dq)
          nlmi  = 1
          if (force ) then
  !           li1 = li + 1
              nlmi1  = 4
          else
              nlmi1  = nlmi
          endif
          
          do  ib = 1, mbc % n(2)
              nlmj = 1
  !     Need to figure out the indexing here
              i = ctsx%struxidx(ib,it)  
              do ilm = 1, nlmi1
                  mbc%vm(ilm,it) = mbc%vm(ilm,it) + &                        
                  &  e2*sum(ctsx%struxd(i+(ilm-1)*nlmj:i+ilm*nlmj-1)*mbc%qmpol(1:nlmj,ib))
              enddo                
              if (mbc % ovl) Ar(it) = Ar(it) + e2*ctsx%struxd(i)*mbc%qmpol(1,ib)
          enddo
          
          if (mbc % ovl) then
              mbc%vu(it) =U*mbc%qmpol(1,it)
              pot0(it) = Ar(it) + mbc%vu(it)
          endif

!            This is only the sumM for the monopolar case     
         sumM = sumM + mbc % qmpol(1,it) * mbc%vm(1,it)
         
         if (mbc%orth(it) == 0) then
         
!          if (.not. quiet) print *, 'ATOM = ', it
          do isp = 1,nsp
!                   if (.not. quiet .and. isp == 1) then
!                    !   print *, 'SPN = ', isp
!                       print *, '  VMAD=', mbc%vm(:,it)
!                       print *, '  Vhub=', vhub(:,it)
!                       print *, '  Vston=', vj(:,it)
!                       
!                   endif
                  
                  x = 1
                  r = 0
                  do ilmp = 1,9 
                      y = 1
                      x = x + r
                      r = 0
                      do ilmpp = 1,9
                          lp = ll(ilmp)
                          lpp = ll(ilmpp)
                          p = 0
                          pp = 0
                          do l = 1, nlt
                              if (llistt(l) == lp) p = 1
                              if (llistt(l) == lpp) pp = 1
                          enddo
                          
                          if (p == 1 .and. pp == 1) then
                              if (ilmpp == ilmp) then
                                    mbc % h(mbc % pth(it)+x,mbc % pth(it)+y,isp) = &
                                &   mbc % h(mbc % pth(it)+x,mbc % pth(it)+y,isp) &
                                                    &  + vhub(ll(ilmp),it) + vj(isp,it)
                                    
                                    
                              endif 
                              do ilm = 1, nlm 
                                  l = ll(ilm)
                                  call getM(zt,l,lp,lpp,M)
                                  mbc % h(mbc % pth(it)+x,mbc % pth(it)+y,isp) = &
                                  & mbc % h(mbc % pth(it)+x,mbc % pth(it)+y,isp) &
                                  &   + mbc%vm(ilm,it) * M  * ctsx%gaunt(ilmpp,ilmp,ilm)   
                              enddo
                              y = y + 1
                              r = 1
                          endif
                      enddo
                  enddo
                
               ! do ilmp = mbc % pth(it) +1,mbc % pth(it+1 )
               !     if (.not. quiet) print *, mbc % h(ilmp,mbc % pth(it) +1:mbc % pth(it+1),isp)
               ! enddo
                
                
          enddo
         else
            do isp =1,nsp
              do ilmp = mbc%pth(it)+1,mbc%pth(it+1)
                  mbc % h(ilmp,ilmp,isp) = mbc % h(ilmp,ilmp,isp) +de(git) +  oppm(isp)*dem(git)
              enddo
            enddo
         endif
      enddo
!    The second order correction to the energy (E_2) from SCTB
      ecorr = 0.5d0*(sumM + sumU)
 ! + sumJ)
!       if (.not. quiet) print *, 'sumM,sumU,sumJ',summ,sumu,sumj
 
      emag = 0.5d0*sumJ
      
     ! call mpi_barrier(mpi_comm_world,ierr)
      
!       if (itt > 1) STOP 'INCTBH'
      
!   including offsite terms (eq 7.81 Finnis' book)
      if (mbc%ovl) then
          do it = 1, mbc % n(2)
!               if (.not. quiet) print *, 'it:',it
              jt = mbc % apt(it)+1
              ib = mbc % bpt(jt)             
              do while (ib /= eol)
!                  if (.not. quiet) print *,'ib:',ib
                  do i = mbc % pth(it)+1,mbc % pth(it+1) 
                      do j = mbc % pth(ib)+1, mbc % pth(ib+1)   
                          Dr = (pot0(it) + pot0(ib))*mbc%s(i,j)*0.5d0
                          do isp = 1, nsp 
                              mbc % h(i,j,isp)=mbc%h(i,j,isp) + Dr
                          enddo
                      end do
                  end do
                  
                  
              !    if (.not. quiet) then
              !        do isp = 1, nsp 
              !            print *, 'OVL SHIFTS',isp
              !            do i = mbc % pth(it)+1,mbc % pth(it+1) 
              !                print *, mbc % h(i, mbc % pth(ib)+1: mbc % pth(ib+1)   ,isp)
              !            enddo
              !        enddo
              !    endif
                  
                  
                  jt = jt + 1
                  ib = mbc % bpt(jt)
              end do
          enddo
      endif
      
      
      if (force) then
          if (.not. quiet) write(6,'("inctbh: calculating contribs to force",/)')
          nlmq1 = 4
          if (.not. allocated(vmg)) allocate(vmg(nlmq,mbc%n(2),3))
          if (.not. allocated(f)) allocate(f(3,mbc%n(2)))
          
          do ib = 1, mbc%n(2)
              nlm1 = 4
              call tbshfl(1,nlmq1,nlm1,1,mbc%vm(1,ib))
          enddo
          
          do  ilm = 1, nlmq
              call scglp1(ilm,kz(ilm),cz(ilm),kx1(ilm),kx2(ilm), &
              &   cx1(ilm),cx2(ilm),ky1(ilm),ky2(ilm),cy1(ilm),cy2(ilm))
          enddo
          
          vmg(1:nlmq,1:mbc%n(2),1:3) = 0d0 
          do ib = 1, mbc%n(2)
              nlm = 1
              do  ilm = 1, nlm
                  vmg(ilm,ib,1) = vmg(ilm,ib,1) + cx1(ilm)*mbc%vm(kx1(ilm),ib) &
                           &       +cx2(ilm)*mbc%vm(kx2(ilm),ib)
                  vmg(ilm,ib,2) = vmg(ilm,ib,2) + cy1(ilm)*mbc%vm(ky1(ilm),ib) &
                           &       +cy2(ilm)*mbc%vm(ky2(ilm),ib)
                  vmg(ilm,ib,3) = vmg(ilm,ib,3) + cz(ilm)*mbc%vm(kz(ilm),ib) 
              enddo
          enddo
          
          do ib = 1, mbc%n(2)
              nlm = 1
              nlm1 = 4
              call tbshfl(0,nlmq1,nlm1,1,mbc%vm(1,ib))
              do i = 1, 3
                  call tbshfl(0,nlmq,nlm,1,vmg(1,ib,i))
              enddo
          enddo   
          
          do il = 0, ll(nlmq) 
              ilm1 = il*il+1
              ilm2 = (il+1)**2
              cc = dsqrt(dfloat((2*il+1)*(2*il+3)))
              vmg(ilm1:ilm2,1:mbc%n(2),1:3) = vmg(ilm1:ilm2,1:mbc%n(2),1:3)*cc
          enddo
          
          f(1:3,:) = 0d0
          do it = 1, mbc%n(2)
              git = mbc % rlst(2) -1 + it
              nlm = 1
              do ilm = 1, nlm 
                  f(1:3,it) = f(1:3,it) - vmg(ilm,it,1:3)*mbc%qmpol(ilm,it)
              enddo
!               if (.not. quiet) print *,it,'Fmad',f(:,it)
              
               fbs(:,git) = fbs(:,git) + f(:,it)
          enddo
      
          deallocate(vmg,f)
          
         
      endif

      
      deallocate(vj,pot0,vhub)
      
      
   end subroutine inctbh   
   
   
   subroutine getM(z,l,lp,lpp,M)
! C ----------------------------------------------------------------------
! Cr Remarks
! Cr  The tight binding parameters in qpol are as follows
! Cr  qpol(1) = M(011) = M(101)
! Cr  qpol(2) = M(112)
! Cr  qpol(3) = M(022) = M(202)
! Cr  qpol(4) = M(121) = M(211)
! Cr  qpol(5) = M(222)
! Cr  qpol(6) = M(123)
! Cr  qpol(7) = M(224)
! Cr  These are converted into Stone's definitions by multiplying the
! Cr  values from ctrl by \sqrt{4\pi/(2l+1)}
! C ----------------------------------------------------------------------
      use mod_conf, only : rtc
      use mod_atom_ar

      implicit none
! C Passed Parameters
      integer          :: ilm,ilmp,ilmpp,n
      real(dp),dimension(:),pointer :: qpl
      real(dp) :: M
! C Local Variables
      integer   :: ll,l,lp,lpp,z
      real(dp)  :: fourpi,sqr4pi,fac,dsqrt,datan
  
      fourpi = 16d0*datan(1d0)
      sqr4pi = dsqrt(fourpi)

!       l = ll(ilm)
!       lp = ll(ilmp)
!       lpp = ll(ilmpp)
      fac = dsqrt(fourpi/(2*l + 1))
      
      qpl => rtc % tbham_conf % a(asort(zatype(z))) % qpol

      M = 0d0
      if (lp .eq. 1 .and. lpp .eq. 0 .and. l .eq. 1) M = qpl(1)
      if (lp .eq. 0 .and. lpp .eq. 1 .and. l .eq. 1) M = qpl(1)
      if (lp .eq. 1 .and. lpp .eq. 1 .and. l .eq. 2) M = qpl(2)
      if (lp .eq. 2 .and. lpp .eq. 0 .and. l .eq. 2) M = qpl(3)
      if (lp .eq. 0 .and. lpp .eq. 2 .and. l .eq. 2) M = qpl(3)
      if (lp .eq. 2 .and. lpp .eq. 1 .and. l .eq. 1) M = qpl(4)
      if (lp .eq. 1 .and. lpp .eq. 2 .and. l .eq. 1) M = qpl(4)
      if (lp .eq. 2 .and. lpp .eq. 2 .and. l .eq. 2) M = qpl(5)
      if (lp .eq. 1 .and. lpp .eq. 2 .and. l .eq. 3) M = qpl(6)
      if (lp .eq. 2 .and. lpp .eq. 1 .and. l .eq. 3) M = qpl(6)
      if (lp .eq. 2 .and. lpp .eq. 2 .and. l .eq. 4) M = qpl(7)
      M = fac*M
      if (l .eq. 0) then
        if (lp .eq. lpp) then
          M = sqr4pi
        endif
      endif
      
   end subroutine getM
   
   
   function ll(ilm) result(l)
!  Function to produce angular momentum from matrix index
      integer, intent(in) :: ilm
      integer :: l
     

     if (ilm == 1) then
        l = 0
     else if (ilm < 5) then
        l = 1
     else if (ilm < 10) then
        l = 2
     else if (ilm < 17) then
        l = 3
     else 
        l = 4
     endif
    
   end function ll
   
   
!    This was started to produce TB multipoles, but was not finished
   
!    subroutine tbmpol(mbc,ef)
!       use ab_io
!       use mod_ham
!       use mod_atom_ar
!       use mod_precision
!       use mod_all_scalar
!       use mod_const
!       
!       type(mbtb_t), intent(inout) :: mbc
!       complex(dp) :: ez, zp
!       integer :: isp
!       real(dp), intent(in) :: ef
!       integer :: it, git
!       complex(dp) :: zp,zfac,zep
!       real(dp) :: phase, w0, rhoc(9,9,mbc%n(2)),rhoc1(9,9),zt
!       integer :: m, p 
!       complex(dp) :: gab(mxnstat, mxnstat)
! 
!       include "../Include/NebList.array"
!       include "../Include/Atom.array"
!       include "../Include/Force.array"
!       include "../Include/PosVel.array"
! 
!       
!       do isp = 1, nsp         
!          mbc % rho = 0.0_dp
!          w0    = kt*(2*m)
!          phase = pi/real(2*m, dp)
!          zp    = cmplx(cos(phase),sin(phase), kind=dp)
!          zfac  = zp*zp
! 
!          do p = 0, m-1
!          
!             zep  = ef+w0*(zp-1)
! 
!             call gtbi(mbc, ez, isp,.false.)
! 
!             do it = 1, mbc % n(2)
!                 gab = 0.0_dp
!                 gab(1 : mbc % pth (it + 1) - mbc % pth(it), 1 : mbc % pth (it + 1) - mbc % pth(it)) =  &
!                     & mbc % gt(mbc % pth(it) + 1 : mbc % pth(it + 1), mbc % pth(it) + 1 : mbc % pth(it + 1))
!                  
!                 rhoc(:,:,it) = rhoc(:,:,it) + real(zp*gab)
!             end do
!       
!             zp   = zfac*zp
!          enddo
!       end do
! 
!       do it = 1, mbc % n(2)
!         git = mbc % rlst(2) -1 + it
!         zt = z(git)
!         
! !    This restructures the rho matrix into full spd format     
!         rhoc1(:,:) = rhoc(:,:,it)
!         rhoc(:,:,it) = 0.0_dp
!         if (mbc % tbc % llist(1,zt) == 2) then
!           rhoc(5:,5:,it) = rhoc1(:5,:5)
!         else if (mbc % tbc % llist(1,zt) == 1) then
!           rhoc(2:,2:,it) = rhoc1(:8,:8)
!         else 
!           rhoc(1,1,it) = rhoc1(1,1)
!           if (mbc % tbc % llist(2,zt) == 1) then
!             rhoc(:,:,it) = rhoc1(:,:)
!           else if (mbc % tbc % llist(2,zt) == 2) then
!             rhoc(1,5:9,it) = rhoc1(1,2:6)
!             rhoc(5:9,1,it) = rhoc1(2:6,1)
!             rhoc(5:9,5:9,it) = rhoc1(2:6,2:6)
!           endif
!         endif
!         
!         enddo
!         
! !         nstt = mbc % pth(it+1) - mbc % pth(it)
!         do ilm = 2,9
!           do ilmp = 1,9
!               do ilmpp = 1,9
! !               still have to bring in the gaunt and deltas
! !                 mbc % qmpol(ilm,it) = mbc % qmpol(ilm,it) + rhoc(ilmp,ilmpp,it)*delta*gaunt(ilmp,ilmpp,ilm)      
!               enddo
!           enddo
!         enddo
!       enddo
!       
!    end subroutine tbmpol    
   
   
   subroutine areduce(n, pos, full, cnds)
      integer, intent(in) :: n, pos(1:)
      real(dp), intent(in) :: full(1:)
      real(dp), intent(out) :: cnds(1:)

      integer :: i

      do i = 1, n
         cnds(i) = cnds(i) + sum(full(pos(i)+1 : pos(i+1)))
      end do

   end subroutine areduce


   subroutine mbldh(mbc,ctsx)
       use ab_io , only : init_ab , free_ab , assoc_ab
      use mod_ham, only : init_ham, free_ham, assoc_ham
      use mod_chi, only : init_chi, free_chi
      use mod_all_scalar, only : nsp, nd, quiet,ninert,forces,rcut,rprune,totnd,mfac
      use topologia, only : mpmap, iproc,sdistrib, nproc, master
      use mod_conf, only : rtc
      use mod_const


      type(mbtb_t), intent(inout) :: mbc
      type(ctsx_t), intent(inout) :: ctsx
      integer :: ia, ja, jb, p, mp, ib, ir,i,j,it
      integer :: loc_clusiz
      
      integer :: info, lwork
      real(dp), allocatable :: rwork(:)
      complex(dp), allocatable :: work(:),si(:,:)
      integer, allocatable :: ipiv(:)
      
      real(dp) :: rcut1,ett(3,3)
      integer  :: rlim1(3),maxneigh,nat,apti(mxtotnd),bpti(mxnbl)
      logical :: stressed
            
      integer :: struxsize, dstrxsize,nlmq1,nlmq,ierr
      integer, allocatable :: indxcg(:),jcg(:)
      real(dp),allocatable :: cg(:),cy(:)
      
!       integer :: map1(mxtotnd), z1(mxtotnd), zinert1(minert), mapi1(mxtotnd)

      
!       double precision, allocatable :: w(:),work(:),sd(:,:),sv(:,:)
!       character :: num,up
!       integer ::  n,lda,lwork,info
!       
      include "../Include/PosVel.array"
      include "../Include/NebList.array"    
      include "../Include/Atom.array"
      include "../Include/ag.conn"

      mbc % nd = mbc % rlst(nr+1)-1



      do ir = 1, nr
         mbc % n(ir) = mbc%rlst(ir+1) - mbc%rlst(ir)
      end do

      if (.not. allocated(mbc % ad)) allocate(mbc % ad(3, mbc % nd))

      if (sum(rlim) == 0 .and. mbc%mxgf) then
        if (ninert == 0) then
            nat = mbc % n(1)
        else
            nat = mbc % nd + ninert
        endif
        if (.not. allocated(mbc % apb)) allocate(mbc % apb(nat+1))
        if (.not. allocated(mbc % bpb)) allocate(mbc % bpb(nat*(mxnnb+1)))
        
      else if (sum(rlim) == 0 .and. .not. mbc%mxgf) then
        nat = mbc % nd + ninert
        if (.not. allocated(mbc % apb)) allocate(mbc % apb(nat+1))
        if (.not. allocated(mbc % bpb)) allocate(mbc % bpb((nat)*(mxnnb+1)))      
      else
        if (.not. allocated(mbc % apb)) allocate(mbc % apb(mxtotnd))
        if (.not. allocated(mbc % bpb)) allocate(mbc % bpb(mxnbl))
        nat = mxtotnd-1
      endif
      
      if (.not. allocated(mbc % apt)) allocate(mbc % apt(mbc % n(2)+1))
      if (.not. allocated(mbc % apd)) allocate(mbc % apd(mbc % n(2)+1))

!       if (.not. mbc%sctb) then
!          if (.not. allocated(mbc % bpt)) allocate(mbc % bpt(mbc % n(2)*(mxnnb+1)))
!          if (.not. allocated(mbc % bpd)) allocate(mbc % bpd(mbc % n(2)*(mxnnb))) ! there will be at least one neighbour in its own region 
!      else
         if (.not. allocated(mbc % bpt)) allocate(mbc % bpt(mbc % n(2)*(80)))
         if (.not. allocated(mbc % bpd)) allocate(mbc % bpd(mbc % n(2)*(80)))
!      endif
      
      mbc % apt = 0
      mbc % bpt = 0
      
      mbc % apd = 0
      mbc % bpd = 0
      
      mbc % ad = ad(:, :mbc % nd)
      
!      if (.not. quiet) then
!          print *, 'in'
!          call print_neblist(1, mbc % nd, aptr, bptr)   
!      endif
!       if (mbc%mxgf) then
!           call extract_boplist(aptr, bptr,mbc % rlst, mbc % apb, mbc % bpb,nat)
!       else
!           mbc % apb = aptr
!           mbc % bpb = bptr
!       endif
      call extract_boplist(aptr, bptr,mbc % rlst, mbc % apb, mbc % bpb,nat,mbc%mxgf)
      
      call extract_tblist(aptr, bptr, 2, 1, mbc % rlst, mbc % apd, mbc % bpd,mbc%mxgf) ! list of atoms from r2 which have neighbours in r1


 !  This allows TB atoms to have a greater range and neighbour list without overflowing BOP's clusters
      rcut1 = rcut
      do i = 1, rtc % tbham_conf % nb
        do j = 0,9
          if (rtc % tbham_conf % b(i-1) % h(j) % tl % rc > rcut1) then
            rcut1 = rtc % tbham_conf % b(i-1) % h(j) % tl % rc
          endif
        enddo
      enddo
!       print *, 'RCUT',rcut1
! 
      if (rcut1 > rcut) then
        rlim1 = 0
        ett = 0.0_dp
        stressed = totstr /= 0.0_dp .and. any(et /= 0.0_dp)      
        if (stressed) ett = et*totstr
        apti = 0
        bpti = 0
        call bldnebt(.false., rlim1, lena, a, rprune, rcut1, stressed, ett, &
                    &  nd, ad, z, ninert, adinert, zinert, &
                    &  map, apti, bpti, maxneigh, totnd, nne, mapi,100)
      else
        apti(:) = aptr(:size(apti))
        bpti(:) = bptr(:size(bpti))
      endif
      

      
      call extract_tblist(apti, bpti, 2, 2, mbc % rlst, mbc % apt, mbc % bpt,mbc%mxgf)
      

     
!      if (.not. quiet) then
!          print *, 'out, 1'
!          call print_neblist(1, mbc % n(1), mbc % apb, mbc % bpb)
!          print *, 'out, 2'
!          call print_neblist(mbc % rlst(2), mbc % n(2), mbc % apt, mbc % bpt )
!          print *, 'out, 2->1'
!          call print_neblist(1, mbc % n(2), mbc % apd, mbc % bpd)
!      endif

      if (mbc % mxgf) then
!           mpmap(iproc) = 0
!           mpmap(iproc+1) = mbc % n(1)
          nat = mbc % n(1)
!           call sdistrib(mbc%n(1), nproc, mpmap, aproc)
      else
!           mpmap(iproc) = 0
!           mpmap(iproc+1) = mbc % nd
          nat = mbc%nd
!           call sdistrib(mbc%nd, nproc, mpmap, aproc)
      endif
      
            

      
!       print *, iproc
!       print *, mpmap
!   ISSUE: THESE ARE GLOBAL
      bptr = 0 
      aptr = 0
      bptr(:size(mbc % bpb))   = mbc % bpb
      aptr(:size(mbc % apb))   = mbc % apb

!       if (iproc == master) then
!           call init_ab (1, nat, nsp)
!           call init_ham(1, nat)
! 
!       endif
!       call mpi_barrier(mpi_comm_world,ierr)

      call init_ab (1, nat, nsp)
      call init_ham(1, nat)
       
      do ia = 1, nat
          if (zc(z(ia)) /= 0) then
              call assoc_ab(ia,1) 
              call assoc_ham(ia)

              call getnch(ia)   !Find number of linear chains per site.
              call bldclus(ia)
              !           loc_clusiz = loc_clusiz + pcluster(nbase+1)-1
              call bldh()   
          endif
      end do

      
      
      call sdistrib(mfac, nproc, mpmap, aproc)
      
   
      if (.not. allocated(mbc % orth)) allocate(mbc % orth(mbc % n(2)))
      mbc % orth(:) = 0 
      if (mbc%ovl) then
          do ia = 1, mbc%n(2)
            if (mbc%bpd(mbc%apd(ia)) /= eol .and. zc(z(ia+mbc%n(1))) /= 0)&
               & mbc%orth(ia) = 1
!              if(.not. quiet) then
!                  write(6,'("Orth:",x,i3,x,i0)') ia, mbc%orth(ia)
!!                   print *, mbc%bpd(mbc%apd(ia):mbc%apd(ia+1)-1)
!              endif
          enddo
      endif

! End BOP r1


!       TB r2

      if (.not. allocated(mbc % pbh)) allocate(mbc % pbh( mbc % n(1)+1))
      call mpos(mbc%rlst(1), mbc%rlst(2)-1, mbc%pbh,mbc)
      mbc % bhsz = mbc % pbh(mbc % n(1)+1)
!       print *, 'pbh', mbc % pbh

      if (.not. allocated(mbc % pth)) allocate(mbc % pth( mbc % n(2)+1))
      call mpos(mbc%rlst(2), mbc%rlst(3)-1, mbc%pth,mbc)
      mbc % thsz = mbc % pth(mbc % n(2)+1)
!        print *, 'pth', mbc % pth

      mbc % ahsz = mbc % bhsz + mbc % thsz
      if (.not. allocated(mbc % pah)) allocate(mbc % pah( mbc % nd+1))
      mbc % pah(:mbc % rlst(2)) = mbc % pbh
      mbc % pah(mbc % rlst(2) : mbc % rlst(3)) = mbc%bhsz + mbc % pth
!        print *, 'pah', mbc % pah
!       print *, 'thsz:', mbc % thsz
        
     

!       print *, 'thsz:', mbc%thsz
      if (.not. allocated(mbc%h0)) allocate(mbc%h0(mbc%thsz,mbc%thsz))
      
      if (.not. allocated(mbc% h)) allocate(mbc%h(mbc%thsz,mbc%thsz,nsp))
   
!       if (.not. allocated(mbc%hp)) allocate(mbc%hp(mbc%thsz,mbc%thsz))
!       print *, 'h'
      if (.not. allocated(mbc%dh)) allocate(mbc%dh(mxnstat,mxnstat,mbc%apd(mbc%n(2)+1)))
!       print *, 'dh'
      if (mbc%ovl) then
          if (.not. allocated(mbc%s))  allocate(mbc%s(mbc%thsz,mbc%thsz))          
!           if (.not. allocated(mbc%hs))  allocate(mbc%hs(mbc%thsz,mbc%thsz))
      endif

      call bldmtbmat(mbc % rlst(2), mbc % rlst(3)-1, mbc % pth, mbc % apt, mbc % bpt, mbc % h0,mbc,.false.)
       
      if (mbc%ovl) call bldmtbmat(mbc % rlst(2), mbc % rlst(3)-1, mbc % pth, mbc % apt, mbc % bpt, mbc % s,mbc,.true.)
       if (mbc%bhsz /= 0) then
          call bldsparsedh(mbc % apd, mbc % bpd, 2, 1, mbc % rlst,  mbc%dh ,mbc)
       endif
       
       
       if (mbc % sctb) then
       
          nlmq = 1
          if (.not. forces) then
              nlmq1 = 1
          else
              nlmq1 = 4
          endif
          
          if (.not. allocated(indxcg)) allocate(indxcg(7400))
          if (.not. allocated(jcg)) allocate(jcg(62200))
          if (.not. allocated(cg)) allocate(cg(62200))
          if (.not. allocated(cy)) allocate(cy(289))
      
          if (.not. allocated(mbc%qmpol)) allocate(mbc%qmpol(9,mbc%n(2)))  
          if (.not. allocated(ctsx%u)) allocate(ctsx%u(mbc%n(2)))
          
          if (.not. allocated(mbc%vm)) allocate(mbc%vm(9,mbc%n(2)))
          if (.not. allocated(mbc%vu)) allocate(mbc%vu(mbc%n(2)))
      
          if (.not. allocated(mbc % dqchi)) allocate(mbc % dqchi(mbc%nd))
          
          call bldU(ctsx,mbc%n(2),mbc%rlst(2))
          
          call sylmnc(cy,16)
          call scg(9,cg,indxcg,jcg)
          call makcgn(9,9,25,indxcg,jcg,cg,ctsx%gaunt)
          
          if (.not. allocated(ctsx%struxidx)) allocate(ctsx%struxidx(mbc%n(2), mbc%n(2)))
          if (.not. allocated(ctsx%dstrxidx)) allocate(ctsx%dstrxidx(1,1))
          call mkstrxidx(forces, mbc%n(2),struxsize,dstrxsize,ctsx%struxidx,ctsx%dstrxidx)
          
!     dstrxd is not currently used, but may be in the future for pressure calcs.
          if (.not. allocated(ctsx%struxd)) allocate(ctsx%struxd(struxsize))
          if (.not. allocated(ctsx%dstrx))allocate(ctsx%dstrxd(dstrxsize))
          
          
          call mkstrxd(forces,mbc%n(2),mbc%rlst(2)-1, 1, nlmq1, &
                &     indxcg,jcg,cg,cy,ctsx%struxd,ctsx%dstrxd, ctsx%struxidx, ctsx%dstrxidx, struxsize, dstrxsize)
               
          deallocate(indxcg,jcg,cg,cy)
          
      endif
 
       
!       if (mbc % thsz /= 0 .and. mbc%ovl) then
! !         if (.not. allocated(si))  allocate(si(mbc%thsz,mbc%thsz))
!         if (.not. allocated(mbc%o))  allocate(mbc%o(mbc%thsz,mbc%thsz))
!         mbc % o(:,:) = mbc%s(:,:)
!         
!         do i = 1, mbc % thsz
!             mbc % o(i,i) =  mbc % o(i,i) - 1
!         enddo
!         allocate(ipiv(mbc % thsz))

!         call zgetrf (mbc % thsz, mbc % thsz,  si, mbc % thsz, ipiv, info)
!         if (info /= 0) stop 'zgetrf info /= 0'
! 
!         allocate(work(1))
!         lwork = -1
! 
!         call zgetri(mbc % thsz, si, mbc % thsz, ipiv, work, lwork, info)
!         if (info /= 0) stop 'zgetri i info /= 0'
! 
!         lwork = work(1)
!         deallocate(work)
!         allocate(work(lwork))
! 
!         call zgetri(mbc % thsz, si, mbc % thsz, ipiv, work, lwork, info)
!         if (info /= 0) stop 'zgetri r info /= 0'
! 
!         deallocate(work,ipiv)
!         
!         do i = 1, mbc % thsz
!           do j = 1, mbc % thsz
!               mbc % hp(i,j) = sum(mbc%h(i,:)*si(:,j))
!           enddo
!         enddo
!       endif
      

       
!      USED TO CHECK EIGVALS OF HAM/OVL  
!        print *, 'HAMILTONIAN'
!        PRINT *, MBC%H    
!        print *, 'OVERLAP'
!        PRINT *, MBC%S    

!        allocate(sv(mbc%thsz,mbc%thsz))
!        
!        sv(:,:) = mbc % s(:,:)
!        num='v'
!        up='l'
!        lda = mbc % thsz
!        n = mbc % thsz
!         allocate(work(1))
!         lwork = -1
!         allocate(W(mbc%thsz))
!       call dsyev(num,up,n,sv,lda,W,WORK,LWORK,INFO)
!        print *, INFO
!        print *, work(1)
!         lwork = int(work(1))
!         deallocate(work)
!         allocate(work(lwork))
!         sv(:,:) = mbc % s(:,:)
!         call dsyev(num,up,n,sv,lda,W,WORK,LWORK,INFO)
!         print *, info
!         print *,'eigen'
!         print *, W
!        deallocate(work)
! 
!        allocate(sd(mbc%thsz,mbc%thsz))
!       
!        
!        sd = 0.0_dp
!        
!        do i = 1,mbc%thsz
!           sd(i,i) = sqrt(W(i))
!        enddo
!        
!        do i = 1,mbc%thsz
!           do j = 1, mbc%thsz
!               mbc % hs(i,j) = sum(mbc % s(i,:) * mbc % h(:,j))
!           enddo
!        enddo
!        
!        deallocate(sd)
!        a(:,:), ac(:,:), aca(:,:), v(:), w(:,:), occ(:)

      if (mbc%mxgf) then
        if (.not. allocated(mbc%a)) allocate(mbc % a(mbc % bhsz, mbc % thsz))
        if (.not. allocated(mbc%ac)) allocate(mbc% ac(mbc % bhsz, mbc % thsz))
      else
        if (.not. allocated(mbc%cn)) allocate(mbc % cn(mbc % thsz, mbc % thsz,nsp))
        if (.not. allocated(mbc%en)) allocate(mbc % en(mbc % thsz,nsp))
      endif
!       if (.not. allocated(mbc%v)) allocate(mbc % v(mbc % thsz))
!       if (.not. allocated(mbc%occ)) allocate(mbc % occ(mbc % thsz))
!       if (.not. allocated(mbc%w)) allocate(mbc % w(mbc % thsz, mbc % thsz))
      if (.not. allocated(mbc%gt)) allocate(mbc % gt(mbc % thsz, mbc % thsz))
      if (.not. allocated(mbc%ggd)) allocate(mbc%ggd(mbc % ahsz))
      if (.not. allocated(mbc%q) .and. .not. mbc%mxgf) allocate(mbc%q(nsp,mbc % n(1)))

!       stop

      

   end subroutine mbldh


   subroutine mpos(bgn, fin, p,mbc)
      type(mbtb_t), intent(inout) :: mbc
      integer, intent(in) :: bgn, fin
      integer, intent(out) :: p(:)
      integer :: ia

!     calculate offsets in matrices using the common nstt. 
!     To embed atoms of the same type but with different number of orbitals it may be better if the conf_t structures are used 

      include "../Include/Atom.array"

      p(1) = 0
      if (bgn < mbc %rlst(2)) then
          do ia = 2, fin - bgn + 2 ! go one extra 
                p(ia) = p(ia-1) + nstt(z( bgn  - 2 + ia)) ! gia = off + ia; gia: global ia
          enddo
      else
          do ia = 2, fin - bgn + 2
              if (mbc%orth(ia-1) == 0) then
                  p(ia) = p(ia-1) + mbc%tbc%nstt(z( bgn  - 2 + ia)) ! gia = off + ia; gia: global ia
              else
                  p(ia) = p(ia-1) + nstt(z( bgn  - 2 + ia))
              endif
          enddo
      endif
      
   end subroutine mpos
   
   subroutine bldU(ctsx,n,off)
      
      use mod_conf, only: rtc
      use mod_atom_ar, only: asort,zatype
      
      type(ctsx_t), intent(inout) :: ctsx
      integer :: it,git,n,off

!     Calculate Hubbard U for sctb case, just to establish a scale of distance from the 
!     BOP region for each TB atom

      include "../Include/Atom.array"
!       include "../Include/PosVel.array"
      
!       dumx = 0.0_dp
!       dumin = 100000_dp
!       du(:) = 0.0_dp
!       do it = 1,mbc%n(2)
!           git = it + mbc % rlst(2) -1
!           do ib = 1,mbc%n(1)     
!               dr = ad(:,git)-ad(:,ib)     
!               magdr = (dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3))
!               if (magdr < 36) du(it) = du(it) + (1/magdr)      
!           enddo
!           if (du(it) > dumx) dumx = du(it)
!           if (du(it) > 10d-8 .and. du(it) < dumin) dumin = du(it)
!       enddo
      
!       print *, 'HUBBARD U'
!       ,dumx
      do it = 1, n
! !           if (mbc%orth(it) == 1) then
! !               mbc%u(it) = 27.211386
!           if (du(it) > 10d-8) then
!               mbc%u(it) = 13.605693 + 1.5*13.605693*((du(it)-dumin)/(dumx-dumin))
!               if (1.5*13.605693*((du(it)-dumin)/(dumx-dumin)) > 13.605693*2) &
!                           &  mbc%u(it) = 13.605693 + 13.605693
!           else 
             git = it + off -1
              ctsx%u(it) = rtc % tbham_conf % a(asort(zatype(z(git)))) % uh
!           endif
!           print *, it, ctsx% u(it)
!           ,du(it), mbc % u(it)
      enddo
      
   end subroutine bldU


   subroutine bldsparsedh(ap, bp, r1, r2, lst, h,mbc)
   
!    build sparse dH using the neighbour table
      type(mbtb_t), intent(inout) :: mbc
      integer, intent(in) :: ap(:), bp(:), r1, r2, lst(:)
      real(dp), intent(out) :: h(:,:,:)

      integer :: i,j,off1,off2, ia
      integer :: za,zb
      integer :: ja,ja0,ib,la,lb,gia,gib
      integer :: nstta,nsttb,n

      real(dp) :: dr(3), dea,magdr,dumx
      real(dp) :: subh(mxnstat,mxnstat)
!       integer  :: llist1(ml),nl1
!       integer  :: llist2(ml),nl2
      real(dp) :: scfcut(14)
      
      
      type(ham_conf_t),pointer :: ham
      
      include "../Include/Atom.array" ! for z mainly
      include "../Include/PosVel.array"

      
      dea = 0.0_dp
      scfcut = 1.0_dp
      h(:,:,:) = 0.0_dp
 
      off1 = lst(r1)-1
      off2 = lst(r2)-1

      ham => rtc % ham_conf 
      
      do gia = lst(r1), lst(r1+1)-1   

         za = z(gia)
         nstta = nstt(za)
!          nl1 = mbc%tbc%nl(za)
!          llist1=mbc%tbc%llist(:,za)
 
         ja = ap(gia-off1) ! ia = gia-off1
         
         if (ja == 0) then
            ja=1
         endif
         ja0 = ja - 1

         do while (bp(ja) /= eol)           
!             print *,ja
!             print *,bp
            ib = bp(ja)
            gib = off2 + ib 

            zb = z(gib)
            nsttb = nstt(zb)
             
            dr = ad(:,gia)-ad(:,gib)     
            
         
!             print *, "ad"
!             print *, gia,ad(:,gia),gib,ad(:,gib)
!             call states(zb,nl2,nsttb,llist2)
!               print *, 'sparse: gia, gib, ja:', gia, gib, ja


            
            call matel(za,nstta,zb,nsttb,dr,subh,dea,scfcut,ham,.false.)     
!             print *, subh
!             if (gia == 41 .and. ib == 40) then
! 	        print *, gia, gib
! 	      print *, subh
! 	    endif
!             subh(1,:) = 0.0_dp
!             n=0
!             if (za == zb .and. mbc%orth(gia-off1) == 0) then
!                 do i = 1, mbc%tbc%nl(za)
!                     if (mbc%tbc%llist(i,za) /= llist(1,za)) then
!                         n = n + 2*mbc%tbc%llist(i,za)+1
!                     else
!                         exit
!                     endif    
!                 enddo
!             endif
        

            h(:,:,ja) = subh(:,:)
!             print *, gia,gib,ja
!             print *,dr
!             do i = 1, nstta
!                 print *, h(i,:,ja)
!             enddo
            

            ja = ja+1
!             print *, gia,gib
!             print *, dr

              
         end do
         
      end do
      
   end subroutine bldsparsedh



   subroutine bldmtbmat(bgn, fin, ph, ap, bp, h,mbc,ovl)
   
!       Build the full size embedded molecular TB Hamiltonian matrix
      type(mbtb_t), intent(inout) :: mbc
      integer, intent(in) :: bgn, fin, ph(:), ap(:), bp(:)
      real(dp), intent(out) :: h(:,:)
  
      integer :: i,j,off, iap, ibp
      integer :: za,zb
      integer :: ja,ja0,ia,ib,la,lb,gia,gib
      integer :: nstta,nsttb
      integer :: llist1(ml),nl1
      integer :: llist2(ml),nl2,off1,off2

      real(dp) :: dr(3), dea
      real(dp) :: subh(mxnstat,mxnstat)
      type(ham_conf_t),pointer :: hamtb
      logical :: ovl

      real(dp) :: scfcut(14)

      include "../Include/Atom.array" ! for z mainly 
      include "../Include/PosVel.array"


      dea = 0.0_dp
      scfcut = 1.0_dp

      h(:,:) = 0.0_dp

      off = bgn - 1
      hamtb => rtc % tbham_conf 

      do ia = 1, fin - off      
         gia = off + ia

         za = z(gia)
         
!          nl1 = mbc%tbc%nl(za)
!          llist1=mbc%tbc%llist(:,za)

         ja = ap(ia)

         ja0 = ja - 1
         do while (bp(ja) /= eol)
            nstta = mbc%tbc%nstt(za)
            ib = bp(ja)

            gib = off + ib
            zb = z(gib)
            nsttb = mbc%tbc%nstt(zb)
!             nl2 = mbc%tbc%nl(zb)
!             llist2=mbc%tbc%llist(:,zb)
            dr = ad(:,gia)-ad(:,gib)

!             print *, "ad"
!             print *, gia,ad(:,gia),gib,ad(:,gib)
             
             
! !              % b(bsort(btype(za,zb)))
            call matel(za,nstta,zb,nsttb,dr,subh,dea,scfcut,hamtb,ovl)
                       
            off1 = 0
            off2 = 0
            if (mbc%orth(ia) + mbc%orth(ib) /= 0) then
                nstta = nstt(za)
                nsttb = nstt(zb)
                off1 = mbc%tbc%nstt(za) - nstta
                off2 = mbc%tbc%nstt(zb) - nsttb
                if (zc(za) == 0 .or. zc(zb) == 0) then
                    ja = ja+1
                    cycle
                endif
            endif

            iap = ph(ia)
            ibp = ph(ib)
            h((ph(ia+1)-nstta)+1:ph(ia+1), (ph(ib+1)-nsttb)+1:ph(ib+1)) =  &
                &         subh(1+off1:mbc%tbc%nstt(za),1+off2:mbc%tbc%nstt(zb))
            
!             if (mbc%orth(ia)+mbc%orth(ib) < 2) then
!                 h(ph(ia)+1+mbc%orth(ib):ph(ia+1), ph(ib)+1+mbc%orth(ia):ph(ib+1)) =  &
!                 &         subh(1+mbc%orth(ia)+mbc%orth(ib):nstta,1+mbc%orth(ib)+mbc%orth(ia):nsttb)
!             else
!                 h(ph(ia)+1:ph(ia+1), ph(ib)+1:ph(ib+1)) = subh(1+mbc%orth(ia):nstta,1+mbc%orth(ib):nsttb)
!             endif
              
!             print *, 'ia,iao,ib,ibo,ovl',ia,mbc%orth(ia),ib,mbc%orth(ib), ovl
!             
! !             print *, 
! !             do i = 1,nstta
! !                 print *,  subh(i, :)
! !             enddo
! !             print *, ' '
!             PRINT *,ia,ib
!             do i = ph(ia)+1,ph(ia+1)
!                 print *,  h(i, ph(ib)+1:ph(ib+1))
!             enddo
!             print *, ph(ia)+1,ph(ia+1),ph(ib)+1,ph(ib+1)
            ja = ja+1
         end do
      end do

   end subroutine bldmtbmat


   subroutine swap(n,a,b)
!    utility routine to test different swap algorithms

      integer, intent(in) :: n
      integer, intent(inout) :: a(n), b(n)

      integer :: i, t

      do i = 1, n
         t = a(i)
         a(i) = b(i)
         b(i) = t
      end do

!      do i = 1, n
!         a(i) = ieor(a(i),b(i))
!         b(i) = ieor(a(i),b(i))
!         a(i) = ieor(a(i),b(i))
!      end do


   end subroutine swap




   subroutine extract_tblist(api, bpi, r1, r2, lst, apo, bpo,mixgf)
         use mod_conf
!    extracts the output position tables apo, bpo from the input ones api, bpi, for the atoms listed in lst,
!       such that the first atom of a bond is bound to be in region r1 and the second one in r2
      include "../Include/NebList.array"
      include "../Include/PosVel.array"
      include "../Include/Atom.array"
      
      integer, intent(in) :: api(:), bpi(:), r1, r2, lst(:)
      integer, intent(out) :: apo(:), bpo(:)

      integer :: ia, ib, ja, jb, p, mp
      logical :: mixgf


      mp = 1
      do ia = lst(r1), lst(r1+1) - 1
         ja = ia - lst(r1) + 1
         apo(ja) = mp
         
!           CHECK THIS WITH HYDROGEN
!  To prevent non-BOP atoms like H from getting BOP neighbours
!          if (r2 == 1 .and. zc(z(ia)) == 0) then
!             bpo(mp) = eol
!             mp = mp +1
!             cycle
! 
!          endif 


         jb = api(ia)
         p = 0
         ib = bpi(jb+p)
         if (.not. (r2 == 1 .and. zc(z(ia)) == 0)) then
            do while (ib /= eol)
                if (rmap(ib, lst) == r2) then
                    bpo(mp) = ib - lst(r2) + 1
                    mp = mp + 1
                elseif(.not. mixgf .and. r2 ==1 .and. map(ib) /= ib) then
                    bpo(mp) = ib - lst(r2) + 1
                    mp = mp + 1
!                 elseif (rmap(ib, lst) == 3 .and. sqrt((ad(1,ia)-ad(1,ib))**2 + &
!                         & (ad(2,ia)-ad(2,ib))**2 + (ad(3,ia)-ad(3,ib))**2) < 4) then
!                         print *, "TB atom ",ia," has neighbour ",ib," at range"
!                         print *, sqrt((ad(1,ia)-ad(1,ib))**2 + (ad(2,ia)-ad(2,ib))**2 + (ad(3,ia)-ad(3,ib))**2)
                end if
                p = p + 1
                ib = bpi(jb+p)
                
            end do
         endif
         bpo(mp) = eol

         if (bpo(max(mp-1,1)) /= eol .or. mp == 1) then
            mp = mp + 1
         else
            apo(ja) = apo(ja)-1
         end if
      end do

      apo(lst(r1+1)- lst(r1) + 1) = mp-1

      
!       if (r2 == 1) then
!           print *, 'apo:', apo
!           print *, 'bpo:', bpo
!       endif
!     stop
   end subroutine extract_tblist

   
   subroutine extract_boplist(api, bpi,lst, apo, bpo,natm,mixgf)
         use mod_conf
         use mod_const
!    extracts the output position tables apo, bpo from the input ones api, bpi, for BOP atoms, 
!   is made so that it can also include the neighbours for atoms in PBC cells for bldclus
      include "../Include/NebList.array"
      include "../Include/PosVel.array"
      include "../Include/Atom.array"
      
      integer, intent(in) :: api(:), bpi(:),lst(:)
      integer, intent(out) :: apo(:), bpo(:)

      integer :: ia, ib, ja, jb, p, mp,natm
      logical :: mixgf

      mp = 1
      do ia = 1,natm
         if ((map(ia) == ia .and. rmap(ia, lst) == 1) &
             & .or. (map(ia) /= ia) .or. ( .not. mixgf)) then
            apo(ia) = mp
            jb = api(ia)
            p = 0
            if (zc(z(ia)) /= 0) then
                ib = bpi(jb)
                if (ib == 0) exit
                do while (ib /= eol)       
    !          Need to test this with multiple atom types
                    if ((map(ib) == ib .and. rmap(ib, lst) == 1) &
                    &      .or. (map(ib) /= ib .and. zc(z(ib)) /= 0) &
                       .or. (.not. mixgf .and. zc(z(ib)) /= 0)) then
                    bpo(mp) = ib 
                    mp = mp + 1
                    end if
                    p = p + 1
                    ib = bpi(jb+p)
                end do
            endif
            bpo(mp) = eol
            if (bpo(max(mp-1,1)) /= eol) then
                mp = mp + 1
            else            
               apo(ia) = apo(ia)-1
            end if
         endif
         if (jb == 0) exit
      end do
    
      apo(ia) = mp-1

!       print *, 'apo:', apo
!       print *, 'bpo:', bpo
! 

   end subroutine extract_boplist


   subroutine print_neblist(s,e,a,b)
!    print the position tables ap, bp from starting atom s to endind atom e
      include "../Include/NebList.array"

      integer, intent(in) :: s,e, a(:), b(:)

      integer :: ia, ja, p, ib, gib

!       print *, s,e
      do ia = 1, e
!         print *, 'ia:', ia
         ja = a(ia)
         p = 0
         ib = b(ja+p)

         do while(ib/=eol)
            gib= ib+s-1
!            print *, '   p, ib, map(ib)', p, ib, map(ib)
            p = p + 1
            ib = b(ja+p)
         end do
      end do

   end subroutine print_neblist
   
   
   subroutine loadtb( tbham_conf,tbem,ovl )
         
      use mod_precision
      use mod_all_scalar
      use mod_const
      use mod_conf
      use mod_pft

      implicit none

      type(tbconf), intent(inout) :: tbem  
      type(ham_conf_t), intent(in), target :: tbham_conf
      type(bond_conf_t), pointer :: b
      type(hop_t), pointer :: h,o
      integer :: zin, i, j,  li, lj, oi, oj, bi
      logical :: ovl

      
      
      do j = 0,natype
         zin =  tbham_conf % a(j) % z
         
         tbem%zc(zin) =  tbham_conf % a(j) % b % cc 
         tbem%nstt(zin) = 0
         
         if (mag) tbem%istn(zin) =  tbham_conf % a(j) % b % st
         
         tbem%nl(zin) =  tbham_conf % a(j) % b % norbs
         tbem%llist(:,zin) =  tbham_conf % a(j) % b % orblist
         
         do i = 1, tbham_conf % a(j) % b % norbs
            tbem%nstt(zin) = tbem%nstt(zin) + 2*(tbham_conf % a(j) % b % orblist(i))+1
         enddo
         
         tbem%es(zin) =  tbham_conf % a(j) % b % es
         tbem%ep(zin) =  tbham_conf % a(j) % b % ep
         tbem%ed(zin) =  tbham_conf % a(j) % b % ed
      end do
      
      do i = 0, tbham_conf % nb - 1
         b => tbham_conf % b(i)
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
      
   end subroutine loadtb


   function rmap(ia, rlst)
!       range map, return the range in which atom ia resides according to the range's list rlst

      integer, intent(in) :: ia, rlst(:)
      integer :: rmap
!       To be generalised in the future. For now there are only 2 blocks

      if (ia < rlst(2)) then
         rmap = 1
      else if (ia >= rlst(2) .and. ia < rlst(3)) then
         rmap = 2
      else
         rmap = 3
      end if

   end function rmap
   
   
      subroutine sylmnc(c,lmx)
! C- Normalization constants for the spherical harmonics
! C ----------------------------------------------------------------
! Ci Inputs
! Ci   lmx
! Co Outputs
! Co   C
! Cr Remarks
! Cr   use together with sylm (from ASW package)
! C ----------------------------------------------------------------
      implicit none
! C ... Passed parameters
      integer  :: lmx
      real(dp) :: c(*)
! C ... Local parameters
      integer  :: i,i1,i2,l,lav,lp1,m,n1,n2,n3
      real(dp) :: fn2,fpi,tlp1,tpi,y0

      tpi = 8d0*datan(1d0)
      fpi = 2d0*tpi
      y0 = 1d0/dsqrt(fpi)
      c(1) = y0
      do  l = 1, lmx
        lp1 = l+1
        tlp1 = l+lp1
        lav = l*lp1 + 1
        c(lav) = dsqrt(tlp1/fpi)
        do  m = 1, l
          n2 = lp1-m
          n1 = n2+1
          n3 = l+m
          fn2 = n2
          do  i = n1, n3
            fn2 = fn2*i
          enddo
          i1 = lav+m
          i2 = lav-m
          c(i1) = dsqrt(tlp1/(fn2*tpi))
          c(i2) = c(i1)
        enddo
      enddo
   end subroutine sylmnc
   

   end module tbbop_emb
