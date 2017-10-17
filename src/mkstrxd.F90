
subroutine mkstrxd(force,nbas,off, nlmq, nlmq1, &
                  & indxcg, jcg, cg, cy, struxd, dstrxd, struxidx, dstrxidx, struxsize, dstrxsize)
!- Make lookup table of strux and radial derivatives
! ----------------------------------------------------------------------
!i Inputs:
!  force   : logical to include higher multipoles for F calc.
!  off     : Offset for TB atoms (emb)
! nlmq and nlmq1 are currently set to 1, could be increased for higher multipoles
!i  ldip  : 3 include 'spherical' dipole correction to Ewald
!i        : 2 include 'slab' dipole correction to Ewald
!i        : any other number - skip the dipole correction
!i   nbas,bas,awld,alat,vol,dlat,nkd,glat,nkg,indxcg,jcg,cg,cy
!i   nlmq1: leading dimension of strx, nlmq1=(ll(nlmq)+1)**2
!i   nlmq : max L-cutoff for multipoles, leading dimensions of dstrx,
!i           second dimension of strx
!i   lmxl : l-cutoff for multipoles for each class of species
!i   pv   : true if radial derivatives required
!i   MOL  : true for molecule (cluster) branch
!i   struxidx    : index for struxd, shall be allocated to (nbas,tbc%escount(pid)),
!i                 tbc%escount(pid) the number of atoms allocated to process 'pid'
!i   dstrxidx    : index for dstrxd
!
!o Outputs:
!o  struxd : coefficients B_LL'(R'-R) (structure constants)
!o  dstrxd : derivatives of structure constants (x dB_LL'(x)/dx) at x = |R'-R|
!o          if pv = F dstrx is not touched
!r Remarks
!r   Instead of making B on each sc iteration, the program prepares the whole B
!r   in the beginning of the self-consistency cycle to be used as a lookup table.
!r   This results in faster calculation at an expence of memory.
!r
!r   Calling mkstrx is set as the default option. To call instead hstr on each
!r   iteration (and therefore to reduce memory but increase computational time),
!r   use switch --sfly
!r
!r   Efficiency issues: there are two symmetry properties of zero energy
!r   structure constants that are used in mkstrx:
!r     B_LL'(R'-R) = B_L'L(R-R')                (1)
!r     B_LL'(R'-R) = (-1)^(l+l') B_L'L(R'-R)    (2)
!r   same properties hold for the radial derivative of B.
!r
!r This symmetry however is not employed because the complexity that a
!r  parallel implementation will require and the diminishing gains + increased overhead with.
!r Curently the matrices and the indices are distributed in contiguous blocks across the long
!r side and no communication whatsoever is needed. For a process with id pid, the block starts
!r at atom tbc%esamap(pid)+1 and ends at tbc%esamap(pid+1) inclusive. For atom ib from this range
!r there is a column of neighbours in struxidx(:,ib-tbc%esamap(pid)). For each of pair (jb,ib)
!r a block of size (nlmj,nlmi1) is preallocated starting at struxd(struxidx(jb,ib-tbc%esamap(pid)))
!r The block size for dstrx is (nlmj,nlmi)
!r
!b Bugs
!b   mkstrx is not written in the most efficient way and could be further refined wrt to
!b   amount of calculated elements of strux. At the moment mkstrx is a result of certain
!b   trade off between performance and clarity of the program.
!b
!u Updates
!u           2013 (DP)  strx moved to indexed linear array with auxiliary index for better compatibility with C
!u           2013 (DP)  Fortran pointer based strx for memory efficiency
!u           2012 (DP)  Parallell distributed strx & dstrx
!u    05 Mar 2010 (SL)  optimization (make strux up to Lmax and for triangle
!u                      jb<=ib only)
!u    19 Jan 2010 (ATP) first written
! ----------------------------------------------------------------------
 
      use mod_atom_ar
      use mod_precision
      use mod_all_scalar
      use mod_const
 
!   need to check how essential these are
!    use tbprl
!    use structures
   implicit none

   integer, intent(in) ::  off, indxcg(7400),jcg(62200),struxsize, dstrxsize
   real(dp), intent(in) :: cy(289),cg(62200)
   integer, intent(in) :: struxidx(nbas,*), dstrxidx(nbas,*),nbas
   real(dp), intent(out) :: struxd(struxsize), dstrxd(dstrxsize)


   integer :: ib,gib,jb,gjb,lmax,lmxf,lmxst,nlmj,nlmi,nlmi1,nlmq,nlmq1
   integer :: li,li1,lj, ilm,jlm, i, cmol, cpv,  u, pi1, pj, sz, dsz
   real(dp) :: dr(3), hl(100),bl(100),struxl(nlmq*nlmq1), dstrxl(nlmq*nlmq)
!    , plat(3,3), qlat(3,3), alat,vol, awald, det, 
   logical :: force
!    procedure(logical) :: bittst
!    integer :: pid
!    logical :: lroot
   
      include "Include/PosVel.array"
      include "Include/NebList.array"    
      include "Include/Atom.array"
      
!    call tcn('mkstrxp')

!    lroot  = tbc % c3d % lrt
!    pid    = tbc % c3d % id
! 
!    pv  = bittst(s_ctrl%ltb,128)
!    mol = bittst(s_ctrl%ltb,2**18)
!    force = bittst(s_ctrl%ltb,16)

!    plat  = s_lat%plat
!    alat  = s_lat%alat
!    vol   = s_lat%vol
!    awald = s_lat%awald


! is this required for the real space calc?
!    call dinv33(s_lat%plat,0,qlat,det)

!    if (.not. pv) 
   pj = 1

   do  ib = 1, nbas
!        print *, 'ATOM 1 = ', ib
       gib = ib+off
!    increase li to include higher multipoles
      li = 0
      nlmi  = (li+1)*(li+1)

      if (force) then
         li1 = li + 1
         nlmi1  = (li1+1)*(li1+1)
      else
         li1=li
         nlmi1  = nlmi
      endif

      do  jb = 1, nbas
!          print *, 'ATOM 2 = ', jb
         gjb = jb+off
         lj = 0
         nlmj = (lj+1)*(lj+1)

         lmxst = li1 + lj
!          print *, 'lmxst',lmxst
!          tau = s_lat%pos(1:3,jb)-s_lat%pos(1:3,ib)
         dr = (ad(1:3,gjb) - ad(1:3,gib))
!          /0.529177 to compare to TBE values
!          print *, gib,gjb,dr   
!          if (.not. mol) then
!             call directshortn(tau,plat,qlat)
!             call rcnsl0(tau,awald,lmxst,alat,glat,nkg,dlat,nkd,vol,cy,hl)
!          else
!             taua = tau * alat

            call soldhj(dr,0d0,0,lmxst,hl,bl)
!          endif

         call hstrd(struxl,dstrxl, nlmi1, nlmi, nlmj, hl, cg,indxcg,jcg,vol)

         pi1 = struxidx(jb,ib)
         sz = nlmj*nlmi1
         struxd(pi1:pi1+sz-1) = struxl(1:sz)
!          if (pv) then
!             pj = dstrxidx(jb,pib)
!             dsz = nlmj*nlmi
!             dstrxd(pj:pj+dsz-1) = dstrxl(1:dsz)
!          end if
!                  print *, struxd(pi1:pi1+sz-1)
!                   print *, ' '
      enddo
   end do
!    stop
!    call tcx('mkstrxp')


end subroutine mkstrxd


subroutine hstrd(strux,dstrx,nlmi1,nlmi,nlmj,hl,cg,indxcg,jcg,vol)
! C- Make structure constants from reduced strux at energy zero
! C ----------------------------------------------------------------------
! Ci Inputs:
! Ci  nlmf,nlm :make B_LL' for L = 1, nlmf, L'= 1, nlm
! Ci  nlmq1 :leading dimensions of B, nlmq1=(ll(nlmq)+1)**2
! Ci  nlmq  :max L-cutoff for multipoles, leading dimension of B'
! Ci  hl    :radial Hankels for given point |tau|
! Ci  cg,indxcg,jcg : Clebsch Gordan coefficients in condensed form
! Ci  vol   :unit cell volume
! Co Outputs: strx,drstrx
! Co  strx  :coefficients B_LL'(tau) (structure constants)
! Co  dstrx :derivatives of structure constants (x dB_LL'/dx) at x = tau
! Co         if pv = F dstrx is not touched
! Cr Remarks: B_LL'(tau) = 4\pi \sum_L" (-1)^l' H_l"(tau) Y_L"(tau)
! Cr         HY are reduced structure constants (see rcnsl0, soldh)
! Cr         If pv is set, returns tau*d/dtau B in drstrx
! C ----------------------------------------------------------------------
   implicit none
   integer nlmi1,nlmi,nlmj,indxcg(*),jcg(*),ldip
   real(8) :: cg(*),strux(nlmj,nlmi1),dstrx(nlmj,nlmi),hl(*),vol
   integer lmxx,icg,icg1,icg2,ii,ilm,indx,ipow,klm,l,lk,ll,llm,lm,lp,lmax,mlm
   parameter (lmxx=12)
   real(8) :: sig(0:lmxx),fourpi,fpibv,sum

   ! !       call tcn('hstr: make strux')
   fourpi = 16d0*datan(1d0)
   lmax = ll(nlmi1) + ll(nlmj)
   if (lmax .gt. lmxx) stop ' change dimensions in hstr'
   ! C --- (-1)^l ---
   sig(0) = 1d0
   if (lmax .gt. 0) then
      do  l = 1, lmax
         sig(l) = -sig(l-1)
      enddo
   endif
   ! C --- add together Gaunt-sums ---
   do  mlm = 1, nlmi1
      lm = ll(mlm)
      do  klm = 1, nlmj
         lk = ll(klm)
         sum = 0d0
!          sumr = 0d0
         ii = max0(mlm,klm)
         indx = (ii*(ii-1))/2+min0(mlm,klm)
         icg1 = indxcg(indx)
         icg2 = indxcg(indx+1)-1
         do  icg = icg1, icg2
            llm = jcg(icg)
            lp = ll(llm)
            ipow = (lm + lk - lp)/2
            if (ipow .eq. 0) then
                sum = sum + cg(icg)*hl(llm)
!                 if (pv) then
!                   sumr = sumr - (lp+1) * cg(icg)*hl(llm)
!                 endif
            endif
         enddo
         strux(klm,mlm) = sum*fourpi*sig(lk)
   ! c need to cut off l+1 terms for B'. Quick fix for now
!          if (pv .and. mlm <= nlmi) dstrx(klm,mlm) = sumr*fourpi*sig(lk)
      enddo
   enddo
   ! C --- the following includes extra p terms 'implicitly' ---
!    if ((ldip .eq. 2 .or. ldip .eq. 3) .and. (.not. MOL)) then
!       if (min(nlmj,nlmi1) > 1) then
!          fpibv = fourpi/vol
!          do  ilm = 2, min(4,nlmj,nlmi1)
!             strux(ilm,ilm) = strux(ilm,ilm) - fpibv
! !             if (pv .and. ilm <= nlmi) dstrx(ilm,ilm) = dstrx(ilm,ilm) + 3d0*fpibv
!          enddo
!       endif
!    endif

   call tbshfl(0,nlmj,nlmj,nlmi1,strux)
   call strfac(0,nlmj,nlmj,nlmi1,strux)

!    if (pv) then
!       call tbshfl(0,nlmj,nlmj,nlmi,dstrx)
!       call strfac(0,nlmj,nlmj,nlmi,dstrx)
!    end if

!       call tcx('hstr: make strux')
end subroutine hstrd

   
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
   