      subroutine bldhdo(isp)
      use mod_precision
      use mod_all_scalar
      use mod_ham
      use mod_const
      use mod_atom_ar, only : demi

!
!    This is a subroutine to set up a tight binding Hamiltonian matrix.
!

      implicit none

      integer, intent(in) :: isp

      
      include "Include/Atom.array"
      include "Include/NebList.array"

      real(dp) :: esa,epa,eda, dea, epde, edde
      real(dp), pointer :: subh(:,:)

      integer :: za
      integer :: ja,ja0,ia,ib
      integer :: nstta,mxh,nla
      integer :: nclus,ic,i,j1, i1,i2, mia
      
      
!       integer :: cnt
!       cnt =0
      nclus = pcluster(nbase+1)-1
      mxh = nclus*mxnstat
      
      do ic = 1, nclus

         ia = cluster(ic)

         za = z(ia)
         mia = map(ia)
         
         if (mia /= 0) then
            dea = de(mia)
            if (mag) dea = dea + oppm(isp)*dem(mia)
         else
            dea = 0.0_dp
!             print *, 'bldhdo:, ic, ia, mapi', ia, ia, mapi(ia)
            if (mag) dea = dea + oppm(isp)*demi(mapi(ia))
         endif
         
       
         
         
         esa = es(za)
         epa = ep(za)
         eda = ed(za)
         nla = nl(za)
         nstta = nstt(za)
         llista = llist(:,za)
         epde = epa+dea
         edde = eda+dea

         ja = aptr(ia)
         ja0 = ja - 1

         do while (bptr(ja) /= ia)
            ja = ja+1
         end do
   
         subh => h(1:mxnstat,1:mxnstat,ja-ja0,ic)

         subh = 0.0_dp
         
         j1 = 1
         do i1 = 1,nla
            select case(llista(i1))
            case (2)
               subh(j1,j1) = edde
               subh(j1+1,j1+1) = edde
               subh(j1+2,j1+2) = edde
               subh(j1+3,j1+3) = edde
               subh(j1+4,j1+4) = edde
               j1 = j1 + 5
            case (1)
               subh(j1,j1) = epde
               subh(j1+1,j1+1) = epde
               subh(j1+2,j1+2) = epde
               j1 = j1 + 3
            case (0)
               subh(j1,j1) = esa+dea
               j1 = j1 + 1
            end select
         enddo
      
      end do

!       print *,'cnt',cnt,nclus
      end subroutine bldhdo
