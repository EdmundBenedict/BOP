
   subroutine bldh()
      use mod_precision
      use mod_all_scalar
      use mod_const
      use mod_conf
      use mod_ham
      use mod_atom_ar
!       , only : btype, bsort
      
!    This is a subroutine to set up a tight binding Hamiltonian matrix.
!

      implicit none


      include "Include/Atom.array"
      include "Include/NebList.array"
      include "Include/PosVel.array"


      real(dp) :: dea

      integer za,zb
      integer ja,ja0,ia,ib
      integer ic,i,j,p
      integer llist1(ml)
      integer llist2(ml)
      integer nstt1,nstt2,nl1,nl2
   !*** Displacement.
      real(dp) :: dr(3)
      
      real(dp) :: scfcut(14)
      
      type(ham_conf_t),pointer :: hambp

      scfcut = 1.0_dp
      
   !
   !    Set up the Hamiltonian.
   !
   !    For each ion (A) do ...
   !

      hambp => rtc % ham_conf 
      
      do ic = 1, pcluster(nbase+1)-1

         ia = cluster(ic)
         za = z(ia)
         nstt1 = nstt(za)
!          call states(za,nl1,nstt1,llist1)
         
         if (map(ia) /= 0) then
            dea = de(map(ia))

         else
            dea = 0.0_dp
         endif
          
         ja0 = aptr(ia) - 1
         p = 1
         ib = bptr(p+ja0)          
         do while (ib /= eol)
            zb = z(ib)
            nstt2 = nstt(zb)
            dr = ad(1:3,ia)-ad(1:3,ib)
               
!             call states(zb,nl2,nstt2,llist2)
            
            
            call matel(za,nstt1,zb,nstt2,dr,h(1,1,p,ic),dea,scfcut,hambp,.false.)

            p = p + 1
            ib = bptr(p+ja0)
         enddo

      enddo
      
      nullify(hambp)

   end subroutine bldh

