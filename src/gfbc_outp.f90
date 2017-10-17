 
      subroutine gfbc_outp()
          use mod_precision


          use mod_all_scalar

          use mod_const
!
!     This subroutine writes the output file required by greenlat.c
!     to apply the elastic and lattice Greens funtions and update
!     the positions of atoms. It is called at the end of relax_ds.f
!     by the GFBCON switch in fort.8.
!
!     At the moment this is for one element only. I guess I'll fix this
!     at some point in the future. 
!
!     M. Cawkwell. 11th September 2003
!
      implicit none
!
!
!    Include constants.
!

!      include "Include/ALL.const"

!
!    Include scalars.
!

!      include "Include/ALL.scalar"

!
!    Include arrays.
!

!      include "Include/Atom.array"
      include "Include/BEC.array"
!      include "Include/BondOrder.array"
      include "Include/Force.array"
!      include "Include/Hamilt.array"
!      include "Include/Misc.array"
!      include "Include/Moment.array"
!      include "Include/NebList.array"
!      include "Include/NotSRT.array"
      include "Include/PosVel.array"
!      include "Include/RecurCoef.array"
!      include "Include/Relax.array"
!      include "Include/SRT.array"
!
      real(dp) :: xdimx, ydimy, zdim1, zdim2, zero
      integer flag, i, one
!
      xdimx = 1.0e4_dp
      ydimy = 1.0e4_dp
      zdim1 = lena(3)
      zdim2 = 0.0_dp
      one = 1
      zero = 0.0_dp

!
!     We need the forces acting on the active atoms. The forces on the 
!     inert atoms will automatically be set to zero.
!     
      flag = 1
      call getetot(flag)
      call erasab()
!
      open (unit=96, status="UNKNOWN", file="reloutput.dat")
!
!     Format of file is as follows
!     line 1 = title
!     line 2 = number of atoms (active and inert)
!     lines 3 and 4 = limits of cell +x, +y, +z; -x, -y, -z
!     line 5 = atomic mass of element (???)
!     REPEATED OVER ALL ATOMS
!     line 6 = x,y,z coordinate of atom in A
!     line 7 = fx, fy, fz in eV / A
!     line 8 = species (1, 2, 3...)
!
      write(96,'("Coords and forces from dislocation relaxation")')
      write(96,20) nd+ninert, one
      write(96,10) xdimx, ydimy, zdim1
      write(96,10) -1.0_dp*xdimx, -1.0_dp*ydimy, zdim2
      write(96,*) "    6.084255408495665E-03        77"
      do i = 1,nd
         write(96,10) ad(1,i), ad(2,i), ad(3,i)
         write(96,10) ftot(1,i), ftot(2,i), ftot(3,i)
         write(96,30) one
      enddo
      do i = 1,ninert
         write(96,10) adinert(1,i), adinert(2,i), adinert(3,i)
         write(96,10) zero, zero, zero
         write(96,30) one
      enddo
!
 20   format(6x,i4,9x,i1)
 10   format(3x,e22.16,3x,e22.16,3x,f22.16)
 30   format(9x,i1)
!      
      return
      end
