 
       function repeng()
          use mod_precision


          use mod_all_scalar

          use mod_const
          use mod_atom_ar
!
!    This is a routine to evaluate the repulsive part of the energy.
!

      implicit none
      real(dp) :: repeng

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

      include "Include/Atom.array"
      include "Include/NebList.array"
      include "Include/PosVel.array"

!
!    Declare the simple variables.
!

      real(dp) :: phi,vscale
      real(dp) :: sum
      real(dp) :: r

      integer i1,i2,i2a
      integer z1,z2,ptr

!
!    Declare the local arrays.
!

!*** The displacement between two atoms in different cells.
      real(dp) :: dr(3)
!
!      EXTERNAL SCALE
!
!
!    Evaluate the energy due to one atom.
!

      repeng = 0.0_dp

      do i1 = 1,nd,1

         sum = 0.0_dp

         z1 = z(i1)
         i2a = aptr(i1)

         do while (bptr(i2a) /= eol)

            i2 = bptr(i2a)
            z2 = z(i2)

            ptr = btype(z1,z2)

            dr(1) = ad(1,i2)-ad(1,i1)
            dr(2) = ad(2,i2)-ad(2,i1)
            dr(3) = ad(3,i2)-ad(3,i1)

            r = sqrt(dr(1)**2+dr(2)**2+dr(3)**2)

            if (r > 1.0e-6_dp) then
               phi = bnddat(15,ptr)* & 
     &         vscale(r,bndscl(1,15,ptr),bndscl(2,15,ptr), & 
     &                bndscl(5,15,ptr),bndscl(6,15,ptr), & 
     &                bndscl(3,15,ptr),bndscl(4,15,ptr), & 
     &                bndscl(7,15,ptr),bndscl(8,15,ptr), & 
     &                bndscl(9,15,ptr),bndscl(10,15,ptr), & 
     &                bndscl(11,15,ptr),bndscl(12,15,ptr))
               sum = sum + phi
            endif

            i2a = i2a+1

         enddo

         repeng = repeng + & 
     &            sum*(embed(1,z1)+sum*(embed(2,z1)+ & 
     &            sum*(embed(3,z1)+sum*embed(4,z1))))

      enddo

      end

