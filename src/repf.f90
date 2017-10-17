 
      subroutine repf()
          use mod_precision


          use mod_all_scalar

          use mod_const
          use mod_atom_ar
!
!    This is a subroutine to evaluate the forces due to the repulsive term.
!

      implicit none

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
      include "Include/Force.array"
      include "Include/NebList.array"
      include "Include/PosVel.array"

!
!    Declare the simple variables.
!

      real(dp) :: dphi,phi,sum,dfi,ff
      real(dp) :: r,vscale,dscale

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
!    Initialise the forces.
!

      do i1 = 1,nd,1
         fpp(1,i1) = 0.0_dp
         fpp(2,i1) = 0.0_dp
         fpp(3,i1) = 0.0_dp
      enddo

      rdfpp = 0.0_dp

!
!    Evaluate the forces.
!

      do i1 = 1,nd,1

         z1 = z(i1)
         i2a = aptr(i1)

         sum = 0.0_dp

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

         dfi = embed(1,z1)+sum*(2.0_dp*embed(2,z1)+ & 
     &         sum*(3.0_dp*embed(3,z1)+sum*4.0_dp*embed(4,z1)))

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

               dphi = bnddat(15,ptr)* & 
     &         vscale(r,bndscl(1,15,ptr),bndscl(2,15,ptr), & 
     &                bndscl(5,15,ptr),bndscl(6,15,ptr), & 
     &                bndscl(3,15,ptr),bndscl(4,15,ptr), & 
     &                bndscl(7,15,ptr),bndscl(8,15,ptr), & 
     &                bndscl(9,15,ptr),bndscl(10,15,ptr), & 
     &                bndscl(11,15,ptr),bndscl(12,15,ptr))

               ff = dfi*dphi/r

               rdfpp = rdfpp + (dr(1)**2+dr(2)**2+dr(3)**2)*ff

               if (i2 <= nd) then

                  fpp(1,i1) = fpp(1,i1) + ff*dr(1)
                  fpp(2,i1) = fpp(2,i1) + ff*dr(2)
                  fpp(3,i1) = fpp(3,i1) + ff*dr(3)

                  fpp(1,i2) = fpp(1,i2) - ff*dr(1)
                  fpp(2,i2) = fpp(2,i2) - ff*dr(2)
                  fpp(3,i2) = fpp(3,i2) - ff*dr(3)

               elseif (map(i2) > 0) then

                  fpp(1,i1) = fpp(1,i1) + ff*dr(1)
                  fpp(2,i1) = fpp(2,i1) + ff*dr(2)
                  fpp(3,i1) = fpp(3,i1) + ff*dr(3)

                  fpp(1,map(i2)) = fpp(1,map(i2)) - ff*dr(1)
                  fpp(2,map(i2)) = fpp(2,map(i2)) - ff*dr(2)
                  fpp(3,map(i2)) = fpp(3,map(i2)) - ff*dr(3)

               else

! Inert atom neighbor. Curses!!!
                  fpp(1,i1) = fpp(1,i1) + 2.0_dp*ff*dr(1)
                  fpp(2,i1) = fpp(2,i1) + 2.0_dp*ff*dr(2)
                  fpp(3,i1) = fpp(3,i1) + 2.0_dp*ff*dr(3)

               endif

            endif

            i2a = i2a+1

         enddo

      enddo

      end

