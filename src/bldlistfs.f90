 
!****************************************************************************
!                                                                           *
!                                                                           *
!                             Alexey Girshick                               *
!              Department of Materials Science and Engineering              *
!                        University of Pennsylvania                         *
!                                                                           *
!                                                                           *
!****************************************************************************
!                                                                           *
!                                                                           *
!              This subroutine builds a complete list of atoms              * 
!                 for the Finnis-Sinclair potentials scheme.                *
!                                                                           *
!                                                                           *
!****************************************************************************


      subroutine bldlistfs(dlim)
          use mod_precision


          use mod_all_scalar

          use mod_const


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
      include "Include/NebList.array"
      include "Include/PosVel.array"

!
!    Common blocks introduced by A.Girshick.
!    
 
      include "Include/ag.conn"

!
!    Declare the simple variables.
!

      integer ia
      integer i1,i2,i3
      integer ns
      integer i,j

!
!    Declare the local arrays.
!

!*** A is the array of primitive translation vectors, and B is the inverse.
      real(dp) :: b(3,3)
!*** The position of one atom.
      real(dp) :: dr(3)
      real(dp) :: pos(3)
      real(dp) :: ld(3)
!*** DLIM sets the limits on the number of neighboring cells to be scanned.
      real(dp) :: dlim(2,3)
      integer  rliml(3)

      real(dp) :: zstr(3)


!
!  If this is stress application, put back the unstressed coordinates.
!

      if (mvflag == 17) then

!         WRITE ( *, '("APPS FS1 = ",F8.4)') TOTSTR

         do i = 1, nd, 1
            do j = 1, 3, 1
               ad( j, i ) = unshad( j, i )
            enddo
         enddo

         do i = 1, ninert, 1
            do j = 1, 3, 1
               adinert( j, i ) = unshind( j, i )
            enddo
         enddo

      endif



      do i = 1, 3
          if ( rlim( i )  ==  0 )   then
               rliml( i ) = 0
          else
               rliml( i ) = int( rcutl/lena(i) ) + 1
          endif
      enddo

!
!    Build list of atoms.
!

!  First include in the list active atoms from the central cell.

      do ia = 1,nd,1
         map(ia) = ia
      enddo


!  Now include in the list inert atoms from the central cell.

      do ia = 1,ninert,1
         ns = nd + ia
         ad( 1, ns ) = adinert( 1, ia )
         ad( 2, ns ) = adinert( 2, ia )
         ad( 3, ns ) = adinert( 3, ia )
         map( ns ) = ns
         z( ns ) = zinert( ia )
      enddo   

      nne = nd + ninert
      ns = nne


!  Include in the list active atoms from the identical cells.

      lena(1) = sqrt(a(1,1)**2 + a(1,2)**2 + a(1,3)**2)
      lena(2) = sqrt(a(2,1)**2 + a(2,2)**2 + a(2,3)**2)
      lena(3) = sqrt(a(3,1)**2 + a(3,2)**2 + a(3,3)**2)

      dlim(1,1) = -rcutl/lena(1)
      dlim(2,1) = 1.0_dp + rcutl/lena(1)

      dlim(1,2) = -rcutl/lena(2)
      dlim(2,2) = 1.0_dp + rcutl/lena(2)

      dlim(1,3) = -rcutl/lena(3)
      dlim(2,3) = 1.0_dp + rcutl/lena(3)

      do i = 1,3,1
         do j = 1,3,1
            b(j,i) = a(j,i)
         enddo
      enddo
      call inv3x3(b)

      do i1 = -rliml(1),rliml(1),1
         do i2 = -rliml(2),rliml(2),1
            do i3 = -rliml(3),rliml(3),1
               if ((i1 /= 0).or.(i2 /= 0).or.(i3 /= 0)) then
                  dr(1) = real(i1, dp)*a(1,1)+real(i2, dp)*a(2,1)+real(i3, dp)*a(3,1)
                  dr(2) = real(i1, dp)*a(1,2)+real(i2, dp)*a(2,2)+real(i3, dp)*a(3,2)
                  dr(3) = real(i1, dp)*a(1,3)+real(i2, dp)*a(2,3)+real(i3, dp)*a(3,3)
                  do ia = 1,nd,1
                     pos(1) = dr(1)+ad(1,ia)
                     pos(2) = dr(2)+ad(2,ia)
                     pos(3) = dr(3)+ad(3,ia)
                     call mul3x3(b,pos,ld)
                     if ( & 
     &                   ((rliml(1) == 0).or. & 
     &                    ((ld(1) >= dlim(1,1)).and. & 
     &                     (ld(1) <= dlim(2,1)))) & 
     &                  .and. & 
     &                   ((rliml(2) == 0).or. & 
     &                    ((ld(2) >= dlim(1,2)).and. & 
     &                     (ld(2) <= dlim(2,2)))) & 
     &                  .and. & 
     &                   ((rliml(3) == 0).or. & 
     &                    ((ld(3) >= dlim(1,3)).and. & 
     &                     (ld(3) <= dlim(2,3)))) & 
     &                  ) then
                        ns = ns + 1
                        ad(1,ns) = pos(1)
                        ad(2,ns) = pos(2)
                        ad(3,ns) = pos(3)
                        map(ns) = ia
                        z(ns) = z(ia)
                        if (ns > mxtotnd) then
                           write(6,'(''Too many atoms. Increase MXTOTND.'')')
                           call panic()
                        endif
                     endif
                  enddo
               endif
            enddo
         enddo
      enddo

!
!    Add inert atoms to the end of the list.
!

      do i1 = -rliml(1),rliml(1),1
         do i2 = -rliml(2),rliml(2),1
            do i3 = -rliml(3),rliml(3),1
               if ((i1 /= 0).or.(i2 /= 0).or.(i3 /= 0)) then
                  dr(1) = real(i1, dp)*a(1,1)+real(i2, dp)*a(2,1)+ real(i3, dp)*a(3,1)
                  dr(2) = real(i1, dp)*a(1,2)+real(i2, dp)*a(2,2)+ real(i3, dp)*a(3,2)
                  dr(3) = real(i1, dp)*a(1,3)+real(i2, dp)*a(2,3)+ real(i3, dp)*a(3,3)
                  do ia = 1,ninert,1
                     pos(1) = dr(1)+adinert(1,ia)
                     pos(2) = dr(2)+adinert(2,ia)
                     pos(3) = dr(3)+adinert(3,ia)
                     call mul3x3(b,pos,ld)
                     if ( & 
     &                   ((rliml(1) == 0).or. & 
     &                    ((ld(1) >= dlim(1,1)).and. & 
     &                     (ld(1) <= dlim(2,1)))) & 
     &                  .and. & 
     &                   ((rliml(2) == 0).or. & 
     &                    ((ld(2) >= dlim(1,2)).and. & 
     &                     (ld(2) <= dlim(2,2)))) & 
     &                  .and. & 
     &                   ((rliml(3) == 0).or. & 
     &                    ((ld(3) >= dlim(1,3)).and. & 
     &                     (ld(3) <= dlim(2,3)))) & 
     &                  ) then
                        ns = ns + 1
                        ad(1,ns) = pos(1)
                        ad(2,ns) = pos(2)
                        ad(3,ns) = pos(3)
                        map(ns) = nd + ia
                        z(ns) = zinert(ia)
                        if (ns > mxtotnd) then
                           write(6,'(''Too many atoms. Increase MXTOTND.'')')
                           call panic()
                        endif
                     endif
                  enddo
               endif
            enddo
         enddo
      enddo

      totnd = ns

!
! All atoms have to be strained when stresses are applied.
!

      if (mvflag == 17) then

!      WRITE ( *, '("APPLS FS2= ",F8.4)') TOTSTR

         do  i = 1, totnd

            do j=1,3
               strad(j,i) = et(j,1) * ad(1,i) &
                        & + et(j,2) * ad(2,i) &
                        & + et(j,3) * ad(3,i)
            enddo
!            PRINT *, I, STRAD(1,I)*TOTSTR, STRAD(2,I)*TOTSTR, 
!     +           STRAD(3,I)*TOTSTR
         
            do  j = 1, 3
               ad(j,i) = ad(j,i) + strad(j,i)*totstr
            enddo
         enddo
!
!     Modified by M. Cawkwell 2/2/2005 
!
!          DO  I = 1, NINERT
!
!            DO J=1,3
!               STRAD(J,I) = ET(J,1) * ADINERT(1,I) +
!     &                      ET(J,2) * ADINERT(2,I) +
!     &                      ET(J,3) * ADINERT(3,I)
!            ENDDO
!         
!            DO  J = 1, 3
!               ADINERT(J,I) = ADINERT(J,I) + STRAD(J,I)*TOTSTR
!            ENDDO
!
!         ENDDO
!
      endif

      return
      end

