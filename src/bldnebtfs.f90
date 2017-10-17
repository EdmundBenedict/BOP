 
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
!         This subroutine builds a "long cutoff" neighboring table          *
!                 for the Finnis-Sinclair potentials scheme.                *
!                                                                           *
!                                                                           *
!****************************************************************************


      subroutine bldnebtfs()
          use mod_precision
          use mod_all_scalar
          use mod_const
          use mod_clock

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

!      include "Include/Atom.array"
      include "Include/NebList.array"
      include "Include/NebListL.array"
      include "Include/PosVel.array"
!      include "Include/RecurCoef.array"

!
!    Common blocks introduced by A.Girshick.
!    
 
      include "Include/ag.conn"

!
!    Declare the simple variables.
!

      real(dp) :: rcutl2
      real(dp) :: ia_pos(3)

      integer inda,indb,totnbox
      integer ia_n(3),ib1,ib2,ib3
      integer ia,ib,ic
      integer ja,jb
      integer nb
      integer ns

!
!    Declare local arrays.
!

!*** Size of periodically repeated unit cell.
      real(dp) :: span(3)
!*** Origin for total cell of atoms.
      real(dp) :: org(3)
!*** Size of one of partitions.
      real(dp) :: dl(3)
!*** DLIM sets the limits on the number of neighboring cells to be scanned.
      real(dp) :: dlim(2,3)
!*** Extremities of atomic coordinates.
      real(dp) :: rmax(3)
      real(dp) :: rmin(3)

!*** NNA is an array of the number of neighbors for each atom
!     INTEGER NNA(MXTOTND)
!*** The number of boxes in a given direction used to partition space.
      integer nbox(3)

!*** The number of atoms in a box.
      integer boxcount(mbox)

!*** The list of atoms for each box.
      integer box(mcount,mbox)

      
      
      
      integer(8) :: c1,c2

!
!    Impose periodic boundary conditions on the atoms.
!
      call system_clock(c1)
      if (mvflag /= 10) call pbc(ad,nd,adinert,ninert,a,rlim)
      call system_clock(c2)
      print *,'pbc',real(c2-c1,dp)/real(cr,dp)
!
!    Build list of atoms.
!

      call system_clock(c1)
      call bldlistfs(dlim)
      call system_clock(c2)
      print *,'bldlistfs',real(c2-c1,dp)/real(cr,dp)

      ns = totnd

      rcutl2 = rcutl**2

!
!    Build the neighbor list - this is simple O(N^2) method.
!
!       print *,'nbtflg',nbtflg
      call system_clock(c1)
      if (nbtflg == 1) then
      
         nb = 1
         do ia = 1,nne,1
            aptrl(ia) = nb
            do ib = 1,ns,1
                
               if (sum((ad(:,ia) - ad(:,ib))**2) <= rcutl2 .and. ia /= ib) then
                  bptrl(nb) = ib
                  nb = nb + 1
               endif
            enddo
            if (nb-aptrl(ia) > mxnnbl) then
                write(6,'(''Number of neighbors for atom # '', I5, '' = '',I3)') ia,nb-aptrl(ia)
                write(6,'(''This is more than the maximum allowed of '', I3)') mxnnbl
                write(6,'(''=> Increase MXNNBL.'')')
               call panic()
            endif

!            WRITE(6,'(2I5)') IA, NB-APTRL(IA)

            bptrl(nb) = eol
            nb = nb + 1

         enddo

!
!    Build the neighbor list - this is the linked list O(N) method.
!     Note: this requires a tetrahedral unit cell.
!

      elseif (nbtflg == 2) then
         if ((rlim(1) == 0).or.(rlim(2) == 0).or.(rlim(3) == 0)) then
            rmin(1) = ad(1,1)
            rmin(2) = ad(2,1)
            rmin(3) = ad(3,1)
            rmax(1) = ad(1,1)
            rmax(2) = ad(2,1)
            rmax(3) = ad(3,1)
            do ia = 2,nd,1
               rmin(1) = min(rmin(1),ad(1,ia))
               rmin(2) = min(rmin(2),ad(2,ia))
               rmin(3) = min(rmin(3),ad(3,ia))
               rmax(1) = max(rmax(1),ad(1,ia))
               rmax(2) = max(rmax(2),ad(2,ia))
               rmax(3) = max(rmax(3),ad(3,ia))
            enddo
            do ia = 1,ninert,1
               rmin(1) = min(rmin(1),adinert(1,ia))
               rmin(2) = min(rmin(2),adinert(2,ia))
               rmin(3) = min(rmin(3),adinert(3,ia))
               rmax(1) = max(rmax(1),adinert(1,ia))
               rmax(2) = max(rmax(2),adinert(2,ia))
               rmax(3) = max(rmax(3),adinert(3,ia))
            enddo
         endif

         if (rlim(1) == 0) then
            span(1) = rmax(1)-rmin(1)+0.001_dp
            org(1) = rmin(1)
         else
            span(1) = lena(1) + 2.0_dp*rcutl
            org(1) = dlim(1,1)*a(1,1)+dlim(1,2)*a(2,1)+dlim(1,3)*a(3,1)
         endif
         if (rlim(2) == 0) then
            span(2) = rmax(2)-rmin(2)+0.001_dp
            org(2) = rmin(2)
         else
            span(2) = lena(2) + 2.0_dp*rcutl
            org(2) = dlim(1,1)*a(1,2)+dlim(1,2)*a(2,2)+dlim(1,3)*a(3,2)
         endif
         if (rlim(3) == 0) then
            span(3) = rmax(3)-rmin(3)+0.001_dp
            org(3) = rmin(3)
         else
            span(3) = lena(3) + 2.0_dp*rcutl
            org(3) = dlim(1,1)*a(1,3)+dlim(1,2)*a(2,3)+dlim(1,3)*a(3,3)
         endif

         nbox(1) = max(1,int(span(1)/rcutl))
         nbox(2) = max(1,int(span(2)/rcutl))
         nbox(3) = max(1,int(span(3)/rcutl))

         if (idebug == 1) then
            write(6,'(''RCUT(L) = '',G13.5)') rcutl
            write(6,'(''SPAN(L) = '',3G13.5)') span
            write(6,'(''ORG(L)  = '',3G13.5)') org
            write(6,*) 'NBOX(L) = ',nbox
         end if

         totnbox = nbox(1)*nbox(2)*nbox(3)
         if (totnbox > mbox) then
               write(6,'(''Too many boxes for neighbor list calculation.'')')
               write(6,'(''Increase MBOX to at least '',I5)') totnbox
            call panic()
         endif

         dl(1) = span(1)/real(nbox(1), dp)
         dl(2) = span(2)/real(nbox(2), dp)
         dl(3) = span(3)/real(nbox(3), dp)
         if (idebug == 1) write(6,'(''DL = '',3G13.5)') dl

         do inda = 1,totnbox,1
            boxcount(inda) = 0
         enddo

         do ia = 1,ns,1

            ia_n = int((ad(:,ia)-org)/dl(:))
            inda = 1 + ia_n(1) + nbox(1)*(ia_n(2) + nbox(2)*ia_n(3))
            if (idebug == 1) then
               if ((inda < 1).or.(inda > totnbox)) then
                     write(6,'(''Attempting to store atom in illegal box.'')')
                     write(6,'(''INDA = '',I5)') inda
                  call panic()
               endif
            endif

            boxcount(inda) = boxcount(inda) + 1
            if (boxcount(inda) > mcount) then
                  write(6,'(''Too many atoms in one box.'')')
                  write(6,'(''Increase MCOUNT.'')')
               call panic()
            endif

            box(boxcount(inda),inda) = ia

         enddo

         nb = 1
         aptrl(1) = nb
         do ia = 1,nne,1
          
            ia_pos = ad(:,ia)

            ia_n = int((ia_pos-org)/dl)

            do ib1 = max(0,ia_n(1)-1),min(nbox(1)-1,ia_n(1)+1),1
               do ib2 = max(0,ia_n(2)-1),min(nbox(2)-1,ia_n(2)+1),1
                  do ib3 = max(0,ia_n(3)-1),min(nbox(3)-1,ia_n(3)+1),1
                     indb = 1 + ib1 + nbox(1)*(ib2 + nbox(2)*ib3)
                     do jb = 1,boxcount(indb),1
                        ib = box(jb,indb)

                        if (sum((ia_pos - ad(:,ib))*(ia_pos - ad(:,ib))) <= rcutl2 .and. ia /= ib) then
                              bptrl(nb) = ib
                              nb = nb + 1
                        endif
                     enddo
                  enddo
               enddo
            enddo

            if (nb-aptrl(ia) > mxnnbl) then
                 write(6,'(''Number of neighbors for atom # '',I5, '' = '', I3)') ia, nb - aptrl(ia) - 2
                 write(6,'(''This is more than the maximum allowed of '', I3)') mxnnbl
                 write(6,'(''=> Increase MXNNBL.'')')
               call panic()
            endif

!                 WRITE(6,'(''Number of neighbors for atom # '',I5,
!     +                  '' = '',I3)') IA, NB - APTRL(IA) - 2

            bptrl(nb) = eol
            nb = nb + 1
            aptrl(ia+1) = nb

         enddo

      endif
      call system_clock(c2)
      print *,'rest',real(c2-c1,dp)/real(cr,dp)


      return
      end

