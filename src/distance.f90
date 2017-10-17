 
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
!          This subroutine calculates the shortest distance between         *
!               the two atoms for a given atomic configuration.             *
!                                                                           *
!                                                                           *
!****************************************************************************

!****************************************************************************
!                                                                           *
!           Rewritten by Matous Mrovec                                      *
!                                                                           *
!****************************************************************************


      subroutine distance
          use mod_precision


          use mod_all_scalar

          use mod_const
          use mod_atom_ar, only : btype
          
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
      include "Include/PosVel.array"
      include "Include/NebListL.array"

!
!    Common blocks introduced by A.Girshick.
!
 
      include "Include/ag.conn" ! for et*totstr only


      integer  i, j, k, l, m, n, ja, ic, bt
      integer  incr(3), nnb(100)
      real(dp) ::  dist, border
      real(dp) ::  rx, ry, rz
      character*5 elm,elm2
      
      integer :: fake(1) ! 'dummy' to pass instead of mapi, it will not be touched
      real(dp) :: ett(3,3)
      
      logical :: stressed
      
      
      write(6,'(/"ENTER CENTRAL ATOM : ",$)')
      read(5,*) ic
      write(6,'(/"ENTER RANGE IN [A] : ",$)')
      read(5,*) border
!      BORDER = 3.8
!      BORDER = LATPAR
!      WRITE(6,'(/"ENTER NAME OF SPECIES : ",$)')
!      READ(5,'(A2)') ELM
      elm='AA'

      stressed = totstr /= 0.0_dp .and. any(et /= 0.0_dp)      
      
      if (stressed) then
         ett = et*totstr
         ad(:,:nd) = unshad(:,:nd)
         adinert(:,:ninert) = unshind(:,:ninert)
      end if
      
      call bldnebt(.true., rlim, lena, a, rprune, rcut, stressed, ett, &
                  & nd, ad, z, ninert, adinert, zinert, &
                  & map, aptr, bptr, maxnumnb, totnd, fake,mxnnb)

      open(unit=28,file="nnb.lst",status='NEW')
      open(unit=29,file="ablock.xbs",status='NEW')
      open(unit=30,file="wblock.xbs",status='NEW')
      write(29,'("* Output file for xbs program")')
      write(30,'("* Output file for xbs program")')

      do i = 1, nd

         write(29,'("atom  ",A5,4X,3(F10.4))') elm, &
     &               ad(1,i),ad(2,i),ad(3,i)
         write(30,'("atom  ",A5,4X,3(F14.8))') elm, &
     &               ad(1,i),ad(2,i),ad(3,i)

      enddo

         write(28,'("* Neighbor list of atom : ",I3)' ) ic
!         WRITE(28,'("* NEIGHBOUR",7X,"X",9X,"Y",9X,"Z",7X,"DIST")')
!         WRITE(28,'("atom",A4,4X,4(F10.4))') ELM, 0.0, 0.0, 0.0, 0.0

         k=0
         ja = aptrl(ic)

         do while ( bptrl( ja )  /=  eol )
            j = bptrl( ja )
            if ( j  <=  9  ) then
               write(elm,'("A",I1)') j
            else if ( j  <=  99 ) then
               write(elm,'("A",I2)') j
            else if ( j  <=  999 ) then 
               write(elm,'("A",I3)') j
            else
               write(elm,'("A",I4)') j
            endif

            rx = ad( 1, ic ) - ad( 1, j )
            ry = ad( 2, ic ) - ad( 2, j )
            rz = ad( 3, ic ) - ad( 3, j )               
            dist = sqrt( rx*rx + ry*ry + rz*rz )

            bt = btype(z(ic),z(j))

            if (dist < border) then
               write(28,'("atom  ",A5,4X,4(F10.4),4X,I2)') elm, &
!     &              AD( 1, J ),AD( 2, J ),AD( 3, J ),DIST,BT &
     &              rx, ry, rz, dist,bt
!            WRITE(28,'("spec  ",A5,4X,"0.500",4X,"Black")') ELM
!            WRITE(28,'(/)')
               k=k+1
            endif
               ja = ja + 1
         enddo

         write(28,'("* No. OF NEIGHBOURS : ",I3)' ) k-1
!         WRITE(28,'("* Atom data")' )
!         WRITE(28,'("spec  A*",4X,"0.500",4X,"Black")')
         write(28,'("* Bond data")' ) 
         write(28,'("bonds  A*  A*",4X,"0.000", &
     &               F10.4,4X,"0.1",4X,"1.00")') border
         write(28,'(/)')


      elm='AA'
      elm2='IA'
      write(30,'("Inert atoms")' )
      do i=1,ninert
         write(30,'("atom  ",A5,4X,3(F14.8))') elm2, &
     &             adinert(1,i),adinert(2,i),adinert(3,i)
      enddo

         write(29,'("* Atom data")' )
         write(29,'("spec  ",A5,4X,"0.50",4X,"Black")') elm
         write(29,'("* Bond data")' )
         write(29,'("bonds  ",2(A5),4X,"0.000", &
     &               F10.4,4X,"0.1",4X,"1.00")') elm,elm,border

         write(30,'("* Atom data")' )
         write(30,'("spec  ",A5,4X,"0.40",4X,"Black")') elm
         write(30,'("spec  ",A5,4X,"0.40",4X,"White")') elm2
         write(30,'("* Bond data")' )
         write(30,'("bonds  ",2(A5),4X,"0.000", &
     &               F10.4,4X,"0.1",4X,"1.00")') elm,elm,border

      close(28)
      close(29)
      close(30)

!  Set up the neighbor lists for the Finnis-Sinclair potential.

!          CALL BLDNEBTFS()

!  Calculate the Finnis-Sinclair energy and the corresponding forces.
      
!          CALL FINNIS()


!      WRITE( 6, '(/F10.7/)' ) DIST0
!      WRITE( 6, '(/I4, 4X, I4)' ) I0, J0
!      WRITE( 6, '(3(/F10.7, 4X, F10.7))' ) 
!     &    AD( 1, I0 ) / LATPAR + XXMIN, AD( 1, J0 ) / LATPAR + XXMIN,
!     &    AD( 2, I0 ) / LATPAR + YYMIN, AD( 2, J0 ) / LATPAR + YYMIN,
!     &    AD( 3, I0 ) / LATPAR + ZZMIN, AD( 3, J0 ) / LATPAR + ZZMIN


!      STOP

      end

