 
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
!              This subroutine manages the relaxation of the                *
!                    grain boundary or a stacking fault.                    *
!                                                                           *
!                                                                           *
!****************************************************************************


      subroutine  grainbound ( )
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
 
!      include "Include/Atom.array"
      include "Include/PosVel.array"

!
!    Common blocks introduced by A.Girshick.
!    
 
      include "Include/ag.conn"


!    Declare the variables.

      integer           i, j, flag, iz
      real(dp) ::  ezero, eshift, surf, zeroheightup
      character*80 filename

      filename = genfile(1:lengfn)//'.gb'
      open(unit = 80,file = filename,status = 'NEW')
      open(unit = 81,file = "shifted.block",status = 'NEW')

!  Memorize the original unshifted coordinates.

      do i = 1, nd, 1
          do j = 1, 3, 1
              unshad( j, i ) = ad( j, i )
          enddo
      enddo

      do i = 1, ninert, 1
          do j = 1, 3, 1
              unshind( j, i ) = adinert( j, i )
          enddo
      enddo
    
      zeroheight = zeroheight * a( 3, 3 ) - 1.d-6
      zeroheightup = zeroheight + 1.d-6
      write(*,'(2(F16.12))') zeroheight, zeroheightup


!     SURF is a conversion factor to express 
!     the gamma-surface energy in mJ/m**2

      surf = a( 1, 1 ) * a( 2, 2 ) / 16020.3


!  Calculate the energy of the ideal crystal EZERO.
 
      flag = 1
      call getetot ( flag )


      ezero = ( eprom + ebond + epair - eatom - uent ) / surf
      write(80,'("Energy of nonshifted and nonrelaxed GB")')
      write(80,'(" E[eV] = ",F16.8/)') ezero*surf
      write(80,'(" Ezero = ",F16.8/)') ezero
      write(80,'(" XSHIFTGB = ",F16.8/)') xshiftgb
      write(80,'(" YSHIFTGB = ",F16.8/)') yshiftgb
      write(80,'(" ZSHIFTGB = ",F16.8/)') zshiftgb

!     Shift the coordinates of the upper part of the block 
!     for active atoms.

      do i = 1,nd,1
          if ( unshad( 3, i )  >=  zeroheight )  then
             if ( unshad( 3, i )  <=  zeroheightup ) then
              ad( 1, i ) = unshad( 1, i ) + xshiftgb/2 * a(1,1)
              if ( ad( 1, i )  >  a(1,1) ) then
                   ad( 1, i ) = ad( 1, i ) - a(1,1)
                   write(80,'(" ATOM = ",I4/)') i
              endif
              ad( 2, i ) = unshad( 2, i ) + yshiftgb/2 * a(2,2)
              if ( ad( 2, i )  >  a(2,2) ) then
                   ad( 2, i ) = ad( 2, i ) - a(2,2)
                   write(80,'(" ATOM = ",I4/)') i
              endif
              ad( 3, i ) = unshad( 3, i ) + zshiftgb/2 * a(3,3)
             else 
              ad( 1, i ) = unshad( 1, i ) + xshiftgb * a(1,1)
              if ( ad( 1, i )  >  a(1,1) ) then
                   ad( 1, i ) = ad( 1, i ) - a(1,1)
                   write(80,'(" ATOM = ",I4/)') i
              endif
              ad( 2, i ) = unshad( 2, i ) + yshiftgb * a(2,2)
              if ( ad( 2, i )  >  a(2,2) ) then
                   ad( 2, i ) = ad( 2, i ) - a(2,2)
                   write(80,'(" ATOM = ",I4/)') i
              endif
              ad( 3, i ) = unshad( 3, i ) + zshiftgb * a(3,3)
             endif
          else
              ad( 1, i ) = unshad( 1, i )
              ad( 2, i ) = unshad( 2, i )
              ad( 3, i ) = unshad( 3, i )
          endif
      write(81,'(3(2X,F16.12))') ad(1,i), ad(2,i), ad(3,i)
      enddo


!     Save the original shiffted coordinates.

      do i = 1,nd,1
          do j = 1, 3, 1
              shcrd( j, i ) = ad( j, i )
          enddo
      enddo           


!     Shift the coordinates of the upper part of the block 
!     for inert atoms.

      do i = 1,ninert,1
          if ( unshind( 3, i )  >=  zeroheight )  then
              adinert( 1, i ) = unshind( 1, i ) + xshiftgb*a(1,1)
              if ( adinert( 1, i )  >  a(1,1) ) then
                   adinert( 1, i ) = adinert( 1, i ) - a(1,1)
                   write(80,'(" IATOM = ",I4/)') i
              endif
              adinert( 2, i ) = unshind( 2, i ) + yshiftgb*a(2,2)
              if ( adinert( 2, i )  >  a(2,2) ) then
                   adinert( 2, i ) = adinert( 2, i ) - a(2,2)
                   write(80,'(" IATOM = ",I4/)') i
              endif
              adinert( 3, i ) = unshind( 3, i ) + zshiftgb * a(3,3)
          else
              adinert( 1, i ) = unshind( 1, i )
              adinert( 2, i ) = unshind( 2, i )
              adinert( 3, i ) = unshind( 3, i )
          endif
      write(81,'(3(2X,F16.12))') adinert(1,i),adinert(2,i),adinert(3,i)
      enddo


!     Call the relaxation subroutine.

      iter = 0
      call relax_sb ( ezero, surf )

! changed temporarily by M.M.

      a(3,3)=a(3,3)+zshiftgb*a(3,3)
      write(*,'(" A(3,3) = ",F16.8/)') a(3,3)
      flag = 1
      call getetot (flag)
      eshift = ( eprom + ebond + epair - eatom - uent ) / surf
      write(80,'(" E[eV]  = ",F16.8)') eshift * surf
      write(80,'(" ESHIFT = ",F16.8)') eshift
      write(80,'(" DELTA Esh-Ez = ",F16.8/)') eshift-ezero
      write(80,'(F16.8,F16.8)') yshiftgb, eshift-ezero

      close(80)
      close(81)

      return
      end

