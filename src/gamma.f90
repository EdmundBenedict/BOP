 
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
!       This subroutine manages the relaxation of the gamma-surface.        *
!                                                                           *
!                                                                           *
!****************************************************************************


      subroutine  gammasurf ( )
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
      include "Include/PosVel.array"

!
!    Common blocks introduced by A.Girshick.
!    
 
      include "Include/ag.conn"


      real(dp) ::  ezero, surf
      integer           i, j, ii
      integer           flag

      character*80 filename


!
!    Declare the arrays.
!

      integer  arraysize
      parameter ( arraysize = 10201 )
      integer  nbr ( arraysize )
      real(dp) ::  xshift ( arraysize ), yshift ( arraysize )
      real(dp) ::  eee ( arraysize )
      

      filename = genfile(1:lengfn)//'.gs'
      open(unit = 80,file = filename,status = 'NEW')
      

!  Calculate the gamma-surface mesh.

      if ( gamxnp  <=  1 ) then  
          gamxnp = 1
          gamxstep = 0.d0
      else
          gamxstep = ( gamxmax - gamxmin ) / real( gamxnp - 1 , dp)
      endif
              
      if ( gamynp  <=  1 ) then
          gamynp = 1
          gamystep = 0.d0
      else
          gamystep = ( gamymax - gamymin ) / real( gamynp - 1 , dp)
      endif
      
      pointsgs = gamxnp * gamynp


      write ( 6, '(//F10.6, 4X, F10.6, 4X, I4, 4X, F10.6 //)' )  &
     &                gamxmin, gamxmax, gamxnp, gamxstep
      write ( 6, '(//F10.6, 4X, F10.6, 4X, I4, 4X, F10.6 //)' )  &
     &                gamymin, gamymax, gamynp, gamystep
      write ( 6,"(//'INDRELAX is ',I1,'   ZEROHEIGHT is ',F16.12//)")  &
     &                indrelax, zeroheight

      gamxmin  = gamxmin  * a( 1, 1 )
      gamxmax  = gamxmax  * a( 1, 1 ) + 1.0e-6_dp
      gamxstep = gamxstep * a( 1, 1 )

      gamymin  = gamymin  * a( 2, 2 )
      gamymax  = gamymax  * a( 2, 2 ) + 1.0e-6_dp
      gamystep = gamystep * a( 2, 2 )

      write ( 80, '( 3(/) )' )
      write ( 80, '("Gamma-surface Calculation (BOP scheme)")')
!      WRITE ( 80, '( " Z - [   0  1 -1  0  ]" )' )
!      WRITE ( 80, '( " X - [   1  0  0  0  ]" )' )
      write ( 80, '( 3( 4X, F7.4 )  )' )  gamxmin, gamxmax, gamxstep
      write ( 80, '( 3( 4X, F7.4 )  )' )  gamymin, gamymax, gamystep
      if (req  ==  11 )     write(80,'("UN")')
      if (req  ==  21 )     write(80,'("MX")')
      if (req  ==  12 )     write(80,'("MY")')
      if (req  ==  22 )     write(80,'("MI")')


!  Create the arrays of displacements in the plane of the boundary.
      
      do i = 1, gamxnp, 1
          do j = 1, gamynp, 1
              ii = ( i - 1 ) * gamynp + j
              xshift( ii ) = ( gamxmin  +  real(i-1, dp) * gamxstep ) 
              yshift( ii ) = ( gamymin  +  real(j-1, dp) * gamystep ) 
              eee( ii ) = 0.d0
          enddo
      enddo


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

!Dimitar Pashov edit
!       write(6,*) zeroheight
!end edit
!     SURF is a conversion factor to express 
!     the gamma-surface energy in mJ/m**2

      surf = a( 1, 1 ) * a( 2, 2 ) / 16020.3

      eflo = 0.d0
      do i = 1, nd, 1
          de( i ) = 0.d0
      enddo

      flag = 1
      call getetot ( flag )

      ezero = ( eprom + ebond + epair - eatom - uent ) / surf

      write ( 6, '(F16.6)' ) eprom + ebond + epair - eatom - uent
      write ( 6, '(F16.6)' ) ezero

!  Cycle on all the points.

      do ii = 1, pointsgs, 1

          write ( 72, "( //80('*')/)" )        
          write ( 72, '(" X = ",F6.4,"   Y = ",F6.4/ )' )  &
     &            xshift( ii ) , yshift( ii )
     
          eflo = 0.d0
          do i = 1, nd, 1
              de( i ) = 0.d0
          enddo

!         Shift the coordinates of the upper part of the block 
!         for active atoms.

          do i = 1,nd,1
              if ( unshad( 3, i )  >=  zeroheight )  then
                  write(6,'("*** SHIFTING ACTIVE ATOM ***")')
                  ad( 1, i ) = unshad( 1, i ) + xshift( ii ) 
                  ad( 2, i ) = unshad( 2, i ) + yshift( ii )
                  ad( 3, i ) = unshad( 3, i )
              else
                  ad( 1, i ) = unshad( 1, i )
                  ad( 2, i ) = unshad( 2, i )   
                  ad( 3, i ) = unshad( 3, i )
              endif
          enddo

!         Save the original unshifted coordinates.
  
          if ( pointsgs  ==  1 ) then
              do i = 1,nd,1
                  do j = 1, 3, 1
                      shcrd( j, i ) = ad( j, i )
                  enddo
              enddo           
          endif

!         Shift the coordinates of the upper part of the block 
!         for inert atoms.

          do i = 1,ninert,1
              if ( unshind( 3, i )  >=  zeroheight )  then
                  adinert( 1, i ) = unshind( 1, i ) + xshift( ii )
                  adinert( 2, i ) = unshind( 2, i ) + yshift( ii )           
                  adinert( 3, i ) = unshind( 3, i )
              else
                  adinert( 1, i ) = unshind( 1, i )
                  adinert( 2, i ) = unshind( 2, i )
                  adinert( 3, i ) = unshind( 3, i )
              endif
          enddo

     
          if ( indrelax  ==  0 )   then

              flag = 1
              call getetot ( flag )

          elseif ( indrelax  ==  1 )   then
      
              iter = 0
              call relax_sb ( ezero, surf )

          endif

          eee(ii) = ( (eprom+ebond+epair-eatom-uent)/surf ) - ezero
          nbr(ii) = maxnumnb
          write(*,'(F16.8)') eprom+ebond+epair-eatom-uent
          write ( 80, '( 6X, F14.8, 4X, F14.8, 4X, F16.8, 3X, I3 )' ) &
     &                    xshift(ii), yshift(ii),  eee(ii),  nbr(ii)
          call flush( 80 )

      enddo

      write ( 80, "( /80('*')//)" )

      call tag( 80 )

      close ( 80 )


      return
      end

