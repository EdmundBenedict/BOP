
!****************************************************************************
!                                                                           *
!                                                                           *
!               This subroutine calculates the energy/volume curves         *
!              at specified c/a ratios.                                     *
!                        Dimitar Pashov                                     *
!                                                                           *
!****************************************************************************



      subroutine getenvol ( )
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
      include "Include/Relax.array"

!      
!    Common blocks introduced by A.Girshick.      
!          
       
      include "Include/ag.conn"


      integer  arraysize
      parameter ( arraysize = 10201 ) 

      character*80 filename
      integer  i,  j,  ii, iilim
      integer  nbr( arraysize )
      integer  flag
      real(dp) ::  a0( 3, 3 )
      real(dp) ::  vstep, rtstep
      real(dp) ::  xxx, yyy, emin
      real(dp) ::  rat0, athrd
      real(dp) ::  eee(arraysize), aaa(arraysize), ccc(arraysize), vvv(arraysize), rat(arraysize)
      real(dp) ::  eeeii, aaaii, cccii



      write ( 6, '(//F12.8, 4X, F12.8, 4X, I4 //)' )  &
     &                                zminv, zmaxv, numpointsv
      write ( 6, '(//F12.8, 4X, F12.8, 4X, I4 //)' )  &
     &                                zminrt, zmaxrt, numpointsrt

      filename = genfile(1:lengfn)//'.ev'
      open(unit = 80,file = filename,status = 'NEW')

!  Determine the steps vSTEP and rtSTEP.

      if ( numpointsv  <=  1 ) then
          numpointsv = 1
          vstep = 0.d0
      else
          vstep = ( zmaxv - zminv ) / real( numpointsv - 1 , dp)
      endif

      if ( numpointsrt  <=  1 ) then
          numpointsrt = 1
          rtstep = 0.d0
      else
          rtstep = ( zmaxrt - zminrt ) / real( numpointsrt - 1 , dp)      
      endif

!  IILIM - total number of points on this energy surface.

      iilim = numpointsv * numpointsrt
      if ( iilim  >  arraysize )   then
          write ( 6, '(//10X,"INCREASE  PARAMETER  ARRAYSIZE"// &
     &               10X,"IN  THE  SUBROUTINE  GETENVOL"//)' )
          stop
      endif

!  Create the arrays of the lattice parameters of distorted lattices.

      rat0 = a(3,3)/a(2,2)
      athrd = 1.0_dp/3.0_dp

      do i = 1, numpointsv, 1
          do j = 1, numpointsrt, 1
              
              ii = ( i - 1 ) * numpointsrt + j
              vvv( ii ) = zminv + real(i-1, dp)*vstep 
              rat( ii ) = zminrt + real(j-1, dp)*rtstep
              
              aaa(ii) = exp(athrd*log(vvv(ii)*rat0/rat(ii)))
              ccc(ii) = exp(athrd*log(vvv(ii)*rat(ii)*rat(ii)/(rat0*rat0)))
              
              
              eee( ii ) = 0.d0
          enddo
      enddo

!  Memorize the original values of the A matrix.

      a0(1,1) = a(1,1)
      a0(1,2) = a(1,2)
      a0(1,3) = a(1,3)

      a0(2,1) = a(2,1)
      a0(2,2) = a(2,2)
      a0(2,3) = a(2,3)
      
      a0(3,1) = a(3,1)
      a0(3,2) = a(3,2)
      a0(3,3) = a(3,3)
 
      flag = 1
      emin = 0.d0
      xxx = 0.d0
      yyy = 0.d0

!  Cycle over all the points.

      do ii = 1, iilim, 1

         eflo = 0.d0
         do i = 1, nd, 1
             de( i ) = 0.d0
         enddo

         do i = 1,3,1
            a( i, 1 ) = a0( i, 1 ) * aaa( ii )
            a( i, 2 ) = a0( i, 2 ) * aaa( ii )
            a( i, 3 ) = a0( i, 3 ) * ccc( ii )
         enddo

         do i = 1,nd,1
            ad(1,i) = a(1,1)*d(1,i)+a(2,1)*d(2,i)+a(3,1)*d(3,i)
            ad(2,i) = a(1,2)*d(1,i)+a(2,2)*d(2,i)+a(3,2)*d(3,i)
            ad(3,i) = a(1,3)*d(1,i)+a(2,3)*d(2,i)+a(3,3)*d(3,i)
             
!       write(*,*)ad(1,i), ad(2,i), ad(3,i) 

         enddo

         do i = 1,ninert,1
            adinert(1,i) = a(1,1)*dinert(1,i)+a(2,1)*dinert(2,i)+ & 
     &                     a(3,1)*dinert(3,i)
            adinert(2,i) = a(1,2)*dinert(1,i)+a(2,2)*dinert(2,i)+ & 
     &                     a(3,2)*dinert(3,i)
            adinert(3,i) = a(1,3)*dinert(1,i)+a(2,3)*dinert(2,i)+ & 
     &                     a(3,3)*dinert(3,i)
         enddo
!
!     Adding this part to perform relaxation of internal
!     degrees of freedom if required when determining lattice
!     parameters
!
!     You need the ELPLUSREL flag in fort.8
!
!     M. Cawkwell 4th November 2004
!
         if( elplusrel  ==  1 ) then
            call relax()
         endif

         flag = 1
         call getetot ( flag )
         
!          EEE(II) = ( EPROM + EBOND + EPAIR - EATOM - UENT )
        eee(ii) = ( eprom + ebond + epair - eatom - uent )/real(nd, dp)
         nbr(ii) = maxnumnb

         if ( eee( ii )  <  emin )   then
             emin = eee( ii )
             xxx  = aaa( ii )
             yyy  = ccc( ii )
         endif



         write ( 6, '(6X, F20.12, 3X, F20.12, 3X, F10.8, 4X, F10.8, 4X, F20.12, 3X, I3)' )  &
     &                 vvv(ii), rat(ii),aaa( ii ), ccc( ii ),  eee( ii ), nbr( ii )


      enddo

!  Print out the results to the channel 80.

      write ( 80, "( /75('*')//)" )
      write ( 80, '( 2(/) )' )
      write ( 80, '( 3 ( 4X, F11.8 ), 4X, I4 )' )   &
     &                           zminv, zmaxv, vstep, numpointsv
      write ( 80, '( 3 ( 4X, F11.8 ), 4X, I4 )' )   &
     &                           zminrt, zmaxrt, rtstep, numpointsrt
      write ( 80, '( 9(/) )' )
      write ( 80, "( /75('*')//)" )
      write ( 80, '(8X,"The lowest energy found:   E = ",F16.12," eV" &
     &  //8X,"at   A = ",F12.8," A0  and   C = ",F12.8," C0" /)' ) &
     &                              emin, xxx, yyy
      write ( 80, "(/75('*')   &
     &   //7X,'Ti cohesive energy as a function of A and C' &
     &   ///12X, 'VOL', 20X, 'RATIO',12x, 'A', 12x, 'C', 17X, 'ENERGY', 15X, 'MAXNUMNB'//)" )              
      do  i = 1, iilim, 1
         write ( 80, '(6X, F20.12, 4X, F20.12, 4X, F10.8, 4X, F10.8, 4X, F20.12, 3X, I3 )' )   &
     &                  vvv(i), rat(i), aaa( i ), ccc( i ),  eee( i ),  nbr( i )


       enddo




      write ( 80, '( /75("*")//)' )

      call tag( 80 )

      close ( 80 )

!  Assign the original values to the A matrix.

      a(1,1) = a0(1,1)
      a(1,2) = a0(1,2)
      a(1,3) = a0(1,3)

      a(2,1) = a0(2,1)
      a(2,2) = a0(2,2)
      a(2,3) = a0(2,3)

      a(3,1) = a0(3,1)
      a(3,2) = a0(3,2)
      a(3,3) = a0(3,3) 



      return
      end

