
   subroutine writeddp( fln )
      use mod_precision
      use mod_const
      use mod_all_scalar
      use mod_io

      implicit none

      character(len=*), intent(in) :: fln
      
      include "Include/PosVel.array"
      include "Include/ag.conn"

      
      real(dp), pointer ::  crd(:,:), ucrd(:,:)
      integer, allocatable :: indpl(:)
      integer i, j, ii, total, num
      character*80  filename
      real(dp) :: strad(3)


!  Let's first transfer coordinates back to the original scale.
!  We have to do it only for the relaxed coordinates, because the unrelaxed 
!  coordinates were never transfered.

      if ( mvflag  /=  17 ) then         
         crd => ad( 1:3, 1:nd) !/ latpar 
      else
!             crd(1:3, :nd) = ( ad( 1:3, :nd) - strad( 1:3, :nd) * totstr )!/ latpar
         crd => unshad(1:3, 1:nd) !/ latpar
      end if
      
!       if (mvflag /= 30) then
         ucrd => unrld(1:3, 1:nd)
!       else
!          ucrd => ad( 1:3, 1:nd)
!       endif
      
!  Added by M.M.
!  Let's also transfer unrelaxed coordinates
!
!      DO I = 1, ND
!          IF ( MVFLAG .EQ. 16 ) THEN
!              UNRLD( 1, I ) = UNRLD( 1, I ) + XXMIN
!              UNRLD( 2, I ) = UNRLD( 2, I ) + YYMIN
!              UNRLD( 3, I ) = UNRLD( 3, I ) + ZZMIN
!          ENDIF
!      ENDDO


!  Now let's pick up those atoms which are within the plotting window.
!  Discrimination is done based on the unrelaxed coordinates.

      allocate(indpl(nd))
      ii = 0
      do i = 1, nd
          if ( (( ucrd( 1, i )  >=  plxmin ) .and. &
     &         ( ucrd( 1, i )  <=  plxmax ) .and. &
     &         ( ucrd( 2, i )  >=  plymin ) .and. &
     &         ( ucrd( 2, i )  <=  plymax )) &
     &     .or.((plxmin == 0.0_dp) .and. &
     &          (plxmax == 0.0_dp) .and. &
     &          (plymin == 0.0_dp) .and. &
     &          (plymax == 0.0_dp))   )    then
              ii = ii + 1
              indpl( ii ) = i
          end if
      end do
      total = ii


!     Open the plot file.

!      IF ( KEY .EQ. 0 ) THEN 
!          IF      ( ITER .LE. 99  ) THEN
!             WRITE(FILENAME,'("st.interm.",I2)') ITER
!          ELSE IF ( ITER .LE. 999 ) THEN
!              WRITE(FILENAME,'("st.interm.",I3)') ITER
!          ELSE
!              WRITE(FILENAME,'("st.interm.",I4)') ITER
!          ENDIF
!          OPEN(UNIT = 87,FILE = FILENAME)
!          NUM = 87
!      ELSE IF ( KEY .EQ. 1 )  THEN
!          NUM = 86
!      ENDIF
!
!     Modified by M. Cawkwell 6th February 2004
!
!       if (mvflag == 16) then
!          if ( key  ==  0 ) then 
!             write(filename,'("ddplot/st.interm.",i0)') iter
!             open(unit = 87,file = filename)
!             num = 87
!          else if ( key  ==  1 )  then
!             num = 86
!             open(unit = 86,file = 'ddplot/dsout.txt')
!          end if
!       end if
! !
!       if (mvflag == 17) then
!          if ( key  ==  0 ) then 
!             write(filename,'("ddplot/st.interm.",i0)') iter
!             open(unit = 87,file = filename)
!             num = 87
!          else if ( key  ==  1 )  then
!             write(filename,'("ddplot/plot.totstr.",F7.5)') totstr
!             open(unit=86, file = filename)
!             num = 86
!          end if
!       end if      
!       



      num = get_new_unit()
      open(num, file=trim(fln), action='write')      

!     And, finally, let's write down the plot file.

      write ( num, '(I4)') total
 
!     Write down the relaxed Z coordinates.

      do i = 1, total
          write ( num, '(F9.5)' )  crd( 3, indpl( i ) )
      end do

!     Write down the relaxed X and Y coordinates and atomic species. 

      do i = 1, total
          ii = indpl( i )
          write ( num, '(F9.5,1X,F9.5,1X,I1)' ) crd( 1, ii ), crd( 2, ii ), ind( ii )
      end do

      write ( num, '(I4)' ) total

!     Write down the unrelaxed Z coordinates.

      do i = 1, total
          write ( num, '(F9.5)' )  ucrd( 3, indpl( i ) )!/latpar
      end do

!     Write down the unrelaxed X and Y coordinates.

      do i = 1, total
          ii = indpl( i )
          write ( num, '(F9.5,1X,F9.5)' )  ucrd(1,ii), ucrd(2,ii)!/latpar /latpar
      end do

!     Finally, this zero is required by the protocol of the plotting routine.

      write ( num, '(I4)') 0
      write ( num, '(E14.7)' ) lena(1)!/latpar
      write ( num, '(E14.7)' ) lena(2)!/latpar
      write ( num, '(E14.7)' ) lena(3)!/latpar

!       call tag( num )

      call close_unit( num )  

      deallocate(indpl)

!       
!       if ( key .eq. 0 ) then
!             write(filename,'("ddplot/st.interm_",f6.4,".",i5.5)')   totstr,iter
!             open(unit = 89,file = filename)
!             num = 89
!       else if ( key .eq. 1 )  then
!          write(filename,'("ddplot/plot_",f6.4,".dat")') totstr
!          open(unit = 88,file = filename)
!          num = 88
!       end if
! 
! 
!    ! *     And, finally, let's write down the plot file.
! 
!       write ( num, '(i4)') total
! 
!    ! *     Write down the relaxed Z coordinates.
! 
!       do i = 1, total
!             write ( num, '(f9.5)' )  crd( 3, indpl( i ) )
!       end do
! 
!    ! *     Write down the relaxed X and Y coordinates and atomic species. 
! 
!       do i = 1, total
!             ii = indpl( i )
!             write ( num, '(f9.5,1x,f9.5,1x,i1)' )  crd( 1:2, ii ), ind( ii )
!       end do
! 
!       write ( num, '(i4)' ) total
! 
! !    *     Write down the unrelaxed Z coordinates.
! 
!       do i = 1, total
!             write ( num, '(f9.5)' ) unrld( 3, indpl( i ))!/latpar
!       end do
! 
! !    *     Write down the unrelaxed X and Y coordinates.
! 
!       do i = 1, total
!             ii = indpl( i )
!             write ( num, '(f9.5,1x,f9.5)' ) unrld(1:2,ii)!/latpar
!       end do
! 
! !    *     Finally, this zero is required by the protocol of the plotting routine.
! 
!       write ( num, '(i4)') 0
!       write ( num, '(e14.7)' ) lena(1)!/latpar
!       write ( num, '(e14.7)' ) lena(2)!/latpar
!       write ( num, '(e14.7)' ) lena(3)!/latpar
! 
!       call tag( num )
! 
!       close( num )  
      
      end subroutine writeddp
