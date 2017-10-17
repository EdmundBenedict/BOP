 
!****************************************************************************
!                                                                           *
!                                                                           *
!                             Alexey Girshick                               *
!              Department of Materials Science and Engineering              *
!                        University of Pennsylvania                         *
!                                                                           *
!                             Rewritten by M.M.                             *
!                             messed up by DLP                              *
!****************************************************************************
!                                                                           *
!                                                                           *
!        This subroutine imposes strain according to the chosen             *
!        applied stress. Usually used for studies of external               *
!                 stresses on dislocations.                                 *
!                                                                           *
!****************************************************************************


   subroutine stress()

      use mod_precision
      use mod_all_scalar
      use mod_const

      implicit none

      include "Include/PosVel.array"
      include "Include/ag.conn"


      integer :: i, j, istr, nstr
      real(dp) :: pmodul, appstr

      character(len=75) :: dummytext
      character(len=38) :: comment
      character(len=80) :: filename
      


!  Read in the stress application information from the channel # 35.

!
!    *****EXAMPLE OF INPUT FILE******
!
!------------------------------------------------------------------
!  APPLIED STRESS: TENSION IN [012] - Mo                   (comment)
!  THE STRESS IS NORMALIZED BY THE SHEAR MODULUS C44=1.089 (comment)
!  1.089                                  (the normalization factor)
!  0.000  0.002  5   (starting stress, increment, no. of increments)
!     0.00000E+00   0.00000E+00   0.00000E+00
!     0.00000E+00   0.40000E+00  -0.48990E+00        (stress tensor)
!     0.00000E+00  -0.48990E+00   0.60000E+00
! 
!    -0.10211E+00   0.29866E-01  -0.12193E-01
!     0.29866E-01   0.64305E-01  -0.18269E+00        (strain tensor)
!    -0.12193E-01  -0.18269E+00   0.16475E+00
!------------------------------------------------------------------


      if (.not. quiet) write ( 6, '(//80("*"))' )
      if (.not. quiet) write ( 6, '(/14X, " --- STRESS APPLICATION --- "/)')

      read ( 35, '(A75)' )  dummytext
      if (.not. quiet) write ( 6, * ) dummytext
      read ( 35, '(A75)' )  dummytext
      if (.not. quiet) write ( 6, * ) dummytext
      read ( 35, * ) pmodul
      if (.not. quiet) write ( 6,'(" The normalization factor = ",F8.4)' ) pmodul
      read ( 35, * ) totstr, appstr, nstr
      
      if (.not. quiet) then
         write ( 6,'(" Starting value of stress = ",F8.4)') totstr
         write ( 6,'(" Stress increment = ",F8.4)') appstr
         write ( 6,'(" Number of stress increments = ",i0)') nstr
      end if

      read ( 35, '(A75)' )  dummytext             
      read ( 35, '(A75)' )  dummytext             
      read ( 35, '(A75)' )  dummytext             
      read ( 35, '(A75)' )  dummytext             

      read ( 35, * )  et( 1, 1), et( 1, 2), et( 1, 3)
      read ( 35, * )  et( 2, 1), et( 2, 2), et( 2, 3)
      read ( 35, * )  et( 3, 1), et( 3, 2), et( 3, 3)

      if (.not. quiet) then
         write ( 6, '(//20X,"Strain tensor in the dislocation basis",/)')
         do i = 1, 3
            write ( 6, '(/11X,3(2X,E16.8))' )  ( et(i,j),j=1,3 )
         enddo
         write ( 6, '(//80("*")/)' )
      end if

      
!  Multiply the strain matrix by the normalization factor.


      et = et * pmodul

!  Memorize the original unshifted coordinates.

      unshad( 1:3, 1:nd ) = ad( 1:3, 1:nd )
      unshind( 1:3, 1:ninert) = adinert( 1:3, 1:ninert)

!
! Start the big cycle over incremental stress application
!

!       nbtflg = 1
!      TOTSTR = 0.0

      do istr = 1, nstr

         iter = 1
         totstr = totstr + appstr
!          the above may introduce error
!         it is safer to use totstr = istr*appstr
         
         if (.not. quiet) then
            write ( 9, '("APPLIED STRESS = ",F8.4)') totstr
            write ( 72, '(/"*****************************"/)')
            write ( 72, '("Step = ",I2)') istr
            write ( 72, '("Applied stress = ",F8.5/)') totstr
         end if

!     Note:
!     All atoms get strained in BLDLIST.F and BLDLISTFS.F 
!     so we don't do anything here.

!     Relax the dislocation.

         call relax_ds()

!         WRITE ( 32, '(/"*****************************"/)')
!         WRITE ( 32, '("Step = ",I2)') ISTR
!         WRITE ( 32, '("Applied stress = ",F8.4/)') TOTSTR

!  Put back the original unshifted coordinates.
!         PRINT *, LENA(1), LENA(2), LENA(3)
!
         if (gfbcon  ==  1) then
            call gfbc_outp()
         endif

!         IF (APPSTR .LT. 0.0D0) THEN
          if (.not. quiet) call outblock(1)            
!         ENDIF


         ad(1:3, 1:nd) = unshad(1:3, 1:nd)
         adinert(1:3, 1:ninert) = unshind(1:3, 1:ninert)

         
         if (.not. quiet) then
            write(filename,'("st.rlx.",f5.3)') totstr
         
            call writeddp('ddplot/'//trim(filename)//'.dat')
            call writexyz('xyz/'//trim(filename)//'.xyz', filename)         
            call writecfg('cfg/'//trim(filename)//'.cfg')
            call writexbs('bs/'//trim(filename)//'.bs')            
            call writecell('cell/'//trim(filename)//'.cell.in')
         end if
         
!         IF (APPSTR.GE.0.0D0) THEN
!            CALL OUTBLOCK(1)
!         ENDIF
!
!         IF (GFBCON .EQ. 1) THEN
!            CALL GFBC_OUTP()
!         ENDIF
!
      enddo

   end subroutine stress


