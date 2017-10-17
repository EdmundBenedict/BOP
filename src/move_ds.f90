
   !****************************************************************************
   !                                                                           *
   !                                                                           *
   !                             Alexey Girshick                               *
   !              Department of Materials Science and Engineering              *
   !                        University of Pennsylvania                         *
   !                                                                           *
   !                                 +DLP                                      *
   !****************************************************************************
   !                                                                           *
   !                                                                           *
   !     This is an atom moving subroutine for the dislocation relaxation.     *
   !                                                                           *
   !                                                                           *
   !****************************************************************************



   subroutine move_ds()
      use mod_precision
      use mod_all_scalar
      use mod_const

      implicit none

      include "Include/PosVel.array"
      include "Include/Force.array"
      include "Include/ag.conn"

      real(dp) ::  ok_new
      real(dp) :: dis( 3 ), totdis( 3 ), displ
      real(dp) :: disorg
      integer  :: i, j, k, seed
      integer  :: stat, kk
      real(dp) :: invet1(3,3), v3(3)
!       real(dp) :: rand
      

      real, parameter :: maxdad = 0.20_dp
      
      if (mvflag == 17) then
!       DLP:  dis corresponds to strained coords. to avoid force oscilation at higher stresses it shall be made consistent by multiplying with (et+I)^(-1)
!       so that  (ad + dis) = (et+I) (unshad + (et+I)^(-1) dis) 
         invet1 = et
         invet1(1,1) = invet1(1,1) + 1.0_dp
         invet1(2,2) = invet1(2,2) + 1.0_dp
         invet1(3,3) = invet1(3,3) + 1.0_dp
         call inv3x3(invet1)
      end if

   !  Added random displacement for first 10 iterations

   !      SEED = ITER
   !      CALL SRAND(SEED)
   !      DO I=1,5
   !       WRITE(*,*) RAND()
   !      ENDDO
   !
   !     Determining whether the GFBCs have been switched on
   !
   !      CALL FNDREC(8,'GFBCON',STAT)
   !      READ(8,*) GFBCON
   !
   !     If they have been switched on, then find the radius of the region
   !     which is to be relaxed atomistically
   !     
      if (gfbcon  ==  1) then
   !     Call subroutine to check the forces in regions I and II
         call gfbcaf(radiusi)
      endif
   !
      
   !  Displace the atoms according to the individual atomic forces
      do  i = 1, nd

   !      IF ( ITER .LT. 10 ) THEN
   !          OK_NEW = OK * (RAND()+1)
   !          WRITE(*,*) ITER, OK_NEW
   !          DO J = 1, 3
   !              DIS( J ) = FTOT( J, I ) * OK_NEW
   !          ENDDO
   !      ELSE
   !          DO J = 1, 3
   !              DIS( J ) = FTOT( J, I ) * OK
   !          ENDDO
   !      ENDIF
   !
   !     Modified by M. Cawkwell 23rd September 2003 for so GFBCs can be
   !     used.
   !     
   !     GFBCs require that only the active atoms in region I are relaxed
   !     using atomistics, so we will set the displacements for active
   !     atoms at distance > R1 from the origin to zero. R1 must be read
   !     in from the input file fort.8 and flag GFBCON must be set = 1.
   !     Make sure that when making blocks for dislocation relaxations 
   !     the radius of your active region is R1+R2, and that R2-R1 is
   !     greater than 2R_cut for your potential
   !
         dis = ftot( 1:3, i ) * ok

            if (gfbcon == 1) then
               if (sum(ad(1:2,i)*ad(1:2,i))  >  radiusi*radiusi) dis = 0.0_dp
            endif

            displ = sqrt(sum(dis*dis))
            if ( displ  >  maxdad ) dis = maxdad * dis / displ

   ! If we applied also stress we have to add the displacement only
   ! to the unstressed coordinates.

            if (mvflag == 17) then
               call mul3x3(invet1,dis,v3)
               unshad( 1:3, i ) = unshad( 1:3, i ) + v3
            else
               ad( 1:3, i ) = ad( 1:3, i ) + dis
            endif
   !           lef = lef + sum(dis * dmdl( 1:3, i ))
      enddo

   !
   !  If this is stress application, put back the unstressed coordinates.
   !

      if (mvflag == 17) then
         ad( :, 1:nd ) = unshad( :, 1:nd )
         adinert( :, 1:ninert ) = unshind( :, 1:ninert )
      endif

   end subroutine move_ds

