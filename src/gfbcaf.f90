 
      subroutine gfbcaf(radiusi)
      use mod_precision
      use mod_all_scalar, only : iter, nd
      use mod_const
      use topologia
          
      implicit none
!
!     This subroutine calculates the average forces in regions I and II
!     during relaxation of a dislocation using Greens function boundary
!     conditions. This should help you see whether the relaxation is
!     proceeding correctly with forces being built up in region II 
!     region I is relaxing. This is also useful to check that the 
!     position update is doing the correct thing.
!
!     M. Cawkwell 25th September 2003
!

      include "Include/PosVel.array"
      include "Include/Force.array"

      integer i, j, regic, regiic
      real(dp) :: radiusi, avfi(3), avfii(3), avfti, avftii
      real(dp) :: dfromo
!
     if (iproc == master) open(unit=22, status="new", file="gfbc.avefi.dat")
     if (iproc == master) open(unit=23, status="new", file="gfbc.avefii.dat")
!
      do j = 1,3
         avfi(j) = 0.0_dp
         avfii(j) = 0.0_dp
      enddo
      avfti = 0.0_dp
      avftii = 0.0_dp
      regic = 0
      regiic = 0
!
      do i = 1,nd
         dfromo = sqrt(ad(1,i)*ad(1,i) + ad(2,i)*ad(2,i))
!
         if (dfromo  <=  radiusi) then
!
            regic = regic + 1
!             do j = 1,3
!                avfi(j) = avfi(j) + sqrt(ftot(j,i)*ftot(j,i)) ! This is fairly dumb
!             enddo
            avfi(1:3) = avfi(1:3) + abs(ftot(1:3,i))
            avfti = avfti + sqrt(sum(ftot(1:3,i)*ftot(1:3,i)))
!
         endif
!
         if (dfromo  >  radiusi) then
!
            regiic = regiic + 1
            avfii(1:3) = avfii(1:3) + abs(ftot(1:3,i))
            avftii = avftii + sqrt(sum(ftot(1:3,i)*ftot(1:3,i)))
            
!             do j = 1,3
!                avfii(j) = avfii(j) + sqrt(ftot(j,i)*ftot(j,i)) 
!             enddo
!             avftii = avftii + sqrt(ftot(1,i)*ftot(1,i) + & 
!      &           ftot(2,i)*ftot(2,i) + ftot(3,i)*ftot(3,i))           
! !
         endif
!
      enddo
!
      avfi(1:3) = avfi(1:3)/real(regic,kind=dp)
      avfii(1:3) = avfii(1:3)/real(regiic,kind=dp)
      avfti = avfti/real(regic,kind=dp)
      avftii = avftii/real(regiic,kind=dp)
!
      if (iproc == master) write(22,10) iter, avfi(1), avfi(2), avfi(3), avfti
      if (iproc == master) write(23,10) iter, avfii(1), avfii(2), avfii(3), avftii
!
 10   format(1x,i4,4(1x,f14.10))
!
      
      end



               
