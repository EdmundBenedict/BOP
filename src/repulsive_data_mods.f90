
      
    module pair_coefficients
        use mod_precision
        use mod_conf
        implicit none
        
!         type(pwp_t), allocatable :: ppot(:)
        type(pwp_t) :: ppot(0:2)
                
    end module pair_coefficients


                  
                  
! !                   Tih%base = (/ 0.5190_dp, 0.5246_dp /)   ! znam's bond part
! !                   Tih%base = (/ 0.1295_dp, 0.1625_dp /)     ! my bond parts
! !                   Tih%base = (/ 0.3595_dp, 0.3922_dp /)
!                   tih%base = (/ 0.1687_dp, 0.2034_dp /)
!                   tih%xp   = (/ 0.1090_dp, 0.2610_dp /)
!                   tih%e    = 0.0_dp
!                   tih%base_e    = 0.0_dp
!                   tih%trgt = tih%xp - tih%base
!                   
! !                   L10%base = (/  0.3388_dp,  0.4684_dp /) ! znam's bond part
! !                   L10%base = (/  0.3097_dp,  0.6361_dp /)   ! my bond parts
! !                   L10%base = (/  0.3605_dp,  0.6619_dp /)
! !                  L10%base = (/  0.3466_dp,  0.5991_dp /)
! !                  L10%base = (/  0.4693_dp, 0.4962_dp /)
! !                  L10%base = (/  0.0350_dp,  0.2515_dp /)
!                   l10%base = (/ 0.961200_dp, 0.836300_dp /)
!                   l10%xp   = (/ -0.2130_dp, -0.0400_dp /)
!                   l10%e    =  -4.52_dp
! !                  L10%base_e    =  -8.263189237_dp
! !                  L10%base_e    =  -7.566772945_dp
!                   l10%base_e    =  -7.680525_dp
!                   l10%trgt = l10%xp - l10%base
!                   
! !                   DO19%base =  (/ 0.5852_dp, 0.7755_dp/)  ! znam's bond part
! !                   DO19%base =  (/ 0.4050_dp, 0.3880_dp/)       ! my bond parts
! !                   DO19%base =  (/ 0.5635_dp, 0.5914_dp/)
! !                   DO19%base =  (/ 0.4284_dp, 0.4093_dp/)
!                   do19%base =  (/ 0.6445_dp, 0.8954_dp /)
!                   do19%xp   =  (/-0.0125_dp, 0.2620_dp/)
!                   do19%e    =  -4.75_dp
!                   do19%base_e = -7.8245_dp
!                   do19%trgt = do19%xp - do19%base
!                   
