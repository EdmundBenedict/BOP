    subroutine elcon_hcp_short()
        use mod_precision
        use mod_const
        use mod_all_scalar
        use mod_conf, only : rtc, ens_t
        
        implicit none
        
!         type(run_conf_t), intent(inout) :: rtconf
!         type(elastic_const_t), pointer :: elcons => null()

        integer, parameter :: strains = 6

        include "Include/Atom.array"
        include "Include/BEC.array"
        include "Include/PosVel.array"
        include "Include/ag.conn"


!        the curvature, slope, becp and the sxx, cxx have I dim = 4 to hold the tot, elc, env, pwp contributions

        real(dp) :: strain((strains),3,3), volume, curvature(4,strains), slope(4,strains)
        real(dp), dimension(4) :: sigma11, sigma33, c11, c12, c13, c33, c44, c66, c1344, c1266
        character(len=20) :: fileext(strains)
        integer :: i, j, k, flag, it
        real(dp), parameter :: zero = 0.0_dp, one = 1.0_dp, half = 0.5_dp
        real(dp) :: ta(3,3), a1(3,3), dalpha, eflo0
        real(dp),allocatable :: de0(:),de1(:)
        type(ens_t) :: en0, aven0
        real(dp) :: becp(4,nt)

!         real(dp) :: ad0(3,size(ad,2))
        real(dp) ::  dlim(2,3)!, unitm(3,3)
        
!         unitm(:,1) = (/ one , zero, zero /)
!         unitm(:,2) = (/ one ,  one, zero /)
!         unitm(:,3) = (/ one , zero,  one /)



!         elcons => rtconf % elastic_const_conf


        fileext( 1 ) = 'c11'
        strain( 1, 1, 1:3 ) = (/ one , zero, zero /)
        strain( 1, 2, 1:3 ) = (/ zero, zero, zero /)
        strain( 1, 3, 1:3 ) = (/ zero, zero, zero /)

        fileext( 2 ) = 'c33'
        strain( 2, 1, 1:3 ) = (/ zero, zero, zero /)
        strain( 2, 2, 1:3 ) = (/ zero, zero, zero /)
        strain( 2, 3, 1:3 ) = (/ zero, zero, one  /)

        fileext( 3 ) = '2c11+2c12'      
        strain( 3, 1, 1:3 ) = (/ one , zero, zero /)
        strain( 3, 2, 1:3 ) = (/ zero, one , zero /)
        strain( 3, 3, 1:3 ) = (/ zero, zero, zero /)

        fileext( 4 ) = 'c11+c33+2c13'
        strain( 4, 1, 1:3 ) = (/ one , zero, zero /)
        strain( 4, 2, 1:3 ) = (/ zero, zero, zero /)
        strain( 4, 3, 1:3 ) = (/ zero, zero,  one /)

        fileext( 5 ) = '4c44'      
        strain( 5, 1, 1:3 ) = (/ zero, zero, one  /)
        strain( 5, 2, 1:3 ) = (/ zero, zero, zero /)
        strain( 5, 3, 1:3 ) = (/ one , zero, zero /)

        fileext( 6 ) = '4c66'
        strain( 6, 1, 1:3 ) = (/ zero, one , zero /)
        strain( 6, 2, 1:3 ) = (/ one , zero, zero /)
        strain( 6, 3, 1:3 ) = (/ zero, zero, zero /)

!         if (.not. quiet) open( unit = 80, file = 'ec_short.out' ) 

        if (nt /= 3 .and. nt /= 5) then
            write(6, '("This is the short version. It is hardcoded to 3 and 5 points only.")')
            call panic()
        end if
        
        
        dalpha = (thi - tlo)/real(nt-1, dp)
        talpha(1) = tlo

        do it=1,nt
            talpha(it+1) = talpha(it) + dalpha
        end do
        
        eflo = 0.0_dp
        de(:) = 0.0_dp
!         call bldnebtfs()
        flag = 1
        call getetot(flag)
        
        allocate(de0(nd),de1(nd))
        de0 = de(:nd)
        eflo0 = eflo

        en0 = rtc % tote
        aven0 = rtc % avre

        a1(:,:) = a(:,:)  

!         ad0(:,:) = ad(:,:)
!          call dumpneighl(i,1,'zero')
        do i = 1, strains
            if (.not. quiet) print '(a)', fileext(i)
            tmatrix(:,:) = strain(i,:,:)
            ta = matmul(a1, tmatrix)
            
            if (nt == 3 .and. talpha(2) == 0.0_dp) then
            
               a(:,:) = a1(:,:) + talpha(1)*ta(:,:)
               ad(:,:nd) = matmul(transpose(a(:,:)), d(:,:nd))
!                write(6,"('original:',16(/,3(x,f12.6)))") ad(:,:nd)
!                ad(:,:) = matmul(unitm + talpha(1)*transpose(tmatrix), ad0(:,:))
!                write(6,"('nov:',32(/,3(x,f12.6)))") ad(:,:nd*2)
!                call bldlistfs(dlim)
!                write(6,"('list:',32(/,3(x,f12.6)))") ad(:,:nd*2)
!                stop
               adinert(:,:ninert) = matmul(transpose(a(:,:)), dinert(:,:ninert))
!                print *, 'ninert:', ninert
               de(:nd) = de0
               eflo = eflo0
               flag = 1

!                call dumpneighl(i,1,'pred')
!                call bldlistfs(dlim)
               call getetot(flag)
!                bec(1) = eprom + ebond + epair - eatom - uent

!                call dumpneighl(i,1,'sled')



               becp(1,1) = rtc % tote % bind
               becp(2,1) = rtc % tote % elct
               becp(3,1) = rtc % tote % env
               becp(4,1) = rtc % tote % pair

               de1 = de(:nd)
               a(:,:) = a1(:,:) + talpha(3)*ta(:,:)
               ad(:,:nd) = matmul(transpose(a(:,:)), d(:,:nd))
               adinert(:,:ninert) = matmul(transpose(a(:,:)), dinert(:,:ninert))
               de(:nd) = de0 + de0 - de1
               eflo = eflo0
               flag = 1
!                call bldlistfs(dlim)
               call getetot(flag)
!                bec(3) = eprom + ebond + epair - eatom - uent
!                bec(2) = etot0

               becp(1,3) = rtc % tote % bind
               becp(2,3) = rtc % tote % elct
               becp(3,3) = rtc % tote % env
               becp(4,3) = rtc % tote % pair

               becp(1,2) = en0 % bind
               becp(2,2) = en0 % elct
               becp(3,2) = en0 % env
               becp(4,2) = en0 % pair                    

                    
            else    
!             print *, redb,'wrong',endc
               do it = 1,nt,1
                  if (talpha(it) /= 0.0_dp) then
                     if (.not. quiet) print '(4x,f12.6)', talpha(it)
                     a(:,:) = a1(:,:) + talpha(it)*ta(:,:)
                     ad(:,:nd) = matmul(transpose(a(:,:)), d(:,:nd))
                     adinert(:,:ninert) = matmul(transpose(a(:,:)), dinert(:,:ninert))
                     
   !                     print *, 'predi de: ', de(:nd)
   !                     eflo = 0.0_dp
   !                     de(:) = 0.0_dp
                     de(:nd) = de0
                     eflo = eflo0
                     flag = 1
                     call getetot(flag)
   !                     print *, 'de: ', de(:nd)
                     
!                      bec(it) = eprom + ebond + epair - eatom - uent
                     becp(1,it) = rtc % tote % bind
                     becp(2,it) = rtc % tote % elct
                     becp(3,it) = rtc % tote % env
                     becp(4,it) = rtc % tote % pair
                  else
!                      bec(it) = etot0
                     becp(1,it) = en0 % bind
                     becp(2,it) = en0 % elct
                     becp(3,it) = en0 % env
                     becp(4,it) = en0 % pair
                  end if
               enddo
            end if

            select case (nt)
                case (3)
                    slope(:,i) = 0.5_dp * (becp(:,3) - becp(:,1)) / dalpha
                    curvature(:,i) = (becp(:,3) - 2*becp(:,2) + becp(:,1))/(dalpha*dalpha)
                case (5)
!                 The linear five-point stencil
                    slope(:,i) = (becp(:,1) - 8*becp(:,2) + 8*becp(:,4) - becp(:,5)) / (12*dalpha)
                    curvature(:,i) = ( -becp(:,1) + 16*becp(:,2) - 30*becp(:,3) + 16*becp(:,4) - becp(:,5) ) &
                                   & / (12 * dalpha * dalpha)
                case default
                    write(6, '("This is the short version. It is hardcoded to 3 and 5 points only.")')
                    call panic()
            end select
            
!             if (.not. quiet) then
!                 write(  6,'(5x, 4(x,f20.12), 8x, a20)' ) curvature(:, i ), fileext( i )
!                 write( 80,'(5x, 4(x,f20.12), 8x, a20)' ) curvature(:, i ), fileext( i )
!             end if
        enddo

        deallocate(de0,de1)
         
        
         
        c11(:) = curvature(:,1)
        c33(:) = curvature(:,2)
        c12(:) = half * curvature(:,3) - curvature(:,1)
        c13(:) = half * ( curvature(:,4) - ( curvature(:,1) + curvature(:,2) ) )
        c44(:) = 0.25_dp * curvature(:,5)
        c66(:) = 0.25_dp * curvature(:,6)      

!         print *,'a:', a(:,:)
        a(:,:) = a1(:,:)        
        ad(:,:nd) = matmul(transpose(a(:,:)), d(:,:nd))
        adinert(:,:ninert) = matmul(transpose(a(:,:)), dinert(:,:ninert))

        volume = (a(1,2)*a(2,3) - a(1,3)*a(2,2))*a(3,1) &
             & + (a(1,3)*a(2,1) - a(1,1)*a(2,3))*a(3,2) &
             & + (a(1,1)*a(2,2) - a(1,2)*a(2,1))*a(3,3)

        sigma11(:) = slope( :, 1 ) / volume 
        sigma33(:) = slope( :, 2 ) / volume

        c1344 = (c13 / volume) - (c44 / volume - ( sigma11+sigma33 )*0.25_dp)
        c1266 = (c12 / volume) - (c66 / volume - sigma11 * 0.5_dp)


!         print *, associated(elcons)
!         print *, elcons % s11, elcons % s33

        rtc % tote = en0
        rtc % avre = aven0

        rtc % elastic_const_conf % s11 % tot = sigma11(1)
        rtc % elastic_const_conf % s11 % elc = sigma11(2)
        rtc % elastic_const_conf % s11 % env = sigma11(3)
        rtc % elastic_const_conf % s11 % pwp = sigma11(4)

        rtc % elastic_const_conf % s33 % tot = sigma33(1)
        rtc % elastic_const_conf % s33 % elc = sigma33(2)
        rtc % elastic_const_conf % s33 % env = sigma33(3)
        rtc % elastic_const_conf % s33 % pwp = sigma33(4)

        rtc % elastic_const_conf % c11 % tot = c11(1)/volume
        rtc % elastic_const_conf % c11 % elc = c11(2)/volume
        rtc % elastic_const_conf % c11 % env = c11(3)/volume
        rtc % elastic_const_conf % c11 % pwp = c11(4)/volume

        rtc % elastic_const_conf % c33 % tot = c33(1)/volume
        rtc % elastic_const_conf % c33 % elc = c33(2)/volume
        rtc % elastic_const_conf % c33 % env = c33(3)/volume
        rtc % elastic_const_conf % c33 % pwp = c33(4)/volume

        rtc % elastic_const_conf % c12 % tot = c12(1)/volume
        rtc % elastic_const_conf % c12 % elc = c12(2)/volume
        rtc % elastic_const_conf % c12 % env = c12(3)/volume
        rtc % elastic_const_conf % c12 % pwp = c12(4)/volume

        rtc % elastic_const_conf % c13 % tot = c13(1)/volume
        rtc % elastic_const_conf % c13 % elc = c13(2)/volume
        rtc % elastic_const_conf % c13 % env = c13(3)/volume
        rtc % elastic_const_conf % c13 % pwp = c13(4)/volume

        rtc % elastic_const_conf % c44 % tot = c44(1)/volume
        rtc % elastic_const_conf % c44 % elc = c44(2)/volume
        rtc % elastic_const_conf % c44 % env = c44(3)/volume
        rtc % elastic_const_conf % c44 % pwp = c44(4)/volume

        rtc % elastic_const_conf % c66 % tot = c66(1)/volume
        rtc % elastic_const_conf % c66 % elc = c66(2)/volume
        rtc % elastic_const_conf % c66 % env = c66(3)/volume
        rtc % elastic_const_conf % c66 % pwp = c66(4)/volume

        rtc % elastic_const_conf % c1344 % tot = c1344(1)
        rtc % elastic_const_conf % c1344 % elc = c1344(2)
        rtc % elastic_const_conf % c1344 % env = c1344(3)
        rtc % elastic_const_conf % c1344 % pwp = c1344(4)

        rtc % elastic_const_conf % c1266 % tot = c1266(1)
        rtc % elastic_const_conf % c1266 % elc = c1266(2)
        rtc % elastic_const_conf % c1266 % env = c1266(3)
        rtc % elastic_const_conf % c1266 % pwp = c1266(4)

! 
!         if (.not. quiet) then
!                 
!             write( 80,'(//"compare the solutions"/)' )
! 
!             write( 80,'(5x, 4(x,f20.12), 8x, 4(x,f20.12))' ) curvature(:, 1 ), c11
!             write( 80,'(5x, 4(x,f20.12), 8x, 4(x,f20.12))' ) curvature(:, 2 ), c33
!             write( 80,'(5x, 4(x,f20.12), 8x, 4(x,f20.12))' ) curvature(:, 3 ), 2.0_dp * ( c12 + c11 )
!             write( 80,'(5x, 4(x,f20.12), 8x, 4(x,f20.12))' ) curvature(:, 4 ), c11 + c33 + 2.0_dp * c13
!             write( 80,'(5x, 4(x,f20.12), 8x, 4(x,f20.12))' ) curvature(:, 5 ), 4.0_dp * c44
!             write( 80,'(5x, 4(x,f20.12), 8x, 4(x,f20.12))' ) curvature(:, 6 ), 4.0_dp * c66
! 
!             write( 80,'(/5x,"volume = ", f10.6/)' )  volume 
! 
!             write( 80,'( 62("_")//19x,"in ev/(a**3)",12x,"in 10**11 pa"// &
!             &             17x,"calc.",10x,"calc.",/, &
!             &             62("_")/ )' )
! 
!             write( 80,'(5x,"c11",8x, 4(x,f20.12), 8x, 4(x,f20.12) )' ) &
!             &        c11/volume, c11*1.602_dp/volume
! 
!             write( 80,'(5x,"c12",8x, 4(x,f20.12), 8x, 4(x,f20.12) )' ) &
!             &        c12/volume, c12*1.602_dp/volume
! 
!             write( 80,'(5x,"c13",8x, 4(x,f20.12), 8x, 4(x,f20.12) )' ) &
!             &        c13/volume, c13*1.602_dp/volume
! 
!             write( 80,'(5x,"c33",8x, 4(x,f20.12), 8x, 4(x,f20.12) )' ) &
!             &        c33/volume, c33*1.602_dp/volume
! 
!             write( 80,'(5x,"c44",8x, 4(x,f20.12), 8x, 4(x,f20.12) )' ) &
!             &        c44/volume, c44*1.602_dp/volume
! 
!             write( 80,'(5x,"c66",8x, 4(x,f20.12), 8x, 4(x,f20.12) )' ) &
!             &        c66/volume, c66*1.602_dp/volume
! 
! ! 
! !             !  This part is for writing out an input file for the fitting program.
! ! 
! !             write ( 80,'(//"sigma11 =",f16.12/)' ) sigma11
! !             write ( 80,'("sigma33 =",f16.12/)' ) sigma33 
! !             write ( 80,'("de/da =",f16.12/)' ) slope( 1 ) / a( 1, 1 ) / 4.d0
! !             write ( 80,'("de/dc =",f16.12/)' ) slope( 2 ) / a( 3, 3 ) / 4.d0
! ! 
! !             write ( 80,'(65("_"))' )
! !             write ( 80,'("ecoh_bond"/)' )
! !             write ( 80,'("deda_bond")' )
! !             write ( 80,'(f18.12)' )  slope( 1 ) / a( 1, 1 ) / 4.d0
! !             write ( 80,'("dedc_bond")' )
! !             write ( 80,'(f18.12)' )  slope( 2 ) / a( 3, 3 ) / 4.d0
! !             write ( 80,'("sig11_bond")' )
! !             write ( 80,'(f18.12)' )  sigma11
! !             write ( 80,'("sig33_bond")' )
! !             write ( 80,'(f18.12)' )  sigma33                      
! !             write ( 80,'("c11_bond")' )
! !             write ( 80,'(f18.12)' )  c11 / volume - sigma11
! !             write ( 80,'("c12_bond")' )
! !             write ( 80,'(f18.12)' )  c12 / volume
! !             write ( 80,'("c13_bond")' )
! !             write ( 80,'(f18.12)' )  c13 / volume
! !             write ( 80,'("c33_bond")' )
! !             write ( 80,'(f18.12)' )  c33 / volume - sigma33
! !             write ( 80,'("c44_bond")' )
! !             write ( 80,'(f18.12)' )  c44 / volume - ( sigma11 + sigma33 )/4.d0
! !             write ( 80,'("c66_bond")' )
! !             write ( 80,'(f18.12)' )  c66 / volume - sigma11 / 2.d0
! ! 
! !             write ( 80,'(///, 16x, "cauchy  pressures", 8x, "calc.", 8x, "exp." )' )
! !             write ( 80,'(/ 20x, "c13 - c44", 10x, f8.4 )' ) &
! !             &   (c13 / volume) - (c44 / volume - ( sigma11+sigma33 )/4.d0)
! ! 
! !             write ( 80,'(/, 20x, "c12 - c66", 10x, f8.4 /)' ) &
! !             &   (c12 / volume) - (c66 / volume - sigma11 / 2.d0)
! ! 
! !             write ( 6,'(///, 16x, "cauchy  pressures", 8x, "calc.", 8x, "exp.")' )
! !             write ( 6,'(/ 20x, "c13 - c44", 10x, f8.4 )' ) &
! !             &   (c13 / volume) - (c44 / volume - ( sigma11+sigma33 )/4.d0)
! !             write ( 6,'(/20x,"c12 - c66",10x,f8.4,7x,f6.3/)' ) &
! !             &   (c12 / volume) - (c66 / volume - sigma11 / 2.d0)
! ! 
! !             write(6, *) 'epair', epair, epair/real(nd,dp)
! ! 
!             call tag( 80 )
! 
!             close( 80 )
!         end if
      contains
         subroutine dumpneighl(s,p,sfx)
               use mod_io, only : get_new_unit, close_unit
               implicit none

               include "Include/NebList.array"
               include "Include/NebListL.array"
               include "Include/ag.conn"
               include "Include/PosVel.array"

               integer, intent(in) :: s,p
               character(len=4), intent(in) :: sfx

               integer :: i, j, ja, u, mj
               character(len=16) :: fn

               mj = 0

               write(fn,"('sdata_s',i1,'_p',i1,'_',a)") s, p, sfx
               u = get_new_unit()
               open(u, file=fn, action='write')
               write(u,'("nne:",i0)') nne
               write(u,'("neighs:")',advance='no')
               do i = 1, nne
                  write(u,'(/,i6,":")',advance='no') i 
                  ja = aptrl(i)
                  j  = bptrl(ja)
                  do while (j /= eol)
                     write(u,'(x,i0)',advance='no') j
                     mj = max(mj,j)
                     ja = ja + 1
                     j = bptrl(ja)
                  end do
               end do          
               
               write(u,'(/,"ad:")')
               do i=1,mj
                  write(u,"(3(x,f12.6))") ad(:,i)
               end do
               
               call close_unit(u)
         end subroutine dumpneighl
         


    end subroutine elcon_hcp_short

