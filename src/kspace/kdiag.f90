
   subroutine kdiag(gpt, knh, enk, kpsi,kovl)
      use mod_precision
      use mod_const
!       use mod_all_scalar
      use mod_io
      use mod_conf
      implicit none

      logical, intent(in) :: gpt
      integer, intent(in) :: knh
      complex(dp), intent(inout) :: kpsi(knh,knh)
      real(dp), intent(out) :: enk(knh)

      integer :: zh_n, zh_lwork, zh_lrwork, zh_liwork, info
      integer, allocatable :: zh_iwork(:)
      real(dp), allocatable :: zh_rwork(:), rkpsi(:,:)
      complex(dp), allocatable :: zh_work(:)
      
      integer, save :: zh_fail = 0
      logical :: ovl
      integer :: ik,i,j
      complex(dp) :: opsi(knh,knh)
      complex(dp), intent(in) :: kovl(knh,knh)
      
!        call print_c(kpsi(:,:),6,'','f8.3','kh')
!        stop
      
      ovl = rtc % ham_conf % ovl
      
      zh_n = knh
      zh_lwork = 1 + 25*knh + 3*knh*knh
      zh_lrwork = 1+4*knh+3*knh*knh+2*knh*20
      zh_liwork = 3 + 5*knh

      allocate(zh_rwork(zh_lrwork), zh_work(zh_lwork), zh_iwork(zh_liwork))
   
   

      if (gpt) then

         allocate(rkpsi(zh_n,zh_n))
         
         rkpsi = real(kpsi)
         
         call dsyevd('V','U',zh_n,rkpsi,knh,enk, & 
&                  zh_work,zh_lwork,zh_iwork,zh_liwork,info)
         if (info /= 0) then
            write(6,'(''Error in DSYEVD: '',I5)') info
            call panic( )
         endif
         kpsi = cmplx(rkpsi,0.0_dp, kind=dp)
         deallocate(rkpsi)
      else
               

!
!       Use ZHEEVD by default but not if it failed last time
!
!          zh_fail = 1
         if (zh_fail == 0) then
!             call print_c(kpsi(:,:),6,'','f8.3','bs')

!             write(1003, "('hr',10(/,10(x,f8.3)))") real(kpsi(:10,:10))
!             write(1003, "('hi',10(/,10(x,f8.3)))") aimag(kpsi(:10,:10))
            if (.not. ovl) then
                call zheevd('V', 'U', zh_n, kpsi, zh_n,enk, zh_work, &
                & zh_lwork, zh_rwork, zh_lrwork, zh_iwork, zh_liwork, info)
            else
                opsi(:,:) = kovl(:,:)
                call zhegvd(1, 'V', 'U', zh_n, kpsi, zh_n, opsi, zh_n, enk, zh_work, &
                & zh_lwork, zh_rwork, zh_lrwork, zh_iwork, zh_liwork, info)
            endif
!                
!             write(1003,"('e',10(x,f8.3))") enk(:10)

!             call print_c(kpsi(:,:),6,'','f8.3','ps')
             call print_c(enk,6,'','f8.3', 'e:') 
!             write(6,'("e:",x,i2,36(x,f6.3))') enk
!             write(6,'(4(/,9(x,f7.4)))') enk
!             stop
         else
            print *, 'zheev'
            call zheev('V','U',zh_n,kpsi,zh_n, & 
   &                 enk,zh_work,zh_lwork,zh_rwork,info)

         endif


   !
   !       Check that the diagonaliser worked.
   !
         if (info > 0) then
            if (zh_fail == 0) then
               write(6,'("ZHEEVD failed to converge; INFO = ", I5)') info
               zh_fail = 1
!                   goto 1
            else
               write(6,'("ZHEEV also failed to converge! INFO = ", I5)') info
               call panic( )
            endif
         elseif (info < 0) then
            write(6,'("ZHEEV(D) was passed an illegal number! INFO = ",I5)') info
            call panic( )
         endif

!
!       Everything worked - reset diagonaliser to ZHEEVD
!
         if (zh_fail == 1) then
            write(6,'("ZHEEV worked.")')
            zh_fail = 0
         endif



      endif
!          write(6,'("done")')


      deallocate(zh_rwork, zh_work, zh_iwork)
      
      
      
      
   end subroutine kdiag
