      
!      based of Matous' routines by the same name
      
      subroutine writecfg(filename)
! c      subroutine writecfg(filename,center)
! c
! c     creates cfg file
      
      use mod_const
      use mod_all_scalar, only : nd, ninert, mvflag
      use mod_atom_ar, only: mass, mg, demi
      
      implicit none
      
!       for istn and zinert only
      include 'Include/Atom.array'
      
      include 'Include/PosVel.array'
      include 'Include/NebList.array'
      include 'Include/Force.array'
      include 'Include/Misc.array'
      include 'Include/ag.conn'
      
      character(len=*) :: filename
      
      real(dp) :: ff, adc, ccmin, ccmax
      
!       real(dp) :: adx,ady,adz
!       real(dp) :: xxmax,yymax,zzmax
!       real(dp) :: xxmin,yymin,zzmin


! c*** a is the array of primitive translation vectors, and b is the inverse.
      real(dp) :: b(3,3),atmp(3,3)

      real(dp) :: pos(3)
      real(dp) :: ld(3),adstra(3),tstra(3)

      integer i,j,num,c
      

! c      real(dp), allocatable ::  tstra(:,:)
! c      real(dp), allocatable ::  tstrin(:,:)
      real(dp), allocatable ::  cxa(:,:)
      real(dp), allocatable ::  cxi(:,:)


! c      allocate (tstra(3,nd))
! c      allocate (tstrin(3,ninert))
      allocate (cxa(3,nd))
      allocate (cxi(3,ninert))


      num=33
      open (unit=num,status='unknown',form='formatted',file=trim(filename))

      xxmax=0.0
      yymax=0.0
      zzmax=0.0
      xxmin=0.0
      yymin=0.0
      zzmin=0.0



! c
! c     handle correctly pbc
! c
! c     get the safe size of the block in the c-direction
! c

      atmp = 0.0_dp
      
      do c = 1, 3
         ccmin = 0.0_dp
         ccmax = 0.0_dp
         
         if (rlim(c) == 0) then
            do i = 1,nd
               if (mvflag == 17) then
                  adc = unshad(c,i) + sum( et(c,1:3) * unshad(1:3,i))*totstr
               else
                  adc = ad(c,i)
               endif

               if (adc < ccmin) ccmin = adc
               if (adc > ccmax) ccmax = adc
            enddo
            
            do i = 1,ninert
               if (mvflag == 17) then
                  adc = unshind(c,i) + sum(et(c,1:3) * unshind(1:3,i))*totstr
               else
                  adc = adinert(c,i)
               endif
               if (adc < ccmin) ccmin = adc
               if (adc > ccmax) ccmax = adc
            enddo
            atmp(c,c)=ccmax-ccmin+10.0_dp
         else
            atmp(c,1:3)=a(c,1:3)
         endif
      end do
      

! c     get inverse matrix of coordinate system
      b = atmp
      call inv3x3(b)

! c
! c     if there is external stress add it to coordinates
! c
      if (mvflag == 17) then

         do i = 1, nd
            do j=1,3
               tstra(j) = sum( et(j,1:3) * unshad(1:3,i)) * totstr
            enddo

            do j=1,3
               adstra(j) = unshad(j,i) + tstra(j)
            enddo

            call mul3x3(b,adstra,ld)

            cxa(1:3,i) = ld
            where (ld < 0.0_dp) cxa(1:3,i) = cxa(1:3,i) + 1.0_dp
!             
!             do j=1,3
!                if (ld(j) < 0.0_dp) then
!                   cxa(j,i)=ld(j)+1.0_dp
!                else
!                   cxa(j,i)=ld(j)
!                endif
!             enddo

! c            cxa(1,i)=(unshad(1,i)+tstra(1,i))/ax1
! c            if (cxa(1,i)<0.0) then
! c               cxa(1,i)=cxa(1,i)+1.0
! c            endif
! c            cxa(2,i)=(unshad(2,i)+tstra(2,i))/ay2
! c            if (cxa(2,i)<0.0) then
! c               cxa(2,i)=cxa(2,i)+1.0
! c            endif
! c            cxa(3,i)=(unshad(3,i)+tstra(3,i))/az3
! c            if (cxa(3,i)<0.0) then
! c               cxa(3,i)=cxa(3,i)+1.0
! c            endif

         enddo

         do i=1,ninert
            do j=1,3
               tstra(j) = sum( et(j,1:3) * unshind(1:3,i)) * totstr
            enddo

            adstra(1:3) = unshind(1:3,i) + tstra(1:3)

            call mul3x3(b,adstra,ld)
            
            cxi(1:3,i) = ld
            where (ld < 0.0_dp) cxi(1:3,i) = cxi(1:3,i) + 1.0_dp
            
!             do j=1,3
!                if (ld(j) < 0.0_dp) then
!                   cxi(j,i)=ld(j)+1.0_dp
!                else
!                   cxi(j,i)=ld(j)
!                endif
!             enddo
! ! 
! c            cxi(1,i)=(unshind(1,i)+tstrin(1,i))/ax1
! c            if (cxi(1,i)<0.0) then
! c               cxi(1,i)=cxi(1,i)+1.0
! c            endif
! c            cxi(2,i)=(unshind(2,i)+tstrin(2,i))/ay2
! c            if (cxi(2,i)<0.0) then
! c               cxi(2,i)=cxi(2,i)+1.0
! c            endif
! c            cxi(3,i)=(unshind(3,i)+tstrin(3,i))/az3
! c            if (cxi(3,i)<0.0) then
! c               cxi(3,i)=cxi(3,i)+1.0
! c            endif

         enddo

      else

         do i = 1, nd

            call mul3x3(b,ad(1,i),ld)

! c            write(*,*) ld(1),ld(2),ld(3)

            cxa(1:3,i) = ld
            where (ld < 0.0_dp) cxa(1:3,i) = cxa(1:3,i) + 1.0_dp
            
!             do j=1,3
!                if (ld(j) < 0.0_dp) then
!                   cxa(j,i)=ld(j)+1.0_dp
!                else
!                   cxa(j,i)=ld(j)
!                endif
!             enddo

! c            cxa(i)=ad(1,i)/ax
! c            if (cxa(1,i)<0.0) then
! c               cxa(1,i)=cxa(1,i)+1.0
! c            endif
! c            cxa(i)=ad(2,i)/ay
! c            if (cxa(2,i)<0.0) then
! c               cxa(2,i)=cxa(2,i)+1.0
! c            endif
! c            cxa(i)=ad(3,i)/az
! c            if (cxa(3,i)<0.0) then
! c               cxa(3,i)=cxa(3,i)+1.0
! c            endif

         enddo

         do i=1,ninert
            call mul3x3(b,adinert(1,i),ld)

            cxi(1:3,i) = ld
            where (ld < 0.0_dp) cxi(1:3,i) = cxi(1:3,i) + 1.0_dp
            
!             do j=1,3
!                if (ld(j) < 0.0_dp) then
!                   cxi(j,i)=ld(j)+1.0_dp
!                else
!                   cxi(j,i)=ld(j)
!                endif
!             enddo

! c            cxi(i)=adinert(1,i)/ax
! c            if (cxi(1,i)<0.0) then
! c               cxi(1,i)=cxi(1,i)+1.0
! c            endif
! c            cxi(i)=adinert(2,i)/ay
! c            if (cxi(2,i)<0.0) then
! c               cxi(2,i)=cxi(2,i)+1.0
! c            endif
! c            cxi(i)=adinert(3,i)/az
! c            if (cxi(3,i)<0.0) then
! c               cxi(3,i)=cxi(3,i)+1.0
! c            endif

         enddo

      endif


      write(num,'("Number of particles = ",i6)') nd+ninert
      write(num,'("A = 1.0 angstrom (basic length-scale)")')

      write(num,'("H0(1,1) = ",f8.3," a")') atmp(1,1)
      write(num,'("H0(1,2) = ",f8.3," a")') atmp(1,2)
      write(num,'("H0(1,3) = ",f8.3," a")') atmp(1,3)

      write(num,'("H0(2,1) = ",f8.3," a")') atmp(2,1)
      write(num,'("H0(2,2) = ",f8.3," a")') atmp(2,2)
      write(num,'("H0(2,3) = ",f8.3," a")') atmp(2,3)

      write(num,'("H0(3,1) = ",f8.3," a")') atmp(3,1)
      write(num,'("H0(3,2) = ",f8.3," a")') atmp(3,2)
      write(num,'("H0(3,3) = ",f8.3," a")') atmp(3,3)

      write(num,'(".NO_VELOCITY.")')
      write(num,'("entry_count = 9")')

      write(num,'("auxiliary[0] = force")')
      write(num,'("auxiliary[1] = force_x")')
      write(num,'("auxiliary[2] = force_y")')
      write(num,'("auxiliary[3] = force_z")')
      write(num,'("auxiliary[4] = magnetic moment")')
      write(num,'("auxiliary[5] = active/inert")')
! c      write(num,'("auxiliary[5] = erep")')
! c      write(num,'("auxiliary[6] = ebond")')
! c      write(num,'("auxiliary[7] = etot")')

! c
! c     write out active atoms
! c
      do i = 1, nd
         write(num,'(f12.5)') mass(i)
         write(num,'(a)') symb(i)
         ff = sqrt(sum(ftot(1:3,i)*ftot(1:3,i)))

         write(num,'(3f12.8,4f12.5,f10.6,i6)') cxa(1:3,i),ff,ftot(1:3,i),mg(i),i
      enddo

! c
! c     write out inert atoms
! c
      do i = 1, ninert
         write(num,'(f12.5)') 1.0_dp
! c         write(num,'("1.00000")')
! c         write(num,'("ar")')
         write(num,'(a)') symi(i)
         write(num,'(3f12.8," 0.0 0.0 0.0 0.0",f10.6," -1000")')  cxi(1:3,i),2*demi(i)/istn(zinert(i))
      enddo

      close(num)

! c      deallocate (tstra)
! c      deallocate (tstrin)
      deallocate (cxa)
      deallocate (cxi)

! c      print *,"cfg written"

      end subroutine writecfg



! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      subroutine writexbs(filename)
! c
! c     creates xbs file
! c
      use mod_const
      use mod_all_scalar, only : nd, ninert

      implicit none

      include 'Include/PosVel.array'
      
      
      integer i,num
      character(len=*) filename

! c
! c     create xbs output file
! c

      num=33
      open (unit=num,status='unknown',form='formatted',file=filename)

      write(33,'("* active atoms")')
      do i = 1, nd
         write(33,'("atom  Ac",4x,3(f14.8))')  ad(1:3,i)
      enddo

      write(33,'("* inert atoms")' )
      do i=1,ninert
         write(33,'("atom  In",4x,3(f14.8))')  adinert(1:3,i)
      enddo

      write(33,'("* atom data")' )
      write(33,'("spec  Ac",4x,"0.40",4x,"blue")')
      write(33,'("spec  In",4x,"0.40",4x,"white")')
      write(33,'("* bond data")' )
      write(33,'("bonds  Ac  Ac  0.000  3.000  0.1  1.0")')

      close(33)

      end subroutine writexbs


