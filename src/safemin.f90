 
   subroutine safemin(oldad,dir,ein,eout,status)
      use mod_precision
      use mod_all_scalar
      use mod_const
          
!
!    This is a subroutine to find the minimum of a function of many variables.
!     This routine only carries out one line minimization. It is designed to
!     pick up errors in forces. Stability is taken as more important than speed
!     here.
!
!    OLDAD(3,ND) : The atomic coordinates on entry
!    DIR(3,ND)   : The search direction
!    EIN         : Energy on entry
!    EOUT        : Energy on exit
!    STATUS      : Error flag. 0 => Successful minimization
!                              1 => Error in forces
!

      implicit none

      
      real(dp), intent(in) :: ein, oldad(3,mxnd), dir(3,mxnd)
      real(dp), intent(out) :: eout
      integer, intent(out) :: status
      
      include "Include/Force.array"
      include "Include/PosVel.array"

      real(dp) :: k0,e0,kmin
      real(dp) :: k1,k2,k3,lambda
      real(dp) :: e1,e2,e3,lef0
      real(dp) :: maxd, maxdis(3)

      integer :: ia,flag, maxidx(1)
      
      character(len=1), parameter :: cdrname(3) = ['x','y','z']

      
      
!       the ef related things may be useful but it seems the time spent 
!   in finding ef is rather small so these are commented out in most of the code
!   flag=2 is pretty much meaningless now
      
!
!    Test forces.
!

      status = 0
!       lef0 = lef

      k0 = 1.0e-6_dp
!       do ia = 1,nd
!          ad(1,ia) = oldad(1,ia) + k0*dir(1,ia)
!          ad(2,ia) = oldad(2,ia) + k0*dir(2,ia)
!          ad(3,ia) = oldad(3,ia) + k0*dir(3,ia)
!          lef = lef + k0*dir(1,ia)*dmdl(1,ia) & 
!      &             + k0*dir(2,ia)*dmdl(2,ia) & 
!      &             + k0*dir(3,ia)*dmdl(3,ia)
!       end do
      
      ad(:,:nd) = oldad(:,:nd) + k0*dir(:,:nd)
      
      flag = 2
      call getetot(flag)
      if (.not. quiet) write(6,'(''LOCC = '',G12.5,'' TOTNIA = '',G12.5)') locc,totnia
      if (abs(locc-totnia)/real(nd, dp) > 0.001_dp) then
         flag = 1
         call getetot(flag)
      endif
      e0 = eprom + ebond + epair - eatom - uent
      if (.not. quiet) then
          write(9,'("flag=",i1,"  e0=",g22.12)') flag,e0
          write(6,'("RELAX flag=",i1,"  e0=",g22.12)') flag,e0
      endif
      iter=iter+1
      
!       transformed from MM's
      if (e0 > ein) then
         if (.not. quiet) write(6,'(/''The forces are inconsistent with the energy.'')')
!   CHANEGE THIS BACK TO 1.0_DP
         if ( (e0-ein) < 1.0_dp ) then
            if (.not. quiet) then
               write(6,'(''but i am going to be nice about it!'')')
               write(6,'(''jumping out of safemin.'')')
               write(6,'(/''* check if the relaxation converges !!! *'')')
            end if
            kmin = 0.001_dp
            ad(:,:nd) = oldad(:,:nd) + kmin*dir(:,:nd)
         else
            if (.not. quiet) write(6,'(''They are making the energy INCREASE.''/)')
            status = 1
         endif
         return
      endif

!
!    Bracket minimum.
!

      k1 = 0.0_dp
      e1 = ein
      k2 = k0
      e2 = e0
      lambda = 2.0_dp

 2    if (k2 == k0) then
         k3 = 0.001_dp
      elseif (k2 == 0.01_dp) then
         k3 = 0.01_dp
      else
         k3 = lambda*k2          
      endif

      
!     until the endif is transformed from Matous'      
      maxdis = k3*maxval(dir(:,:nd), dim=2)
      
      if (any(maxdis > 0.1_dp)) then
         if (.not. quiet) write(*,'(" *** WARNING ***")')
         if (.not. quiet) write(*,'(" Displacement > 0.1 ")')

         maxidx = maxloc(maxdis)
!          maxd = maxdis(maxidx(1))
         maxd   = maxval(maxdis)
         if (.not. quiet) write(*,'("maxd_",a1,"=",f12.3)') cdrname(maxidx(1)), maxd
         
         k3 = k3/(int(maxd/0.1_dp)+1.0_dp) ! not quite sure what the int in doing
         if (.not. quiet) write(*,'("Setting k3 to ",e10.3)') k3
         
         ad(:,:nd) = oldad(:,:nd) + k3*dir(:,:nd)

         flag = 1
         call getetot(flag)
! c         write(6,'(''locc = '',g12.5,'' totnia = '',g12.5)') locc,totnia
         eout = eprom + ebond + epair + eenv - eatom - uent
         if (.not. quiet) then
            write(9,'("flag=",i1,"  eout=",g22.12)') flag,eout
            write(6,'("RELAX flag=",i1,"  eout=",g22.12)') flag,eout
          endif
         iter=iter+1

         return
     end if
      
!       lef = lef0
!       do ia = 1,nd
!          ad(1,ia) = oldad(1,ia) + k3*dir(1,ia)
!          ad(2,ia) = oldad(2,ia) + k3*dir(2,ia)
!          ad(3,ia) = oldad(3,ia) + k3*dir(3,ia)
!          lef = lef + k3*dir(1,ia)*dmdl(1,ia) & 
!      &             + k3*dir(2,ia)*dmdl(2,ia) & 
!      &             + k3*dir(3,ia)*dmdl(3,ia)
!       end do

      ad(:,:nd) = oldad(:,:nd) + k3*dir(:,:nd)
      
      
      flag = 2
      call getetot(flag)
      if (.not. quiet) write(6,'(''LOCC = '',G12.5,'' TOTNIA = '',G12.5)') locc,totnia
      if (abs(locc-totnia)/real(nd, dp) > 0.001_dp) then
         flag = 1
         call getetot(flag)
      endif
      e3 = eprom + ebond + epair - eatom - uent
      
      if (.not. quiet) then
          write(9,'("flag=",i1,"  e3=",g22.12)') flag,e3
          write(6,'("RELAX flag=",i1,"  e3=",g22.12)') flag,e3
      endif
      iter=iter+1
      
      if (e3 < e2) then
         k1 = k2
         e1 = e2
         k2 = k3
         e2 = e3
         goto 2
      endif

!
!    Interpolate to find minimum.
!

      kmin = 0.5_dp*(((k1*k1-k3*k3)*(e1-e2)-(k1*k1-k2*k2)*(e1-e3)) & 
     &              /((k1-k3)*      (e1-e2)-(k1-k2)*      (e1-e3)))

!       lef = lef0
!       do ia = 1,nd,1
!          ad(1,ia) = oldad(1,ia) + kmin*dir(1,ia)
!          ad(2,ia) = oldad(2,ia) + kmin*dir(2,ia)
!          ad(3,ia) = oldad(3,ia) + kmin*dir(3,ia)
!          lef = lef + kmin*dir(1,ia)*dmdl(1,ia) & 
!      &             + kmin*dir(2,ia)*dmdl(2,ia) & 
!      &             + kmin*dir(3,ia)*dmdl(3,ia)
!       end do
      
      ad(:,:nd) = oldad(:,:nd) + kmin*dir(:,:nd)
      
      flag = 2
      call getetot(flag)
      if (.not. quiet) write(6,'(''LOCC = '',G12.5,'' TOTNIA = '',G12.5)') locc,totnia
      if (abs(locc-totnia)/real(nd, dp) > 0.001_dp) then
         flag = 1
         call getetot(flag)
      endif
      eout = eprom + ebond + epair - eatom - uent
      if (.not. quiet) THEN
          write(9,'("flag=",i1,"  eout=",g22.12)') flag,eout
          write(6,'("RELAX flag=",i1,"  eout=",g22.12)') flag,eout
      ENDIF
      iter=iter+1
      
   end subroutine safemin

