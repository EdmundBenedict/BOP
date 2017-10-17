 
      subroutine lgfsafemin(oldad,dir,ein,eout,status,testf,morig, mfcomp, mcad)
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

!      include "Include/Atom.array"
      include "Include/Force.array"
!      include "Include/Hamilt.array"
!      include "Include/Misc.array"
!      include "Include/NebList.array"
      include "Include/PosVel.array"
!      include "Include/Relax.array"

!
!    Declare the simple variables.
!

      real(dp) :: ein,eout,k0,e0,kmin
      real(dp) :: k1,k2,k3,lambda
      real(dp) :: e1,e2,e3,lef0
      real(dp) :: testf, tmpdis, meadj, mcad(3,nd)

      integer status,ia,flag
      integer morig, mfcomp

!
!    Declare the arrays.
!

      real(dp) :: oldad(3,mxnd)
      real(dp) :: dir(3,mxnd)

!
!    Test forces.
!

      status = 0
      lef0 = lef

      k0 = 1.0e-6_dp
      do 1 ia = 1,nd,1
         ad(1,ia) = oldad(1,ia) + k0*dir(1,ia)
         ad(2,ia) = oldad(2,ia) + k0*dir(2,ia)
         ad(3,ia) = oldad(3,ia) + k0*dir(3,ia)
         lef = lef + k0*dir(1,ia)*dmdl(1,ia) & 
     &             + k0*dir(2,ia)*dmdl(2,ia) & 
     &             + k0*dir(3,ia)*dmdl(3,ia)
 1    continue
      flag = 2
      call getetot(flag)
      write(6,'(''LOCC = '',G12.5,'' TOTNIA = '',G12.5)') locc,totnia
      if (abs(locc-totnia)/real(nd, dp) > 0.001_dp) then
         flag = 1
         call getetot(flag)
      endif
!
!     Adjusting energy due to work done on atom at origin
!
      tmpdis=ad(mfcomp,morig)-mcad(mfcomp,morig)
      meadj = -1.0_dp*testf*tmpdis
      e0 = eprom + ebond + epair - eatom - uent + meadj
      ftot(mfcomp,morig) = ftot(mfcomp,morig) + testf
!
!      IF (E0.GT.EIN) THEN
!         WRITE(6,'(/''The forces are inconsistent with the energy.'')')
!         WRITE(6,'(''They are making the energy INCREASE.''/)')
!         STATUS = 1
!         RETURN
!      ENDIF

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

      lef = lef0
      do 3 ia = 1,nd,1
         ad(1,ia) = oldad(1,ia) + k3*dir(1,ia)
         ad(2,ia) = oldad(2,ia) + k3*dir(2,ia)
         ad(3,ia) = oldad(3,ia) + k3*dir(3,ia)
         lef = lef + k3*dir(1,ia)*dmdl(1,ia) & 
     &             + k3*dir(2,ia)*dmdl(2,ia) & 
     &             + k3*dir(3,ia)*dmdl(3,ia)
 3    continue
      flag = 2
      call getetot(flag)
      write(6,'(''LOCC = '',G12.5,'' TOTNIA = '',G12.5)') locc,totnia
      if (abs(locc-totnia)/real(nd, dp) > 0.001_dp) then
         flag = 1
         call getetot(flag)
      endif
      tmpdis=ad(mfcomp,morig)-mcad(mfcomp,morig)
      meadj = -1.0_dp*testf*tmpdis
      e3 = eprom + ebond + epair - eatom - uent + meadj
      ftot(mfcomp,morig) = ftot(mfcomp,morig) + testf      

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

      kmin = 0.5_dp*(((k1*k1-k3*k3)*(e1-e2)-(k1*k1-k2*k2)*(e1-e3))/ & 
     &              ((k1-k3)*      (e1-e2)-(k1-k2)*      (e1-e3)))

      lef = lef0
      do 4 ia = 1,nd,1
         ad(1,ia) = oldad(1,ia) + kmin*dir(1,ia)
         ad(2,ia) = oldad(2,ia) + kmin*dir(2,ia)
         ad(3,ia) = oldad(3,ia) + kmin*dir(3,ia)
         lef = lef + kmin*dir(1,ia)*dmdl(1,ia) & 
     &             + kmin*dir(2,ia)*dmdl(2,ia) & 
     &             + kmin*dir(3,ia)*dmdl(3,ia)
 4    continue
      flag = 2
      call getetot(flag)
      write(6,'(''LOCC = '',G12.5,'' TOTNIA = '',G12.5)') locc,totnia
      if (abs(locc-totnia)/real(nd, dp) > 0.001_dp) then
         flag = 1
         call getetot(flag)
      endif
      tmpdis=ad(mfcomp,morig)-mcad(mfcomp,morig)
      meadj = -1.0_dp*testf*tmpdis
      eout = eprom + ebond + epair - eatom - uent + meadj
      ftot(mfcomp,morig) = ftot(mfcomp,morig) + testf

      end

