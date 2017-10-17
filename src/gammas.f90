 
      subroutine gammas()
          use mod_precision


          use mod_all_scalar

          use mod_const

          use mod_srt

!
!     Calculates the gamma-surface. AMB 10.1.95
!
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
      include "Include/BEC.array"
!      include "Include/BondOrder.array"
      include "Include/Force.array"
!      include "Include/Hamilt.array"
!      include "Include/Misc.array"
!      include "Include/Moment.array"
!      include "Include/NebList.array"
!      include "Include/NotSRT.array"
      include "Include/PosVel.array"
!      include "Include/RecurCoef.array"
      include "Include/Relax.array"
!      include "Include/SRT.array"

!
!    Declare the simple variables.
!

      real(dp) :: ein,eout,ered,u0
      real(dp) :: alpha,maxf,magf,dmu
      real(dp) :: gd1,alpha1,alpha2,alphau,alphal
      real(dp) :: xstart,ystart

      integer ng1,ng2
      integer usrexit
      integer i,j,flag,iat
      integer icom,nsrch
      common /gamcom/ xstart,ystart,ng1,ng2

      character*80 filename
      integer ix,iy
      integer ngmax
      parameter (ngmax=21)
      real(dp) :: gs(ngmax,ngmax),acell(3,3)
      real(dp) :: step0,shift(2),offset
      real(dp) :: veca(3),vecb(3),vecc(3),area,len1,len2
      real(dp) :: adw(3,mxtotnd)
      real(dp) :: adinw(3,minert)
      real(dp) :: len1a,len1b,len2a,len2b,zero

!
!     SAVE AN INITIAL STEP
!
      zero = 0.0_dp

      step0=step

!     print*,'<gammas>:xstart,ystart=',xstart,ystart

      if(ng1 > ngmax)stop'<GAMMA>: NG1 v NGMAX MISMATCH'
      if(ng2 > ngmax)stop'<GAMMA>: NG2 v NGMAX MISMATCH'

      if(mxiter == 0)then
         write(6,*)' UNRELAXED GAMMA-SURFACE '
      else
         write(6,*)'   RELAXED GAMMA-SURFACE '
      endif
!
!     SAVE A COPY OF UNDEFORMED UNIT CELL
!
      do i=1,3
         do j=1,3
            acell(i,j)=a(i,j)
         enddo
      enddo

      do iat=1,nd
         do i=1,3
         adw(i,iat)=ad(i,iat)
         enddo
      enddo
      do iat = 1, ninert
         do i=1,3
         adinw(i,iat)=adinert(i,iat)
         enddo
      enddo


!     MAKE A NORMALISATION FOR THE SURFACE OF A FAULT

      do i=1,3
         veca(i)=a(1,i)
         vecb(i)=a(2,i)
      enddo

      call cross(veca,vecb,vecc)
      area=sqrt( vecc(1)**2 + vecc(2)**2 + vecc(3)**2 )
!    EV/A^2 --> mJ/m^2
      area=area/16021.3_dp

      len1=0e0_dp
      len2=0e0_dp
      len1=sqrt(veca(1)**2+veca(2)**2+veca(3)**2)
      len2=sqrt(vecb(1)**2+vecb(2)**2+vecb(3)**2)

      filename = genfile(1:lengfn)//'.gs'
      open(unit = 70,file = filename,status = 'NEW')
      open(unit = 71,file = 'gs.gs',status='unknown')
      open(unit = 60,file = 'convrg.gs',status='unknown')

      write(6,'(2(1X,I3))')ng1,ng2
      len1a=zero
      len1b=len1
      len2a=zero
      len2b=len1
      if(ng1 <= 1)then
         len1a=xstart
         len1b=xstart
      endif
      if(ng2 <= 1)then
         len2a=ystart
         len2b=ystart
      endif
      write(70,'(4(1X,G10.3))')len1a,len1b,len2a,len2b
      write(71,'(4(1X,G10.3))')len1a,len1b,len2a,len2b
      write(70,'(2(1X,I3))')ng1,ng2
      write(71,'(2(1X,I3))')ng1,ng2

!    LOOP OVER THE GAMMA-SURFACE (P1,P2) <<-------------------------------
!
      offset=1e30_dp
      do 100 ix=1,ng1,1
         do 100 iy=1,ng2,1
!D       DO 100 IX=2,2
!D       DO 100 IY=2,2

            step=step0

!
!     Retrieve the unit cell
!

            do i = 1,3,1
               do j=1,2,1
                  a(j,i) = acell(j,i)
               enddo
            enddo

      do iat=1,nd
         do i=1,3
         ad(i,iat)=adw(i,iat)
         enddo
      enddo
      do iat = 1, ninert
         do i=1,3
         adinert(i,iat)=adinw(i,iat)
         enddo
      enddo

!
!     SHEAR THE UNIT CELL
!
      if(ng1 > 1)then
      shift(1)=xstart+((ix-1)*acell(1,1)+ & 
     &(iy-1)*acell(2,1))/(ng1-1.)
      else
         shift(1)=xstart
      endif
      if(ng2 > 1)then
         shift(2)=ystart+((ix-1)*acell(1,2)+ & 
     &(iy-1)*acell(2,2))/(ng2-1.)
      else
         shift(2)=ystart
      endif 
 
            do i = 1,2,1
               a(3,i) = acell(3,i) + shift(i)
            enddo

        do i=1,3,1            
        write(6,'(''A = '',3F13.5)') (a(i,j),j=1,3)
        write(60,'(''A = '',3F13.5)') (a(i,j),j=1,3)
        enddo
        print*
        
!
!    Initialize the Hessian matrix.
!

      if (hessmx < 3*nd) then
         write(6,'(''Increase the size of HESSMX to at least '', & 
     &             I5)') 3*nd
         call panic()
      endif

      do i = 1,3*nd,1
         do j = 1,3*nd,1
            hess(j,i) = 0.0_dp
         enddo
         hess(i,i) = 1.0_dp
      enddo

!
!    Relax structure.
!
      iter = 1
      icom = 0

      write(60,'(/''======================================'')')
      write(60,'(''Starting relaxation with:'')')
      write(60,'(''FTOL   = '',G12.5)') ftol
      write(60,'(''MXITER = '',I5)') mxiter
      write(60,'(''STEP,STPMIN = '',2G12.5)')step, stpmin

      flag = 1
      call getetot(flag)
      call erasab()

      ein  = eprom + ebond + epair - eatom - uent
      eout = ein
      u0   = eout + 0.5_dp*uent 
      ered = 1.0e30_dp

      do i = 1, nd, 1
         write(60,'(''IAT='',I4,1X,''AD='',3(1X,F10.5))') & 
     &        i,ad(1,i),ad(2,i),ad(3,i)
      enddo
      write(60,*)
      do i = 1, nd, 1
         write(60,'(''IAT='',I4,1X,''FORCE='',3(1X,F10.5))') & 
     &i,ftot(1,i),ftot(2,i),ftot(3,i)
      enddo
      call constr(ftot,nd,cnst_a,cnst_v,cnst_n,mcnst_n)

      maxf = 0.0_dp
      do i = 1,nd,1
         newftot(1,i) = -ftot(1,i)
         newftot(2,i) = -ftot(2,i)
         newftot(3,i) = -ftot(3,i)
         magf = sqrt(ftot(1,i)**2+ftot(2,i)**2+ftot(3,i)**2)

         if (magf > maxf) maxf = magf
      enddo

      do while ((iter <= mxiter).and.(maxf > ftol).and. & 
     &          (icom /= 1).and.(usrexit() == 0).and. & 
     &          (step > stpmin))

         write(6,'(/''Iteration #             '',I5)') iter
         write(6,'(''Free energy,U(T=0)(start)= '',2G22.15)')ein,u0
!        WRITE(6,'(''Fermi energy (start)     = '',G12.5)') LEF
         write(6,'(''Maximum force (start)    = '',G12.5)') maxf

         call mstmin(3*nd,3*mxnd,hess,newftot,ftot,cg,newad, & 
     &               step,0.0_dp,alpha,icom, & 
     &               gd1,alpha1,alpha2,alphau,alphal,nsrch,6)

         dmu = (locc-totnia)/dndm
         call move(ad,nd,cg,alpha,dmdl,dmu)
         lef = lef + dmu

         flag = 2
         call getetot(flag)

         eout = eprom + ebond + epair - eatom - uent
         u0   = eout + 0.5_dp*uent

         call constr(ftot,nd,cnst_a,cnst_v,cnst_n,mcnst_n)

         maxf = 0.0_dp
         do i = 1,nd,1
            ftot(1,i) = -ftot(1,i)
            ftot(2,i) = -ftot(2,i)
            ftot(3,i) = -ftot(3,i)
            magf = sqrt(ftot(1,i)**2+ftot(2,i)**2+ftot(3,i)**2)
            if (magf > maxf) maxf = magf
         enddo

         write(6,'(''Free energy, U(T=0)(finish) = '',2G22.15)')eout,u0
!        WRITE(6,'(''Fermi energy (finish)     = '',G12.5)') LEF
         write(6,'(''Maximum force (finish)    = '',G12.5)') maxf

         ered = ein-eout
         if (ered < 0.0_dp) then
            write(6,'(/''Energy increasing => Reducing step.'')')
            write(6,'(''Old step = '',G12.5)') step
            step = step/2.0_dp
            write(6,'(''New step = '',G12.5)') step
         endif

         ein = eout
         iter = iter + 1

      enddo

      do i = 1,nd,1
         ftot(1,i) = -ftot(1,i)
         ftot(2,i) = -ftot(2,i)
         ftot(3,i) = -ftot(3,i)
      enddo

      write(60,'(''Stopping relaxation after ITER: '',I4)')iter
      if ((icom == 1).or.(maxf < ftol)) then
         write(60,'(''Relaxation successfully completed.'')')
      elseif (iter > mxiter) then
         write(60,'(''Maximum number of iterations exceeded.'')')
      elseif (step < stpmin) then
         write(60,'(''Minimum step size reached.'')')
      endif
      write(60,'(''Maximum force = '',G12.5,''eV/A'')') maxf

      do i = 1, nd, 1
         write(60,'(''IAT='',I4,1X,''AD='',3(1X,F10.5))') & 
     &        i,ad(1,i),ad(2,i),ad(3,i)
      enddo
      write(60,*)
      do i = 1, nd, 1
         write(60,'(''IAT='',I4,1X,''F='',3(1X,F10.5))') & 
     &        i,ftot(1,i),ftot(2,i),ftot(3,i)
      enddo
      write(60,'(''======================================''/)')


!
!     PRINT-OUT THE GAMMA-SURFACES
!
      gs(ix,iy) = u0
      offset = min(offset,u0)

! --> DUMP THE G-SURFACE (NO OFFSET !!) TO CHANNEL 71:

      write(71,'(4(1x,g13.7))')shift(1),shift(2),gs(ix,iy)
      
 100  continue

! -->  PRINT OUT THE G-SURFACE (WITH OFFSET, IN mJ/m^2  !!!)

      do iy=1,ng2,1
         do ix=1,ng1,1
            if(ng1 > 1.or.ng2 > 1)then
               gs(ix,iy)=(gs(ix,iy)-offset)/area
            else
               gs(ix,iy)=gs(ix,iy)/area
            endif
            write(70,'(4(1x,g13.7))')gs(ix,iy)
         enddo
      enddo
      close(70)
      close(71)
      close(60)

      end

      subroutine cross(a,b,c)
          use mod_precision


!
!    Evaluates cross product of two vectors. C = A x B
!

      real(dp) :: a(3),b(3),c(3)

      c(1) =  a(2)*b(3)-a(3)*b(2)
      c(2) = -a(1)*b(3)+a(3)*b(1)
      c(3) =  a(1)*b(2)-a(2)*b(1)
      end

