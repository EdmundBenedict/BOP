 
      subroutine getvib()
          use mod_precision
          use mod_all_scalar
          use mod_const
          use mod_atom_ar
          
          
!** This, together with bldlcell has to be reviewed and either removed or adapted to the parallel distribution and dynamic arrays
          
!
!    This is a routine to evaluate the total energy and forces for a configuration.
!
!     FLAG = 0 => Read the coefficients from disk, and find Fermi energy.
!     FLAG = 1 => Evaluate the coefficients, write them to disk, and find Fermi energy.
!     FLAG = 2 => Evaluate the coefficients, do not write them to disk, and do not find Fermi energy.
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

      include "Include/Force.array"
      include "Include/Atom.array"
      include "Include/PosVel.array"
      include "Include/Vib.array"

!
!    Declare the simple variables.
!


      integer i,j,ii,jj,ia,flag,nds

      real(dp) :: de_o(mxnd),dq_o(mxnd)
      real(dp) :: dmu,dmdl_0(3,mxnd),dndm_0,lef_0,totnia_0

      character*80 filename
!
!    Set DELTA by hand in this version
!
!     DELTA=0.01D0
!
!    Create the Large Unit cell and Mappings
!
      nds=nd
      call bldlcell()
!CC
      call getocc(1,locc)
      flag = 1
      call getetot(flag)
      call erasab()
!CC
      do i=1,nd
       dmdl_0(1,i)=dmdl(1,i)
       dmdl_0(2,i)=dmdl(3,i)
       dmdl_0(3,i)=dmdl(2,i)
      enddo
      dndm_0=dndm
      lef_0=lef
      totnia_0=totnia
!CC
      do i=1,nd
       dq_o(i)=dq(i)
       de_o(i)=de(i)
      enddo
!
!
!     Loops here to calculate the Force constant matrix
!
      do ii=1,3
       do jj=1,nds
        write(6,'(''II='',I5,'' JJ='',I5)')ii,jj
!       CALL FLUSH(6)
        do ia=1,2
         if(ia == 1)ad(ii,jj)=ad(ii,jj)+delta
         if(ia == 2)ad(ii,jj)=ad(ii,jj)-delta
!
!    Evaluate the band structure energy and forces.
!
!CC
         if (dndm_0 > 1.0e-6_dp) then
            dmu = (locc-totnia_0)/dndm_0
            if(ia == 1)dmu = dmu + dmdl_0(ii,jj)*delta
            if(ia == 2)dmu = dmu - dmdl_0(ii,jj)*delta
            lef = lef_0 + dmu
         endif
!CC
!
         flag=2
         call getetot(flag)
!
!    Evaluate the total forces.
!
         write(*,'("zacatek")')
         do i=1,3
          do j = 1,nd,1
           if(ia == 1) dynr(ii,i,jj,j)=-0.5*(fbs(i,j)+fpp(i,j))/delta
           if(ia == 2) dynr(ii,i,jj,j)=dynr(ii,i,jj,j)+0.5*(fbs(i,j)+fpp(i,j))/delta
         write(*,'(''dynr(1,'',i3,'',1,'',i3,'')='',f15.6)') & 
     &        i,j,dynr(ii,i,jj,j)
         write(*,'("hola")')
          enddo
         enddo
!
         do i=1,nd
          dq(i)=dq_o(i)
          de(i)=de_o(i)
         enddo
         if(ia == 1)ad(ii,jj)=ad(ii,jj)-delta
         if(ia == 2)ad(ii,jj)=ad(ii,jj)+delta
!
        enddo
       enddo
      enddo

      filename = genfile(1:lengfn)//'.afcm'
      open(unit=65,file=filename,form='unformatted',status='unknown')
      write(65) dynr,mapr,mapx
      close(65)

      end

