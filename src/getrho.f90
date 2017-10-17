 
      subroutine getrho(elo,ehi)
          use mod_precision
          use mod_all_scalar
          use mod_const
          use ab_io
          use mod_funptr
!
!    This is a routine to evaluate the electron density at each active
!     atom, for a given range of energies.
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

      include "Include/Atom.array"
!      include "Include/NotSRT.array"
      include "Include/PosVel.array"
!      include "Include/RecurCoef.array"
!      include "Include/SRT.array"

!
!    Declare the simple variables.
!

      real(dp) :: elo,ehi,rho,nlo,nhi
      real(dp) :: numelnull

      integer la
      integer ia
      integer wrflag

      character*80 filename

!
!    Open the output file.
!

      if (idebug == 1) write(6,'(''LENGFN = '',I5)') lengfn

      filename = genfile(1:lengfn)//'.rho'
      open(unit = 99,file = filename,status = 'NEW')

!
!    Evaluate the densities.
!

      wrflag = 1
      do ia = 1,nd,1
         call rdab(wrflag,ia)
         rho = 0.0_dp
         do la = 1,nchain,1
            if (term == 1) then
               nlo = numelsrt(elo,la)
               nhi = numelsrt(ehi,la)
            elseif ((term == 2).or.(term == 3)) then
               nlo = numelnull(elo,lchain(la),mrec,diag(1,la), & 
     &                         eigvec(1,1,la),kt)
               nhi = numelnull(ehi,lchain(la),mrec,diag(1,la), & 
     &                         eigvec(1,1,la),kt)
            endif
            rho = rho + wt(la)*(nhi-nlo)
         enddo
         write(99,'(''['',F9.3,'','',F9.3,'','',F9.3,''] '',G22.15)') & 
     &         ad(1,ia),ad(2,ia),ad(3,ia),rho
      enddo

!
!    Close the output file.
!

      close(50)
      close(99)

      end

