 
       function getefnull(nelec)
          use mod_precision


          use mod_all_scalar

          use mod_const

!
!    This is a function to find the fermi energy for the null terminator.
!

      implicit none
      real(dp) :: getefnull


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
!      include "Include/NotSRT.array"
!      include "Include/RecurCoef.array"


!
!    Common blocks introduced by A.Girshick.
!
 
      include "Include/ag.conn"   
     

!
!    Declare the simple variables.
!

      real(dp) :: nelec,efhi,nf,nflo,nfhi
      real(dp) :: ef,numel,def

      integer niter


!
!    Find limits within which fermi energy lies.
!

!     WRITE(6,'(''NELEC = '',G12.5)') NELEC

      nflo = numel(eflo)

!     WRITE(6,'(''EFLO = '',G12.5,'' NFLO = '',G12.5)') EFLO,NFLO

      def = 1.0_dp
      if (nflo > nelec) then
         do while(nflo > nelec)
            efhi = eflo
            nfhi = nflo
            eflo = eflo - def
            nflo = numel(eflo)
!           WRITE(6,'(''EFLO = '',G12.5,'' NFLO = '',G12.5)') EFLO,NFLO
         enddo
      else
         efhi = eflo
         nfhi = nflo
         do while(nfhi < nelec)
            eflo = efhi
            nflo = nfhi
            efhi = efhi + def
            nfhi = numel(efhi)
!           WRITE(6,'(''EFHI = '',G12.5,'' NFHI = '',G12.5)') EFHI,NFHI
         enddo
      endif

!
!    Refine estimate of the Fermi energy.
!

      nf = -1.0_dp
      niter = 1

      do while ((abs(nf-nelec) > 1.0e-8_dp*real(nd, dp)).and. & 
     &          (niter <= 100))

!
!       Bisection step.
!

         ef = 0.5_dp*(efhi+eflo)
         nf = numel(ef)

         if (nf > nelec) then
            efhi = ef
         else
            eflo = ef
         endif

!
!       Linear interpolation step.
!

         ef = (efhi*(nelec-nflo)+eflo*(nfhi-nelec))/(nfhi-nflo)
         if ((ef < eflo).or.(ef > efhi)) ef = 0.5_dp*(efhi+eflo)
         nf = numel(ef)

         if (nf > nelec) then
            efhi = ef
         else
            eflo = ef
         endif

!        WRITE(6,'(''EFLO = '',G12.5,'' EFHI = '',G12.5,
!    +             '' NF = '',G12.5)') EFLO,EFHI,NF

         niter = niter + 1

      enddo

      getefnull = ef
      totnia = nf

!     WRITE(6,'(''NITER = '',I4)') NITER

      end

