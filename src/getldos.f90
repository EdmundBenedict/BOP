 
      subroutine getldos(ia,la,ne,flag)
          use mod_precision
          use mod_all_scalar
          use mod_const
!           use mod_ham
          use ab_io
          use mod_funptr
          use mod_g0n, only : g00
!
!    This is a routine to calculate the density of states
!     and the integerated density of states.
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
!      include "Include/RecurCoef.array"
!      include "Include/SRT.array"

!
!    Declare the simple variables.
!


!       complex(dp) :: g00

      real(dp) :: e,estep,idos,dos,numelnull
      real(dp) :: emin,emax

      integer ne,flag
      integer i,wrflag
      integer ia,la

      character*80 filename
      character*1 lnum
      character*4 anum

!
!    Read in data.
!

      write(6,'(''LA = '',I5,'' WT = '',F12.5)') la,wt(la)

!       wrflag = 1
!       do i = 1,ia,1
!          call rdab(wrflag,i)
!       enddo
!       close(50)


        call assoc_ab(ia,1)

!
!    Find the limits of the density of states.
!

      if (flag == 1) then
         if (term == 1) then
            emin = (lainf(la)-2.0_dp*lbinf(la))
            emax = (lainf(la)+2.0_dp*lbinf(la))
         elseif ((term == 2).or.(term == 3)) then
            emin = diag(1,la) - 10.0_dp*kt
            emax = diag(lchain(la)+1,la) + 10.0_dp*kt
         endif
      elseif (flag == 2) then
         emin = 1.0e30_dp
         emax = -emin
         if (term == 1) then
            do la = 1,nchain,1
               emin = min((lainf(la)-2.0_dp*lbinf(la)),emin)
               emax = max((lainf(la)+2.0_dp*lbinf(la)),emax)
            enddo
         elseif ((term == 2).or.(term == 3)) then
            do la = 1,nchain,1
               emin = min(diag(1,la),emin)
               emax = max(diag(lchain(la)+1,la),emax)
            enddo
            emin = emin - etail
            emax = emax + etail
         endif
      endif

!
!    Open the output files.
!

      if (flag == 1) then

         write(lnum,'(I1)') la

         write(anum,'(I4)') ia
         if (ia < 10) then
            anum(1:3) = '000'
         elseif (ia < 100) then
            anum(1:2) = '00'
         elseif (ia < 1000) then
            anum(1:1) = '0'
         endif

         filename = genfile(1:lengfn)//'_A'//anum//'_L'//lnum//'.ldos'
         open(unit = 10,file = filename,status = 'NEW')

         filename = genfile(1:lengfn)//'_A'//anum//'_L'//lnum//'.ildos'
         open(unit = 11,file = filename,status = 'NEW')

      elseif (flag == 2) then

         write(anum,'(I4)') ia
         if (ia < 10) then
            anum(1:3) = '000'
         elseif (ia < 100) then
            anum(1:2) = '00'
         elseif (ia < 1000) then
            anum(1:1) = '0'
         endif

         filename = genfile(1:lengfn)//'_A'//anum//'.ldos'
         open(unit = 10,file = filename,status = 'NEW')

         filename = genfile(1:lengfn)//'_A'//anum//'.ildos'
         open(unit = 11,file = filename,status = 'NEW')

      endif

!
!    Write out the density of states and the integrated density of states.
!

      if (flag == 1) then

         estep = (emax-emin)/real(ne-1, dp)
         e = emin
         do i = 1,ne,1
            if (term == 1) then
               dos = -wt(la)*aimag(g00(cmplx(e, kind=dp), & 
     &               arec(0:lchain(la),la),brec(0:lchain(la)+1,la), & 
     &               lchain(la),lainf(la),lbinf(la)))/pi
               idos = numelsrt(e,la)*wt(la)
            elseif ((term == 2).or.(term == 3)) then
               dos = -wt(la)*aimag(g00(cmplx(e,kt, kind=dp), & 
     &               arec(0:lchain(la),la),brec(0:lchain(la)+1,la), & 
     &               lchain(la),lainf(la),lbinf(la)))/pi
               idos = numelnull(e,lchain(la),mrec,diag(1,la), & 
     &                          eigvec(1,1,la),kt)*wt(la)
            endif
            write(10,*) e,dos/wt(la)
            write(11,*) e,idos/wt(la)
            e = e+ estep
         enddo

      elseif (flag == 2) then

         estep = (emax-emin)/real(ne-1, dp)
         e = emin
         do i = 1,ne,1
            if (term == 1) then
               dos = 0.0_dp
               idos = 0.0_dp
               do la = 1,nchain,1
                  dos = dos-wt(la)*aimag(g00(cmplx(e, kind=dp), & 
     &                  arec(0:lchain(la),la),brec(0:lchain(la)+1,la), & 
     &                  lchain(la),lainf(la),lbinf(la)))/pi
                  idos = idos + & 
     &                   numelsrt(e,la)*wt(la)
               enddo
            elseif ((term == 2).or.(term == 3)) then
               dos = 0.0_dp
               idos = 0.0_dp
               do la = 1,nchain,1
                  dos = dos-wt(la)*aimag(g00(cmplx(e,kt, kind=dp), & 
     &                  arec(0:lchain(la),la),brec(0:lchain(la)+1,la), & 
     &                  lchain(la),lainf(la),lbinf(la)))/pi
                  idos = idos + & 
     &                   numelnull(e,lchain(la),mrec,diag(1,la), & 
     &                   eigvec(1,1,la),kt)*wt(la)
               enddo
            endif
            write(10,*) e,dos/wt(la)
            write(11,*) e,idos/wt(la)
            e = e+ estep
         enddo

      endif

!
!    Close the output files.
!

      close(10)
      close(11)

      end

