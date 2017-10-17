 
      subroutine getdos(ne)
          use mod_precision
          use mod_all_scalar
          use mod_const
          use ab_io
!           use mod_ham
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

      real(dp) :: e,delte,idos,dos,numelnull
      real(dp) :: emin,emax,maxdos

      integer i,wrflag,ne
      integer ia,la

      character*80 filename

!
!    Find the limits of the density of states.
!
      
      wrflag = 1
      if (term == 1) then
         emin = 1.0e30_dp
         emax = -emin
         do ia = 1,nd,1
            print *, 'ia,nchain',ia, nchain
!             call assoc_ham(ia)
            call assoc_ab(ia,1)
            print *, 'ia,nchain',ia, nchain
            do la = 1,nchain,1
               if (emin > (lainf(la)-2.0_dp*lbinf(la))) emin = (lainf(la)-2.0_dp*lbinf(la))
               if (emax < (lainf(la)+2.0_dp*lbinf(la))) emax = (lainf(la)+2.0_dp*lbinf(la))
            enddo
         enddo
      elseif ((term == 2).or.(term == 3)) then
         emin = 1.0e30_dp
         emax = -emin
         do ia = 1,nd,1
!             call assoc_ham(ia)
            call assoc_ab(ia,1)
            do la = 1,nchain,1
               if (emin > diag(1,la)) emin = diag(1,la)
               if (emax < diag(lchain(la)+1,la)) emax = diag(lchain(la)+1,la)
            enddo
         enddo
         emin = emin - etail
         emax = emax + etail
      endif
      close(50)

!
!    Open the output files.
!

      filename = genfile(1:lengfn)//'.dos'
      open(unit = 10,file = filename,status = 'NEW')

      filename = genfile(1:lengfn)//'.idos'
      open(unit = 11,file = filename,status = 'NEW')

!
!    Write out the density of states and the integrated density of states.
!
      write(10,'(5x,"efermi",/,g13.5,/,5x,"e",13x,"states")') lef
      write(11,'(5x,"efermi",/,g13.5,/,5x,"e",13x,"states")') lef
      
      maxdos = -999.0_dp
      delte = (emax-emin)/real(ne-1, dp)
      e = emin
      do i = 1,ne,1
         wrflag = 1
         idos = 0.0_dp
         dos = 0.0_dp
         do ia = 1,nd,1
!             call assoc_ham(ia)
            call assoc_ab(ia,1)

            do la = 1,nchain,1
               if (term == 1) then
                  dos = dos - wt(la)*aimag(g00(cmplx(e, kind=dp), & 
     &                        arec(0:lchain(la),la),brec(0:lchain(la)+1,la), & 
     &                        lchain(la),lainf(la), & 
     &                        lbinf(la)))/pi
                  idos = idos + numelsrt(e,la)*wt(la)
               elseif ((term == 2).or.(term == 3)) then
                  dos = dos - wt(la)*aimag(g00(cmplx(e,kt, kind=dp), & 
     &                        arec(0:lchain(la),la),brec(0:lchain(la)+1,la), & 
     &                        lchain(la),lainf(la), & 
     &                        lbinf(la)))/pi
                  idos = idos + numelnull(e,lchain(la),mrec, & 
     &                                    diag(1,la), & 
     &                                    eigvec(1,1,la),kt) & 
     &                                   *wt(la)
               endif
            enddo
         enddo
!          write(10,'(2(x,G13.5))') e,0.5d0*dos/dble(nd)
!          write(11,'(2(x,G13.5))') e,0.5d0*idos/dble(nd)
         write(10,'(2(x,G13.5))') e,dos/real(nd, dp)
         write(11,'(2(x,G13.5))') e,idos/real(nd, dp)
         
         e = e + delte
         if (dos > maxdos) maxdos=dos
         close(50)
      enddo

! !
! !     Write Fermi energy to .dos file.
! !
!       write(10,'(A)') '&'
!       write(10,'(2G13.5)') lef, 0.0d0
!       write(10,'(2G13.5)') lef, maxdos + 1.0d0


!
!    Close the output files.
!

      close(10)
      close(11)

      end

