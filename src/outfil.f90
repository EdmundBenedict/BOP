 

      subroutine outfil(flag)
          use mod_precision
          use mod_all_scalar
          use mod_const
          use mod_srt
          use ab_io
          use mod_g0n, only : g00
          use mod_ham
          use mod_chi
          use mod_kspace
          use mod_funptr
          use mod_atom_ar, only: dq, btype
          
!
!    This is a routine to write data to the output file.
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
      include "Include/BEC.array"
!       include "Include/BondOrder.array"
      include "Include/Force.array"
!       include "Include/Hamilt.array"
!       include "Include/KHamilt.array"
!      include "Include/Misc.array"
!      include "Include/Moment.array"
      include "Include/NebList.array"
!      include "Include/NotSRT.array"
      include "Include/PosVel.array"
!      include "Include/RecurCoef.array"
!      include "Include/Relax.array"
!      include "Include/SRT.array"

!
!    Common blocks introduced by A.Girshick.
!    
 
      include "Include/ag.conn"


!
!    Declare the simple variables.
!

!       complex(dp) :: g00

      real(dp) :: dossr,dosef
      real(dp) :: avennb
      real(dp) :: bsons
      real(dp) :: niter
      real(dp) :: totocc
      real(dp) :: u0
      real(dp) :: tret
      real(dp) :: sumenk,sumea


      integer la,n,nmax
      integer i,j,k,ia,za,ib
      integer flag,wrflag
      integer stat, fcfiton

      character*80 filename

!
!    Setup output file.
!

      if (flag == 1) then

!
!       Write summary of system to output file.
!
         write(9,'(/''Energy shifts.'')')
         write(9,'(8G13.5)') (de(i),i=1,nd,1)
         if (potflg == 2) then
            write(9,'(''*******************************'')')
            write(9,'(''*          B  O  P            *'')')
            write(9,'(''*******************************'')')
         endif

         if (potflg == 4) then
            write(9,'(''*******************************'')')
            write(9,'(''*      K - S  P  A  C  E      *'')')
            write(9,'(''*******************************'')')
         endif

         if (potflg == 2) then

         if (term == 1) then
            write(9,'(''*******************************'')')
            if (momflg == 1) then
               write(9,'(''*USING: AVERAGED MOMENTS      *'')')
               write(9,'(''*       SQUARE ROOT TERMINATOR*'')')
            elseif (momflg == 2) then
               write(9,'(''*USING: LOCALLY ORIENTED ORBITALS*'')')
               write(9,'(''*       SQUARE ROOT TERMINATOR   *'')')
            elseif (momflg == 3) then
               write(9,'(''*USING: SEPERATED ORBITALS    *'')')
               write(9,'(''*       SQUARE ROOT TERMINATOR*'')')
            endif
            if (chi_meth == 1) then
               write(9,'(''*       ANALYTIC INTEGRATION  *'')')
            else
               write(9,'(''*       NUMERICAL INTEGRATION *'')')
            endif
            write(9,'(''*******************************'')')
         elseif (term == 2) then
            if (momflg == 1) then
               write(9,'(''**************************'')')
               write(9,'(''*USING: AVERAGED MOMENTS *'')')
               write(9,'(''*       FINITE TERMINATOR*'')')
               write(9,'(''**************************'')')
            elseif (momflg == 2) then
               write(9,'(''**********************************'')')
               write(9,'(''*USING: LOCALLY ORIENTED ORBITALS*'')')
               write(9,'(''*       FINITE TERMINATOR        *'')')
               write(9,'(''**********************************'')')
            elseif (momflg == 3) then
               write(9,'(''***************************'')')
               write(9,'(''*USING: SEPERATED ORBITALS*'')')
               write(9,'(''*       FINITE TERMINATOR *'')')
               write(9,'(''***************************'')')
            endif
         elseif (term == 3) then
            if (momflg == 1) then
               write(9,'(''*************************'')')
               write(9,'(''*USING: AVERAGED MOMENTS*'')')
               write(9,'(''*       NO TERMINATOR   *'')')
               write(9,'(''*************************'')')
            elseif (momflg == 2) then
               write(9,'(''**********************************'')')
               write(9,'(''*USING: LOCALLY ORIENTED ORBITALS*'')')
               write(9,'(''*       NO TERMINATOR            *'')')
               write(9,'(''**********************************'')')
            elseif (momflg == 3) then
               write(9,'(''***************************'')')
               write(9,'(''*USING: SEPERATED ORBITALS*'')')
               write(9,'(''*       NO TERMINATOR     *'')')
               write(9,'(''***************************'')')
            endif
         endif
   
         write(9,'(/''========================================'')')
         write(9,'(''Number of exact recursion levels   = '',I3)') nbase
         write(9,'(''Total number of recursion levels   = '',I3)') nrec
         write(9,'(''Maximum number of recursion levels = '',I3)') mrec
         write(9,'(''========================================'')')

         write(9,'(/''=========================================='', & 
     &             ''============'')')
         write(9,'(''Number of atoms in central cell         = '',I5)') & 
     &                                                               nd
         write(9,'(''=========================================='', & 
     &             ''============'')')

         endif


      elseif (flag == 2.and.fs_only == 0) then

!
!       Write out electronic information.
!

         call getocc(1,totocc)

         write(9,'(/''=========================================='', & 
     &             ''============'')')
         write(9,'(''Total number of atoms                   = '',I5)') & 
     &                                                            totnd
         write(9,'(''Average number of nearest neighbours    = '', & 
     &         G12.5)') avennb(aptr,nd)
         if (potflg == 2) then
         write(9,'(''Average recursion cluster size          = '', & 
     &         G12.5)') aveclusiz
         endif
         write(9,'(''Total number of electrons per unit cell = '', & 
     &             G12.5)') totocc
         write(9,'(''=========================================='', & 
     &             ''============'')')

!
!       Write out the energies.
!
! 
!          if (potflg == 4) then
!             call iktran(kk,nk,a)
!             do j = 1,nk,1
!                write(9,'(/''K # '',I3,'' : ('',G12.5,'','',G12.5, & 
!      &                    '','',G12.5,'')'')') j,kk(1,j),kk(2,j),kk(3,j)
!                write(9,'(''Eigenvalues and occupancies :'')')
!                do i = 1, kmxh, 3
!                   if (i == kmxh) then
!                      write(9,'(2G12.5,2X)') enk(i,j),occ(i,j)
!                   else if (i ==  kmxh-1) then
!                      write(9,'(2(2G12.5,2X))') enk(i,j),occ(i,j), & 
!      &                    enk(i+1,j),occ(i+1,j)
!                   else
!                      write(9,'(3(2G12.5,2X))') enk(i,j),occ(i,j), & 
!      &                    enk(i+1,j),occ(i+1,j), & 
!      &                    enk(i+2,j),occ(i+2,j)
!                   endif
!                enddo
!             enddo
!          endif

!
!       Calculate the density of states at the Fermi level, and the onsite
!        band energy.
!
         
         if (potflg == 2) then

         bsons = 0.0_dp
         dosef = 0.0_dp
         dossr = 0.0_dp
         wrflag = 1
         do ia = 1,nd,1
            call rdab(wrflag,ia)
            call assoc_chi(ia,1)
            bsons = bsons + ebsos(ia,lef)
            do la = 1,nchain,1
               if (term == 1) then
                  dosef = dosef - wt(la)*aimag(g00(cmplx(lef, kind=dp), & 
     &                        arec(0:lchain(la),la),brec(0:lchain(la)+1,la), & 
     &                     lchain(la),lainf(la),lbinf(la)))/pi
               elseif ((term == 2).or.(term == 3)) then
                  dosef = dosef - wt(la)*aimag(g00(cmplx(lef,kt, kind=dp), & 
     &                        arec(0:lchain(la),la),brec(0:lchain(la)+1,la), & 
     &                     lchain(la),lainf(la),lbinf(la)))/pi
               endif
               nmax = lchain(la)
               if (term == 1) then
                  call getchisrt(lef,nmax+1,la)
               elseif ((term == 2).or.(term == 3)) then
                  call getchinull(lef,nmax,mrec,chia(0,la), & 
     &                            chib(0,la),diag(1,la), & 
     &                            eigvec(1,1,la),kt)
               endif
               do n = 0,nmax,1
                  dossr = dossr + wt(la)*chia(n,la)
               enddo
            enddo
         enddo
         dosef = dosef/real(nd, dp)
         dossr = dossr/real(nd, dp)
         close(50)

         endif

!
!       Write out the energies.
!

         write(9,'(/''========================================='', & 
     &              ''==============='')')
         write(9,'(''Fermi energy                      = '',G12.5, & 
     &             '' eV'')') lef
         write(9,'(''Density of states at Fermi energy = '',G12.5, & 
     &             ''/eV/atom'')') dosef
         write(9,'(''Electron temperature              = '',G12.5, & 
     &             '' eV'')') kt
         write(9,'(''========================================='', & 
     &              ''==============='')')
         write(9,'(/''========================================='', & 
     &              ''==============='')')
         write(9,'(''Occupancy             = '',G12.5)') locc
         write(9,'(''On site energy        = '',G22.15,'' eV'')') & 
     &                                                            eprom
         write(9,'(''Bond energy           = '',G22.15,'' eV'')') & 
     &                                                            ebond
         write(9,'(''Band structure energy = '',G22.15,'' eV'')') & 
     &                                                      eprom+ebond
         write(9,'(''Repulsive energy      = '',G22.15,'' eV'')') & 
     &                                                            epair
         write(9,'(''Atomic band energy    = '',G22.15,'' eV'')') & 
     &                                                           eatom
         write(9,'(''Electron entropy      = '', & 
     &             G22.15,'' eV'')') uent
         if ((mvflag >= 2).and.(mvflag <= 5)) then
            ke = 1.5_dp*kb*temp*real(nd-1, dp)
            write(9,'(''Kinetic energy        = '', & 
     &                G22.15,'' eV'')') ke
            etot = eprom+ebond+epair-eatom+ke-uent
         else
            etot = eprom+ebond+epair-eatom-uent
            u0   = etot + 0.5_dp*uent
         endif
         write(9,'(''Total free     energy    = '',G22.15,'' eV'')') & 
     &                                                             etot
         write(9,'(''Total energy (kT=0)      = '',G22.15,'' eV'')') & 
     &                                                             u0
         write(9,'(''Average free     energy  = '',G22.15,'' eV'')') & 
     &                                                    etot/real(nd, dp)
         write(9,'(''Average energy (kT=0)    = '',G22.15,'' eV'')') & 
     &                                                    u0/real(nd, dp)
         write(9,'(''========================================='', & 
     &              ''==============='')')

!
!       Check consistency between on-site recursion method and
!        bond order potential method.
!
         if (potflg == 2) then

         write(9,'(/''========================='', & 
     &              ''======================'')')
         write(9,'(''Consistency check (matrix only):'')')
         write(9,'(''On-site band energy    = '',G22.15)') bsons
         write(9,'(''Inter-site band energy = '',G22.15)') eprom+ebond
         write(9,'(''Number of electrons    = '',G22.15)') locc
         write(9,'(''Occupancy              = '',G22.15)') totnia
         write(9,'(''========================='', & 
     &             ''======================'')')

!
!       Convergence test. The density of states at the fermi energy is
!        calculated two ways, and the results compared.
!

         write(9,'(/''============================================='', & 
     &              ''===='')')
         write(9,'(''Convergence test:'')')
         write(9,'(''DOS at Fermi energy      = '',G22.15)') dosef
         write(9,'(''Sum of susceptibilities  = '',G22.15)') dossr
         write(9,'(''============================================='', & 
     &             ''===='')')


         elseif (potflg == 4) then
!             sumenk = 0.0_dp
!             do j = 1,nk,1
!                do i = 1,kmxh,1
!                   sumenk = sumenk + wtk(j)*enk(i,j)
!                enddo
!             enddo
!             sumea = 0.0_dp
!             do ia = 1,nd,1
!                za = z(ia)
!                do la = 1,nl(za),1
!                   if (llist(la,za) == 0) then
!                      sumea = sumea + (es(za) + de(ia))
!                   elseif (llist(la,za) == 1) then
!                      sumea = sumea + 3.0_dp*(ep(za) + de(ia))
!                   elseif (llist(la,za) == 2) then
!                      sumea = sumea + 5.0_dp*(ed(za) + de(ia))
!                   endif
!                enddo
!             enddo
!             write(9,'(/''========================='', & 
!      &                 ''======================'')')
!             write(9,'(''Consistency check:'')')
!             write(9,'(''Sum of eigenvalues     = '',G22.15)') sumenk
!             write(9,'(''Sum of onsite energies = '',G22.15)') sumea
!             write(9,'(''Number of electrons    = '',G22.15)') locc
!             write(9,'(''Occupancy              = '',G22.15)') totnia
!             write(9,'(''========================='', & 
!      &                ''======================'')')
            continue
         endif



!
!       Write out atomic parameters used.
!
         write(9,'(A)')' '
         write(9,'(A)')'Atomic Parameters.'
         write(9,'(A)')'=================='
         write(9,'(A)') ' Z   Charge     Es       Ep       Ed       C1 & 
     &     C2       C3       C4'
         do i = 0, natype
           za = atype(i)
           write(9,'(A3,2X,F4.1,2X,7(F9.4))')listsymb(za),zc(za),es(za), & 
     &    ep(za),ed(za),embed(1,za),embed(2,za),embed(3,za),embed(4,za)
         enddo

!
!     Write out bond energy parameters.
!
         write(9,'(A)') ' '
         write(9,'(A)') 'Bond energy parameters.'
         write(9,'(A)') '======================='
         write(9,'(A)') '  Z1  Z2  Vsss   Vsps   Vpss   Vpps   Vppp   Vs & 
     &ds   Vdss   Vpds   Vdps   Vpdp   Vdpp   Vdds   Vddp   Vddd   Phi0 & 
     &'
         do i = 0, mxz
            do j = 0, mxz
               if (btype(i,j) /= 0) then
                  ib = btype(i,j)
                  write(9,'(A3,1X,A3,1X,15(F6.3,1X))')  & 
     &                 listsymb(i),listsymb(j),(bnddat(ia,ib),ia=1,15)
               endif
            enddo
         enddo

!
!     Write out bond scaling parameters.
!
         write(9,'(A)') ' '
         write(9,'(A)') 'Bond scaling parameters.'
         write(9,'(A)') '======================='
         write(9,'(A)') '  Z1  Z2   r0     rc     r1    rcut    n    nc & 
     &     d0     dc     d1    dcut    m     mc'
         do i = 0, mxz
            do j = 0, mxz
               if (btype(i,j) /= 0) then
                  ib = btype(i,j)
                  do k=1,15
                  if (bndscl(1,k,ib) /= 0.0) then
                  write(9,'(a, 1x, a, 1x, 12(f6.3, 1x) )')  & 
     &                 listsymb(i),listsymb(j),(bndscl(ia,k,ib),ia=1,12)
                  endif
                  enddo
               write(9,'(A)') '------------------------'
               endif
            enddo
         enddo


!
!       Write out energy shifts.
!

         write(9,'(/''Energy shifts.'')')
         write(9,'(8G13.5)') (de(i),i=1,nd,1)

!
!       Write out errors in local charge.
!

         write(9,'(/''Errors in local charge.'')')
         write(9,'(8G13.5)') (dq(i),i=1,nd,1)


!       Write out contributions to the cohesive energy.

         write( 9, '(/"EBOND = ", F12.6 )' )    ebond / nd
         write( 9, '( "EPAIR = ", F9.6 )' )    epair / nd            
         write( 9, '( "EPROM = ", F9.6 )' )    eprom / nd        
         write( 9, '( "EATOM = ", F9.6 )' )  - eatom / nd         
         write( 9, '( "UENT  = ", F9.6 )' )  - uent  / nd        


      elseif ( flag  ==  2   .and.  fs_only  ==  1 ) then

         write(9,'(''Repulsive energy  = '',G22.15,'' eV'')') & 
     &                                           epair/real(nd, dp)
      

         continue


      elseif (flag == 3) then

!
!       Write out forces.
!

         write(9,'(/''Band structure forces (eV/A):'')')
         write(9,'(3G13.5)') (fbs(1,j),fbs(2,j),fbs(3,j),j=1,nd,1)

         write(9,'(/''Pair potential forces (eV/A):'')')
         write(9,'(3G13.5)') (fpp(1,j),fpp(2,j),fpp(3,j),j=1,nd,1)

         write(9,'(/''Total forces (eV/A).'')')
         write(9,'(3G13.5)') (ftot(1,j),ftot(2,j),ftot(3,j),j=1,nd,1)

!
!     Adding this subroutine to help with the fitting of force constants
!     
!     M. Cawkwell 15th April 2004
!
!         CALL FNDREC(8,'FCFITON',STAT)
!         READ(8,*) FCFITON
!
!         IF (FCFITON.EQ.1) THEN
!            CALL FCFIT()
!         ENDIF
!
!       Write out the block data.
!

         call outblock( 1 )


!
!       Write out thermal averages.
!

         if ((mvflag >= 2).and.(mvflag <= 4).and.(iter > nequil)) then

            niter = real(iter-nequil, dp)
            avge = sume/niter
            sde = sqrt(sumee/niter - avge**2)
            avgt = sumt/niter
            sdt = sqrt(sumtt/niter - avgt**2)
            avgp = sump/niter
            sdp = sqrt(sumpp/niter - avgp**2)
            avgv = sumv/niter
            sdv = sqrt(abs(sumvv/niter - avgv**2))

            write(9,'(/''======================================'')')
            write(9,'(''Thermal averages:'')')
            write(9,'(''Energy      = '',G12.5,'' ('',G9.2,'')'')') & 
     &             avge,sde
            write(9,'(''Temperature = '',G12.5,'' ('',G9.2,'')'')') & 
     &             avgt,sdt
            write(9,'(''Pressure    = '',G12.5,'' ('',G9.2,'')'')') & 
     &             avgp,sdp
            write(9,'(''Volume      = '',G12.5,'' ('',G9.2,'')'')') & 
     &             avgv,sdv
            write(9,'(''======================================'')')

         endif

! !
! !    Write out the elapsed time.
! !
!          call pclock(tret)
!          write(9,'(/''================================='')')
!          write(9,'(''Total CPU time    = '',G12.5,''s'')') tret
!          if ((mvflag >= 1).and.(mvflag <= 5)) then
!             write(9,'(''Time/atom/MD step = '',G12.5,''s'')')  & 
!      &                tret/real(nd, dp)/real(iter-1, dp)
!          endif
!          write(9,'(''================================='')')
! 
! !
! !    Close the output files.
! !


      endif

      end

