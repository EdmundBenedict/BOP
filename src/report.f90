 
      subroutine report(flag)
          use mod_precision
          use mod_all_scalar
          use mod_const
          use ab_io
          use mod_chi
          use mod_funptr
!
!    This is a routine to summarize the coefficients evaluated.
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
!       include "Include/BondOrder.array"
!      include "Include/Force.array"
!      include "Include/Hamilt.array"
!      include "Include/Moment.array"
      include "Include/NebList.array"
!      include "Include/NotSRT.array"
      include "Include/PosVel.array"
!      include "Include/RecurCoef.array"
!      include "Include/SRT.array"

!
!    Declare the simple variables.
!

      integer ib,lb,nsttb,nlb,nstta
      integer ja,ja0,nla
      integer i,j,nmax
      integer ia,la
      integer flag,wrflag
      real(dp) :: mom2, mom3, mom4

!
!    Test to see what is requested.
!

      if (flag == 0) then

!
!       Do nothing.
!

      elseif (flag == 1) then

!
!       Ensure that the onsite and intersite recursion coefficients
!        obey the sum rule.
!

         call sumrule()

      elseif (flag == 2) then

!
!       Write out recursion coefficients.
!

         wrflag = 1
         do ia = 1,nd,1
            call rdab(wrflag,ia)
            write(9,'(''Atom # '',I4)') ia
            call states(z(ia),nla,nstta,llista)
            do la = 1,nchain,1
               nmax = lchain(la)
               write(9,'(''A('',I2,'') = '',G22.15,'' B('',I2,'') = '', & 
     &                   G22.15)')(i,arec(i,la),i,brec(i,la), & 
     &                   i=0,nmax,1)
               write(9,'(''AINF = '',G22.15,'' BINF = '',G22.15)') & 
     &                   lainf(la),lbinf(la)
            enddo
            call states(z(ia),nla,nstta,llista)
            do la = 1,nstta,1
               ja = aptr(ia)
               ja0 = ja
               do while(bptr(ja) /= eol)
                  ib = bptr(ja)
                  if (ib /= ia) then
                     write(9,'(''Neighbor Atom # '',I4)') ib
                     call states(z(ib),nlb,nsttb,llistb)
                     do lb = 1,nsttb,1
                        write(9,'(''LA = '',I1,'' LB = '',I1)') la,lb
                        write(9,'(''DA('',I2,'') = '',G22.15, & 
     &                            '' DB('',I2,'') = '',G22.15)') & 
     &                        (i,darec(i,la,lb,ja-ja0+1), & 
     &                         i,dbrec(i,la,lb,ja-ja0+1), & 
     &                                  i=0,nmax,1)
                     enddo
                  endif
                  ja = ja + 1
               enddo
            enddo
         enddo
         close(50)

      elseif (flag == 3) then

!
!       Write out eigenvalues and vectors for no terminator.
!

         if ((term == 2).or.(term == 3)) then
            
            do ia = 1,nd,1
               do la = 1,nchain,1
                  nmax = lchain(la)
                  write(9,'(''Eigenvalues:'')')
                  write(9,'(6G13.5)') (diag(i,la),i=1,nmax+1,1)
                  write(9,'(''Eigenvectors:'')')
                  write(9,'(6G13.5)') ((eigvec(j,i,la), & 
     &                                 j=1,nmax+1,1),i=1,nmax+1,1)
               enddo
            enddo
         endif

      elseif (flag == 4) then

!
!       Write out susceptibilities.
!

         wrflag = 1
         do ia = 1,nd,1
            call rdab(wrflag,ia)
            call assoc_chi(ia,1)
            write(9,'(''Atom # '',I4)') ia
            do la = 1,nchain,1
               nmax = lchain(la)
               if (term == 1) then
                  if (chi_meth == 1) then
                     call getchisrt(lef,nmax+1,la)
                  else
                     call getchisrt_c(lef,nmax+1,la)
                  endif
               elseif ((term == 2).or.(term == 3)) then
                  call getchinull(lef,nmax,mrec,chia(0,la), & 
     &                            chib(0,la),diag(1,la), & 
     &                            eigvec(1,1,la),kt)
               endif
               write(9,'(''L = '',I4)') la
               write(9,'(''CHIA('',I2,'') = '',G12.5,''  CHIB('',I2, & 
     &                   '') = '',G12.5)')(j,chia(j,la),j, & 
     &                      chib(j,la),j=0,nmax,1)
            enddo
         enddo
         close(50)


      elseif (flag == 5) then

!
!       Write out some recursion coefficients.
!

         wrflag = 1
         do ia = 1,nd,1
            call rdab(wrflag,ia)
            write(9,'(''Atom # '',I4)') ia
            call states(z(ia),nla,nstta,llista)
            write(6,'('' NCHAIN '',I4)') nchain
            do la = 1,nchain,1
               nmax = lchain(la)
               write(9,'(''A('',I2,'') = '',G22.15,'' B('',I2,'') = '', & 
     &                   G22.15)')(i,arec(i,la),i,brec(i,la), & 
     &                   i=0,nmax,1)
               write(9,'(''AINF = '',G22.15,'' BINF = '',G22.15)') & 
     &                   lainf(la),lbinf(la)
            
            mom2 = arec(0,la)**2 + brec(1,la)**2
            mom3 = arec(0,la)**3 + 2*arec(0,la)*brec(1,la)**2 + &
     &             arec(1,la)*brec(1,la)**2 
            mom4 = arec(0,la)**4 + 3*arec(0,la)**2*brec(1,la)**2 + &
     &             2*arec(0,la)*arec(1,la)*brec(1,la)**2 + &
     &             arec(1,la)**2*brec(1,la)**2 +  &
     &             brec(1,la)**2*brec(2,la)**2 + brec(1,la)**4
            write(19,'(" Moments ")')
            write(19,'(''Atom # '',I4)') ia
            write(19,'("2nd, 3rd, 4th ",3G22.15)') mom2,mom3,mom4

            enddo
         enddo
         close(50)


      endif

      end

