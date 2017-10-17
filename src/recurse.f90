

   subroutine recurse(near,lim)
      use mod_precision
      use mod_const
      use mod_all_scalar, only : momflg, term, chi_meth, nsp, nbase, quiet
      use mod_srt
      use mod_ham
      use ab_io
      use topologia, only : iproc, mpmap
      use mod_clock

      implicit none
            
      logical, intent(in) :: near ! calculate only nearest neighbour zeta or full cluster    
      real(dp) :: zeta(0:2*mrec+2,mxnstat,mxcls,mxnstat)
      integer :: nclus,lim ! lim allows recurse to run independently of mpmap
      
      include "Include/Atom.array"
      
      integer :: ia, la, nmax, i, j, ifail,isp,in,fn
      
      integer(8) :: c1, c2
      
      call system_clock(c1)

      if (lim == 0) then
          in = mpmap(iproc)+1
          fn = mpmap(iproc+1)
      else
          in = 1
          fn = lim
      endif
      
      
      do isp = 1,nsp
         do ia = in, fn
   !          if (ia>2) exit
            call assoc_ab(ia,isp)
            call assoc_ham(ia)
            
            if (near) then
               nclus = pcluster(2)-1
            else
               nclus = pcluster(nbase+1)-1
            end if
      
!             print *, 'nclus', nclus
!             print *, 'pcluster', pcluster
!             print *, 'cluster', cluster(:nclus)
            
      
            call bldhdo(isp)

   !          Evaluate the recursion coefficients and their derivatives.
   !

            if (momflg == 1) then
               call precab(ia, zeta, near)               
            elseif ((momflg == 2).or.(momflg == 3)) then
               call  recab(ia, zeta, near)
            endif

            call dabdl(ia, zeta, nclus)
   !
   !          Find the band limits for square root terminator.
   !
   !       WRITE(6,'("Find the band limits for square root term")')

            if ((term == 1).or.(term == 2)) then
               do la = 1,nchain
                  call locabinf(arec(0,la), brec(0,la), lchain(la), lainf(la), lbinf(la))
               enddo
   !
   !          Pad out recursion coefficients for finite terminator.
   !
   !        WRITE(6,'("Pad out recursion coefficients for finite")')

               if (term == 2) then
                  call fintrm(ia, nclus)
                  lbinf(:nchain) = 0.0_dp
               endif

            elseif (term == 3) then
               lbinf(la:nchain) = 0.0_dp
            endif



   !
   !          Put in truncator.
   !

            if (momflg == 1) then
               call trncav(ia, nclus)
            elseif ((momflg == 2).or.(momflg == 3)) then
               call trncno(ia, nclus)
            endif

   !
   !          Find roots of polynomials.
   !
   !       WRITE(6,'("Find roots of polynomials")')

            if (term == 1) then
               if (chi_meth == 1) then
                  do la = 1,nchain
                     call getroots(arec(0,la),brec(0,la), & 
      &                             lchain(la),lchain(la)+1, & 
      &                             lainf(la),lbinf(la), & 
      &                             root(1,la),f1(1,la), & 
      &                             f2(1,la),forder(la), & 
      &                             salpha,sbeta,sp,f,fnag,work1, & 
      &                             mrec)
                  enddo
               endif
            elseif ((term == 2).or.(term == 3)) then
               do la = 1,nchain
                  nmax = lchain(la)
                  do i = 0,nmax
                     diag(i+1,la) = arec(i,la)
                     work1(i+1) = brec(i,la)
                  enddo
                  do i = 1,nmax+1
                     do j = 1,nmax+1
                        eigvec(j,i,la) = 0.0_dp
                     enddo
                     eigvec(i,i,la) = 1.0_dp
                  enddo
                  ifail = 0
   !                 CALL F02AMF(NMAX+1,MACPREC,DIAG(1,LA),WORK1,
   !    +                        EIGVEC(1,1,LA),MREC+2,IFAIL)
                  call tql2(mrec+2,nmax+1,diag(1,la),work1, eigvec(1,1,la),ifail)
               enddo
            endif
         end do
      end do
      
      call system_clock(c2)
      
      if (.not. quiet) print *, 't(recurse@root):', real(c2-c1,dp)/real(cr,dp)

   end subroutine recurse
