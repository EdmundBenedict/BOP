 
   function eff(v,fun,ef0,tol,mxit)
      use mod_precision
      use mod_all_scalar
      use mod_const
      
      implicit none
 
      abstract interface
         function fun_i(x) result (r)
            use mod_precision, only : dp
            real(dp), intent(in) :: x
            real(dp) :: r
         end function fun_i
      end interface
      
      real(dp), intent(in) :: v, ef0, tol
      integer, intent(in) :: mxit
      procedure(fun_i) :: fun
      real(dp) :: eff
      
      
      real(dp) :: efhi,eflo,nf,nflo,nfhi
      real(dp) :: ef,def,efmid,nfmid
      real(qp) :: c,d,dsc
      

      integer :: niter


      real(qp) :: p(3)


      integer :: it
      real(dp) :: ef_prev, dnf_prev, dnf
      real(qp) :: x0,x1,x2,y0,y1,y2
!     integer, save :: counter = 0
      integer :: rfmethod

      

      
      
      !
!    Find limits within which fermi energy lies.
!
!      WRITE(*,'("STARTING eff")')
!      WRITE(*,'("FINDING FERMI")')
      eflo = ef0
      nflo = fun(eflo)
!        print *, 'nflo, eflo ', nflo, eflo
!       
      if (abs(nflo - v) < 1e-10_dp) then
        eff = eflo
        totnia = nflo
        return
! !       else
! !         print *, eflo,nflo, nflo - v
      endif

      def = 1.0_dp
      if (nflo > v) then
         do while(nflo > v)
            efhi = eflo
            nfhi = nflo
            eflo = eflo - def
            nflo = fun(eflo)
!             print *,'going down',eflo,efhi
         enddo
      else
         efhi = eflo
         nfhi = nflo
         do while(nfhi < v)
            eflo = efhi
            nflo = nfhi
            efhi = efhi + def
            nfhi = fun(efhi)
!             print *,'going up',eflo,efhi
         enddo
      endif






      
      
!       print *,grn,eflo,(efhi-eflo),(nflo-v),(nfhi-nflo),endc
      efmid = (efhi*(v-nflo)+eflo*(nfhi-v))/(nfhi-nflo)
!       efmid = eflo - (efhi-eflo)*(nflo-v)/(nfhi-nflo)
!       print *, 'unsafe:', nfhi, nflo, nfhi-nflo, eflo,efhi, efmid
      nfmid = fun(efmid)
!       counter = counter + 1
!         print *, 'end unsafe:' ,nfhi, nflo, nfhi-nflo, efmid
!      PRINT *,'NFMID = ',NFMID

      nf = -1.0_dp
      niter = 1

!        print *,'v:', v

      do while ((abs(nf-v) > nd*tol) .and. (niter <= mxit)) !nd*1.0e-10_dp
!            write(6,'("efi:",x,i3,3(x,g22.14))') niter, eflo,efmid,efhi
!            write(6,'("nfi:",x,i3,3(x,g22.14))') niter, nflo,nfmid,nfhi
!     print *, '      nit:', niter
!          x(1,:) = (/ 1.0_dp, eflo, eflo*eflo /)
!          x(2,:) = (/ 1.0_dp, efmid, efmid*efmid /)
!          x(3,:) = (/ 1.0_dp, efhi, efhi*efhi /)
! 
!          call inv3x3(x)
! 
!          p(1) = x(1,1)*nflo + x(1,2)*nfmid + x(1,3)*nfhi
!          p(2) = x(2,1)*nflo + x(2,2)*nfmid + x(2,3)*nfhi
!          p(3) = x(3,1)*nflo + x(3,2)*nfmid + x(3,3)*nfhi
!          

         
         x0 = real(eflo         , qp)
         x1 = real(efmid - x0   , qp)
         x2 = real(efhi - x0    , qp)
         y0 = real(nflo         , qp)
         y1 = real(nfmid - y0   , qp)
         y2 = real(nfhi - y0    , qp)
         
            
            
         p(3) = x2*y1 - x1*y2               ! x**2
         p(2) = x1*x1*y2 - x2*x2*y1         ! x
         p(1) = (y0-real(v,qp))*x1*x2*(x1-x2)    ! 1
         
         
         if ( abs(x1*x2*(x1-x2)) < 1e-10_qp .and. .not. quiet ) then
            print *, 'sungular mat, det:', x1*x2*(x1-x2)
            write(6,'("xi:",3(x,g22.14))') x0,x1,x2
            write(6,'("yi:",3(x,g22.14))') y0,y1,y2
            write(6,"('p_sing:',3(x,g22.14))") p
         end if
         
!          if (abs(p(3)) < (1.0e-10_dp) .and. ((1000.0_dp*abs(p(3)) < abs(p(2))) .or. (1000.0_dp*abs(p(3)) < abs(p(1)))) ) then
        if ((abs(p(3)) < (1.0e-30_qp) .and. ((abs(p(3)) < abs(p(2))) .or. (abs(p(3)) < abs(p(1)))) ) .or. all(p == 0.0_qp)) then
!             write(6,"('p:',3(x,f22.14))") p
!             print *,viob,'sing->bisect',endc
            if (nfmid > v) then
               ef = 0.5_dp*(eflo+efmid)
            else
               ef = 0.5_dp*(efhi+efmid)
            endif
         else
         
!             c = -p(2)/(2.0d0*p(3))
            c = -0.5_qp*p(2)/p(3) 
            dsc = c*c-p(1)/p(3)
!             write(6,"('c,dsc:',2(x,f22.14))") c,dsc
            if (dsc < 0.0_qp) then
!                print *,viob,'imag->bisect',endc
               if (nfmid > v) then
                  ef = 0.5_dp*(eflo+efmid)
               else
                  ef = 0.5_dp*(efhi+efmid)
               endif
            else
               d = sqrt(dsc)
               ef = real(c+d + x0, dp)
               if ((ef < eflo).or.(ef > efhi)) then
                  ef = real(c-d + x0, dp)
                  if ((ef < eflo).or.(ef > efhi)) then
!                      print *,viob,'nons->bisect',endc
                     if (nfmid > v) then
                        ef = 0.5_dp*(eflo+efmid)
                     else
                        ef = 0.5_dp*(efhi+efmid)
                     endif
                  endif
               endif
            endif
         endif

!         if ((nfhi < v) .or. (nflo > v)) then
!             write(6,'("STOP")')
!             write(6,'("efi:",x,i3,3(x,g22.14))') niter, eflo,efmid,efhi
!             write(6,'("nfi:",x,i3,3(x,g22.14))') niter, nflo,nfmid,nfhi
!             write(6,'("xi:",3(x,g22.14))') x0,x1,x2
!             write(6,'("yi:",3(x,g22.14))') y0,y1,y2
!             write(6,"('p_sing:',3(x,g22.14))") p
!             stop 'cant go on'
!         end if


         nf = fun(ef)
         if (.not. quiet) write (6,"('i,ef,nf:',x,i3,3(x,g22.14))") niter, ef, nf, v
!          counter = counter + 1



        if (nf /= nfmid) then                                            
            if (nf > nfmid .and. nf > v ) then
                efhi = ef
                nfhi = nf
            else if (nfmid > nf .and. nfmid > v ) then
                efhi  = efmid
                efmid = ef
                nfhi  = nfmid
                nfmid = nf
            else if (v > nf .and. v > nfmid ) then
                if (nfmid > nf) then
                    eflo = ef
                    nflo = nf
                else
                    eflo  = efmid
                    efmid = ef
                    nflo  = nfmid
                    nfmid = nf
                end if
            end if
        
        else if ( nf < v) then
            eflo  = ef
            efmid = 0.5_dp*(ef + efhi)
            nflo  = nf
            nfmid = fun(efmid)
        else 
            efhi = ef
            efmid = 0.5_dp*(ef + eflo)
            nfhi  = nf
            nfmid = fun(efmid)
        end if
        
        
!         idebug = 1
!          if (idebug == 1) then
!             idebug = 0
! !             write(6,'("     i,eflo,efmid,efhi,nf:",i0,6(x,f22.14))') niter,eflo,efmid,efhi
! !             write(6,'("     i,nflo,nfmid,nfhi,nf:",i0,6(x,f22.14))') niter,eflo,efmid,efhinf, nf-v
!             write(6,"(a,'i,nf:',x,i3,x,f22.14,a)") grnb, niter, nf,endc
!          endif
         niter = niter + 1

      enddo
      
      eflo = ef
      nflo = nf

      eff = ef
      totnia = nf

!     print *,'counter:', counter
!     end select

      end function eff

