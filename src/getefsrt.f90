 
       function getefsrt(nelec)
          use mod_precision


          use mod_all_scalar

          use mod_const

!
!    This is a function to find the fermi energy for a square root terminator.
!

      implicit none
        
        
        interface  
            function zbrent(x1,x2,f1,f2,tol)
                use mod_precision
                implicit none
                real(dp), intent(in) :: x1,x2,f1,f2,tol
                real(dp) :: zbrent
            end function zbrent
            
            function zriddr(x1,x2,f1,f2,xacc)
                use mod_precision
                implicit none
                real(dp), intent(in) :: x1,x2,f1,f2,xacc
                real(dp) :: zriddr
            end function zriddr
            
            function numeldif(ef)
                use mod_precision
                implicit none
                real(dp), intent(in) :: ef
                real(dp) :: numeldif
            end function numeldif
        end interface  
        
        
      
      
      
      real(dp) :: getefsrt


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
!      include "Include/RecurCoef.array"
!      include "Include/SRT.array"
!      include "Include/Force.array"

!
!    Common blocks introduced by A.Girshick.
!    
 
      include "Include/ag.conn"


!
!    Declare the simple variables.
!

      real(dp) :: nelec,efhi,nf,nflo,nfhi
      real(dp) :: ef,def,efmid,nfmid
      real(qp) :: c,d,dsc
      real(dp) :: numel

      integer niter
      

!
!    Declare local arrays.
!

!       real(dp) :: x(3,3)
      
      real(qp) :: p(3)


      integer :: it
      real(dp) :: ef_prev, dnf_prev, dnf
      real(qp) :: x0,x1,x2,y0,y1,y2
!     integer, save :: counter = 0
    integer :: rfmethod
!
!    Find limits within which fermi energy lies.
!
!      WRITE(*,'("STARTING GETEFSRT")')
!      WRITE(*,'("FINDING FERMI")')

      nflo = numel(eflo)
!        print *, 'nflo, eflo ', nflo, eflo
!       
      if (abs(nflo - nelec) < 1e-10_dp) then
        getefsrt = eflo
        totnia = nflo
        return
! !       else
! !         print *, eflo,nflo, nflo - nelec
      endif

      def = 1.0_dp
      if (nflo > nelec) then
         do while(nflo > nelec)
            efhi = eflo
            nfhi = nflo
            eflo = eflo - def
            nflo = numel(eflo)
!             print *,'going down',eflo,efhi
         enddo
      else
         efhi = eflo
         nfhi = nflo
         do while(nfhi < nelec)
            eflo = efhi
            nflo = nfhi
            efhi = efhi + def
            nfhi = numel(efhi)
!             print *,'going up',eflo,efhi
         enddo
      endif






!     print *, blu,'eflo,efhi,nflo,nfhi:',eflo,efhi,nflo,nfhi,endc

!     rfmethod = 3
! 
!     select case (rfmethod)
!         case(0)
! 
!       ef = zriddr(eflo,efhi,nflo-nelec,nfhi-nelec,1.0e-8_dp)
!       nf = numel(ef)
! 
!       eflo = ef
!       nflo = nf
! 
!       getefsrt = ef
!       totnia = nf
! 
! 
!     case(1)
! 
! 
! 
! 
!       ef = zbrent(eflo,efhi,nflo-nelec,nfhi-nelec,1.0e-8_dp)
!       nf = numel(ef)
! 
!       eflo = ef
!       nflo = nf
! 
!       getefsrt = ef
!       totnia = nf
! 
!     case(2)
! 
! 
! 
! 
! 
! 
! 
! 
! 
! 
!       
!       if (abs(nflo - nelec) < abs(nfhi - nelec)) then
!         ef_prev = efhi
!         dnf_prev = nfhi - nelec
!         ef = eflo
!         dnf = nflo - nelec      
!       else      
!         ef_prev = eflo
!         dnf_prev = nflo - nelec
!         ef = efhi
!         dnf = nfhi - nelec  
!       end if
! 
!       
!       
!       
!       
! !       call panic()
!       
!       do it = 1, 100
!           if (it > 1) dnf = numeldif(ef)
!                     
!           if (.not. quiet) then
!             write(6,'(/,"it:",i3)') it
!             write(6,'("ef,nf:",4(x,f22.14))') ef_prev, ef, dnf_prev, dnf
!           end if
!           
!           if (abs(dnf) <= 1.0e-10_dp) exit
!           def = dnf * (ef - ef_prev) / ( dnf - dnf_prev )
!           
!           ef_prev = ef
!           dnf_prev = dnf
!           ef = ef - def
!         
!       end do
! 
!       eflo = ef
!       nflo = dnf + nelec
! 
!       getefsrt = ef
!       totnia = dnf + nelec
! 
!     return
! 
! 
! 
! 
! 
!     case(3)



      
      
!       print *,grn,eflo,(efhi-eflo),(nflo-nelec),(nfhi-nflo),endc
      efmid = (efhi*(nelec-nflo)+eflo*(nfhi-nelec))/(nfhi-nflo)
!       efmid = eflo - (efhi-eflo)*(nflo-nelec)/(nfhi-nflo)
!       print *, 'unsafe:', nfhi, nflo, nfhi-nflo, eflo,efhi, efmid
      nfmid = numel(efmid)
!       counter = counter + 1
!         print *, 'end unsafe:' ,nfhi, nflo, nfhi-nflo, efmid
!      PRINT *,'NFMID = ',NFMID

      nf = -1.0_dp
      niter = 1

!        print *,'nelec:', nelec

      do while ((abs(nf-nelec) > nd*1.0e-10_dp) .and. (niter <= 100))
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
         p(1) = (y0-real(nelec,qp))*x1*x2*(x1-x2)    ! 1
         
         
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
            if (nfmid > nelec) then
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
               if (nfmid > nelec) then
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
                     if (nfmid > nelec) then
                        ef = 0.5_dp*(eflo+efmid)
                     else
                        ef = 0.5_dp*(efhi+efmid)
                     endif
                  endif
               endif
            endif
         endif

!         if ((nfhi < nelec) .or. (nflo > nelec)) then
!             write(6,'("STOP")')
!             write(6,'("efi:",x,i3,3(x,g22.14))') niter, eflo,efmid,efhi
!             write(6,'("nfi:",x,i3,3(x,g22.14))') niter, nflo,nfmid,nfhi
!             write(6,'("xi:",3(x,g22.14))') x0,x1,x2
!             write(6,'("yi:",3(x,g22.14))') y0,y1,y2
!             write(6,"('p_sing:',3(x,g22.14))") p
!             stop 'cant go on'
!         end if


         nf = numel(ef)
!         write (6,"('i,ef,nf:',x,i3,2(x,g22.14))") niter, ef, nf
!          counter = counter + 1



        if (nf /= nfmid) then                                            
            if (nf > nfmid .and. nf > nelec ) then
                efhi = ef
                nfhi = nf
            else if (nfmid > nf .and. nfmid > nelec ) then
                efhi  = efmid
                efmid = ef
                nfhi  = nfmid
                nfmid = nf
            else if (nelec > nf .and. nelec > nfmid ) then
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
        
        else if ( nf < nelec) then
            eflo  = ef
            efmid = 0.5_dp*(ef + efhi)
            nflo  = nf
            nfmid = numel(efmid)
        else 
            efhi = ef
            efmid = 0.5_dp*(ef + eflo)
            nfhi  = nf
            nfmid = numel(efmid)
        end if
        
        
!         idebug = 1
!          if (idebug == 1) then
!             idebug = 0
! !             write(6,'("     i,eflo,efmid,efhi,nf:",i0,6(x,f22.14))') niter,eflo,efmid,efhi
! !             write(6,'("     i,nflo,nfmid,nfhi,nf:",i0,6(x,f22.14))') niter,eflo,efmid,efhinf, nf-nelec
!             write(6,"(a,'i,nf:',x,i3,x,f22.14,a)") grnb, niter, nf,endc
!          endif
         niter = niter + 1

      enddo
      
      eflo = ef
      nflo = nf

      getefsrt = ef
      totnia = nf

!     print *,'counter:', counter
!     end select

      end

