
      function vpotenv( x, env, del, delam )
          use mod_precision
!  this function evaluates the environmentally dependent pair interaction
      implicit none
      include 'Include/Envir.array'
      real(dp) :: vpotenv, x, env(4), delam, del(4),d

      if(x  >=  env(4)) then
        vpotenv = 0.0_dp
      else if ( x > del(4) ) then      
        d = (x - del(4))/(env(4) - del(4))
        
        vpotenv = (env(1)/x) * exp(-(env(2) + delam) * (x - 2*env(3))) &
              & * ((6.0_dp*d + 3.0_dp)*d + 1.0_dp) * (1.0_dp - d)**3
!                 & * (((-6.0_dp*d + 15.0_dp)*d -10.0_dp)*d*d*d + 1.0_dp  )
         print "('env: (',f8.4,'/x)*exp(-(',f8.4,'+',f8.4,')*(x-2*',f8.4,'))  x=',f8.4)", env(1),env(2),delam,env(3),x
      else
        vpotenv = (env(1)/x) * exp(-(env(2) + delam) * (x - 2*env(3)))
      end if

      end function vpotenv


!  this function evaluates the force for the spline part 
!  of the environmental term
     
      function fes( x, env, del, delam, ddelam )
          use mod_precision

!
      implicit none
!
      include 'Include/Envir.array'
!
      real(dp) :: env(4), del(4), x, delam, ddelam
      real(dp) :: fes, rpom, dpom, fpom, apom, bpom, fdpom
      real(dp) :: d_fpom, d_fdpom, d_apom, d_bpom
       
      fes = 0.d0
       
      rpom = del(4)
      dpom = rpom - env(4)
      
      fpom = (env(1)/rpom) * exp( - (env(2)+delam) * &
     &     (rpom - 2*env(3)))
      fdpom = - fpom * ( (env(2)+delam) + 1/rpom)
         
      apom = ( fdpom * dpom - 2 * fpom)/(dpom**3.0_dp)
      bpom = ( 3 * fpom - fdpom * dpom)/(dpom**2.0_dp)
    

      d_fpom = (env(1)/rpom) * exp( - (env(2)+delam) * &
     &     (rpom - 2*env(3)))*(-ddelam)*(env(2)+delam)
      d_fdpom = - d_fpom * ( (env(2)+delam) + 1/rpom) &
     &     - fpom * ddelam
 
      d_apom = (d_fdpom * dpom - 2 * d_fpom)/(dpom**3.0_dp)
      d_bpom = ( 3 * d_fpom - d_fdpom * dpom)/(dpom**2.0_dp)      

      fes = d_apom * ( x - env(4))**3.0_dp + &
     &     3 * apom * ( x - env(4))**2.0_dp + &
     &     d_bpom * ( x - env(4))**2.0_dp + &
     &     2 * bpom * ( x - env(4))
      
      return
      end


!******************************************************************************
!
!  The functions below added to the package to create some kind of uniform
!  interface for using these potetnial functions, like one which is used
!  in the relaxation program: every function has 2 parameters, and it
!  decides what to do based on IPICK, which has to be supplied as an input.
!  The reason we are using these functions above is that they work in units
!  of angstroms.
!
!******************************************************************************


      function  chi( r, ipick )
          use mod_precision

!
      implicit none
!
      include 'Include/Envir.array'
!
      real(dp), intent(in) :: r
      integer, intent(in) :: ipick
      real(dp) :: chi
      
      real(dp) :: d, dc, echi
!       real(dp), parameter :: cs(3) = (/ -6.0_dp, 15.0_dp, -10.0_dp /)
      

      select case(ipick)
        case (0)
            if( r >= aenv(4)) then
                chi = 0.0_dp
            else if ( r > adel(4) ) then
                d = (r - adel(4))/(aenv(4) - adel(4))
                chi = adel(1) * exp( - r * adel(2) ) * ((6.0_dp*d + 3.0_dp)*d + 1.0_dp) * (1.0_dp - d)**3
                !(((-6.0_dp*d + 15.0_dp)*d -10.0_dp)*d*d*d + 1.0_dp  )
                print "('lama: ',f8.4,'*exp(',f8.4,'*x)  x=',f8.4)", adel(1), -adel(2), r
            else
                chi = adel(1) * exp( - r * adel(2) ) 
               print '(i3,2(t5,4(x,f12.8),/))', 0, adel, aenv
            end if
        case(1)
            if( r >= benv(4)) then
                chi = 0.0_dp
            else if ( r > bdel(4) ) then
                d = (r - bdel(4))/(benv(4) - bdel(4))
                chi = bdel(1) * exp( - r * bdel(2) ) * ((6.0_dp*d + 3.0_dp)*d + 1.0_dp) * (1.0_dp - d)**3
!                 (((-6.0_dp*d + 15.0_dp)*d -10.0_dp)*d*d*d + 1.0_dp  )
                print "('lamb: ',f8.4,'*exp(',f8.4,'*x)  x=',f8.4)", bdel(1), -bdel(2), r
            else
                chi = bdel(1) * exp( - r * bdel(2) ) 
                  print '(i3,2(t5,4(x,f12.8),/))', 1, bdel, benv
            end if
            
      case default
         write ( 6, '("chi: ERROR IN IPICK CHOICE")' )
         call panic()
      end select
 
      end function chi
      








! first derivative of CHI with respect to R

      function  dchi( r, ipick )
          use mod_precision

!
      implicit none
!
      include 'Include/Envir.array'
!
      integer ipick
      real(dp) :: dchi, r, rpom, dpom, fpom, fdpom, fd2pom
      real(dp) :: apom, bpom, cpom

      dchi = 0.0_dp
       
      if ( ipick  ==  0 )  then
         
         rpom = adel(4)
         if(r >  aenv(4)) goto 100
         if(r >  rpom) then
            dpom = rpom - aenv(4)
            fpom = adel(1)*exp( - rpom * adel(2) )
            fdpom = - adel(2) * fpom
            fd2pom = -adel(2) * fdpom
            
            apom = ( fd2pom * dpom * dpom - 6 * dpom * fdpom + 12 * fpom ) / (2*dpom**5)
            bpom = ( 7 * dpom * fdpom - fd2pom * dpom * dpom - 15 * fpom ) / (dpom**4)
            cpom = ( fd2pom * dpom * dpom - 8 * dpom * fdpom + 20 * fpom) / (2*dpom**3)
            
            dchi  = 5 * apom * ( r - aenv(4) )**4 + &
        &           4 * bpom * ( r - aenv(4) )**3 + &
        &           3 * cpom * ( r - aenv(4) )**2
            
            return
!             GOTO 100
         endif
         dchi = -adel(2)*adel(1)*exp( - r * adel(2) )
         
      else if ( ipick  ==  2 )  then
     
         rpom = bdel(4)
         if(r >  benv(4)) goto 100  
         if(r >  rpom) then
            dpom = rpom - benv(4)
            fpom = bdel(1)*exp( - rpom * bdel(2) )
            fdpom = - bdel(2) * fpom  
            fd2pom = -bdel(2) * fdpom
            
            apom = ( fd2pom * dpom * dpom - 6 * dpom * fdpom + 12 * fpom ) / (2*dpom**5)
            bpom = ( 7 * dpom * fdpom - fd2pom * dpom * dpom - 15 * fpom ) / (dpom**4)
            cpom = ( fd2pom * dpom * dpom - 8 * dpom * fdpom + 20 * fpom) / (2*dpom**3)
            
            dchi  = 5 * apom * ( r - benv(4) )**4 + &
        &           4 * bpom * ( r - benv(4) )**3 + &
        &           3 * cpom * ( r - benv(4) )**2  
            
            return
!             GOTO 100
         endif
         dchi = -bdel(2)*bdel(1)*exp( - r * bdel(2) )
      
      
      else if ( ipick  ==  1 )  then
     
         rpom = abdel(4)
         if(r >  abenv(4)) goto 100  
         if(r >  rpom) then
            dpom = rpom - abenv(4)
            fpom = abdel(1)*exp( - rpom * abdel(2) )
            fdpom = - abdel(2) * fpom  
            fd2pom = -abdel(2) * fdpom
            
            apom = ( fd2pom * dpom * dpom - 6 * dpom * fdpom + 12 * fpom ) / (2*dpom**5)
            bpom = ( 7 * dpom * fdpom - fd2pom * dpom * dpom - 15 * fpom ) / (dpom**4)
            cpom = ( fd2pom * dpom * dpom - 8 * dpom * fdpom + 20 * fpom) / (2*dpom**3)
            
            dchi  = 5 * apom * ( r - abenv(4) )**4 + &
        &           4 * bpom * ( r - abenv(4) )**3 + &
        &           3 * cpom * ( r - abenv(4) )**2
            
            return
!             GOTO 100
         endif
         dchi = -abdel(2)*abdel(1)*exp( - r * abdel(2) )      
      
      else
         write ( 6, '("ERROR IN IPICK CHOICE")' )
         stop
      endif
 100  continue
      
      return
      end



! this function returns the value of parameter lambda_zero

      function getlamzero( ipick )
          use mod_precision

!
      implicit none
!      
      include 'Include/Envir.array'
!
      integer ipick
      real(dp) :: getlamzero

      if ( ipick  ==  0 ) then
         getlamzero = aenv(2)
      elseif ( ipick  ==  1 ) then
         getlamzero = benv(2)
      else
         write ( 6, '("getlamzero: ERROR IN IPICK CHOICE")' )
         call panic()
      endif
      
      return
      end

! this function returns the value of parameter Rcut
        
      function getrcut( ipick )
          use mod_precision

!
      implicit none
!
      include 'Include/Envir.array'
!
      integer ipick
      real(dp) :: getrcut
        
      if ( ipick  ==  0 ) then
         getrcut = aenv(4)
      elseif ( ipick  ==  1 ) then
         getrcut = benv(4)  
      else
         write ( 6, '("getrcut: ERROR IN IPICK CHOICE")' )
         call panic()
      endif
      
      return
      end
        
      function getrc( ipick )
          use mod_precision

!
      implicit none
!
      include 'Include/Envir.array'
!
      integer ipick
      real(dp) :: getrc
        
      if ( ipick  ==  0 ) then
         getrc = aenv(3)
      elseif ( ipick  ==  1 ) then
         getrc = benv(3)  
      else
         write ( 6, '("getrc: ERROR IN IPICK CHOICE")' )
         call panic()
      endif
      
      return
      end
      
      
      
      

      function  venv( r, ipick, delam )
          use mod_precision

!     
      implicit none
!
      include 'Include/Envir.array'
!
      integer ipick
      real(dp) :: r, venv, delam, vpotenv
      external vpotenv
      
!       write(6,'("r = ", f12.6)') r
!       write(6,'("aenv = ", f12.6)') aenv
!       write(6,'("abenv = ", f12.6)') abenv
!       write(6,'("benv = ", f12.6)') benv
!       write(6,'("adel = ", f12.6)') adel
!       write(6,'("abdel = ", f12.6)') abdel
!       write(6,'("bdel = ", f12.6)') bdel
      
      
      if      ( ipick  ==  0 )  then
         venv = vpotenv( r, aenv,  adel, delam )
      else if ( ipick  ==  1 )  then
         venv = vpotenv( r, abenv, abdel , delam )
      else if ( ipick  ==  2 )  then
         venv = vpotenv( r, benv,  bdel, delam )
      else 
         write ( 6, '("ERROR IN IPICK CHOICE")' )
         stop
      endif
      

      return
      end
      
      
      
      
      
      
      
      
      
      function  fenvspl( r, ipick, delam, ddelam )
          use mod_precision

!     
      implicit none
!
      include 'Include/Envir.array'
!
      integer ipick
      real(dp) :: r, fenvspl, delam, ddelam, fes
      external fes
      
      if      ( ipick  ==  0 )  then
         fenvspl = fes( r, aenv,  adel, delam, ddelam ) 
      else if ( ipick  ==  1 )  then
         fenvspl = fes( r, abenv, adel , delam, ddelam )
      else if ( ipick  ==  2 )  then
         fenvspl = fes( r, benv, bdel, delam, ddelam )
      else
         write ( 6, '("ERROR IN IPICK CHOICE")' )
         stop
      endif
      
      return
      end 







     function gsppot(r, bt)
         use mod_precision

        use pair_coefficients
        use mod_gsp

        implicit none

        real(dp), intent(in) :: r
        integer, intent(in) :: bt
        real(dp) :: gsppot
        real(dp) :: y
        integer :: i
        
!         do i=0,2
!         print '(i0,8(x,g12.6))',  i,   ppot(i) % v,       ppot(i) % r0, ppot(i) % rc, ppot(i) % n, ppot(i) % nc, ppot(i) % r1, ppot(i) % rcut
!         end do
        
        if (r <= ppot(bt) % r1) then
            
            gsppot = ppot(bt) % v * gsp(r, ppot(bt) % r0, ppot(bt) % rc, ppot(bt) % n, ppot(bt) % nc)
        
        else if (r >= ppot(bt) % rcut) then
            
            gsppot = 0.0_dp
        
        else
            y = (r - ppot(bt) % r1)/(ppot(bt) % rcut - ppot(bt) % r1)
            gsppot = ppot(bt) % v * gsp(r, ppot(bt) % r0, ppot(bt) % rc, ppot(bt) % n, ppot(bt) % nc) &
                    &  * ((6.0_dp*y + 3.0_dp)*y + 1.0_dp) * (1.0_dp - y)**3
            
        end if    
            
    end function gsppot



    

    function  mpw( ipick )
          use mod_precision

!
      implicit none
!
      include 'Include/Envir.array'
!
      integer ipick
      real(dp) :: mpw

!       MPW = ADEL(3)
      
      select case(ipick)
        case(0)
            mpw = adel(3)
        case(1)
            mpw = bdel(3)
        case default
            write(6,'("Wrong ipick pased to mpw(). Expected 0 or 1, got ",i0,". Exitting!")') ipick
            call panic()
      end select
      
     
      return
      end



