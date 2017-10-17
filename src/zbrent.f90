    
    
    
    function numeldif(ef)
        use mod_precision
        use mod_all_scalar, only : locc
        implicit none
        
        integer :: flag
        
        real(dp), intent(in) :: ef
        real(dp) :: numeldif
        
        integer, save :: counter = 0
        
        interface
            function numel(flag,ef)
                use mod_precision
                implicit none
                integer, intent(inout) :: flag
                real(dp), intent(in) :: ef
                real(dp) :: numel
            end function numel
            
        end interface
        
        
        flag = 1
        numeldif = numel(flag,ef) - locc
        counter = counter + 1
        
        print *, 'counter:', counter
    
    end function numeldif
    
    
    
    function zbrent(x1,x2,f1,f2,tol)
    use mod_precision
    implicit none

    interface
        function numeldif(ef)
            use mod_precision
            implicit none
            real(dp), intent(in) :: ef
            real(dp) :: numeldif
        end function numeldif
    end interface


    real(dp), intent(in) :: x1,x2,f1,f2,tol
    real(dp) :: zbrent



    integer, parameter :: itmax=100
    real(dp), parameter :: eps=epsilon(x1)
    integer :: iter
    real(dp) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
    a=x1
    b=x2
    fa = f1
    fb = f2
    c=b
    fc=fb
    do iter=1,itmax
        if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
            c=a
            fc=fa
            d=b-a
            e=d
        end if
        if (abs(fc) < abs(fb)) then
            a=b
            b=c
            c=a
            fa=fb
            fb=fc
            fc=fa
        end if
        tol1=2.0_dp*eps*abs(b)+0.5_dp*tol
        xm=0.5_dp*(c-b)
        if (abs(xm) <= tol1 .or. abs(fb) < tol) then
            zbrent=b
            return
        end if
        if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
            s=fb/fa
            if (a == c) then
                p=2.0_dp*xm*s
                q=1.0_dp-s
            else
                q=fa/fc
                r=fb/fc
                p=s*(2.0_dp*xm*q*(q-r)-(b-a)*(r-1.0_dp))
                q=(q-1.0_dp)*(r-1.0_dp)*(s-1.0_dp)
            end if
            if (p > 0.0) q=-q
            p=abs(p)
            if (2.0_dp*p  <  min(3.0_dp*xm*q-abs(tol1*q),abs(e*q))) then
                e=d
                d=p/q
            else
                d=xm
                e=d
            end if
        else
            d=xm
            e=d
        end if
        a=b
        fa=fb
        b=b+merge(d,sign(tol1,xm), abs(d) > tol1 )
        fb=numeldif(b)
        
        write(6,'(/,"it:",i3)') iter
        write(6,'("ef,nf-ne:",2(x,f22.14))') b, fb
    end do
    write(6,"('zbrent: exceeded maximum iterations')")
    zbrent=b
    end function zbrent
