    function zriddr(x1,x2,f1,f2,xacc)
    
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


    real(dp), intent(in) :: x1,x2,f1,f2,xacc
    real(dp) :: zriddr
    
    integer, parameter :: maxit=60
    real(dp), parameter :: unused=-1.11e30_dp
    integer :: j
    real(dp) :: fh,fl,fm,fnew,s,xh,xl,xm,xnew
    fl=f1
    fh=f2
    if ((fl > 0.0 .and. fh < 0.0) .or. (fl < 0.0 .and. fh > 0.0)) then
        xl=x1
        xh=x2
        zriddr=unused
        do j=1,maxit
            xm=0.5_dp*(xl+xh)
            fm=numeldif(xm)
            s=sqrt(fm**2-fl*fh)
            if (s == 0.0) return
            xnew=xm+(xm-xl)*(sign(1.0_dp,fl-fh)*fm/s)
            if (abs(xnew-zriddr) <= xacc) return
            zriddr=xnew
            fnew=numeldif(zriddr)
            if (fnew == 0.0) return
            if (sign(fm,fnew) /= fm) then
                xl=xm
                fl=fm
                xh=zriddr
                fh=fnew
            else if (sign(fl,fnew) /= fl) then
                xh=zriddr
                fh=fnew
            else if (sign(fh,fnew) /= fh) then
                xl=zriddr
                fl=fnew
            else
                write(6,"('zriddr: never get here')")
            end if
            if (abs(xh-xl) <= xacc) return
            write(6,'(/,"   j:",i3)') j
!             write(6,'("   ef,nf-ne:",2(x,f22.14))') zriddr, fnew
        end do
        write(6,"('zriddr: exceeded maximum iterations')")
    else if (fl == 0.0) then
        zriddr=x1
    else if (fh == 0.0) then
        zriddr=x2
    else
        write(6,"('zriddr: root must be bracketed')")
    end if
    end function
