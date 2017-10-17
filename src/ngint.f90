x0 = max((lainf(la)-2*lbinf(la)-ef)/kt,-10.0_dp)
x1 = min((lainf(la)+2*lbinf(la)-ef)/kt, 10.0_dp)
m = min(mx,max(11,nrec*nint(50*kt)))
m = 2*(m/2)+1
x = x0
dx = (x1-x0)/real(m-1, dp)
do i = 1,m
    zi = cmplx(ef + kt*x, kind=dp)
    y(i) = aimag(g00(zi,arec(0:lchain(la),la), &
                        & brec(0:lchain(la)+1,la),lchain(la), & 
                        & lainf(la),lbinf(la)))*delta(kt*x,kt)
    x = x + dx
enddo





x0 = max((lainf(la)-2*lbinf(la)-ef)/kt,-10.0_dp)
x1 = min((lainf(la)+2*lbinf(la)-ef)/kt, 10.0_dp)
m = min(mx,max(11,nrec*nint(50*kt)))
m = 2*(m/2)+1
x = x0
dx = (x1-x0)/real(m-1, dp)
do i = 1,m
    zi = cmplx(ef + kt*x, kind=dp)
    y(i) = aimag(g00(zi,arec(0:lchain(la),la), &
                    & brec(0:lchain(la)+1,la),lchain(la), & 
                    & lainf(la),lbinf(la)))*entden(kt*x,kt)
    x = x + dx
enddo