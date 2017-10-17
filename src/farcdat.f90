 
      subroutine farcdat(p,mu,kt,epout,spout)
          use mod_precision


!
!******************************************************************
!     I=\int P(E)f(E) =  -T \sum_i P( cxp(i)*T+mu )*sp(i)
!******************************************************************
!
!    This routine passes back the values of the poles and the factors
!     needed to carry out the integral shown above.
!

      implicit none

!
!    Declare the simple variables.
!

      complex(dp) :: epout

      real(dp) :: mu,kt,spout

      integer p
      integer nmx

!
!    Define parameters.
!

      parameter (nmx= 15)

!
!    Declare the arrays.
!

      complex(dp) :: cxp(nmx)

      real(dp) :: sp(nmx)

!
!    Assign data values.
!
!    Fifteen point approximation.
!

      data cxp( 1)/ (-4.58902802629    ,  1.35980605083 )   /
      data cxp( 2)/ (-4.56715604516    ,  1.34275624602 )  /
      data cxp( 3)/ (-2.99891216014    ,  2.97176715264 )   /
      data cxp( 4)/ (-2.89092370441    ,  3.00341937426 )   /
      data cxp( 5)/ (-1.47681348666    ,  3.26093144117 )   /
      data cxp( 6)/ (-.352368246614    ,  2.47906680939 )   /
      data cxp( 7)/ ( .000000000000e+00,  2.39188426776 )   /
      data cxp( 8)/ ( .000000000000e+00,  3.87129391754 )   /
      data cxp( 9)/ ( .000000000000e+00,  3.01389125368 )   /
      data cxp(10)/ ( .352368246614    ,  2.47906680939 )   /
      data cxp(11)/ ( 1.47681348666    ,  3.26093144117 )   /
      data cxp(12)/ ( 2.89092370441    ,  3.00341937426 )   /
      data cxp(13)/ ( 2.99891216014    ,  2.97176715264 )   /
      data cxp(14)/ ( 4.56715604514    ,  1.34275624603 )   /
      data cxp(15)/ ( 4.58902802631    ,  1.35980605082 )   /

      data sp( 1)/-1.000000    /
      data sp( 2)/ 1.000000    /
      data sp( 3)/-1.000000    /
      data sp( 4)/ 1.000000    /
      data sp( 5)/-1.000000    /
      data sp( 6)/ 1.000000    /
      data sp( 7)/ 1.000000    /
      data sp( 8)/-1.000000    /
      data sp( 9)/ 1.000000    /
      data sp(10)/ 1.000000    /
      data sp(11)/-1.000000    /
      data sp(12)/ 1.000000    /
      data sp(13)/-1.000000    /
      data sp(14)/ 1.000000    /
      data sp(15)/-1.000000    /

!
!    Seven point approximation.
!

!     DATA CXP(1)/(-1.771028710319321D0,3.55783255876122D0)/
!     DATA CXP(2)/(-0.3006180464402131D0, 2.323529036976914D0)/
!     DATA CXP(3)/(0.D0, 7.170464262981526D0)/
!     DATA CXP(4)/(0.D0, 2.583135129624091D0)/
!     DATA CXP(5)/(0.D0, 7.284992346676701D0)/
!     DATA CXP(6)/(0.3006180464402222D0, 2.323529036976921D0)/
!     DATA CXP(7)/(1.771028710319322D0,  3.557832558761221D0)/

!     DATA SP(1)/-1D0/,SP(2)/1D0/,SP(3)/1D0/,SP(4)/1D0/,SP(5)/-1D0/
!     DATA SP(6)/1D0/,SP(7)/-1D0/

!
!    Set up values.
!

      epout = cmplx(kt, kind=dp)*cxp(p)+cmplx(mu, kind=dp)
      spout = sp(p)

      end

