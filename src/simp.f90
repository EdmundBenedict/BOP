 
       function simp(h,y,n)
          use mod_precision

!****************************************************************
!	THE SUBROUTINE CALCULATES NDIM INTEGRALS
!	BETWEEN X1 AND XN OF FUNCTION Y ( N.LE.NDIM
!	 XN=X1+(N-1)*H ) FOR THE SPECIAL CASE:
!	NDIM=2*ND+1 Y(1)=0 . Z(2*I+1) IS CALCULATED
!	BY THE SIMPSON'S METHOD 
!	     B
!	SIMP=$  dX*Y(X),
!	     A
!****************************************************************
!     implicit real*8(a-h,o-z)
    implicit none
    real(dp) :: simp
    integer, intent(in) :: n
    real(8), intent(in) :: h, y(n)
    
    
    integer :: k
    real(8) :: c3
     
!
    if (.not. ((n/2+n/2) < n)) then
        write(0,*)'<SIMP> WORKS ONLY IF N=2*ND+1, I.E. IS ODD'
        stop
    end if
!
    c3=h/3.0_dp
    simp = 0.0_dp
    if(n /= 1) then
        do k=2,n-1,2
            simp = simp + 4*y(k) + 2*y(k+1)
        end do
        simp=(simp+y(1)-y(n))*c3
    end if
    end
