 

   subroutine mul3x3(m,v1,v2)
      use mod_precision
      
!       perform    m^t  v1 -> v2
!    v1 shalt not be the same adress as v2


      implicit none

      integer i,j

      real(dp), intent(in) :: m(3,3)
      real(dp), intent(in) :: v1(3)
      real(dp), intent(out) :: v2(3)


      do i = 1, 3
         v2(i) = sum(m(1:3,i)*v1(1:3))
      end do

   end subroutine mul3x3
   
   
   
!    
!    subroutine mul3x3i(t, a, b)
!       
!       use mod_precision
!       
! !       regular matrix multiplication assumed
! !  a b -> b
!    
!    
!    end subroutine mul3x3i

   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
