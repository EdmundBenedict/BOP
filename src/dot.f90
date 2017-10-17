 
       function dot(v1,v2,n)
          use mod_precision


!
!    This is a routine to find the dot product between two vectors.
!

      
      implicit none
      real(dp) :: dot

!
!    Declare the simple variables.
!
      
      
      integer i,n
      real(8) :: ddot
!       real    :: sdot

!
!    Declare the arrays.
!

      real(dp) :: v1(n)
      real(dp) :: v2(n)

!
!    Evaluate the dot product.
!

!       dot = 0.0d0
!       do i = 1,n,1
!          dot = dot + v1(i)*v2(i)
!       enddo
!        print *,'n: ',n
      if (n<512) then
        
        dot = dot_product(v1,v2)
!         dot = sum(v1*v2)
      else
        dot = ddot(n,v1,1,v2,1)
      endif
      
      end

