 
       function avennb(aptr,nd)
          use mod_precision


!
!    This a routine to evaluate the average number of nearest neighbours.
!

      implicit none
      real(dp) :: avennb

!
!    Declare the simple variables.
!

      integer nd,i,sum

!
!    Declare the arrays.
!

      integer aptr(nd+1)

!
!    Evaluate the average cluster size.
!

      sum = 0
      do i = 1,nd,1
         sum = sum + aptr(i+1)-aptr(i)-1
      enddo
      avennb = real(sum, dp)/real(nd, dp)

      end

