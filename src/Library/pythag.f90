 
       function pythag(a,b)
          use mod_precision

           real(dp) :: pythag
!       implicit none
      real(dp) :: a,b
!
!     finds sqrt(a**2+b**2) without overflow or destructive underflow
!
      real(dp) :: p,r,s,t,u
      p = max(abs(a),abs(b))
      if (p  ==  0.0_dp) go to 20
      r = (min(abs(a),abs(b))/p)**2
!      write(6,'("a = ",g14.7,"  b = ",g14.7)') a, b
   10 continue
         t = 4.0_dp + r
         if (t  ==  4.0_dp) go to 20
         s = r/t
         u = 1.0_dp + 2.0_dp*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
!         pythag = sqrt(a*a+b*b)
      end
