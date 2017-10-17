
   module mod_tail
      
      use mod_precision
      use mod_const
      
      implicit none
      
      private
            
      type :: tail_t

   ! t:    
   !     0: no tail
   !     1: aug poly
   !     2: mult type I:  f^{\{n,m\}}(x) = (1-x)^{m+1} \sum_{i=0}^n \binomial{m+i}{m} x^i
   !        n: derivatives = 0 on the left side (r1) <-- not yet, 2,2 is hardcoed for now
   !        m: derivatives = 0 on the right side (rc) <-|
   !     3: mult type II: f^{\{n,m\}}(x) = (1-x^n)^m
   !       ~ n-1 devivatives = 0 on left side (r1)
   !         m: derivatives = 0 on the right side (rc)

         integer :: t = 0
         real(dp) :: r1, rc, n, m
         real(dp), allocatable :: cache(:)
         logical :: mlt = .false., vltl = .false.
         procedure(tail_i), pointer, nopass :: f=>null()
      end type tail_t
      
      
      abstract interface
         subroutine tail_i(r,t,d, f)
            import dp, tail_t
            real(dp), intent(in) :: r
            type(tail_t), intent(in) :: t
            integer, intent(in) :: d
            real(dp), intent(out) :: f(0:d)
         end subroutine tail_i
      end interface
   
      public :: tail_t, tail_i, polytail, binmtail, dstptail, nulltail
   contains
   
   subroutine polytail(r,t,d, f)
      real(dp), intent(in) :: r
      type(tail_t), intent(in) :: t
      integer, intent(in) :: d
      real(dp), intent(out) :: f(0:d)
      real(dp) :: x, omx, omx2, omx3, pl, pl1

      if (d > 2) then
         print *, 'Required derivative higher than supported',d
         stop
      end if
      
      x = (r - t % r1) * t % cache(1)
      
      omx = 1 - x
      omx2 = omx*omx
      omx3 = omx2*omx
      pl = (t % cache(3)*x + t % cache(4))*x + t % cache (5)
      f(0) = omx3 * pl
      
      if (d>0) then
         pl1 = 2*t % cache(3)*x + t % cache(4)
         f(1) = (-3*omx2*pl + omx3*pl1) * t % cache(1)
         if (d>1) f(2) = (6*omx*pl - 6*omx2*pl1 + omx3 * 2 * t % cache(3)) * t % cache(1)**2
      end if
   end subroutine polytail   
   
   
   subroutine binmtail(r,t,d, f)
      real(dp), intent(in) :: r
      type(tail_t), intent(in) :: t
      integer, intent(in) :: d
      real(dp), intent(out) :: f(0:d)
      real(dp) :: x, omx, omx2, omx3
      
      x = (r - t % r1) * t % cache(1)
      omx = 1-x
!       omx  = (t % rc - r) * t % cache(1)
      omx2 = omx*omx
      omx3 = omx2*omx
      
      if (d > 2) then
         print *, 'Required derivative higher than supported',d
         stop
      end if

      f(0) = (6*x*x + 3*x + 1)*omx3
      if (d>0) then
         f(1) = -30*x*x*omx2*t % cache(1)
         if (d>1) f(2) = -60*(omx - x)*x*omx*t % cache(1)*t % cache(1)
      end if
      
   end subroutine binmtail   
   
   
   
   subroutine dstptail(r,t,d, f)
      real(dp), intent(in) :: r
      type(tail_t), intent(in) :: t
      integer, intent(in) :: d
      real(dp), intent(out) :: f(0:d)
      real(dp) :: x, dx, r1, rc, n, m, ircmr1, xn, oxn, ixoxn
      
      r1 = t%r1
      rc = t%rc
      n  = t%n
      m  = t%m
      ircmr1 = t % cache(1)
      x = (r - r1)*ircmr1
      
      if (d > 2) then
         print *, 'Required derivative higher than supported',d
         stop
      end if
      
      xn = x**n
      oxn = 1-xn
      
      f(0) = oxn**m
      
      if (d>0) then
         ixoxn = 1.0_dp/(x*oxn)
         f(1) = -n*m*(f(0)*xn)*ixoxn*ircmr1         
         if (d>1) f(2) = f(1)*ixoxn*(n*(1-m*xn) - oxn)*ircmr1 ! or the commented out:
!        if (d>1) f(2) = f(1)*ixoxn*(n - 1 + (1-m*n)*xn)*ircmr1 ! or the one above
      end if
      
   end subroutine dstptail
   
   
   subroutine nulltail(r,t,d, f)
      real(dp), intent(in) :: r
      type(tail_t), intent(in) :: t
      integer, intent(in) :: d
      real(dp), intent(out) :: f(0:d)

      print *, 'nulltail. exitting..'
      print *, 'r, t % t:', r, t%t
      
      
      f(0) = 1.0_dp
      f(1:d) = 0.0_dp
      stop
      
   end subroutine nulltail
      
   
   
   end module mod_tail