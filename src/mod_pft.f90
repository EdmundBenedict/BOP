
   module mod_pft
   
      use mod_precision
      use mod_const
      
      implicit none
      
      private
      
      type :: pft_t
         integer :: t = 0, n = 0
         real(dp), allocatable :: a(:), cache(:) 
         procedure(pft_i), pointer, nopass :: f
      end type pft_t


      abstract interface
         subroutine pft_i(r,p,d, f)
            import dp, pft_t
            real(dp), intent(in) :: r
            type(pft_t), intent(in) :: p
            integer, intent(in) :: d
            real(dp), intent(out) :: f(0:d)
         end subroutine pft_i
      end interface
      
      public :: pft_t, tailed_fun, init_pfts
      
   
   contains
   
     
      subroutine gsph(r,p,d, f)
      
!       This is in need of sanitisation and optimisation
      
         real(dp), intent(in) :: r
         type(pft_t), intent(in) :: p
         integer, intent(in) :: d
         real(dp), intent(out) :: f(0:d)
         
         real(dp) :: r0, rd, n, nc, f1, f2
         
!          scl = p % cache(1)
         r0  = p % a(1)
         rd  = p % a(2)
         n   = p % a(3)
         nc  = p % a(4)
         
!          r**-n*exp(-n*(r/rd)**nc)
         
!          -n*r**-n*exp(-n*(r/rd)**nc)*(1 + nc*(r/rd)**nc)/r
         
!          print *, n,nc
         if (n == 0.0_dp) then
            f1 = 1.0_dp
            if (nc == 0.0_dp) then
               f2 = 1.0_dp
            elseif (nc == 1.0_dp) then
               f2 = exp((r0/rd)-(r/rd))
!                print *, 'r0-r',r0-r
!                print *, "fv",(r0/rd)-(r/rd)
            else
               f2 = exp((r0/rd)**nc-(r/rd)**nc)
            endif
         else
   !          f1 = exp(n*log(r0/r))
            f1 = (r0/r)**n
            if (nc == 0.0_dp) then
               f2 = 1.0_dp
            elseif (nc == 1.0_dp) then
               f2 = exp(n*((r0/rd)-(r/rd)))
            else
               f2 = exp(n*((r0/rd)**nc-(r/rd)**nc))
            endif
         endif

         f(0) = f1*f2

         if (d < 1) return
         
         f1 = f(0)

         if (n == 0.0_dp) then
            if (nc == 0.0_dp) then
               f2 = 0.0_dp
            elseif (nc == 1.0_dp) then
               f2 = -1.0_dp/rd
            else
               f2 = -(nc/r)*((r/rd)**nc)
            endif
         else
            if (nc == 0.0_dp) then
               f2 = -n/r
            elseif (nc == 1.0_dp) then
               f2 = -n*(1.0_dp/r + 1.0_dp/rd)
            else
               f2 = -(n/r)*(1.0_dp+nc*(r/rd)**nc)
            endif
         endif

         f(1) = f1*f2
         
         if (d < 2) return
         
         f1 = f(0)
         
         if (n == 0.0_dp) then
            if (nc == 0.0_dp) then
               f2 = 0.0_dp
            elseif (nc == 1.0_dp) then
               f2 = 1.0_dp/rd**2
            else
               f2 = (nc/r**2)*(r/rd)**nc*(nc*(r/rd)**nc-nc+1)
            endif
         else
            if (nc == 0.0_dp) then
               f2 = n*(n+1)/r**2
            else

               f2 = n/r**2 * ( 1 + n*(1+nc*(r/rd)**nc)**2 + &
      &                     nc*(r/rd)**nc - nc**2*(r/rd)**nc ) 

            endif
         endif

         f(2) = f1*f2
         
         
         
      end subroutine gsph
      
   
   
      subroutine spln(r,p,d, f)
         real(dp), intent(in) :: r
         type(pft_t), intent(in) :: p
         integer, intent(in) :: d
         real(dp), intent(out) :: f(0:d)
! There is room of optimisation here too. may be precompute the diffs then power separately?!         
         integer :: i, n
         real(dp) :: w(0:2), v(0:2), rkr, invr, g
         
         n = p % n / 2 ! number of pairs
!          order of decreasing Rk assumed
         
         f = 0.0_dp
            
         if (p%t == 5 .and. r < p % a(2*n)) then
            rkr = p % a(2*n) - r
            v(2) = p % a(2*n-1) * rkr
            v(1) = v(2) * rkr
            v(0) = v(1) * rkr
            if (d > 0) v(1:2) = v(1:2) * [-3, 6]
            
            invr = 1.0_dp/r
            w(0) = exp(rkr) * invr
            g = invr + 1
            if (d > 0) w(1) = -w(0)*g
            if (d > 1) w(2) = w(0)*(2*invr*g+1.0_dp) 
            
            f(0) = v(0)*w(0)
            if (d > 0) f(1) = -(v(1)*w(0) + v(0)*w(1))/3.0_dp
            if (d > 1) f(2) = (v(2)*w(0) + 2*v(1)*w(1) + v(0)*w(2))/6.0_dp
            
            n = n-1
         end if
         
         do i = 1, n
            if (r >= p % a(2*i)) exit
            f(0) = f(0) + p % a(2*i-1) * (p % a(2*i) - r) ** 3
         end do
         
         if (d < 1) return

         do i = 1, n
            if (r >= p % a(2*i)) exit
            f(1) = f(1) + p % a(2*i-1) * (p % a(2*i) - r) ** 2
         end do
         
         f(1) = -3*f(1)
         
         if (d < 2) return

         do i = 1, n
            if (r >= p % a(2*i)) exit
            f(2) = f(2) + p % a(2*i-1) * (p % a(2*i) - r)
         end do
         
         f(2) = 6*f(2)
         
      end subroutine spln
      
         
   
      subroutine rexs(r,p,d, f)
!       The form is a*r**(-n)*exp(-m*(r-rd))
!      the params series should be supplied in quads:
!      a1 n1 m1 rd1   a2 n2 m2 rd2 ....
      
         real(dp), intent(in) :: r
         type(pft_t), intent(in) :: p
         integer, intent(in) :: d
         real(dp), intent(out) :: f(0:d)
         
         integer :: i, n
         real(dp) :: ir, t1, t2
         
         
         n = (p % n / 4) * 4
         
         ir = 1.0_dp/r
         
         f = 0.0_dp
         do i = 1 , n, 4
            t1 = p % a(i) * ir ** p % a(i+1) * exp(- p % a(i+2) * (r - p % a(i+3)))
            f(0) = f(0) + t1
            if (d>0) then
               t2 = p % a(i+1)*ir + p % a(i+2)
               f(1) = f(1) + t1*t2
               if (d>1) f(2) = f(2) + t1 * (t2*t2 + p % a(i+1) * ir * ir)
            end if
         end do
         
         if (d > 0) f(1) = -f(1)
         
      end subroutine rexs
      
      
      subroutine plnm(r,p,d, f)
         real(dp), intent(in) :: r
         type(pft_t), intent(in) :: p
         integer, intent(in) :: d
         real(dp), intent(out) :: f(0:d)
         
         f = 0.0_dp
         stop 'plain polynomial potential not implemented yet'
         
      end subroutine plnm
      
      subroutine nullfun(r,p,d, f)
         real(dp), intent(in) :: r
         type(pft_t), intent(in) :: p
         integer, intent(in) :: d
         real(dp), intent(out) :: f(0:d)

         print *, 'nullfun. exitting..'
         print *, 'r, p % t:', r, p%t
         
         f(0) = 1.0_dp
         f(1:d) = 0.0_dp
         stop
         
      end subroutine nullfun
      
         
      subroutine init_pfts(p,t)
         use mod_tail
         
         type(pft_t), intent(inout) :: p
         type(tail_t), intent(inout) :: t
         real(dp) :: fv(0:2), rcmr1
         integer :: sz, i
         
         
   !       Still not sure whether it is not better to leave these selections in the reading part in mod_conf
         select case(t%t) 
            case (0)
               t%f => nulltail
            case (1)  
               t%f   => polytail
            case (2)
               t%mlt = .true.

               t%f   => binmtail
            case (3)   
               t%mlt = .true.
               t%f   => dstptail
            case default
               print *,'invalid tail choice ', t%t
               stop 
         end select
         
         
         select case (p%t)
            case (0)
               p%f => nullfun
            case (1)
               p%f => gsph
            case (2)
               p%f => rexs
            case (3)
               p%f => plnm
            case (4)
               p%f => spln
            case (5)
               p%f => spln
            case default
               print *,'invalid fun choice ', p%t
               stop    
         end select
         
         
         
   !       if (p % t == 1) then
   !          if (.not. allocated(p%cache)) allocate (p % cache()) !
   !       end if
         
         rcmr1 = t % rc - t % r1
         sz = 2
         if (.not. t % mlt) sz = 5
         if (.not. allocated(t%cache)) allocate(t%cache(sz))         
         t%cache(1) = 1.0_dp/rcmr1
         t%cache(2) = rcmr1
         if (.not. t % mlt) then
            t % vltl = .true.
            call tailed_fun(0.5_dp*(t % r1 + t % rc), p, t, 0, fv) ! the r = halfway bwt r1 and rc here is just so we are surely in the tail. other than this it is irrelevant.
            t % vltl = .false.
         end if
         
         if (p%t == 4 .or. p%t == 5 ) then
!          check whether the spline pairs are in order of descending Rk
            do i = 2, (p % n / 2)
               if (p % a(2*i) > p % a(2*(i-1))) stop 'Order the spline parameters in pairs Ak,Rk by decreasing Rk and try again ..'
            end do
         end if
         
         
      end subroutine init_pfts
                  
      subroutine tailed_fun(r,p,t,d, f)
         use mod_tail
         
         real(dp), intent(in) :: r
         type(pft_t), intent(in) :: p         
         type(tail_t), intent(inout) :: t
         integer, intent(in) :: d
         real(dp), intent(out) :: f(0:d)
         real(dp) :: tv(0:d)
         real(dp) :: fv(0:2)  ! these are not up to d as the tail need the high derivatives on the boundaties in any case
         

         procedure(pft_i), pointer :: pf
         procedure(tail_i), pointer :: tf
         pf => p % f
         tf => t % f

         if (r < t % rc) then
!             print *, 'tail_fun r < rc'
            if (r <= t % r1 .or. t % t == 0) then
!                print *, 'tail_fun r <= r1 or t == 0 ', r, t % r1 , t % t, p % t
!                print *, 'tail_fun assoc pf', associated(p%f)
               call pf(r, p, d, f)
            else
!                print *, 'tail_fun r > r1 .and. t /= 0 ', r , t % r1 , t % t
               if (t % mlt) then
!                   print *, 'tail_fun mlt' 
                  call pf(r, p, d, fv)
                  call tf(r, t, d, tv)
                  f(0) = fv(0)*tv(0)
                  if (d>0) f(1) = fv(1)*tv(0) + fv(0)*tv(1)
                  if (d>1) f(2) = fv(2)*tv(0) + 2*fv(1)*tv(1) + fv(0)*tv(2)
               else
!                   print *, 'tail_fun nomlt'
                  if (t % vltl) then
!                      print *, 'poly tail calc'
                     call pf(t % r1, p, 2, fv)
                     t % cache(3) = (0.5_dp * fv(2) * t % cache(2) + 3 * fv(1)) * t % cache(2) + 6 * fv(0)
                     t % cache(4) = fv(1) * t % cache(2) + 3 * fv(0)
                     t % cache(5) = fv(0)
                  end if
                  call tf(r, t, d, f)
 
               end if
            end if
         else if (t % rc <= r ) then
            f = 0.0_dp
         end if
      end subroutine tailed_fun
      
      
    
   
   end module mod_pft
   
   
   
   
   
