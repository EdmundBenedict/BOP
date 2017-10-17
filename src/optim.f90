    
    
    module opti
    
    use mod_precision
    use mod_const
    
    implicit none
    
    contains
    
         
         
         
         
         
    subroutine nrm1(n, x, dist, fun, a1, )     
         
         integer, intent(in) :: n
         real(dp), intent(inout), dimension(n) :: x, dist, a1
         
         
         do                         
            x(:n) = x(:n) - sum((/ (de(ia)*zc(z(ia)), ia = 1,nd) /))/sumz
            
            call spnbop(n, x, dq_o, dqde_o, nsp, eelct, mg_i, mg_o)
            
            if (dqmax <= qerr .or. it == mxit) exit

            de(:nd) = de(:nd) - dq/chia(:nd)

            it = it + 1
         end do
    
    
    
    
    
    
    
    
    end module opti