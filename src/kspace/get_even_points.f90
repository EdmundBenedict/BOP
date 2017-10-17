    subroutine get_even_points(nv, oc, n,  p)
        use mod_precision

        
        implicit none
        
        integer, intent(in) :: nv(3)
        logical, intent(in) :: oc(3)
        integer, intent(in) :: n
        
        real(8), intent(out) :: p(3,n)
        
        real(8), dimension(3) :: d, s
        integer :: i, j, k
        
        
        if (product(nv) /= n) then
            write(6,'("The passed size does not match the requrements!",/ &
                      & "There is a bug in whatever call get_even_points")')
            call panic()
        endif
        
        d = 1.0_dp/nv
        
        where(oc) 
            s = d*0.5_dp
        elsewhere 
            s = 0.0_dp
        end where
        
        do k = 0, nv(3)-1
            do j = 0, nv(2)-1
                do i = 0, nv(1)-1
                    p(:, k*nv(1)*nv(2) + j*nv(1) + i + 1) = (/ i, j, k /)*d + s
                end do
            end do
        end do
        
        
        
        do i = 1,n
            print '("ek ",3(x,f8.5))', p(:,i)
        end do
        
    
    
    end subroutine get_even_points