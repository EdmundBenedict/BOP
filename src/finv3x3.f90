    function finv3x3(m)
        use mod_precision

        
!         Just a functional interface to inv3x3
    
        real(8), intent(in) :: m(3,3)
        
        real(8) :: finv3x3(3,3)
    
        finv3x3(:,:) = m(:,:)
    
        call inv3x3(finv3x3)
    
    
    end function finv3x3