    
    function femag()
        
        use mod_precision
        use mod_const        
        use topologia
        use mod_atom_ar, only : mg
        
        implicit none
        
        real(dp) :: femag
        integer :: ia
        
        include 'Include/Atom.array'
        
        
        femag = 0.0_dp
        
        do ia = mpmap(iproc)+1, mpmap(iproc+1)
            femag = femag + istn(z(ia))*(mg(ia)*mg(ia))
        end do
        
        femag = 0.25_dp*femag

    end function femag