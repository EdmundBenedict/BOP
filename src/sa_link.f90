    module sa_link
        use mod_precision
        use mod_conf
        implicit none
        
        type cell_patch
            type(cell_t) :: cell
            type(ens_t) :: en
            type(elastic_const_t) :: elcon
        end type cell_patch
        
        
        type(run_conf_t), private, pointer :: current_rtconf => null()
        type(cell_t), private, pointer :: current_cell => null()
        
    contains
    
    
        subroutine bop_en(rtconf, cellpkg)
            use mod_all_scalar, only : eprom, ebond, epair, eatom, uent,  etot
            type(run_conf_t), intent(inout), target :: rtconf
            type(cell_patch), intent(inout), target :: cellpkg
            real(dp) :: real_num_atoms
            
            
!             print *,' calculating ', cellpkg % cell % name
            
            
            rtconf % mvflag = -1
            rtconf % cell => cellpkg % cell
            
            current_rtconf => rtconf
            current_cell => rtconf % cell
            
            call bop(rtconf)
            
            cellpkg % en  = rtconf % avre
            
        end subroutine bop_en
        
        
        
        subroutine bop_el(rtconf, cellpkg)
            type(run_conf_t), intent(inout), target :: rtconf
            type(cell_patch), intent(inout), target :: cellpkg
            real(dp) :: real_num_atoms
            
            rtconf % mvflag = 12
            rtconf % cell => cellpkg % cell
            
            
            current_rtconf => rtconf
            current_cell => rtconf % cell
            
            call bop(rtconf)
            
            cellpkg % elcon % s11 = rtconf % elastic_const_conf % s11
            cellpkg % elcon % s33 = rtconf % elastic_const_conf % s33
            cellpkg % elcon % c11 = rtconf % elastic_const_conf % c11
            cellpkg % elcon % c33 = rtconf % elastic_const_conf % c33
            cellpkg % elcon % c12 = rtconf % elastic_const_conf % c12
            cellpkg % elcon % c13 = rtconf % elastic_const_conf % c13
            cellpkg % elcon % c44 = rtconf % elastic_const_conf % c44
            cellpkg % elcon % c66 = rtconf % elastic_const_conf % c66
            cellpkg % elcon % c1344 = rtconf % elastic_const_conf % c1344
            cellpkg % elcon % c1266 = rtconf % elastic_const_conf % c1266
            
            cellpkg % en  = rtconf % avre        
        
        end subroutine bop_el


        subroutine bop_st(rtconf, cellpkg)
            type(run_conf_t), intent(inout), target :: rtconf
            type(cell_patch), intent(inout), target :: cellpkg
            real(dp) :: real_num_atoms
            
            rtconf % mvflag = 36
            rtconf % cell => cellpkg % cell
            
            
            current_rtconf => rtconf
            current_cell => rtconf % cell
            
            call bop(rtconf)
            
            cellpkg % elcon % s11 = rtconf % elastic_const_conf % s11
            cellpkg % elcon % s33 = rtconf % elastic_const_conf % s33
            
            cellpkg % en = rtconf % avre

        end subroutine bop_st
        
        
        

    
    
    
    end module sa_link
