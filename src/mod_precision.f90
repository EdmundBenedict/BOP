
    module mod_precision
        implicit none

    !       integer, parameter :: sp = kind(1.0)
    integer, parameter :: dp = 8, qp = 16    !kind(1.0_dp)
!     integer, parameter :: dp = 8, qp = 8    !kind(1.0_dp)

    
    
    end module mod_precision    


    module mod_clock
        use mod_precision
        implicit none
        
        integer(8) :: cr,cm
    contains

        subroutine get_clock_prec()
        call system_clock(count_rate = cr, count_max = cm )
        end subroutine get_clock_prec

    end module mod_clock

