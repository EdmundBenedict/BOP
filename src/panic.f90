 
      subroutine panic()
          use mod_precision


!
!    This is a routine to leave the program in the event of a major error.
!

      implicit none

      write(6,'(/''Fatal error found. Emergency exit requested.'')')

      call dump()

      call erasab()

      call close_all_files()

      stop

      end



      subroutine close_all_files()
          use mod_precision

          implicit none
        
          integer :: u
        
          do u = 1,4,1
              call close_unit_safely(u)
          end do
        
          do u = 7,100,1
              call close_unit_safely(u)
          end do
        

      contains
            subroutine close_unit_safely(u)
                use mod_precision

                implicit none
                integer :: u
                logical :: ostat
                
                inquire(unit=u,opened=ostat) 
                if (ostat) close(u)
            end subroutine close_unit_safely
      end subroutine close_all_files
      

      
      
      
      
      
      
      
      