
   subroutine writexyz(fln,nm)
   
      use mod_const
      use mod_all_scalar
      use mod_io
      
      implicit none
   
      character(len=*), intent(in) :: nm, fln
   
      integer :: i, u
   
      include 'Include/Atom.array'
      include 'Include/PosVel.array'

      
      
      u = get_new_unit()
      open(u, file=trim(fln))
            
      write(u, '(i7)')  nd + ninert
      write(u, '(a)') nm
      
      do i = 1, nd
!          write(u, '(a,x,3(x,es12.5)))') symb(z(i)), ad(1:3,i)
         write(u, '(a,x,3(x,es12.5))') 'Ac', ad(1:3,i)
      end do
      
      do i = 1, ninert
!          write(u, '(a,x,3(x,es12.5)))') symbi(zinert(i)), adinert(1:3,i)
         write(u, '(a,x,3(x,es12.5))') 'In', adinert(1:3,i)
      end do

   
      call close_unit(u)
   
   
   
   
   end subroutine writexyz