
   subroutine print_forces(l, u)
   
   use mod_precision
   use mod_const
   use mod_all_scalar, only : nd
   
   implicit none
   
   character(len=*), intent(in) :: l
   integer, intent(in) :: u
   
   include 'Include/PosVel.array'
   include 'Include/Force.array'
   
   
   integer :: i
   
   write(u, "(/,a,4x,'total force:',28x,'bandstructure',27x,'classical')") l
   do i = 1, nd
      write(u, '(a,3(x,3(x,es12.5)))') l, ftot(:,i), fbs(:,i), fpp(:,i)
   end do
   
   
   end subroutine print_forces
   
   
   
   