

   subroutine writecell(fln)
   
   use mod_precision
   use mod_const
   use mod_all_scalar
   use mod_conf 
   use mod_atom_ar
   use mod_io
   
   implicit none
   
   include 'Include/PosVel.array'
   include 'Include/Atom.array'
   
   character(len=*), intent(in) :: fln
   
   type(cell_t), pointer :: cell
   
   real(dp) :: inva(3,3)
   
   integer :: i, u
   
   
   cell => rtc % cell
   
   do i = 1,3
      cell % a(1:3,i) = a(1:3,i)/lena
   end do
   
   cell % lena = lena
   
   inva = a
   call inv3x3(inva)
   
   do i = 1, nd
      call mul3x3(inva, ad(1,i), cell % d(i) % crds)
      cell % d(i) % de = de(i)
   end do

   if (mag) then
      do i = 1, nd
         cell % d(i) % mg = mg(i)   
      end do
   end if
   
   do i = 1, ninert
      call mul3x3(inva, adinert(1,i), cell % dinert(i) % crds)
   end do
   
   
   
   u = get_new_unit()
   open(u,file=trim(fln),action='write')
   call print_cell(cell,'',u)
   call close_unit(u)
   
   
   end subroutine writecell