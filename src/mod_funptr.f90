
   module mod_funptr

   use mod_precision
   
   implicit none
   
   public

   
   abstract interface
      subroutine s_r_i_i(ef,nmax,la)
         import dp
         real(dp), intent(in) :: ef
         integer, intent(in) :: nmax, la
      end subroutine s_r_i_i
      
      function f_i_r(ia,ef)
         import dp
         integer, intent(in) :: ia
         real(dp), intent(in) :: ef
         real(dp) :: f_i_r
      end function f_i_r

      function f_r_i(ef,la)
         import dp
         integer, intent(in) :: la
         real(dp), intent(in) :: ef
         real(dp) :: f_r_i
      end function f_r_i

      function f_i_r_i(ia,ef,isp)
         import dp
         real(dp), intent(in) :: ef
         integer, intent(in) :: ia, isp
         real(dp) :: f_i_r_i
      end function f_i_r_i
   end interface
   
   procedure(s_r_i_i), pointer :: getchisrt => null(), dchisr => null()
   procedure(s_r_i_i), private :: getchisrt_a, getchisrt_c, dchisr_a, dchisr_c
   
   procedure(f_i_r_i), pointer :: feprom => null()
   procedure(f_i_r_i), private :: epromavg, eprommix, epromnoavg
   
   procedure(f_i_r  ), pointer :: ebsos     => null()
   procedure(f_i_r  ), private :: ebsos_a, ebsos_c
   
   procedure(f_r_i  ), pointer :: numelsrt  => null()
   procedure(f_r_i  ), private :: numelsrt_a, cnel

   contains
            
   subroutine select_feprom(momflg)
      integer, intent(in) :: momflg
      
      if      (momflg == 1) then
            feprom => epromavg        !          Averaged moments.
      else if (momflg == 2) then
            feprom => eprommix        !          Mixed basis.
      else if (momflg == 3) then
            feprom => epromnoavg       !          Non-averaged basis.
      end if
      
   end subroutine select_feprom
   
   
   subroutine select_sri(chi_meth)
      integer, intent(in) :: chi_meth
      
      if (chi_meth == 1) then
            getchisrt => getchisrt_a
            ebsos     => ebsos_a
            numelsrt  => numelsrt_a
            dchisr    => dchisr_a
      else
            getchisrt => getchisrt_c
            ebsos     => ebsos_c
            numelsrt  => cnel
            dchisr    => dchisr_c
      end if
      
   end subroutine select_sri
   
   
   
   end module mod_funptr