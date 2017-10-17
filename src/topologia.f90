   module topologia

   implicit none

   integer :: nproc, iproc, master, aproc ! aproc : number of active processors in cases when there are more than necessary

   integer, allocatable :: mpmap(:)
   
   contains

   subroutine sdistrib(n, p, m, a)
      integer, intent(in) :: n, p
      integer, intent(out) :: m(0:p), a
      integer :: i, mdl

      m(0) = 0
      m(1:p) = n/p
      mdl = mod(n,p)
      if (mdl>0) m(2:mdl+1) = m(2:mdl+1) + 1
      a = min(p,n+1)
      do i=1,p
          m(i) = m(i) + m(i-1)
      end do
!       print *, 'mpmap:', mpmap
!       print *, 'aproc:', a
   end subroutine sdistrib

   end module topologia

! 
! program pr
!    use distrib
!    implicit none
! 
!    integer :: n=13,p=8
!    integer :: m(0:19)
! 
! !    m = 0
!    call sdistrib(n,p,m)
! 
!    print *, n,p
!    print *,m(0:p)
!    
! 
! 
! end program pr

