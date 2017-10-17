 
      subroutine erasab()
          use mod_precision


!
!    This is a routine to erase the temporary file used to store data
!     used for finding the fermi energy.
!
      use ab_io
      implicit none

!
!    Include constants.
!

!      include "Include/ALL.const"

!
!    Include scalars.
!

!      include "Include/ALL.scalar"

!
!    Declare the simple variables.
!
!       logical, optional, intent(in) :: hstore
      character*80 filename 
      character*1024 scratch

      logical ex

! !
! !    Erase the file.
! !
! !       if (present(hstore)) then
! !         if (.not. hstore) then
! !             call free_ab()
! !         else 
! !             call get_environment_variable('BOPS_SCRATCH', scratch)
! !             if (len(trim(scratch))<1) scratch = 'tmp'
! !             
! !             filename = genfile(1:lengfn)//'.ab'
! !             inquire(file=trim(scratch)//filename,exist=ex)
! !             if (ex) then
! !                 open(unit=50,file=trim(scratch)//'/'//filename,form='UNFORMATTED', & 
! !             &        status='OLD')
! !                 close(50,status='DELETE')
! !             endif
! !             
! !         end if
! ! 
! !       else
! !         print *,'debug 8'
!         call free_ab()
! !       end if



      end

