   module mod_ham

   use mod_precision, only : dp
   use mod_const
   
   implicit none
 
   private
!
!    Hamiltonian information.
!

!*** Hamiltonian matrix.
      real(dp), public, pointer :: h(:,:,:,:) => null()

!*** PCLUSTER is a set of pointers to the start of the cluster lists.
      integer, public, pointer :: pcluster(:) => null()

!*** CLUSTER is the array of indices for the atoms in each cluster.
      integer, public, pointer :: cluster(:) => null()

!*** Index to position in cluster list from atom index.
      integer, public, pointer :: decipher(:) => null()

! !
! !    Declare the common blocks.
! !

!       common /hamint/pcluster,cluster,decipher
!       common /hamdp/h,subh,grad
      public :: free_ham, assoc_ham, init_ham

! To not initialise any of the huge arrays from the type because gfortran 4.7 frowns duting runtime
      type, private :: ham_arch_t
         integer :: pcluster(0:mrec+3), cluster(mxcls), decipher(mxtotnd)
         real(dp) :: h(mxnstat,mxnstat,mxnnb,mxcls)
      end type


      type(ham_arch_t), private, target, allocatable :: ham_arch(:)

      contains
    
      subroutine init_ham(a1,a2)
         integer, intent(in) :: a1,a2
         integer :: i

!          if (.not. allocated(ham_arch) ) then
!                allocate(ham_arch(mpmap(iproc)+1:mpmap(iproc+1)))
!          else
!                if ((mpmap(iproc+1) - mpmap(iproc)) /= size(ham_arch)) then
!                   deallocate(ham_arch)
!                   allocate(ham_arch(mpmap(iproc)+1:mpmap(iproc+1)))
!                end if
!          end if
         
         call free_ham()
         allocate(ham_arch(a1:a2))
         
         do i = a1, a2
            ham_arch(i) % decipher = 0
         end do
         
         

      end subroutine init_ham
         

      subroutine assoc_ham(ia)

         
         integer,intent(in) :: ia
         
         pcluster => ham_arch(ia) % pcluster
         cluster  => ham_arch(ia) % cluster
         decipher => ham_arch(ia) % decipher
         h        => ham_arch(ia) % h


      end subroutine assoc_ham

         
      
      subroutine free_ham()
         if (allocated(ham_arch)) deallocate(ham_arch)
      end subroutine free_ham
      
      


   end module mod_ham
