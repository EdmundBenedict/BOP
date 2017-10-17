   module mod_chi

   use mod_precision, only : dp
   use mod_const
   use mod_all_scalar, only : forces
   
   implicit none
 
   private

!
!    General Bond Order expansion information.
!

!*** The bond order
      real(dp), public :: bo(mxnstat,mxnstat), dbo(mxnstat,mxnstat)

!*** The susceptibilities. !*** The derivatives of the susceptibilities.
      real(dp), pointer, public :: chia(:,:)=>null(),  chib(:,:)=>null(), &
      &                           dchia(:,:)=>null(), dchib(:,:)=>null()

      real(dp), public, allocatable :: dq2chia(:) ! chiaia(mxnd)



      public :: free_chi, assoc_chi, init_chi


      type chi_arch_t
         real(dp) :: chia(0:mterm,mchain), chib(1:mterm,mchain), dchia(0:mterm,mchain), dchib(1:mterm,mchain)
      end type


      type(chi_arch_t), private, target, allocatable :: chi_arch(:,:)

      contains
    
      subroutine init_chi(a1,a2,nsp)
         integer, intent(in) :: a1,a2,nsp
!          integer :: i
!          
!          if (.not. allocated(chi_arch) ) then
!                allocate(chi_arch(mpmap(iproc)+1:mpmap(iproc+1),nsp))
!          else
!                if ((mpmap(iproc+1) - mpmap(iproc)) /= size(chi_arch,1)) then
!                   deallocate(chi_arch)
!                   allocate(chi_arch(mpmap(iproc)+1:mpmap(iproc+1),nsp))
!                end if
!          end if
!          
!          if (.not. allocated(dq2chia)) then
!             allocate(dq2chia(mpmap(iproc)+1:mpmap(iproc+1)))
!          else
!             if ((mpmap(iproc+1) - mpmap(iproc)) /= size(dq2chia)) then
!                 deallocate(dq2chia)
!                 allocate(dq2chia(mpmap(iproc)+1:mpmap(iproc+1)))
!             end if
!          end if
!   

          call free_chi()
          allocate(chi_arch(a1:a2,nsp))
          allocate(dq2chia (a1:a2))
         
      end subroutine init_chi
         

      subroutine assoc_chi(ia,isp)

         
         integer,intent(in) :: ia, isp
         
         chia  => chi_arch(ia,isp) % chia  
         chib  => chi_arch(ia,isp) % chib  
         if (forces) then
            dchia => chi_arch(ia,isp) % dchia 
            dchib => chi_arch(ia,isp) % dchib 
         end if
!          if(ia>2)then
!          chia  => chi_arch(2) % chia  
!          chib  => chi_arch(2) % chib  
!          dchia => chi_arch(2) % dchia 
!          dchib => chi_arch(2) % dchib 
!          end if
!          

      end subroutine assoc_chi

         
      
      subroutine free_chi()

         if (allocated(chi_arch)) deallocate(chi_arch)
         if (allocated(dq2chia )) deallocate(dq2chia)
      end subroutine free_chi
      


   end module mod_chi