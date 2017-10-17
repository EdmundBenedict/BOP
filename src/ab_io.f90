
    module ab_io
    use mod_precision, only : dp
    use mod_const
    use mod_srt
    
    implicit none
    
    private
    
      include "Include/Atom.array"

!    Recursion coefficient information.


!*** The recursion coefficients.
      real(dp), public, pointer :: arec(:,:)=>null(), brec(:,:)=>null()

!*** Derivatives of recursion coefficients.
      real(dp), public, pointer :: darec(:,:,:,:)=>null(), dbrec(:,:,:,:)=>null()


!*** Local values of termination coefficients. !*** The number of states associated with a linear chain
      real(dp), public, pointer :: lainf(:)=>null(), lbinf(:)=>null(), wt(:)=>null()

!*** The number of linear chains associated with a given site.
      integer, public, pointer :: nchain => null()

!*** The length of the linear chains.
      integer, public, pointer :: lchain(:) => null()


!
!    Non-square root terminator integration information.
!

!*** The parameters used for the null terminator integration scheme.
      real(dp), public, pointer :: diag(:,:)=>null()
      real(dp), public, pointer :: eigvec(:,:,:)=>null()
      


    
    public :: assoc_ab, free_ab,  init_ab

    
    type, private :: cn_arch_t
        integer :: nchain ! quirky name for number of orbitals, mchain is the max allowed
        integer, dimension(mchain) :: lchain, forder
        real(dp), dimension(mchain) :: wt
    end type
    
    type, private :: ab_arch_t
        real(dp), dimension(mchain) :: lainf, lbinf
        real(dp) :: arec(0:2*mrec,mchain),               brec(0:2*mrec+1,mchain), &
                 & darec(0:mrec,mxnstat,mxnstat,mxcls), dbrec(0:mrec+1,mxnstat,mxnstat,mxcls)
        real(dp) :: diag(mrec+2,mchain), eigvec(mrec+2,mrec+2,mchain)
        complex(dp), dimension(2*mterm+2,mchain) :: root, f1, f2
    end type
    
    
    type( ab_arch_t), private, target, allocatable :: ab_arch(:,:)
    type(cn_arch_t), private, target, allocatable :: cn_arch(:)
    
    
    
    contains
    
    subroutine init_ab(a1,a2,nsp)
        
        
        integer, intent(in) :: a1, a2, nsp
!         
!         if (.not. allocated(ab_arch) ) then
!             allocate(ab_arch(mpmap(iproc)+1:mpmap(iproc+1),nsp))
!         else
!             if ((mpmap(iproc+1) - mpmap(iproc)) /= size(ab_arch,1)) then
!                 deallocate(ab_arch)
!                 allocate(ab_arch(mpmap(iproc)+1:mpmap(iproc+1),nsp))
!             end if
!         end if
!         
!         if (.not. allocated(cn_arch) ) then
!             allocate(cn_arch(mpmap(iproc)+1:mpmap(iproc+1)))
!         else
!             if ((mpmap(iproc+1) - mpmap(iproc)) /= size(cn_arch)) then
!                 deallocate(cn_arch)
!                 allocate(cn_arch(mpmap(iproc)+1:mpmap(iproc+1)))
!             end if
!         end if
        
        call free_ab()
        allocate(ab_arch(a1:a2,nsp))
        allocate(cn_arch(a1:a2))
        
    end subroutine init_ab
        

    subroutine assoc_ab(ia,isp)
        use mod_all_scalar, only : term, chi_meth
        
        integer,intent(in) :: ia, isp
        
        
        lainf  => ab_arch(ia,isp) % lainf
        lbinf  => ab_arch(ia,isp) % lbinf
        arec   => ab_arch(ia,isp) % arec
        brec   => ab_arch(ia,isp) % brec
        darec  => ab_arch(ia,isp) % darec
        dbrec  => ab_arch(ia,isp) % dbrec
        
        nchain => cn_arch(ia) % nchain
        lchain => cn_arch(ia) % lchain
        wt     => cn_arch(ia) % wt
        forder => cn_arch(ia) % forder

        if (term == 1 .and. chi_meth == 1) then
            root => ab_arch(ia,isp) % root
            f1   => ab_arch(ia,isp) % f1
            f2   => ab_arch(ia,isp) % f2
        else if ((term == 2).or.(term == 3)) then
            diag   => ab_arch(ia,isp) % diag 
            eigvec => ab_arch(ia,isp) % eigvec
        end if                          
    end subroutine assoc_ab
    
        
    subroutine free_ab()
        if (allocated(ab_arch)) deallocate(ab_arch)
        if (allocated(cn_arch)) deallocate(cn_arch)        
    end subroutine free_ab
    
    end module ab_io