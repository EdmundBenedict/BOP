    
    module mod_srt
!
!    Masato's square root terminator integration scheme information.
!
    
!*** The parameters for the polynomial used to evaluate integrals analytically
    
    use mod_precision
    use mod_const
    
    implicit none 
    
    private
    
!     include "Include/ALL.const"
    
!     integer, parameter :: dp = 8
    
    complex(dp), public ::  srtin(0:2*(mrec+3))
    complex(dp), public ::  srtjn(0:2*(mrec+3))
    real(dp), public :: srtrn(0:2*(mrec+3)+2*mfdpol+2)
    complex(dp), public ::  srtvn(0:4*mrec+6)
    complex(dp), public ::  srtvn0(0:4*mrec+6+2*mfdpol+2)
    complex(dp), public ::  srtvn1(0:4*mrec+6+2*mfdpol+2)
    real(dp), public :: srtun(0:2*mfdpol+2)
    real(dp), public :: srtun0(0:2*(mrec+3)+2*mfdpol+2)
    real(dp), public :: srtun1(0:2*(mrec+3)+2*mfdpol+2)

!*** First and second derivatives of polynomial evaluated at roots.
    complex(dp), public, pointer, dimension(:,:) ::  f1=>null(),f2=>null()     !(2*mterm+2,mchain)

!*** Roots.
    complex(dp), public, pointer ::  root(:,:)=>null() !(2*mterm+2,mchain)
    complex(dp), public ::  xp(matrts)

!*** Arrays with information for Masato's integrals
    real(dp), public :: pd(0:2*mrec+6,0:mrec+2)
    real(dp), public :: pd1(0:2*mrec+6,0:mrec+3)
    real(dp), public :: salpha(0:mrec+3)
    real(dp), public :: sbeta(0:mrec+3)
    real(dp), public :: sp(-1:mrec+4,-1:mrec+4)
    real(dp), public :: spa(0:mrec+4)
    real(dp), public :: spb(0:mrec+4)
    real(dp), public :: wn(-1:mrec+2,0:mrec+3)
    real(dp), public :: wn1(-1:mrec+2,0:mrec+3)
    real(dp), public :: sd(0:2*(mrec+3),0:mrec+2)
    real(dp), public :: f(0:2*(mrec+3))
    real(dp), public :: fnag(0:2*(mrec+3))
    real(dp), public :: work1(4*(mrec+3))
    real(dp), public :: np(matrts)

!*** Fermi Dirac polynomial coefficients.
    real(dp), public :: fdpol(0:mfdpol)

!*** Pascal's triangle.
    real(dp), public :: triang(2*mfdpol+3,2*mfdpol+3)

!*** Order of polynomials
    integer,public , pointer :: forder(:)=>null()

! !
! !    Declare the common blocks.
! !

! common /srtint/forder
! common /srtdp/pd,pd1,salpha,sbeta,sp,spa,spb, & 
! &              wn,wn1,sd,f,fnag,work1,srtrn,np, & 
! &              triang,fdpol,srtun0,srtun1,srtun
! common /srtdc/srtin,srtjn,srtvn,f1,f2,root,xp, & 
! &              srtvn0,srtvn1

    end module mod_srt
