 
      subroutine locsym(h,transm,work,eval,mxnstat,nstat)
          use mod_precision


!
!    This is a routine to find starting states for the recursion
!     which refelect the local symmetry.
!

      implicit none

!
!    Declare the simple variables.
!

      integer i,j,ifail
      integer nstat
      integer mxnstat

!
!    Declare the arrays.
!

!*** Hamiltonian matrix.
      real(dp) :: h(mxnstat,mxnstat)
!*** Work space
      real(dp) :: work(mxnstat,mxnstat)
!*** Transformation matrix. Indices: Basis state index, Transformed state index, Atom index.
      real(dp) :: transm(mxnstat,mxnstat)
!*** Eigenvalues
      real(dp) :: eval(mxnstat)

!
!    Force nearly zero components of the hamiltonian type matrix to zero.
!

      do j = 1,nstat,1
         do i = 1,nstat,1
            if (abs(h(i,j)) < 1.0e-6_dp) h(i,j)=0.0_dp
         enddo
      enddo

!
!    Diagonalize the Hamiltonian.
!

!     CALL F02ABF(H,MXNSTAT,NSTAT,EVAL,TRANSM,MXNSTAT,WORK,IFAIL)
      call diagsym(mxnstat,nstat,h,eval,transm,work,ifail)

      end

