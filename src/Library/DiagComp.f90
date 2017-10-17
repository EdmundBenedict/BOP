 
 
      subroutine diagcomp(nm,n,ar,ai,w,zr,zi,fv1,fv2,fm1,ierr)
          use mod_precision

!
      integer i,j,n,nm,ierr
      real(dp) :: ar(nm,n),ai(nm,n),w(n),zr(nm,n),zi(nm,n), & 
     &       fv1(n),fv2(n),fm1(2,n)
!
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find the eigenvalues and eigenvectors 
!     of a complex hermitian matrix.
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrix  a=(ar,ai).
!
!        ar  and  ai  contain the real and imaginary parts,
!        respectively, of the complex hermitian matrix.
!
!
!     on output
!
!        w  contains the eigenvalues in ascending order.
!
!        zr  and  zi  contain the real and imaginary parts,
!        respectively, of the eigenvectors.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for tqlrat
!           and tql2.  the normal completion code is zero.
!
!        fv1, fv2, and  fm1  are temporary storage arrays.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      if (n  <=  nm) go to 10
      ierr = 10 * n
      go to 50
!
   10 call  htridi(nm,n,ar,ai,w,fv1,fv2,fm1)
!     .......... find both eigenvalues and eigenvectors ..........
   20 do 40 i = 1, n
         do 30 j = 1, n
            zr(j,i) = 0.0_dp
   30    continue
         zr(i,i) = 1.0_dp
   40 continue

      call  tql2(nm,n,w,fv1,zr,ierr)
      if (ierr  /=  0) go to 50
      call  htribk(nm,n,ar,ai,fm1,n,zr,zi)
   50 return
      end

