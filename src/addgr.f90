 
      subroutine addgr(nd,totnd,ad,bins,nbins,rmax,dr)
          use mod_precision


!
!    This procedure updates the g(r) summation.
!

      implicit none

!
!    Declare the simple variables.
!

      real(dp) :: rmax,dr

      integer ind
      integer ia,ib
      integer nd,totnd,nbins

!
!    Declare the arrays.
!

      real(dp) :: dad(3)
      real(dp) :: ad(3,totnd)

      integer bins(nbins)

!
!    Update the histogram.
!

      do ia = 1,nd,1
         do ib = 1,totnd,1
            dad = ad(:,ia) - ad(:,ib)
            ind = int(sum(dad*dad)/dr)
            if ((ind >= 1).and.(ind <= nbins)) bins(ind) = bins(ind) + 1
         enddo
      enddo

      end
