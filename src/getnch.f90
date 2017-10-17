 
      subroutine getnch(ia)
          use mod_precision
          use mod_all_scalar
          use mod_const
          use ab_io

!
!    This is a routine to set up the number of chains with their weights
!     for all the active atoms.
!

      implicit none

      include "Include/Atom.array"

      integer ia,j

!    Set up the number of linear chains associated with each site.
!
      if (momflg == 1) then

!       Averaged Moments.

         nchain = nl(z(ia))
         
         if (nchain > mchain) then
            write(6,'(''Too many chains called for, for atom '', I4)') ia
            write(6,'(''Increase MXNSTAT to at least '',I2)') nchain
            call panic()
         endif
         
         do j = 1,nchain
            wt(j) = real(2*llist(j,z(ia))+1, dp)

         enddo

      elseif ((momflg == 2).or.(momflg == 3)) then

!       No averaging of moments.

         nchain = nstt(z(ia))
         
         if (nchain > mchain) then
            write(6,'(''Too many chains called for, for atom '', & 
     &                I4)') ia
            write(6,'(''Increase MXNSTAT to at least '',I2)') & 
     &                nchain
            call panic()
         endif
        
         wt(1:nchain) = 1.0_dp
      else
         write(6,'(''Illegal value of MOMFLG'')')
         call panic()
      endif

      end

