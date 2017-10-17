 
      subroutine fcfit()
          use mod_precision


          use mod_all_scalar

          use mod_const
!
!     This subroutine should help with the fitting of force constants
!     and help to show the differences between the force constants from
!     BOP and values calculated ab initio
!
!     M. Cawkwell. 15th April 2004
!
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
!    Include arrays.
!
      include "Include/Force.array"
      include "Include/ag.conn"
!
      integer i, j
      real(dp) :: abinitf(3,8)
      character*40 flenm
!
      write(flenm,'("fcfit.A",F3.1,".mu",F4.2,".dat")') aescra, aescrmu
      open (unit=38, status="UNKNOWN", file=flenm)
!
!     Now we list the forces calculated ab initio (make sure you get
!     the atoms in the right order in the list)
!
!     In this case we're doing it for bcc Mo:
!
      do i = 1,3
         do j = 1,8
            abinitf(i,j) = 0.0_dp
         enddo
      enddo
!
      abinitf(1,1) = -3.882e-3_dp
      abinitf(1,2) = 2.026e-2_dp
      abinitf(2,2) = 3.414e-2_dp
      abinitf(3,2) = abinitf(2,2)
      abinitf(1,3) = 5.194e-2_dp
      abinitf(1,4) = 3.093e-2_dp
      abinitf(2,4) = -4.055e-2_dp
      abinitf(3,4) = abinitf(2,4)
      abinitf(1,5) = 4.078e-2_dp
      abinitf(1,6) = -1.543e-2_dp
      abinitf(1,7) = 3.051e-1_dp
      abinitf(1,8) = -6.087e-1_dp
!
      write(38,'("FORCES FROM SBOP")')
      write(38,'("__________________________________________________")')
      write(38,10) (ftot(1,i), ftot(2,i), ftot(3,i), i = 1,8)
      write(38,'("__________________________________________________")')
!
      write(38,'("DIFF. BETWEEN BOP AND AB INITIO FORCES: bcc Mo")')
      write(38,'("__________________________________________________")')
!
      write(38,10) (ftot(1,i)-abinitf(1,i), ftot(2,i)-abinitf(2,i), & 
     &     ftot(3,i)-abinitf(3,i), i=1,8)
!
      write(38,'("__________________________________________________")')      
      write(38,'("RELATIVE ERROR IN FORCES: (FTOT - ABINITF)/ABINITF")')
      write(38,'("__________________________________________________")')
      write(38,11) "ORIGIN (1/2, 1/2, 1/2) XX",  & 
     &     (ftot(1,8)-abinitf(1,8))/abinitf(1,8)
      write(38,11) "1ST NN (1/4, 1/4, 1/4) XX", & 
     &     (ftot(1,2)-abinitf(1,2))/abinitf(1,2)
      write(38,11) "1ST NN (1/4, 1/4, 1/4) XY", & 
     &     (ftot(2,2)-abinitf(2,2))/abinitf(2,2)
      write(38,11) "1ST NN (3/4, 1/4, 1/4) XX", & 
     &     (ftot(1,4)-abinitf(1,4))/abinitf(1,4)
      write(38,11) "1ST NN (3/4, 1/4, 1/4) XY", & 
     &     (ftot(2,4)-abinitf(2,4))/abinitf(2,4)
      write(38,11) "2ND NN (O, 1/2, 1/2) XX  ", & 
     &     (ftot(1,7)-abinitf(1,7))/abinitf(1,7)
      write(38,11) "2ND NN (1/2, 1/2, 0) XX  ", & 
     &     (ftot(1,6)-abinitf(1,6))/abinitf(1,6)
      write(38,11) "3RD NN (1/2, 0, 0) XX    ", & 
     &     (ftot(1,3)-abinitf(1,3))/abinitf(1,3)
      write(38,11) "3RD NN (0, 1/2, 0) XX    ", & 
     &     (ftot(1,5)-abinitf(1,5))/abinitf(1,5)
      write(38,11) "5TH NN (0, 0, 0) XX      ", & 
     &     (ftot(1,1)-abinitf(1,1))/abinitf(1,1)
!
 10   format(1x,f6.4,1x,f6.4,1x,f6.4)
 11   format(1x,a25,15x,f9.5)
!
      close(38)
!
      return
!
      end
