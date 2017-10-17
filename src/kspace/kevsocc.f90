 
      subroutine kevsocc()
          use mod_precision


          use mod_all_scalar

          use mod_const

!
!    Calculate the energy vs occupancy curves.
!

      implicit none

!
!    Include constants.
!

!      include "../Include/ALL.const"

!
!    Include scalars.
!

!      include "../Include/ALL.scalar"

!
!    Include arrays.
!

      include "../Include/Atom.array"
      include "../Include/KHamilt.array"
!      include "../Include/Misc.array"
!      include "../Include/NebList.array"

!
!    Declare the simple variables.
!

      real(dp) :: sum,nocc,ef,bseng,keprom,etotal

      integer n

      character*80 filename

!
!     Open the output files.
!
      
      filename = genfile(1:lengfn)//'.eprom'
      open(unit = 1,file = filename,status = 'NEW')
      
      filename = genfile(1:lengfn)//'.ebond'
      open(unit = 2,file = filename,status = 'NEW')
      
      filename = genfile(1:lengfn)//'.eband'
      open(unit = 3,file = filename,status = 'NEW')
 
      write(1,'(2G23.15)') 0.0_dp,0.0_dp
      write(2,'(2G23.15)') 0.0_dp,0.0_dp
      write(3,'(2G23.15)') 0.0_dp,0.0_dp

!     
!     Evaluate and write out the energies.
!     

!      ETOT = ETOTAL(KPSI,1)

      do n = 1,2*kmxh

         nocc = real(n, dp)
         call fndocc(nk,mxnk,nd,mxnd,kmxh,mxnstat,enk, & 
     &               wtk,nocc,occ,ef,kt,sum)
         eband = bseng(kmxh,enk,wtk,occ)
         eprom = keprom(kpsi,occ,khpos, & 
     &                  wtk,z,llista)
         ebond = eband - eprom

         write(1,'(2G23.15)') nocc/real(nd, dp),eprom/real(nd, dp)
         write(2,'(2G23.15)') nocc/real(nd, dp),ebond/real(nd, dp)
         write(3,'(2G23.15)') nocc/real(nd, dp),eband/real(nd, dp)

      enddo
      
!     
!     Close the output files.
!     
      
      close(1)
      close(2)
      close(3)
      
      end
