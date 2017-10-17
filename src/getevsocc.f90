 
      subroutine getevsocc()
          use mod_precision
          use mod_all_scalar
          use mod_const
          use ab_io
          use mod_chi

!
!    Calculate the energy vs occupancy curves.
!
! This is not up to date, do not use before fixing it!

      implicit none


      include "Include/Atom.array"
      include "Include/Misc.array"

      real(dp) :: getefsrt,getefnull

      integer i,flag,wrflag
      integer ia,la

      character*80 filename

!     
!     Find the occupancies and fermi energies.
!     
      flag = 0
      call getocc(nfermi,multocc)
      if (term == 1) then
         do i = 1,nfermi,1
            multef(i) = getefsrt(multocc(i),flag)
         enddo
      elseif ((term == 2).or.(term == 3)) then
         do i = 1,nfermi,1
            multef(i) = getefnull(multocc(i),flag)
         enddo
      endif
      
!
!     Open the output files.
!
      
      filename = genfile(1:lengfn)//'.eprom'
      open(unit = 1,file = filename,status = 'NEW')
      
      filename = genfile(1:lengfn)//'.ebond'
      open(unit = 2,file = filename,status = 'NEW')
      
      filename = genfile(1:lengfn)//'.eband'
      open(unit = 3,file = filename,status = 'NEW')
      
!     
!     Evaluate and write out the energies.
!     
      
      do i = 1,nfermi,1

         eprom = 0.0_dp
         ebond = 0.0_dp

!     
!     Calculate the susceptibilities.
!     

         wrflag = 1
         do ia = 1,nd,1
            call rdab(wrflag,ia)
            do la = 1,nchain,1
               if (term == 1) then
                  if (chi_meth == 1) then
                     call getchisrt(multef(i), lchain(la)+1,la)
                  else
                     call getchisrt_c(multef(i), lchain(la)+1,la)
                  endif
               elseif ((term == 2).or.(term == 3)) then
                  call getchinull(multef(i),lchain(la),mrec,chia(0,la),chib(0,la),diag(1,la),eigvec(1,1,la),kt)
               endif
            enddo
         
!     
!          The energy.
!     
         
            if (momflg == 1) then
            
!     
!          Averaged moments.
               call epromavg(ia,multef(i))
               call febond(ia)
            
            elseif (momflg == 2) then
 
!     
!          Mixed basis.
!     
 
               call eprommix(ia,multef(i))
               call febond(ia)

            elseif (momflg == 3) then
 
!     
!          Non-averaged basis.
!     
 
               call epromnoavg(ia,multef(i))
               call febond(ia)
 
            endif
         
         enddo

         write(1,'(2G23.15)') multocc(i)/real(nd, dp),eprom/real(nd, dp)
         write(2,'(2G23.15)') multocc(i)/real(nd, dp),ebond/real(nd, dp)
         write(3,'(2G23.15)') multocc(i)/real(nd, dp), & 
     &        (eprom+ebond)/real(nd, dp)

         close(50)
 
      enddo
      
!     
!     Close the output files.
!     
      
      close(1)
      close(2)
      close(3)
      
      end
