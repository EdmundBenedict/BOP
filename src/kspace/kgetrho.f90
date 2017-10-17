 
      subroutine kgetrho(elo,ehi)
          use mod_precision


          use mod_all_scalar

          use mod_const

!
!    This is a routine to evaluate the electron density at each active
!     atom, for a given range of energies.
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
      include "../Include/PosVel.array"
      include "../Include/KHamilt.array"


!
!    Declare the simple variables.
!

      real(dp) :: elo,ehi
      real(dp) :: wi, wisum

      integer la,za, nla, nstta
      integer ia, ik, n

      character*80 filename

!
!    Open the output file.
!

      if (idebug == 1) write(6,'(''LENGFN = '',I5)') lengfn

      filename = genfile(1:lengfn)//'.rho'
      open(unit = 99,file = filename,status = 'NEW')

      write(99,'(I6)') nd
      write(99,'(3F12.5)') lena(1),lena(2),lena(3)

      wisum = 0.0_dp
!
!    Evaluate the densities.
!
      do ia = 1, nd, 1

         za = z(ia)
         wi = 0.0_dp
         
         call states(za,nla,nstta,llista)


         do la = 1, nstta, 1

            do ik = 1, nk, 1
               
               do n = 1, khpos(nd+1), 1

                  if ((enk(n,ik) < ehi).and. & 
     &                (enk(n,ik) > elo)) then

                     wi = wi + wtk(ik)* & 
     &                   (real(kpsi(khpos(ia)+la,n,ik), dp)**2 + & 
     &                    aimag(kpsi(khpos(ia)+la,n,ik))**2)

                  endif

               enddo

            enddo

         enddo

         wisum = wisum + wi

         write(99,'(3F9.3,G22.15)') & 
     &         ad(1,ia),ad(2,ia),ad(3,ia),wi         


      enddo

      write(99,'(A,F9.3,A,F9.3,A,F9.3)') 'NO. OF STATES IN RANGE ', & 
     &     elo,' TO ',ehi,' IS ',wisum

!
!    Close the output file.
!

      close(99)



      end

