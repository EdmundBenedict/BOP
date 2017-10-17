 
       subroutine rhoh()
          use mod_precision
          use mod_all_scalar
          use mod_const
          use mod_kspace


      implicit none

      include "../Include/Atom.array"

      
      real(dp) :: el

      integer nstta,nla
      integer ia,la,ma,i,ip,ik

      logical usepot
      real(dp) :: oe(0:2), epromes(0:2)




      do ik = 1,nk

         do ih = 





      do ik = 1,nk


         i = 1
         do while (i <= knh)


            do ia = 1,nd,1

               call states(z(ia),nla,nstta,llista)
               call onsite(z(ia),oe(0),oe(1),oe(2))

               ip = khpos(ia)
               do la = 1,nla,1
                  el = oe(llista(la))
                  
                  do ma = 1,2*llista(la)+1
                     ip = ip + 1                     
                     epromes(llista(la)) = epromes(llista(la)) &
     &                         + occ(i,ik) * (real(kpsi(ip,i,ik))**2 & 
     &                                     + aimag(kpsi(ip,i,ik))**2) * el
                  enddo
               enddo
            enddo
            i = i + 1
         enddo

!          endif

      enddo

!      IF (USEPOT(P_KSPACE)) CALL GANDSDP(704,SUM)

      eproms = epromes(0)
      epromp = epromes(1)
      epromd = epromes(2)

      keprom = sum(epromes)

      end subroutine rhoh

