
   function keprom()
      use mod_precision
      use mod_all_scalar, only : mag, nd, nsp
      use mod_const
      use mod_kspace, only : khpos, krhod

!    This is a routine to evaluate the promotion energy.

      implicit none

      integer :: isp
      
      include "../Include/Atom.array"

      real(dp) :: keprom
      

      real(dp) :: el

      integer nstta,nla,ll
      integer ia,la,ma,i,ip,ik,ierr

      real(dp) :: oe(0:2), epromes(0:2), deia


      epromes = 0.0_dp
      
      do isp = 1, nsp
         do ia = 1, nd
            deia = de(ia)
            if (mag) deia = deia + oppm(isp)*dem(ia)
         
            oe(0) = es(z(ia)) + deia
            oe(1) = ep(z(ia)) + deia
            oe(2) = ed(z(ia)) + deia

            ip = khpos(ia) + 1
            do la = 1,nl(z(ia))
               ll = llist(la,z(ia))
               ma = 2*ll
               epromes(ll) = epromes(ll) + oe(ll)*sum(krhod(ip:ip+ma,isp))
               ip = ip + ma + 1
            enddo
         enddo
      end do
      
!       if (.not.mag) epromes = 2*epromes
      
!       eproms = epromes(0)
!       epromp = epromes(1)
!       epromd = epromes(2)

      keprom = sum(epromes)
      if (.not. mag) keprom = 2*keprom

   end function keprom

