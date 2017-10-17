 
      subroutine shiftons(isp, knh, kham,  kpsi)

      use mod_precision
      use mod_all_scalar, only : nd, mag
      use mod_const
      use mod_kspace, only : khpos
      use topologia, only : iproc, mpmap
      use mod_io
      
      implicit none

      integer, intent(in) :: isp, knh
      complex(dp), intent(in) :: kham(knh,knh)
      complex(dp), intent(out) :: kpsi(knh,knh)
      
      
      include "../Include/Atom.array"
      
      
      real(dp) :: dea,phase,rfac,ifac
      
      integer :: i,j,info
      integer :: za,zb,mapb
      integer :: ja,ja0,ia,la,nla,nma,ma,nstta,ip
      
      
      real(dp) :: oe(0:2)
     

      kpsi(:,:) = kham(:,:)
      

      do ia = 1,nd
         dea = de(ia)
         if (mag) dea = dea + oppm(isp)*dem(ia)
         
         call states(z(ia),nla,nstta,llista)
         call onsite(z(ia),oe(0),oe(1),oe(2))
!          print '("ia,nla,llista,oe",5(x,i0),3(x,f6.3))', ia,nla,llista,oe
         j = khpos(ia)
         do la = 1,nla
            nma = 2*llista(la)+1
            do ma = 1, nma
               ip = j + ma
               kpsi(ip, ip) = kpsi(ip,ip) + dea
            end do
            j = j + nma
         end do

      end do
      

   end subroutine shiftons