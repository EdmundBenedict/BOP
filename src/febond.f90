 
      function febond(ia,ef)
          use mod_precision
          use mod_all_scalar
          use mod_const
          use ab_io
          use mod_chi
          use mod_ham

!
!    This is a subroutine to evaluate the bond contribution to the
!     band structure energy
!

      implicit none


      include "Include/Atom.array"

      integer, intent(in) :: ia
      real(dp), intent(in) :: ef
      real(dp) :: febond

      real(dp) :: hbo
      
      integer nmax
      integer nla
      integer la,lb
      integer nstta,nsttb
      integer za
      integer jb
      integer lla0,nma,ma

      febond = 0.0_dp

      za = z(ia)
      call states(za,nla,nstta,llista)

      if (momflg == 1) then
         lla0 = 0
         do la = 1,nla
            nma = 2*llista(la) + 1
            mlista(lla0+1:lla0+nma) = la
            lla0 = lla0 + nma
         enddo
      else
         do la = 1, nstta
            mlista(la) = la
         end do
      endif

!    For each neighbor ....
      do jb = pcluster(1), pcluster(2)-1 ! starts from pcls(1) to avoid onsite (ia==ib)
        
         nsttb = nstt(z(cluster(jb)))
         do la = 1,nstta
            ma = mlista(la)
            nmax = lchain(ma)
            do lb = 1, nsttb
               hbo =  sum(chia(0:nmax,ma) * darec(0:nmax,la,lb,jb)) &
                & + 2*sum(chib(1:nmax,ma) * dbrec(1:nmax,la,lb,jb)) 
!                   bo(la,lb) =   dot_product(chia(0:nmax,ma), darec(0:nmax,la,lb,jb)) &
!                           & + 2*dot_product(chib(1:nmax,ma), dbrec(1:nmax,la,lb,jb)) 

               if (term == 1) hbo = hbo + 2*chib(nmax+1,ma) * dbrec(nmax+1,la,lb,jb)
               
               hbo = -2*hbo
               febond = febond + hbo*darec(0,la,lb,jb)

!                print *, 'ia, la, ib, lb, dar, h',  ia, la, cluster(jb), lb, darec(0,la,lb,jb), h(la,lb,jb,1)          
            enddo
         enddo
         
        
      enddo
      
      end function febond
      
      
      
      
