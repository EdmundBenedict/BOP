 
      subroutine dabdl(ia, zeta, nclus)
          use mod_precision
          use mod_all_scalar
          use mod_const
          use ab_io
          use mod_ham
!
!    This is a routine to evaluate the derivatives of the recursion coefficients.
!

      implicit none
      
      integer, intent(in) :: ia, nclus
      real(dp), intent(in) :: zeta(0:2*mrec+2,mxnstat,mxcls,mxnstat)

      real(dp) :: o(0:mrec,0:2*mrec)

      include "Include/Atom.array"

      integer :: jb
      integer la,lb,n,m
      integer nla,nstta,nlb,nsttb
      integer nmax,ma,lla,lla0,nma
      

!
!    Build O matrix and derivatives of the recursion coefficients.
!
      o = 0.0_dp

      call states(z(ia),nla,nstta,llista)

      if (momflg == 1) then
         lla0 = 0
         do la = 1,nla
            nma = 2*llista(la) + 1
            do ma = 1,nma
               lla = lla0 + ma
               mlista(lla) = la
            enddo
            lla0 = lla0 + nma
         enddo
      else
         do la = 1,nstta
            mlista(la) = la
         enddo
      endif

      do jb = 1, nclus ! start from pcluster(1) == 2 to avoid the onsite (ia==ib) or 1 to include it
         nsttb = nstt(z(cluster(jb)))
         do la = 1,nstta
            ma = mlista(la)
            if (lchain(ma) == nrec) then
               do lb = 1,nsttb
               
                  o(0,0:2*nrec) = zeta(0:2*nrec,lb,jb,la)
                  
!                   print *, 'zeta:', zeta(0:2*nrec,lb,jb,la)
!                   print *, 'brec(1,ma)', brec(1,ma)
!                   print *, 'o(0,1)', o(0,1)
!                   print *, 'brec(1,ma)*o(0,1)', brec(1,ma)*o(0,1)
                  
                  darec(0,la,lb,jb) = brec(1,ma)*o(0,1)
                  dbrec(0,la,lb,jb) = 0.0_dp

                  m = 1
                  do n = m,2*nrec-m
                     o(m,n) = ((arec(n,ma)-arec(m-1,ma))*o(m-1,n) & 
  &                           + brec(n,ma)*o(m-1,n-1) & 
  &                           + brec(n+1,ma)*o(m-1,n+1))/brec(m,ma)
                  enddo
                  darec(m,la,lb,jb) = brec(m+1,ma)*o(m,m+1) - brec(m,ma)*o(m-1,m)
                  dbrec(m,la,lb,jb) = 0.5_dp*brec(m,ma)*(o(m,m)-o(m-1,m-1))

                  do m = 2,nrec
                     do n = m,2*nrec-m
                        o(m,n) = ((arec(n,ma)-arec(m-1,ma))*o(m-1,n) & 
  &                              - brec(m-1,ma)*o(m-2,n) & 
  &                              + brec(n,ma)*o(m-1,n-1) & 
  &                              + brec(n+1,ma)*o(m-1,n+1))/brec(m,ma)
                     enddo
                     
!                      Some elements of o are not set but then used to form d[ab]rec in the next 2 lines. the corresponding d[ab]rec does not seem to be ever used
!                      print *, '++++++++++++++++++++'
!                      print *, 'brec(m,ma)*o(m-1,m)'  , brec(m,ma)*o(m-1,m)
!                      print *, 'brec(m+1,ma)*o(m,m+1)',brec(m+1,ma)*o(m,m+1)
!                      print *, 'o(m,m+1)', o(m,m+1)
!                      print *, 'brec(m+1,ma)', brec(m+1,ma)
                     
                     
!                      print *, '+', brec(m+1,ma)*o(m,m+1) - brec(m,ma)*o(m-1,m)
                     
                     darec(m,la,lb,jb) = brec(m+1,ma)*o(m,m+1) - brec(m,ma)*o(m-1,m)
                     dbrec(m,la,lb,jb) = 0.5_dp*brec(m,ma)*(o(m,m)-o(m-1,m-1))
                  enddo

               enddo

            else

               nmax = lchain(ma)

               do lb = 1,nsttb

                  o(0,0:nmax) = zeta(0:nmax,lb,jb,la)

                  darec(0,la,lb,jb) = brec(1,ma)*o(0,1)
                  dbrec(0,la,lb,jb) = 0.0_dp

                  if (nmax > 0) then
                     m = 1
                     do n = m, nmax
                        o(m,n) = ((arec(n,ma)-arec(m-1,ma))*o(m-1,n) & 
  &                              + brec(n,ma)*o(m-1,n-1) & 
  &                              + brec(n+1,ma)*o(m-1,n+1))/brec(m,ma)
                     enddo
                     darec(m,la,lb,jb) = brec(m+1,ma)*o(m,m+1) - brec(m,ma)*o(m-1,m)
                     dbrec(m,la,lb,jb) = 0.5_dp*brec(m,ma)*(o(m,m)-o(m-1,m-1))
                  endif

                  if (nmax > 1) then
                     do m = 2, nmax
                        do n = m, nmax
                           o(m,n) = ((arec(n,ma)-arec(m-1,ma))*o(m-1,n) & 
  &                                 - brec(m-1,ma)*o(m-2,n) & 
  &                                 + brec(n,ma)*o(m-1,n-1) & 
  &                                 + brec(n+1,ma)*o(m-1,n+1))/brec(m,ma)
                        enddo
                        darec(m,la,lb,jb) = brec(m+1,ma)*o(m,m+1) - brec(m,ma)*o(m-1,m)
                        dbrec(m,la,lb,jb) = 0.5_dp*brec(m,ma)*(o(m,m)-o(m-1,m-1))
                     enddo
                  endif
               enddo
            endif
         enddo
      enddo

      end subroutine dabdl

