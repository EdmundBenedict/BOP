 
   subroutine trncav(ia, nclus)
      use mod_precision
      use mod_all_scalar
      use mod_const
      use ab_io
      use mod_ham
!
!    This is a routine to avaluate the truncator for averaged moments.
!

      implicit none

      integer, intent(in) :: ia, nclus
      
      include "Include/Atom.array"

      
      real(dp) :: taua,taub

      integer i
      integer nmax
      integer nla
      integer la,lb,lla
      integer nstta,nsttb
      integer za
      integer jla,jb

!
!    For each bond do ..
!

      za = z(ia)
      call states(za,nla,nstta,llista)
      jla = 1
      do la = 1,nla
         do lla = 1,2*llista(la)+1
            mlista(jla) = la
            jla = jla + 1
         enddo
      enddo

!
!    For each neighbor ....
!
      do jb = 1, nclus ! start from pcluster(1)==2 to avoid the onsite or 1 to include it
         nsttb = nstt(z(cluster(jb)))
         do la = 1,nstta
            lla = mlista(la)
            nmax = lchain(lla)
            do lb = 1,nsttb
               if (term == 1) then
!old
!                 TAUA = ((LBINF(LLA)/BREC(1,LLA))**2)*DAREC(0,LA,LB,jb)
!                 DO I = 0,NMAX-1,1
!                    TAUA = TAUA - DAREC(I,LA,LB,jb)
!                 ENDDO
!                 DAREC(NMAX,LA,LB,jb) = TAUA
!
!                 TAUB = 0.5D0*(LAINF(LLA)-AREC(0,LLA))*
!     +                  DAREC(0,LA,LB,jb)*
!     +                  LBINF(LLA)/(BREC(1,LLA)**2)
!                 DO I = 1,NMAX,1
!                    TAUB = TAUB - DBREC(I,LA,LB,jb)*LBINF(LLA)/BREC(I,LLA)
!                 ENDDO
!old              DBREC(NMAX+1,LA,LB,jb) = TAUB
                  taua = sum(darec(0:nmax-1,la,lb,jb))
                  
                  darec(nmax,la,lb,jb) = taua*((lbinf(lla)/brec(nmax,lla))**2 - 1.0_dp)
                  dbrec(nmax+1,la,lb,jb) = (0.5_dp*(lainf(lla)-arec(nmax,lla))*taua/brec(nmax,lla)**2)*lbinf(lla)

               elseif ((term == 2).or.(term == 3)) then
                  if (nmax >= nrec) darec(nmax,la,lb,jb) = -sum(darec(0:nmax-1,la,lb,jb))
               endif
            enddo
         enddo
      enddo

      end subroutine trncav

