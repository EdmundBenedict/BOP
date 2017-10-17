 
      subroutine fintrm(ia, nclus)
          use mod_precision
          use mod_all_scalar
          use mod_const
          use ab_io
          use mod_ham, only : cluster

!
!    This is a routine to add in a finite terminator.
!

      implicit none

      integer, intent(in) :: ia, nclus
      
      include "Include/Atom.array"
      include "Include/NebList.array"

      real(dp) :: rulea,ruleb,ra,rb
      real(dp) :: g,phi,w

      integer nma
      integer n
      integer nmax
      integer nla,nlb
      integer ib,la,lb,lla
      integer nstta,nsttb
      integer za
      integer jla,jla0,jb

!
!    Include termination a's and b's, and estimate the corresponding
!    da's and db's from the sum rule.
!

      if (momflg == 1) then

         za = z(ia)
         call states(za,nla,nstta,llista)
         jla = 1
         do la = 1,nla
            do lla = 1,2*llista(la)+1
               mlista(jla) = la
               jla = jla + 1
            enddo
         enddo

         jla0 = 1
         do la = 1,nla
            nma = 2*llista(la)+1
            nmax = lchain(la)

            if ((nmax == nrec).and.(nmax < mrec)) then

               if (nmax >= 1) then
                  call getptm(arec(nmax-1,la),arec(nmax,la), & 
     &                        brec(nmax,la),lainf(la), & 
     &                        lbinf(la),g,phi,w)
                  do n = nmax+1,mrec
                     call evlptm(arec(n,la),brec(n,la), & 
     &                           n-nmax+1,lainf(la), & 
     &                           lbinf(la),g,phi,w)
                  enddo
               else
                  do n = nmax+1,mrec
                     arec(n,la) = lainf(la)
                     brec(n,la) = lbinf(la)
                  enddo
               endif

               do jb = 1, nclus ! start from pcluster(1) == 2 to avoid the onsite (ia==ib) or 1 to include it
                  nsttb = nstt(z(cluster(jb)))
                  jla = jla0
                  do lla = 1,nma
                     do lb = 1,nsttb

                        rulea = brec(nrec+1,la)**2 - brec(nrec,la)**2
                        ra = rulea/(brec(1,la)**2)
                        darec(nrec,jla,lb,jb) = ra*darec(0,jla,lb,jb)

                        do n = nrec+1,mrec

                           rulea = brec(n+1,la)**2 - brec(n,la)**2
                           ra = rulea/(brec(1,la)**2)
                           darec(n,jla,lb,jb) = ra*darec(0,jla,lb,jb)

                           ruleb = 0.5_dp*brec(n,la)*(arec(n,la)-arec(n-1,la))
                           rb = ruleb/(brec(1,la)**2)
                           dbrec(n,jla,lb,jb) = rb*darec(0,jla,lb,jb)
                        enddo
                     enddo
                     jla = jla + 1
                  enddo
               enddo
               lchain(la) = mrec
            endif
            jla0 = jla0 + nma
         enddo
      elseif (momflg == 2) then
         za = z(ia)
         call states(za,nla,nstta,llista)

         do la = 1,nstta
            nmax = lchain(la)
            if ((nmax == nrec).and.(nmax < mrec)) then
               if (nmax >= 1) then
                  call getptm(arec(nmax-1,la),arec(nmax,la), & 
     &                        brec(nmax,la),lainf(la), & 
     &                        lbinf(la),g,phi,w)
                  do n = nmax+1,mrec,1
                     call evlptm(arec(n,la),brec(n,la), & 
     &                           n-nmax+1,lainf(la), & 
     &                           lbinf(la),g,phi,w)
                  enddo
               else
                  do n = nmax+1,mrec,1
                     arec(n,la) = lainf(la)
                     brec(n,la) = lbinf(la)
                  enddo
               endif

               do jb = 1, nclus ! start from pcluster(0) to include ib == ia 
                  ib = cluster(jb)
                  nsttb = nstt(z(ib))
                  do lb = 1,nsttb
                     if ((ia /= ib).or.(la /= lb)) then
                        rulea = brec(nrec+1,la)**2 - brec(nrec,la)**2
                        ra = rulea/(brec(1,la)**2)
                        darec(nrec,la,lb,jb) = ra*darec(0,la,lb,jb)

                        do n = nrec+1,mrec
                           rulea = brec(n+1,la)**2 - brec(n,la)**2
                           ra = rulea/(brec(1,la)**2)
                           darec(n,la,lb,jb) = ra*darec(0,la,lb,jb)

                           ruleb = 0.5_dp*brec(n,la)*(arec(n,la)-arec(n-1,la))
                           rb = ruleb/(brec(1,la)**2)
                           dbrec(n,la,lb,jb) = rb*darec(0,la,lb,jb)
                        enddo
                     endif
                  enddo
               enddo
               lchain(la) = mrec
            endif
         enddo
      elseif (momflg == 3) then
         za = z(ia)
         call states(za,nla,nstta,llista)

         do la = 1,nstta

            nmax = lchain(la)

            if ((nmax == nrec).and.(nmax < mrec)) then

               if (nmax >= 1) then
                  call getptm(arec(nmax-1,la),arec(nmax,la), & 
     &                        brec(nmax,la),lainf(la), & 
     &                        lbinf(la),g,phi,w)
                  do n = nmax+1,mrec
                     call evlptm(arec(n,la),brec(n,la), & 
     &                           n-nmax+1,lainf(la), & 
     &                           lbinf(la),g,phi,w)
                  enddo
               else
                  do n = nmax+1,mrec
                     arec(n,la) = lainf(la)
                     brec(n,la) = lbinf(la)
                  enddo
               endif

               do jb = 1, nclus ! start from pcluster(1) == 2 to avoid the onsite (ia==ib) or 1 to include
                  nsttb = nstt(z(cluster(jb)))
                  do lb = 1,nsttb
                     rulea = brec(nrec+1,la)**2 - brec(nrec,la)**2
                     ra = rulea/(brec(1,la)**2)
                     darec(nrec,la,lb,jb) = ra*darec(0,la,lb,jb)
                     do n = nrec+1,mrec,1
                        rulea = brec(n+1,la)**2 - brec(n,la)**2
                        ra = rulea/(brec(1,la)**2)
                        darec(n,la,lb,jb) = ra*darec(0,la,lb,jb)
                        
                        ruleb = 0.5_dp*brec(n,la)*(arec(n,la)-arec(n-1,la))
                        rb = ruleb/(brec(1,la)**2)
                        dbrec(n,la,lb,jb) = rb*darec(0,la,lb,jb)
                     enddo
                  enddo
               enddo
               lchain(la) = mrec
            endif
         enddo
      endif
      
      end subroutine fintrm

