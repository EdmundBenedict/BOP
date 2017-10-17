
   subroutine trncno(ia, nclus)
   
   use mod_precision
   use mod_all_scalar
   use mod_const
   use ab_io
   use mod_ham
   
   !
   !    This is a routine to evaluate the truncator for
   !     no averaging of moments.
   !

   implicit none

   integer, intent(in) :: ia, nclus

   include "Include/Atom.array"

   real(dp) :: taua,taub

   integer i,nmax
   integer nla,nlb
   integer ib,la,lb
   integer nstta,nsttb
   integer za,zb
   integer ja,ja0,jb

!
!    For each bond do ..
!

   nstta = nstt(z(ia))
!
!    For each neighbor ....
!


   do jb = 1, nclus
      nsttb = nstt(z(cluster(jb)))
      do la = 1,nstta
         nmax = lchain(la)
         do lb = 1,nsttb

            if (term == 1) then

               taua = ((lbinf(la)/brec(1,la))**2)* darec(0,la,lb,jb)
               do i = 0,nmax-1
                  taua = taua - darec(i,la,lb,jb)
               enddo
               darec(nmax,la,lb,jb) = taua

               taub = 0.5_dp*(lainf(la)-arec(0,la))* & 
   &                      darec(0,la,lb,jb)* & 
   &                      lbinf(la)/(brec(1,la)**2)
               do i = 1,nmax
                  taub = taub - dbrec(i,la,lb,jb)* & 
   &                                 lbinf(la)/brec(i,la)
               enddo
               dbrec(nmax+1,la,lb,jb) = taub

            elseif ((term == 2).or.(term == 3)) then

               if (nmax >= nrec) then

                  taua = 0.0_dp
                  do i = 0,nmax-1
                     taua = taua - darec(i,la,lb,jb)
                  enddo
                  darec(nmax,la,lb,jb) = taua
               endif
            endif
         enddo
      enddo
   enddo

   end subroutine trncno

