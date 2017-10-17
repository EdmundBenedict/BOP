 
      subroutine precab(ia)
          use mod_precision
          use mod_all_scalar
          use mod_const
          use ab_io
          use mod_ham
          use mod_clock
!
!    This is a routine to build the parallel recursion coefficients,
!     for averaged moments.
!

      implicit none

      include "Include/Atom.array"
      include "Include/Moment.array"
      include "Include/NebList.array"


      real(dp) :: dot,lng,inv_lng,expfac,fla

      integer n0,n1,n2,n3,t1,t2,t3,tmp
      integer ja,ja0,jb
      integer ncross,m
      integer lla,ma,lla0,nma,la,lb,n,ia,ib,i
      integer nla,nstta,nlb,nsttb
      integer findj
      integer dj
      character(len=100) :: precsform,pform1
      real(dp) :: t2lng
      integer(8) :: c1,c2,nc1,nc2,nc
      integer :: cnt

      real(dp),parameter :: zero = 1.0e-3_dp

      cnt = 0
      call states(z(ia),nla,nstta,llista)
      ja = findj(bptr(aptr(ia)),mxnb,ia)
      lla0 = 0
      do la = 1,nla,1
         brec(0,la) = 0.0_dp
         lchain(la) = nrec
         nma = 2*llista(la)+1
         fla = sqrt(real(nma, dp))
         do lla = lla0+1,lla0+nma,1
            ! lla = lla0 + ma
            do lb = 1,mxnstat,1
               chebt(lb,1,lla,1) = 0.0_dp
               chebt(lb,1,lla,2) = 0.0_dp
               chebt(lb,1,lla,3) = 0.0_dp
               do jb = 1,mxnnb,1
                  zeta(0,lb,jb,lla) = 0.0_dp
               enddo
            enddo
            chebt(lla,1,lla,2) = 1.0_dp/fla
            zeta(0,lla,ja,lla) = 1.0_dp
         enddo
         
         lla0 = lla0 + nma
      enddo

      decipher(ia) = 1

!
!    Part 1: Recursion with growing cluster.
!

      t1 = 1
      t2 = 2
      t3 = 3
      if (ia == 1) call system_clock(c1)
      nc = 0
      do n = 1,nbase
      print *,'np1:',n
         n3 = pcluster(n-1)-1
         n2 = pcluster(n)-1
         n1 = pcluster(n+1)-1
        
!          print *,(/n1,n2,n3/)*mxnstat

         do i = n2+1,n1
            decipher(cluster(i)) = i
         enddo
!          if (ia == 1) print '("prec arec:",6(x,f16.12))',arec(:5,1)
!          if (ia == 1) print '("prec brec:",5(x,f16.12))',brec(1:5,1)
         lla0 = 0
         do la = 1,nla
!          print *, chebt(1:mxnstat,1,la,t2)
            nma = 2*llista(la)+1
! print *,'faza 1', nla, lchain(la)
            if (lchain(la) == nrec) then
! print *, 'faza 1'
               arec(n-1,la) = 0.0_dp
               t2lng = 0.0_dp
               do lla = lla0+1,lla0+nma,1
!                t2lng = t2lng + dot(chebt(1,1,lla,t2),chebt(1,1,lla,t2),n2*mxnstat)
!                   lla = lla0 + ma
!                   write(pform1,'(a,i0,a)') '(64x,x,f16.12,"|",',n2*mxnstat,'(x,f16.12))'
!                   print trim(pform1), sqrt(dot(chebt(1,1,lla,t2),chebt(1,1,lla,t2),n2*mxnstat)), chebt(1:mxnstat,1:n2,lla,t2)
                  if (ia == 1) call system_clock(nc1)
!                   call gethu(chebt(1,1,lla,t2),n2,chebt(1,1,lla,t1),n1)
                  call gethu(chebt(1,1,lla,t2),n2,hup(1,1,lla,n-1),n1))
                  if (ia == 1) then
                     call system_clock(nc2)
                     nc = nc + nc2 - nc1
                     cnt = cnt + 1
                  end if
!                if (ia == 1 .and. la==1) then
!                   write(precsform,'(a,i0,a)') '("prec recs:",2(x,i0),3(x,f16.12),"|",',n2*mxnstat,'(x,f16.12))'
                  
!                   print trim(precsform),n-1,la,arec(n-1,la),brec(n-1,la),sqrt(dot(chebt(1,1,lla,t1),chebt(1,1,lla,t1),n2*mxnstat)), chebt(1:mxnstat,1:n2,lla,t1)
                  
!                endif

                  arec(n-1,la) = arec(n-1,la) + dot(chebt(1,1,lla,t2), chebt(1,1,lla,t1), n2*mxnstat)
                  
               enddo
!                print *,redb,'t2lng:',t2lng,endc
               lng = 0.0_dp
               do lla = lla0+1,lla0+nma,1
                  call daxpy(n2*mxnstat, -arec(n-1,la), chebt(1,1,lla,t2), 1, chebt(1,1,lla,t1), 1 )
                  call daxpy(n3*mxnstat, -brec(n-1,la), chebt(1,1,lla,t3), 1, chebt(1,1,lla,t1), 1 )
                  lng = lng + dot(chebt(1,1,lla,t1), chebt(1,1,lla,t1),n1*mxnstat)
               end do
               


!                lng = 0.0d0
!                print *,nma
!                do lla = lla0+1,lla0+nma,1
!                   lla = lla0 + ma
!                   lng = lng + dot(chebt(1,1,lla,t1), & 
!      &                            chebt(1,1,lla,t1),n1*mxnstat)
!                enddo
               lng = sqrt(lng)
               inv_lng = 1.0_dp/lng
               brec(n,la) = lng

               if (brec(n,la) <= zero) then
                  lchain(la) = n-1
!                   m = n
                  do m = n,2*nrec,1
!                   do while (m.le.2*nrec)
                     arec(m,la) = 0.0_dp
                     brec(m,la) = 0.0_dp
                     ja = aptr(ia)
                     ja0 = ja
                     ib = bptr(ja)
                     do while (ib /= eol)
                        dj = ja-ja0+1
                        call states(z(ib),nlb,nsttb,llistb)
                        do lla = lla0+1,lla0+nma,1
!                            lla = lla0 + ma
                           do lb = 1,nsttb,1
                              zeta(m,lb,dj,lla) = 0.0_dp
                           enddo
                        enddo
                        ja = ja + 1
                        ib = bptr(ja)
                     enddo
!                      m = m + 1
                  enddo
               else
! print *, 'faza 1 else'
                  do lla = lla0+1,lla0+nma,1
!                      lla = lla0 + ma
                     do i = 1,n1,1
                        do lb = 1,nstt(z(cluster(i))),1
!                            chebt(lb,i,lla,t1) = chebt(lb,i,lla,t1)/lng
                           chebt(lb,i,lla,t1) = inv_lng*chebt(lb,i,lla,t1)
                        enddo
                     enddo
                     
                     
                  enddo
               endif

            endif

            lla0 = lla0 + nma

         enddo

         ja = aptr(ia)
         ja0 = ja
         ib = bptr(ja)
         do while (ib /= eol)
            dj = ja-ja0+1
            jb = decipher(ib)
            call states(z(ib),nlb,nsttb,llistb)
            lla0 = 0
            do la = 1,nla,1
               nma = 2*llista(la)+1
               fla = sqrt(real(nma, dp))
               if (lchain(la) == nrec) then
                  do lla = lla0+1,lla0+nma,1
!                      lla = lla0 + ma
                     do lb = 1,nsttb,1
                        zeta(n,lb,dj,lla) = chebt(lb,jb,lla,t1)*fla
                     enddo
                  enddo
               endif
               lla0 = lla0 + nma
            enddo
            ja = ja + 1
            ib = bptr(ja)
         enddo
 
         tmp = t3
         t3 = t2
         t2 = t1
         t1 = tmp

      enddo
      if (ia == 1) then 
         call system_clock(c2)
         print '("part1:",3(x,f16.12),x,i0)', real(c2-c1,dp)/real(cr,dp), real(nc,dp)/real(cr,dp), real(c2-c1-nc,dp)/real(cr,dp),cnt
      end if


!
!   Part 2: Recursion at full cluster size.
!

      ncross = 2*nrec-nbase+1

      n3 = pcluster(nbase)-1
      n2 = pcluster(nbase+1)-1
      n1 = pcluster(nbase+1)-1
      if (ia == 1) call system_clock(c1)
      nc = 0
      do n = nbase+1,nrec+1,1
      print *,'np2:',n
         lla0 = 0
         do la = 1,nla,1

            nma = 2*llista(la)+1

            if (lchain(la) == nrec) then
! print *, 'faza 2'
               arec(n-1,la) = 0.0_dp
               do lla = lla0+1,lla0+nma,1
!                   lla = lla0 + ma
                  if (ia == 1) call system_clock(nc1)
                  call gethu(chebt(1,1,lla,t2),n2,chebt(1,1,lla,t1),n1)
                  if (ia == 1) then
                     call system_clock(nc2)
                     nc = nc + nc2 - nc1
                     cnt = cnt + 1
                  end if
                  arec(n-1,la) = arec(n-1,la) + dot(chebt(1,1,lla,t2), chebt(1,1,lla,t1),n2*mxnstat)
               enddo

               do lla = lla0+1,lla0+nma,1
!                   lla = lla0 + ma
                  
                  if (n2*mxnstat < 100) then
                    do i = 1,n2,1
                        do lb = 1,nstt(z(cluster(i))),1
                            chebt(lb,i,lla,t1) = chebt(lb,i,lla,t1) - arec(n-1,la)*chebt(lb,i,lla,t2)
                        enddo
                    enddo
                  else  
                    call daxpy(n2*mxnstat, -arec(n-1,la), chebt(1,1,lla,t2), 1, chebt(1,1,lla,t1), 1 )                  
                  endif
                  
                  if (n3*mxnstat < 100) then
                    do i = 1,n3,1
                        do lb = 1,nstt(z(cluster(i))),1
                            chebt(lb,i,lla,t1) = chebt(lb,i,lla,t1) - brec(n-1,la)*chebt(lb,i,lla,t3)
                        enddo
                    enddo
                  else
                    call daxpy(n3*mxnstat, -brec(n-1,la), chebt(1,1,lla,t3), 1, chebt(1,1,lla,t1), 1 )                  
                  endif
               enddo

               lng = 0.0_dp
               do lla = lla0+1,lla0+nma,1
!                   lla = lla0 + ma
                  lng = lng + dot(chebt(1,1,lla,t1), & 
     &                            chebt(1,1,lla,t1),n1*mxnstat)
               enddo
               lng = sqrt(lng)
               brec(n,la) = lng

               if (brec(n,la) <= zero) then
                  lchain(la) = n-1
!                   m = n
                  do m=n,2*nrec,1
                     arec(m,la) = 0.0_dp
                     brec(m,la) = 0.0_dp
                     ja = aptr(ia)
                     ja0 = ja
                     ib = bptr(ja)
                     do while (ib /= eol)
                        dj = ja-ja0+1
                        call states(z(ib),nlb,nsttb,llistb)
                        do lla = lla0+1,lla0+nma,1
!                            lla = lla0 + ma
                           do lb = 1,nsttb,1
                              zeta(m,lb,dj,lla) = 0.0_dp
                           enddo
                        enddo
                        ja = ja + 1
                        ib = bptr(ja)
                     enddo
!                      m = m + 1
                  enddo
               else
! print *, 'faza 2 else'
                  do lla = lla0+1,lla0+nma,1
                     ! lla = lla0 + ma
                     do i = 1,n1,1
                        do lb = 1,nstt(z(cluster(i))),1
                           chebt(lb,i,lla,t1) = chebt(lb,i,lla,t1)/lng
                        enddo
                     enddo
                  enddo
               endif

            endif

            lla0 = lla0 + nma

         enddo

         ja = aptr(ia)
         ja0 = ja
         ib = bptr(ja)
         do while (ib /= eol)
            
            dj = ja-ja0+1
            jb = decipher(ib)
            call states(z(ib),nlb,nsttb,llistb)
            lla0 = 0
            do la = 1,nla,1
               nma = 2*llista(la)+1
               fla = sqrt(real(nma, dp))
               if (brec(nrec+1,la) > zero) then
                  do lla = lla0+1,lla0+nma,1
                     ! lla = lla0 + ma
                     do lb = 1,nsttb,1
                        zeta(n,lb,dj,lla) = chebt(lb,jb,lla,t1)*fla
                     enddo
                  enddo
               endif
               lla0 = lla0 + nma
            enddo
            ja = ja + 1
            ib = bptr(ja)
         enddo
 
         tmp = t3
         t3 = t2
         t2 = t1
         t1 = tmp

         n3 = pcluster(nbase+1)-1

      enddo
      if (ia == 1) then 
         call system_clock(c2)
         print '("part2:",3(x,f16.12),x,i0)', real(c2-c1,dp)/real(cr,dp), real(nc,dp)/real(cr,dp), real(c2-c1-nc,dp)/real(cr,dp),cnt
      end if

!
!   Part 3: Extension at constant cluster size.
!

      expfac = 1.0_dp
      if (ia == 1) call system_clock(c1)
      nc=0
      do n = nrec+2,ncross,1
      print *,'np3:',n
         lla0 = 0
         do la = 1,nla,1

            nma = 2*llista(la)+1

            if (brec(nrec+1,la) > zero) then
! print *, 'faza 3'
               arec(n-1,la) = arec(nrec,la)
               brec(n,la) = brec(nrec+1,la)*expfac

               do lla = lla0+1,lla0+nma,1
                  ! lla = lla0 + ma
                  if (ia == 1) call system_clock(nc1)
                  call gethu(chebt(1,1,lla,t2),n2,chebt(1,1,lla,t1),n1)
                  if (ia == 1) then
                     call system_clock(nc2)
                     nc = nc + nc2 - nc1
                     cnt = cnt + 1
                  end if
                  do i = 1,n1,1
                     do lb = 1,nstt(z(cluster(i))),1
                        chebt(lb,i,lla,t1) = (chebt(lb,i,lla,t1) -  arec(n-1,la)*chebt(lb,i,lla,t2) -  brec(n-1,la)*chebt(lb,i,lla,t3))/ brec(n,la)
                     enddo
                  enddo
               enddo

            endif

            lla0 = lla0 + nma

         enddo

         ja = aptr(ia)
         ja0 = ja
         ib = bptr(ja)
         do while (ib /= eol)
            dj = ja-ja0+1
            jb = decipher(ib)
            call states(z(ib),nlb,nsttb,llistb)
            lla0 = 0
            do la = 1,nla,1
               nma = 2*llista(la)+1
               fla = sqrt(real(nma, dp))
               if (brec(nrec+1,la) > zero) then
                  do lla = lla0+1,lla0+nma,1
                     ! lla = lla0 + ma
                     do lb = 1,nsttb,1
                        zeta(n,lb,dj,lla) = & 
     &                           chebt(lb,jb,lla,t1)*fla
                     enddo
                  enddo
               endif
               lla0 = lla0 + nma
            enddo
            ja = ja + 1
            ib = bptr(ja)
         enddo

         tmp = t3
         t3 = t2
         t2 = t1
         t1 = tmp

      enddo
      if (ia == 1) then 
         call system_clock(c2)
         print '("part3:",3(x,f16.12),x,i0)', real(c2-c1,dp)/real(cr,dp), real(nc,dp)/real(cr,dp), real(c2-c1-nc,dp)/real(cr,dp),cnt
      end if

!
!   Part 4: Extension with shrinking cluster.
!

      n0 = nbase
      if (ia == 1) call system_clock(c1)
      nc=0
      do n = ncross+1,2*nrec,1
      print *,'np4:',n
         n3 = pcluster(n0)-1
         n2 = pcluster(n0+1)-1
         n1 = pcluster(n0)-1

         do i = n1+1,n2,1
            decipher(cluster(i)) = 0
         enddo

         lla0 = 0
         do la = 1,nla,1

            nma = 2*llista(la)+1

            if (brec(nrec+1,la) > zero) then
!         print *, 'faza 4'
               arec(n-1,la) = arec(nrec,la)
               brec(n,la) = brec(nrec+1,la)*expfac

               do lla = lla0+1,lla0+nma,1
!                   lla = lla0 + ma
                  if (ia == 1) call system_clock(nc1)
                  call gethu(chebt(1,1,lla,t2),n2,chebt(1,1,lla,t1),n1)
                  if (ia == 1) then
                     call system_clock(nc2)
                     nc = nc + nc2 - nc1
                     cnt = cnt + 1
                  end if
                  do i = 1,n1,1
                     do lb = 1,nstt(z(cluster(i))),1
                        chebt(lb,i,lla,t1) = (chebt(lb,i,lla,t1) & 
     &                        -  arec(n-1,la)*chebt(lb,i,lla,t2) & 
     &                        -  brec(n-1,la)*chebt(lb,i,lla,t3))/ & 
     &                           brec(n,la)
                     enddo
                  enddo
               enddo

            endif

            lla0 = lla0 + nma

         enddo

         ja = aptr(ia)
         ja0 = ja
         ib = bptr(ja)
         do while (ib /= eol)
            dj = ja-ja0+1
            jb = decipher(ib)
            call states(z(ib),nlb,nsttb,llistb)
            lla0 = 0
            do la = 1,nla,1
               nma = 2*llista(la)+1
               fla = sqrt(real(nma, dp))
               if (brec(nrec+1,la) > zero) then
                  do lla = lla0+1,lla0+nma,1
                     ! lla = lla0 + ma
                     do lb = 1,nsttb,1
                        zeta(n,lb,dj,lla) = & 
     &                           chebt(lb,jb,lla,t1)*fla
                     enddo
                  enddo
               endif
               lla0 = lla0 + nma
            enddo
            ja = ja + 1
            ib = bptr(ja)
         enddo

         tmp = t3
         t3 = t2
         t2 = t1
         t1 = tmp

         n0 = n0 - 1

      enddo
      if (ia == 1) then 
         call system_clock(c2)
         print '("part4:",3(x,f16.12),x,i0)', real(c2-c1,dp)/real(cr,dp), real(nc,dp)/real(cr,dp), real(c2-c1-nc,dp)/real(cr,dp), cnt
      end if
      print *,'precab cnt:', cnt
      do i = 1,n1,1
         decipher(cluster(i)) = 0
      enddo

      end

s