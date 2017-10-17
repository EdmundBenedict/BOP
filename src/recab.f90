 
      subroutine recab(ia, zeta, near)
          use mod_precision
          use mod_all_scalar
          use mod_const
          use ab_io
          use mod_ham
!
!    This is a routine to build the recursion coefficients.
!

      implicit none

      integer, intent(in) :: ia
      logical, intent(in) :: near ! calculate only nearest neighbour zeta or full cluster
      real(dp), intent(out) :: zeta(0:2*mrec+2,mxnstat,mxcls,mxnstat)
   
!*** Vectors constructed from powers of the Hamiltonian.
      real(dp) :: chebt(mxnstat,mxcls,mxnstat,3)
      
      include "Include/Atom.array"

      real(dp) :: dot,lng,expfac

      integer n0,n2,n3,t1,t2,t3,tmp
      integer, target :: n1
      integer jb
      integer ncross,m
      integer la,lb,n,i
      integer nla,nstta,nlb,nsttb

      real(dp), parameter :: zero = 1.0e-3_dp

      integer, pointer :: nclus

      
      if (near) then
         allocate(nclus)
         nclus = pcluster(2)-1
      else
         nclus => n1
      end if
      
      chebt(:,1,:,:) = 0.0_dp
      zeta(0,:,:,:) = 0.0_dp
      brec(0,:) = 0.0_dp
      lchain(:) = nrec
      
      do la = 1,mxnstat
         chebt(la,1,la,2) = 1.0_dp
         zeta(0,la,1,la) = 1.0_dp
      enddo

      call states(z(ia),nla,nstta,llista)

      decipher(ia) = 1

!
!    Part 1: Recursion with growing cluster.
!

      t1 = 1
      t2 = 2
      t3 = 3

      do n = 1,nbase,1

         n3 = pcluster(n-1)-1
         n2 = pcluster(n)-1
         n1 = pcluster(n+1)-1

         do i = n2+1,n1,1
            decipher(cluster(i)) = i
         enddo

         do la = 1,nstta,1

            if (lchain(la) == nrec) then

               call gethu(chebt(1,1,la,t2),n2,chebt(1,1,la,t1),n1,1)
               arec(n-1,la) = dot(chebt(1,1,la,t2), chebt(1,1,la,t1),n2*mxnstat)
               do i = 1,n2,1
                  do lb = 1,nstt(z(cluster(i))),1
                     chebt(lb,i,la,t1) = chebt(lb,i,la,t1) - arec(n-1,la)*chebt(lb,i,la,t2)
                  enddo
               enddo
               do i = 1,n3,1
                  do lb = 1,nstt(z(cluster(i))),1
                     chebt(lb,i,la,t1) = chebt(lb,i,la,t1) - brec(n-1,la)*chebt(lb,i,la,t3)
                  enddo
               enddo
               lng = dot(chebt(1,1,la,t1), chebt(1,1,la,t1), n1*mxnstat)
               lng = sqrt(lng)
               brec(n,la) = lng
               if (brec(n,la) <= zero) then
                  lchain(la) = n-1
                  m = n
                  do while (m <= 2*nrec)
                     arec(m,la) = 0.0_dp
                     brec(m,la) = 0.0_dp
                     do jb = 1, nclus
                        do lb = 1,nstt(z(cluster(jb)))
                           zeta(m,lb,jb,la) = 0.0_dp
                        enddo
                     enddo
                     m = m + 1
                  enddo
               else
                  do i = 1,n1,1
                     do lb = 1,nstt(z(cluster(i)))
                        chebt(lb,i,la,t1) = chebt(lb,i,la,t1)/lng
                     enddo
                  enddo
               endif

            endif

         enddo

         do jb = 1, nclus
            nsttb = nstt(z(cluster(jb)))
            do la = 1,nstta
               if (lchain(la) == nrec) then
                  do lb = 1,nsttb
                     zeta(n,lb,jb,la) = chebt(lb,jb,la,t1)
                  enddo
               endif
            enddo
         enddo
 
         tmp = t3
         t3 = t2
         t2 = t1
         t1 = tmp

      enddo

!
!   Part 2: Recursion at full cluster size.
!

      ncross = 2*nrec-nbase+1

      n3 = pcluster(nbase)-1
      n2 = pcluster(nbase+1)-1
      n1 = pcluster(nbase+1)-1

      do n = nbase+1,nrec+1

         do la = 1,nstta

            if (lchain(la) == nrec) then

               call gethu(chebt(1,1,la,t2),n2,chebt(1,1,la,t1),n1,1)
               arec(n-1,la) = dot(chebt(1,1,la,t2), chebt(1,1,la,t1), n2*mxnstat)
               do i = 1,n2
                  do lb = 1,nstt(z(cluster(i)))
                     chebt(lb,i,la,t1) = chebt(lb,i,la,t1) - arec(n-1,la)*chebt(lb,i,la,t2)
                  enddo
               enddo
               do i = 1,n3
                  do lb = 1,nstt(z(cluster(i)))
                     chebt(lb,i,la,t1) = chebt(lb,i,la,t1) - brec(n-1,la)*chebt(lb,i,la,t3)
                  enddo
               enddo
               lng = dot(chebt(1,1,la,t1), chebt(1,1,la,t1),n1*mxnstat)
               lng = sqrt(lng)
               brec(n,la) = lng
               if (brec(n,la) <= zero) then
                  lchain(la) = n-1
                  m = n
                  do while (m <= 2*nrec)
                     arec(m,la) = 0.0_dp
                     brec(m,la) = 0.0_dp
                     do jb = 1, nclus
                        nsttb = nstt(z(cluster(jb)))
                        do lb = 1,nsttb
                           zeta(m,lb,jb,la) = 0.0_dp
                        enddo
                     enddo
                     m = m + 1
                  enddo
               else
                  do i = 1,n1
                     do lb = 1,nstt(z(cluster(i)))
                        chebt(lb,i,la,t1) = chebt(lb,i,la,t1)/lng
                     enddo
                  enddo
               endif

           endif

         enddo
 
         do jb = 1, nclus
            nsttb = nstt(z(cluster(jb)))
            do la = 1,nstta
               if (lchain(la) == nrec) then
                  do lb = 1,nsttb
                     zeta(n,lb,jb,la) = chebt(lb,jb,la,t1)
                  enddo
               endif
            enddo
         enddo

         tmp = t3
         t3 = t2
         t2 = t1
         t1 = tmp

         n3 = pcluster(nbase+1)-1

      enddo

!
!   Part 3: Extension at constant cluster size.
!

      expfac = 1.0_dp

      do n = nrec+2,ncross

         do la = 1,nstta
            if (brec(nrec+1,la) > zero) then
               arec(n-1,la) = arec(nrec,la)
               brec(n,la) = brec(nrec+1,la)*expfac
               call gethu(chebt(1,1,la,t2),n2,chebt(1,1,la,t1),n1,1)
               do i = 1,n1
                  do lb = 1,nstt(z(cluster(i)))
                     chebt(lb,i,la,t1) = (chebt(lb,i,la,t1) & 
     &                       -  arec(n-1,la)*chebt(lb,i,la,t2) & 
     &                       -  brec(n-1,la)*chebt(lb,i,la,t3))/ & 
     &                      brec(n,la)
                  enddo
               enddo
            endif
         enddo

         do jb = 1, nclus
            nsttb = nstt(z(cluster(jb)))
            do la = 1,nstta
               if (brec(nrec+1,la) > zero) then
                  do lb = 1,nsttb,1
                     zeta(n,lb,jb,la) = chebt(lb,jb,la,t1)
                  enddo
               endif
            enddo
         enddo

         tmp = t3
         t3 = t2
         t2 = t1
         t1 = tmp

      enddo

!
!   Part 4: Extension with shrinking cluster.
!

      n0 = nbase

      do n = ncross+1,2*nrec

         n3 = pcluster(n0)-1
         n2 = pcluster(n0+1)-1
         n1 = pcluster(n0)-1

         do i = n1+1,n2
            decipher(cluster(i)) = 0
         enddo

         do la = 1,nstta
            if (brec(nrec+1,la) > zero) then
               arec(n-1,la) = arec(nrec,la)
               brec(n,la) = brec(nrec+1,la)*expfac
               call gethu(chebt(1,1,la,t2),n2,chebt(1,1,la,t1),n1,1)
               do i = 1,n1
                  do lb = 1,nstt(z(cluster(i)))
                     chebt(lb,i,la,t1) = (chebt(lb,i,la,t1) & 
     &                        -  arec(n-1,la)*chebt(lb,i,la,t2) & 
     &                        -  brec(n-1,la)*chebt(lb,i,la,t3))/ & 
     &                      brec(n,la)
                  enddo
               enddo
            endif
         enddo

         do jb = 1, nclus
            nsttb = nstt(z(cluster(jb)))
            do la = 1,nstta
               if (brec(nrec+1,la) > zero) then
                  do lb = 1,nsttb
                     zeta(n,lb,jb,la) = chebt(lb,jb,la,t1)
                  enddo
               endif
            enddo
         enddo

         tmp = t3
         t3 = t2
         t2 = t1
         t1 = tmp

         n0 = n0 - 1

      enddo

      do i = 1,n1
         decipher(cluster(i)) = 0
      enddo

      end

