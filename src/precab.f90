
   subroutine precab(ia, zeta, near)

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
   
   integer, intent(in) :: ia
   logical, intent(in) :: near ! calculate only nearest neighbour zeta or full cluster
   real(dp), intent(out) :: zeta(0:2*mrec+2,mxnstat,mxcls,mxnstat)
   
   
   
!*** Vectors constructed from powers of the Hamiltonian.
   real(dp) :: chebt(mxnstat,mxcls,mxnstat,3)
      
   include "Include/Atom.array"


   real(dp) :: dot,lng,inv_lng,expfac,fla

   integer n0,n2,n3,t1,t2,t3,tmp
   integer, target :: n1, tnclus
   integer jb
   integer ncross,m
   integer lla,ma,lla0,nma,la,lb,n,ib,i
   integer nla,nstta,nsttb
!    integer(8) :: c1,c2,nc1,nc2,nc

   real(dp), parameter :: zero = 1.0e-3_dp
   integer, pointer :: nclus
   
   if (near) then
      tnclus = pcluster(2)-1 
      nclus => tnclus
   else
      nclus => n1
   end if
   
   call states(z(ia),nla,nstta,llista)
   
   lla0 = 0
   chebt(:,1,:,:) = 0.0_dp
!    chebt = 0.0_dp

!    zeta(0,:,:,:) = 0.0_dp
   zeta = 0.0_dp
   
   do la = 1,nla
      brec(0,la) = 0.0_dp
      lchain(la) = nrec
      nma = 2*llista(la)+1
      fla = sqrt(real(nma, dp))
      do lla = lla0+1,lla0+nma
         ! lla = lla0 + ma
         chebt(lla,1,lla,2) = 1.0_dp/fla
         zeta(0,lla,1,lla) = 1.0_dp
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
   do n = 1,nbase
      n3 = pcluster(n-1)-1
      n2 = pcluster(n)-1
      n1 = pcluster(n+1)-1
     
      do i = n2+1,n1
         decipher(cluster(i)) = i
      enddo

      lla0 = 0
      do la = 1,nla
         nma = 2*llista(la)+1
         if (lchain(la) == nrec) then
            call gethu(chebt(1,1,lla0+1,t2),n2,chebt(1,1,lla0+1,t1),n1,nma)

            arec(n-1,la) = 0.0_dp
            do lla = lla0+1,lla0+nma
               arec(n-1,la) = arec(n-1,la) + dot(chebt(1,1,lla,t2), chebt(1,1,lla,t1), n2*mxnstat)
            enddo

            lng = 0.0_dp
            do lla = lla0+1,lla0+nma
               call daxpy(n2*mxnstat, -arec(n-1,la), chebt(1,1,lla,t2), 1, chebt(1,1,lla,t1), 1 )
               call daxpy(n3*mxnstat, -brec(n-1,la), chebt(1,1,lla,t3), 1, chebt(1,1,lla,t1), 1 )
               lng = lng + dot(chebt(1,1,lla,t1), chebt(1,1,lla,t1),n1*mxnstat)
            end do
            
            lng = sqrt(lng)
            brec(n,la) = lng
            
            if (brec(n,la) <= zero) then
               lchain(la) = n-1
               do m = n,2*nrec
                  arec(m,la) = 0.0_dp
                  brec(m,la) = 0.0_dp
                  do jb = 1, nclus
                     nsttb = nstt(z(cluster(jb)))
                     zeta(m,1:nsttb,jb,lla0+1:lla0+nma) = 0.0_dp
                  enddo
               enddo
            else
               inv_lng = 1.0_dp/lng
               do lla = lla0+1,lla0+nma
                  do i = 1,n1
                     do lb = 1,nstt(z(cluster(i)))
                        chebt(lb,i,lla,t1) = inv_lng*chebt(lb,i,lla,t1)
                     enddo
                  enddo
               enddo
            endif
         endif
         lla0 = lla0 + nma
      enddo

      
      
      do jb = 1, nclus
         nsttb = nstt(z(cluster(jb)))
         lla0 = 0
         do la = 1,nla
            nma = 2*llista(la)+1
            if (lchain(la) == nrec) then
               fla = sqrt(real(nma, dp))
               zeta(n,1:nsttb,jb,lla0+1:lla0+nma) = chebt(1:nsttb,jb,lla0+1:lla0+nma,t1)*fla
            end if
            lla0 = lla0 + nma
         enddo
      end do

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
      lla0 = 0
      do la = 1,nla

         nma = 2*llista(la)+1

         if (lchain(la) == nrec) then
            arec(n-1,la) = 0.0_dp
            call gethu(chebt(1,1,lla0+1,t2),n2,chebt(1,1,lla0+1,t1),n1,nma)
            do lla = lla0+1,lla0+nma,1
               arec(n-1,la) = arec(n-1,la) + dot(chebt(1,1,lla,t2), chebt(1,1,lla,t1),n2*mxnstat)
            enddo

            do lla = lla0+1,lla0+nma               
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
                 do i = 1,n3
                     do lb = 1,nstt(z(cluster(i)))
                         chebt(lb,i,lla,t1) = chebt(lb,i,lla,t1) - brec(n-1,la)*chebt(lb,i,lla,t3)
                     enddo
                 enddo
               else
                 call daxpy(n3*mxnstat, -brec(n-1,la), chebt(1,1,lla,t3), 1, chebt(1,1,lla,t1), 1 )                  
               endif
            enddo

            lng = 0.0_dp
            do lla = lla0+1,lla0+nma
               lng = lng + dot(chebt(1,1,lla,t1), chebt(1,1,lla,t1),n1*mxnstat)
            enddo
            lng = sqrt(lng)
            brec(n,la) = lng

            if (brec(n,la) <= zero) then
               lchain(la) = n-1
               do m=n,2*nrec
                  arec(m,la) = 0.0_dp
                  brec(m,la) = 0.0_dp
                  do jb = 1, nclus
                     nsttb = nstt(z(cluster(jb)))
                     zeta(m,1:nsttb,jb,lla0+1:lla0+nma) = 0.0_dp
                  enddo
               enddo
            else
               do lla = lla0+1,lla0+nma
                  do i = 1,n1
                     do lb = 1,nstt(z(cluster(i)))
                        chebt(lb,i,lla,t1) = chebt(lb,i,lla,t1)/lng
                     enddo
                  enddo
               enddo
            endif
         endif
         lla0 = lla0 + nma
      enddo

            
      do jb = 1, nclus
         nsttb = nstt(z(cluster(jb)))
         lla0 = 0
         do la = 1,nla
            nma = 2*llista(la)+1
            if (brec(nrec+1,la) > zero) then
               fla = sqrt(real(nma, dp))
               zeta(n,1:nsttb,jb,lla0+1:lla0+nma) = chebt(1:nsttb,jb,lla0+1:lla0+nma,t1)*fla
            end if
            lla0 = lla0 + nma
         enddo
      end do
      
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
   do n = nrec+2,ncross,1
      lla0 = 0
      do la = 1,nla

         nma = 2*llista(la)+1

         if (brec(nrec+1,la) > zero) then
            arec(n-1,la) = arec(nrec,la)
            brec(n,la) = brec(nrec+1,la)*expfac
            call gethu(chebt(1,1,lla0+1,t2),n2,chebt(1,1,lla0+1,t1),n1,nma)
            do lla = lla0+1,lla0+nma
               do i = 1,n1
                  do lb = 1,nstt(z(cluster(i)))
                     chebt(lb,i,lla,t1) = (chebt(lb,i,lla,t1) -  arec(n-1,la)*chebt(lb,i,lla,t2) &
                                     & -  brec(n-1,la)*chebt(lb,i,lla,t3))/ brec(n,la) 
                  enddo
               enddo
            enddo

         endif

         lla0 = lla0 + nma

      enddo
            
      do jb = 1, nclus
         nsttb = nstt(z(cluster(jb)))
         lla0 = 0
         do la = 1,nla
            nma = 2*llista(la)+1
            if (brec(nrec+1,la) > zero) then
               fla = sqrt(real(nma, dp))
               zeta(n,1:nsttb,jb,lla0+1:lla0+nma) = chebt(1:nsttb,jb,lla0+1:lla0+nma,t1)*fla
            end if
            lla0 = lla0 + nma
         enddo
      end do
      

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

      do i = n1+1,n2,1
         decipher(cluster(i)) = 0
      enddo

      lla0 = 0
      do la = 1,nla
         nma = 2*llista(la)+1
         if (brec(nrec+1,la) > zero) then
            arec(n-1,la) = arec(nrec,la)
            brec(n,la)   = brec(nrec+1,la)*expfac
            call gethu(chebt(1,1,lla0+1,t2),n2,chebt(1,1,lla0+1,t1),n1,nma)
            do lla = lla0+1,lla0+nma
               do i = 1,n1
                  nsttb = nstt(z(cluster(i)))
                  chebt(1:nsttb,i,lla,t1) = (chebt(1:nsttb,i,lla,t1) &
                            & - arec(n-1,la)*chebt(1:nsttb,i,lla,t2) &
                            & - brec(n-1,la)*chebt(1:nsttb,i,lla,t3))/brec(n,la)
               enddo
            enddo
         endif
         lla0 = lla0 + nma
      enddo

      do jb = 1, nclus
         nsttb = nstt(z(cluster(jb)))
         lla0 = 0
         do la = 1,nla
            nma = 2*llista(la)+1
            if (brec(nrec+1,la) > zero) then
               fla = sqrt(real(nma, dp))
               zeta(n,1:nsttb,jb,lla0+1:lla0+nma) = chebt(1:nsttb,jb,lla0+1:lla0+nma,t1)*fla
            end if
            lla0 = lla0 + nma
         enddo
      end do      
      
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

! 
!    subroutine daxpyl(n, a, x, sx, y, sy)
!       implicit none
!       integer, intent(in) :: n, sx, sy
!       real(8), intent(in) :: a, x(0:n-1)
!       real(8), intent(inout) :: y(0:n-1)
!       
!       integer :: i
!       
!       do i = 0, n-1
!          y(i*sy) = y(i*sy) + a*x(i*sx)
!       end do
!       
!       
!    end subroutine daxpyl
!    