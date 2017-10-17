
   subroutine forcebop()

      use mod_precision
      use mod_all_scalar
      use mod_const
      use ab_io
      use mod_ham, only : cluster, pcluster, assoc_ham ! assoc_ham for clustr and pcluster
      use mod_chi
      use mod_funptr
      use topologia, only : iproc, mpmap, aproc, master, nproc
      use mod_conf
!       use mod_atom_ar, only : btype

      use mpi
   !
   !    This is a subroutine to evaluate the bond contribution to the
   !     band structure energy, and the band structure forces.
   !

      implicit none


      include "Include/Atom.array"
      include "Include/Force.array"
      include "Include/Moment.array"
      include "Include/NebList.array"
      include "Include/PosVel.array"

      real(dp) :: bondeng,hopint

      integer i,j,nmax
      integer nla,nlb
      integer ia,ib,la,lb, mib
      integer nstta,nsttb
      integer za,zb
      integer jb
      integer lla0,nma,ma
      integer :: nmb, llb0, mb, nmaxb, ja, mlistb(mxnstat)
      type(ham_conf_t),pointer :: ham

   !
   !    Declare the local arrays.
   !

   !*** A displacement vector.
      real(dp) :: dr(3)
   !*** Forces.
      real(dp) :: fab(3), dfab(3)

      real(dp) :: grad(3,mxnstat,mxnstat)

      integer :: ierr, isp
      real(dp) :: hbo, dhbo, hbor

      real(dp) :: scfcut(14),dsf(14,3)

      scfcut = 1.0_dp
      dsf = 0.0_dp

      fbs(:,:nd) = 0.0_dp
!       dmdl(:,:nd) = 0.0_dp

!       rdfbs_tot = 0.0_dp
      rdfbs = 0.0_dp

      
            ham => rtc % ham_conf 
      do isp = 1, nsp
         do ia = mpmap(iproc)+1, mpmap(iproc+1)

            call assoc_ab(ia,isp)
            call assoc_chi(ia,isp)
            call assoc_ham(ia)

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

!             ! Calculate the chi derivatives, needed for dmdl only
!             do la = 1,nchain
!                nmax = lchain(la)
!                if (term == 1) then
!                   call dchisr(lef,nmax+1,la)
!                else if ((term == 2).or.(term == 3)) then
!                   call dchinl(lef,nmax,mrec,dchia(0,la), &
!                               & dchib(0,la),diag(1,la), &
!                               & eigvec(1,1,la),kt)
!                end if
!             end do
!

         !    For each neighbor ....
            do jb = pcluster(1), pcluster(2)-1 ! avoid onsite by starting at pcuster(1)
               ib = cluster(jb)

               mib = map(ib)
               zb = z(ib)
               call states(zb,nlb,nsttb,llistb)



               dr = ad(1:3,ia) - ad(1:3,ib)

!     it may be possible to avoid computing grad 2 times in the magnetic case
               call grdmat(grad,dr,za,nstta,zb,nsttb,scfcut,dsf,ham,.false.)

               if (momflg == 2) call utranv(transm(1,1,ia),grad,transm(1,1,mib),trwork,mxnstat,nstta,nsttb)



               fab  = 0.0_dp
               dfab = 0.0_dp

               do la = 1,nstta
                  ma = mlista(la)
                  nmax = lchain(ma)
                  do lb = 1, nsttb

! may be done with dot_product instead sum, not sure which is faster
                     hbo =    sum(chia(0:nmax,ma) * darec(0:nmax,la,lb,jb)) &
                        & + 2*sum(chib(1:nmax,ma) * dbrec(1:nmax,la,lb,jb))
                     if (term == 1) hbo = hbo + 2*chib(nmax+1,ma)*dbrec(nmax+1,la,lb,jb)


                     fab = fab + hbo * grad(1:3,la,lb)

! ! !                   Needed for dmdl only
!                      dhbo =  sum(dchia(0:nmax,ma) * darec(0:nmax,la,lb,jb)) &
!                        & + 2*sum(dchib(1:nmax,ma) * dbrec(1:nmax,la,lb,jb))
!                      if (term == 1) dhbo = dhbo + 2*dchib(nmax+1,ma)*dbrec(nmax+1,la,lb,jb)

!                      print *, darec(0,la,lb,jb),  dhbo
! ! ! !                      fab = fab + darec(0,la,lb,jb) * dhbo
!
!
!                      dfab = dfab + dhbo * grad(:,la,lb)

                  enddo
               enddo



               if (.not. mag) then
                  fab = 2*fab
!                   dfab = 2*dfab
               end if

               if (mib /= 0) then
! Symmetrisation
                  fbs(:,ia) = fbs(:,ia) + fab
                  fbs(:,mib) = fbs(:,mib) - fab
!   Commented because it is not in the version from Matous but it looks like a very useful thing so I am keeping it around
!                      dmdl(:,ia) = dmdl(:,ia) + dfab
!                      dmdl(:,mib) = dmdl(:,mib) - dfab
               else
                  fbs(:,ia) = fbs(:,ia) + 2*fab

!                      dmdl(:,ia) = dmdl(:,ia) + 2*dfab
               endif

!                print *, "ia,ib",ia,ib
!                print *, fab
               
               rdfbs = rdfbs - sum(dr*fab)

!                rdfbs_tot = rdfbs_tot - sum(dr*fab)
!
!                rdfbs(1,1) = rdfbs(1,1) - fab(1)*dr(1)
!                rdfbs(2,2) = rdfbs(2,2) - fab(2)*dr(2)
!                rdfbs(3,3) = rdfbs(3,3) - fab(3)*dr(3)
!
!                rdfbs(1,2) = rdfbs(1,2) - fab(1)*dr(2)
!                rdfbs(2,3) = rdfbs(2,3) - fab(2)*dr(3)
!                rdfbs(3,1) = rdfbs(3,1) - fab(3)*dr(1)
            enddo ! ib
         end do ! ia
      end do ! isp

!       rdfbs(2,1) = rdfbs(1,2)
!       rdfbs(3,2) = rdfbs(2,3)
!       rdfbs(1,3) = rdfbs(3,1)


      if (nproc > 1) then
          call mpi_allreduce(mpi_in_place,  fbs, 3*nd, mpi_real8, mpi_sum, mpi_comm_world, ierr)
    !       call mpi_allreduce(mpi_in_place, dmdl, 3*nd, mpi_real8, mpi_sum, mpi_comm_world, ierr)
          call mpi_allreduce(mpi_in_place, rdfbs,   1, mpi_real8, mpi_sum, mpi_comm_world, ierr)
    !       call mpi_allreduce(mpi_in_place, rdfbs_tot,   1, mpi_real8, mpi_sum, mpi_comm_world, ierr)
    !       call mpi_allreduce(mpi_in_place, rdfbs,   9, mpi_real8, mpi_sum, mpi_comm_world, ierr)
      end if
      
      

!       print *, 'SUM'
!       do ia = mpmap(iproc)+1, mpmap(iproc+1)
!           print *, ia,fbs(:,ia)
!       enddo
 
   end subroutine forcebop





!
!
!                if (momflg == 1) then
!                   llb0 = 0
!                   do lb = 1,nlb
!                      nmb = 2*llistb(lb) + 1
!                      mlistb(llb0+1:llb0+nmb) = lb
!                      llb0 = llb0 + nmb
!                   enddo
!                else
!                   do lb = 1, nsttb
!                      mlistb(lb) = lb
!                   end do
!                endif
!
















!                      write (568,'("ia,ib,la,lb,bo,bo:",4(x,i5),x,f20.16)',advance='no') ia, ib, la, lb, hbo
!
!                      call assoc_ab(ib,isp)
!                      call assoc_chi(ib,isp)
!                      call assoc_ham(ib)
! !
! !
! !
!                      do ja = pcluster(1), pcluster(2)-1 ! avoid onsite by starting at pcuster(1)
!                         if (ia /= cluster(ja)) cycle
!                         mb = mlistb(lb)
!                         nmaxb = lchain(mb)
!
!                         hbor =  sum(chia(0:nmaxb,mb) * darec(0:nmaxb,lb,la,ja)) &
!                            & + 2*sum(chib(1:nmaxb,mb) * dbrec(1:nmaxb,lb,la,ja))
!                         if (term == 1) hbor = hbor + 2*chib(nmaxb+1,mb)*dbrec(nmaxb+1,lb,la,ja)
!
! !                         write (567,*) hbo, hbor
!                         hbo = 0.5_dp *( hbo + hbor)
!                      end do
!
!                      call assoc_ab(ia,isp)
!                      call assoc_chi(ia,isp)
!                      call assoc_ham(ia)

!                      write (568,'(x,f20.16)',advance='yes') hbo


!                mib = 0
