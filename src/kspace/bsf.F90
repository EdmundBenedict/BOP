! 
   subroutine bsf()
      use mod_precision
      use mod_all_scalar
      use mod_const
      use mod_kspace
      use mod_atom_ar, only : btype
      use topologia
      use mod_conf
      
#ifdef MPI
      use mpi
#endif

!    evaluate the band structure energy contribution to the forces.


      implicit none

      include "../Include/Atom.array"
      include "../Include/Force.array"
      include "../Include/NebList.array"
      include "../Include/PosVel.array"

      real(dp) :: rcc,icc,rbo
      real(dp) :: rfac,ifac,phase,rphs,iphs

      integer ia,ib,ja,mapb,i,j
      integer ik,in
      integer iha,ihb, iha0, ihb0
      integer la,lb,nstta,nsttb
      integer za,zb
      integer ptr, count1, test_atom
      type(ham_conf_t),pointer :: ham

      logical usepot

      real(dp) :: fab(3),sgrad(3,mxnstat,mxnstat)
! 
   !*** The difference between two atomic positions.
      real(dp) :: dr(3)
   !*** Screening function
      real(dp) :: scfcut(14),dsf(3,14)

      complex(dp) :: cc,ccj
      
      integer :: ierr, isp
     

      
      fbs(:,:nd) = 0.0_dp
      rdfbs = 0.0_dp
      band_forces = 0.0_dp

      scfcut = 1.0_dp
      dsf = 0.0_dp
      
!       this is slow
                     
      ham => rtc % ham_conf 
      
      do ia = 1,nd

         za = z(ia)
         nstta = nstt(za)
!          call states(za,nla,nstta,llista)
         ja = aptr(ia)
         ptr = 1

         iha0 = khpos(ia)
         do while(bptr(ja) /= eol)

            ib = bptr(ja)
            mapb = map(ib)
            zb = z(mapb)
            nsttb = nstt(zb)
!             call states(zb,nlb,nsttb,llistb)
            
            ihb0 = khpos(mapb)
            
            dr = ad(1:3,ia) - ad(1:3,ib)

            call grdmat(kgrad,dr,za,nstta,zb,nsttb,scfcut,dsf,ham,.false.)
            
            if (ham % ovl) call grdmat(sgrad,dr,za,nstta,zb,nsttb,scfcut,dsf,ham,.true.)

            fab = 0.0_dp

            do isp = 1, nsp
               do ik = mpmap(iproc)+1, mpmap(iproc+1)
                  phase = sum(kk(1:3,ik)*dr)
                  rphs = cos(phase)
                  iphs = sin(phase)
                  do in = 1, knh                  
!                      if (occ(in,ik) < 1.0e-6_dp) exit
                     do lb = 1, nsttb
                        ccj = dconjg(kpsi(ihb0+lb,in,ik,isp))
                        do la = 1, nstta
                           cc = kpsi(iha0+la,in,ik,isp)*ccj                                    
                           rbo = occ(in,ik,isp) * (rphs*real(cc) + iphs*aimag(cc)) ! f(n) Re( e^{-ikr} C_{ni} C^*_{nj}  )
                           if (.not. ham % ovl) then
                              fab = fab + kgrad(:,la,lb)*rbo
                           else 
                              fab = fab + (kgrad(:,la,lb) - sgrad(:,la,lb)*enk(in,ik,isp))*rbo
                           endif 
!                            print *, kgrad(:,la,lb), rbo
                        enddo
                     enddo
                  end do
               end do
            end do

            fab = -fab
            if (.not.mag) fab = 2 * fab
            
#ifdef MPI            
            call mpi_allreduce(mpi_in_place, fab, 3, mpi_real8, mpi_sum, mpi_comm_world, ierr)
#endif            

!             print *, fab
            fbs(:,ia) = fbs(:,ia) + fab

            fbs(:,mapb) = fbs(:,mapb) - fab

            rdfbs = rdfbs - sum(dr*fab)

!             band_forces(:,ptr,ia) = band_forces(:,ptr,ia) + fab
!             
!             do count1 = 0, mxnnb-1
!                test_atom = bptr(aptr(mapb)+count1)
!                if (((test_atom <= nd).and.(test_atom == ia)).or. &
!    &                (test_atom > nd).and.(map(test_atom) == ia)) then
!                   band_forces(:,count1+1,mapb) = band_forces(:,count1+1,mapb) - fab(:)
!                   exit
!                endif
!             enddo

            ptr = ptr + 1

            ja = ja + 1

         enddo

      enddo

      
      do ia =1,nd
          print *, fbs(:,ia)
      enddo
      
!          call system_clock(c2)
      
!          write(6,'("bsf wallclock time: ",g,"s")'), real(c2 - c1,8)/real(cr,8)

   end subroutine bsf
! 
