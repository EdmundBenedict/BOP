 
      subroutine matel(z1,nstt1,z2,nstt2,dr,subh,de,scf,ham,ovl)

!
!    This is a routine to evaluate the elements of the tight-binding matrix.
!
!=============================================================================
!    Note that the angular momentum contributions are indexed as follows:
!
!		s		0
!
!		x		0
!		y		1
!		z		2
!
!		3zz-rr		0
!		xx-yy		1
!		xy		2
!		xz		3
!		yz		4
!
!=============================================================================
!
          use mod_precision
          use mod_const, only: mxnstat,ml,sqrt3
          use mod_conf, only: rtc,ham_conf_t,bond_conf_t
          use mod_atom_ar,only:zatype,asort,bsort,btype

      implicit none
      
      integer, intent(in) :: z1,z2
      real(dp), intent(in) :: dr(3), de, scf(14)
      real(dp), intent(out) :: subh(mxnstat,mxnstat)
      
      type(ham_conf_t),target :: ham
      type(bond_conf_t),pointer:: bnd
      real(dp) :: es,ep,ed
      real(dp) :: vsss,vsps,vpss,vpps,vppp
      real(dp) :: vsds,vdss,vpds,vdps,vpdp,vdpp
      real(dp) :: vdds,vddp,vddd
      real(dp) :: magdr,l,m,n

!       integer ia,ib
      integer nstt1,nstt2,nl1,nl2
      integer i1,i2,j1,j2,nstt2b,i
      logical :: ovl
! This switch is for embedding, ovl to switch to the overlap matrix.
      
      
!       real(dp) :: bnddat(15,mbtype)
!       real(dp) :: bndscl(14,15,mbtype)
      real(dp) :: rfun(0:0,14)




      integer :: llist1(ml), llist2(ml)

      real(8) :: ll, mm, nn, inv_magdr, epde, edde, &
      & vpps_sub_vppp, ll_add_mm, ll_sub_mm, &
      & sqrt3_x_vsds, sqrt3_x_vpds, sqrt3_x_vdss, neg_sqrt3_x_vdps, &
      & xyz, nn_sub_half__ll_add_mm, one_sub_two_ll, one_sub_two_mm, one_sub_two_nn, &
      & term1, term2ll, term2mm, term2nn, term3, &
      & sqrt3_x_vdds, sqrt3_x_vddp, sqrt3_x_vddd, term2, term4, term5

!         integer, save :: counter=0

!
!    Find some parameters.
!
!           call states(z1,nl1,nstt1,llist1)
!           call states(z2,nl2,nstt2,llist2)

      magdr = sqrt(dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3))

      subh = 0.0_dp

      
      nl1 = ham % a(asort(zatype(z1))) %b % norbs
      llist1 = ham % a(asort(zatype(z1))) %b % orblist
!       if (.not. dh) then
        nl2 = ham % a(asort(zatype(z2))) %b % norbs
        llist2 = ham % a(asort(zatype(z2))) %b % orblist
!       else
!         call states(z2,nl2,nstt2,llist2)
!       endif
      

      if (magdr > 1.0e-6_dp) then
         
         bnd => ham % b(bsort(btype(z1,z2)))
!         print '("intersite: z1: ",i0,"  z2: ",i0,"  dr:",3(x,f12.6))', z1,z2,dr
!         counter = counter + 1
!         print '("c: ", i8, 3x, "z1: ", i8, 3x, "z2: ", i8, 3x, "dr: ", 3f12.6)',  counter, z1, z2, dr

!
!       Find an off-diagonal sub-matrix.
!
!       Evaluate the radial functions.
!
! 
!          call radf(magdr,vsss,vsps,vpss,vpps,vppp, & 
!      &             vsds,vdss,vpds,vdps,vpdp,vdpp, & 
!      &             vdds,vddp,vddd,bnddat(1,btype(z1,z2)), & 
!      &             bndscl(1,1,btype(z1,z2)))
! 
! 
! ! Apply the screening function
! !       if (z1 /= z2) then
! !          print *, 'z12sss', z1,z2
! !          do i =1,15
! !             print *,bnddat(i,btype(z1,z2)),"|",bndscl(:,i,btype(z1,z2))
! !          end do
! !          
! !       end if
!          vsss = vsss * scf(1)
!          vsps = vsps * scf(2)
!          vpss = vpss * scf(3)
!          vpps = vpps * scf(4)
!          vppp = vppp * scf(5)
!          vsds = vsds * scf(6)
!          vdss = vdss * scf(7)
!          vpds = vpds * scf(8)
!          vdps = vdps * scf(9)
!          vpdp = vpdp * scf(10)
!          vdpp = vdpp * scf(11)
!          vdds = vdds * scf(12)
!          vddp = vddp * scf(13)
!          vddd = vddd * scf(14)

         call eval_rfun(z1, z2, magdr, 0, rfun,bnd,ovl)
         
         
         nullify(bnd)
         


         vsss = rfun(0, 1)
         vsps = rfun(0, 2)
         vpss = rfun(0, 3)
         vpps = rfun(0, 4)
         vppp = rfun(0, 5)
         vsds = rfun(0, 6)
         vdss = rfun(0, 7)
         vpds = rfun(0, 8)
         vdps = rfun(0, 9)
         vpdp = rfun(0,10)
         vdpp = rfun(0,11)
         vdds = rfun(0,12)
         vddp = rfun(0,13)
         vddd = rfun(0,14)
         
         if (rtc%scf /= 0) then
            vsss = vsss * scf( 1)
            vsps = vsps * scf( 2)
            vpss = vpss * scf( 3)
            vpps = vpps * scf( 4)
            vppp = vppp * scf( 5)
            vsds = vsds * scf( 6)
            vdss = vdss * scf( 7)
            vpds = vpds * scf( 8)
            vdps = vdps * scf( 9)
            vpdp = vpdp * scf(10)
            vdpp = vdpp * scf(11)
            vdds = vdds * scf(12)
            vddp = vddp * scf(13)
            vddd = vddd * scf(14)
         end if
         
         
         
         
      
      
      
      
      
      
      
! 
! 
!       print '("vsss: ", f12.6)', vsss
!       print '("vsps: ", f12.6)', vsps
!       print '("vpss: ", f12.6)', vpss
!       print '("vpps: ", f12.6)', vpps
!       print '("vppp: ", f12.6)', vppp
!       print '("vsds: ", f12.6)', vsds
!       print '("vdss: ", f12.6)', vdss
!       print '("vpds: ", f12.6)', vpds
!       print '("vdps: ", f12.6)', vdps
!       print '("vpdp: ", f12.6)', vpdp
!       print '("vdpp: ", f12.6)', vdpp
!       print '("vdds: ", f12.6)', vdds
!       print '("vddp: ", f12.6)', vddp
!       print '("vddd: ", f12.6)', vddd
! 
! 
! 
! 
! !
!       Set up the off-diagonal sub-matrix.
!
         inv_magdr = 1.0_dp/magdr
         l = dr(1)*inv_magdr
         m = dr(2)*inv_magdr
         n = dr(3)*inv_magdr

         ll = l*l
         mm = m*m
         nn = n*n

         ll_add_mm = ll+mm
         ll_sub_mm = ll-mm
        
         j1 = 1
         do i1 = 1,nl1

            select case (llist1(i1))
            
            case (1)

               j2 = 1
               do i2 = 1,nl2,1

                  select case (llist2(i2))
                  case (1)
! pp
                     vpps_sub_vppp = vpps-vppp
                     
                     subh(j1,j2)     = ll*vpps_sub_vppp + vppp
                     subh(j1,j2+1)   = l*m*vpps_sub_vppp
                     subh(j1,j2+2)   = l*n*vpps_sub_vppp
                     subh(j1+1,j2)   = m*l*vpps_sub_vppp
                     subh(j1+1,j2+1) = mm*vpps_sub_vppp + vppp
                     subh(j1+1,j2+2) = m*n*vpps_sub_vppp
                     subh(j1+2,j2)   = n*l*vpps_sub_vppp
                     subh(j1+2,j2+1) = n*m*vpps_sub_vppp
                     subh(j1+2,j2+2) = nn*vpps_sub_vppp + vppp
                     j2 = j2 + 3

                  case (2)
! pd
                     sqrt3_x_vpds = sqrt3*vpds
                     xyz = l*m*n*(sqrt3_x_vpds - 2*vpdp)
                     nn_sub_half__ll_add_mm = nn-0.5_dp*ll_add_mm
                     one_sub_two_ll = 1-2*ll
                     one_sub_two_mm = 1-2*mm
                     one_sub_two_nn = 1-2*nn
                     term1 = (nn_sub_half__ll_add_mm*vpds - sqrt3*nn*vpdp)
                     term2ll = (ll*sqrt3_x_vpds + one_sub_two_ll*vpdp)
                     term2mm = (mm*sqrt3_x_vpds + one_sub_two_mm*vpdp)
                     term2nn = (nn*sqrt3_x_vpds + one_sub_two_nn*vpdp)
                     term3 = ll_sub_mm*(0.5_dp*sqrt3_x_vpds - vpdp)
                     
                     subh(j1,j2)     = l*term1
                     subh(j1,j2+1)   = l*(term3 + vpdp)
                     subh(j1,j2+2)   = m*term2ll
                     subh(j1,j2+3)   = n*term2ll
                     subh(j1,j2+4)   = xyz
                     subh(j1+1,j2)   = m*term1
                     subh(j1+1,j2+1) = m*(term3 - vpdp )
                     subh(j1+1,j2+2) = l*term2mm
                     subh(j1+1,j2+3) = xyz
                     subh(j1+1,j2+4) = n*term2mm
                     subh(j1+2,j2)   = n*(nn_sub_half__ll_add_mm*vpds + sqrt3*ll_add_mm*vpdp)
                     subh(j1+2,j2+1) = n*term3
                     subh(j1+2,j2+2) = xyz
                     subh(j1+2,j2+3) = l*term2nn
                     subh(j1+2,j2+4) = m*term2nn
                     j2 = j2 + 5
                  
                  case (0)
! ps
                     subh(j1,j2)   = -l*vpss
                     subh(j1+1,j2) = -m*vpss
                     subh(j1+2,j2) = -n*vpss
                     j2 = j2 + 1


                  end select

               enddo

               j1 = j1 + 3

            case (2)

               j2 = 1
               do i2 = 1,nl2,1

                  select case (llist2(i2))
                  
                  case (1)
! dp
                     neg_sqrt3_x_vdps = -sqrt3*vdps
                     
                     xyz = l*m*n*(neg_sqrt3_x_vdps + 2*vdpp)
                     
                     nn_sub_half__ll_add_mm = nn-0.5_dp*ll_add_mm
                     one_sub_two_ll = 1-2*ll
                     one_sub_two_mm = 1-2*mm
                     one_sub_two_nn = 1-2*nn
                     
                     term1 = (-nn_sub_half__ll_add_mm*vdps + sqrt3*nn*vdpp)
                     term2ll = (ll*neg_sqrt3_x_vdps - one_sub_two_ll*vdpp)
                     term2mm = (mm*neg_sqrt3_x_vdps - one_sub_two_mm*vdpp)
                     term2nn = (nn*neg_sqrt3_x_vdps - one_sub_two_nn*vdpp)
                     term3 = ll_sub_mm*(0.5_dp*neg_sqrt3_x_vdps + vdpp)
                     
                     
                     subh(j1,j2)     = l*term1
                     subh(j1,j2+1)   = m*term1
                     subh(j1,j2+2)   = -n*(nn_sub_half__ll_add_mm*vdps + sqrt3*ll_add_mm*vdpp)
                     subh(j1+1,j2)   = l*(term3 - vdpp )
                     subh(j1+1,j2+1) = m*(term3 + vdpp)
                     subh(j1+1,j2+2) = n*term3
                     subh(j1+2,j2)   = m*term2ll
                     subh(j1+2,j2+1) = l*term2mm
                     subh(j1+2,j2+2) = xyz
                     subh(j1+3,j2)   = n*term2ll
                     subh(j1+3,j2+1) = xyz
                     subh(j1+3,j2+2) = l*term2nn
                     subh(j1+4,j2)   = xyz
                     subh(j1+4,j2+1) = n*term2mm
                     subh(j1+4,j2+2) = m*term2nn
                     
                     j2 = j2 + 3

                  case (2)
! dd
                     nn_sub_half__ll_add_mm = nn-0.5_dp*ll_add_mm
                     
                     sqrt3_x_vdds = sqrt3 * vdds
                     sqrt3_x_vddp = sqrt3 * vddp
                     sqrt3_x_vddd = sqrt3 * vddd
                     
                     term1 = (0.5_dp*nn_sub_half__ll_add_mm*sqrt3_x_vdds & 
     &                               - nn*sqrt3_x_vddp + 0.25_dp*(1+nn)*sqrt3_x_vddd)
                     term2 = (nn_sub_half__ll_add_mm*sqrt3_x_vdds & 
     &                               + (ll_add_mm-nn)*sqrt3_x_vddp - 0.5_dp*ll_add_mm*sqrt3_x_vddd)
                     term3 = ll_sub_mm*(1.5_dp*vdds - 2*vddp + 0.5_dp*vddd)
                     term4 = (3.0_dp*vdds - 4.0_dp*vddp + vddd)
                     term5 =  vddp -vddd
                     
                     
                     subh(j1,j2)     = nn_sub_half__ll_add_mm * nn_sub_half__ll_add_mm * vdds & 
     &                               + ll_add_mm * (3.0_dp * nn *  vddp + 0.75_dp * ll_add_mm * vddd)
!                                       ll_add_mm * 0.75d0 * (nn *  (4.0d0 * vddp - vddd)+ vddd  )
     
                     
                     subh(j1,j2+1)   = ll_sub_mm* term1
                     subh(j1,j2+2)   = 2*l*m*term1
                     subh(j1,j2+3)   = l*n*term2
                     subh(j1,j2+4)   = m*n*term2


!                     WRITE(79,'(2I3,F10.6)') J1,J2,SUBH(J1,J2)
!                     WRITE(79,'(2I3,F10.6)') J1,J2+1,SUBH(J1,J2+1)
!                     WRITE(79,'(2I3,F10.6)') J1,J2+2,SUBH(J1,J2+2)
!                     WRITE(79,'(2I3,F10.6)') J1,J2+3,SUBH(J1,J2+3)
!                     WRITE(79,'(2I3,F10.6)') J1,J2+4,SUBH(J1,J2+4)


                     subh(j1+1,j2)   = subh(j1,j2+1)
                     subh(j1+1,j2+1) = ll_sub_mm*ll_sub_mm*term4*0.25_dp + nn*vddd + ll_add_mm*vddp
                     subh(j1+1,j2+2) = l*m*term3
                     subh(j1+1,j2+3) = n*l*(term3 + term5)
                     subh(j1+1,j2+4) = m*n*(term3 - term5)
                     subh(j1+2,j2)   = subh(j1,j2+2)
                     subh(j1+2,j2+1) = subh(j1+1,j2+2)
                     subh(j1+2,j2+2) = ll*mm*term4 + nn*vddd + ll_add_mm*vddp
                     subh(j1+2,j2+3) = m*n*( ll*term4 + term5)
                     subh(j1+2,j2+4) = l*n*( mm* term4 + term5)
                     subh(j1+3,j2)   = subh(j1,j2+3)
                     subh(j1+3,j2+1) = subh(j1+1,j2+3)
                     subh(j1+3,j2+2) = subh(j1+2,j2+3)
                     subh(j1+3,j2+3) = ll*nn*term4 + mm*vddd+ (ll+nn)*vddp 
                     subh(j1+3,j2+4) = l*m*(nn*term4 + term5)
                     subh(j1+4,j2)   = subh(j1,j2+4)
                     subh(j1+4,j2+1) = subh(j1+1,j2+4)
                     subh(j1+4,j2+2) = subh(j1+2,j2+4)
                     subh(j1+4,j2+3) = subh(j1+3,j2+4)
                     subh(j1+4,j2+4) = nn*mm*term4 + ll*vddd + (nn+mm)*vddp
                     
                     j2 = j2 + 5
                  
                  case (0)
! ds
                     sqrt3_x_vdss = sqrt3 * vdss
                     subh(j1,j2)   = (nn-0.5_dp*ll_add_mm)*vdss
                     subh(j1+1,j2) = 0.5_dp*ll_sub_mm*sqrt3_x_vdss
                     subh(j1+2,j2) = l*m*sqrt3_x_vdss
                     subh(j1+3,j2) = l*n*sqrt3_x_vdss
                     subh(j1+4,j2) = m*n*sqrt3_x_vdss
                     j2 = j2 + 1

                  
                  end select

               enddo

               j1 = j1 + 3

            case (0)

               j2 = 1
               do i2 = 1,nl2,1
                  
                  select case (llist2(i2))
                  
                  case (1)
! sp
                     subh(j1,j2)   = l*vsps
                     subh(j1,j2+1) = m*vsps
                     subh(j1,j2+2) = n*vsps
                     j2 = j2 + 3

                  case (2)
! sd
                     sqrt3_x_vsds = sqrt3*vsds
                     
                     subh(j1,j2)   = (nn-0.5_dp*ll_add_mm)*vsds
                     subh(j1,j2+1) = 0.5_dp*ll_sub_mm*sqrt3_x_vsds
                     subh(j1,j2+2) = l*m*sqrt3_x_vsds
                     subh(j1,j2+3) = l*n*sqrt3_x_vsds
                     subh(j1,j2+4) = m*n*sqrt3_x_vsds
                     j2 = j2 + 5
                  case (0)
! ss
!                   if (z1 /= z2) print *, 'z12sss', z1,z2, vsss
                     subh(j1,j2) = vsss
                     j2 = j2 + 1


                  end select

               enddo

               j1 = j1 + 1



            end select

         enddo
        
      else
!
!    If this is on-site sub-matrix, then find the on-site terms.
!
!      
!           print '("onsite: z1: ",i0,"  z2: ",i0,"  dr:",3(x,f12.6))', z1,z2,dr
!          call onsite(z1,es,ep,ed)

        if (.not. ovl) then

          es = ham % a(asort(zatype(z1))) % b % es
          ep = ham % a(asort(zatype(z1))) % b % ep
          ed = ham % a(asort(zatype(z1))) % b % ed
          
          epde = ep+de
          edde = ed+de
          do i1 = 1,nstt1
              do i2 = i1+1,nstt2
                subh(i1,i2) = 0.0_dp
                subh(i2,i1) = 0.0_dp
              enddo
          enddo
          j1 = 1
          do i1 = 1,nl1,1
              select case(llist1(i1))
              
              case (1)
                subh(j1,j1) = epde
                subh(j1+1,j1+1) = epde
                subh(j1+2,j1+2) = epde
                j1 = j1 + 3
              case (2)
                subh(j1,j1) = edde
                subh(j1+1,j1+1) = edde
                subh(j1+2,j1+2) = edde
                subh(j1+3,j1+3) = edde
                subh(j1+4,j1+4) = edde
                j1 = j1 + 5
              case (0)
                subh(j1,j1) = es+de
                j1 = j1 + 1
              end select
          enddo
          
        else
              do i1 = 1,nstt1
                subh(i1,i1) = 1.0_dp
              enddo
        endif
         
      
      endif

!      DO J1=1,5
!         WRITE(79,'(5F10.6)') SUBH(J1,1),SUBH(J1,2),SUBH(J1,3),
!     +            SUBH(J1,4),SUBH(J1,5)
!      ENDDO

      end

