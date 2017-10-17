 
!****************************************************************************
!                                                                           *
!                              Alexey Girshick                              *
!                                                                           *
!              Department of Materials Science and Engineering,             *
!                      University of Pennsylvania, USA                      *
!                                   and                                     *
!                         Department of Materials                           *
!                        University of Oxford,  UK                          *
!                                                                           *
!                                June, 1995                                 *
!                                                                           *
!****************************************************************************
!                                                                           *
!   This subroutine calculates energy and forces on atoms for the           *
!   Finnis-Sinclair potential scheme. Potential reading part and all the    *
!   potentials functions are located in the file "pot.f". This subroutine   *
!   can be used with any functional form of pair potentials in the          *
!   Finnis-Sinclair scheme, while "pot.f" is specific for the cubic         *
!   splines.                                                                *
!                                                                           *
!****************************************************************************
!     
!     Modified to include screening function for many-body repulsive
!     term by M. Cawkwell, January 2004 whilst visiting DGP and DNM
!     at the Department of Materials, University of Oxford.
!

      subroutine finnis ( )
          use mod_precision
          use mod_all_scalar
          use mod_const
          use mod_clock
          use mod_conf, only : rtc

      implicit none

!
!    Include constants.
!

!      include "Include/ALL.const"

!
!    Include scalars.
!

!      include "Include/ALL.scalar"

!
!    Include arrays.
!

!      include "Include/Atom.array"
      include "Include/Force.array"
      include "Include/PosVel.array"
      include "Include/NebList.array"
      include "Include/NebListL.array"

!
!    Common blocks introduced by A.Girshick.
!    
 
      include "Include/ag.conn"
!      include 'Include/Envir.array'

!
!    Declare the variables.
!

      integer  i, j, k, p, l, nsb, mjc
      integer  ja, la, ipick
      integer  ipjl, ippl
      integer laptrl(mxtotnd)
      integer lbptrl(mxtotnd)
      integer stat
      integer lpcnt

      real(dp) ::  rx, ry, rz, rpj, rh , mpw, getlamzero
      real(dp) ::  xpl , ypl , zpl, rpl
      real(dp) ::  xjl , yjl , zjl, rjl
      real(dp) ::  rcutenv, rcenv, r1env, fenvspl
      real(dp) ::  getrcut, getrc, getr1, dchi
      real(dp) ::  e_atomp, rho_sum, del_sum ,fab( 3 )
      real(dp) ::  dvr, dpr
      real(dp) ::  expot, dexpot, venv, gsppot
      real(dp) ::  chi,pom_del, pom_ddel, pom_delbig
      real(dp) ::  rho( mxtotnd )
      real(dp) ::  dellam ( mxtotnd ), ddellam( mxtotnd)
      real(dp) ::  derlamx ( mxtotnd ), derlamy( mxtotnd )
      real(dp) ::  derlamz ( mxtotnd )
      real(dp) ::  env, dprx, dpry, dprz
      real(dp) ::  pomx, pomy, pomz
      real(dp) ::  ejl, ftbx, ftby, ftbz, test
      real(dp) ::  delrmax
      
      integer(8) :: c1, c2

!       print *, 'nne', nne
      epair = 0.0_dp
      eenv = 0.0_dp
      eclas = 0.0_dp
!       call system_clock(c1)
      if( env_include  ==  1 )  then

         do  p = 1, nne
            del_sum = 0.0_dp
            ja = aptrl( p )
            do while ( bptrl( ja )  /=  eol )
               
               j = bptrl( ja )
               
               if ( j  /=  p ) then

                  ipick = ind( j )

                  rx = ad( 1, p ) - ad( 1, j )
                  ry = ad( 2, p ) - ad( 2, j )
                  rz = ad( 3, p ) - ad( 3, j )
                  rpj = sqrt( rx*rx + ry*ry + rz*rz )
                  
                  del_sum = del_sum + chi( rpj, ipick)
                  
               endif
               ja = ja + 1
            enddo
            
            dellam( p ) = exp( (1.0_dp / mpw(ind(p))) * log( del_sum ) )
            
         enddo
         
!        Now fill out the rest of DELLAM array using the periodicity.
                
         do  j = nne+1, totnd
            dellam( j )  = dellam( map( j ) )
         enddo
      endif


!       call system_clock(c2)
!       print *,'dl',real(c2-c1,dp)/real(cr,dp)


!
!  Main cycle - over all atoms.
!       write(85,*) 'n it jt bt rx ry rz r'

!       call system_clock(c1)
      do p = 1, nd     
         e_atomp = 0.d0         
         ja = aptrl( p )
         j = bptrl( ja )
         print *, p, ind(p), rtc % ham_conf % a(ind(p)) % symb
         do while ( j  /=  eol )
            
!             print *, '    ', j
            if ( j  /=  p ) then

               

               ipick = ind( p ) + ind( j )

               rx = ad( 1, p ) - ad( 1, j )
               ry = ad( 2, p ) - ad( 2, j )
               rz = ad( 3, p ) - ad( 3, j )
               rpj = sqrt( rx*rx + ry*ry + rz*rz )

!                print *, p,j, ind(p),ind(j), ipick, rtc % ham_conf % b(ipick) % name
            
!                write(85, '(4(x,i3), 4(x,f20.14))') &
!                 & p, ind(p), ind(j), ipick, rx ,ry, rz, rpj 
!                 flush(85)
!                 
               if (env_include == 1) then
                  pom_del = 0.5_dp * ( dellam(p) + dellam(j) )
                  eenv = eenv + 0.5_dp * venv(rpj,ipick,pom_del)
               end if
               
               if (vpair_include == 1) then
!                  e_atomp = e_atomp + 0.5_dp * expot( rpj, ipick)
                 epair = epair + 0.5_dp * gsppot( rpj, ipick)
               end if
            endif

            ja = ja + 1
            j = bptrl( ja )
            
         enddo

         
!          write(6,'("EPAIR = ",f12.6)') epair
      enddo
      
!       call system_clock(c2)
!       print *,'rep',real(c2-c1,dp)/real(cr,dp)

      eclas = epair + eenv


!       
!       
!       
!     if (forces) then
! 
!       dprx = 0.0_dp
!       dpry = 0.0_dp
!       dprz = 0.0_dp
!       
!       ftbx = 0.0_dp
!       ftby = 0.0_dp
!       ftbz = 0.0_dp
! 
! !end edit
! 
! 
! 
!       epair = 0.0_dp
!       fpp(:,:nd) = 0.0_dp
!       rho = 0.0_dp
!       
! 
! 
! !  Get RHO for all atoms.
!  
!  
! !Dimitar Pashov comment
! !       IF( EMB_INCLUDE .EQ. 1 )  THEN
! ! 
! ! *        First get RHO for all atoms in the original block.
! ! 
! !          DO  P = 1, NNE
! !             RHO_SUM = 0.D0
! !             JA = APTRL( P )
! ! 
! !             DO WHILE ( BPTRL( JA ) .NE. EOL )
! !                 J = BPTRL( JA )
! !                 IF ( J .NE. P ) THEN
! !                     IPICK = IND( P ) + IND( J )
! !                     RX = AD( 1, P ) - AD( 1, J )
! !                     RY = AD( 2, P ) - AD( 2, J )
! !                     RZ = AD( 3, P ) - AD( 3, J )               
! !                     RPJ = DSQRT( RX*RX + RY*RY + RZ*RZ )
! !                     RHO_SUM = RHO_SUM + PHI( RPJ, IPICK )
! !                 ENDIF
! !                 JA = JA + 1
! !             ENDDO     
! ! 
! !             RHO( P ) = DSQRT( RHO_SUM )
! !          ENDDO
! ! 
! ! *        Now fill out the rest of RHO array using the periodicity.
! ! 
! !          DO  J = NNE+1, TOTND
! !              RHO( J ) = RHO( MAP( J ) )
! !          ENDDO
! !       ENDIF
! !end comment
! 
! !
! !    Additional environment dependent term
! !     Introduced by Stefan Znam
! !
!       
! 
!       if( env_include  ==  1 )  then
! 
! !        First get DEL_LAM for all atoms in the original block.
! 
! ! !Dimitar Pashov edit
! !          write(6,*) 'NNE = ', nne
! !          write(6,*) 'AD = ', ad
! !          do while (aptrl(i) .ne. eol)
! !             
! !         write(6,*) 'BPTRL = ', bptrl
! !         
! ! !end edit Dimitar Pashov
! 
!          do  p = 1, nne
! !Dimitar Pashov
! !         write(6,*) 'P = ', P
! !end edit
! 
!             del_sum = 0.d0
!             derlamx( p ) = 0.d0
!             derlamz( p ) = 0.d0
!             derlamy( p ) = 0.d0
! 
!             ja = aptrl( p )
! !Dimitar Pashov
! !         write(6,*) 'JA = APTRL( P ) = ', JA
! !end edit
! !             write(6,*) BPTRL
! 
!             do while ( bptrl( ja )  /=  eol )
!                j = bptrl( ja )
!                
! !Dimitar Pashov
! !                write(6,*) 'JA = ', JA, ';   J = BPTRL( JA ) = ', J
! !                write(6,*) 'P, J = ', P, J
! !end edit
!                
!                
!                if ( j  /=  p ) then
! !Dimitar Pashov
! !                   IPICK = IND( J )
!                   ipick = ind(j) + ind(p)
! !Dimitar Pashov
!                   
!                   rx = ad( 1, p ) - ad( 1, j )
!                   ry = ad( 2, p ) - ad( 2, j )
!                   rz = ad( 3, p ) - ad( 3, j )
!                   rpj = sqrt( rx*rx + ry*ry + rz*rz )
! !Dimitar Pashov
! !                write(6,*) 'R(P,J) = ', rpj
! !end edit
!                   del_sum = del_sum + chi( rpj, ipick)
!                   derlamx( p ) = derlamx( p ) + dchi( rpj, ipick) * rx / rpj
!                   derlamy( p ) = derlamy( p ) + dchi( rpj, ipick) * ry / rpj
!                   derlamz( p ) = derlamz( p ) + dchi( rpj, ipick) * rz / rpj
!                endif
!                ja = ja + 1
!             enddo
! !  
! !             DELLAM( P ) = DEXP( (1.0D0 / MPW(IPICK)) * DLOG( DEL_SUM ) )
! !             DDELLAM( P ) = 1.0D0/MPW(IPICK) * DEXP( (1.D0/MPW(IPICK) &
! !      &                      - 1.0D0 ) * DLOG(DEL_SUM) )
!             
!             dellam( p ) = exp( (1.0_dp / mpw(ind(p))) * log( del_sum ) )
!             ddellam( p ) = 1.0_dp/mpw(ind(p)) * exp( (1.0_dp/mpw(ind(p)) - 1.0_dp ) * log(del_sum) )
!             
! 
! !Dimitar Pashov
! !             write(6,'("delsum(",i0,") = ",f12.6)') P, del_sum
! !             write(6,'("DELLAM(",i0,") = ",f12.6)') P, DELLAM(P)
! !             write(6,'(///)')
! !Dimitar Pashov
! 
! 
!          enddo
! 
! !        Now fill out the rest of DELLAM array using the periodicity.
!                 
!          do  j = nne+1, totnd
!             dellam( j )  = dellam( map( j ) )
!             ddellam( j ) = ddellam( map( j ) )
!             derlamx( j ) = derlamx( map( j ) )
!             derlamy( j ) = derlamy( map( j ) )
!             derlamz( j ) = derlamz( map( j ) )
!          enddo
!       endif
! !
! !  Main cycle - over all atoms.
! 
!       do p = 1, nd     
!          e_atomp = 0.d0
!          
!          ja = aptrl( p )
! 
!          do while ( bptrl( ja )  /=  eol )
!             j = bptrl( ja )
! 
! !     Add the contributions from the P - J interaction.
!  
!             if ( j  /=  p ) then
! !
!                ipick = ind( p ) + ind( j )
!                rx = ad( 1, p ) - ad( 1, j )
!                ry = ad( 2, p ) - ad( 2, j )
!                rz = ad( 3, p ) - ad( 3, j )
!                rpj = sqrt( rx*rx + ry*ry + rz*rz )
!                 
!                
! 
!                select case(env_include)
!                   
!                   case(0)
!                   
!                   dpr = 0.d0
!                   
!                   case(1)
!                   
!                   rcutenv = getrcut(ipick)
!                   
!                   pom_del = 0.5_dp * ( dellam(p) + dellam(j) )
! 
!                   pom_ddel = 0.5_dp * dchi(rpj,ipick) * ddellam(j) 
! 
!                   pomx = 0.5_dp * derlamx(p) * ddellam(p) 
!                   pomy = 0.5_dp * derlamy(p) * ddellam(p) 
!                   pomz = 0.5_dp * derlamz(p) * ddellam(p) 
! 
!                   env = venv(rpj,ipick,pom_del)
! !                   write(6,'("ENV = ", f12.6)') env
! !                   write(6,'("RPJ = ", f12.6)') rpj
! !                   write(6,'("pom_del = ", f12.6)') pom_del
!                   e_atomp = e_atomp + 0.5_dp * env
!                   
!                   rcenv = getrc(ipick)   
!                                         
!                   pom_delbig = getlamzero(ipick) + 0.5*(dellam(p)+dellam(j))
! 
!                   rh = 1.d0/rpj + pom_delbig + (rpj-2*rcenv)*pom_ddel
! 
!                   dpr = env * rh
! 
!                   dprx = env * (rpj-2*rcenv)*pomx
!                   dpry = env * (rpj-2*rcenv)*pomy
!                   dprz = env * (rpj-2*rcenv)*pomz
! 
! !     Adding the three-body force contributions of the envir. term
! 
!                   ftbx = 0.0_dp
!                   ftby = 0.0_dp
!                   ftbz = 0.0_dp
! 
!                   la = aptrl( p )
! !
!                   lpcnt = 0
!                   delrmax = 1000.0_dp
!                   omsenvij = 0.0_dp
!                   do mjc = 1,3
!                      senvfc(mjc) = 0.0_dp
!                   enddo
! !
! !     Loop over atom k
! !
!                   do while ( bptrl( la )  /=  eol )
!                      l = bptrl( la )
!                      if ( (l  /=  p) .and. (l  /=  j) ) then
! 
!                         xpl = ad( 1, p ) - ad( 1, l )
!                         ypl = ad( 2, p ) - ad( 2, l )
!                         zpl = ad( 3, p ) - ad( 3, l )
!                         rpl = sqrt( xpl*xpl + ypl*ypl + zpl*zpl )
!      
!                         xjl = ad( 1, j ) - ad( 1, l )
!                         yjl = ad( 2, j ) - ad( 2, l )
!                         zjl = ad( 3, j ) - ad( 3, l )
!                         rjl = sqrt( xjl*xjl + yjl*yjl + zjl*zjl )
! 
!                         pom_del = 0.5_dp * ( dellam(j) + dellam(l) )
!                         
!                         ipjl = ind( j ) + ind( l )
!                         ippl = ind(p) + ind(l)
!                         
!                         ejl = 0.25_dp * venv (rjl, ipjl , pom_del ) * ( rjl - 2 * getrc ( ipjl ) )
!                         
!                         pomx =  ddellam( j ) * dchi ( rpj, ipick ) * rx/rpj + xpl/rpl * ddellam( l ) * dchi ( rpl, ippl )
! 
!                         pomy =  ddellam( j ) * dchi ( rpj, ipick ) * ry/rpj + ypl/rpl * ddellam( l ) * dchi ( rpl, ippl )   
! 
!                         pomz =  ddellam( j ) * dchi ( rpj, ipick ) * rz/rpj + zpl/rpl * ddellam( l ) * dchi ( rpl, ippl )
!      
!                         ftbx = ftbx + pomx * ejl
!                         ftby = ftby + pomy * ejl
!                         ftbz = ftbz + pomz * ejl
!                       
!                      endif
! 
!                      la = la + 1
!                   enddo
!                   
!                   case default
!                   write(0,'("Error in env_include. Expected 0 or 1, got ", i0, ".")') env_include
!                end select
! 
!                select case(vpair_include)
!                     case(0)
!                         dvr = 0.0_dp
!                     case(1)
! !                         DVR = DVEE( RPJ, IPICK )
! !                         E_ATOMP = E_ATOMP + 0.5D0 * VEE( RPJ, IPICK)
!                         dvr = dexpot( rpj, ipick )
!                         e_atomp = e_atomp + 0.5_dp * expot( rpj, ipick)
!                     case default
!                         write(0,'("Error in vpair_include. Expected 0 or 1, got ", i0, ".")') vpair_include
!                         call panic()
!                end select
! 
!                fab( 1 ) =  dprx + ftbx + (rx / rpj)*(dpr - dvr )
!                fab( 2 ) =  dpry + ftby + (ry / rpj)*(dpr - dvr )
!                fab( 3 ) =  dprz + ftbz + (rz / rpj)*(dpr - dvr )
! 
!                fpp( 1, p ) = fpp( 1, p ) + fab( 1 )
!                fpp( 2, p ) = fpp( 2, p ) + fab( 2 )
!                fpp( 3, p ) = fpp( 3, p ) + fab( 3 )
! 
!             endif
!             ja = ja + 1
!          enddo
! 
!          epair = epair + e_atomp
! !          write(6,'("EPAIR = ",f12.6)') epair
!       enddo
!       
!       end if
      
      end
