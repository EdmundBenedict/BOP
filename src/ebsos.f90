 
      function ebsos_a(ia,ef)
          use mod_precision
          use mod_all_scalar
          use mod_const
          use mod_srt
          use ab_io

!    This is a routine to calculate the band structure energy and
!     number of electrons using the onsite expressions.

      implicit none

      integer, intent(in) :: ia
      real(dp), intent(in) :: ef
      
      include "Include/Atom.array"

!    Declare the simple variables.

      real(dp) :: u0,u1
      real(dp) :: prod,ebsos_a
      real(dp) :: getebsnull

      integer la
      integer n

!
!    Calculate the band structure energy and number of electrons.
!
      ebsos_a = 0.0_dp
      do la = 1,nchain
         if (term == 1) then
            prod = 0.5_dp
            do n = 1,lchain(la),1
               prod = prod*(brec(n,la)/(2.0_dp*lbinf(la)))
            enddo
            prod = -prod*prod/pi

            call get_srtun(1,ef,la)
            u0 = srtun(0)
            u1 = srtun(1)

            ebsos_a = ebsos_a + wt(la)*prod*(lainf(la)*u0-2.0_dp*lbinf(la)*u1)
         elseif ((term == 2).or.(term == 3)) then
            ebsos_a = ebsos_a + wt(la)*getebsnull(ef,lchain(la),mrec,diag(1,la),eigvec(1,1,la),kt)
         endif

      enddo

      end

!-------------------------------------------------------------------------

      function ebsos_c(ia,ef)
          use mod_precision
          use mod_all_scalar
          use mod_const
          use mod_srt
          use ab_io
          use mod_g0n, only : g00
!
!     AMB: complex plane integration for onsite estimate of energy.
!     This is a routine to calculate the band structure energy and
!     number of electrons using the onsite expressions.
!

      implicit none

      integer, intent(in) :: ia
      real(dp), intent(in) :: ef
      
      include "Include/Atom.array"

      complex(dp) :: zp,zfac,zep,w0
      real(dp) :: phase,ebsos_c,getebsnull
      
      integer la,m,p

!
!    Evaluate the number of electrons.
!
      ebsos_c = 0.0_dp
      do la = 1,nchain

         if (term == 1) then

            m     = nint(mfac*0.25_dp*(ef-lainf(la) +2*lbinf(la))/kt)
!            m     = nint(mfac*0.25_dp*1.8_dp/kt)
            w0    = kt*(2*m)
            phase = pi/real(2*m, dp)
            zp    = cmplx(cos(phase),sin(phase), kind=dp)
            zfac  = zp*zp
            
            do p = 0, m-1
               zep = ef + w0*(zp-1)
               ebsos_c = ebsos_c + real(zep*zp*g00(zep,arec(0:lchain(la),la), &
                                                     & brec(0:lchain(la)+1,la), & 
                                                     & lchain(la),lainf(la),lbinf(la)))
               zp  = zp*zfac               
            enddo
            
            ebsos_c = 4*kt*wt(la)*ebsos_c 
            
         elseif ((term == 2).or.(term == 3)) then
            ebsos_c = ebsos_c + wt(la)*getebsnull(ef,lchain(la),mrec,diag(1,la),eigvec(1,1,la),kt)
         endif

      enddo

      end

