 
       function numelsrt_a(ef,la)
          use mod_precision
          use mod_all_scalar
          use mod_const
          use mod_srt
          use ab_io


!
!    This is a function to evaluate the number of electrons for the
!     square root terminator for a given state.
!

      implicit none
      
      real(dp), intent(in) :: ef
      integer, intent(in) :: la
      
      real(dp) :: numelsrt_a
      
      real(dp) :: prod
      
      integer i

!
!    Evaluate the number of electrons.
!
         prod = 0.5_dp
         do i = 1,lchain(la)
            prod = prod*(brec(i,la)/(2*lbinf(la)))
         enddo
         prod = -prod*prod/pi
         call get_srtun(0,ef,la)
         numelsrt_a = srtun(0)*prod

      end

!-----------------------------------------------------------------

       function cnel(ef,la)
          use mod_precision
          use mod_all_scalar
          use mod_const
          use ab_io
          use mod_g0n, only : g00
!           use topologia, only : iproc

!     AMB: complex plane integration for number of electrons
!
!     This is a function to evaluate the number of electrons for the
!     square root terminator for a given state.


      implicit none
      
      
      real(dp), intent(in) :: ef
      integer, intent(in) :: la
      
      
      real(dp) :: cnel

      complex(dp) :: zp,zfac,zep
      real(dp) :: phase,w0
      integer :: m,p


!    Evaluate the number of electrons.


       cnel  = 0.0_dp
!      M     = NINT(MFAC*0.25D0*(EF-LAINF(LA)+2.0D0*LBINF(LA))/KT)
       m     = nint(mfac*0.25_dp*(ef-lainf(la)+2*lbinf(la))/kt)
!       print *,ef, lainf(la), lbinf(la), (ef-lainf(la)+2.0_dp*lbinf(la))
!       m     = nint(mfac*0.25_dp*1.8_dp/kt)

!       write(100+iproc,*) ef, lainf(la), lbinf(la)
!       write(100+iproc,*) '      ',ef, arec(0:lchain(la),la), brec(0:lchain(la)+1,la)
      
      w0    = kt*(2*m)
      phase = pi/real(2*m, dp)
      zp    = cmplx(cos(phase),sin(phase), kind=dp)
      zfac  = zp*zp
!       print *, la, lchain(la)
      do p = 0, m-1
         zep  = ef+w0*(zp-1)
         cnel = cnel + real(zp*g00(zep,arec(0:lchain(la),la),brec(0:lchain(la)+1,la), &
                                   & lchain(la),lainf(la),lbinf(la)))
         zp   = zfac*zp
      enddo
      
      cnel = 4*kt*cnel

      end function cnel

