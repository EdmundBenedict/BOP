 
       function entropy(ef)
          use mod_precision
          use mod_all_scalar
          use mod_const
          use mod_srt
          use ab_io
          use mod_g0n, only : g00
!
!    This is a function to evaluate the entropic contribution to
!     the free energy of the electrons.
!

      implicit none
      real(dp) :: entropy


      include "Include/Atom.array"

      complex(dp) :: zi

      real(dp) :: nu,etau
      real(dp) :: entden,ef
      real(dp) :: prod,sum
      real(dp) :: taur,nur,nurs
      real(dp) :: e0,e1,i0,x0r
      real(dp) :: simp,x0,x1,x,dx

      integer :: i, la, r, s, m

      integer, parameter :: mx = 1111

      real(dp) :: y(mx)
      
      
      complex(dp) :: zp,zfac,zep,w0
      real(dp) :: phase
      
      integer :: p

!    Evaluate the entropy. for the atom to which the recursion coefficients are associated
      
      entropy = 0.0_dp

      if (term == 1) then

         if (chi_meth == 1) then

            do la = 1,nchain

               nu = -cmplx((ef-lainf(la))/(2*lbinf(la)), kind=dp)
               etau = -cmplx(kt/(2*lbinf(la)), kind=dp)

               e0 = nu - fdx0*etau
               e1 = nu + fdx0*etau

               prod = 0.5_dp
               do i = 1,lchain(la)
                  prod = prod*(brec(i,la)/(2*lbinf(la)))
               enddo
               prod = -wt(la)*kt*prod*prod/pi

               call get_ai_un(e0,mrec,root(1,la),f1(1,la),f2(1,la), & 
     &                  forder(la),2*mfdpol+2,srtin,srtrn,srtun0)

               call get_ai_un(e1,mrec,root(1,la),f1(1,la),f2(1,la), & 
     &                  forder(la),2*mfdpol+2,srtin,srtrn,srtun1)

               i0 = srtun1(0)-srtun0(0)
               taur = etau*etau
               nur = nu*nu
               x0r = fdx0*fdx0
               do r = 0,mfdpol,1
                  nurs = nur
                  sum = 0.0_dp
                  do s = 0,2*r+2
                     sum = sum + triang(2*r+3,s+1)*nurs*(srtun1(s)-srtun0(s))
                     nurs = -nurs/nu
                  enddo
                  entropy = entropy + (2*r+1)*fdpol(r)*(sum/taur-x0r*i0)*prod/real(2*r+2, dp)
                  taur = taur*etau*etau
                  nur = nur*nu*nu
                  x0r = x0r*fdx0*fdx0
               enddo
            enddo

         else

            do la = 1,nchain
            
! !             this is wrong, what should happen is to find the residuals of sigma(f(z))) instead of only f(z)
!             f(z) is dealth with in the bratkovski horsfield 1995/6 paper appendix a. 
!             if residuals are found for ln(z) = n*(z**(1/n)-1) for n->inf and sigma it may be possible to sync 
!             the bond and entropic energy
! 
!                m     = nint(mfac*0.25_dp*(ef-lainf(la) +2*lbinf(la))/kt)
!    !            m     = nint(mfac*0.25_dp*1.8_dp/kt)
!                w0    = kt*(2*m)
!                phase = pi/real(2*m, dp)
!                zp    = cmplx(cos(phase),sin(phase), kind=dp)
!                zfac  = zp*zp
!                
!                do p = 0, m-1
!                   zep = ef + w0*(zp-1)
!                   entropy = entropy &
!                                & + real(zp*zentden(zep)*g00(zep,arec(0:lchain(la),la), &
!                                                  & brec(0:lchain(la)+1,la), & 
!                                                  & lchain(la),lainf(la),lbinf(la)))
!                   zp  = zp*zfac                                    
!                end do
!                
!                entropy = -4*kt*wt(la)*entropy
! !             
            
!             
               x0 = max((lainf(la)-2*lbinf(la)-ef)/kt,-10.0_dp)
               x1 = min((lainf(la)+2*lbinf(la)-ef)/kt, 10.0_dp)
!                m = min(mx,max(11,nrec*nint(50*kt)))
               m = min(mx,max(111,nrec*nint(500*kt)))
               m = 2*(m/2)+1
               x = x0
               dx = (x1-x0)/real(m-1, dp)
               do i = 1,m
                  zi = cmplx(ef + kt*x, kind=dp)
                  y(i) = aimag(g00(zi,arec(0:lchain(la),la), &
                                    & brec(0:lchain(la)+1,la),lchain(la), & 
                                    & lainf(la),lbinf(la)))*entden(kt*x,kt)
                  x = x0 + i*dx
               enddo
               entropy = entropy + wt(la)*simp(dx,y,m)
            enddo
            entropy = 2*kt*kt*entropy/pi
         endif

      elseif ((term == 2).or.(term == 3)) then

         do la = 1,nchain
            do i = 1,lchain(la)+1
               entropy = entropy + eigvec(1,i,la)*eigvec(1,i,la)*wt(la)*entden(diag(i,la)-ef,kt)
            enddo
         enddo         
         entropy = -2*kt*entropy
      endif

!       contains
!       
!          function zentden(x)
!             use mod_precision
! 
!             ! evluate complex entropy density function.
! 
!             implicit none
!             complex(dp), intent(in) :: x
!             complex(dp) :: zentden
!             complex(dp) :: n,m
! 
!             n = 1.0_dp/(1.0_dp+exp(x))
!             m = 1.0_dp - n
! 
!             zentden = n*log(n) + m*log(m)
!          end function zentden
      end function entropy

