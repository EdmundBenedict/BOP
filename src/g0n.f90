 
    module mod_g0n
    
    use mod_precision
    implicit none
    
    contains
 
     pure function g00(z,a,b,nrec,ainf,binf)


!
!    This is the Green's function.
!
      
      integer, intent(in) :: nrec
      real(dp), intent(in) :: ainf,binf,a(0:nrec), b(0:nrec+1)
      complex(dp), intent(in) :: z
      
      
      complex(dp) :: g00

      complex(dp) :: v,w,t
      real(dp) :: binf2
      integer :: i

!    Evaluate the termination function.
!
      if (binf /= 0.0_dp) then
         v = z-ainf
         binf2 = binf+binf
         w = sqrt((v - binf2)*(v + binf2))
         if (imag(w) > 0.0_dp) then
            t = 0.5_dp*(v - w)
         else
            t = 0.5_dp*(v + w)
         endif
!          t = 0.5_dp*(v - sign(1.0_dp,imag(w))*w)
      else         
        t = 0.0_dp
      endif

!
!    Evaluate the continued fraction.
!

      g00 = 1.0_dp/(z-(a(nrec)+t))
      do i = nrec,1,-1
        g00 = 1.0_dp/(z-(a(i-1)+b(i)*b(i)*g00))
      enddo
      


      end function g00

!------------------------------------------------------------------------------

      pure subroutine getg0n(z,a,b,g0n,nmax,nrec,ainf,binf)


!
!    This is a Green's function.
!

      integer, intent(in) :: nmax,nrec
      real(dp), intent(in) :: a(0:nrec), b(0:nrec+1), ainf,binf
      complex(dp), intent(in) :: z
      complex(dp), intent(out) :: g0n(0:nmax)
      
      integer :: i
      real(8) :: inv_binf
      complex(dp) :: x

      



      
!
!    Evaluate the Green's function.
!

      g0n(0) = g00(z,a,b,nrec,ainf,binf)
      g0n(1) = ((z-a(0))*g0n(0)-1.0_dp)/b(1)
      
      if (nmax <= nrec) then
         do i = 2,nmax
            g0n(i) = ((z-a(i-1))*g0n(i-1)-b(i-1)*g0n(i-2))/b(i)         
         enddo
      else
         do i = 2,nrec
            g0n(i) = ((z-a(i-1))*g0n(i-1)-b(i-1)*g0n(i-2))/b(i)
         enddo

         inv_binf = 1.0_dp/binf

         g0n(nrec+1) = ((z-a(nrec))*g0n(nrec)-b(nrec)*g0n(nrec-1))*inv_binf
     
         x = (z-ainf)*inv_binf
         
         do i = nrec+2,nmax
            g0n(i) = x*g0n(i-1)-g0n(i-2)
         enddo
         
      endif

      end subroutine getg0n

    end module mod_g0n