 
      function epromnoavg(ia,ef,isp)
          use mod_precision
          use mod_all_scalar
          use mod_const
          use mod_funptr
          use ab_io, only : diag, eigvec, lchain


!    This is a subroutine to evaluate the promotion energy.

      implicit none


      include "Include/Atom.array"

      integer, intent(in) :: ia,isp
      real(dp), intent(in) :: ef
      real(dp) :: epromnoavg
      procedure(real(dp)) :: numelnull
      real(dp) :: ela, deia, nia

      integer lla,ila,ima
      integer nstta,nla,la
      integer za


!    Evaluate the on site energies.

      epromnoavg = 0.0_dp
      za = z(ia)
      deia = de(ia)
      if (mag) deia = deia + oppm(isp)*dem(ia)
      call states(za,nla,nstta,llista)
      la = 1
      do ila = 1,nla
         
         lla = llista(ila)
         
         if (lla == 0) then
            ela = es(za) + deia
         elseif (lla == 1) then
            ela = ep(za) + deia
         elseif (lla == 2) then
            ela = ed(za) + deia
         endif
         
         do ima = 1,2*lla+1
            if (term == 1) then
               nia = numelsrt(ef,la)
            elseif ((term == 2).or.(term == 3)) then
               nia = numelnull(ef,lchain(la),mrec, diag(1,la), eigvec(1,1,la),kt)
            endif
            la = la + 1
            epromnoavg = epromnoavg + nia*ela
         enddo
      enddo

      end function epromnoavg

