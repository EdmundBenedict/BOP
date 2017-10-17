 
      function epromavg(ia,ef,isp)
          use mod_precision
          use mod_all_scalar
          use mod_const
          use mod_funptr
          use ab_io, only : diag, eigvec, lchain

!
!    This is a subroutine to evaluate the promotion energy.
!

     
      implicit none

      include "Include/Atom.array"
      
      integer, intent(in) :: ia,isp
      real(dp), intent(in) :: ef
      
      real(dp) :: epromavg
      real(dp) :: nia
      procedure(real(dp)) :: numelnull
      real(dp) :: ela, deia

      integer lla
      integer nstta,nla,la
      integer za

      epromavg = 0.0_dp
      
      za = z(ia)
      
      deia = de(ia)
      if (mag) deia = deia + oppm(isp)*dem(ia)
      
      call states(za,nla,nstta,llista)
      
      do la = 1,nla
         lla = llista(la)
         if (lla == 0) then
            ela = es(za) + deia
         elseif (lla == 1) then
            ela = ep(za) + deia
         elseif (lla == 2) then
            ela = ed(za) + deia
         endif
         if (term == 1) then
            nia = numelsrt(ef,la)
         elseif ((term == 2).or.(term == 3)) then
            nia = numelnull(ef,lchain(la),mrec,diag(1,la), eigvec(1,1,la),kt)
         endif
         
         print *,'1st nia',nia
         
         nia = nia*(2*lla+1)
         
         if (lla  ==  0) eproms = nia*ela
         if (lla  ==  1) epromp = nia*ela
         if (lla  ==  2) epromd = nia*ela
         
         epromavg = epromavg + nia*ela
        print *, 'ia,la,nia,ela', ia, la, nia, ela
      enddo

      end function epromavg

