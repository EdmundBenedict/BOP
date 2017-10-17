 
      subroutine envscr(i, j, k, rij, rik, rjk, lpcnt)
          use mod_precision


          use mod_const
!
!     This subroutine is called from FINNIS() and calculates the
!     screening function S_pj for the many-body repulsive term. This is
!     switched on using ENVSCRON in your input file fort.8. 
!     
!     The form for the screening functions is as follows:
!
!     S_ij^k = O_ik * O_jk / O_ij
!
!     where O_ij = A_ij * exp (-mu_ij * R_ij)
!
!     The parameters A should be equal to 1 but have been included
!     here to give a little more flexibility during fitting if this
!     is needed. TRY WITH ALL THE A'S EQUAL TO 1 FIRST!!! - this is 
!     it should be 'in theory'. The A's and mu's are read from fort.12
!     fort.13 and fort.14.
!
!     The [1 - S_ij]'s that will be used to screen the many-body 
!     repulsive term can be built up using either the sum or product
!     (Baskes) methods. See which works better for you although the
!     product might be faster and will not suffer such strong 
!     oscillations. Somethings to play with...
!
!     Anyway, that is the general idea...
!
!     Work started on this code in January 2004 with the help of DGP 
!     and DNM at the Department of Materials, University of Oxford.
!
!     M. Cawkwell 14/1/2004. 
!
      implicit none
!
!
!    Include constants.
!

!     include "Include/ALL.const"

!
!    Include scalars.
!

!      include "Include/ALL.scalar"

!
!    Include arrays.
!

!      include "Include/Atom.array"
      include "Include/PosVel.array"
!      include "Include/NebList.array"
!      include "Include/NebListL.array"

!
!    Common blocks introduced by A.Girshick.
!    
 
      include "Include/ag.conn"
!      include 'Include/Envir.array'

!
!    Declare the variables.
!
      integer i, j, k, m, n, ipickij, ipickik, ipickjk
      integer ka, lpcnt 
      real(dp) :: sijk, rij, rjk, rik
      real(dp) :: rcij(3), rcik(3), dsijk(3)
      real(dp) :: aik, ajk, aij, muik, mujk, muij
!
      lpcnt = lpcnt + 1
!
!     I is the label of atom number I, J is the label of atom J etc...
!
!     Determining species of atoms at sites I and J:
!     E.G. Mo = 0, Si = 1 so Mo-Mo: IPICK  = 0, Mo-Si: IPICK = 1
!     Si-Si: IPICK = 2
!
      ipickij = ind(i) + ind(j)
      ipickik = ind(i) + ind(k)
      ipickjk = ind(j) + ind(k)
!
      do m = 1,3
         rcij(m) = ad(m,i) - ad(m,j)
         rcik(m) = ad(m,i) - ad(m,k)
      enddo
!
      if (ipickij  ==  0) then
         aij = aescra
         muij = aescrmu
      else if (ipickij  ==  1) then
         aij = abescra
         muij = abescrmu
      else if (ipickij  ==  2) then
         aij = bescra
         muij = bescrmu
      endif
!
      if (ipickik  ==  0) then
         aik = aescra
         muik = aescrmu
      else if (ipickik  ==  1) then
         aik = abescra
         muik = abescrmu
      else if (ipickik  ==  2) then
         aik = bescra
         muik = bescrmu
      endif            
!     
      if (ipickjk  ==  0) then
         ajk = aescra
         mujk = aescrmu
      else if (ipickjk  ==  1) then
         ajk = abescra
         mujk = abescrmu
      else if (ipickjk  ==  2) then
         ajk = bescra
         mujk = bescrmu
      endif
!     
!     Now calculating S_ij^k for given the constants and the distances
!
!      SIJK = DEXP( -1.0D0 * MUIJ * (RIK+RJK-RIJ))
      sijk = (aik*ajk / aij) * exp( (-1.0_dp * muik*rik) + & 
     &     (-1.0_dp * mujk*rjk) + (muij*rij))
!     
!     Now Sij is constructed from the SIJK's using either the sum
!     or product forms as selected in fort.8
!
!     Now, doing the product form (from Baskes)
!
      if (senvprod  ==  1) then
!     
         if (lpcnt  ==  1) then
            omsenvij = (1.0_dp - sijk)
         else if (lpcnt  >  1) then
!     
            omsenvij = omsenvij * ( 1.0_dp - sijk )
!     
         endif
!     
      endif
!     
!     Now doing the sum form: S_ij = Sum( k .ne. i,j) S_ij^k
!     
      if (senvsum  ==  1) then
!
         if (lpcnt  ==  1) then
            omsenvij = 1.0_dp - sijk
         else if (lpcnt  >  1) then
            omsenvij = omsenvij - sijk
         endif
!
      endif
!  
!     Now preparing the derivative of S_ij^k for the calculation of
!     forces using the screened many-body repulsive term.
!
      do m = 1,3
         dsijk(m) = sijk * (((muij*rcij(m))/rij) - ((muik*rcik(m))/rik))
         senvfc(m) = senvfc(m) + dsijk(m)
      enddo
!
      return
      end
