 
       function bseatm(es,ep,ed,nl,llist,zc,z)
          use mod_precision


!
!    This is a routine to evaluate the bandstructure energy for
!     an isolated atom.
!

      implicit none
      real(dp) :: bseatm

!
!    Declare the simple variables.
!

      real(dp) :: es,ep,ed,zc
      real(dp) :: temp,sum,occ

      integer nl,z,i

!
!    Declare the arrays.
!

      real(dp) :: e(3)
      real(dp) :: n(3)

      integer llist(nl)

!
!    Store onsite energies and maximum occupancies.
!
      sum = 0.0_dp
      do i = 1,nl,1
         if (llist(i) == 0) then
            e(i) = es
            n(i) = 2.0_dp
         elseif (llist(i) == 1) then
            e(i) = ep
            n(i) = 6.0_dp
         elseif (llist(i) == 2) then
            e(i) = ed
            n(i) = 10.0_dp
         endif
         sum = sum + n(i)
      enddo
      
      


      if (sum < zc) then
         write(6,'(''Too many electrons for atom type '',I3)') z
         write(6,'(''ZC = '',G12.5)') zc
         write(6,'(''Number of states = '',G12.5)') sum
         call panic()
      endif

!
!    Sort the energies.
!

      if (nl == 2) then
         if (e(1) > e(2)) then
            temp = e(1)
            e(1) = e(2)
            e(2) = temp
            temp = n(1)
            n(1) = n(2)
            n(2) = temp
         endif
      elseif (nl == 3) then
         if (e(1) > e(2)) then
            temp = e(1)
            e(1) = e(2)
            e(2) = temp
            temp = n(1)
            n(1) = n(2)
            n(2) = temp
         endif
         if (e(2) > e(3)) then
            temp = e(3)
            e(3) = e(2)
            e(2) = temp
            temp = n(3)
            n(3) = n(2)
            n(2) = temp
            if (e(1) > e(2)) then
               temp = e(1)
               e(1) = e(2)
               e(2) = temp
               temp = n(1)
               n(1) = n(2)
               n(2) = temp
            endif
         endif
      endif

!
!    Calculate the total energy.
!

      bseatm = 0.0_dp
      sum = zc
      do i = 1,nl,1
         occ = min(n(i),sum)
         sum = sum - occ
         bseatm = bseatm + occ*e(i)
!         print *, 'i,occ,bseatm:',i,occ,bseatm
      enddo

      end

