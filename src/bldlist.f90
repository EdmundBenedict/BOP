 
      subroutine bldlist( n, rlim, a, dlim, ad, ntot, map, z, lmapi, nd, mapi)
          use mod_precision
          use mod_const

!
!    This is a routine to build the complete list of atoms.
!

      implicit none
      
      logical, intent(in) :: lmapi ! whether we shalt fill mapi
      
!*** The number of unit cells to include in each dimension.
      integer, intent(in) ::  rlim(3), n, nd ! n: total in ucell, nd: dynamic

      real(dp), intent(in) :: a(3,3), dlim(2,3)

      
!*** Basis set in Cartesian coordinates.
      real(dp), intent(inout) :: ad(3,mxtotnd)
      integer, intent(inout), dimension(mxtotnd) :: map, z, mapi
      
!*** MAP is an array linking the central cell atoms to their periodically
!    repeated neighbors.
      integer, intent(out) :: ntot
      
!       include "Include/ag.conn"


!
!    Declare the simple variables.
!

      integer ia
      integer i1,i2,i3
      integer ns
      integer i,j

!
!    Declare the local arrays.
!

!*** A is the array of primitive translation vectors, and B is the inverse.
      real(dp) :: b(3,3)
!*** The position of one atom.
      real(dp) :: dr(3)
      real(dp) :: pos(3)
      real(dp) :: ld(3)

      real(dp) :: zstr(3)

!
!    Build list of atoms.
!


      ns = n

      b = a
      call inv3x3(b)

      do i1 = -rlim(1),rlim(1)
         do i2 = -rlim(2),rlim(2)
            do i3 = -rlim(3),rlim(3)
               if ((i1 /= 0).or.(i2 /= 0).or.(i3 /= 0)) then
                  dr = real(i1, dp)*a(1, 1:3)+real(i2, dp)*a(2, 1:3)+real(i3, dp)*a(3, 1:3)
                  do ia = 1,n
                     pos = dr + ad(1:3, ia)
                     call mul3x3(b,pos,ld)
                     if ( all((rlim(1:3) == 0).or. & 
     &                    ((ld(1:3) >= dlim(1,1:3)).and. & 
     &                     (ld(1:3) <= dlim(2,1:3)))) ) then
                        ns = ns + 1
                        ad(1:3, ns) = pos                        
                        map(ns) = map(ia)
                        z(ns) = z(ia)
                        if (lmapi .and. ia > nd) then
                           mapi(ns) = ia - nd ! ad hoc fix for the magnetic case taken from Matous
!                            print *, 'ns, mapi', ns, mapi(ns), map(ns)
                        end if   
                        if (ns > mxtotnd) then
                           write(6,'(''Too many atoms. Increase MXTOTND.'')')
                           call panic()
                        endif
                     endif
                  enddo
               endif
            enddo
         enddo
      enddo

!     ntot = n*(2*rlim(1)+1)*(2*rlim(2)+1)*(2*rlim(3)+1)
!     write(6,'(''Number of atoms before pruning = '',I6)') ntot
!     write(6,'(''Number of atoms after pruning = '',I6)') ns
    ntot = ns

! !
! !    Add inert atoms to the end of the list.
! !
! 
!       do i1 = -rlim(1),rlim(1)
!          do i2 = -rlim(2),rlim(2)
!             do i3 = -rlim(3),rlim(3)
!                dr(:) = real(i1, dp)*a(1,:)+real(i2, dp)*a(2,:)+real(i3, dp)*a(3,:)
!                do ia = 1,ninert
!                   pos(:) = dr(:)+adinert(:,ia)
!                   call mul3x3(b,pos,ld)
!                   if (all((rlim(:) == 0).or. & 
!      &                 ((ld(:) >= dlim(1,:)).and. & 
!      &                  (ld(:) <= dlim(2,:))))) then
!                      ns = ns + 1
!                      ad(:,ns) = pos(:)
!                      map(ns) = 0
!                      z(ns) = zinert(ia)
! 
!                      if (ns > mxtotnd) then
!                         write(6,'(''Too many atoms. Increase MXTOTND.'')')
!                         call panic()
!                      endif
!                   endif
!                enddo
!             enddo
!          enddo
!       enddo
! 
!       ntot = ns


      end subroutine bldlist

