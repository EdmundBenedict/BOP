
   subroutine gethu(u,nu,hu,nhu,nma)
   
   use mod_precision
   use mod_const
   use mod_ham
   
!    This is a routine to multiply the Hamiltonian into a vector.
   
   implicit none

   include "Include/Atom.array"
   include "Include/NebList.array"

   integer, intent(in) :: nu,nhu,nma
   real(dp), intent(in) :: u(mxnstat,mxcls,nma)
   real(dp), intent(inout) :: hu(mxnstat,mxcls,nma)
   
   integer :: nsttb,nstta,ia,ib,ja,ja0,ic,lb,i,dj,lla
!       real(dp) :: dot
!       character(len=1) :: tmpchr
   
   hu(1:mxnstat,1:nhu,1:nma) = 0.0_dp 

!       write(601,'(/,"s:     i     ia     dj     ib     ic",/,37("-"))')
   do i = 1,nu
      ia = cluster(i)
      nstta = nstt(z(ia))
      ja = aptr(ia)
      ja0 = ja -1
      ib = bptr(ja)
      do while (ib /= eol)
!             tmpchr = 'o'
         ic = decipher(ib)
         if (ic /= 0) then
            dj = ja-ja0
            nsttb = nstt(z(ib))

            do lla = 1, nma
               do lb = 1, nsttb
! if issues occur try putting the limit up to nstta
! mxnstat is used because it constant (f90 parameter) and allows the compiler to vectorise the dot_product 
!                   hu(lb,ic,lla) = hu(lb,ic,lla) + dot_product(h(:,lb,dj,i),u(:,i,lla))
                  hu(lb,ic,lla) = hu(lb,ic,lla) + sum(h(:,lb,dj,i)*u(:,i,lla))
!                      hu(lb,ic,lla) = hu(lb,ic,lla) + dot(h(1,lb,dj,i),u(1,i,lla),mxnstat)
               end do
            end do
!                 tmpchr = 'i'
         endif
!             write(601,'(a,5(x,i6))') tmpchr,i,ia,ja-ja0,ib,ic
         ja = ja + 1
         ib = bptr(ja)
      enddo
   enddo
   
   end subroutine gethu
