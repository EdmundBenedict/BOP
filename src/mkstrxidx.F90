subroutine mkstrxidx(force, nbas, ssize, dsize, struxidx, dstrxidx)
!- Make structure constants index arrays
!--------------------------------------------------------------------------------
!i Inputs:
!i   force       : whether to go to a higher harmonic in strux
!i   nbas        : number of atoms
!
!o Outputs:
!o   ssize       : size of the compact structure constants array in units of real(8)
!o   dsize       : size of the compact derivatives of the structure constants in units of real(8)
!o   struxidx    : index for struxd, shall be allocated to (nbas,tbc%escount(pid)),
!o                 tbc%escount(pid) the number of atoms allocated to process 'pid'
!o   dstrxidx    : index for dstrxd
!
!r Remark:
!r The sizes of each block are calculated and local offsets generated.
!r The indexing is transposed for cache friendliness and also ease of distribution in the parallel case
!r Neither the indices nor the structure constants are ever shared.
!r There is a 'real and present' danger the default integer(4) overflows for largeish systems,
!r  If this happens switching to integer(8) shall help but I guess there
!r  is a lot more to wory about in other places so I am leaving to int4 for now.


!    use tbprl   WHAT IS THIS FOR?
   implicit none

   integer, intent(in) :: nbas
   integer, intent(out) :: struxidx(nbas,*), dstrxidx(nbas,*), ssize, dsize
   logical :: force
   integer :: sidx, didx, pib, ib, jb, li, lj, nlmi, nlmi1, nlmj

!    call tcn('mkstrxidx')

   sidx = 1
   didx = 1

   do  ib = 1, nbas
! This needs to be altered for higher multipoles
      nlmi  = 1

      if (force ) then
         nlmi1  = 4
      else
         nlmi1  = 1
      endif

      do  jb = 1, nbas
         nlmj = 1

         struxidx(jb,ib) = sidx
         sidx = sidx + nlmj*nlmi1
         if (sidx < struxidx(jb,ib)) call panic()
!          rx('MKSTRIDX: Integer overflow. Convert to to int8 or address type')

!          if (pv) then
!             dstrxidx(jb,pib) = didx
!             didx = didx + nlmj*nlmi
!          end if
      end do
   end do

   ssize = sidx
   dsize = didx

!    call tcx('mkstrxidx')

end subroutine mkstrxidx
