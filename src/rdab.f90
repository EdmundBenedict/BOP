 
      subroutine rdab(flag,ia)
          use mod_precision


!
!    This is a routine to read in the recursion coefficients from
!     high-speed disk in the temp directory.
!
      use ab_io
      implicit none

!
!    Include constants.
!

!      include "Include/ALL.const"

!
!    Include scalars.
!

!      include "Include/ALL.scalar"

!
!    Include arrays.
!

!      include "Include/Atom.array"
!      include "Include/Moment.array"
!      include "Include/NebList.array"
!      include "Include/NotSRT.array"
!      include "Include/RecurCoef.array"
!      include "Include/SRT.array"

!
!    Declare the simple variables.
!
!     logical, optional, intent(in) :: hstore
      integer flag,i,j
      integer ia,la,ib,lb,ja,ja0,jb,n
      integer ma,nma,lla0,nla,nstta

      character*80 filename 
      character*1024 scratch
! 
!       integer(8) :: cr, cm, c1, c2
!       integer(8), save :: cc=0
! 
! !       write(6,*) ia
! 
!       call system_clock(count_rate = cr, count_max = cm )
! 
! 
!       
!        call system_clock(c1)
!        
!     if (present(hstore)) then 
!         if (.not. hstore) then
!             call assoc_ab(ia)
!         else
!       
! 
!         !
!         !    Test flag to see if need to open file.
!         !
! 
!             if (flag.eq.1) then
!                 filename = genfile(1:lengfn)//'.ab'
!                 call get_environment_variable('BOPS_SCRATCH', scratch)
!                 if (len(trim(scratch))<1) scratch = 'tmp'
!                 open(unit=50,file=trim(scratch)//'/'//filename,form='UNFORMATTED', & 
!             &        status='OLD')
!                 flag = 0
!             endif
! 
!         !
!         !    Read in data.
!         !
!             call states(z(ia),nla,nstta,llista)
!             lla0 = 0
!             read(50) nchain
!             do la = 1,nchain,1
!                 if (momflg.eq.1) nma = 2*llista(la) + 1
!                 read(50) lchain(la),wt(la)
!                 read(50) (arec(n,la),brec(n,la),n=0,lchain(la),1)
!                 read(50) lainf(la),lbinf(la)
!                 ja = aptr(ia)
!                 ja0 = ja
!                 do while(bptr(ja).ne.eol)
!                     jb = ja-ja0+1
!                     ib = bptr(ja)
!                     if (momflg.eq.1) then
!                     do ma = 1,nma,1
!                         do lb = 1,nstt(z(ib)),1
!                             read(50) (darec(n,lla0+ma,lb,jb), dbrec(n,lla0+ma,lb,jb), n=0,lchain(la)+1,1)
!                         enddo
!                     enddo
!                     elseif ((momflg.eq.2).or.(momflg.eq.3)) then
!                     do lb = 1,nstt(z(ib)),1
!                         read(50) (darec(n,la,lb,jb),dbrec(n,la,lb,jb), n=0,lchain(la)+1,1)
!                     enddo
!                     endif
!                     ja = ja + 1
!                 enddo
!                 if (term.eq.1) then
!                     read(50) forder(la)
!                     read(50) (root(i,la),f1(i,la),f2(i,la), i=1,forder(la),1)
!                 elseif ((term.eq.2).or.(term.eq.3)) then
!                     read(50) (diag(i,la),i=1,lchain(la)+1,1)
!                     read(50) ((eigvec(i,j,la), i=1,lchain(la)+1,1), j=1,lchain(la)+1,1)
!                 endif
!                 if (momflg.eq.1) lla0 = lla0 + nma
!             enddo
!             
!           end if
!       else
        if (flag==1) flag = 0
        call assoc_ab(ia,1)
!       end if

!         call system_clock(c2)
!        
!        cc = cc + c2 - c1
!        write(6,'("rdab wallclock time: ",g,"s cc:",g,"s")'), real(c2 - c1,8)/real(cr,8), cc/real(cr,8)
!        
      end

