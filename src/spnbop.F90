   subroutine spnbop(act)
    
!     Optionally spin polarised BOP calculation of electronic energy, charges and 
!     derivatives of the charges with respect to the onsite energies.

!         dmdn may be useful as derivative in the fermi finding routine
!   the fermi finders have to be replaced with some ignorant rootfinding
!  function only aware that the numel function only increases


      use mod_precision
      use mod_all_scalar
      use mod_const
      use ab_io
      use mod_chi
      use mod_clock
      use topologia, only : iproc, mpmap, aproc, master
      use mod_atom_ar
      
      implicit none
      
      integer, intent(in) :: act
      
      real(dp) :: ef
      character(len=25) :: ffmt
      integer :: ia
      
      integer(8) :: c1, c2
      
      procedure(real(dp)) :: getefsrt, getefnull, numel, eff
      
      include 'Include/Atom.array'
      
      if ((.not. quiet) .and. (nd < 100)) then
         write(ffmt,"(a,i0,a)") '(a,":",', nd, '(x,f22.14))'
         write(6, ffmt) 'de', de(:nd)
         if (mag) write(6, ffmt) 'dm', dem(:nd)
   !                if (mag) write(6, ffmt) 'mg', mg(:nd)
      end if
      
      
      call recurse(.true.,0) ! near = .true. instead of full cluster which is used in embedding
      
      call system_clock(c1)
      ef = lef
      if (term == 1) then
   !                   ef = getefsrt(locc)
            ef = eff(locc,numel,ef,eftol,100)
      elseif ((term == 2).or.(term == 3)) then
            ef = getefnull(locc)
      end if
      call system_clock(c2)
      
      if (.not. quiet) print *, 't(eff@root): ', real(c2-c1,dp)/real(cr,dp)
      
      do ia = mpmap(iproc)+1, mpmap(iproc+1)
         dq(ia) = dq(ia) - zc(z(ia))
         if (abs(dq(ia)) > dqmax) dqmax = abs(dq(ia))
      end do
      
      if ((.not. quiet) .and. nd < 100) then
                  write(6, ffmt) 'dq', dq(:nd)
                  if (mag) write(6, ffmt) 'mg', mg(:nd)
         write(6,'("ef:",x,f22.14)') ef    
      end if
      
      if (act > 1) then
   !             call makechi(ef)
   !             dqde(:nd) = -2*chiaia(:nd)
            call get_dq2chia(ef)
            
            
   !                 if (.not. quiet) then
   !                 write(6, ffmt) 'ch', dq2chia(:nd)
   !                 end if
            
   !                 stop 
            
            if (act>2) then
               call eval_bsens(ef)
            end if
      end if    
      
   !             if (.not. quiet) then
   !                etott = eprom + ebond + epair - eatom - uent
   !                write(6,'(/,"it:",i3)') it
   !                write(6,'("ef:",x,f22.14)') lef
   !                write(6, ffmt) 'de', de(:nd)
   !                write(6, ffmt) 'dq', dq(:nd)
   !                write(6, ffmt) 'ch', dq2chia(:nd)
   !    !             write(6,"('de,dq:',2(x,f22.14))") de(9),dq(9)
   !                
   !                if (abs(dqmax) <= qerr) then 
   !                   write(9,"('c')",advance='no')
   !                else if (it == 1) then
   !                   write(9,"('s')",advance='no')
   !                else if (it == mxit) then
   !                   write(9,"('x')",advance='no')
   !                else
   !                   write(9,"('i')",advance='no')
   !                end if
   !                
   ! !                write(9,"(2(x,f22.14),x,i3)") etott, dqmax, it
   !                write(9,"(1(x,f22.14),x,i3)") dqmax, it
   !                
   !             end if

      lef = ef
      
      
      eelct = ((ebond + eprom) - uent) - eatom + emag
     
   end subroutine spnbop
