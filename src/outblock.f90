

   subroutine outblock( key )
      
      use mod_precision
      use mod_all_scalar
      use mod_const
      use mod_atom_ar, only: mg,  demi
      implicit none

      integer key

      include "Include/Atom.array"
      include "Include/Misc.array"
      include "Include/PosVel.array"

      include "Include/ag.conn"


      character*80 filename
      integer  i, num, j


      if ( key  ==  0 )  then
         open(unit = 87,file = "block.interm")
         num = 87          
      else if ( key  ==  1 )  then
         num = 32
      else if ( key  ==  2 )  then
         write(filename,'("rad.block.",i0)') klp
         open(unit = 87,file = filename)
         num = 87
      endif

      if ( mvflag  ==  18 )  then
         write(num,'(''LOG_RAD'')')
!
!     Modified by M. Cawkwell: ln(R/b), ie, r_0 = b = LENA(3)
!
         write(num,'(F18.12)') log( radius / lena(3) )   
      endif

!
!      Write out the primitive translation vectors.
!
         
      write(num,'(''IF_RELAXED'')')
      write(num,'("1")')
      write(num,'(''ITER'')')
      write(num,'(I4/)') iter
      
      if (mvflag == 18) then
         write(num,'(''E_COH'')')
         write(num,'(F16.12)') e_coh
      endif
      
      write(num,'(''A'')')
      write(num,'(3(3(x,F8.5),/))') 1.d0, 0.d0, 0.d0, &
                              & 0.d0, 1.d0, 0.d0, &
                              & 0.d0, 0.d0, 1.d0
      
      write(num,'(''LEN'')')
      write(num,'(3(F16.12,2X))') a(1,1),a(2,2),a(3,3)     
      write(num,'(''LATPAR'')')
      write(num,'(F18.12)')       latpar
!      WRITE(NUM,'(''LIMITS'')')
!      WRITE(NUM,'(2(F16.12,2X))') XXMIN, XXMAX
!      WRITE(NUM,'(2(F16.12,2X))') YYMIN, YYMAX        
!      WRITE(NUM,'(2(F16.12,2X))') ZZMIN, ZZMAX        

!      WRITE(NUM,'(''NULLHEIGHT'')')
!      WRITE(NUM,'(F18.12)')       ZEROHEIGHT

!
!    Write out the number of active atoms.
!
   
      write(num,'(''ND'')')
      write(num,'(I6)') nd

!
!    Write out the final atomic positions and velocities of active atoms.
!
                                               
      write(num,'(''D'')')
      do i = 1, nd
         write(num,'(3(2x,f16.12),3x,a2,3x,f16.12)') (ad(j,i)/a(j,j), j = 1, 3), symb(i), mg(i)
      end do

!    
!    Write out the number of inert atoms.
!

      write(num,'(''NINERT'')')
      write(num,'(I6)') ninert

!     
!    Write out the positions of the inert atoms.
!

      write(num,'(''DINERT'')')
      do i = 1, ninert
         write(num,'(3(2X,F16.12),3X,A2,3x,f16.12)') &
            & (adinert(j,i)/a(j,j), j = 1, 3), symi(i), (2*demi(i))/istn(zinert(i))
      end do

!    
!    Write out the unrelaxed positions of active atoms.
!           
         
      write(num,'(''UNRLD'')')
      do i = 1, nd
         write(num,'(3(2X,F16.12),3X,A2)') (unrld(j,i)/a(j,j), j =1, 3), symb(i)
      end do

!
!    Write out the energy shifts and EF level.
!    

      write(num, '(''EFLO'')')
      write(num, '(F16.12)') eflo
      write(num, '(''DE'')')
      
      write(num, *) de(:nd)

!      CLOSE( NUM )

   end subroutine outblock

