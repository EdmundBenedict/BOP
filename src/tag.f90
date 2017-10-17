 
!****************************************************************************
!                                                                           *
!                                                                           *
!                             Alexey Girshick                               *
!              Department of Materials Science and Engineering              *
!                        University of Pennsylvania                         *
!                                                                           *
!                                                                           *
!****************************************************************************
!                                                                           *
!                                                                           *
!             This subroutine puts a tag on the output file.                *
!                                                                           *
!                                                                           *
!****************************************************************************


      subroutine tag( num )
          use mod_precision


          use mod_all_scalar

          use mod_const

          use mod_srt
          use mod_atom_ar

!
!    This is a routine to write data to the output file.
!

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

      include "Include/Atom.array"
      include "Include/BEC.array"
!      include "Include/BondOrder.array"
!      include "Include/Force.array"
!      include "Include/Hamilt.array"
!      include "Include/Misc.array"
!      include "Include/Moment.array"
!      include "Include/NebList.array"
!      include "Include/NotSRT.array"
      include "Include/PosVel.array"
!      include "Include/RecurCoef.array"
!      include "Include/Relax.array"
!      include "Include/SRT.array"

!
!    Common blocks introduced by A.Girshick.
!    
 
      include "Include/ag.conn"
      include 'Include/Envir.array'

!
!    Declare the simple variables.
!

      integer           i, j, k, ia, ib, num, za


      write ( num, '( /80("_") )' )
      write ( num, '( /6X,"The results above have been obtained using" &
     &                    " the following parameters:")' )
      write ( num, '(  80("_")// )' )  

      if ( scf_include  ==  1 ) then
        write (num,'(/22X,"SCREENING OF BOND INTEGRALS INCLUDED"//)')
      else
        write ( num, '( /25X,"NO SCREENING OF BOND INTEGRALS"//)' )
      endif

      if      ( fs_include  ==  0 ) then
          write ( num, '( /27X,"ONLY BONDING PART INCLUDED"//)' )
      else if ( fs_only  ==  0  .and.  emb_include  ==  0 ) then   
          write ( num, '( /20X,"BONDING PART AND PAIR POTENTIAL" &
     &                         " INCLUDED"//)' )
      else if ( fs_only  ==  0  .and.  emb_include  ==  1 ) then
          write ( num, '( /20X,"BONDING PART AND FINNIS-SINCLAIR" &
     &                         " INCLUDED"//)' )
      else if ( fs_only  ==  1  .and.  emb_include  ==  0 ) then
          write ( num, '( /26X,"ONLY PAIR POTENTIAL INCLUDED"//)' )
      else if ( fs_only  ==  1  .and.  emb_include  ==  1 ) then
          write ( num, '( /26X,"ONLY FINNIS-SINCLAIR INCLUDED"//)' )
      endif

      if ( momflg  ==  1 )  then
          write( num,'(/18X, "MOMFLG:   ",I2,"     (using averaged" &
     &                       " moments)" )' )  momflg
      else 
          write( num,'(/18X, "MOMFLG:   ",I2 )' )  momflg
      endif
      if ( term  ==  1 )  then
          write( num,'(/18X, "TERM:     ",I2,"     (using square root" &
     &                       " terminator)"  )' )  term
      else
          write( num,'(/18X, "TERM:     ",I2 )' )  term
      endif

      write( num,'(/18X, "NREC:     ", I2, 5X, &
     &           "(total number of recursion levels)" )' )  nrec
      write( num,'(/18X, "NBASE:    ", I2, 5X, &
     &           "(number of exact recursion levels)" )' )  nbase

      if      ( chi_meth  ==  1 )  then
          write( num,'(/18X, "CHI_METH: ",I2,"     (using analytical" &
     &                       " integration)" )' )  chi_meth
      else if ( chi_meth  ==  2 )  then
          write( num,'(/18X, "CHI_METH: ",I2,"     (using numerical" &
     &                       " integration)" )' )  chi_meth
          write( num,'(/18X, "MFAC:    ", I3, 5X, &
     &               "(number of points prefactor)" )' )  mfac
      endif
      write( num,'(//18X, "KT:     ",F4.2 , 5X, &
     &           "(electronic temperature)" )' )  kt
      if ( qerr  >  0.d0  )  then
          write( num,'(/18X, "QERR:   ",F4.2"     (using local " &
     &                       " charge neutrality)" )' )  qerr
      else
          write( num,'(/18X, "QERR:    .00     (without local"  &
     &                       " charge neutrality)" )' )
      endif
      write( num,'(/18X, "RCUT:   ",F4.2 , 5X, &
     &           "(bond part cut-off)" )' )  rcut
      write( num,'(/18X, "RCUTL:  ",F4.2 , 5X, &
     &           "(pair potential cut-off)" )' )  rcutl
      write( num,'(/18X, "RPRUNE: ",F4.1 , 5X, &
     &           "(width of inert atoms zone)" )' )  rprune

      if ( mvflag  ==  12 )  then
          write( num, '( //18X, "TRANGE:   ",2(F7.5,2X),I2 )' )  &
     &                   tlo, thi, nt
          write( num,'(/18X, "POLYORD:  ", I1 )' ) polyord
          write( num,'(/18X, "STRTYPE:  ", A3 )' ) strtyp
      endif

      if ( mvflag  ==  14 )  then
          if ( indrelax  ==  0 )  then
              write( num,'(//18X,"INDREL:    0",5X,"(gamma-surface" &
     &                       " not relaxed)"  )' ) 
          else if ( indrelax  ==  1 )  then
              write( num,'(//18X,"INDREL:    1",5X,"(gamma-surface"     &
     &                       " relaxed)"  )' )     
          endif
      endif

      if ( mvflag  ==  15 )  then
          write( num,'(//18X, "XSHIFTGB: ",F6.4, &
     &           " (shift in the x-direction)" )' )  xshiftgb
          write( num,'(/ 18X, "YSHIFTGB: ",F6.4,      &
     &           " (shift in the y-direction)" )' )  yshiftgb
          write( num,'(/ 18X, "ZSHIFTGB: ",F6.4, &
     &           " (shift in the z-direction)" )' )  zshiftgb
      endif

      if ( ( mvflag  >=  14 ) .and.  ( mvflag  <=  17 ) )  then
          write( num,'(//18X, "OK:   ",F6.4, 5X, &
     &           "(relaxation constant)" )' )  ok
          write( num,'(/ 18X, "FTOL: ",F6.4, 5X, &
     &           "(required force tolerance)" )' )  ftol           
          write( num,'(/ 18X, "MXITER: ",I7, 5X, &
     &           "(maximum number of iterations)" )' )  mxiter
          if ( const_volume  ==  0  .and.  mvflag  <=  15 )  then
              write( num,'(/ 18X,"RELTYPE:   1",5X,"(relaxation" &
     &                       " at constant pressure)"  )' )
          else if ( const_volume  ==  1 .and. mvflag  <=  15 ) then
              write( num,'(/18X, "RELTYPE:   2",5X,"(relaxation" &
     &                       " at constant volume)"  )' )
          endif
      endif


      if ( fs_only  ==  0 ) then
         write( num,'(///3X, "Atomic data:",11X,"NEL",7X,"ES", &
     &                   7X,"EP",7X,"ED")' )
         do i = 0, natype
           za = atype(i)
           write( num, '(/ 19X, A3, 2X, F6.3, 3(2X, F7.4) )' )  &
     &                listsymb(za), zc(za), es(za), ep(za), ed(za)
         enddo

         write( num,'(//3X, "Bond scaling:",8X,"R0",6X,"RC", &
     &                 6X,"R1",5X,"RCUT",6X,"N",6X,"NC")' )
         do i = 0, mxz
           do j = 0, mxz
             if ( btype(i,j)  /=  0 ) then
                ib = btype(i,j)
                do k=1,15
                if (bndscl(1,k,ib) /= 0.0) then
                write( num,'(/12X,A2," - ",A2,1X,5(2X,F6.4),F7.4)' )  &
     &             listsymb(i),listsymb(j),(bndscl(ia,k,ib),ia=1,6)
                endif
                enddo
                write(9,'(A)') '------------------------'
             endif
           enddo
         enddo

         write( num,'(//"Hopping integrals:")' )
         do i = 0, mxz
           do j = 0, mxz
             if ( btype(i,j)  /=  0 ) then
               ib = btype(i,j)          
               do ia = 1, 14   
                 if ( bnddat(ia,ib)  /=  0.0_dp ) then
                   if ( ia  ==   1 ) then
                        write( num, '(25X, A2," - ",A2,2X,"V(ss_s):", &
     &                  F9.4/)') listsymb(i), listsymb(j), bnddat(ia,ib)
                   endif
                   if ( ia  ==   2 ) then
                        write( num, '(25X, A2," - ",A2,2X,"V(sp_s):", &
     &                  F9.4/)') listsymb(i), listsymb(j), bnddat(ia,ib)
                   endif
                   if ( ia  ==   3 ) then
                        write( num, '(25X, A2," - ",A2,2X,"V(ps_s):", &
     &                  F9.4/)') listsymb(i), listsymb(j), bnddat(ia,ib)
                   endif
                   if ( ia  ==   4 ) then
                        write( num, '(25X, A2," - ",A2,2X,"V(pp_s):", &
     &                  F9.4/)') listsymb(i), listsymb(j), bnddat(ia,ib)
                   endif
                   if ( ia  ==   5 ) then
                        write( num, '(25X, A2," - ",A2,2X,"V(pp_p):", &
     &                  F9.4/)') listsymb(i), listsymb(j), bnddat(ia,ib)
                   endif
                   if ( ia  ==   6 ) then
                        write( num, '(25X, A2," - ",A2,2X,"V(sd_s):", &
     &                  F9.4/)') listsymb(i), listsymb(j), bnddat(ia,ib)
                   endif
                   if ( ia  ==   7 ) then
                        write( num, '(25X, A2," - ",A2,2X,"V(ds_s):", &
     &                  F9.4/)') listsymb(i), listsymb(j), bnddat(ia,ib)
                   endif
                   if ( ia  ==   8 ) then
                        write( num, '(25X, A2," - ",A2,2X,"V(pd_s):", &
     &                  F9.4/)') listsymb(i), listsymb(j), bnddat(ia,ib)
                   endif
                   if ( ia  ==   9 ) then
                        write( num, '(25X, A2," - ",A2,2X,"V(dp_s):", &
     &                  F9.4/)') listsymb(i), listsymb(j), bnddat(ia,ib)
                   endif
                   if ( ia  ==  10 ) then
                        write( num, '(25X, A2," - ",A2,2X,"V(pd_p):", &
     &                  F9.4/)') listsymb(i), listsymb(j), bnddat(ia,ib)
                   endif
                   if ( ia  ==  11 ) then
                        write( num, '(25X, A2," - ",A2,2X,"V(dp_p):", &
     &                  F9.4/)') listsymb(i), listsymb(j), bnddat(ia,ib)
                   endif
                   if ( ia  ==  12 ) then
                        write( num, '(25X, A2," - ",A2,2X,"V(dd_s):",  &
     &                  F9.4/)') listsymb(i), listsymb(j), bnddat(ia,ib)
                   endif
                   if ( ia  ==  13 ) then
                        write( num, '(25X, A2," - ",A2,2X,"V(dd_p):",  &
     &                  F9.4/)') listsymb(i), listsymb(j), bnddat(ia,ib)
                   endif
                   if ( ia  ==  14 ) then
                     write( num, '(25X, A2," - ",A2,2X,"V(dd_d):", &
     &               F9.4/)') listsymb(i), listsymb(j), bnddat(ia,ib)
                   endif
                 endif
               enddo
             endif
           enddo
         enddo

      endif

!  The following piece was written for prinitng out the Finnis-Sinclair
!  potential for pure element. It has to be modified for the case of binary 
!  alloy  potential.

      if ( fs_include  ==  1 ) then
         write( num,'(//26X, "Pair potential  " &
     &            A2," - ",A2// )' ) listsymb( 22 ), listsymb( 22 ) 
          do i = 1, 6
              write( num,'(16X,"A\(",I1," ) = ",F18.14,6X,"R\(", &
     &             I1," ) = ",F6.4/)' ) i,arep( i, 1 ),i,arep( i, 2 )
          enddo
         write( num,'(//26X, "Pair potential  " &
     &            A2," - ",A2// )' ) listsymb( 13 ), listsymb( 13 ) 
          do i = 1, 6
              write( num,'(16X,"A\(",I1," ) = ",F18.14,6X,"R\(", &
     &             I1," ) = ",F6.4/)' ) i,brep( i, 1 ),i,brep( i, 2 )
          enddo
         write( num,'(//26X, "Pair potential  " &
     &            A2," - ",A2// )' ) listsymb( 13 ), listsymb( 22 ) 
          do i = 1, 6
              write( num,'(16X,"A\(",I1," ) = ",F18.14,6X,"R\(", &
     &             I1," ) = ",F6.4/)' ) i,abrep( i, 1 ),i,abrep( i, 2 )
          enddo


          if ( emb_include  ==  1 ) then

              do i = 1, mxz       
                  if ( btype ( i, i )  /=  0 )  then
                     write( num,'(//24X, "Embedding potential  " &
     &                A2," - ",A2// )' ) listsymb( i ), listsymb( i )
                  endif
              enddo
              do i = 1, 2
                  write( num,'(16X,"a\(",I1," ) = ",F18.14,6X,"r\(", &
     &                 I1," ) = ",F6.4/ )' ) i,acoh(i,1),i,acoh(i,2)    
              enddo
          endif

          write( num,'(/31X, "Repulsive core ")' )
          write( num,'(/7X,"R_0 = ",F4.2,3X,"A = ",F6.2,3X, &
     &          "Alpha = ",F6.2,3X,"NR = ",F4.2,3X,"N = ",F3.1//)' ) &
     &           (acore( i ), i=1,5 )


          if ( env_include  ==  1 ) then


!              DO I = 1, MXZ
!                  IF ( BTYPE ( I, I ) .NE. 0 )  THEN
!                      WRITE( NUM,'(/24X, "Environmental term  "
!     &                A2," - ",A2/ )' ) LISTSYMB( I ), LISTSYMB( I )



!                 ENDIF

!             ENDDO  

!           WRITE( NUM,'(7X,"Aenv = ",F11.5,13X,"lam0 = ",F9.6//,7X,
!     &          "d0   = ",F5.3,19X,"Rcut = ",F5.2)' )
!     &           (AENV( I ), I=1,4 )


!           WRITE( NUM,'(/7X,"Cenv = ",F12.4,17X,"ny   = ",F9.6//,7X,
!     &          "n    = ",F4.2,20X,"Rtail = ",f5.2//)' )
!     &           (ADEL( I ), I=1,4 )             
! 

!             ENDIF  
!         
!             ENDDO  
                     write( num,'(/24X, "Environmental term  " &
     &                A2," - ",A2/ )' ) listsymb( 22 ), listsymb( 22 )

           write( num,'(7X,"Aenv = ",F11.5,13X,"lam0 = ",F9.6//,7X, &
     &          "d0   = ",F5.3,19X,"Rcut = ",F5.2)' ) &
     &           (aenv( i ), i=1,4 )


           write( num,'(/7X,"Cenv = ",F12.4,17X,"ny   = ",F9.6//,7X, &
     &          "n    = ",F4.2,20X,"Rtail = ",f5.2//)' ) &
     &           (adel( i ), i=1,4 )             

                     write( num,'(/24X, "Environmental term  " &
     &                A2," - ",A2/ )' ) listsymb( 13 ), listsymb( 13 )

           write( num,'(7X,"Aenv = ",F11.5,13X,"lam0 = ",F9.6//,7X, &
     &          "d0   = ",F5.3,19X,"Rcut = ",F5.2)' ) &
     &           (benv( i ), i=1,4 )


           write( num,'(/7X,"Cenv = ",F12.4,17X,"ny   = ",F9.6//,7X, &
     &          "n    = ",F4.2,20X,"Rtail = ",f5.2//)' ) &
     &           (bdel( i ), i=1,4 )             

                     write( num,'(/24X, "Environmental term  " &
     &                A2," - ",A2/ )' ) listsymb( 13 ), listsymb( 22 )

           write( num,'(7X,"Aenv = ",F11.5,13X,"lam0 = ",F9.6//,7X, &
     &          "d0   = ",F5.3,19X,"Rcut = ",F5.2)' ) &
     &           (abenv( i ), i=1,4 )


           write( num,'(/ 7x, "Cenv = ", f12.4, 17x, "ny   = ", f9.6, 2/, 7x, &
     &          "n    = ", f4.2, 20x, "Rtail = ", f5.2, /, /)' ) &
     &           (abdel( i ), i=1,4 )             

          endif

      endif


      return
      end



