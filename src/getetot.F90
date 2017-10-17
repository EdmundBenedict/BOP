
   subroutine getetot(flag)
      use mod_precision
      use mod_all_scalar
      use mod_const
      use mod_conf, only : rtc
      use tbbop_emb, only : emb_sc
      
      use mod_forcedetails, only: print_forcedetails
      
      
!       REMOVE THIS
      use topologia, only: iproc,master

   !
   !    This is a routine to evaluate the total energy and forces for a configuration.
   !
   !     FLAG = 0 => Read the coefficients from disk, and find Fermi energy.
   !     FLAG = 1 => Evaluate the coefficients, write them to disk, and find Fermi energy.
   !     FLAG = 2 => Evaluate the coefficients, do not write them to disk, and do not find Fermi energy.
   !

      implicit none

      include "Include/Atom.array"
      include "Include/Force.array"
      include "Include/PosVel.array"
      include "Include/NebList.array"
      include "Include/NebListL.array"
      include "Include/ag.conn"


      integer, intent(in) :: flag
      real(dp) ::   repeng, u1, rxab,ryab,rzab
      real(dp) ::   scf(14),scfcut(14)

      integer ::  i, j, ii
      integer ::  ia, ib, ja, bn, ianb
      real(dp) :: ndr, ett(3,3)
      logical :: stressed

      if (forces) then
         ftot(:,:nd) = 0.0_dp
         fatom(:,:nd) = 0.0_dp
         fbs(:,:nd) = 0.0_dp
         fpp(:,:nd) = 0.0_dp
         fblock(:,:nd+1) = 0.0_dp
      end if
         
      ebond = 0.0_dp     
      eprom = 0.0_dp     
      eband = 0.0_dp     
      emag  = 0.0_dp      
      uent  = 0.0_dp     
      eelct = 0.0_dp     
      eenv  = 0.0_dp     
      epair = 0.0_dp     
      eclas = 0.0_dp    
      ebind = 0.0_dp                              
      
      ett = 0.0_dp
      
      stressed = totstr /= 0.0_dp .and. any(et /= 0.0_dp)      
      if (stressed) ett = et*totstr
      
      if ( fs_only  ==  0 )  then

         call bldnebt(.false., rlim, lena, a, rprune, rcut, stressed, ett, &
                   &  nd, ad, z, ninert, adinert, zinert, &
                   &  map, aptr, bptr, maxnumnb, totnd, nne, mapi,mxnnb)
         
         if (potflg == sbop) then
!             if (.not. quiet) then
!                write(9,*) "****************************************"
!                write(9,*) "*   Real-space calculation of energy   *"
!                write(9,*) "****************************************"
!             end if

            call getebsfbs(flag)

#ifdef DISLDBG
            call print_forcedetails()
#endif
               
         elseif (potflg == kspace) then
!             write(*,*) "****************************************"
!             write(*,*) "*     K-space calculation of energy    *"
!             write(*,*) "****************************************"

            call kebsfbs()
            
         endif

         if (potflg == tbbop .or. potflg == 10) then
            call emb_sc()
         
         end if
         eelct = eelct + ((ebond + eprom) - uent) - eatom + emag
         

      endif
      
      
      if ((fs_include  ==  1 .or. env_include  ==  1 .or. vpair_include  ==  1)   &
               &   .and. (potflg /= 9 .and. potflg /=10)) then     
         if (stressed) then
            ad(:,:nd) = unshad(:,:nd)
            adinert(:,:ninert) = unshind(:,:ninert)
         end if

! embedding potentials need proper map() while the rest of the code uses map=0 for inert atoms.
         call bldnebt(.true., rlim, lena, a, rprune, rcutl, stressed, ett, &
                   &  nd, ad, z, ninert, adinert, zinert, &
                   &  map, aptrl, bptrl, maxnumnb, totnd, nne, mapi,mxnnb)  ! mapi will not be accessed  
         
         call classic()
         
         !    Add in one body energy and forces.
         
         if (iproc == master) then
            print *, 'FORCES'
            do i = 1, nd
                print *, 'ATOM:',i
                print *, 'FBS:',fbs(:,i)
                print *, 'FPP',fpp(:,i)
                print *, 'FTOT:',fbs(:,i)+fpp(:,i)
                print *,' '
            enddo
         endif
         if (v1flag > 0) then
            call onebdy(u1,fpp)
            epair = epair + u1
         endif      
      endif
      

      
      
      ebind = eelct + eclas
      
      if (iproc == master) then
          print *, 'eelct',eelct
          print *,'ebond',ebond
          print *,'eprom',eprom
          print *,'uent',uent
          print *,'eatom',eatom
          print *,'emag',emag
          print *,'eclas',eclas
          print *,'ebind',ebind
      endif
      
      rtc % tote % bond = ebond
      rtc % tote % prom = eprom
      rtc % tote % band = eband
      rtc % tote % atom = eatom
      rtc % tote % mag  = emag
      rtc % tote % ent  = uent 
      rtc % tote % elct = eelct
      rtc % tote % env  = eenv 
      rtc % tote % pair = epair
      rtc % tote % clas = eclas

      rtc % tote % bind = ebind

      ndr = real(nd,dp)  
         
      rtc % avre % bond = ebond/ndr
      rtc % avre % prom = eprom/ndr
      rtc % avre % band = eband/ndr
      rtc % avre % atom = eatom/ndr
      rtc % avre % mag  = emag /ndr
      rtc % avre % ent  = uent /ndr 
      rtc % avre % elct = eelct/ndr
      rtc % avre % env  = eenv /ndr
      rtc % avre % pair = epair/ndr
      rtc % avre % clas = eclas/ndr
      
      rtc % avre % bind = ebind/ndr


!       This is for all the practical routines which do their own energy sum and thus rely on ebond.
      if (mag) ebond = ebond + emag
      
      

      if (forces) then
         ftot(:,:nd) = fbs(:,:nd) + fpp(:,:nd)
!                   
!          fblock(:,nd+1) = 0.0_dp
!          do i = nd, 1, -1
!             fblock(:,i) = fblock(:,i+1) + ftot(:,i)
!          end do
         
!          if (dndm > 1.0e-6_dp) then
!             dmdl(:,:nd) = dmdl(:,:nd)/dndm
!          else
!             dndm = 0.0_dp
!             dmdl(3,:nd) = 0.0_dp
!          endif
!          
         dndm = 0.0_dp
         dmdl = 0.0_dp
!         if (.not. quiet) call print_forces('', 6)
      end if
      
   end subroutine getetot





