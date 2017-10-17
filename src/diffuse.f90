 
 
      subroutine diffuse()
          use mod_precision


          use mod_all_scalar

          use mod_const

!     This routine calculates a diffusion barrier.
!
!     The method employed is to push a chosen atom (D_ATOM)
!     in a certain direction (DISP_VEC) in a series of steps (N_DISP).
!
!     At each point on the diffusion path the diffusing atom is
!     constrained to lie in the plane perpendicular to the direction
!     of motion. All other atoms are relaxed unconstrained using a
!     conjugate gradients algorithm.
!
!     In order to prevent the whole unit cell just translating
!     with the diffusing atom a sphere of a given radius (RC_DIFF) is
!     constructed around the diffusing atom. All atoms which lie outside
!     this sphere are transformed into `inert' atoms and do not contribute
!     to the total energy.
!
!     The routine creates a file with extension '.diff' which stores
!     details of the distance down the diffusion path (diffusion co-ordinate),
!     energy difference (relative to the relaxed initial location of the
!     diffusing atom) along with the diffusing atom's position and force upon it.
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
      include "Include/Force.array"
!      include "Include/Hamilt.array"
      include "Include/Misc.array"
!      include "Include/Moment.array"
!      include "Include/NebList.array"
!      include "Include/NotSRT.array"
      include "Include/PosVel.array"
!      include "Include/RecurCoef.array"
      include "Include/Relax.array"
!      include "Include/SRT.array"

!
!    Declare the simple variables.
!

      real(dp) :: dist,fac,rxqd,ryqd,rzqd
      real(dp) :: d_dist,eref,geteatom

      integer q,qq
      integer nactive

      character*80 dum
      character*80 filename

!
!    Declare local arrays.
!

!*** The initial and final co-ordinates of the diffusing atom.
      real(dp) :: d_start(3), d_end(3)
!*** The unit vector of diffusion.
      real(dp) :: d_unit_vec(3)

!
!    Trap constrained motion.
!

      if (cnst_n > 0) then
         print *,'[DIFFUSE] THIS ROUTINE CAN NOT HANDLE &  
     & EXPLICITLY CONSTRAINED ATOMS AT THE MOMENT'
         call panic()
      endif
         
!   Start and end points of diffusion path.
      do q = 1,3,1
         d_start(q) = ad(q,d_atom)
         d_end(q)   = d_start(q) + disp_vec(q)
      enddo

!   Length of diffusion path.
      d_dist = sqrt( (d_end(1)-d_start(1))**2 + &  
     &     (d_end(2)-d_start(2))**2 +  &  
     &     (d_end(3)-d_start(3))**2 )
      
!   Unit vector along diffusion path.
      do q=1,3
         d_unit_vec(q) = (d_end(q)-d_start(q))/d_dist
      enddo

      ninert  = 0
      nactive = 0

!     Copy diffusing atom into temporary storage
!     so the diffusing atom is always at position 1
!     in the co-ordinate data structure.
      nactive = nactive + 1
      newad(1,nactive) = ad(1,d_atom)
      newad(2,nactive) = ad(2,d_atom)
      newad(3,nactive) = ad(3,d_atom)
      newz(nactive)    = z(d_atom)

      do q = 1, nd
         if (q /= d_atom) then
            rxqd = ad(1,q) - ad(1,d_atom)
            ryqd = ad(2,q) - ad(2,d_atom)
            rzqd = ad(3,q) - ad(3,d_atom)
            if (rxqd > lena(1)/2.0_dp) then
               rxqd = rxqd - lena(1)
            elseif (rxqd < -lena(1)/2.0_dp) then
               rxqd = rxqd + lena(1)
            endif
            if (ryqd > lena(2)/2.0_dp) then
               ryqd = ryqd - lena(2)
            elseif (ryqd < -lena(2)/2.0_dp) then
               ryqd = ryqd + lena(2)
            endif
            if (rzqd > lena(3)/2.0_dp) then
               rzqd = rzqd - lena(3)
            elseif (rzqd < -lena(3)/2.0_dp) then
               rzqd = rzqd + lena(3)
            endif
            
            dist = sqrt(rxqd*rxqd + ryqd*ryqd + rzqd*rzqd)

            if (dist > rc_diff) then
! Atom will be fixed. Copy into inert atom array.

               write(6,'(A,I4,4F12.5)')'FIXING ATOM ',q, &  
     &              ad(1,q),ad(2,q),ad(3,q),dist
               ninert            = ninert + 1
               if (ninert > minert) then
                  print *,'Too many inert atoms. Increase MINERT'
                  call panic()
               endif
               adinert(1,ninert) = ad(1,q)
               adinert(2,ninert) = ad(2,q)
               adinert(3,ninert) = ad(3,q)
               zinert(ninert)    = z(q)

            else
! Atom will be active. Copy into temporary storage.

               nactive = nactive + 1
               newad(1,nactive) = ad(1,q)
               newad(2,nactive) = ad(2,q)
               newad(3,nactive) = ad(3,q)
               newz(nactive)    = z(q)

            endif
         endif
      enddo

      write(6,'(A,I4)') 'NUMBER OF ATOMS FIXED FOR DIFFUSION =',ninert

!     Form new shorter AD array for active atoms.
      do q = 1, nactive
         ad(1,q) = newad(1,q)
         ad(2,q) = newad(2,q)
         ad(3,q) = newad(3,q)
         z(q)    = newz(q)
      enddo
      nd = nactive

      print *,'SO  NUMBER OF ACTIVE ATOMS = ',nd
      print *,'AND NUMBER OF INERT  ATOMS = ',ninert

!     Reset chain characteristics and occupancies.
      call getocc(1,locc)
      eatom = geteatom(z,nd)

      if (idebug == 1) then
         print *,'***** ACTIVE ATOMS *****'
         do q = 1, nd
            write(6,'(i4,3f12.5,i4)') q,ad(1,q),ad(2,q),ad(3,q),z(q)
         enddo
         print *,'***** INERT ATOMS ******'
         do q = 1, ninert
            write(6,'(i4,3f12.5,i4)') q,adinert(1,q),adinert(2,q), &
     &           adinert(3,q),zinert(q)
         enddo
      endif

      write(6,*)
      write(6,*) '------------------- DIFFUSION RUN PARAMETERS ----'
      write(6,'(A,3F12.5)') 'DIFFUSING ATOM STARTS AT : ', & 
     &     d_start(1),d_start(2),d_start(3)
      write(6,'(A,3F12.5)') 'AND FINISHS  AT : ',d_end(1),d_end(2), & 
     &     d_end(3)
      write(6,'(A,F12.5)')  'SO THE DIFFUSION DISTANCE IS ',d_dist
      write(6,'(A,3F12.5)') 'DIFFUSION VECTOR AND PLANE OF CONSTRAINT=', & 
     &     (d_unit_vec(q),q=1,3)
      write(6,'(A,I5)') 'NUMBER OF DIFFUSION STEPS BEING PERFORMED', & 
     &     n_disp
      write(6,'(A)') '--------------------------------------------'
      write(6,*)
      write(6,*)


!     Constrain diffusing atom to lie in a plane perpendicular
!     to its diffusion direction.
      d_atom = 1
      cnst_n = 0
      cnst_n = cnst_n + 1
      cnst_a(cnst_n) = d_atom
      do q = 1,3
         cnst_v(q,cnst_n) = d_unit_vec(q)
      enddo
      
      print *,'NUMBER OF CONSTRAINTS = ',cnst_n
      do q = 1, cnst_n
         write(6,'(2(A,I4),A,3F12.5)') 'INDEX = ',q, &
     &        ' ATOM = ',cnst_a(q), &
     &        ' CONSTRAINT = ',(cnst_v(qq,q),qq=1,3)
      enddo
      
      if (n_disp /= 0) then
         fac = d_dist/real(n_disp, dp)
      else
         fac = 0.0_dp
      endif
      
      do q = 0, n_disp
         
         if (q /= 0) then
            ad(1,d_atom) = ad(1,d_atom) + fac*d_unit_vec(1)
            ad(2,d_atom) = ad(2,d_atom) + fac*d_unit_vec(2)
            ad(3,d_atom) = ad(3,d_atom) + fac*d_unit_vec(3)
         endif
         
         write(6,'(A,3F12.5)') 'DIFFUSING ATOM NOW AT ', &
     &        (ad(qq,d_atom),qq=1,3)
         
         iter = 1

         call relax()

         etot = eprom+ebond+epair-eatom-uent
         
         filename =  genfile(1:lengfn)//'.diff'
         if (q == 0) then
            open(43,file=filename,status='NEW')
         else
            open(43,file=filename,status='OLD')
            do while (.true.)
               read(43,'(A)',end=1) dum
            enddo
         endif
         
 1       if (q == 0) then
            eref = etot
            write(43,'(A,I5)')     '# DIFFUSING ATOM NUMBER : ', & 
     &           d_atom
            write(43,'(A,3F12.5)') '# WHICH STARTS AT : ', & 
     &           d_start(1),d_start(2),d_start(3)
            write(43,'(A,3F12.5)') '# AND FINISHS  AT : ', & 
     &           d_end(1),d_end(2),d_end(3)
            write(43,'(A,F12.5)')  '# SO THE DIFFUSION DISTANCE IS' & 
     &           ,d_dist
            write(43,'(A,I5)')     '# NUMBER OF DIFFUSION STEPS = ' & 
     &           ,n_disp    
            write(43,'(A,3F12.5)') '# DIFFUSION VECTOR AND PLANE OF & 
     &CONSTRAINT = ',(d_unit_vec(qq),qq=1,3)
            write(43,'(A,F12.5)')  '# REFERENCE ENERGY = ',eref
            write(43,'(A)') '#'
            write(43,'(A)') '#     D_R         DELTA_E      X_D  & 
     &Y_D         Z_D         FX_D        FY_D        FZ_D'

         endif
         
         write(6,'(11F12.5)') q*fac, etot - eref, ad(1,d_atom),  & 
     &        ad(2,d_atom), ad(3,d_atom), & 
     &        ftot(1,d_atom),ftot(2,d_atom),ftot(3,d_atom)
         write(43,'(11F12.5)') q*fac, etot - eref,ad(1,d_atom), & 
     &        ad(2,d_atom), ad(3,d_atom), & 
     &        ftot(1,d_atom),ftot(2,d_atom),ftot(3,d_atom)
         
         close(43)
         
      enddo                     ! DO Q
      
      
      return
      end

