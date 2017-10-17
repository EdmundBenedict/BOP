 
!   SUBROUTINE TO CALCULATE AND WRITE OUT ATOMIC STRESS TENSORS
      subroutine katstress()
          use mod_precision


          use mod_all_scalar

          use mod_const
      implicit  none
!      include   "../Include/ALL.const"
!      include   "../Include/ALL.scalar"
      include   "../Include/Atom.array"
      include   "../Include/Force.array"
      include   "../Include/NebList.array"
      include   "../Include/PosVel.array"
!
      integer   counter, ia, ja, i, j
      real*8    f_at_stress(6,mxnd), f_tot_stress(6)
      real*8    kin_tot_stress(6)
      real*8    kin_at_stress(6,mxnd),dr(3),at_vol
      real*8    ev_gpa, kin_ev
      character*80 filename
!
      parameter(ev_gpa = 160.2189,kin_ev = 103.6413)
!
      filename = genfile(1:lengfn)//'.atstress'
      open(unit = 196,file = filename,status = 'NEW')
!     STRESS TENSOR IN ORDER xx,xy,xz,yy,yz,zz
      do 50 i = 1, mxnd
         do 40 j = 1, 6
            f_at_stress(j,i) = 0.0
            kin_at_stress(j,i) = 0.0
 40      continue
 50   continue
      do 55 j = 1, 6
         f_tot_stress(j) = 0.0
         kin_tot_stress(j) = 0.0
 55   continue
      at_vol = vol/real(nd, dp)
!
      do 70 ia = 1, nd
         ja = aptr(ia)
         counter = 1
         do while(bptr(ja) /= eol)
            dr(1) = ad(1,ia) - ad(1,bptr(ja))
            dr(2) = ad(2,ia) - ad(2,bptr(ja))
            dr(3) = ad(3,ia) - ad(3,bptr(ja))
!
            f_at_stress(1,ia) = f_at_stress(1,ia) - (dr(1) * &
     &                          (band_forces(1,counter,ia) + &
     &                           rep_forces(1,counter,ia)))/2.0

            f_at_stress(2,ia) = f_at_stress(2,ia) - (dr(1) * &
     &                          (band_forces(2,counter,ia) + &
     &                           rep_forces(2,counter,ia)))/2.0

            f_at_stress(3,ia) = f_at_stress(3,ia) - (dr(1) * &
     &                          (band_forces(3,counter,ia) + &
     &                           rep_forces(3,counter,ia)))/2.0

            f_at_stress(4,ia) = f_at_stress(4,ia) - (dr(2) * &
     &                          (band_forces(2,counter,ia) + &
     &                           rep_forces(2,counter,ia)))/2.0

            f_at_stress(5,ia) = f_at_stress(5,ia) - (dr(2) * &
     &                          (band_forces(3,counter,ia) + &
     &                           rep_forces(3,counter,ia)))/2.0

            f_at_stress(6,ia) = f_at_stress(6,ia) - (dr(3) * &
     &                          (band_forces(3,counter,ia) + &
     &                           rep_forces(3,counter,ia)))/2.0
!
            ja = ja + 1
            counter = counter + 1
         enddo
         kin_at_stress(1,ia) = kin_at_stress(1,ia) -  &
     &                          (mass(ia) * vel(1,ia) * vel(1,ia))
         kin_at_stress(2,ia) = kin_at_stress(2,ia) -  &
     &                          (mass(ia) * vel(1,ia) * vel(2,ia))
         kin_at_stress(3,ia) = kin_at_stress(3,ia) -  &
     &                          (mass(ia) * vel(1,ia) * vel(3,ia))
         kin_at_stress(4,ia) = kin_at_stress(4,ia) -  &
     &                          (mass(ia) * vel(2,ia) * vel(2,ia))
         kin_at_stress(5,ia) = kin_at_stress(5,ia) -  &
     &                          (mass(ia) * vel(2,ia) * vel(3,ia))
         kin_at_stress(6,ia) = kin_at_stress(6,ia) -  &
     &                          (mass(ia) * vel(3,ia) * vel(3,ia))
!
         do 60 j = 1, 6
            kin_at_stress(j,ia) = kin_at_stress(j,ia)*kin_ev
            f_tot_stress(j) = f_tot_stress(j) + f_at_stress(j,ia)
            kin_tot_stress(j) = kin_tot_stress(j) + kin_at_stress(j,ia)
 60      continue
 70   continue
!
      do 80 j = 1, 6
         f_tot_stress(j) = (f_tot_stress(j)/vol)*ev_gpa
         kin_tot_stress(j) = (kin_tot_stress(j)/vol)*ev_gpa
 80   continue
!
      do 90 ia = 1, nd
         do 85, j = 1, 6
            f_at_stress(j,ia) = (f_at_stress(j,ia)/at_vol)*ev_gpa
            kin_at_stress(j,ia) = (kin_at_stress(j,ia)/at_vol)*ev_gpa
 85      continue
 90   continue
!
      write(196,110) nd, (f_tot_stress(j) + kin_tot_stress(j),j=1,6)
      do 100 ia = 1, nd
         write(196,120) ia, listsymb(z(ia)), &
     &                (f_at_stress(j,ia) + kin_at_stress(j,ia),j=1,6)
 100  continue
 110  format(i5,2x,6e15.7)
 120  format(i5,2x,a2,2x,6e15.7)
      return
      end
