 
      subroutine bldlcell()
          use mod_precision
          use mod_all_scalar
          use mod_const
          use mod_atom_ar
          
!** This, together with getvib has to be reviewed and either removed or adapted to the parallel distribution and dynamic arrays          
!   
                    
          
      implicit none
!
!      include "Include/ALL.const"
!      include "Include/ALL.scalar"
      include "Include/Atom.array"
      include "Include/PosVel.array"
      include "Include/NebList.array"
      include "Include/Vib.array"
!
      integer i,i1,i2,i3,ia,ix,iy,iz,ntot
!
!    routine to calculate large cell for phonon calculation
!    note - a tetragonal cell is assumed here!!
!    In the mapping of atom to lattice vector the average
!    coordinate is assumed to be 0.5,0.5,0.5
!
      if(nper(1)*nper(2)*nper(3)*nd > mxnd)then
         write(6,'(''NEED TO INCREASE MXND FOR LARGE CELL'')')
         write(6,'(''STOPPING BOP!!!'')')
         stop
      endif
!
      lena(1) = sqrt(a(1,1)**2+a(1,2)**2+a(1,3)**2)
      lena(2) = sqrt(a(2,1)**2+a(2,2)**2+a(2,3)**2)
      lena(3) = sqrt(a(3,1)**2+a(3,2)**2+a(3,3)**2)
!     


      ntot=0
      do i1=1,nper(1)
         do i2=1,nper(2)
            do i3=1,nper(3)
               do ia=1,nd
                  ntot=ntot+1
                  ad(1,ntot) = (i1-1)*a(1,1) + (i2-1)*a(2,1) + (i3-1)*a(3,1) + ad(1,ia)
                  ad(2,ntot) = (i1-1)*a(1,2) + (i2-1)*a(2,2) + (i3-1)*a(3,2) + ad(2,ia)
                  ad(3,ntot) = (i1-1)*a(1,3) + (i2-1)*a(2,3) + (i3-1)*a(3,3) + ad(3,ia)
                  z(ntot)=z(ia)
                  mapx(ntot)=ia
                  de(ntot)=de(ia)
                  dq(ntot)=dq(ia)
!     
                  if(mod(nper(1),2) == 1)then
                     if((i1-1) <= nper(1)/2)then
                        ix=i1-1
                     else
                        ix=(i1-1)-nper(1)
                     endif
                  elseif(mod(nper(1),2) == 0)then
                     if((i1-1) < nper(1)/2)then
                        ix=i1-1
                     elseif((i1-1) > nper(1)/2)then
                        ix=(i1-1)-nper(1)
                     elseif(ad(1,ia)/lena(1) < 0.5_dp)then
                        ix=i1-1
                     else
                        ix=(i1-1)-nper(1)
                     endif
                  endif
                  if(mod(nper(2),2) == 1)then
                     if((i2-1) <= nper(2)/2)then
                        iy=i2-1
                     else
                        iy=(i2-1)-nper(2)
                     endif
                  elseif(mod(nper(2),2) == 0)then
                     if((i2-1) < nper(2)/2)then
                        iy=i2-1
                     elseif((i2-1) > nper(2)/2)then
                        iy=(i2-1)-nper(2)
                     elseif(ad(2,ia)/lena(2) < 0.5_dp)then
                        iy=i2-1
                     else
                        iy=(i2-1)-nper(2)
                     endif
                  endif
                  if(mod(nper(3),2) == 1)then
                     if((i3-1) <= nper(3)/2)then
                        iz=i3-1
                     else
                        iz=(i3-1)-nper(3)
                     endif
                  elseif(mod(nper(3),2) == 0)then
                     if((i3-1) < nper(3)/2)then
                        iz=i3-1
                     elseif((i3-1) > nper(3)/2)then
                        iz=(i3-1)-nper(3)
                     elseif(ad(3,ia)/lena(3) < 0.5_dp)then
                        iz=i3-1
                     else
                        iz=(i3-1)-nper(3)
                     endif
                  endif
                  
                  mapr(1,ntot)=ix
                  mapr(2,ntot)=iy
                  mapr(3,ntot)=iz
                  
               enddo
            enddo
         enddo
      enddo
!
      nd=ntot
      a(1,1)=a(1,1)*real(nper(1), dp)
      a(1,2)=a(1,2)*real(nper(1), dp)
      a(1,3)=a(1,3)*real(nper(1), dp)
      a(2,1)=a(2,1)*real(nper(2), dp)
      a(2,2)=a(2,2)*real(nper(2), dp)
      a(2,3)=a(2,3)*real(nper(2), dp)
      a(3,1)=a(3,1)*real(nper(3), dp)
      a(3,2)=a(3,2)*real(nper(3), dp)
      a(3,3)=a(3,3)*real(nper(3), dp)

      lena(1) = sqrt(a(1,1)**2+a(1,2)**2+a(1,3)**2)
      lena(2) = sqrt(a(2,1)**2+a(2,2)**2+a(2,3)**2)
      lena(3) = sqrt(a(3,1)**2+a(3,2)**2+a(3,3)**2)

      eatom=eatom*nper(1)*nper(2)*nper(3)
!
!    Rescale RLIM according to RPRUNE.
!
      do i = 1,3,1
         if (rlim(i) /= 0) then
            rlim(i) = int(rprune/lena(i))+1
         endif
      enddo
!     
      return
      end
      
