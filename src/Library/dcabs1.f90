 
       function dcabs1(zr,zi)
          use mod_precision

           real(dp) :: dcabs1
      real(dp) :: zr,zi
      dcabs1 = abs(zr) + abs(zi)
      end
      subroutine cpy(tr,ti,dr,di,t1,t2)
          use mod_precision

!- complex multiply (t1,t2) = (tr,ti) * (dr,di) 
      real(dp) :: tr,ti,dr,di,t1,t2      
      real(dp) :: tmp
      tmp = tr*dr - ti*di
      t2  = tr*di + ti*dr
      t1 = tmp
      end
