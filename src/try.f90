 
      character function getname(z)
          use mod_precision


!
!    This is a function to return the chemical symbol corresponding
!     to an atomic number.
!

      implicit none

      integer z
      character listname(103)

!
!    Assign data values.
!

      data listname/'H','He','Li','Be','B','C','N','O','F','Ne', & 
     &              'Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca', & 
     &              'Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn', & 
     &              'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr', & 
     &              'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn', & 
     &              'Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd', & 
     &              'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb', & 
     &              'Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg', & 
     &              'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th', & 
     &              'Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm', & 
     &              'Md','No','Lw'/


      getname = listname(z)
      if ((z < 0) .or. (z > 103)) then
            write(6,'(I3,'' is not yet in the periodic table.'')') z
!            CALL PANIC()
      endif

      end
