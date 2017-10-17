 
       function getmass(z)
          use mod_precision


!
!    This is a function to return the mass of an ion core
!     in atomic mass units.
!

      implicit none
      real(dp) :: getmass

!
!    Declare the simple variables.
!

      integer z

!
!    Declare the arrays.
!

      real(dp) :: mass(103)

!
!    Assign values to the mass.
!

      data mass/  1.0_dp,  4.0_dp,  6.9_dp,  9.0_dp, 10.8_dp, 12.0_dp, & 
     &           14.0_dp, 16.0_dp, 19.0_dp, 20.2_dp, 23.0_dp, 24.3_dp, & 
     &           27.0_dp, 28.1_dp, 31.0_dp, 32.1_dp, 35.5_dp, 39.9_dp, & 
     &           39.1_dp, 40.1_dp, 45.0_dp, 47.9_dp, 50.9_dp, 52.0_dp, & 
     &           54.9_dp, 55.9_dp, 58.9_dp, 58.7_dp, 63.6_dp, 65.4_dp, & 
     &           69.7_dp, 72.6_dp, 74.9_dp, 79.0_dp, 79.9_dp, 83.8_dp, & 
     &           85.5_dp, 87.6_dp, 88.9_dp, 91.2_dp, 92.9_dp, 95.9_dp, & 
     &           98.9_dp,101.1_dp,102.9_dp,106.4_dp,107.9_dp,112.4_dp, & 
     &          114.8_dp,118.7_dp,121.8_dp,127.6_dp,126.9_dp,131.3_dp, & 
     &          132.9_dp,137.3_dp,138.9_dp,140.1_dp,140.9_dp,144.2_dp, & 
     &          145.0_dp,150.4_dp,152.0_dp,157.3_dp,158.9_dp,162.5_dp, & 
     &          164.9_dp,167.3_dp,168.9_dp,173.0_dp,175.0_dp,178.5_dp, & 
     &          181.0_dp,183.9_dp,186.2_dp,190.2_dp,192.2_dp,195.1_dp, & 
     &          197.0_dp,200.6_dp,204.4_dp,207.2_dp,209.0_dp,210.0_dp, & 
     &          210.0_dp,222.0_dp,223.0_dp,226.0_dp,227.0_dp,232.0_dp, & 
     &          231.0_dp,238.0_dp,237.1_dp,244.0_dp,243.0_dp,247.0_dp, & 
     &          247.0_dp,251.0_dp,254.0_dp,257.0_dp,256.0_dp,254.0_dp, & 
     &          257.0_dp/

!
!    Return the value of the mass.
!

      if (z == 0) then
         getmass = 28.1_dp ! Silicon for some reason
      else
         getmass = mass(z)
      endif

      end

