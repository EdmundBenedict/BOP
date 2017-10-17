 
       function geteatom(z,nd)
          use mod_precision


!
!    This is a function to evaluate the bandstructure energies
!     of the sum of isolated atoms.
!

      implicit none
      real(dp) :: geteatom

!
!    Declare the simple variables.
!

      real(dp) :: bseatm

      integer nd,i,zi
      integer mxz,ml

!
!    Define the parameters.
!

      parameter(mxz = 103)
      parameter(ml =3)

!
!    Declare the arrays.
!

      real(dp) :: zc(0:mxz)
      real(dp) :: es(0:mxz)
      real(dp) :: ep(0:mxz)
      real(dp) :: ed(0:mxz)

      integer nl(0:mxz)
      integer nstt(0:mxz)
      integer llist(ml,0:mxz)
      integer z(nd)

!
!    Declare the common blocks.
!

      common /zcom/zc
      common /lcom/nl,nstt,llist
      common /ecom/es,ep,ed

!
!    Return the value of the energy.
!
      geteatom = 0.0_dp
      do i = 1,nd
         zi = z(i)
         geteatom = geteatom + bseatm(es(zi),ep(zi),ed(zi), nl(zi),llist(1,zi),zc(zi),zi)
      enddo
      end

