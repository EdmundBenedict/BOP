    module mod_const

    use mod_precision

    implicit none

!    precision

!     integer, parameter :: dp = 4
!
!    Array dimensions.
!


!*** MXND is the maximum number of atoms allowed in the unit cell.
!      integer, parameter :: MXND = 2000
!       integer, parameter :: mxnd = 128
      integer, parameter :: mxnd = 4500

!*** KMXND is the maximum number of atoms allowed in the unit cell
!***  for k-space calculations.
      integer, parameter :: kmxnd = 32

!*** MXS is the maximum number of atoms allowed in the small unit cell.
!*** this is relevant to calculating the vibrational frequencies.
      integer, parameter :: mxs = 8

!*** MXTOTND is the maximum total number of atoms. It includes the
!    unit cell and all the other atoms needed to build the recursion
!*** clusters.
!       integer, parameter :: mxtotnd = 55000
      integer, parameter :: mxtotnd =  32*1024 !80*1024 !32*1024 !16354 !7168

!*** MXNNB is the maximum number of nearest neighbours for any given
!*** atom. Note that an atom is always a neighbour to itself.
!       integer, parameter :: mxnnb = 32, mxnnbl = 32
      integer, parameter :: mxnnb = 32, mxnnbl = 32

!*** When building moments, a cluster is built for the first moment,
!    another for the second moment etc. MXCLS is the maximum allowed
!*** sum of all these cluster sizes for one recursion chain.
!       integer, parameter :: mxcls = 5000
      integer, parameter :: mxcls = 1024 !8*1024 !8192 !2048

!*** The maximum number of allowed inert atoms.
      integer, parameter :: minert = 1000

!*** MXNSTAT is the maximum number of states assigned to any one atom.
     integer, parameter :: mxnstat = 6
!       integer, parameter :: mxnstat = 9


      integer, parameter :: kmxh = kmxnd*mxnstat ! max dimension of the tb kspace hamiltonian

!*** The maximum number of atoms allowed to form an isolated cluster.
!     integer, parameter :: MXISOCLS = 1

!*** The maximum number of allowed recursion levels.
    integer, parameter :: mrec = 4
!      integer, parameter :: mrec = 5
!       integer, parameter :: mrec = 16
!      integer, parameter :: mrec = 35


!*** The maximum number of allowed types of bond.
      integer, parameter :: matype = 2, mbtype = 4

!*** Maximum number of local densities of states that can be calculated.
!      integer, parameter :: mxldos = 100

!*** Maximum number of histogram bins.
      integer, parameter :: mbins = 100

!*** Maximum number of levels padded for a terminator.
      integer, parameter :: mterm = mrec+2

!*** Maximum length of neighbor table.
      integer, parameter :: mxnb = mxtotnd*(mxnnb+1), mxnbl = mxtotnd*(mxnnbl+1)

!*** Maximum number of angular momenta for orbitals.
      integer, parameter :: ml = 3

!*** Maximum number of recursion chains per site.
      integer, parameter :: mchain = mxnstat

!*** Maximum number of constraints.
      integer, parameter :: mcnst_n = 100

!*** Maximum number of fermi levels.
      integer, parameter :: mfermi = 100

!*** Maximum number of isolated clusters.
!     integer, parameter :: MISOL = 3

!*** Maximum number of atoms in neighbor list box.
      integer, parameter :: mcount = 500

!*** Maximum number of boxes used to generate neighbor list.
      integer, parameter :: mbox = mxtotnd

!*** Maximum number of points in binding energy curve.
      integer, parameter :: mxnt = 11

!*** Maximum order of polynomial that can be fitted to binding energy curve.
      integer, parameter :: pmax = mxnt-1

!*** Maximum number of atom types.
      integer, parameter :: mxz = 103

!*** Maximum number of roots for arctan approximation to fermi function.
      integer, parameter :: matrts = 3

!*** Maximum order of the polynomial approximation to the Fermi-Dirac function.
      integer, parameter :: mfdpol = 1

!*** The maximum dimension of the Hessian matrix. This should be 3*MXND, but
!*** to save memory it might be set smaller when not carrying out a relaxation.
      integer, parameter :: hessmx = 3*mxnd

!*** The maximum number of K points.
!      integer, parameter :: mxnk = 256

!*** The maximum number of bins for the density of states.
      integer, parameter :: mxne = 1000

!*** The maximum number of fitting parameters.
      integer, parameter :: mftpar = 100

!*** The maximum number of different displacements for a configuration.
      integer, parameter :: mlist = 100

!*** The maximum number of points in any dimension for the energy grid.
!*** This should be a power of two.
      integer, parameter :: mgrid = 16

!*** The maximum number of points in a curve to be fitted.
      integer, parameter :: mcrv = 100

!*** The maximum number of atoms in cluster used during a scan.
      integer, parameter :: s_md = 20

!*** The maximum size for the matrix recursion matrices.
      integer, parameter :: matsiz = mxnstat*mxnnb

!----------------------------------------------------

!
!    System constants.
!


!*** End of naighbor list marker.
      integer, parameter :: eol = -1000000

!*** Word size for integer.
!      integer, parameter :: intwd = 4

!*** Word size for double precision.
!      integer, parameter :: dprwd = 8

!*** Word size for double precision complex.
!      integer, parameter :: dcmwd = 2*dprwd

!*** Machine precision.
!      real(dp), parameter :: macprec = 1.0e-15_dp

!----------------------------------------------------

!
!    Program tolerances.
!


!*** Maximum error in position of minimum of binding energy curve.
!      real(dp), parameter :: xtol = 1.0e-6_dp

!*** Minimum step size for relaxation.
      real(dp), parameter :: stpmin = 1.0e-4_dp

!*** Maximum allowed inconsistency in band structure energy.
      real(dp), parameter :: maxerr = 0.01_dp

!*** The width of the Fermi-Dirac function to be approximated by
!***  a polynomial.
      real(dp), parameter :: fdx0 = 3.0_dp

      integer, private :: i,j,k
      character(len=1), parameter :: orbs(0:2) = ['s', 'p', 'd']
      character(len=3), parameter :: hnam(0:9) = [(((orbs(j)//orbs(i)//orbs(k), k = 0,j), j = 0,i), i = 0,2)]
      character(len=3), parameter :: hlnam(14) &
      & = (/ 'sss', 'sps',  'pss',  'pps',  'ppp',  'sds',  'dss',  'pds', 'dps', 'pdp', 'dpp', 'dds', 'ddp', 'ddd' /)

      character(len=1) :: crdnam(3) = ['x', 'y', 'z']

!----------------------------------------------------



      real(dp), parameter :: oppm(2) = [-1.0_dp, 1.0_dp] ! operation plus/minus


      integer, dimension(0:2,0:2,0:2) :: hlidx, hsidx


!
!    Numerical constants.
!


!*** Boltzmann's constant.
      real(dp), parameter :: kb = 8.617e-5_dp

!*** Pi.
      real(dp), parameter :: pi = 3.141592653589_dp

!
      real(dp), parameter :: sqrt3 = sqrt(3.0_dp)
      real(dp), parameter :: sqrt2 = sqrt(2.0_dp)


!----------------------------------------------------

    character(len=2), parameter :: esc = achar(27)//'['
    character(len=4), parameter :: endc = esc//'0m'  !use this thing to stop a colour

    character(len=5), parameter ::  gry = esc//'30m', &
                                 &  red = esc//'31m', &
                                 &  grn = esc//'32m', &
                                 &  yel = esc//'33m', &
                                 &  blu = esc//'34m', &
                                 &  vio = esc//'35m', &
                                 &  cyn = esc//'36m', &
                                 &  wht = esc//'37m'
    ! bold colours
    character(len=7), parameter ::  gryb = esc//'1;30m', &
                                 &  redb = esc//'1;31m', &
                                 &  grnb = esc//'1;32m', &
                                 &  yelb = esc//'1;33m', &
                                 &  blub = esc//'1;34m', &
                                 &  viob = esc//'1;35m', &
                                 &  cynb = esc//'1;36m', &
                                 &  whtb = esc//'1;37m'



    contains

        subroutine init_messy_consts()

        hlidx = huge(0)
        hlidx(0,0,0) = 1
        hlidx(0,1,0) = 2
        hlidx(1,0,0) = 3
        hlidx(1,1,0) = 4
        hlidx(1,1,1) = 5
        hlidx(0,2,0) = 6
        hlidx(2,0,0) = 7
        hlidx(1,2,0) = 8
        hlidx(2,1,0) = 9
        hlidx(1,2,1) = 10
        hlidx(2,1,1) = 11
        hlidx(2,2,0) = 12
        hlidx(2,2,1) = 13
        hlidx(2,2,2) = 14

        hsidx = huge(0)
        hsidx(0,0,0) = 0
        hsidx(0,1,0) = 1
        hsidx(1,1,0) = 2
        hsidx(1,1,1) = 3
        hsidx(0,2,0) = 4
        hsidx(1,2,0) = 5
        hsidx(1,2,1) = 6
        hsidx(2,2,0) = 7
        hsidx(2,2,1) = 8
        hsidx(2,2,2) = 9
        end subroutine init_messy_consts

    end module mod_const

