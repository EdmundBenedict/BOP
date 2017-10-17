    module mod_all_scalar
      use mod_const
      implicit none
!*** Timer values.
      real(dp) time
!*** Fermi energy, electron temperature, number of electrons
      real(dp) :: lef,kt,locc
!*** Promotion, bond, repulsive, atomic, entropic, kinetic  and total energy
      real(dp) :: etot, ebind, eelct, eband, eprom, ebond, eatom, uent, ke, eclas, epair, eenv, emag
      real(dp) :: eproms, epromp, epromd
!*** Energy of isolated cluster.
!     real(dp) :: ECLUSTERS
!*** Atomic and total cluster cutoffs.
      real(dp) :: rcut,rprune,rcutl
!*** Average cluster size, and average number of neighbors.
      real(dp) :: aveclusiz
!*** Allowed tolerance on ionic charge.
      real(dp) :: qerr
!*** The MD step size, temperature, volume, pressure and piston mass.
      real(dp) :: dt,temp,vol,press,mpiston
!*** The contribution to the pressure from the pair potential and band structure force.
      real(dp) :: rdfpp, rdfbs
!*** The temperature limits and step size for simulated annealing.
      real(dp) :: tmin,tmax,dtemp
!*** The force tolerance and step size for atomic relaxation.
      real(dp) :: ftol,step
!*** The sum, sum of squares, mean and standard deviation of the
!*** energy, pressure, temperature and volume.
      real(dp) :: sume,sump,sumt,sumv
      real(dp) :: sumee,sumpp,sumtt,sumvv
      real(dp) :: avge,avgp,avgt,avgv
      real(dp) :: sde,sdp,sdt,sdv
!*** Range and interval for pair correlation function.
      real(dp) :: mda_gr_r,mda_gr_dr
!*** Range for angular correlation function.
      real(dp) :: mda_g3_r
!*** Range and interval for Stillinger-Weber three body energy distribution.
      real(dp) :: mda_sw_e,mda_sw_de
!*** Electron charge distribution energy range.
!     real(dp) :: XA_RHO_LO,XA_RHO_HI
!*** The range of displacement for binding energy curves.
      real(dp) :: tlo,thi
!*** Diffusion cluster size.
      real(dp) :: rc_diff
!*** The rate of change of number of electrons with chemical potential.
      real(dp) :: dndm
!*** The total number of electrons.
      real(dp) :: totnia
!*** The maximum excess charge.
      real(dp) :: dqmax
! !*** The energy shift parameters.
!       real(dp) :: sumz,sumzde,aa
!*** The tolerance in temperature.
      real(dp) :: ttol
!*** Fitting tolerance for z dependance of energy surface.
      real(dp) :: s_tol
!*** Separation between planes, and external force on upper plane.
      real(dp) :: dzpln,fxtpln,vxpln,vypln,fplnx,fplny,fplnz
!*** The radius and external pressure for the Stoneham bomb.
      real(dp) :: bombr0,bombpx,bombp
!*** The extra width added to the density of states.
      real(dp) :: etail
!*** The onsite energy
      real(dp) :: eos
!*** Broadening for matrix recursion terminator.
      real(dp) :: mbinf
!*** The band limits used by matrix recursion.
      real(dp) :: melo,mehi
!*** Thermostat relaxation time.
      real(dp) :: tau
!*** The prefactor for calculating the number of poles for Alex's
!*** integration scheme
      integer :: mfac
!*** Maximum number of iterations.
      integer :: mxiter
!*** Order of polynomial and number of points for binding energy
!*** curve fit.
      integer :: polyord,nt
!*** Number of inert atoms.
      integer :: ninert
!*** Flags and number of bins for pair correlation function, angular
!*** correlation function, and Stillinger-Weber three body energy distribtion.
!*** Plus time-averaged atomic positions and order parameter.
      integer :: mda_gr,mda_gr_n
      integer :: mda_g3,mda_g3_n
      integer :: mda_sw,mda_sw_n
      integer :: mda_adav, mda_op
!*** Flags for mean squared displacement calculation.
      integer :: mda_msd, msd_atom
!*** Number of equilibration steps, monitor output frequency,
!*** configuration write frequency and number of iterations for molecular dynamics.
      integer :: nequil,monitorper,writeper,iter, datform
!*** The autosave frequency, and the trace flag.
      integer :: autosave,trace
!*** The number of recursion levels, and number of levels used to
!*** define cluster size.
      integer :: nrec,nbase
!*** The total number of atoms in the whole cluster and the unit cell.
      integer :: totnd,nd
!*** The flag for the kind of atomic movement.
      integer :: mvflag
!*** The kind of terminator and moments used.
      integer :: term,momflg
!*** The number of fermi energies.
      integer :: nfermi
!*** The number of constraints.
      integer :: cnst_n
!*** The number of isolated clusters.
!     integer :: NISOL
!*** The length of the file name
      integer :: lengfn
!*** The order of polynomials for Masato's integration scheme.
      integer :: pdord,pd1ord
!*** The density of states flag and number of bins.
!     integer :: XA_DOS,XA_DOS_N
!*** The local density of states flag, number of bins, and number of states.
!     integer :: XA_LDOS,XA_LDOS_NE,XA_LDOS_NS
!*** The charge density flag.
!     integer :: XA_RHO
!*** The energy vs. occupancy flag.
!     integer :: XA_EVSOCC
!*** The calibration flag and atom.
      integer :: calibrate,calibatom
!*** The diffusion atom, and number of steps.
      integer :: d_atom,n_disp
!*** The number of atom and bond types.
      integer :: natype,nbtype
!*** The one body potential flag. If this is set to 1, then the one body potential is included.
      integer :: v1flag
!!*** The size of the K-space Hamiltonian.
!      integer :: kmxh
! !*** The number of K points
!       integer :: nk
!*** Flag indicating the k-point basis.
      integer :: kbase
!*** The number of fitting parameters.
      integer :: nftpar
!*** The debug verbosity flag.
      integer :: idebug
!*** The number of displacements for a given configuration.
      integer :: nlist
!*** The number of atoms in cluster to be used for the energy scan.
      integer :: s_nd
!*** Flag indicating whether FFT of the energy in a scan is to be evaluated.
      integer :: s_fft
!*** The number of points in a curve that is to be fitted
      integer :: ncrv
!*** The type of relaxation algorithm used. 1 => Variable Metric
!***                                        2 => Steepest Descent
      integer :: rlxflg
!*** The type of method used to calculate repsonse functions.
!***        1 => Masato's analytic method.
!***        2 => AMB's complex plane integration.
      integer :: chi_meth
!*** The zapper flag and atom type.
      integer :: zapflg,zapz
!*** The size of the matrix recursion matrices.
      integer :: matn
!*** The length of the matrix recursion chain
      integer :: mlchain
!*** XFRCFLG is set to 1 if extra forces are to be used.
      integer :: xfrcflg

!*** POTFLG is the flag that determines how the energy is to be
!***  evaluated.
      integer :: potflg
      integer, parameter :: mixed=-1, mtbeam=0, tersoff=1, sbop=2, clust=3, &
     &  kspace=4, glob = 5,  foe = 7, densel = 8, tbbop = 9

!*** NBLOCK is the number of blocks.
      integer :: nblock


!*** The generic file name.
      character(len=80) genfile

!*** Verbosity *** D.P. *** to be changed to a variety of levels later
      logical :: quiet = .false.

! only rebuild the diagonal part of the hamiltonian
      logical :: onsonly = .false.

      logical :: forces = .false.

! magnetic? a convenience switch saves typing nsp == 1
      logical :: mag = .false.

! number of spins
      integer :: nsp = 1

! maximum magnetic mome      
      real(dp) :: mgmax
      
! radius of i, used in strain
      real(dp) :: radiusi
      
! occupancy tolerance for the fermi level
      real(dp) :: eftol
      
      
             
    end module mod_all_scalar
