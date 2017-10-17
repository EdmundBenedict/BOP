
#ifdef MPI         
    module mod_bandfit_mpi
    use mpi
#else
    module mod_bandfit
#endif    
    use mod_precision
    use mod_all_scalar
    use mod_const
    use mod_conf
    use topologia, only : iproc, mpmap
    use mod_atom_ar, only : btype
    implicit none
    
    integer, allocatable, private :: pidx(:,:) ! hops, btype 
    real(dp), allocatable, private :: p0(:,:) ! hops, btype 
    
    
    
    contains
    
    subroutine init_pars(par, npar)
        integer, intent(out) :: npar
        real(dp), intent(out), allocatable :: par(:)
        integer :: i, j, k, oj, lj, oi, li, bi, nb , idx, hidx
        type(bond_conf_t), pointer :: b
        real(dp), allocatable :: tmp_par(:)
        
        idx = 0
        nb = rtc%ham_conf%nb - 1
        allocate(pidx(0:9, 0:nb), p0(0:9, 0:nb), tmp_par(10*(nb+1)))

        pidx = 0
        
        do i = 0, nb
            b => rtc % ham_conf % b(i)
            do oj = 1, b % atom1 % b % norbs
                lj = b % atom1 % b % orblist(oj)
                do oi = 1, b % atom0 % b % norbs
                    li = b % atom0 % b % orblist(oi)
                    if (li > lj .or. (b % atom0 % z > b % atom1 % z .and. li == lj)) cycle
!                     print *, b% name,' ', orbs(li), orbs(lj)
                    do bi = 0, li ! <=
                        hidx = hsidx(li,lj,bi)
                        if (b % h(hidx) % fit) then
                            idx = idx + 1
!                             print *, b % name,' ', hnam(hidx), b % h(hidx) % v
                            tmp_par(idx) = b % h(hidx) % v
                            pidx(hidx,i) = idx
                            if (b % atom0 % z < b % atom1 % z .and. li == lj) pidx(hidx, b % rev % id) = idx
                        end if
                    end do
                end do
            end do 
        end do
        
        npar = idx
        allocate(par(npar))
        par = tmp_par(:npar)
        deallocate(tmp_par)

        do i = 0, nb
            do j = 0,9
                p0(j,i) = rtc % ham_conf % b(i) % h(j) % v
            end do
        end do
        
        
        print *, par(:npar)
        do i = 0, nb
            print '(10(x,i4))', pidx(:,i)
        end do
        

        
!         stop        
    end subroutine init_pars
    
    subroutine eval_evals(par, npar,  evals)
        use mod_kspace
        integer :: npar
        real(dp), intent(in) :: par(npar)
        real(dp), intent(out) :: evals
    
        qerr = -1.0_dp
     
        call kebsfbs()
        
        

    
        
    
    end subroutine eval_evals
    
#ifdef MPI         
    end module mod_bandfit_mpi
#else
    end module mod_bandfit
#endif 











    subroutine bandfit()

#ifdef MPI         
    use mod_bandfit_mpi
#else
    use mod_bandfit
#endif
    
    implicit none
    
    real(dp), allocatable :: par(:)
    integer :: npar
    
    call init_pars(par, npar)
    
    call bldnebt()
    
    
    
! 
!     
!     print *, enk(:,mpmap(iproc)+1)
!     
    stop
!     
    
    end subroutine bandfit
    
    
    
    