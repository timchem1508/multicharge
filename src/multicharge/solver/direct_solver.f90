module direct_solver
    use mctc_env, only: error_type, fatal_error, wp
    use multicharge_blas, only: symv
    use multicharge_lapack, only: sytrf, sytrs, sytri
    use solver_type_cache, only: cache_container
    use solver, only: mchrg_solver_type
    implicit none
    private

    public :: direct_solver_type

    !> Direct solver using LAPACK
    type, extends(mchrg_solver_type) :: direct_solver_type
    contains
       procedure :: solve => solve_direct
       procedure :: update => update_direct
    end type direct_solver_type

contains

    !> Update method for direct solver
    subroutine update_direct(self, cache, amat, xvec, vrhs, ainv, cpq, info)
        class(direct_solver_type), intent(in) :: self
        type(cache_container), intent(inout) :: cache
        real(wp), intent(in)  :: amat(:, :)
        real(wp), intent(in)  :: xvec(:)
        real(wp), allocatable, intent(out) :: vrhs(:)
        real(wp), allocatable, intent(out), optional :: ainv(:, :)
        logical, intent(in), optional :: cpq
        integer, intent(out), optional :: info
        
        integer :: ndim
        
        ndim = size(xvec)
        
        ! Allocate and initialize vrhs
        allocate(vrhs(ndim))
        vrhs = xvec
        
        ! Allocate and initialize ainv if present
        if (present(ainv)) then
            allocate(ainv(ndim, ndim))
            ainv = amat
        end if
        
        if (present(info)) info = 0
    end subroutine update_direct

    !> Solve method for direct solver
    subroutine solve_direct(self, amat, xvec, vrhs, ainv, cpq, error, info)
        class(direct_solver_type), intent(in) :: self
        type(error_type), allocatable, intent(out) :: error
        real(wp), intent(in)  :: amat(:, :)
        real(wp), intent(in)  :: xvec(:)
        real(wp), allocatable, intent(out) :: vrhs(:)
        real(wp), allocatable, intent(out), optional :: ainv(:, :)
        logical, intent(in), optional :: cpq
        integer, intent(out), optional :: info
    
        integer  :: local_info
        integer :: ndim, ic, jc
        integer, allocatable :: ipiv(:)
        logical :: want_cpq
        type(cache_container), allocatable :: cache
    
        ! Dimensions match check
        ndim = size(xvec)
        if (size(amat,1) /= ndim .or. size(amat,2) /= ndim) then
            call fatal_error(error, "solve_direct: dimension mismatch.")
            if (present(info)) info = -1 
            return
        end if
        
        ! Update cache and prepare vrhs and ainv
        allocate(cache)
        call self%update(cache, amat, xvec, vrhs, ainv, cpq, info)
    
        ! Logical: solve coupled-perturbed equations flag
        want_cpq = .false.
        if (present(cpq)) want_cpq = cpq
    
        allocate(ipiv(ndim))
        call sytrf(ainv, ipiv, info=local_info, uplo='l')
        if (local_info /= 0) then
            call fatal_error(error, "solve_direct: Bunch-Kaufman factorization failed.")
            if (present(info)) info = local_info
            return
        end if
    
        if (want_cpq) then
            call sytri(ainv, ipiv, info=local_info, uplo='l')
            if (local_info /= 0) then
                call fatal_error(error, "solve_direct: Inversion of factorized matrix failed.")
                if (present(info)) info = local_info
                return
            end if
            call symv(ainv, xvec, vrhs, uplo='l')
            do ic = 1, ndim
                do jc = ic + 1, ndim
                    ainv(ic, jc) = ainv(jc, ic)
                end do
            end do
        else
            call sytrs(ainv, vrhs, ipiv, info=local_info, uplo='l')
            if (local_info /= 0) then
                call fatal_error(error, "solve_direct: Solving factorized system failed.")
                if (present(info)) info = local_info
                return
            end if
        end if
    
        if (present(info)) info = local_info
    end subroutine solve_direct

end module direct_solver