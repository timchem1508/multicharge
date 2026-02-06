module direct_solver
    use mctc_env, only: error_type, fatal_error, wp
    use multicharge_blas, only: symv
    use multicharge_lapack, only: sytrf, sytrs, sytri
    use solver_type_cache, only: cache_container
    use print_matrix, only: write_vector, write_matrix
    use solver, only: mchrg_solver_type
    implicit none
    private

    public :: direct_solver_type

    !> Direct solver using LAPACK
    type, extends(mchrg_solver_type) :: direct_solver_type
    contains
       procedure :: solve
       procedure :: update 
    end type direct_solver_type

contains

    !> Update method for direct solver
    subroutine update(self, cache, vrhs, ainv, cpq)
        class(direct_solver_type), intent(in) :: self
        type(cache_container), intent(inout) :: cache
        real(wp), intent(inout) :: vrhs(:)
        real(wp), intent(out) :: ainv(:, :)
        logical, intent(in), optional :: cpq
        
    end subroutine update

    !> Solve method for direct solver
    subroutine solve(self, amat, xvec, vrhs, ainv, cpq, error)
        class(direct_solver_type), intent(in) :: self
        real(wp), intent(in)  :: amat(:, :)
        real(wp), intent(in)  :: xvec(:)
        real(wp), intent(inout) :: vrhs(:)
        real(wp), intent(out) :: ainv(:, :)
        logical, intent(in), optional :: cpq
        type(error_type), allocatable, intent(out) :: error
    
        integer  :: local_info
        integer :: ndim, ic, jc
        integer, allocatable :: ipiv(:)
        logical :: want_cpq
        type(cache_container), allocatable :: cache
    
        ! Dimensions match check
        ndim = size(xvec)
        if (size(amat,1) /= ndim .or. size(amat,2) /= ndim) then
            call fatal_error(error, "solve_direct: dimension mismatch.")
            return
        end if

        ainv = amat
        vrhs = xvec
        
        ! Update cache and prepare vrhs and ainv
        allocate(cache)
        call self%update(cache, vrhs, ainv, cpq)
    
        ! Logical: solve coupled-perturbed equations flag
        want_cpq = .false.
        if (present(cpq)) want_cpq = cpq
    
        allocate(ipiv(ndim))
        call sytrf(ainv, ipiv, info=local_info, uplo='l')
        if (local_info /= 0) then
            call fatal_error(error, "solve_direct: Bunch-Kaufman factorization failed.")
            return
        end if
    
        if (want_cpq) then
            call sytri(ainv, ipiv, info=local_info, uplo='l')
            if (local_info /= 0) then
                call fatal_error(error, "solve_direct: Inversion of factorized matrix failed.")
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
                return
            end if
        end if
    
    end subroutine solve

end module direct_solver