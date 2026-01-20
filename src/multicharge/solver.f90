module solver
    use mctc_env, only: error_type, fatal_error, wp
    use multicharge_blas, only: symv, gemv
    use multicharge_lapack, only: sytrf, sytrs, sytri
    use multicharge_model_cache, only: model_cache, cache_container
    use print_matrix, only: write_matrix, write_vector
    implicit none
    private

    public :: mchrg_solver_type, new_mchrg_solver

    type, abstract :: mchrg_solver_type
    contains
       procedure(solve_if), deferred :: solve
       procedure(update_if), deferred :: update
    end type mchrg_solver_type

    abstract interface
        subroutine solve_if(self, amat, xvec, vrhs, ainv, cpq, error, info)
            import :: mchrg_solver_type, error_type, wp
            class(mchrg_solver_type), intent(in) :: self
            real(wp), intent(in)  :: amat(:, :)
            real(wp), intent(in)  :: xvec(:)
            real(wp), allocatable, intent(out) :: vrhs(:)
            real(wp), allocatable, intent(out), optional :: ainv(:, :)
            logical, intent(in), optional :: cpq
            type(error_type), allocatable, intent(out) :: error
            integer , intent(out), optional :: info
        end subroutine solve_if

        subroutine update_if(self, cache, amat, xvec, vrhs, ainv, cpq, info)
            import :: mchrg_solver_type, cache_container, wp
            class(mchrg_solver_type), intent(in) :: self
            type(cache_container), intent(inout) :: cache
            real(wp), intent(in)  :: amat(:, :)
            real(wp), intent(in)  :: xvec(:)
            real(wp), allocatable, intent(out) :: vrhs(:)
            real(wp), allocatable, intent(out), optional :: ainv(:, :)
            logical, intent(in), optional :: cpq
            integer , intent(out), optional :: info
        end subroutine update_if  

    end interface

    ! Direct solver using LAPACK
    type, extends(mchrg_solver_type) :: direct_solver_type
    contains
       procedure :: solve => solve_direct
       procedure :: update => update_direct
    end type direct_solver_type

    ! CG solver with Jacobi preconditioner
    type, extends(mchrg_solver_type) :: cg_solver_type
!     integer :: max_iter = 1000
 !    real(wp) :: tol = 1.0e-8_wp
    contains
       procedure :: solve => solve_cg
       procedure :: update => update_cg
    end type cg_solver_type

contains

! Update method for direct solver (does nothing in this implementation)
subroutine update_direct(self, cache, amat, xvec, vrhs, ainv, cpq, info)
    class(direct_solver_type), intent(in) :: self
    type(cache_container), intent(inout) :: cache
    real(wp), intent(in)  :: amat(:, :)
    real(wp), intent(in)  :: xvec(:)
    real(wp), allocatable, intent(out) :: vrhs(:)
    real(wp), allocatable, intent(out), optional :: ainv(:, :)
    logical, intent(in), optional :: cpq
    integer , intent(out), optional :: info
    
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

! Update method for CG solver (does nothing in this implementation)
subroutine update_cg(self, cache, amat, xvec, vrhs, ainv, cpq, info)
    class(cg_solver_type), intent(in) :: self
    type(cache_container), intent(inout) :: cache
    real(wp), intent(in)  :: amat(:, :)
    real(wp), intent(in)  :: xvec(:)
    real(wp), allocatable, intent(out) :: vrhs(:)
    real(wp), allocatable, intent(out), optional :: ainv(:, :)
    logical, intent(in), optional :: cpq
    integer , intent(out), optional :: info
    
    integer :: ndim
    
    ndim = size(xvec)
    
    ! Allocate and initialize vrhs (initial guess for CG is zero)
    allocate(vrhs(ndim))
    vrhs = 0.0_wp ! we can start from zero guess, then switch to the local charges
    
    ! Allocate ainv if present (though CG doesn't use it)
    if (present(ainv)) then
        allocate(ainv(ndim, ndim))
    end if
    
    if (present(info)) info = 0
end subroutine update_cg

subroutine solve_direct(self, amat, xvec, vrhs, ainv, cpq, error, info)
    class(direct_solver_type), intent(in) :: self
    type(error_type), allocatable, intent(out) :: error
    real(wp), intent(in)  :: amat(:, :)
    real(wp), intent(in)  :: xvec(:)
    real(wp), allocatable, intent(out) :: vrhs(:)
    real(wp), allocatable, intent(out), optional :: ainv(:, :)
    logical, intent(in), optional :: cpq
    integer , intent(out), optional :: info

    integer  :: local_info
    integer :: ndim, ic, jc
    integer , allocatable :: ipiv(:)
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
        call sytri(ainv, ipiv, info=local_info, uplo='l') !ipiv into a cache?
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

subroutine solve_cg(self, amat, xvec, vrhs, ainv, cpq, error, info)
    class(cg_solver_type), intent(in) :: self
    type(error_type), allocatable, intent(out) :: error
    real(wp), intent(in)  :: amat(:, :)
    real(wp), intent(in)  :: xvec(:)
    real(wp), allocatable, intent(out) :: vrhs(:)
    real(wp), intent(out), allocatable, optional :: ainv(:, :) ! not used in CG
    logical, intent(in), optional :: cpq
    integer , intent(out), optional :: info
    
    integer :: ndim, it, maxit
    real(wp) :: tol, bnorm, rnorm, alpha, beta, denom
    real(wp), allocatable :: r(:), p(:), z(:), Ap(:), Mdiag(:)
    integer  :: local_info
    real(wp) :: rz_old, rz_new
    type(cache_container), allocatable :: cache
    
    ! Dimensions match check
    ndim = size(xvec)
    if (size(amat,1) /= ndim .or. size(amat,2) /= ndim) then  ! REMOVED: .or. size(vrhs) /= ndim
        call fatal_error(error, "solve_cg: dimension mismatch.")
        if (present(info)) info = -1 
        return
    end if

    ! Prepare/cache - update will allocate and initialize vrhs
    allocate(cache)
    call self%update(cache, amat, xvec, vrhs, ainv, cpq, info)
 
    ! Global thresholds
    tol = 1.0e-11_wp
    maxit = max(10, ndim*20)

    allocate(r(ndim), p(ndim), z(ndim), Ap(ndim), Mdiag(ndim))

    ! Jacobi preconditioner (inverse of diagonal)
    do it = 1, ndim
        Mdiag(it)=amat(it,it)
        if (abs(Mdiag(it)) < tol**2) Mdiag(it) = tol**2
        Mdiag(it) = 1.0_wp / Mdiag(it)
    end do

    ! Initial residual r = b - A*x (x=vrhs, which is 0.0_wp from update_cg)
    call symv(amat, vrhs, Ap, alpha=1.0_wp, beta=0.0_wp)
    ! Residual compute
    r = xvec - Ap
    ! Apply preconditioner z = M * r
    z = r * Mdiag
    ! Initial search direction
    p = z
    ! Initial direction update factor
    bnorm = sum(xvec*xvec)
    if (bnorm < tol**2) bnorm = 1.0_wp
    rnorm = sum(r*r)
    if (rnorm / bnorm <= tol**2) then
        if (present(info)) info = 0 
        return
    end if

    ! Dynamical residual
    rz_old = sum(r*z)
    local_info = -1 

    ! Conjugate Gradient iterations
    do it = 1, maxit
        call gemv(amat, p, Ap, alpha=1.0_wp, beta=0.0_wp, trans='n')
        ! Compute step size alpha
        denom = sum(p * Ap) + tiny(1.0_wp)
        write(*,*) "Iteration ", it, " Residual norm: ", sqrt(rnorm / bnorm)
        if (abs(denom) < tol**2) then
            local_info = 0 
            exit
        end if
        alpha = rz_old / denom
        write(*,*) " Alpha: ", alpha
        ! Update solution and residual
        vrhs = vrhs + alpha * p
        write(*,*) " Max abs(vrhs): ", maxval(abs(vrhs))
        r = r - alpha * Ap
        write(*,*) " Max abs(residual): ", maxval(abs(r))
        ! Check convergence
        rnorm = sqrt(sum(r*r))
        if (rnorm / bnorm <= tol) then
            write(*,*) "CG converged in ", it, " iterations."
            call write_vector(vrhs, "CG Solution Vector")
            local_info = 0 
            exit
        end if
        ! Apply preconditioner z = M * r
        z = r * Mdiag
        rz_new = sum(r * z)
        beta = rz_new / rz_old
        p = z + beta * p
        rz_old = rz_new
        if (it == maxit) then
            local_info = 1   ! did not converge
            call fatal_error(error, "solve_cg: CG did not converge within max iterations.")
        end if
    end do

    if (present(info)) info = local_info
end subroutine solve_cg


function new_mchrg_solver(use_cg) result(solver)
    logical, intent(in), optional :: use_cg
    class(mchrg_solver_type), allocatable :: solver
    logical :: cg
    character(len=32) :: env

    cg = .true.
    if (present(use_cg)) then
       cg = use_cg
    else
       call get_environment_variable("MCHARGE_SOLVER", env)
       if (trim(env) == "DIRECT") cg = .false.
    end if

    if (cg) then
       allocate(cg_solver_type :: solver)
    else
       allocate(direct_solver_type :: solver)
    end if
 end function new_mchrg_solver

end module solver