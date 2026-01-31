module iterative_solver
    use mctc_env, only: error_type, fatal_error, wp
    use multicharge_blas, only: symv, gemv
    use solver_type_cache, only: cache_container
    use print_matrix, only: write_vector, write_matrix
    use solver, only: mchrg_solver_type
    implicit none
    private

    public :: cg_solver_type

    !> CG solver with Jacobi preconditioner
    type, extends(mchrg_solver_type) :: cg_solver_type
    contains
       procedure :: solve => solve_cg
       procedure :: update => update_cg
    end type cg_solver_type

contains

    !> Update method for CG solver
    subroutine update_cg(self, cache, amat, xvec, vrhs, ainv, cpq)
        class(cg_solver_type), intent(in) :: self
        type(cache_container), intent(inout) :: cache
        real(wp), intent(in)  :: amat(:, :)
        real(wp), intent(in)  :: xvec(:)
        real(wp), intent(inout) :: vrhs(:)
        real(wp), intent(out) :: ainv(:, :)
        logical, intent(in), optional :: cpq

    end subroutine update_cg

    !> Solve method for CG solver
    subroutine solve_cg(self, amat, xvec, vrhs, ainv, cpq, error)
        class(cg_solver_type), intent(in) :: self
        type(error_type), allocatable, intent(out) :: error
        real(wp), intent(in)  :: amat(:, :)
        real(wp), intent(in)  :: xvec(:)
        real(wp), intent(inout) :: vrhs(:)
        real(wp), intent(out) :: ainv(:, :)
        logical, intent(in), optional :: cpq
        
        integer :: ndim, it, maxit
        real(wp) :: tol, tol_square, bnorm, rnorm, alpha, beta, denom
        real(wp), allocatable :: r(:), p(:), z(:), Ap(:), Mdiag(:)
        real(wp) :: rz_old, rz_new
        type(cache_container), allocatable :: cache
        
        if (present(cpq) .and. cpq .eqv. .true.) then
            call fatal_error(error, "solve_cg: The inverse matrix cannot be calculated using an iterative solver.")
            return
        end if 

        ! Dimensions match check
        ndim = size(xvec)
        if (size(amat,1) /= ndim .or. size(amat,2) /= ndim) then
            call fatal_error(error, "solve_cg: dimension mismatch.")
            return
        end if

        ainv = amat
        !call write_vector(vrhs, "Initial VRHS Vector")
        ! Prepare/cache
        !allocate(cache)
        !call self%update(cache, amat, xvec, vrhs, ainv, cpq)
     
        ! Global thresholds
        tol = 1.0e-11_wp
        tol_square = tol**2
        maxit = max(10, ndim*20)
        write(*,*) "CG Solver: max iterations = ", maxit
    
        allocate(r(ndim), p(ndim), z(ndim), Ap(ndim), Mdiag(ndim))
    
        ! Jacobi preconditioner (inverse of diagonal)
        do it = 1, ndim
            Mdiag(it) = amat(it,it)
            if (abs(Mdiag(it)) < tol_square) Mdiag(it) = tol_square
            Mdiag(it) = 1.0_wp / Mdiag(it)
        end do
    
        ! Initial residual r = b - A*x
        call symv(amat, vrhs, Ap, alpha=1.0_wp, beta=0.0_wp)
        r = xvec - Ap
        
        ! Apply preconditioner z = M * r
        z = r * Mdiag
        
        ! Initial search direction
        p = z
        
        ! Initial direction update factor
        bnorm = dot_product(xvec,xvec)
        if (bnorm < tol_square) bnorm = 1.0_wp
        rnorm = dot_product(r,r)
        
        ! Dynamical residual
        rz_old = dot_product(r,z)
    
        ! Conjugate Gradient iterations
        do it = 1, maxit
            call symv(amat, p, Ap, alpha=1.0_wp, beta=0.0_wp)
            
            ! Compute step size alpha
            denom = dot_product(p,Ap) + tiny(1.0_wp)
            
            if (abs(denom) < tol_square) then
                exit
            end if
            
            alpha = rz_old / denom
            
            ! Update solution and residual
            vrhs = vrhs + alpha * p
            !call write_vector(vrhs, "CG Solution Vector")
            !write(*,*) " iteration ", it
            r = r - alpha * Ap
            
            ! Check convergence
            rnorm = dot_product(r,r)
            if (rnorm / bnorm <= tol_square) then
                write(*,*) "CG converged in ", it, " iterations."
                !call write_vector(vrhs, "CG Solution Vector")
                exit
            end if
            
            ! Apply preconditioner z = M * r
            z = r * Mdiag
            rz_new = dot_product(r, z)
            beta = rz_new / rz_old
            p = z + beta * p
            rz_old = rz_new
            
            if (it == maxit) then
                call fatal_error(error, "solve_cg: CG did not converge within max iterations.")
            end if
        end do
    
    end subroutine solve_cg

end module iterative_solver