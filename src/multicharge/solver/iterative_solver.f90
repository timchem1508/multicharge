module iterative_solver
    use mctc_env, only: error_type, fatal_error, wp
    use multicharge_blas, only: symv, gemv
    use solver_type_cache, only: cache_container
    use print_matrix, only: write_vector
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
    subroutine update_cg(self, cache, amat, xvec, vrhs, ainv, cpq, info)
        class(cg_solver_type), intent(in) :: self
        type(cache_container), intent(inout) :: cache
        real(wp), intent(in)  :: amat(:, :)
        real(wp), intent(in)  :: xvec(:)
        real(wp), allocatable, intent(out) :: vrhs(:)
        real(wp), allocatable, intent(out), optional :: ainv(:, :)
        logical, intent(in), optional :: cpq
        integer, intent(out), optional :: info
        
        integer :: ndim
        
        ndim = size(xvec)
        
        ! Allocate and initialize vrhs (initial guess for CG is zero)
        allocate(vrhs(ndim))
        vrhs = 0.0_wp 
        
        ! Allocate ainv if present (though CG doesn't use it)
        if (present(ainv)) then
            allocate(ainv(ndim, ndim))
        end if
        
        if (present(info)) info = 0
    end subroutine update_cg

    !> Solve method for CG solver
    subroutine solve_cg(self, amat, xvec, vrhs, ainv, cpq, error, info)
        class(cg_solver_type), intent(in) :: self
        type(error_type), allocatable, intent(out) :: error
        real(wp), intent(in)  :: amat(:, :)
        real(wp), intent(in)  :: xvec(:)
        real(wp), allocatable, intent(out) :: vrhs(:)
        real(wp), intent(out), allocatable, optional :: ainv(:, :)
        logical, intent(in), optional :: cpq
        integer, intent(out), optional :: info
        
        integer :: ndim, it, maxit
        real(wp) :: tol, tol_square, bnorm, rnorm, alpha, beta, denom
        real(wp), allocatable :: r(:), p(:), z(:), Ap(:), Mdiag(:)
        integer  :: local_info
        real(wp) :: rz_old, rz_new
        type(cache_container), allocatable :: cache
        
        ! Dimensions match check
        ndim = size(xvec)
        if (size(amat,1) /= ndim .or. size(amat,2) /= ndim) then
            call fatal_error(error, "solve_cg: dimension mismatch.")
            if (present(info)) info = -1 
            return
        end if
    
        ! Prepare/cache
        allocate(cache)
        call self%update(cache, amat, xvec, vrhs, ainv, cpq, info)
     
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
    
        ! Initial residual r = b - A*x (x=vrhs, which is 0.0_wp)
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
        
        if (rnorm / bnorm <= tol_square) then
            if (present(info)) info = 0 
            return
        end if
    
        ! Dynamical residual
        rz_old = dot_product(r,z)
        local_info = -1 
    
        ! Conjugate Gradient iterations
        do it = 1, maxit
            call symv(amat, p, Ap, alpha=1.0_wp, beta=0.0_wp)
            
            ! Compute step size alpha
            denom = dot_product(p,Ap) + tiny(1.0_wp)
            
            if (abs(denom) < tol_square) then
                local_info = 0 
                exit
            end if
            
            alpha = rz_old / denom
            
            ! Update solution and residual
            vrhs = vrhs + alpha * p
            r = r - alpha * Ap
            
            ! Check convergence
            rnorm = dot_product(r,r)
            if (rnorm / bnorm <= tol_square) then
                write(*,*) "CG converged in ", it, " iterations."
                call write_vector(vrhs, "CG Solution Vector")
                local_info = 0 
                exit
            end if
            
            ! Apply preconditioner z = M * r
            z = r * Mdiag
            rz_new = dot_product(r, z)
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

end module iterative_solver