module multicharge_solver_cg
    use mctc_env, only: error_type, fatal_error, wp
    use multicharge_blas, only: symv, gemv
    use multicharge_solver_type, only: mchrg_solver_type, mchrg_solver_input
    use multicharge_solver_cache, only: cache_container, mchrg_solver_cache
    use print_matrix, only: write_vector, write_matrix
    implicit none
    private

    public :: mchrg_solver_cg, new_cg_solver, cg_input

    type, extends(mchrg_solver_cache), public :: cg_cache
    end type cg_cache

    !> Input for CG solver
    type, extends(mchrg_solver_input) :: cg_input
        !> Maximal number of iterations
        integer, allocatable :: cgmiter
        !> Convergence tolerance
        real(wp), allocatable :: cgtol 
        !> Use iterative CG solver
        logical :: cg = .true.
   end type cg_input

    !> CG solver with Jacobi preconditioner
    type, extends(mchrg_solver_type) :: mchrg_solver_cg
        integer, allocatable :: cgmiter
        real(wp), allocatable :: cgtol
    contains
        procedure :: solve
        procedure :: update
    end type mchrg_solver_cg

contains

    subroutine new_cg_solver(self, input)
        class(mchrg_solver_type), intent(out) :: self
        type(cg_input), intent(in) :: input     
                
        real(wp), parameter :: cgmiter_def = 1000
        real(wp), parameter :: cgtol_def = 1.0e-15_wp

        select type (self)
        type is (mchrg_solver_cg)
            
            self%need_pos_def = .true.

            if (allocated(input%cgmiter)) then
                self%cgmiter = input%cgmiter
                write(*,*) "Maximum iterations:", self%cgmiter 
            else
                write(*,*) "Default maximum number of iterations is used: 1000 it."
                self%cgmiter = int(cgmiter_def)
            end if
            
            if (allocated(input%cgtol)) then
                self%cgtol = input%cgtol
                write(*,*) "Tolerance:", self%cgtol
            else 
                write(*,*) "Default tolerance is used: 1.0e-15"
                self%cgtol = cgtol_def
            end if
            
        end select

    end subroutine new_cg_solver

    !> Update method for CG solver
    subroutine update(self, cache, vrhs, ainv, cpq)
        class(mchrg_solver_cg), intent(in) :: self
        type(cache_container), intent(inout) :: cache
        real(wp), intent(inout) :: vrhs(:)
        real(wp), intent(out) :: ainv(:, :)
        logical, intent(in), optional :: cpq

    end subroutine update

    !> Solve method for CG solver
    subroutine solve(self, amat, xvec, vrhs, ainv, cpq, error)
        class(mchrg_solver_cg), intent(in) :: self
        !> A matrix of Ax=b system
        real(wp), intent(in)  :: amat(:, :)
        !> Initial search direction (b)
        real(wp), intent(in)  :: xvec(:)
        !> Initial guess and solution
        real(wp), intent(inout) :: vrhs(:)
        real(wp), intent(out) :: ainv(:, :)
        logical, intent(in), optional :: cpq
        type(error_type), allocatable, intent(out) :: error
        
        !> Maximal number of iterations
        integer :: maxit
        !> Tolerance of the solver
        real(wp) :: tol, tol_square
        
        !> Iterations counter
        integer :: it
        !> Size of the xvec
        integer :: ndim

        !> Search direction (p)
        real(wp), allocatable :: direction(:)
        !> Initial direction norm 
        real(wp) :: bnorm
        !> Residual r = b - A*p
        real(wp), allocatable :: residual(:)
        !> Residual norm 
        real(wp) :: rnorm
        !> Diagonal preconditioner M^-1
        real(wp), allocatable :: Mdiag(:)
        !> Preconditioned residual z=M^-1*r
        real(wp), allocatable ::  zres(:)

        !> Matrix-vector product A*p
        real(wp), allocatable :: Ap(:)
        !> Denominator of the step p^T*Ap
        real(wp) :: denom
        !> Step length
        real(wp) :: alpha

        !> Direction update factor p_new/p
        real(wp) :: beta
        !> Dynamical residuals
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
        allocate(cache)
        call self%update(cache, vrhs, ainv, cpq)
     
        ! Global thresholds
        tol = self%cgtol
        tol_square = tol**2
        maxit = self%cgmiter

        !write(*,*) "CG Solver: max iterations = ", maxit
    
        allocate(residual(ndim), direction(ndim), zres(ndim), Ap(ndim), Mdiag(ndim))
    
        ! Jacobi preconditioner (inverse of diagonal)
        !$omp parallel default(none) &
        !$omp shared(Mdiag, amat, ndim, tol_square) private(it) 
        !$omp do schedule(runtime)
        do it = 1, ndim
            Mdiag(it) = amat(it,it)
            if (abs(Mdiag(it)) < tol_square) Mdiag(it) = tol_square
            Mdiag(it) = 1.0_wp / Mdiag(it)
        end do
        !$omp end do
        !$omp critical (solve_cg_)
    
        ! Initial residual r = b - A*x
        call symv(amat, vrhs, Ap, alpha=1.0_wp, beta=0.0_wp)
        residual = xvec - Ap
        
        ! Apply preconditioner z = M * r
        zres = residual * Mdiag
        
        ! Initial search direction
        direction = zres
        
        ! Initial direction update factor
        bnorm = dot_product(xvec,xvec)
        if (bnorm < tol_square) bnorm = 1.0_wp
        rnorm = dot_product(residual,residual)
        
        
        ! Dynamical residual
        rz_old = dot_product(residual,zres)
        !$omp end critical (solve_cg_)
    
        ! Conjugate Gradient iterations
        !$omp shared(Mdiag, amat, ndim, tol_square, vrhs, maxit) private(it) 
        !$omp do schedule(runtime)
        do it = 1, maxit
            call symv(amat, direction, Ap, alpha=1.0_wp, beta=0.0_wp)
            
            ! Compute step size alpha
            denom = dot_product(direction,Ap) + tiny(1.0_wp)
            
            if (abs(denom) < tol_square) then
                exit
            end if
            
            ! Step update
            alpha = rz_old / denom
            
            ! Update solution and residual
            vrhs = vrhs + alpha * direction
            !call write_vector(vrhs, "CG Solution Vector")
            !write(*,*) " iteration ", it
            residual = residual - alpha * Ap
            
            ! Check convergence
            rnorm = dot_product(residual,residual)
            if (rnorm / bnorm <= tol_square) then
                write(*,*) "CG converged in ", it, " iterations."
                !call write_vector(vrhs, "CG Solution Vector")
                exit
            end if
            
            ! Apply preconditioner z = M * r
            zres = residual * Mdiag
            rz_new = dot_product(residual, zres)
            beta = rz_new / rz_old
            direction = zres + beta * direction
            rz_old = rz_new
            
            if (it == maxit) then
                call fatal_error(error, "solve_cg: CG did not converge within max iterations.")
            end if
        end do
        !$omp end do
        !$omp end parallel
    
    end subroutine solve

end module multicharge_solver_cg