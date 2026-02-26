module multicharge_solver_cg
    use mctc_env, only: error_type, fatal_error, wp
    use mctc_env_timer, only : timer_type, format_time
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
        ! Maximal number of iterations
        integer, allocatable :: cgmiter
        ! Convergence tolerance
        real(wp), allocatable :: cgtol 
        ! Output verbose
        integer, allocatable :: verbose 
        ! Use iterative CG solver
        logical :: cg = .true.
   end type cg_input

    !> CG solver with Jacobi preconditioner
    type, extends(mchrg_solver_type) :: mchrg_solver_cg
        integer, allocatable :: cgmiter
        real(wp), allocatable :: cgtol
        integer, allocatable :: verbose
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
            else
                self%cgmiter = int(cgmiter_def)
            end if
            
            if (allocated(input%cgtol)) then
                self%cgtol = input%cgtol
            else 
                self%cgtol = cgtol_def
            end if

            if (allocated(input%verbose)) then
                self%verbose = input%verbose
            else 
                self%verbose = 0
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
        ! A matrix of Ax=b system
        real(wp), intent(in)  :: amat(:, :)
        ! Initial search direction (b)
        real(wp), intent(in)  :: xvec(:)
        ! Initial guess and solution
        real(wp), intent(inout) :: vrhs(:)
        ! Inverse matrix and coupled perturbed logical 
        ! not used in CG but required by the interface
        real(wp), intent(out) :: ainv(:, :)
        logical, intent(in), optional :: cpq
        type(error_type), allocatable, intent(out) :: error
        
        ! Maximal number of iterations
        integer :: maxit
        ! Tolerance of the solver
        real(wp) :: tol, tol_square
        
        ! Iterations counter
        integer :: it
        ! Size of the xvec
        integer :: ndim

        ! Search direction (p)
        real(wp), allocatable :: direction(:)
        ! Initial direction norm 
        real(wp) :: bnorm
        ! Residual r = b - A*p
        real(wp), allocatable :: residual(:)
        ! Residual norm 
        real(wp) :: rnorm
        ! Diagonal preconditioner M^-1
        real(wp), allocatable :: Mdiag(:)
        ! Preconditioned residual z=M^-1*r
        real(wp), allocatable ::  zres(:)

        ! Matrix-vector product A*p
        real(wp), allocatable :: Ap(:)
        ! Denominator of the step p^T*Ap
        real(wp) :: denom
        ! Step length
        real(wp) :: alpha

        ! Direction update factor p_new/p
        real(wp) :: beta
        ! Dynamical residuals
        real(wp) :: rz_old, rz_new
        ! Relative residuals test rnorm/bnorm
        real(wp) :: rel_res

        type(cache_container), allocatable :: cache
        type(timer_type) :: timer
        
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

        ! Prepare/cache
        allocate(cache)
        call self%update(cache, vrhs, ainv, cpq)
     
        ! Global thresholds
        tol = self%cgtol
        tol_square = tol**2
        maxit = self%cgmiter
    
        allocate(residual(ndim), direction(ndim), zres(ndim), Ap(ndim), Mdiag(ndim))

        call timer%push("total")
        call timer%push("initialization")   
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

        call timer%pop
        if (self%verbose > 1) then
            write(*, '(a, 1x, a)') "Initialisation time:", format_time(timer%get("initialization"))
        end if

        if (self%verbose > 1) then
            write(*,*)
            write(*,*) ' iter      |residual|        step      relative residual', &
                      &'    Time / s'
        end if

        if (self%verbose == 1) then
            write(*,*)
            write(*,*) ' iter      |residual|        step      relative residual'
        end if
    
        ! Conjugate Gradient iterations
        !$omp shared(Mdiag, amat, ndim, tol_square, vrhs, maxit) private(it) 
        !$omp do schedule(runtime)
        do it = 1, maxit
            call timer%push("iteration")
            !$omp critical (solve_cg_)
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
            residual = residual - alpha * Ap
            
            ! Check convergence
            rnorm = dot_product(residual,residual)
            rel_res = rnorm / bnorm

            if (rel_res <= tol_square) then
                if (self%verbose > 0) then
                    write(*,*)
                    write(*,'(a, i0, a, es15.5)') "CG converged in ", it, " iterations with residual norm ", sqrt(rnorm)
                end if
                exit
            end if
            
            ! Apply preconditioner z = M * r
            zres = residual * Mdiag
            rz_new = dot_product(residual, zres)
            beta = rz_new / rz_old
            direction = zres + beta * direction
            rz_old = rz_new

            call timer%pop

            if (self%verbose == 1) then
                write(*, '(i6,*(1x, es15.5))') it, sqrt(rnorm), alpha, sqrt(rel_res)
            end if

            if (self%verbose > 1) then
                write(*, '(i6,*(1x, es15.5))') it, sqrt(rnorm), alpha, sqrt(rel_res), timer%get("iteration")
            end if
            
            if (it == maxit) then
                call fatal_error(error, "solve_cg: CG did not converge within max iterations.")
            end if
        end do
        !$omp end do
        !$omp end parallel
        call timer%pop
        if (self%verbose > 1) then
            write(*, '(a, 1x, a)') "CG total time : ", format_time(timer%get("total"))
            write(*,*)
        end if 
    
    end subroutine solve

end module multicharge_solver_cg