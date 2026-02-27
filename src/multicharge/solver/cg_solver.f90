module multicharge_solver_cg
    use mctc_env, only: error_type, fatal_error, wp, timer_type, format_time
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

    !> Update method for CG solver (not used, but required by interface)
    subroutine update(self, cache, vrhs, ainv, cpq)
        class(mchrg_solver_cg), intent(in) :: self
        type(cache_container), intent(inout) :: cache
        real(wp), intent(inout) :: vrhs(:)
        real(wp), intent(out) :: ainv(:, :)
        logical, intent(in), optional :: cpq
        ! This routine intentionally left empty.
        ! The dummy arguments are unused, but we must set ainv to avoid -Wunused-dummy-argument.
        ainv = 0.0_wp
    end subroutine update

    !> Solve method for CG solver
    subroutine solve(self, amat, xvec, vrhs, ainv, cpq, error)
        class(mchrg_solver_cg), intent(in) :: self
        ! A matrix of Ax=b system
        real(wp), intent(in)  :: amat(:, :)
        ! Right-hand side vector (b)
        real(wp), intent(in)  :: xvec(:)
        ! On input: initial guess; on output: solution
        real(wp), intent(inout) :: vrhs(:)
        ! Inverse matrix – not computed by CG, but required by interface
        real(wp), intent(out) :: ainv(:, :)
        ! Flag for coupled-perturbed equations (should always be .false. for CG)
        logical, intent(in), optional :: cpq
        ! Error handling
        type(error_type), allocatable, intent(out) :: error
        
        ! Maximal number of iterations
        integer :: maxit
        ! Tolerance of the solver
        real(wp) :: tol, tol_square
        
        ! Iterations counter
        integer :: it
        ! Size of the system
        integer :: ndim

        ! Search direction (p)
        real(wp), allocatable :: direction(:)
        ! Norm of the RHS
        real(wp) :: bnorm
        ! Residual r = b - A*x
        real(wp), allocatable :: residual(:)
        ! Residual norm 
        real(wp) :: rnorm
        ! Diagonal preconditioner M^-1
        real(wp), allocatable :: Mdiag(:)
        ! Preconditioned residual z = M^-1 * r
        real(wp), allocatable ::  zres(:)

        ! Matrix-vector product A*p
        real(wp), allocatable :: Ap(:)
        ! Denominator p^T * A * p
        real(wp) :: denom
        ! Step length
        real(wp) :: alpha

        ! Update factor for search direction
        real(wp) :: beta
        ! Old and new values of r^T * z
        real(wp) :: rz_old, rz_new
        ! Relative residual norm (|r| / |b|)
        real(wp) :: rel_res

        type(cache_container), allocatable :: cache
        type(timer_type) :: timer

        ! --------------------------------------------------------------------
        ! 1.  Basic checks and early returns
        ! --------------------------------------------------------------------

        ! Ensure ainv is always defined (required by intent(out))
        ainv = amat

        ! CG cannot compute the inverse matrix
        if (present(cpq) .and. cpq) then
            call fatal_error(error, "solve_cg: The inverse matrix cannot be calculated using an iterative solver.")
            return
        end if 

        ! Dimensions must match
        ndim = size(xvec)
        if (size(amat,1) /= ndim .or. size(amat,2) /= ndim .or. size(vrhs) /= ndim &
                .or. size(ainv,1) /= ndim .or. size(ainv,2) /= ndim) then
            call fatal_error(error, "solve_cg: dimension mismatch.")
            return
        end if

        ! Prepare cache (unused, but required by interface)
        allocate(cache)
        call self%update(cache, vrhs, ainv, cpq)
     
        ! Global thresholds
        tol = self%cgtol
        tol_square = tol**2
        maxit = self%cgmiter
    
        allocate(residual(ndim), direction(ndim), zres(ndim), Ap(ndim), Mdiag(ndim))

        ! --------------------------------------------------------------------
        ! 2.  Timer start (only if verbose > 0)
        ! --------------------------------------------------------------------
        if (self%verbose > 0) call timer%push("total")
        if (self%verbose > 0) call timer%push("initialization")

        ! --------------------------------------------------------------------
        ! 3.  Jacobi preconditioner (inverse of diagonal)
        ! --------------------------------------------------------------------
        !$omp parallel do default(none) shared(Mdiag, amat, ndim, tol_square) private(it)
        do it = 1, ndim
            Mdiag(it) = amat(it,it)
            if (abs(Mdiag(it)) < tol_square) Mdiag(it) = tol_square
            Mdiag(it) = 1.0_wp / Mdiag(it)
        end do
        !$omp end parallel do
    
        ! --------------------------------------------------------------------
        ! 4.  Initial residual and search direction
        ! --------------------------------------------------------------------
        call symv(amat, vrhs, Ap, alpha=1.0_wp, beta=0.0_wp)   ! Ap = A * x0
        residual = xvec - Ap                                    ! r0 = b - A*x0
        
        zres = residual * Mdiag                                 ! z0 = M^{-1} * r0
        direction = zres                                        ! p0 = z0
        
        bnorm = dot_product(xvec, xvec)
        if (bnorm < tol_square) bnorm = 1.0_wp
        rnorm = dot_product(residual, residual)                 ! |r0|^2
        rz_old = dot_product(residual, zres)                    ! r0^T * z0

        if (self%verbose > 0) call timer%pop   ! pop "initialization"

        ! --------------------------------------------------------------------
        ! 5.  Optional output header
        ! --------------------------------------------------------------------
        if (self%verbose > 1) then
            write(*, '(a, 1x, a)') "Initialisation time:", format_time(timer%get("initialization"))
            write(*,*)
            write(*,*) ' iter      |residual|        step      relative residual    Time / s'
        end if
        if (self%verbose == 1) then
            write(*,*)
            write(*,*) ' iter      |residual|        step      relative residual'
        end if

        ! --------------------------------------------------------------------
        ! 6.  Main CG loop
        ! --------------------------------------------------------------------
        do it = 1, maxit
            if (self%verbose > 0) call timer%push("iteration")

            call symv(amat, direction, Ap, alpha=1.0_wp, beta=0.0_wp)   ! Ap = A * p
            
            denom = dot_product(direction, Ap) + tiny(1.0_wp)
            if (abs(denom) < tol_square) then
                if (self%verbose > 0) call timer%pop   ! pop "iteration"
                exit
            end if
            
            alpha = rz_old / denom                                      ! α = r·z / (p·Ap)
            
            vrhs = vrhs + alpha * direction                            ! x = x + α p
            residual = residual - alpha * Ap                           ! r = r - α Ap
            
            rnorm = dot_product(residual, residual)                    ! |r|^2
            rel_res = rnorm / bnorm                                    ! |r| / |b|

            if (rel_res <= tol_square) then
                if (self%verbose > 0) then
                    write(*,*)
                    write(*,'(a, i0, a, es15.5)') "CG converged in ", it, &
                        & " iterations with residual norm ", sqrt(rnorm)
                end if
                if (self%verbose > 0) call timer%pop   ! pop "iteration"
                exit
            end if
            
            zres = residual * Mdiag                                    ! z = M^{-1} * r
            rz_new = dot_product(residual, zres)                       ! new r·z
            beta = rz_new / rz_old                                     ! β = (r·z)_new / (r·z)_old
            direction = zres + beta * direction                        ! p = z + β p
            rz_old = rz_new

            if (self%verbose > 0) call timer%pop   ! pop "iteration"

            ! Verbose output during iterations
            if (self%verbose == 1) then
                write(*, '(i6,*(1x, es15.5))') it, sqrt(rnorm), alpha, sqrt(rel_res)
            end if
            if (self%verbose > 1) then
                write(*, '(i6,*(1x, es15.5))') it, sqrt(rnorm), alpha, sqrt(rel_res), timer%get("iteration")
            end if
            
            ! Non‑convergence after maxit
            if (it == maxit) then
                if (self%verbose > 0) then
                    call timer%pop   ! pop "iteration"
                    call timer%pop   ! pop "total"
                end if
                call fatal_error(error, "solve_cg: CG did not converge within max iterations.")
                return
            end if
        end do

        ! --------------------------------------------------------------------
        ! 7.  Finalise timer and return
        ! --------------------------------------------------------------------
        if (self%verbose > 0) call timer%pop   ! pop "total"

        if (self%verbose > 1) then
            write(*, '(a, 1x, a)') "CG total time : ", format_time(timer%get("total"))
            write(*,*)
        end if 
    
    end subroutine solve

end module multicharge_solver_cg