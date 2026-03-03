! This file is part of multicharge.
! SPDX-Identifier: Apache-2.0
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

!> @file multicharge/solver/cg.f90
!> Provides implementation of the conjugate gradient solver for linear systems of equations.

module cg_solver
    use iso_fortran_env, only : output_unit
    use mctc_env, only: error_type, fatal_error, wp, timer_type, format_time
    use multicharge_blas, only: symv, gemv
    use solver_type, only: mchrg_solver_type, mchrg_solver_input
    use solver_cache, only: cache_container, mchrg_solver_cache
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

        ! Default values        
        integer, parameter :: cgmiter_def = 1000
        real(wp), parameter :: cgtol_def = 1.0e-15_wp
        integer, parameter :: verbose_def = 0

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
                self%verbose = verbose_def
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

        ainv = 0.0_wp
    end subroutine update

    !> Solve method for CG solver
    subroutine solve(self, amat, xvec, vrhs, ainv, cpq, new_unit, error)
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
        ! Output unit
        integer, intent(in), optional :: new_unit
        !> Error handling
        type(error_type), allocatable, intent(out) :: error
        
        ! Maximal number of iterations
        integer :: maxit
        ! Tolerance of the solver
        real(wp) :: tol, tol_square
        
        ! Counters
        integer :: it, iat
        ! Size of the system
        integer :: ndim

        ! Search direction (p)
        real(wp), allocatable :: direction(:)
        ! Norm of the RHS
        real(wp) :: bnorm
        ! Residual 
        real(wp), allocatable :: residual(:)
        ! Residual norm 
        real(wp) :: rnorm
        ! Diagonal preconditioner 
        real(wp), allocatable :: Mdiag(:)
        ! Preconditioned residual 
        real(wp), allocatable ::  zres(:)

        ! Matrix-vector product 
        real(wp), allocatable :: Ap(:)
        ! Denominator of the step length
        real(wp) :: denom
        ! Step length
        real(wp) :: alpha

        ! Update factor for search direction
        real(wp) :: beta
        real(wp) :: rz_old, rz_new
        ! Relative residual norm (|r| / |b|)
        real(wp) :: rel_res

        type(cache_container), allocatable :: cache
        type(timer_type) :: timer
        integer :: unit

        if (present(new_unit)) then
            unit = new_unit
        else
            unit = output_unit
        end if

        ainv = amat

        ! CG cannot compute the inverse matrix
        if (present(cpq) .and. cpq) then
            call fatal_error(error, "solve_cg: The inverse matrix cannot be calculated using an iterative solver.")
            return
        end if 

        ! Dimensions check
        ndim = size(xvec)
        if (size(amat,1) /= ndim .or. size(amat,2) /= ndim .or. size(vrhs) /= ndim &
                .or. size(ainv,1) /= ndim .or. size(ainv,2) /= ndim) then
            call fatal_error(error, "solve_cg: dimension mismatch.")
            return
        end if

        allocate(cache)
        call self%update(cache, vrhs, ainv, cpq)

        tol = self%cgtol
        tol_square = tol**2
        maxit = self%cgmiter
    
        allocate(residual(ndim), direction(ndim), zres(ndim), Ap(ndim), Mdiag(ndim))

        if (self%verbose > 1) call timer%push("total")
        if (self%verbose > 1) call timer%push("initialization")

        ! Diagonal preconditioner M^-1 (Jacobi preconditioner)
        !$omp parallel do default(none) shared(Mdiag, amat, ndim, tol_square) private(iat)
        do iat = 1, ndim
            Mdiag(iat) = amat(iat,iat)
            if (abs(Mdiag(iat)) < tol_square) Mdiag(iat) = tol_square
            Mdiag(iat) = 1.0_wp / Mdiag(iat)
        end do
        !$omp end parallel do
    
        ! Initial residual r = b - A*x
        call symv(amat, vrhs, Ap, alpha=1.0_wp, beta=0.0_wp)   
        residual = xvec - Ap                                    
        
        ! Initial preconditioned residual z = M^-1 * r
        zres = residual * Mdiag                                 
        direction = zres                                    
        
        ! Initial norm of the right-hand side (b)
        bnorm = dot_product(xvec, xvec)
        if (bnorm < tol_square) bnorm = 1.0_wp

        ! Initial residual norm and r^T * z
        rnorm = dot_product(residual, residual)                
        rz_old = dot_product(residual, zres)                    

        if (self%verbose > 1) call timer%pop  ! initialization timer stop

        ! Print header
        call print_cg_header(unit, self%verbose, timer)

        ! Main CG iteration loop
        do it = 1, maxit

            if (self%verbose > 1) call timer%push("iteration")

            ! Matrix-vector product Ap = A * p
            call symv(amat, direction, Ap, alpha=1.0_wp, beta=0.0_wp)

            denom = 0.0_wp
            !$omp parallel do reduction(+:denom) default(none) &
            !$omp shared(ndim,direction,Ap) private(iat)
            do iat = 1, ndim
                denom = denom + direction(iat) * Ap(iat)
            end do
            !$omp end parallel do

            denom = denom + tiny(1.0_wp)

            if (abs(denom) < tol_square) then
                if (self%verbose > 0) call timer%pop
                exit
            end if

            ! Step length alpha = (r^T * z) / (p^T * A * p)
            alpha = rz_old / denom

            ! Update solution x = x + alpha * p and residual r = r - alpha * A*p
            !$omp parallel do default(none) shared(ndim,vrhs,residual,alpha,direction,Ap) private(iat)
            do iat = 1, ndim
                vrhs(iat)     = vrhs(iat)     + alpha * direction(iat)
                residual(iat) = residual(iat) - alpha * Ap(iat)
            end do
            !$omp end parallel do

            ! Compute the new residual norm
            rnorm = 0.0_wp
            !$omp parallel do reduction(+:rnorm) default(none) &
            !$omp shared(ndim,residual) private(iat)
            do iat = 1, ndim
                rnorm = rnorm + residual(iat) * residual(iat)
            end do
            !$omp end parallel do

            ! Relative residual norm to check convergence
            rel_res = rnorm / bnorm

            if (rel_res <= tol_square) then
                if (self%verbose > 0) then
                    call print_cg_convergence(unit, it, sqrt(rnorm), self%verbose)
                    call timer%pop
                end if
                exit
            end if

            ! Preconditioned updated residual z = M^-1 * r
            !$omp parallel do default(none) shared(ndim,zres,residual,Mdiag) private(iat)
            do iat = 1, ndim
                zres(iat) = residual(iat) * Mdiag(iat)
            end do
            !$omp end parallel do

            rz_new = 0.0_wp
            !$omp parallel do reduction(+:rz_new) default(none) &
            !$omp shared(ndim,residual,zres) private(iat)
            do iat = 1, ndim
                rz_new = rz_new + residual(iat) * zres(iat)
            end do
            !$omp end parallel do

            ! Update search direction p = z + beta * p
            beta   = rz_new / rz_old
            rz_old = rz_new

            !$omp parallel do default(none) shared(ndim,direction,zres,beta) private(iat)
            do iat = 1, ndim
                direction(iat) = zres(iat) + beta * direction(iat)
            end do
            !$omp end parallel do

            if (self%verbose > 1) call timer%pop ! iteration timer stop

            ! Print iteration progress
            call print_cg_iteration(unit, it, sqrt(rnorm), alpha, sqrt(rel_res), self%verbose, timer)

            if (it == maxit) then
                if (self%verbose > 1) then
                    call timer%pop   ! pop "iteration"
                    call timer%pop   ! pop "total"
                end if
                call fatal_error(error, "solve_cg: CG did not converge within max iterations.")
                return
            end if

        end do

        if (self%verbose > 1) call timer%pop   ! pop total

        ! Print final summary
        call print_cg_final(unit, timer, self%verbose)
    
    end subroutine solve

    !> Print header for CG solver
    subroutine print_cg_header(unit, verbose, timer)
        integer, intent(in) :: unit, verbose
        type(timer_type), intent(in), optional :: timer

        if (verbose > 1) then
            write(unit, '(a, 1x, a)') "Initialisation time:", format_time(timer%get("initialization"))
            write(unit, '(a)') ''
            write(unit, '(2X,A,6X,A,8X,A,6X,A,4X,A)') &
                'iter', '|residual|', 'step', 'relative residual', 'Time / s'
        else if (verbose == 1) then
            write(unit, '(a)') ''
            write(unit, '(2X,A,6X,A,8X,A,6X,A)') &
                'iter', '|residual|', 'step', 'relative residual'
        end if
    end subroutine print_cg_header

    !> Print convergence message
    subroutine print_cg_convergence(unit, iter, res_norm, verbose)
        integer, intent(in) :: unit, iter, verbose
        real(wp), intent(in) :: res_norm

        if (verbose > 0) then
            write(unit, '(a)') ''
            write(unit, '(a, i0, a, es15.5)') &
                "CG converged in ", iter, " iterations with residual norm ", res_norm
        end if
    end subroutine print_cg_convergence

    !> Print iteration progress
    subroutine print_cg_iteration(unit, iter, res_norm, alpha, rel_res, verbose, timer)
        integer, intent(in) :: unit, iter, verbose
        real(wp), intent(in) :: res_norm, alpha, rel_res
        type(timer_type), intent(in), optional :: timer

        if (verbose == 1) then
            write(unit, '(i6,*(1x, es15.5))') iter, res_norm, alpha, rel_res
        else if (verbose > 1) then
            write(unit, '(i6,*(1x, es15.5))') iter, res_norm, alpha, rel_res, timer%get("iteration")
        end if
    end subroutine print_cg_iteration

    !> Print final summary
    subroutine print_cg_final(unit, timer, verbose)
        integer, intent(in) :: unit, verbose
        type(timer_type), intent(in) :: timer

        if (verbose > 1) then
            write(unit, '(a, 1x, a)') "CG total time : ", format_time(timer%get("total"))
            write(unit, '(a)') ''
        end if
    end subroutine print_cg_final

end module cg_solver