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

module multicharge_solver_cg
    use iso_fortran_env, only : output_unit
    use mctc_env, only: error_type, fatal_error, wp, timer_type, format_time
    use multicharge_blas, only: axpy, scal, dot, symv, gemv
    use multicharge_lapack, only: sytrf, sytrs
    use multicharge_solver_type, only: mchrg_solver_type, mchrg_solver_input
    implicit none
    private

    public :: cg_solver, new_cg_solver, cg_input

    !> Input for CG solver
    type, extends(mchrg_solver_input) :: cg_input
        !> Maximal number of iterations
        integer, allocatable :: cgmiter 
        !> Convergence tolerance
        real(wp), allocatable :: cgtol 
        !> Output verbosity
        integer, allocatable :: verbosity
        !> Use iterative CG solver
        logical :: cg = .true.
   end type cg_input

    !> CG solver with Jacobi preconditioner
    type, extends(mchrg_solver_type) :: cg_solver
        integer, allocatable :: cgmiter
        real(wp), allocatable :: cgtol
        integer, allocatable :: verbosity
    contains
        procedure :: solve
    end type cg_solver

    real(wp), parameter :: eps = tiny(1.0_wp)

    ! Default values        
    integer, parameter :: cgmiter_def = 1000
    real(wp), parameter :: cgtol_def = 1.0e-15_wp
    integer, parameter :: verbosity_def = 0

contains

    subroutine new_cg_solver(self, input)
        class(cg_solver), intent(out) :: self
        type(cg_input), intent(in) :: input     

        self%need_pos_def = .true.
        if (allocated(input%cgmiter)) then
            self%cgmiter = input%cgmiter
        else
            self%cgmiter = cgmiter_def
        end if
        
        if (allocated(input%cgtol)) then
            self%cgtol = input%cgtol
        else 
            self%cgtol = cgtol_def
        end if
        if (allocated(input%verbosity)) then
            self%verbosity = input%verbosity
        else 
            self%verbosity = verbosity_def
        end if


    end subroutine new_cg_solver

    !> Solve procedure for the classical cg and block-cg

    subroutine solve(self, amat, xvec, vrhs, ainv, cpq, new_unit, error)
        class(cg_solver), intent(in) :: self
        !> A matrix of Ax=b system
        real(wp), intent(in)  :: amat(:, :)
        !> Right-hand side vector (b)
        real(wp), intent(in)  :: xvec(:)
        !> On input: initial guess; on output: solution
        real(wp), intent(inout), contiguous :: vrhs(:)
        !> Inverse matrix – not computed by CG, but required by interface
        real(wp), intent(out), optional :: ainv(:, :)
        !> Flag for coupled-perturbed equations
        logical, intent(in), optional :: cpq
        !> Output unit
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
        ! Search direction
        real(wp), allocatable :: dir(:)
        ! Norm of the RHS
        real(wp) :: xvecnorm
        ! Residual 
        real(wp), allocatable :: res(:)
        ! Residual norm 
        real(wp) :: resnorm
        ! Diagonal preconditioner 
        real(wp), allocatable :: prec(:)
        ! Preconditioned residual 
        real(wp), allocatable ::  precres(:)
        ! amat-dir product 
        real(wp), allocatable :: Adir(:)
        ! Denominator of the step length
        real(wp) :: denom
        ! Step length
        real(wp) :: step
        ! Update factor for search direction
        real(wp) :: updfact
        ! Profection of preconditioned residual and an original one
        real(wp) :: resdot_old, resdot_new
        ! Relative residual norm (|resnorm| / |vrhs|)
        real(wp) :: rel_resnorm

        type(timer_type) :: timer
        integer :: unit

        ! CG cannot compute the inverse matrix
        if (present(ainv) .or. present(cpq)) then
            call fatal_error(error, "The inverse matrix cannot be calculated using an iterative solver.")
            return
        end if 

        if (present(new_unit)) then
            unit = new_unit
        else
            unit = output_unit
        end if

        ! Dimensions check
        ndim = size(xvec)
        if (size(amat,1) /= ndim .or. size(amat,2) /= ndim .or. size(vrhs) /= ndim) then
            call fatal_error(error, "dimension mismatch.")
            return
        end if

        tol = self%cgtol
        tol_square = tol**2
        maxit = self%cgmiter
    
        allocate(res(ndim), dir(ndim), precres(ndim), Adir(ndim), prec(ndim))

        if (self%verbosity > 1) call timer%push("total")
        if (self%verbosity > 1) call timer%push("initialization")

        ! Diagonal preconditioner 
        do iat = 1, ndim
            prec(iat) = 1.0_wp / (amat(iat,iat) + eps)
        end do
    
        ! Initial residual
        call symv(amat, vrhs, Adir, alpha=1.0_wp, beta=0.0_wp)   
        res(:) = xvec(:) - Adir(:)                                                                
        
        ! Initial preconditioned residual precres = M^-1 * res
        precres(:) = res(:) * prec(:)                                 
        dir(:) = precres(:)                                    
        
        ! Initial norm 
        xvecnorm = dot(xvec, xvec)
        if (xvecnorm < tol_square) xvecnorm = 1.0_wp

        ! Initial resnorm and res^T * precres
        resnorm = dot(res, res)                
        resdot_old = dot(res, precres)                    

        if (self%verbosity > 1) call timer%pop 

        ! Print header
        call print_cg_header(unit, self%verbosity, maxit, tol, timer)

        ! Main CG iteration loop

        do it = 1, maxit

            if (self%verbosity > 1) call timer%push("iteration")

            ! Matrix-vector product
            call symv(amat, dir, Adir, alpha=1.0_wp, beta=0.0_wp)

            denom = dot(dir, Adir)
            if (abs(denom) < tol_square) then
                if (self%verbosity > 0) call timer%pop
                exit
            end if

            ! Step length step = (res^T * precres) / (dir^T * amat * dir)
            step = resdot_old / (denom + eps)

            ! Update solution vrhs = vrhs + step * dir 
            call axpy(xvec=dir, yvec=vrhs, alpha=step)
            !Update residual res = res - step * amat * dir
            call axpy(xvec=Adir, yvec=res, alpha=-step)

            ! Compute the new residual norm
            resnorm = dot(res, res)

            ! Relative residual norm to check convergence
            rel_resnorm = resnorm / xvecnorm

            if (rel_resnorm <= tol_square) then
                if (self%verbosity > 0) then
                    call print_cg_convergence(unit, it, sqrt(resnorm), self%verbosity)
                    call timer%pop
                end if
                exit
            end if

            ! Updated preconditioned residual
            precres(:) = prec(:) * res(:)

            ! Update search direction
            resdot_new = dot(res, precres)           
            updfact = resdot_new / (resdot_old + eps)
            resdot_old = resdot_new
            call scal(alpha=updfact, xvec=dir)
            call axpy(xvec=precres, yvec=dir, alpha=1.0_wp)

            ! iteration timer pop
            if (self%verbosity > 1) call timer%pop 

            ! Print iteration progress
            call print_cg_iteration(unit, it, sqrt(resnorm), step, sqrt(rel_resnorm), self%verbosity, timer)

            if (it == maxit) then
                if (self%verbosity > 1) then
                    ! pop "iteration"
                    call timer%pop
                    ! pop "total"
                    call timer%pop   
                end if
                call fatal_error(error, "CG did not converge within max iterations.")
                return
            end if

        end do

        ! pop total
        if (self%verbosity > 1) call timer%pop   

        ! Print final summary
        call print_cg_final(unit, timer, self%verbosity)
    
    end subroutine solve


    !> Print header for CG solver
    subroutine print_cg_header(unit, verbosity, maxit, tol, timer)
        integer, intent(in) :: unit, verbosity, maxit
        real(wp), intent(in) :: tol
        type(timer_type), intent(in), optional :: timer

        if (verbosity > 1) then
            write(unit, '(a)') "Using Conjugate Gradient Solver"
            write(unit, '(a)')
            write(unit, '(a, 1x, i6)') "Max iterations : ", maxit
            write(unit, '(a, 1x, es10.2)') "Tolerance      : ", tol
            write(unit, '(a)') "Preconditioner : Jacobi (Diagonal)"
            write(unit, '(a, 1x, a)') "Initialisation time:", format_time(timer%get("initialization"))
            write(unit, '(a)') ''
            write(unit, '(2X,A,6X,A,8X,A,6X,A,4X,A)') &
                'iter', '|residual|', 'step', 'relative residual', 'Time / s'
        else if (verbosity == 1) then
            write(unit, '(a)') "Using Conjugate Gradient Solver"
            write(unit, '(a)')
            write(unit, '(a, 1x, i6)') "Max iterations : ", maxit
            write(unit, '(a, 1x, es10.2)') "Tolerance      : ", tol
            write(unit, '(a)') "Preconditioner : Jacobi (Diagonal)"
            write(unit, '(a)') ''
            write(unit, '(2X,A,6X,A,8X,A,6X,A)') &
                'iter', '|residual|', 'step', 'relative residual'
        end if
    end subroutine print_cg_header

    !> Print convergence message
    subroutine print_cg_convergence(unit, iter, res_norm, verbosity)
        integer, intent(in) :: unit, iter, verbosity
        real(wp), intent(in) :: res_norm

        if (verbosity > 0) then
            write(unit, '(a)') ''
            write(unit, '(a, i0, a, es15.5)') &
                "CG converged in ", iter, " iterations with residual norm ", res_norm
            write(unit, '(a)') ''
        end if
    end subroutine print_cg_convergence

    !> Print iteration progress
    subroutine print_cg_iteration(unit, iter, res_norm, step, rel_resnorm, verbosity, timer)
        integer, intent(in) :: unit, iter, verbosity
        real(wp), intent(in) :: res_norm, step, rel_resnorm
        type(timer_type), intent(in), optional :: timer

        if (verbosity == 1) then
            write(unit, '(i6,*(1x, es15.5))') iter, res_norm, step, rel_resnorm
        else if (verbosity > 1) then
            write(unit, '(i6,*(1x, es15.5))') iter, res_norm, step, rel_resnorm, timer%get("iteration")
        end if
    end subroutine print_cg_iteration

    !> Print final summary
    subroutine print_cg_final(unit, timer, verbosity)
        integer, intent(in) :: unit, verbosity
        type(timer_type), intent(in) :: timer

        if (verbosity > 1) then
            write(unit, '(a, 1x, a)') "CG total time : ", format_time(timer%get("total"))
            write(unit, '(a)') ''
        end if
    end subroutine print_cg_final

end module multicharge_solver_cg