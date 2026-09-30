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
!> Provides implementation of the conjugate gradient solver for linear systems
!> of equations

!> Conjugate gradient solver
module multicharge_solver_cg
   use iso_fortran_env, only : output_unit
   use mctc_env, only : error_type, fatal_error, format_time, timer_type, wp
   use multicharge_blas, only : axpy, dot, scal, symv
   use multicharge_solver_type, only : mchrg_solver_input, mchrg_solver_type
   implicit none
   private

   public :: cg_solver, new_cg_solver, cg_input


   !> Input configuration for the conjugate-gradient solver
   type, extends(mchrg_solver_input) :: cg_input

      !> Maximum number of iterations
      integer, allocatable :: cgmiter

      !> Convergence tolerance
      real(wp), allocatable :: cgtol

      !> Output verbosity
      integer, allocatable :: verbosity

      !> Whether to use the iterative conjugate-gradient solver
      logical :: cg = .true.

   end type cg_input

   !> Conjugate-gradient solver with a Jacobi preconditioner
   type, extends(mchrg_solver_type) :: cg_solver

      !> Maximum number of iterations
      integer, allocatable :: cgmiter

      !> Convergence tolerance
      real(wp), allocatable :: cgtol

      !> Output verbosity
      integer, allocatable :: verbosity
   contains

      !> Solve the linear system iteratively
      procedure :: solve

   end type cg_solver

   !> Positive number used to prevent division by zero
   real(wp), parameter :: eps = tiny(1.0_wp)

   !> Default maximum number of iterations
   integer, parameter :: cgmiter_def = 1000

   !> Default convergence tolerance
   real(wp), parameter :: cgtol_def = 1.0e-15_wp

   !> Default output verbosity
   integer, parameter :: verbosity_def = 0


contains


!> Construct a conjugate-gradient solver from its input configuration
subroutine new_cg_solver(self, input)

   !> Conjugate-gradient solver instance
   class(cg_solver), intent(out) :: self

   !> Conjugate-gradient solver configuration
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


!> Solve a linear system with a diagonally preconditioned CG method
subroutine solve(self, amat, xvec, vrhs, ainv, cpq, new_unit, error)

   !> Conjugate-gradient solver instance
   class(cg_solver), intent(in) :: self

   !> Dense coefficient matrix of the linear system
   real(wp), intent(in) :: amat(:, :)

   !> Right-hand side vector
   real(wp), intent(in) :: xvec(:)

   !> On input: initial guess; on output: solution
   real(wp), intent(inout), contiguous :: vrhs(:)

   !> Inverse matrix
   real(wp), intent(out), optional :: ainv(:, :)

   !> Whether to solve coupled-perturbed equations
   logical, intent(in), optional :: cpq

   !> Output unit
   integer, intent(in), optional :: new_unit

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: maxit, it, iat, ndim, unit
   real(wp) :: tol, tol_square, xvecnorm, resnorm, rel_resnorm
   real(wp) :: denom, step, updfact, resdot_old, resdot_new
   real(wp), allocatable :: dir(:), res(:), prec(:), precres(:), adir(:)
   type(timer_type) :: timer

   ! CG cannot compute the inverse matrix
   if (present(ainv) .or. present(cpq)) then
      call fatal_error(error, &
         & "The inverse matrix cannot be calculated using an iterative solver.")
      return
   end if

   if (present(new_unit)) then
      unit = new_unit
   else
      unit = output_unit
   end if

   ! Dimensions check
   ndim = size(xvec)
   if (size(vrhs) /= ndim) then
      call fatal_error(error, "Dimension mismatch between xvec and vrhs.")
      return
   end if

   if (size(amat, 1) /= ndim .or. size(amat, 2) /= ndim) then
      call fatal_error(error, "Dimension mismatch of the coefficient matrix.")
      return
   end if

   tol = self%cgtol
   tol_square = tol**2
   maxit = self%cgmiter

   allocate(res(ndim), dir(ndim), precres(ndim), adir(ndim), prec(ndim))

   if (self%verbosity > 1) call timer%push("total")
   if (self%verbosity > 1) call timer%push("initialization")

   ! Diagonal preconditioner
   do iat = 1, ndim
      prec(iat) = 1.0_wp / (amat(iat, iat) + eps)
   end do

   ! Initial residual
   call symv(amat, vrhs, adir, alpha=1.0_wp, beta=0.0_wp)
   res(:) = xvec(:) - adir(:)

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
      call symv(amat, dir, adir, alpha=1.0_wp, beta=0.0_wp)

      denom = dot(dir, adir)
      if (abs(denom) < tol_square) then
         if (self%verbosity > 0) call timer%pop
         exit
      end if

      ! Step length step = (res^T * precres) / (dir^T * amat * dir)
      step = resdot_old / (denom + eps)

      ! Update solution vrhs = vrhs + step * dir
      call axpy(xvec=dir, yvec=vrhs, alpha=step)
      ! Update residual res = res - step * amat * dir
      call axpy(xvec=adir, yvec=res, alpha=-step)

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

      ! Stop the iteration timer
      if (self%verbosity > 1) call timer%pop

      ! Print iteration progress
      call print_cg_iteration(unit, it, sqrt(resnorm), step, sqrt(rel_resnorm), &
         & self%verbosity, timer)

      if (it == maxit) then
         if (self%verbosity > 1) then
            ! Stop the iteration timer
            call timer%pop
            ! Stop the total timer
            call timer%pop
         end if
         call fatal_error(error, "CG did not converge within max iterations.")
         return
      end if

   end do

   ! Stop the total timer
   if (self%verbosity > 1) call timer%pop

   ! Print final summary
   call print_cg_final(unit, timer, self%verbosity)

end subroutine solve


!> Print header for CG solver
subroutine print_cg_header(unit, verbosity, maxit, tol, timer)

   !> Output unit
   integer, intent(in) :: unit

   !> Verbosity level
   integer, intent(in) :: verbosity

   !> Maximum number of iterations
   integer, intent(in) :: maxit

   !> Convergence tolerance
   real(wp), intent(in) :: tol

   !> Timer holding the accumulated initialization time
   type(timer_type), intent(in), optional :: timer

   if (verbosity > 1) then
      write(unit, '(a)') "Using Conjugate Gradient Solver"
      write(unit, '(a)')
      write(unit, '(a, 1x, i6)') "Max iterations : ", maxit
      write(unit, '(a, 1x, es10.2)') "Tolerance      : ", tol
      write(unit, '(a)') "Preconditioner : Jacobi (Diagonal)"
      write(unit, '(a, 1x, a)') "Initialisation time:", &
         & format_time(timer%get("initialization"))
      write(unit, '(a)') ''
      write(unit, '(2X,A,6X,A,8X,A,6X,A,4X,A)') &
         & 'iter', '|residual|', 'step', 'relative residual', 'Time / s'
   else if (verbosity == 1) then
      write(unit, '(a)') "Using Conjugate Gradient Solver"
      write(unit, '(a)')
      write(unit, '(a, 1x, i6)') "Max iterations : ", maxit
      write(unit, '(a, 1x, es10.2)') "Tolerance      : ", tol
      write(unit, '(a)') "Preconditioner : Jacobi (Diagonal)"
      write(unit, '(a)') ''
      write(unit, '(2X,A,6X,A,8X,A,6X,A)') &
         & 'iter', '|residual|', 'step', 'relative residual'
   end if

end subroutine print_cg_header


!> Print convergence message
subroutine print_cg_convergence(unit, iter, res_norm, verbosity)

   !> Output unit
   integer, intent(in) :: unit

   !> Iteration at which convergence was reached
   integer, intent(in) :: iter

   !> Converged residual norm
   real(wp), intent(in) :: res_norm

   !> Verbosity level
   integer, intent(in) :: verbosity

   if (verbosity > 0) then
      write(unit, '(a)') ''
      write(unit, '(a, i0, a, es15.5)') &
         "CG converged in ", iter, " iterations with residual norm ", res_norm
      write(unit, '(a)') ''
   end if

end subroutine print_cg_convergence


!> Print iteration progress
subroutine print_cg_iteration(unit, iter, res_norm, step, rel_resnorm, &
      & verbosity, timer)

   !> Output unit
   integer, intent(in) :: unit

   !> Current iteration number
   integer, intent(in) :: iter

   !> Current residual norm
   real(wp), intent(in) :: res_norm

   !> Current step length
   real(wp), intent(in) :: step

   !> Current relative residual norm
   real(wp), intent(in) :: rel_resnorm

   !> Verbosity level
   integer, intent(in) :: verbosity

   !> Timer holding the accumulated iteration time
   type(timer_type), intent(in), optional :: timer

   if (verbosity == 1) then
      write(unit, '(i6,*(1x, es15.5))') iter, res_norm, step, rel_resnorm
   else if (verbosity > 1) then
      write(unit, '(i6,*(1x, es15.5))') iter, res_norm, step, rel_resnorm, &
         & timer%get("iteration")
   end if

end subroutine print_cg_iteration


!> Print the final timing summary
subroutine print_cg_final(unit, timer, verbosity)

   !> Output unit
   integer, intent(in) :: unit

   !> Timer holding the accumulated execution time
   type(timer_type), intent(in) :: timer

   !> Verbosity level
   integer, intent(in) :: verbosity

   if (verbosity > 1) then
      write(unit, '(a, 1x, a)') "CG total time : ", &
         & format_time(timer%get("total"))
      write(unit, '(a)') ''
   end if

end subroutine print_cg_final


end module multicharge_solver_cg
