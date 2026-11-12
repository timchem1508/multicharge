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
   use, intrinsic :: iso_fortran_env, only : output_unit
   use mctc_env, only : error_type, fatal_error, timer_type, wp
   use multicharge_blas, only : axpy, dot, scal
   use multicharge_output, only : print_cg_convergence, print_cg_final, &
      & print_cg_header, print_cg_iteration
   use multicharge_solver_precond, only : jacobi_precond, new_jacobi_precond, &
      & precond_type
   use multicharge_solver_trialprod, only : trialprod
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

      !> Enforce the total charge by projection in a single CG solve
      logical, allocatable :: projected

      !> Preconditioner built by the caller, see precond_kind
      integer, allocatable :: precond

      !> Neighbour cutoff of the sparse approximate inverse in Bohr
      real(wp), allocatable :: precond_cutoff

   end type cg_input

   !> Preconditioned conjugate-gradient solver
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
   if (allocated(input%projected)) self%projected = input%projected
   if (allocated(input%precond)) self%precond = input%precond
   if (allocated(input%precond_cutoff)) self%precond_cutoff = input%precond_cutoff

end subroutine new_cg_solver


!> Solve a linear system with a preconditioned CG method.
!>
!> With the total sum of the solution given, the system A v + lambda 1 = b with
!> sum(v) = qtot is solved instead. The iterates are kept on the constraint
!> plane and the preconditioner is projected onto it. The constant part of the
!> residual is irrelevant for the projected iteration, it is removed after each
!> update to avoid cancellation and accumulated into the Lagrange multiplier.
subroutine solve(self, tprod, xvec, vrhs, ainv, qtot, lambda, prec, new_unit, &
   & error)

   !> Conjugate-gradient solver instance
   class(cg_solver), intent(in) :: self

   !> Trial product operator representing the linear system matrix
   class(trialprod), intent(in) :: tprod

   !> Right-hand side vector
   real(wp), intent(in) :: xvec(:)

   !> On input: initial guess; on output: solution
   real(wp), intent(inout), contiguous :: vrhs(:)

   !> Inverse matrix, not available for the iterative solver
   real(wp), intent(out), optional :: ainv(:, :)

   !> Total sum of the solution, enforced by projection
   real(wp), intent(in), optional :: qtot

   !> Lagrange multiplier of the projected system
   real(wp), intent(out), optional :: lambda

   !> Preconditioner, the Jacobi preconditioner is used if absent
   class(precond_type), intent(in), optional, target :: prec

   !> Output unit
   integer, intent(in), optional :: new_unit

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: maxit, it, ndim, unit
   logical :: converged, projected
   real(wp) :: tol, tol_square, xvecnorm, resnorm, rel_resnorm
   real(wp) :: denom, step, updfact, resdot_old, resdot_new, precsum, lamacc
   real(wp), allocatable :: dir(:), res(:), precres(:), adir(:), pone(:)
   type(jacobi_precond), target :: jacobi
   class(precond_type), pointer :: mprec
   type(timer_type) :: timer

   ! CG cannot compute the inverse matrix
   if (present(ainv)) then
      call fatal_error(error, &
         & "The inverse matrix cannot be calculated using an iterative solver.")
      return
   end if

   if (self%cgmiter <= 0) then
      call fatal_error(error, "Maximal number of CG iterations must be positive.")
      return
   end if

   if (.not. self%cgtol > 0.0_wp) then
      call fatal_error(error, "CG convergence tolerance must be positive.")
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

   projected = present(qtot)

   tol = self%cgtol
   tol_square = tol**2
   maxit = self%cgmiter

   allocate(res(ndim), dir(ndim), precres(ndim), adir(ndim), pone(ndim))

   if (self%verbosity > 1) call timer%push("total")
   if (self%verbosity > 1) call timer%push("initialization")

   ! Preconditioner, the Jacobi preconditioner is built if none is given
   if (present(prec)) then
      mprec => prec
   else
      call new_jacobi_precond(jacobi, tprod, ndim)
      mprec => jacobi
   end if

   ! Preconditioned unit vector spans the correction to the constraint plane
   precres(:) = 1.0_wp
   call mprec%apply(precres, pone)
   precsum = sum(pone)

   ! Move the initial guess onto the constraint plane along the preconditioner
   if (projected) then
      vrhs(:) = vrhs + (qtot - sum(vrhs)) / precsum * pone
   end if

   ! Initial residual
   call tprod%multiply(vrhs, adir)
   res(:) = xvec(:) - adir(:)
   lamacc = 0.0_wp
   call remove_constant(res)

   ! Initial preconditioned residual precres = M^-1 * res
   call precondition(res, precres)
   dir(:) = precres(:)

   ! Initial norm
   xvecnorm = dot(xvec, xvec)
   if (xvecnorm < tol_square) xvecnorm = 1.0_wp

   ! Initial resnorm and res^T * precres
   resnorm = residual_norm(res)
   resdot_old = dot(res, precres)
   converged = resnorm / xvecnorm <= tol_square

   if (self%verbosity > 1) call timer%pop

   ! Print header
   call print_cg_header(unit, self%verbosity, maxit, tol, mprec%get_name(), timer)

   ! Main CG iteration loop
   it = 0
   if (.not. converged) then
      do it = 1, maxit

         if (self%verbosity > 1) call timer%push("iteration")

         ! Matrix-vector product
         call tprod%multiply(dir, adir)

         ! Non-positive curvature breaks down the iteration, unless the
         ! search direction vanished for an already converged residual
         denom = dot(dir, adir)
         if (.not. denom > 0.0_wp) then
            if (self%verbosity > 1) call timer%pop
            exit
         end if

         ! Step length step = (res^T * precres) / (dir^T * amat * dir)
         step = resdot_old / (denom + eps)

         ! Update solution vrhs = vrhs + step * dir
         call axpy(xvec=dir, yvec=vrhs, alpha=step)
         ! Update residual res = res - step * amat * dir
         call axpy(xvec=adir, yvec=res, alpha=-step)
         call remove_constant(res)

         ! Compute the new residual norm
         resnorm = residual_norm(res)

         ! Relative residual norm to check convergence
         rel_resnorm = resnorm / xvecnorm

         converged = rel_resnorm <= tol_square
         if (converged) then
            if (self%verbosity > 1) call timer%pop
            call print_cg_convergence(unit, it, sqrt(resnorm), self%verbosity)
            exit
         end if

         ! Updated preconditioned residual
         call precondition(res, precres)

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

      end do
   end if

   ! Stop the total timer
   if (self%verbosity > 1) call timer%pop

   if (.not. converged) then
      if (it > maxit) then
         call fatal_error(error, "CG did not converge within max iterations.")
      else
         call fatal_error(error, "CG breakdown: non-positive curvature, "//&
            & "the coefficient matrix is not positive definite.")
      end if
      return
   end if

   ! The converged residual is the Lagrange multiplier times the unit vector
   if (present(lambda)) lambda = lamacc

   ! Print final summary
   call print_cg_final(unit, timer, self%verbosity)

contains

   !> Apply the preconditioner, projected onto the constraint plane if requested
   subroutine precondition(rv, zv)
      !> Residual vector
      real(wp), intent(in) :: rv(:)
      !> Preconditioned residual
      real(wp), intent(out) :: zv(:)

      call mprec%apply(rv, zv)
      if (projected) zv(:) = zv - sum(zv) / precsum * pone
   end subroutine precondition

   !> Squared norm of the residual
   function residual_norm(rv) result(rnorm)
      !> Residual vector
      real(wp), intent(in) :: rv(:)
      !> Squared residual norm
      real(wp) :: rnorm

      rnorm = dot(rv, rv)
   end function residual_norm

   !> Move the constant part of the residual into the Lagrange multiplier
   subroutine remove_constant(rv)
      !> Residual vector
      real(wp), intent(inout) :: rv(:)

      real(wp) :: shift

      if (.not. projected) return
      shift = sum(rv) / real(ndim, wp)
      rv(:) = rv - shift
      lamacc = lamacc + shift
   end subroutine remove_constant

end subroutine solve


end module multicharge_solver_cg
