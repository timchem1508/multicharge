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
!> Provides implementation of the (block) conjugate gradient solver

module multicharge_solver_cg
   use iso_fortran_env, only : output_unit
   use mctc_env, only : error_type, fatal_error, format_time, timer_type, wp, i8
   use mctc_csrlist, only : csr_list, find_lower_start, spsymv_csr, spmm_csr
   use mctc_io_structure_atomsort, only : get_spatial_order
   use multicharge_blas, only : axpy, dot, gemm, scal, symv
   use multicharge_lapack, only : potrf, potrs, syevd
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

      !> Convergence tolerance of the block CG for the inverse on a neighborlist
      real(wp), allocatable :: invtol

      !> Output verbosity
      integer, allocatable :: verbosity

      !> Whether to use a neighborlist representation
      logical, allocatable :: use_nlist
   end type cg_input

   !> Conjugate-gradient solver with a Jacobi preconditioner
   type, extends(mchrg_solver_type) :: cg_solver
      !> Maximum number of iterations
      integer, allocatable :: cgmiter

      !> Convergence tolerance
      real(wp), allocatable :: cgtol

      !> Convergence tolerance of the block CG for the inverse on a neighborlist
      real(wp), allocatable :: invtol

      !> Output verbosity
      integer, allocatable :: verbosity

      !> Whether to use a neighborlist representation
      logical, allocatable :: use_nlist
   contains
      !> Solve the linear system iteratively
      procedure :: solve

      !> Solve the linear system for a block of right-hand sides
      procedure :: solve_block

      !> Invert the dense coefficient matrix block-wise
      procedure :: invert

      !> Invert the coefficient matrix on the pattern of its neighborlist
      procedure :: invert_list
   end type cg_solver

   !> Positive number used to prevent division by zero
   real(wp), parameter :: eps = tiny(1.0_wp)

   !> Default maximum number of iterations
   integer, parameter :: cgmiter_def = 1000

   !> Default convergence tolerance
   real(wp), parameter :: cgtol_def = 1.0e-15_wp

   !> Default convergence tolerance of the inverse on a neighborlist
   real(wp), parameter :: invtol_def = 1.0e-11_wp

   !> Default output verbosity
   integer, parameter :: verbosity_def = 0

   !> Default neighborlist usage
   logical, parameter :: use_nlist_def = .false.

   !> Default number of columns per block of right-hand sides
   integer, parameter :: block_size_def = 16

   !> Relative eigenvalue threshold of the Gram matrix
   real(wp), parameter :: orth_thr = 1.0e3_wp * epsilon(1.0_wp)


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
   if (allocated(input%invtol)) then
      self%invtol = input%invtol
   else
      self%invtol = invtol_def
   end if
   if (allocated(input%verbosity)) then
      self%verbosity = input%verbosity
   else
      self%verbosity = verbosity_def
   end if
   if (allocated(input%use_nlist)) then
      self%use_nlist = input%use_nlist
   else
      self%use_nlist = use_nlist_def
   end if

end subroutine new_cg_solver


!> Solve a linear system with a diagonally preconditioned CG method
subroutine solve(self, amat, alist, xvec, vrhs, ainv, cpq, list, new_unit, error)
   !> Conjugate-gradient solver instance
   class(cg_solver), intent(in) :: self

   !> Dense coefficient matrix of the linear system
   real(wp), intent(in), optional :: amat(:, :)

   !> Coefficient matrix values on the upper triangle of the neighborlist
   real(wp), intent(in), optional :: alist(:)

   !> Right-hand side vector
   real(wp), intent(in) :: xvec(:)

   !> On input: initial guess; on output: solution
   real(wp), intent(inout), contiguous :: vrhs(:)

   !> Inverse matrix
   real(wp), intent(out), optional :: ainv(:, :)

   !> Whether to solve coupled-perturbed equations
   logical, intent(in), optional :: cpq

   !> Optional sorted complete neighborlist representation of the matrix
   type(csr_list), intent(in), optional :: list

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
   real(wp), allocatable :: precres(:)
   ! amat-dir product
   real(wp), allocatable :: Adir(:)
   ! Denominator of the step length
   real(wp) :: denom
   ! Step length
   real(wp) :: step
   ! Update factor for search direction
   real(wp) :: updfact
   ! Projection of preconditioned residual and an original one
   real(wp) :: resdot_old, resdot_new
   ! Relative residual norm (|resnorm| / |vrhs|)
   real(wp) :: rel_resnorm

   type(timer_type) :: timer
   integer :: unit
   logical :: nlist

   nlist = self%use_nlist .and. .not. present(amat) .and. &
      & present(list) .and. present(alist)

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
   if (.not. nlist) then
      if (.not. present(amat)) then
         call fatal_error(error, "No coefficient matrix provided.")
         return
      end if
      if (size(amat, 1) /= ndim .or. size(amat, 2) /= ndim) then
         call fatal_error(error, "dimension mismatch.")
         return
      end if
   end if

   tol = self%cgtol
   tol_square = tol**2
   maxit = self%cgmiter

   allocate(res(ndim), dir(ndim), precres(ndim), Adir(ndim), prec(ndim))

   if (self%verbosity > 1) call timer%push("total")
   if (self%verbosity > 1) call timer%push("initialization")

   ! Diagonal preconditioner
   if (nlist) then
      do iat = 1, ndim
         prec(iat) = 1.0_wp / (alist(list%inl(iat)) + eps)
      end do
   else
      do iat = 1, ndim
         prec(iat) = 1.0_wp / (amat(iat, iat) + eps)
      end do
   end if

   ! Initial residual
   if (nlist) then
      call spsymv_csr(ndim, alist, list%inl, list%nlat, vrhs, Adir, uplo='U')
   else
      call symv(amat, vrhs, Adir, alpha=1.0_wp, beta=0.0_wp)
   end if
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
      if (nlist) then
         call spsymv_csr(ndim, alist, list%inl, list%nlat, dir, Adir, uplo='U')
      else
         call symv(amat, dir, Adir, alpha=1.0_wp, beta=0.0_wp)
      end if

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
      call print_cg_iteration(unit, it, sqrt(resnorm), step, sqrt(rel_resnorm), &
         & self%verbosity, timer)

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


!> Solve a linear system with multiple right-hand sides by block CG
subroutine solve_block(self, amat, alist, bmat, xmat, list, new_unit, error, &
   & niter)
   !> Conjugate-gradient solver instance
   class(cg_solver), intent(in) :: self

   !> Dense coefficient matrix of the linear system
   real(wp), intent(in), optional :: amat(:, :)

   !> Coefficient matrix values on the upper triangle of the neighborlist
   real(wp), intent(in), optional :: alist(:)

   !> Block of right-hand sides
   real(wp), intent(in) :: bmat(:, :)

   !> On input: initial guess; on output: solution
   real(wp), intent(inout) :: xmat(:, :)

   !> Optional sorted complete neighborlist representation of the matrix
   type(csr_list), intent(in), optional :: list

   !> Output unit
   integer, intent(in), optional :: new_unit

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Number of iterations until convergence
   integer, intent(out), optional :: niter

   ! Maximal number of iterations
   integer :: maxit
   ! Tolerance of the solver
   real(wp) :: tol, tol_square
   ! Counters
   integer :: it, iat, ivec
   ! Size of the system, number of right-hand sides, rank of the search space
   integer :: ndim, nrhs, nrank
   ! Search directions P
   real(wp), allocatable :: dir(:, :)
   ! Squared norms of the right-hand sides
   real(wp), allocatable :: bnorm(:)
   ! Residuals R
   real(wp), allocatable :: res(:, :)
   ! Squared residual norms
   real(wp), allocatable :: resnorm(:)
   ! Diagonal preconditioner
   real(wp), allocatable :: prec(:)
   ! Preconditioned residuals Z
   real(wp), allocatable :: precres(:, :)
   ! amat-dir product Q = A P
   real(wp), allocatable :: adir(:, :)
   ! Cholesky factor of P^T A P
   real(wp), allocatable :: dtad(:, :)
   ! Step lengths alpha = (P^T Q)^-1 P^T R
   real(wp), allocatable :: step(:, :)
   ! Update factors beta = (P^T Q)^-1 Q^T Z
   real(wp), allocatable :: updfact(:, :)
   ! Largest relative residual norm over all right-hand sides
   real(wp) :: rel_resnorm

   type(timer_type) :: timer
   integer :: unit, info
   logical :: nlist

   nlist = self%use_nlist .and. .not. present(amat) .and. &
      & present(list) .and. present(alist)
   if (present(niter)) niter = 0

   if (present(new_unit)) then
      unit = new_unit
   else
      unit = output_unit
   end if

   ! Dimensions check
   ndim = size(bmat, 1)
   nrhs = size(bmat, 2)
   if (size(xmat, 1) /= ndim .or. size(xmat, 2) /= nrhs) then
      call fatal_error(error, "Dimension mismatch between bmat and xmat.")
      return
   end if
   if (.not. nlist) then
      if (.not. present(amat)) then
         call fatal_error(error, "No coefficient matrix provided.")
         return
      end if
      if (size(amat, 1) /= ndim .or. size(amat, 2) /= ndim) then
         call fatal_error(error, "dimension mismatch.")
         return
      end if
   end if
   if (nrhs == 0) return

   tol = self%cgtol
   tol_square = tol**2
   maxit = self%cgmiter

   allocate(res(ndim, nrhs), dir(ndim, nrhs), precres(ndim, nrhs), &
      & adir(ndim, nrhs), prec(ndim), bnorm(nrhs), resnorm(nrhs), &
      & dtad(nrhs, nrhs), step(nrhs, nrhs), updfact(nrhs, nrhs))

   if (self%verbosity > 1) call timer%push("total")
   if (self%verbosity > 1) call timer%push("initialization")

   ! Diagonal preconditioner
   if (nlist) then
      do iat = 1, ndim
         prec(iat) = 1.0_wp / (alist(list%inl(iat)) + eps)
      end do
   else
      do iat = 1, ndim
         prec(iat) = 1.0_wp / (amat(iat, iat) + eps)
      end do
   end if

   ! Initial residual R = B - A X
   call block_matmul(nlist, xmat, res, amat, alist, list)
   res(:, :) = bmat - res

   ! Initial norms
   do ivec = 1, nrhs
      bnorm(ivec) = dot(bmat(:, ivec), bmat(:, ivec))
      if (bnorm(ivec) < tol_square) bnorm(ivec) = 1.0_wp
      resnorm(ivec) = dot(res(:, ivec), res(:, ivec))
   end do
   rel_resnorm = maxval(resnorm / bnorm)

   ! Initial search directions P = orth(M^-1 R)
   do ivec = 1, nrhs
      precres(:, ivec) = prec * res(:, ivec)
   end do
   call orthonormalize(precres, dir, nrank, error)
   if (allocated(error)) return

   if (self%verbosity > 1) call timer%pop

   ! Print header
   call print_block_cg_header(unit, self%verbosity, maxit, tol, nrhs, timer)

   ! Initial guess is already converged
   if (rel_resnorm <= tol_square .or. nrank == 0) then
      if (self%verbosity > 0) then
         call print_cg_convergence(unit, 0, sqrt(maxval(resnorm)), &
            & self%verbosity)
      end if
      if (self%verbosity > 1) call timer%pop
      call print_cg_final(unit, timer, self%verbosity)
      return
   end if

   ! Main block CG iteration loop

   do it = 1, maxit

      if (self%verbosity > 1) call timer%push("iteration")

      ! Matrix-block product Q = A P
      call block_matmul(nlist, dir(:, :nrank), adir(:, :nrank), amat, alist, &
         & list)

      ! Shrink the small matrices to the current rank to keep them contiguous
      if (size(dtad, 1) /= nrank) then
         deallocate(dtad, step, updfact)
         allocate(dtad(nrank, nrank), step(nrank, nrhs), updfact(nrank, nrhs))
      end if

      ! Cholesky factorization of P^T Q
      call gemm(dir(:, :nrank), adir(:, :nrank), dtad, transa='t')
      call potrf(dtad, info=info)
      if (info /= 0) then
         if (self%verbosity > 1) call timer%pop
         if (self%verbosity > 1) call timer%pop
         call fatal_error(error, &
            & "Block CG: P^T A P is not positive definite.")
         return
      end if

      ! Step lengths alpha = (P^T Q)^-1 (P^T R)
      call gemm(dir(:, :nrank), res, step, transa='t')
      call potrs(dtad, step)

      ! Update solution X = X + P alpha
      call gemm(dir(:, :nrank), step, xmat, beta=1.0_wp)
      ! Update residual R = R - Q alpha
      call gemm(adir(:, :nrank), step, res, alpha=-1.0_wp, beta=1.0_wp)

      ! Compute the new residual norms
      do ivec = 1, nrhs
         resnorm(ivec) = dot(res(:, ivec), res(:, ivec))
      end do

      ! Largest relative residual norm to check convergence
      rel_resnorm = maxval(resnorm / bnorm)

      if (rel_resnorm <= tol_square) then
         if (self%verbosity > 1) call timer%pop
         if (self%verbosity > 0) then
            call print_cg_convergence(unit, it, sqrt(maxval(resnorm)), &
               & self%verbosity)
         end if
         if (present(niter)) niter = it
         exit
      end if

      ! Updated preconditioned residuals Z = M^-1 R
      do ivec = 1, nrhs
         precres(:, ivec) = prec * res(:, ivec)
      end do

      ! Update factors (P^T Q)^-1 (Q^T Z), beta is their negative
      call gemm(adir(:, :nrank), precres, updfact, transa='t')
      call potrs(dtad, updfact)

      ! Update search directions P = orth(Z + P beta)
      call gemm(dir(:, :nrank), updfact, precres, alpha=-1.0_wp, beta=1.0_wp)
      call orthonormalize(precres, dir, nrank, error)
      if (allocated(error)) then
         if (self%verbosity > 1) call timer%pop
         if (self%verbosity > 1) call timer%pop
         return
      end if

      ! iteration timer pop
      if (self%verbosity > 1) call timer%pop

      ! Print iteration progress
      call print_block_cg_iteration(unit, it, sqrt(maxval(resnorm)), nrank, &
         & sqrt(rel_resnorm), self%verbosity, timer)

      if (nrank == 0) then
         if (self%verbosity > 1) call timer%pop
         call fatal_error(error, "Block CG: search space collapsed.")
         return
      end if

      if (it == maxit) then
         if (self%verbosity > 1) call timer%pop
         call fatal_error(error, &
            & "Block CG did not converge within max iterations.")
         return
      end if

   end do

   ! pop total
   if (self%verbosity > 1) call timer%pop

   ! Print final summary
   call print_cg_final(unit, timer, self%verbosity)

end subroutine solve_block


!> Invert a dense symmetric positive definite matrix by block CG
subroutine invert(self, amat, ainv, new_unit, error)
   !> Conjugate-gradient solver instance
   class(cg_solver), intent(in) :: self

   !> Dense coefficient matrix of the linear system
   real(wp), intent(in) :: amat(:, :)

   !> Inverse of the coefficient matrix
   real(wp), intent(out) :: ainv(:, :)

   !> Output unit
   integer, intent(in), optional :: new_unit

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   ! Size of the system, offset and number of columns of the current block
   integer :: ndim, ioff, ncol
   ! Number of trailing rows
   integer :: nsub
   ! Counters
   integer :: iat, ivec
   ! Diagonal of the coefficient matrix
   real(wp), allocatable :: diag(:)
   ! Known rows X_K of the block columns
   real(wp), allocatable :: xknown(:, :)
   ! Right-hand sides and solutions of the trailing system
   real(wp), allocatable :: bsub(:, :), xsub(:, :)
   ! Storage of the dense trailing matrix A_UU
   real(wp), allocatable, target :: abuf(:)
   real(wp), pointer :: asub(:, :)

   ! Dimensions check
   ndim = size(ainv, 1)
   if (size(ainv, 2) /= ndim) then
      call fatal_error(error, "Inverse matrix must be square.")
      return
   end if
   if (size(amat, 1) /= ndim .or. size(amat, 2) /= ndim) then
      call fatal_error(error, "dimension mismatch.")
      return
   end if
   if (ndim == 0) return

   allocate(diag(ndim), xknown(ndim, block_size_def), &
      & abuf(int(ndim, i8) * ndim))
   do iat = 1, ndim
      diag(iat) = amat(iat, iat) + eps
   end do

   do ioff = 0, ndim - 1, block_size_def
      ncol = min(block_size_def, ndim - ioff)
      nsub = ndim - ioff

      ! Rows of the block columns known from earlier columns, A^-1 = A^-T
      do ivec = 1, ncol
         xknown(:ioff, ivec) = ainv(ioff + ivec, :ioff)
         ainv(:ioff, ioff + ivec) = xknown(:ioff, ivec)
      end do

      ! Trailing system A_UU X_U = I_U - A_UK X_K
      allocate(bsub(nsub, ncol), source=0.0_wp)
      do ivec = 1, ncol
         bsub(ivec, ivec) = 1.0_wp
      end do
      asub(1:nsub, 1:nsub) => abuf(1:int(nsub, i8) * nsub)
      asub(:, :) = amat(ioff + 1:, ioff + 1:)
      if (ioff > 0) then
         call gemm(amat(ioff + 1:, :ioff), xknown(:ioff, :ncol), bsub, &
            & alpha=-1.0_wp, beta=1.0_wp)
      end if

      ! Initial guess from the diagonal of the trailing matrix
      allocate(xsub(nsub, ncol))
      do ivec = 1, ncol
         xsub(:, ivec) = bsub(:, ivec) / diag(ioff + 1:)
      end do

      call self%solve_block(amat=asub, bmat=bsub, xmat=xsub, &
         & new_unit=new_unit, error=error)
      if (allocated(error)) return

      ainv(ioff + 1:, ioff + 1:ioff + ncol) = xsub
      deallocate(bsub, xsub)
   end do

end subroutine invert


!> Invert a symmetric positive definite matrix on a neighborlist
subroutine invert_list(self, alist, list, ainvlist, xyz, new_unit, error, half)
   !> Conjugate-gradient solver instance
   class(cg_solver), intent(in) :: self

   !> Coefficient matrix values on the upper triangle of the list of the values
   real(wp), intent(in) :: alist(:)

   !> Sorted complete neighborlist, its complete rows select the subsystems
   type(csr_list), intent(in), target :: list

   !> Inverse on the upper triangle of the list of the values
   real(wp), intent(out) :: ainvlist(:)

   !> Cartesian coordinates for spatially compact blocks
   real(wp), intent(in), optional :: xyz(:, :)

   !> Output unit
   integer, intent(in), optional :: new_unit

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Sorted list of the values, each row starts with the upper triangle of the
   !> same row of the complete list, absent for the complete list itself
   type(csr_list), intent(in), optional, target :: half

   ! Size of the system, number of blocks and counters
   integer :: ndim, nblk, iblk, iat, max_nnz
   ! Atoms ordered into blocks of columns
   integer, allocatable :: order(:)
   ! List of the values and the end of the upper triangle in each of its rows
   type(csr_list), pointer :: hlist
   integer(i8), allocatable :: hfin(:)
   ! Per-thread workspace for the subsystem matrix and block CG vectors
   integer, allocatable :: loc(:), col_id(:), sidx(:), subinl(:), subjat(:)
   real(wp), allocatable :: subdiag(:), prec(:), subval(:)
   real(wp), allocatable :: work_x(:), work_r(:), work_z(:), work_p(:), work_q(:)
   ! Size of the subsystem and number of iterations of a block
   integer :: nsub, niter, maxsub, maxiter
   integer(i8) :: sumsub, sumiter
   real(wp) :: tblk(3), sumtime(3)
   ! Solver for the blocks without output
   type(cg_solver) :: blksolver
   integer :: unit
   logical :: failed, skip

   if (.not. self%use_nlist) then
      call fatal_error(error, "Inversion on a neighborlist requires use_nlist.")
      return
   end if
   if (.not. list%complete) then
      call fatal_error(error, &
         & "Inversion on a neighborlist requires a sorted complete list.")
      return
   end if
   if (allocated(list%nltr)) then
      call fatal_error(error, &
         & "Inversion on a neighborlist requires one entry per atom pair.")
      return
   end if
   if (present(half)) then
      hlist => half
   else
      hlist => list
   end if
   ndim = size(list%inl) - 1
   if (size(hlist%inl) /= ndim + 1) then
      call fatal_error(error, "Dimension mismatch between the lists of the inverse.")
      return
   end if
   if (size(ainvlist, kind=i8) /= size(hlist%nlat, kind=i8) &
      & .or. size(alist, kind=i8) < size(hlist%nlat, kind=i8)) then
      call fatal_error(error, "Dimension mismatch between list and ainvlist.")
      return
   end if
   if (present(xyz)) then
      if (size(xyz, 1) /= 3 .or. size(xyz, 2) < ndim) then
         call fatal_error(error, "Dimension mismatch between list and xyz.")
         return
      end if
   end if
   if (ndim <= 0) return
   ainvlist(:) = 0.0_wp

   if (present(new_unit)) then
      unit = new_unit
   else
      unit = output_unit
   end if

   ! Blocks of spatially close atoms
   if (present(xyz)) then
      call get_spatial_order(xyz(:, :ndim), order, block_size_def)
   else
      order = [(iat, iat = 1, ndim)]
   end if
   nblk = (ndim + block_size_def - 1) / block_size_def

   ! Precompute upper-triangle row bounds of the values once
   allocate(hfin(ndim))
   max_nnz = 0
   do iat = 1, ndim
      hfin(iat) = find_lower_start(hlist, iat) - 1_i8
      max_nnz = max_nnz + int(hfin(iat) - hlist%inl(iat))
   end do

   call new_cg_solver(blksolver, cg_input(cgmiter=self%cgmiter, &
      & cgtol=self%invtol, verbosity=0, use_nlist=.true.))

   failed = .false.
   sumsub = 0_i8
   sumiter = 0_i8
   sumtime(:) = 0.0_wp
   maxsub = 0
   maxiter = 0
   !$omp parallel default(none) &
   !$omp shared(ndim, nblk, max_nnz, order, hfin, alist, list, hlist, ainvlist) &
   !$omp shared(blksolver, error, failed) &
   !$omp private(iblk, loc, col_id, sidx, subdiag, prec, subinl, subjat, subval) &
   !$omp private(work_x, work_r, work_z, work_p, work_q, nsub, niter, skip, tblk) &
   !$omp reduction(+:sumsub, sumiter, sumtime) reduction(max:maxsub, maxiter)
   allocate(loc(ndim), col_id(ndim), source=0)
   allocate(sidx(ndim), subdiag(ndim), prec(ndim), subinl(ndim + 1))
   allocate(subjat(max(1, max_nnz)), subval(max(1, max_nnz)))
   allocate(work_x(block_size_def * ndim), work_r(block_size_def * ndim), &
      & work_z(block_size_def * ndim), work_p(block_size_def * ndim), &
      & work_q(block_size_def * ndim))
   !$omp do schedule(dynamic)
   do iblk = 1, nblk
      !$omp atomic read
      skip = failed
      if (skip) cycle
      call invert_block(blksolver, alist, list, hlist, hfin, &
         & order((iblk - 1) * block_size_def + 1:min(iblk * block_size_def, ndim)), &
         & loc, col_id, sidx, subdiag, prec, subinl, subjat, subval, &
         & work_x, work_r, work_z, work_p, work_q, ainvlist, nsub, niter, &
         & tblk, error, failed)
      sumsub = sumsub + nsub
      sumtime(:) = sumtime + tblk
      sumiter = sumiter + niter
      maxsub = max(maxsub, nsub)
      maxiter = max(maxiter, niter)
   end do
   !$omp end do
   deallocate(loc, col_id, sidx, subdiag, prec, subinl, subjat, subval, &
      & work_x, work_r, work_z, work_p, work_q)
   !$omp end parallel
   if (allocated(error)) return

   if (self%verbosity > 0) then
      write(unit, '(a, 1x, i0, a, i0, a)') "Inverse on neighborlist blocks   :", &
         & nblk, " of ", min(block_size_def, ndim), " columns"
      write(unit, '(a, 1x, f0.1, a, i0)') "Subsystem atoms per block        :", &
         & real(sumsub, wp) / nblk, ", max ", maxsub
      write(unit, '(a, 1x, f0.1, a, i0)') "Block CG iterations per block    :", &
         & real(sumiter, wp) / nblk, ", max ", maxiter
      write(unit, '(a, 1x, es9.2)') "Block CG tolerance               :", &
         & self%invtol
      write(unit, '(a, 3(1x, f0.3))') "Select/CG/store time (s)         :", &
         & sumtime
      write(unit, '(a)') ''
   end if

end subroutine invert_list


!> Solve a block of columns of the inverse on its subsystem
subroutine invert_block(solver, alist, list, hlist, hfin, cols, loc, col_id, sidx, &
   & subdiag, prec, subinl, subjat, subval, work_x, work_r, work_z, work_p, &
   & work_q, ainvlist, nsub, niter, times, error, failed)
   !> Conjugate-gradient solver for the block
   type(cg_solver), intent(in) :: solver

   !> Coefficient matrix values on the upper triangle of the list of the values
   real(wp), intent(in) :: alist(:)

   !> Sorted complete neighborlist, its complete rows select the subsystem
   type(csr_list), intent(in) :: list

   !> Sorted list of the values, each row starts with the upper triangle of the
   !> same row of the complete list
   type(csr_list), intent(in) :: hlist

   !> End of the upper triangle in each row of the list of the values
   integer(i8), intent(in) :: hfin(:)

   !> Atoms of the block columns
   integer, intent(in) :: cols(:)

   !> Local index of the atoms in the subsystem, zero on entry and exit
   integer, intent(inout) :: loc(:)

   !> Column index (1..ncol) of subsystem atoms, zero on entry and exit
   integer, intent(inout) :: col_id(:)

   !> Atoms of the subsystem
   integer, intent(inout) :: sidx(:)

   !> Preallocated subsystem diagonal buffer
   real(wp), intent(inout) :: subdiag(:)

   !> Preallocated diagonal preconditioner buffer
   real(wp), intent(inout) :: prec(:)

   !> Preallocated subsystem row-offset buffer
   integer, intent(inout) :: subinl(:)

   !> Preallocated subsystem column-index buffer
   integer, intent(inout) :: subjat(:)

   !> Preallocated subsystem matrix-value buffer
   real(wp), intent(inout) :: subval(:)

   !> Preallocated block-CG vector buffers
   real(wp), target, intent(inout) :: work_x(:), work_r(:), work_z(:)
   real(wp), target, intent(inout) :: work_p(:), work_q(:)

   !> Inverse on the list of the values, the rows of the block atoms are set
   real(wp), intent(inout) :: ainvlist(:)

   !> Number of atoms in the subsystem
   integer, intent(out) :: nsub

   !> Number of block CG iterations
   integer, intent(out) :: niter

   !> Times of the subsystem selection, block CG and storage
   real(wp), intent(out) :: times(3)

   !> Error handling, shared by all blocks
   type(error_type), allocatable, intent(inout) :: error

   !> Whether any block failed, shared by all blocks
   logical, intent(inout) :: failed

   integer :: ncol, isub, iat, jat, ivec
   integer(i8) :: kat, tick(4), rate
   real(wp), pointer :: xsub(:, :), rsub(:, :)
   type(error_type), allocatable :: blkerror

   call system_clock(tick(1), rate)
   ncol = size(cols)

   ! Collect unique neighbors of the block atoms and sort in ascending order
   nsub = 0
   do ivec = 1, ncol
      iat = cols(ivec)
      do kat = list%inl(iat), list%inl(iat + 1) - 1
         jat = list%nlat(kat)
         if (loc(jat) /= 0) cycle
         nsub = nsub + 1
         loc(jat) = 1
         sidx(nsub) = jat
      end do
   end do
   call sort_int(sidx(:nsub))
   do isub = 1, nsub
      loc(sidx(isub)) = isub
   end do

   call system_clock(tick(2))

   xsub(1:ncol, 1:nsub) => work_x(1:ncol * nsub)
   rsub(1:ncol, 1:nsub) => work_r(1:ncol * nsub)

   call solve_sub_block(solver, alist, hlist, hfin, cols, nsub, sidx, loc, &
      & col_id, subdiag, prec, subinl, subjat, subval, xsub, rsub, work_z, &
      & work_p, work_q, niter, blkerror)

   call system_clock(tick(3))
   if (allocated(blkerror)) then
      !$omp critical (invert_list_error)
      if (.not. allocated(error)) call move_alloc(blkerror, error)
      !$omp end critical (invert_list_error)
      !$omp atomic write
      failed = .true.
   else
      do ivec = 1, ncol
         iat = cols(ivec)
         do kat = hlist%inl(iat), hfin(iat)
            ainvlist(kat) = xsub(ivec, loc(hlist%nlat(kat)))
         end do
      end do
   end if

   loc(sidx(:nsub)) = 0
   call system_clock(tick(4))
   times(:) = real(tick(2:4) - tick(1:3), wp) / real(rate, wp)

end subroutine invert_block


!> Solve the subsystem block linear system A_SS X_S = I_S in transposed (ncol, nsub) layout
subroutine solve_sub_block(self, alist, list, ufin, cols, nsub, sidx, loc, &
   & col_id, subdiag, prec, subinl, subjat, subval, xmat, res, work_z, &
   & work_p, work_q, niter, error)
   !> Conjugate-gradient solver instance
   type(cg_solver), intent(in) :: self

   !> Coefficient matrix values on the upper triangle of the neighborlist
   real(wp), intent(in) :: alist(:)

   !> Sorted list with the matrix on the upper triangle of its rows
   type(csr_list), intent(in) :: list

   !> End of the upper triangle in each row of the list
   integer(i8), intent(in) :: ufin(:)

   !> Atoms of the block columns
   integer, intent(in) :: cols(:)

   !> Number of atoms in the subsystem
   integer, intent(in) :: nsub

   !> Sorted atoms of the subsystem
   integer, intent(in) :: sidx(:)

   !> Index of the atoms of the neighborlist in the subsystem, zero outside
   integer, intent(in) :: loc(:)

   !> Column index (1..ncol) of subsystem atoms, zero on entry and exit
   integer, intent(inout) :: col_id(:)

   !> Subsystem diagonal buffer
   real(wp), intent(inout) :: subdiag(:)

   !> Diagonal preconditioner buffer
   real(wp), intent(inout) :: prec(:)

   !> Subsystem row-offset buffer
   integer, intent(inout) :: subinl(:)

   !> Subsystem column-index buffer
   integer, intent(inout) :: subjat(:)

   !> Subsystem matrix-value buffer
   real(wp), intent(inout) :: subval(:)

   !> Solution vectors in (ncol, nsub) layout
   real(wp), contiguous, intent(inout) :: xmat(:, :)

   !> Residual vectors in (ncol, nsub) layout
   real(wp), contiguous, intent(inout) :: res(:, :)

   !> Preallocated workspace for preconditioned residuals, search directions and A*P
   real(wp), target, intent(inout) :: work_z(:), work_p(:), work_q(:)

   !> Number of iterations until convergence
   integer, intent(out) :: niter

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: ncol, nrank, maxit, it, isub, ivec, info
   real(wp) :: tol_square, rel_resnorm
   real(wp) :: resnorm(block_size_def)
   real(wp), target :: dtad_buf(block_size_def * block_size_def)
   real(wp), target :: step_buf(block_size_def * block_size_def)
   real(wp), target :: updfact_buf(block_size_def * block_size_def)
   real(wp), pointer :: precres(:, :), dir(:, :), adir(:, :)
   real(wp), pointer :: dtad(:, :), step(:, :), updfact(:, :)

   niter = 0
   ncol = size(cols)
   if (ncol == 0 .or. nsub == 0) return

   tol_square = self%cgtol**2
   maxit = self%cgmiter

   precres(1:ncol, 1:nsub) => work_z(1:ncol * nsub)

   ! Extract subsystem CSR matrix and form initial guess X_0 and residual R_0 = I_S - A_SS X_0
   call get_sub_matrix_and_res(alist, list, ufin, cols, nsub, sidx, loc, &
      & col_id, subdiag, prec, subinl, subjat, subval, xmat, res)

   ! Initial squared residual norms (right-hand sides are unit vectors so ||b_i||^2 = 1)
   resnorm(:ncol) = 0.0_wp
   do isub = 1, nsub
      do ivec = 1, ncol
         resnorm(ivec) = resnorm(ivec) + res(ivec, isub)**2
      end do
      precres(:, isub) = prec(isub) * res(:, isub)
   end do
   rel_resnorm = maxval(resnorm(:ncol))

   ! Initial search directions P = orth(M^-1 R)
   call orthonormalize_rows(ncol, nsub, precres, work_p, nrank, error)
   if (allocated(error)) return
   if (rel_resnorm <= tol_square .or. nrank == 0) return

   dir(1:nrank, 1:nsub) => work_p(1:nrank * nsub)
   adir(1:nrank, 1:nsub) => work_q(1:nrank * nsub)

   do it = 1, maxit

      ! Matrix-block product Q = A_SS P
      call sub_block_matmul(nsub, nrank, subdiag, subinl, subjat, subval, &
         & dir, adir)

      dtad(1:nrank, 1:nrank) => dtad_buf(1:nrank * nrank)
      step(1:nrank, 1:ncol) => step_buf(1:nrank * ncol)
      updfact(1:nrank, 1:ncol) => updfact_buf(1:nrank * ncol)

      ! Cholesky factorization of P^T Q
      call gemm(dir, adir, dtad, transb='t')
      call potrf(dtad, info=info)
      if (info /= 0) then
         call fatal_error(error, &
            & "Block CG: P^T A P is not positive definite.")
         return
      end if

      ! Step lengths alpha = (P^T Q)^-1 (P^T R)
      call gemm(dir, res, step, transb='t')
      call potrs(dtad, step)

      ! Update solution X = X + P alpha and residual R = R - Q alpha
      call gemm(step, dir, xmat, transa='t', beta=1.0_wp)
      call gemm(step, adir, res, transa='t', alpha=-1.0_wp, beta=1.0_wp)

      ! Compute new residual norms
      resnorm(:ncol) = 0.0_wp
      do isub = 1, nsub
         do ivec = 1, ncol
            resnorm(ivec) = resnorm(ivec) + res(ivec, isub)**2
         end do
      end do
      rel_resnorm = maxval(resnorm(:ncol))

      if (rel_resnorm <= tol_square) then
         niter = it
         exit
      end if

      ! Updated preconditioned residuals Z = M^-1 R
      do isub = 1, nsub
         precres(:, isub) = prec(isub) * res(:, isub)
      end do

      ! Update factors (P^T Q)^-1 (Q^T Z), beta is their negative
      call gemm(adir, precres, updfact, transb='t')
      call potrs(dtad, updfact)

      ! Update search directions P = orth(Z + P beta)
      call gemm(updfact, dir, precres, transa='t', alpha=-1.0_wp, beta=1.0_wp)
      call orthonormalize_rows(ncol, nsub, precres, work_p, nrank, error)
      if (allocated(error)) return

      if (nrank == 0) then
         call fatal_error(error, "Block CG: search space collapsed.")
         return
      end if

      if (it == maxit) then
         call fatal_error(error, &
            & "Block CG did not converge within max iterations.")
         return
      end if

      dir(1:nrank, 1:nsub) => work_p(1:nrank * nsub)
      adir(1:nrank, 1:nsub) => work_q(1:nrank * nsub)

   end do

end subroutine solve_sub_block


!> Multiply the coefficient matrix with a block of vectors
subroutine block_matmul(nlist, vec, avec, amat, alist, list)
   !> Whether to use the compressed-row representation
   logical, intent(in) :: nlist

   !> Block of input vectors
   real(wp), intent(in) :: vec(:, :)

   !> Block of output vectors
   real(wp), intent(inout) :: avec(:, :)

   !> Dense coefficient matrix of the linear system
   real(wp), intent(in), optional :: amat(:, :)

   !> Coefficient matrix values on the upper triangle of the list
   real(wp), intent(in), optional :: alist(:)

   !> Optional sorted complete neighborlist representation of the matrix
   type(csr_list), intent(in), optional :: list

   integer :: ndim, nvec

   if (nlist) then
      ndim = size(vec, 1)
      nvec = size(vec, 2)
      call spmm_csr('N', ndim, nvec, ndim, 1.0_wp, 'SU', alist, list%nlat, &
         & list%inl(:ndim), list%inl(2:ndim + 1), vec, ndim, 0.0_wp, avec, ndim)
   else
      call gemm(amat, vec, avec)
   end if

end subroutine block_matmul


!> Extract the subsystem matrix A_SS and form initial diagonal guess X_0 and residual R_0
subroutine get_sub_matrix_and_res(alist, list, ufin, cols, nsub, sidx, loc, &
   & col_id, subdiag, prec, subinl, subjat, subval, xmat, res)
   !> Matrix values on the upper triangle of the list
   real(wp), intent(in) :: alist(:)

   !> Sorted list with the matrix on the upper triangle of its rows
   type(csr_list), intent(in) :: list

   !> End of the upper triangle in each row of the list
   integer(i8), intent(in) :: ufin(:)

   !> Atoms of the block columns
   integer, intent(in) :: cols(:)

   !> Number of atoms in the subsystem
   integer, intent(in) :: nsub

   !> Sorted atoms of the subsystem in the neighborlist
   integer, intent(in) :: sidx(:)

   !> Index of the atoms of the neighborlist in the subsystem, zero outside
   integer, intent(in) :: loc(:)

   !> Column index (1..ncol) of subsystem atoms, zero on entry and exit
   integer, intent(inout) :: col_id(:)

   !> Diagonal elements of the subsystem matrix
   real(wp), intent(inout) :: subdiag(:)

   !> Diagonal preconditioner of the subsystem matrix
   real(wp), intent(inout) :: prec(:)

   !> Row offsets of the off-diagonal upper triangle of the subsystem matrix
   integer, intent(inout) :: subinl(:)

   !> Subsystem column indices of the off-diagonal upper triangle
   integer, intent(inout) :: subjat(:)

   !> Matrix values of the off-diagonal upper triangle
   real(wp), intent(inout) :: subval(:)

   !> Initial solution guess in (ncol, nsub) layout
   real(wp), contiguous, intent(out) :: xmat(:, :)

   !> Initial residual R_0 = I_S - A_SS X_0 in (ncol, nsub) layout
   real(wp), contiguous, intent(out) :: res(:, :)

   integer :: ncol, isub, jsub, iat, ivec, jvec, nnz
   integer(i8) :: kat
   real(wp) :: val, diag_val, x0(block_size_def)

   ncol = size(cols)
   xmat(:, :nsub) = 0.0_wp
   res(:, :nsub) = 0.0_wp

   do ivec = 1, ncol
      iat = cols(ivec)
      isub = loc(iat)
      col_id(isub) = ivec
      diag_val = alist(list%inl(iat))
      x0(ivec) = 1.0_wp / (diag_val + eps)
      xmat(ivec, isub) = x0(ivec)
      res(ivec, isub) = 1.0_wp - diag_val * x0(ivec)
   end do

   nnz = 0
   subinl(1) = 1
   do isub = 1, nsub
      iat = sidx(isub)
      diag_val = alist(list%inl(iat))
      subdiag(isub) = diag_val
      prec(isub) = 1.0_wp / (diag_val + eps)
      ivec = col_id(isub)
      ! Traverse upper-triangle neighbors in ascending atom order so subjat is increasing
      do kat = ufin(iat), list%inl(iat) + 1, -1
         jsub = loc(list%nlat(kat))
         if (jsub == 0) cycle
         val = alist(kat)
         nnz = nnz + 1
         subjat(nnz) = jsub
         subval(nnz) = val
         if (ivec /= 0) res(ivec, jsub) = -val * x0(ivec)
         jvec = col_id(jsub)
         if (jvec /= 0) res(jvec, isub) = -val * x0(jvec)
      end do
      subinl(isub + 1) = nnz + 1
   end do

   do ivec = 1, ncol
      col_id(loc(cols(ivec))) = 0
   end do

end subroutine get_sub_matrix_and_res


!> Multiply the subsystem matrix with a block of vectors in (nrank, nsub) layout
pure subroutine sub_block_matmul(nsub, nrank, subdiag, subinl, subjat, subval, &
   & vec, avec)
   !> Number of atoms in the subsystem
   integer, intent(in) :: nsub

   !> Number of active vectors in the block
   integer, intent(in) :: nrank

   !> Diagonal elements of the subsystem matrix
   real(wp), intent(in) :: subdiag(:)

   !> Row offsets of the off-diagonal upper triangle of the subsystem matrix
   integer, intent(in) :: subinl(:)

   !> Subsystem column indices of the off-diagonal upper triangle
   integer, intent(in) :: subjat(:)

   !> Matrix values of the off-diagonal upper triangle
   real(wp), intent(in) :: subval(:)

   !> Block of input vectors of shape (nrank, nsub)
   real(wp), contiguous, intent(in) :: vec(:, :)

   !> Block of output vectors of shape (nrank, nsub)
   real(wp), contiguous, intent(out) :: avec(:, :)

   integer :: isub, jsub, k
   real(wp) :: val, vi(block_size_def), ai(block_size_def)

   avec(:, :nsub) = 0.0_wp
   if (nrank == block_size_def) then
      do isub = 1, nsub
         vi(:) = vec(:, isub)
         ai(:) = avec(:, isub) + subdiag(isub) * vi(:)
         do k = subinl(isub), subinl(isub + 1) - 1
            jsub = subjat(k)
            val = subval(k)
            ai(:) = ai(:) + val * vec(:, jsub)
            avec(:, jsub) = avec(:, jsub) + val * vi(:)
         end do
         avec(:, isub) = ai(:)
      end do
   else
      do isub = 1, nsub
         avec(:, isub) = avec(:, isub) + subdiag(isub) * vec(:, isub)
         do k = subinl(isub), subinl(isub + 1) - 1
            jsub = subjat(k)
            val = subval(k)
            avec(:, isub) = avec(:, isub) + val * vec(:, jsub)
            avec(:, jsub) = avec(:, jsub) + val * vec(:, isub)
         end do
      end do
   end if

end subroutine sub_block_matmul


!> Orthonormal basis of the row space of vec in (ncol, nsub) layout
subroutine orthonormalize_rows(ncol, nsub, vec, out_buf, nrank, error)
   !> Number of input vectors (rows)
   integer, intent(in) :: ncol

   !> Length of each vector (columns)
   integer, intent(in) :: nsub

   !> Input block of vectors of shape (ncol, nsub)
   real(wp), contiguous, intent(in) :: vec(:, :)

   !> Output buffer for the orthonormal basis of shape (nrank, nsub)
   real(wp), target, intent(inout) :: out_buf(:)

   !> Rank of the block of vectors
   integer, intent(out) :: nrank

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: ivec, isub, info
   real(wp) :: vnorm(block_size_def)
   real(wp), target :: gram_buf(block_size_def * block_size_def)
   real(wp), target :: eval_buf(block_size_def)
   real(wp), target :: evec_buf(block_size_def * block_size_def)
   real(wp), pointer :: gram(:, :), eval(:), evec(:, :), outvec(:, :)

   gram(1:ncol, 1:ncol) => gram_buf(1:ncol * ncol)
   eval(1:ncol) => eval_buf(1:ncol)

   ! Gram matrix and its eigendecomposition in ascending order
   call gemm(vec, vec, gram, transb='t')
   call syevd(gram, eval, info=info)
   if (info /= 0) then
      call fatal_error(error, "Block CG: eigendecomposition failed.")
      return
   end if

   ! Keep the dominant eigenvectors, largest eigenvalue first
   nrank = 0
   if (eval(ncol) > eps) then
      nrank = count(eval > orth_thr * eval(ncol))
   end if
   if (nrank == 0) return

   evec(1:nrank, 1:ncol) => evec_buf(1:nrank * ncol)
   outvec(1:nrank, 1:nsub) => out_buf(1:nrank * nsub)
   do ivec = 1, nrank
      evec(ivec, :) = gram(:, ncol - ivec + 1)
   end do

   ! New basis E^T vec, normalized row by row
   call gemm(evec, vec, outvec)
   vnorm(:nrank) = 0.0_wp
   do isub = 1, nsub
      vnorm(:nrank) = vnorm(:nrank) + outvec(:, isub)**2
   end do
   vnorm(:nrank) = 1.0_wp / sqrt(vnorm(:nrank))
   do isub = 1, nsub
      outvec(:, isub) = outvec(:, isub) * vnorm(:nrank)
   end do

end subroutine orthonormalize_rows


!> Sort an integer array in ascending order in place (Shell sort)
pure subroutine sort_int(arr)
   !> Array to sort in ascending order
   integer, intent(inout) :: arr(:)

   integer, parameter :: gaps(9) = [1750, 701, 301, 132, 57, 23, 10, 4, 1]
   integer :: n, igap, gap, i, j, val

   n = size(arr)
   do igap = 1, size(gaps)
      gap = gaps(igap)
      if (gap >= n) cycle
      do i = gap + 1, n
         val = arr(i)
         j = i
         do while (j > gap)
            if (arr(j - gap) <= val) exit
            arr(j) = arr(j - gap)
            j = j - gap
         end do
         arr(j) = val
      end do
   end do

end subroutine sort_int


!> Orthonormal basis of the column space of vec
subroutine orthonormalize(vec, outvec, nrank, error)
   !> Input block of vectors
   real(wp), contiguous, intent(in) :: vec(:, :)

   !> Output orthonormal basis in the first nrank columns
   real(wp), contiguous, intent(out) :: outvec(:, :)

   !> Rank of the block of vectors
   integer, intent(out) :: nrank

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: nvec, ivec, jvec, info
   real(wp) :: vnorm
   real(wp), allocatable :: gram(:, :), eval(:), evec(:, :)

   nvec = size(vec, 2)
   allocate(gram(nvec, nvec), eval(nvec))

   ! Gram matrix and its eigendecomposition in ascending order
   call gemm(vec, vec, gram, transa='t')
   call syevd(gram, eval, info=info)
   if (info /= 0) then
      call fatal_error(error, "Block CG: eigendecomposition failed.")
      return
   end if

   ! Keep the dominant eigenvectors, largest eigenvalue first
   nrank = 0
   if (eval(nvec) > eps) then
      nrank = count(eval > orth_thr * eval(nvec))
   end if
   if (nrank == 0) then
      outvec(:, :) = 0.0_wp
      return
   end if
   allocate(evec(nvec, nrank))
   do ivec = 1, nrank
      evec(:, ivec) = gram(:, nvec - ivec + 1)
   end do

   ! New basis vec V, normalized column by column
   call gemm(vec, evec, outvec(:, :nrank))
   do ivec = 1, nrank
      vnorm = sqrt(dot(outvec(:, ivec), outvec(:, ivec)))
      call scal(alpha=1.0_wp / vnorm, xvec=outvec(:, ivec))
   end do
   do jvec = nrank + 1, nvec
      outvec(:, jvec) = 0.0_wp
   end do

end subroutine orthonormalize


!> Print header for block CG solver
subroutine print_block_cg_header(unit, verbosity, maxit, tol, nrhs, timer)
   !> Output unit
   integer, intent(in) :: unit

   !> Verbosity level
   integer, intent(in) :: verbosity

   !> Maximum number of iterations
   integer, intent(in) :: maxit

   !> Convergence tolerance
   real(wp), intent(in) :: tol

   !> Number of right-hand sides
   integer, intent(in) :: nrhs

   !> Timer holding the accumulated initialization time
   type(timer_type), intent(in), optional :: timer

   if (verbosity > 0) then
      write(unit, '(a)') "Using Block Conjugate Gradient Solver"
      write(unit, '(a)')
      write(unit, '(a, 1x, i6)') "Max iterations : ", maxit
      write(unit, '(a, 1x, es10.2)') "Tolerance      : ", tol
      write(unit, '(a, 1x, i6)') "Right-hand sides:", nrhs
      write(unit, '(a)') "Preconditioner : Jacobi (Diagonal)"
   end if
   if (verbosity > 1) then
      write(unit, '(a, 1x, a)') "Initialisation time:", &
         & format_time(timer%get("initialization"))
      write(unit, '(a)') ''
      write(unit, '(2X,A,6X,A,8X,A,2X,A,4X,A)') &
         'iter', '|residual|', 'rank', 'relative residual', 'Time / s'
   else if (verbosity == 1) then
      write(unit, '(a)') ''
      write(unit, '(2X,A,6X,A,8X,A,2X,A)') &
         'iter', '|residual|', 'rank', 'relative residual'
   end if
end subroutine print_block_cg_header


!> Print block CG iteration progress
subroutine print_block_cg_iteration(unit, iter, res_norm, nrank, rel_resnorm, &
   & verbosity, timer)
   !> Output unit
   integer, intent(in) :: unit

   !> Current iteration number
   integer, intent(in) :: iter

   !> Largest residual norm over all right-hand sides
   real(wp), intent(in) :: res_norm

   !> Current rank of the search space
   integer, intent(in) :: nrank

   !> Largest relative residual norm over all right-hand sides
   real(wp), intent(in) :: rel_resnorm

   !> Verbosity level
   integer, intent(in) :: verbosity

   !> Timer holding the accumulated iteration time
   type(timer_type), intent(in), optional :: timer

   if (verbosity == 1) then
      write(unit, '(i6, 1x, es15.5, 1x, i11, 1x, es15.5)') &
         & iter, res_norm, nrank, rel_resnorm
   else if (verbosity > 1) then
      write(unit, '(i6, 1x, es15.5, 1x, i11, *(1x, es15.5))') &
         & iter, res_norm, nrank, rel_resnorm, timer%get("iteration")
   end if
end subroutine print_block_cg_iteration


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


!> Print final summary
subroutine print_cg_final(unit, timer, verbosity)
   !> Output unit
   integer, intent(in) :: unit

   !> Timer holding the accumulated execution time
   type(timer_type), intent(in) :: timer

   !> Verbosity level
   integer, intent(in) :: verbosity

   if (verbosity > 1) then
      write(unit, '(a, 1x, a)') "CG total time : ", format_time(timer%get("total"))
      write(unit, '(a)') ''
   end if
end subroutine print_cg_final

end module multicharge_solver_cg
