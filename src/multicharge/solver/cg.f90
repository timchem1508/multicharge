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
!> of equations, including a block variant for multiple right-hand sides and
!> the block-wise inversion of the coefficient matrix.

module multicharge_solver_cg
   use iso_fortran_env, only : output_unit
   use mctc_env, only : error_type, fatal_error, format_time, timer_type, wp, i8
   use mctc_csrlist, only : csr_list, spgemv_csr, spsymv_csr, spmm_csr
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

      !> Output verbosity
      integer, allocatable :: verbosity

      !> Whether to use a neighborlist representation
      logical, allocatable :: use_nlist

      !> Threshold for dropping elements of the inverse in the charge
      !> derivatives on a neighborlist, relative to its largest diagonal
      !> element, which makes the compressed derivatives sparser
      real(wp), allocatable :: ainvthr

      !> Buffer distance in Bohr of the subsystems for the inverse on the
      !> neighborlist beyond the rows of the block atoms, at most the cutoff of
      !> the list is effective
      real(wp), allocatable :: ainvbuf
   end type cg_input

   !> Conjugate-gradient solver with a Jacobi preconditioner
   type, extends(mchrg_solver_type) :: cg_solver
      !> Maximum number of iterations
      integer, allocatable :: cgmiter

      !> Convergence tolerance
      real(wp), allocatable :: cgtol

      !> Output verbosity
      integer, allocatable :: verbosity

      !> Whether to use a neighborlist representation
      logical, allocatable :: use_nlist

      !> Threshold for dropping elements of the inverse in the charge
      !> derivatives on a neighborlist, relative to its largest diagonal
      !> element, which makes the compressed derivatives sparser
      real(wp), allocatable :: ainvthr

      !> Buffer distance in Bohr of the subsystems for the inverse on the
      !> neighborlist beyond the rows of the block atoms, at most the cutoff of
      !> the list is effective
      real(wp), allocatable :: ainvbuf
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

   !> Default output verbosity
   integer, parameter :: verbosity_def = 0

   !> Default neighborlist usage
   logical, parameter :: use_nlist_def = .false.

   !> Default threshold for skipping elements of the inverse
   real(wp), parameter :: ainvthr_def = 0.0_wp

   !> Default buffer distance around the subsystems of the inverse
   real(wp), parameter :: ainvbuf_def = 0.0_wp

   !> Default number of columns per block of right-hand sides
   integer, parameter :: block_size_def = 16

   !> Relative eigenvalue threshold of the Gram matrix for dropping linearly
   !> dependent search directions
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
   if (allocated(input%ainvthr)) then
      self%ainvthr = input%ainvthr
   else
      self%ainvthr = ainvthr_def
   end if
   if (allocated(input%ainvbuf)) then
      self%ainvbuf = input%ainvbuf
   else
      self%ainvbuf = ainvbuf_def
   end if

end subroutine new_cg_solver


!> Solve a linear system with a diagonally preconditioned CG method
subroutine solve(self, amat, alist, xvec, vrhs, ainv, cpq, list, new_unit, error)
   !> Conjugate-gradient solver instance
   class(cg_solver), intent(in) :: self

   !> Dense coefficient matrix of the linear system
   real(wp), intent(in), optional :: amat(:, :)

   !> Coefficient matrix values in compressed-row storage
   real(wp), intent(in), optional :: alist(:)

   !> Right-hand side vector
   real(wp), intent(in) :: xvec(:)

   !> On input: initial guess; on output: solution
   real(wp), intent(inout), contiguous :: vrhs(:)

   !> Inverse matrix
   real(wp), intent(out), optional :: ainv(:, :)

   !> Whether to solve coupled-perturbed equations
   logical, intent(in), optional :: cpq

   !> Optional neighborlist representation of the matrix
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
      if (list%complete) then
         call spgemv_csr(ndim, alist, list%inl, list%nlat, vrhs, Adir)
      else
         call spsymv_csr(ndim, alist, list%inl, list%nlat, vrhs, Adir)
      end if
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
         if (list%complete) then
            call spgemv_csr(ndim, alist, list%inl, list%nlat, dir, Adir)
         else
            call spsymv_csr(ndim, alist, list%inl, list%nlat, dir, Adir)
         end if
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


!> Solve a linear system with multiple right-hand sides using the
!> breakdown-free block conjugate gradient method with a Jacobi
!> preconditioner (Ji, Sosonkina, Li, Co-HPC 2014). New search directions are
!> orthonormalized by an eigendecomposition of their Gram matrix.
subroutine solve_block(self, amat, alist, bmat, xmat, list, new_unit, error, &
   & niter)
   !> Conjugate-gradient solver instance
   class(cg_solver), intent(in) :: self

   !> Dense coefficient matrix of the linear system
   real(wp), intent(in), optional :: amat(:, :)

   !> Coefficient matrix values in compressed-row storage
   real(wp), intent(in), optional :: alist(:)

   !> Block of right-hand sides
   real(wp), intent(in) :: bmat(:, :)

   !> On input: initial guess; on output: solution
   real(wp), intent(inout) :: xmat(:, :)

   !> Optional neighborlist representation of the matrix
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
   !$omp parallel do private(ivec) &
   !$omp shared(bmat, bnorm, res, resnorm, tol_square, nrhs)
   do ivec = 1, nrhs
      bnorm(ivec) = dot(bmat(:, ivec), bmat(:, ivec))
      if (bnorm(ivec) < tol_square) bnorm(ivec) = 1.0_wp
      resnorm(ivec) = dot(res(:, ivec), res(:, ivec))
   end do
   !$omp end parallel do
   rel_resnorm = maxval(resnorm / bnorm)

   ! Initial search directions P = orth(M^-1 R)
   !$omp parallel do private(ivec) shared(dir, prec, res, nrhs)
   do ivec = 1, nrhs
      dir(:, ivec) = prec * res(:, ivec)
   end do
   !$omp end parallel do
   call orthonormalize(dir, nrank, error)
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
      !$omp parallel do private(ivec) shared(res, resnorm, nrhs)
      do ivec = 1, nrhs
         resnorm(ivec) = dot(res(:, ivec), res(:, ivec))
      end do
      !$omp end parallel do

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
      dir(:, :) = precres
      call orthonormalize(dir, nrank, error)
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


!> Invert a dense symmetric positive definite matrix with the block CG solver.
!> The columns of the inverse are solved in blocks. Since the inverse is
!> symmetric, the rows K of a block that belong to earlier columns are already
!> known and only the trailing rows U are solved from A_UU X_U = I_U - A_UK X_K,
!> which shrinks with every block.
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


!> Invert a symmetric positive definite matrix in complete compressed-row
!> storage with the block CG solver and keep only the elements on the pattern
!> of the matrix. The columns are solved in blocks of spatially close atoms,
!> each from A_SS X_S = I_S on the subsystem S of the atoms in the rows of the
!> block atoms, extended by the atoms within the buffer distance beyond the
!> radius of these rows. The inverse decays fast with the distance, without a
!> buffer the error of the subsystem is of the order of the elements outside
!> of the pattern, which are dropped. The blocks are independent and solved in
!> parallel.
subroutine invert_list(self, alist, list, ainvlist, xyz, new_unit, error)
   !> Conjugate-gradient solver instance
   class(cg_solver), intent(in) :: self

   !> Coefficient matrix values in compressed-row storage
   real(wp), intent(in) :: alist(:)

   !> Complete neighborlist representation of the matrix
   type(csr_list), intent(in) :: list

   !> Inverse of the coefficient matrix on the pattern of the list
   real(wp), intent(out) :: ainvlist(:)

   !> Cartesian coordinates for spatially compact blocks and the buffer of the
   !> subsystems, the blocks follow the atom order if absent
   real(wp), intent(in), optional :: xyz(:, :)

   !> Output unit
   integer, intent(in), optional :: new_unit

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   ! Size of the system, number of blocks and counters
   integer :: ndim, nblk, iblk, iat
   ! Atoms ordered into blocks of columns
   integer, allocatable :: order(:)
   ! Local index of the atoms in the subsystem of a block, zero outside
   integer, allocatable :: loc(:)
   ! Size of the subsystem and number of iterations of a block
   integer :: nsub, niter, maxsub, maxiter
   integer(i8) :: sumsub, sumiter
   ! Coordinates and buffer distance of the subsystems
   real(wp), allocatable :: pos(:, :)
   real(wp) :: buffer
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
         & "Inversion on a neighborlist requires a complete list.")
      return
   end if
   ! Periodic lists with translation indices hold an atom pair once per image
   if (allocated(list%nltr)) then
      call fatal_error(error, &
         & "Inversion on a neighborlist requires one entry per atom pair.")
      return
   end if
   if (size(ainvlist, kind=i8) /= size(list%nlat, kind=i8) &
      & .or. size(alist, kind=i8) < size(list%nlat, kind=i8)) then
      call fatal_error(error, "Dimension mismatch between list and ainvlist.")
      return
   end if
   ndim = size(list%inl) - 1
   if (present(xyz)) then
      if (size(xyz, 1) /= 3 .or. size(xyz, 2) < ndim) then
         call fatal_error(error, "Dimension mismatch between list and xyz.")
         return
      end if
   end if
   buffer = self%ainvbuf
   if (buffer > 0.0_wp .and. .not. present(xyz)) then
      call fatal_error(error, &
         & "Buffer of the inversion subsystems requires coordinates.")
      return
   end if
   if (ndim <= 0) return

   if (present(new_unit)) then
      unit = new_unit
   else
      unit = output_unit
   end if

   ! Blocks of spatially close atoms
   allocate(order(ndim))
   do iat = 1, ndim
      order(iat) = iat
   end do
   if (present(xyz)) then
      pos = xyz(:, :ndim)
      call get_spatial_order(pos, order, block_size_def)
   else
      allocate(pos(3, 0))
   end if
   nblk = (ndim + block_size_def - 1) / block_size_def

   call new_cg_solver(blksolver, cg_input(cgmiter=self%cgmiter, &
      & cgtol=self%cgtol, verbosity=0, use_nlist=.true.))

   failed = .false.
   sumsub = 0_i8
   sumiter = 0_i8
   maxsub = 0
   maxiter = 0
   !$omp parallel default(none) &
   !$omp shared(ndim, nblk, order, alist, list, ainvlist, pos, buffer, &
   !$omp& blksolver, error, failed) &
   !$omp private(iblk, loc, nsub, niter, skip) &
   !$omp reduction(+:sumsub, sumiter) reduction(max:maxsub, maxiter)
   allocate(loc(ndim), source=0)
   !$omp do schedule(dynamic)
   do iblk = 1, nblk
      !$omp atomic read
      skip = failed
      if (skip) cycle
      call invert_block(blksolver, alist, list, &
         & order((iblk - 1) * block_size_def + 1:min(iblk * block_size_def, ndim)), &
         & pos, buffer, loc, ainvlist, nsub, niter, error, failed)
      sumsub = sumsub + nsub
      sumiter = sumiter + niter
      maxsub = max(maxsub, nsub)
      maxiter = max(maxiter, niter)
   end do
   !$omp end do
   deallocate(loc)
   !$omp end parallel
   if (allocated(error)) return

   if (self%verbosity > 0) then
      write(unit, '(a, 1x, i0, a, i0, a)') "Inverse on neighborlist blocks   :", &
         & nblk, " of ", min(block_size_def, ndim), " columns"
      write(unit, '(a, 1x, f0.1, a, i0)') "Subsystem atoms per block        :", &
         & real(sumsub, wp) / nblk, ", max ", maxsub
      write(unit, '(a, 1x, f0.1, a, i0)') "Block CG iterations per block    :", &
         & real(sumiter, wp) / nblk, ", max ", maxiter
      write(unit, '(a)') ''
   end if

end subroutine invert_list


!> Solve a block of columns of the inverse on its subsystem of the atoms in
!> the rows of the block atoms and the neighbors of these atoms within the
!> buffer distance beyond the row radius of a block atom. The columns are
!> stored in the rows of the block atoms.
subroutine invert_block(solver, alist, list, cols, xyz, buffer, loc, ainvlist, &
   & nsub, niter, error, failed)
   !> Conjugate-gradient solver for the block
   type(cg_solver), intent(in) :: solver

   !> Coefficient matrix values in compressed-row storage
   real(wp), intent(in) :: alist(:)

   !> Complete neighborlist representation of the matrix
   type(csr_list), intent(in) :: list

   !> Atoms of the block columns
   integer, intent(in) :: cols(:)

   !> Cartesian coordinates, only used with a buffer
   real(wp), intent(in) :: xyz(:, :)

   !> Buffer distance beyond the radius of the rows of the block atoms
   real(wp), intent(in) :: buffer

   !> Local index of the atoms in the subsystem, zero on entry and exit
   integer, intent(inout) :: loc(:)

   !> Inverse on the pattern of the list, the rows of the block atoms are set
   real(wp), intent(inout) :: ainvlist(:)

   !> Number of atoms in the subsystem
   integer, intent(out) :: nsub

   !> Number of block CG iterations
   integer, intent(out) :: niter

   !> Error handling, shared by all blocks
   type(error_type), allocatable, intent(inout) :: error

   !> Whether any block failed, shared by all blocks
   logical, intent(inout) :: failed

   integer :: ncol, ncore, nrej, isub, iat, jat, ivec
   integer(i8) :: kat, pos
   real(wp) :: dist2
   logical :: inbuf
   ! Atoms of the subsystem and neighbors rejected for the buffer
   integer, allocatable :: sidx(:), rej(:)
   ! Squared radius of the rows of the block atoms extended by the buffer
   real(wp), allocatable :: rad2(:)
   ! Subsystem matrix A_SS in compressed-row storage
   real(wp), allocatable :: alsub(:)
   type(csr_list) :: lsub
   ! Right-hand sides and solutions of the subsystem
   real(wp), allocatable :: bsub(:, :), xsub(:, :)
   type(error_type), allocatable :: blkerror

   ncol = size(cols)
   allocate(sidx(size(loc)))

   ! Atoms in the rows of the block atoms, the block atoms come first
   nsub = 0
   do ivec = 1, ncol
      iat = cols(ivec)
      do kat = list%inl(iat), list%inl(iat + 1) - 1
         jat = list%nlat(kat)
         if (loc(jat) /= 0) cycle
         nsub = nsub + 1
         loc(jat) = nsub
         sidx(nsub) = jat
      end do
   end do

   ! Neighbors of these atoms within the buffer distance beyond the radius of
   ! the row of any block atom, rejected neighbors are marked as visited
   if (buffer > 0.0_wp) then
      allocate(rad2(ncol), rej(size(loc)))
      do ivec = 1, ncol
         iat = cols(ivec)
         dist2 = 0.0_wp
         do kat = list%inl(iat), list%inl(iat + 1) - 1
            dist2 = max(dist2, sum((xyz(:, list%nlat(kat)) - xyz(:, iat))**2))
         end do
         rad2(ivec) = (sqrt(dist2) + buffer)**2
      end do
      ncore = nsub
      nrej = 0
      do isub = 1, ncore
         iat = sidx(isub)
         do kat = list%inl(iat), list%inl(iat + 1) - 1
            jat = list%nlat(kat)
            if (loc(jat) /= 0) cycle
            inbuf = .false.
            do ivec = 1, ncol
               dist2 = sum((xyz(:, jat) - xyz(:, cols(ivec)))**2)
               if (dist2 <= rad2(ivec)) then
                  inbuf = .true.
                  exit
               end if
            end do
            if (inbuf) then
               nsub = nsub + 1
               loc(jat) = nsub
               sidx(nsub) = jat
            else
               nrej = nrej + 1
               loc(jat) = -1
               rej(nrej) = jat
            end if
         end do
      end do
      loc(rej(:nrej)) = 0
   end if

   ! Subsystem matrix with the elements in the list order, the diagonal stays
   ! first in each row
   lsub%complete = .true.
   allocate(lsub%inl(nsub + 1))
   lsub%inl(1) = 1
   do isub = 1, nsub
      iat = sidx(isub)
      pos = lsub%inl(isub)
      do kat = list%inl(iat), list%inl(iat + 1) - 1
         if (loc(list%nlat(kat)) /= 0) pos = pos + 1
      end do
      lsub%inl(isub + 1) = pos
   end do
   allocate(lsub%nlat(lsub%inl(nsub + 1) - 1), alsub(lsub%inl(nsub + 1) - 1))
   do isub = 1, nsub
      iat = sidx(isub)
      pos = lsub%inl(isub)
      do kat = list%inl(iat), list%inl(iat + 1) - 1
         jat = loc(list%nlat(kat))
         if (jat == 0) cycle
         lsub%nlat(pos) = jat
         alsub(pos) = alist(kat)
         pos = pos + 1
      end do
   end do

   ! Unit vectors of the block columns, initial guess from the diagonal
   allocate(bsub(nsub, ncol), xsub(nsub, ncol), source=0.0_wp)
   do ivec = 1, ncol
      isub = loc(cols(ivec))
      bsub(isub, ivec) = 1.0_wp
      xsub(isub, ivec) = 1.0_wp / (alsub(lsub%inl(isub)) + eps)
   end do

   call solver%solve_block(alist=alsub, bmat=bsub, xmat=xsub, list=lsub, &
      & error=blkerror, niter=niter)
   if (allocated(blkerror)) then
      !$omp critical (invert_list_error)
      if (.not. allocated(error)) call move_alloc(blkerror, error)
      !$omp end critical (invert_list_error)
      !$omp atomic write
      failed = .true.
   else
      ! Columns of the symmetric inverse in the rows of the block atoms
      do ivec = 1, ncol
         iat = cols(ivec)
         do kat = list%inl(iat), list%inl(iat + 1) - 1
            ainvlist(kat) = xsub(loc(list%nlat(kat)), ivec)
         end do
      end do
   end if

   loc(sidx(:nsub)) = 0

end subroutine invert_block


!> Order the atoms into spatially compact blocks by recursive coordinate
!> bisection along the largest extent, all blocks except for the last one
!> hold exactly bsize atoms
recursive subroutine get_spatial_order(xyz, order, bsize)
   !> Cartesian coordinates of all atoms
   real(wp), intent(in) :: xyz(:, :)

   !> Atoms to order, reordered on output
   integer, intent(inout) :: order(:)

   !> Number of atoms per block
   integer, intent(in) :: bsize

   integer :: nat, nleft, idir
   real(wp) :: extent(3)

   nat = size(order)
   if (nat <= bsize) return

   do idir = 1, 3
      extent(idir) = maxval(xyz(idir, order)) - minval(xyz(idir, order))
   end do
   idir = maxloc(extent, 1)
   call sort_by_key(xyz(idir, order), order)

   ! Half of the blocks to each side, the incomplete block goes to the right
   nleft = ((nat + bsize - 1) / bsize / 2) * bsize
   call get_spatial_order(xyz, order(:nleft), bsize)
   call get_spatial_order(xyz, order(nleft + 1:), bsize)

end subroutine get_spatial_order


!> Stable sort of an index array by ascending keys with a bottom-up merge sort
subroutine sort_by_key(key, idx)
   !> Keys of the entries of the index array
   real(wp), intent(in) :: key(:)

   !> Index array, sorted on output
   integer, intent(inout) :: idx(:)

   integer :: n, width, lo, mid, hi, i, j, k
   integer, allocatable :: perm(:), tmp(:)

   n = size(idx)
   allocate(perm(n), tmp(n))
   do i = 1, n
      perm(i) = i
   end do

   width = 1
   do while (width < n)
      do lo = 1, n, 2 * width
         mid = min(lo + width - 1, n)
         hi = min(lo + 2 * width - 1, n)
         i = lo
         j = mid + 1
         k = lo
         do while (i <= mid .and. j <= hi)
            if (key(perm(j)) < key(perm(i))) then
               tmp(k) = perm(j)
               j = j + 1
            else
               tmp(k) = perm(i)
               i = i + 1
            end if
            k = k + 1
         end do
         tmp(k:k + mid - i) = perm(i:mid)
         k = k + mid - i + 1
         tmp(k:k + hi - j) = perm(j:hi)
      end do
      perm(:) = tmp
      width = 2 * width
   end do

   idx(:) = idx(perm)

end subroutine sort_by_key


!> Multiply the coefficient matrix with a block of vectors, avec = A vec,
!> using either the dense matrix or its compressed-row representation
subroutine block_matmul(nlist, vec, avec, amat, alist, list)
   !> Whether to use the compressed-row representation
   logical, intent(in) :: nlist

   !> Block of input vectors
   real(wp), intent(in) :: vec(:, :)

   !> Block of output vectors
   real(wp), intent(inout) :: avec(:, :)

   !> Dense coefficient matrix of the linear system
   real(wp), intent(in), optional :: amat(:, :)

   !> Coefficient matrix values in compressed-row storage
   real(wp), intent(in), optional :: alist(:)

   !> Optional neighborlist representation of the matrix
   type(csr_list), intent(in), optional :: list

   integer :: ndim, nvec

   ndim = size(vec, 1)
   nvec = size(vec, 2)

   if (nlist) then
      if (list%complete) then
         call csr_block_matmul(alist, list, vec, avec)
      else
         call spmm_csr("N", ndim, nvec, ndim, 1.0_wp, "SU", alist, &
            & list%nlat, list%inl(1:ndim), list%inl(2:ndim+1), vec, ndim, &
            & 0.0_wp, avec, ndim)
      end if
   else
      call gemm(amat, vec, avec)
   end if

end subroutine block_matmul


!> Multiply a matrix in complete compressed-row storage with a block of
!> vectors, avec = A vec. The block is transposed first, the rows of the block
!> gathered for each matrix element are then contiguous in memory.
subroutine csr_block_matmul(alist, list, vec, avec)
   !> Matrix values in complete compressed-row storage
   real(wp), intent(in) :: alist(:)

   !> Complete neighborlist representation of the matrix
   type(csr_list), intent(in) :: list

   !> Block of input vectors
   real(wp), intent(in) :: vec(:, :)

   !> Block of output vectors
   real(wp), intent(inout) :: avec(:, :)

   integer :: ndim, nvec, iat
   integer(i8) :: kat
   real(wp), allocatable :: vect(:, :), acc(:)

   ndim = size(vec, 1)
   nvec = size(vec, 2)
   allocate(vect(nvec, ndim))

   !$omp parallel default(none) shared(ndim, nvec, alist, list, vec, avec, vect) &
   !$omp private(iat, kat, acc)
   allocate(acc(nvec))
   !$omp do schedule(static)
   do iat = 1, ndim
      vect(:, iat) = vec(iat, :)
   end do
   !$omp end do
   !$omp do schedule(guided)
   do iat = 1, ndim
      acc(:) = 0.0_wp
      do kat = list%inl(iat), list%inl(iat + 1) - 1
         acc(:) = acc + alist(kat) * vect(:, list%nlat(kat))
      end do
      avec(iat, :) = acc
   end do
   !$omp end do
   deallocate(acc)
   !$omp end parallel

end subroutine csr_block_matmul


!> Replace a block of vectors by an orthonormal basis of its column space.
!> The basis is obtained from the eigendecomposition of the Gram matrix
!> vec^T vec, eigenvectors with negligible eigenvalues are dropped and the
!> remaining columns of vec are zeroed.
subroutine orthonormalize(vec, nrank, error)
   !> On input: block of vectors; on output: orthonormal basis in the first
   !> nrank columns
   real(wp), contiguous, intent(inout) :: vec(:, :)

   !> Rank of the block of vectors
   integer, intent(out) :: nrank

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: nvec, ivec, jvec, info
   real(wp) :: vnorm
   real(wp), allocatable :: gram(:, :), eval(:), evec(:, :), tmp(:, :)

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
      vec(:, :) = 0.0_wp
      return
   end if
   allocate(evec(nvec, nrank))
   do ivec = 1, nrank
      evec(:, ivec) = gram(:, nvec - ivec + 1)
   end do

   ! New basis vec V, normalized column by column
   allocate(tmp, source=vec)
   call gemm(tmp, evec, vec(:, :nrank))
   do ivec = 1, nrank
      vnorm = sqrt(dot(vec(:, ivec), vec(:, ivec)))
      call scal(alpha=1.0_wp / vnorm, xvec=vec(:, ivec))
   end do
   do jvec = nrank + 1, nvec
      vec(:, jvec) = 0.0_wp
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
