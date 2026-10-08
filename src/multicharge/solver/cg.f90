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

   public :: cg_solver, new_cg_solver, cg_input, get_blocks

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

      !> Threshold for skipping elements of the inverse in the charge
      !> derivatives, relative to its largest diagonal element
      real(wp), allocatable :: ainvthr

      !> Whether to contract the charge derivatives with the sparse inverse
      !> instead of dense slabs, chosen from an operation count if not set
      logical, allocatable :: sparse_qgrad

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

      !> Whether to use a neighborlist representation
      logical, allocatable :: use_nlist

      !> Threshold for skipping elements of the inverse in the charge
      !> derivatives, relative to its largest diagonal element
      real(wp), allocatable :: ainvthr

      !> Whether to contract the charge derivatives with the sparse inverse
      !> instead of dense slabs, chosen from an operation count if not set
      logical, allocatable :: sparse_qgrad
   contains
      !> Solve the linear system iteratively
      procedure :: solve

      !> Solve the linear system for a block of right-hand sides
      procedure :: solve_block

      !> Invert the coefficient matrix block-wise
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
   if (allocated(input%sparse_qgrad)) self%sparse_qgrad = input%sparse_qgrad

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


!> Split the columns of a matrix into blocks for the block CG solver. Each
!> block holds at most block_size columns, the last block is zero-padded.
subroutine get_blocks(mat, box, ncol, block_size)
   !> Matrix to split, typically a nat x nat matrix of right-hand sides
   real(wp), intent(in) :: mat(:, :)

   !> Column blocks of the matrix, dimension (size(mat, 1), block_size, nblk)
   real(wp), allocatable, intent(out) :: box(:, :, :)

   !> Number of occupied columns in each block
   integer, allocatable, intent(out) :: ncol(:)

   !> Maximum number of columns per block, default is 16
   integer, intent(in), optional :: block_size

   integer :: bsize, nrow, nvec, nblk, iblk, ivec

   if (present(block_size)) then
      bsize = max(1, block_size)
   else
      bsize = block_size_def
   end if

   nrow = size(mat, 1)
   nvec = size(mat, 2)
   bsize = max(1, min(bsize, nvec))
   nblk = (nvec + bsize - 1) / bsize

   allocate(box(nrow, bsize, nblk), source=0.0_wp)
   allocate(ncol(nblk))
   do iblk = 1, nblk
      ivec = (iblk - 1) * bsize
      ncol(iblk) = min(bsize, nvec - ivec)
      box(:, :ncol(iblk), iblk) = mat(:, ivec+1:ivec+ncol(iblk))
   end do

end subroutine get_blocks


!> Solve a linear system with multiple right-hand sides using the
!> breakdown-free block conjugate gradient method with a Jacobi
!> preconditioner (Ji, Sosonkina, Li, Co-HPC 2014). New search directions are
!> orthonormalized by an eigendecomposition of their Gram matrix.
subroutine solve_block(self, amat, alist, bmat, xmat, list, new_unit, error)
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


!> Invert a symmetric positive definite matrix with the block CG solver. The
!> columns of the inverse are solved in blocks. Since the inverse is symmetric,
!> the rows K of a block that belong to earlier columns are already known and
!> only the trailing rows U are solved from A_UU X_U = I_U - A_UK X_K, which
!> shrinks with every block.
subroutine invert(self, amat, alist, ainv, list, new_unit, error)
   !> Conjugate-gradient solver instance
   class(cg_solver), intent(in) :: self

   !> Dense coefficient matrix of the linear system
   real(wp), intent(in), optional :: amat(:, :)

   !> Coefficient matrix values in compressed-row storage
   real(wp), intent(in), optional :: alist(:)

   !> Inverse of the coefficient matrix
   real(wp), intent(out) :: ainv(:, :)

   !> Optional neighborlist representation of the matrix
   type(csr_list), intent(in), optional :: list

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
   ! Trailing matrix A_UU in compressed-row storage
   real(wp), allocatable :: alsub(:)
   type(csr_list) :: lsub

   logical :: nlist

   nlist = self%use_nlist .and. .not. present(amat) .and. &
      & present(list) .and. present(alist)

   ! Dimensions check
   ndim = size(ainv, 1)
   if (size(ainv, 2) /= ndim) then
      call fatal_error(error, "Inverse matrix must be square.")
      return
   end if
   if (nlist) then
      if (size(list%inl) < ndim + 1) then
         call fatal_error(error, "Dimension mismatch between list and ainv.")
         return
      end if
   else
      if (.not. present(amat)) then
         call fatal_error(error, "No coefficient matrix provided.")
         return
      end if
      if (size(amat, 1) /= ndim .or. size(amat, 2) /= ndim) then
         call fatal_error(error, "dimension mismatch.")
         return
      end if
   end if
   if (ndim == 0) return

   allocate(diag(ndim), xknown(ndim, block_size_def))
   if (nlist) then
      do iat = 1, ndim
         diag(iat) = alist(list%inl(iat)) + eps
      end do
      lsub%complete = list%complete
      allocate(lsub%inl(ndim + 1), lsub%nlat(size(list%nlat)), &
         & alsub(size(alist)))
   else
      do iat = 1, ndim
         diag(iat) = amat(iat, iat) + eps
      end do
      allocate(abuf(int(ndim, i8) * ndim))
   end if

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
      if (nlist) then
         call get_trailing_csr(ioff, alist, list, xknown(:ioff, :ncol), alsub, &
            & lsub, bsub)
      else
         asub(1:nsub, 1:nsub) => abuf(1:int(nsub, i8) * nsub)
         asub(:, :) = amat(ioff + 1:, ioff + 1:)
         if (ioff > 0) then
            call gemm(amat(ioff + 1:, :ioff), xknown(:ioff, :ncol), bsub, &
               & alpha=-1.0_wp, beta=1.0_wp)
         end if
      end if

      ! Initial guess from the diagonal of the trailing matrix
      allocate(xsub(nsub, ncol))
      do ivec = 1, ncol
         xsub(:, ivec) = bsub(:, ivec) / diag(ioff + 1:)
      end do

      if (nlist) then
         call self%solve_block(alist=alsub, bmat=bsub, xmat=xsub, list=lsub, &
            & new_unit=new_unit, error=error)
      else
         call self%solve_block(amat=asub, bmat=bsub, xmat=xsub, &
            & new_unit=new_unit, error=error)
      end if
      if (allocated(error)) return

      ainv(ioff + 1:, ioff + 1:ioff + ncol) = xsub
      deallocate(bsub, xsub)
   end do

end subroutine invert


!> Invert a symmetric positive definite matrix in complete compressed-row
!> storage with the block CG solver and keep only the elements on the pattern
!> of the matrix. The block columns are solved from the trailing systems as in
!> invert, but the known rows X_K are taken from the rows of the block columns
!> on the pattern. This is accurate if the inverse decays faster than the
!> pattern extends, elements outside of it are dropped.
subroutine invert_list(self, alist, list, ainvlist, new_unit, error)
   !> Conjugate-gradient solver instance
   class(cg_solver), intent(in) :: self

   !> Coefficient matrix values in compressed-row storage
   real(wp), intent(in) :: alist(:)

   !> Complete neighborlist representation of the matrix
   type(csr_list), intent(in) :: list

   !> Inverse of the coefficient matrix on the pattern of the list
   real(wp), intent(out) :: ainvlist(:)

   !> Output unit
   integer, intent(in), optional :: new_unit

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   ! Size of the system, offset and number of columns of the current block
   integer :: ndim, ioff, ncol
   ! Number of trailing rows
   integer :: nsub
   ! Counters
   integer :: iat, jat, ivec
   integer(i8) :: kat
   ! Diagonal of the coefficient matrix
   real(wp), allocatable :: diag(:)
   ! Known rows X_K of the block columns
   real(wp), allocatable :: xknown(:, :)
   ! Right-hand sides and solutions of the trailing system
   real(wp), allocatable :: bsub(:, :), xsub(:, :)
   ! Trailing matrix A_UU in compressed-row storage
   real(wp), allocatable :: alsub(:)
   type(csr_list) :: lsub

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
   if (ndim <= 0) return

   allocate(diag(ndim), xknown(ndim, block_size_def))
   do iat = 1, ndim
      diag(iat) = alist(list%inl(iat)) + eps
   end do
   lsub%complete = .true.
   allocate(lsub%inl(ndim + 1), lsub%nlat(size(list%nlat)), alsub(size(alist)))
   ainvlist(:) = 0.0_wp

   do ioff = 0, ndim - 1, block_size_def
      ncol = min(block_size_def, ndim - ioff)
      nsub = ndim - ioff

      ! Rows of the block columns known from earlier columns, A^-1 = A^-T,
      ! restricted to the pattern of the rows
      xknown(:ioff, :ncol) = 0.0_wp
      do ivec = 1, ncol
         iat = ioff + ivec
         do kat = list%inl(iat), list%inl(iat + 1) - 1
            jat = list%nlat(kat)
            if (jat <= ioff) xknown(jat, ivec) = ainvlist(kat)
         end do
      end do

      ! Trailing system A_UU X_U = I_U - A_UK X_K
      allocate(bsub(nsub, ncol), source=0.0_wp)
      do ivec = 1, ncol
         bsub(ivec, ivec) = 1.0_wp
      end do
      call get_trailing_csr(ioff, alist, list, xknown(:ioff, :ncol), alsub, &
         & lsub, bsub)

      ! Initial guess from the diagonal of the trailing matrix
      allocate(xsub(nsub, ncol))
      do ivec = 1, ncol
         xsub(:, ivec) = bsub(:, ivec) / diag(ioff + 1:)
      end do

      call self%solve_block(alist=alsub, bmat=bsub, xmat=xsub, list=lsub, &
         & new_unit=new_unit, error=error)
      if (allocated(error)) return

      ! Block columns into every row holding one of them, the rows K of the
      ! block columns are the known elements
      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(ndim, ioff, ncol, list, ainvlist, xsub, xknown) &
      !$omp private(iat, jat, kat)
      do iat = 1, ndim
         do kat = list%inl(iat), list%inl(iat + 1) - 1
            jat = list%nlat(kat) - ioff
            if (jat < 1 .or. jat > ncol) cycle
            if (iat > ioff) then
               ainvlist(kat) = xsub(iat - ioff, jat)
            else
               ainvlist(kat) = xknown(iat, jat)
            end if
         end do
      end do
      deallocate(bsub, xsub)
   end do

end subroutine invert_list


!> Extract the trailing block A_UU, rows and columns ioff+1:n, of a matrix in
!> compressed-row storage and subtract A_UK X_K from the right-hand sides
subroutine get_trailing_csr(ioff, alist, list, xknown, alsub, lsub, bsub)
   !> Number of leading rows and columns K
   integer, intent(in) :: ioff

   !> Coefficient matrix values in compressed-row storage
   real(wp), intent(in) :: alist(:)

   !> Neighborlist representation of the matrix
   type(csr_list), intent(in) :: list

   !> Known rows X_K of the block of solutions
   real(wp), intent(in) :: xknown(:, :)

   !> Values of the trailing matrix A_UU
   real(wp), intent(inout) :: alsub(:)

   !> Compressed-row index of the trailing matrix A_UU
   type(csr_list), intent(inout) :: lsub

   !> Right-hand sides, I_U on input and I_U - A_UK X_K on output
   real(wp), intent(inout) :: bsub(:, :)

   integer :: nsub, ncol, isub, iat, jat, ivec
   integer(i8) :: kat, pos
   real(wp), allocatable :: rsum(:)

   nsub = size(bsub, 1)
   ncol = size(bsub, 2)

   ! Number of trailing columns in each trailing row, the diagonal stays first
   lsub%inl(1) = 1
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(nsub, ioff, list, lsub) private(isub, iat, kat)
   do isub = 1, nsub
      iat = ioff + isub
      lsub%inl(isub + 1) = 0
      do kat = list%inl(iat), list%inl(iat + 1) - 1
         if (list%nlat(kat) > ioff) lsub%inl(isub + 1) = lsub%inl(isub + 1) + 1
      end do
   end do
   do isub = 1, nsub
      lsub%inl(isub + 1) = lsub%inl(isub + 1) + lsub%inl(isub)
   end do

   ! Trailing matrix and the A_UK X_K contributions held by the rows U
   !$omp parallel default(none) &
   !$omp shared(nsub, ncol, ioff, alist, list, xknown, alsub, lsub, bsub) &
   !$omp private(isub, iat, jat, kat, pos, rsum)
   allocate(rsum(ncol))
   !$omp do schedule(runtime)
   do isub = 1, nsub
      iat = ioff + isub
      pos = lsub%inl(isub)
      rsum(:) = 0.0_wp
      do kat = list%inl(iat), list%inl(iat + 1) - 1
         jat = list%nlat(kat)
         if (jat > ioff) then
            lsub%nlat(pos) = jat - ioff
            alsub(pos) = alist(kat)
            pos = pos + 1
         else
            rsum(:) = rsum + alist(kat) * xknown(jat, :)
         end if
      end do
      bsub(isub, :) = bsub(isub, :) - rsum
   end do
   !$omp end do
   deallocate(rsum)
   !$omp end parallel

   ! In upper triangular storage A_UK is held by the rows K as A_KU
   if (.not. list%complete) then
      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(ncol, ioff, alist, list, xknown, bsub) &
      !$omp private(ivec, iat, jat, kat)
      do ivec = 1, ncol
         do jat = 1, ioff
            do kat = list%inl(jat), list%inl(jat + 1) - 1
               iat = list%nlat(kat)
               if (iat > ioff) then
                  bsub(iat - ioff, ivec) = bsub(iat - ioff, ivec) &
                     & - alist(kat) * xknown(jat, ivec)
               end if
            end do
         end do
      end do
   end if

end subroutine get_trailing_csr


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

   character(len=2) :: matdescra
   integer :: ndim, nvec

   ndim = size(vec, 1)
   nvec = size(vec, 2)

   if (nlist) then
      if (list%complete) then
         matdescra = "G "
      else
         matdescra = "SU"
      end if
      call spmm_csr("N", ndim, nvec, ndim, 1.0_wp, matdescra, alist, &
         & list%nlat, list%inl(1:ndim), list%inl(2:ndim+1), vec, ndim, &
         & 0.0_wp, avec, ndim)
   else
      call gemm(amat, vec, avec)
   end if

end subroutine block_matmul


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
