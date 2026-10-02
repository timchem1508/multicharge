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
!> of equations.

module multicharge_solver_cg
   use iso_fortran_env, only : output_unit
   use mctc_env, only : error_type, fatal_error, format_time, timer_type, wp, i8
   use mctc_csrlist, only : csr_list, spgemv_csr, spsymv_csr
   use multicharge_blas, only : axpy, dot, scal, symv
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

   !> Default neighborlist usage
   logical, parameter :: use_nlist_def = .false.


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


!> Separate a sparse matrix into rectangular column blocks
!>
!> The columns of an nrow x ncol matrix in compressed-row storage are cut into
!> contiguous groups of bsize columns, the last group may be narrower. Every
!> block spans all rows and stores only its non-zero elements, sorted by row,
!> in a single pass over the matrix (stable counting sort, O(nnz)). The input
!> rows are distributed over the threads, which fill disjoint segments.
subroutine get_blocks(error, bsize, nrow, ncol, ia, ja, a, bptr, brow, bcol, bval)
   !$ use omp_lib, only : omp_get_max_threads, omp_get_thread_num

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Number of columns per block
   integer, intent(in) :: bsize

   !> Number of rows of the matrix
   integer, intent(in) :: nrow

   !> Number of columns of the matrix
   integer, intent(in) :: ncol

   !> Row offsets into ja and a, size nrow + 1
   integer(i8), intent(in) :: ia(:)

   !> Column indices of the non-zero elements
   integer, intent(in) :: ja(:)

   !> Non-zero elements of the matrix
   real(wp), intent(in) :: a(:)

   !> Offset of every block into brow, bcol and bval, size number of blocks + 1
   integer(i8), allocatable, intent(out) :: bptr(:)

   !> Row index of each non-zero element, sorted by row inside a block
   integer, allocatable, intent(out) :: brow(:)

   !> Column index of each non-zero element local to its block, 1 to bsize
   integer, allocatable, intent(out) :: bcol(:)

   !> Value of each non-zero element
   real(wp), allocatable, intent(out) :: bval(:)

   integer :: nblk, nthr, tid, ilo, ihi, i, ib, it
   integer(i8) :: k, pos, total
   integer(i8), allocatable :: cnt(:, :)
   logical :: bad

   if (bsize < 1) then
      call fatal_error(error, "Block size must be positive.")
      return
   end if
   if (size(ia) /= nrow + 1) then
      call fatal_error(error, "Row offsets do not match the number of rows.")
      return
   end if

   nblk = (ncol + bsize - 1) / bsize

   nthr = 1
   !$ nthr = omp_get_max_threads()
   allocate(cnt(nblk, 0:nthr-1), source=0_i8)
   bad = .false.

   ! Count the elements of every block per thread
   !$omp parallel default(none) num_threads(nthr) &
   !$omp& shared(nrow, ncol, bsize, nthr, ia, ja, cnt) reduction(.or.: bad) &
   !$omp& private(tid, ilo, ihi, i, k, ib)
   tid = 0
   !$ tid = omp_get_thread_num()
   ilo = int((int(tid, i8) * nrow) / nthr) + 1
   ihi = int((int(tid + 1, i8) * nrow) / nthr)
   do i = ilo, ihi
      do k = ia(i), ia(i+1) - 1
         if (ja(k) < 1 .or. ja(k) > ncol) then
            bad = .true.
            cycle
         end if
         ib = (ja(k) - 1) / bsize + 1
         cnt(ib, tid) = cnt(ib, tid) + 1
      end do
   end do
   !$omp end parallel
   if (bad) then
      call fatal_error(error, "Column index out of range.")
      return
   end if

   ! Block offsets and the first write position of every thread in each block
   allocate(bptr(nblk + 1))
   total = 0
   do ib = 1, nblk
      bptr(ib) = total + 1
      do it = 0, nthr - 1
         k = cnt(ib, it)
         cnt(ib, it) = total
         total = total + k
      end do
   end do
   bptr(nblk + 1) = total + 1

   allocate(brow(total), bcol(total), bval(total))

   ! Scatter the elements, every thread keeps the row order inside a block
   !$omp parallel default(none) num_threads(nthr) &
   !$omp& shared(nrow, bsize, nthr, ia, ja, a, cnt, brow, bcol, bval) &
   !$omp& private(tid, ilo, ihi, i, k, ib, pos)
   tid = 0
   !$ tid = omp_get_thread_num()
   ilo = int((int(tid, i8) * nrow) / nthr) + 1
   ihi = int((int(tid + 1, i8) * nrow) / nthr)
   do i = ilo, ihi
      do k = ia(i), ia(i+1) - 1
         ib = (ja(k) - 1) / bsize + 1
         pos = cnt(ib, tid) + 1
         cnt(ib, tid) = pos
         brow(pos) = i
         bcol(pos) = ja(k) - (ib - 1) * bsize
         bval(pos) = a(k)
      end do
   end do
   !$omp end parallel

end subroutine get_blocks


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
