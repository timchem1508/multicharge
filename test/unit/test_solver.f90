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

module test_solver
   use mctc_env, only: wp, i8
   use mctc_env_testing, only: new_unittest, unittest_type, error_type, test_failed
   use mctc_io_structure, only: structure_type, new
   use mstore, only: get_structure
   use multicharge_model_type, only: mchrg_model_type
   use multicharge_model_eeqbc, only: eeqbc_model
   use multicharge_param, only: new_eeq2019_model, new_eeqbc2025_model
   use multicharge_charge, only: get_charges, get_eeq_charges, get_eeqbc_charges
   use multicharge_solver_type, only: mchrg_solver_type, mchrg_solver_input
   use multicharge_solver_direct, only : direct_solver, new_direct_solver, direct_input
   use multicharge_solver_cg, only : cg_solver, new_cg_solver, cg_input
   use mctc_csrlist, only : csr_list
   implicit none
   private

   public :: collect_solver

   real(wp), parameter :: thr = 100 * epsilon(1.0_wp)
   real(wp), parameter :: thr1 = 1.0e5_wp*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))
   real(wp), parameter :: thr_rel = 1.0e-6_wp


contains


!> Collect all unit tests for the CG solver
subroutine collect_solver(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
   & new_unittest("cg-identity-2x2", test_cg_identity_2x2), &
   & new_unittest("cg-diagonal-5x5", test_cg_diagonal_5x5), &
   & new_unittest("cg-spd-small", test_cg_spd_small), &
   & new_unittest("cg-spd-large", test_cg_spd_large), &
   & new_unittest("cg-ill-conditioned", test_cg_ill_conditioned), &
   & new_unittest("cg-zero-rhs", test_cg_zero_rhs), &
   & new_unittest("cg-random-spd", test_cg_random_spd), &
   & new_unittest("block-cg-dense", test_block_cg_dense), &
   & new_unittest("block-cg-sparse", test_block_cg_sparse), &
   & new_unittest("cg-invert-dense", test_cg_invert_dense), &
   & new_unittest("cg-invert-list", test_cg_invert_list), &
   & new_unittest("cg-invert-list-upper", test_cg_invert_list_upper), &
   & new_unittest("cg-invert-list-local", test_cg_invert_list_local), &
   & new_unittest("cg-invert-list-buffer", test_cg_invert_list_buffer) &
   & ]

end subroutine collect_solver

subroutine solver_maker(solver, input, error)
   !> Solver type
   class(mchrg_solver_type), intent(out), allocatable :: solver
   !> Solver input
   class(mchrg_solver_input), intent(in) :: input
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   select type (input)
   type is (cg_input)
      block
         class(cg_solver), allocatable :: tmp
         allocate(tmp)
         call new_cg_solver(tmp, input)
         call move_alloc(tmp, solver)
      end block
   type is (direct_input)
      block
         class(direct_solver), allocatable :: tmp
         allocate(tmp)
         call new_direct_solver(tmp, input)
         call move_alloc(tmp, solver)
      end block
   class default
      allocate(error)
      return
   end select

end subroutine solver_maker

!> Test: Identity matrix 2x2
subroutine test_cg_identity_2x2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 2
   real(wp) :: amat(n, n), xvec(n), vrhs(n)
   real(wp) :: expected(n)

   ! Solver variables
   class(mchrg_solver_type), allocatable :: solver
   class(mchrg_solver_input), allocatable :: solver_input
   real(wp) :: tol = 1.0e-15_wp
   integer :: maxiter = 1000
   integer :: verbosity = 0

   allocate(cg_input :: solver_input)
   select type (solver_input)
   type is (cg_input)
      solver_input%cgtol = tol
      solver_input%cgmiter = maxiter
      solver_input%verbosity = verbosity
   end select
   call solver_maker(solver, solver_input,  error)

   ! Identity matrix
   amat = 0.0_wp
   amat(1,1) = 1.0_wp
   amat(2,2) = 1.0_wp

   ! RHS vector
   xvec = [1.0_wp, 2.0_wp]

   ! Initial guess
   vrhs = [0.0_wp, 0.0_wp]

   ! Solve iteratively
   call solver%solve(amat=amat, xvec=xvec, vrhs=vrhs, error=error)
   if (allocated(error)) return

   ! Reference direct solver
   deallocate(solver_input)
   deallocate(solver)
   allocate(direct_input :: solver_input)
   select type(solver_input)
   type is (direct_input)
      solver_input%verbosity = verbosity
   end select
   call solver_maker(solver, solver_input,  error)
   expected = [0.0_wp, 0.0_wp]
   call solver%solve(amat=amat, xvec=xvec, vrhs=expected,  error=error)
   if (allocated(error)) return

   ! Check solution
   if (any(abs(vrhs - expected) > thr)) then
      call test_failed(error, "CG solver failed for identity matrix")
      print'(a)', "Solution:"
      print'(3es21.14)', vrhs
      print'(a)', "Expected:"
      print'(3es21.14)', expected
   end if

end subroutine test_cg_identity_2x2

!> Test: Diagonal matrix 5x5
subroutine test_cg_diagonal_5x5(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 5
   real(wp) :: amat(n, n), xvec(n), vrhs(n)
   real(wp) :: expected(n), diag(n)
   integer :: i

   ! Solver variables
   class(mchrg_solver_type), allocatable :: solver
   class(mchrg_solver_input), allocatable :: solver_input
   real(wp) :: tol = 1.0e-15_wp
   integer :: maxiter = 1000
   integer :: verbosity = 0

   allocate(cg_input :: solver_input)
   select type (solver_input)
   type is (cg_input)
      solver_input%cgtol = tol
      solver_input%cgmiter = maxiter
      solver_input%verbosity = verbosity
   end select
   call solver_maker(solver, solver_input,  error)

   ! Diagonal matrix with increasing values
   amat = 0.0_wp
   diag = [1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp]
   do i = 1, n
      amat(i,i) = diag(i)
   end do

   ! RHS vector
   xvec = [1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp]

   ! Initial guess
   vrhs = [0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]

   ! Solve
   call solver%solve(amat=amat, xvec=xvec, vrhs=vrhs, error=error)
   if (allocated(error)) return

   ! Reference direct solver
   deallocate(solver_input)
   deallocate(solver)
   allocate(direct_input :: solver_input)
   select type(solver_input)
   type is (direct_input)
      solver_input%verbosity = verbosity
   end select
   call solver_maker(solver, solver_input,  error)
   expected = [0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
   call solver%solve(amat=amat, xvec=xvec, vrhs=expected, error=error)
   if (allocated(error)) return

   ! Check solution
   if (any(abs(vrhs - expected) > thr)) then
      call test_failed(error, "CG solver failed for diagonal matrix")
      print'(a)', "Solution:"
      print'(3es21.14)', vrhs
      print'(a)', "Expected:"
      print'(3es21.14)', expected
   end if

end subroutine test_cg_diagonal_5x5

!> Test: Small SPD matrix
subroutine test_cg_spd_small(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 3
   real(wp) :: amat(n, n), xvec(n), vrhs(n)
   real(wp) :: expected(n)

   ! Solver variables
   class(mchrg_solver_type), allocatable :: solver
   class(mchrg_solver_input), allocatable :: solver_input
   real(wp) :: tol = 1.0e-15_wp
   integer :: maxiter = 1000
   integer :: verbosity = 0

   allocate(cg_input :: solver_input)
   select type (solver_input)
   type is (cg_input)
      solver_input%cgtol = tol
      solver_input%cgmiter = maxiter
      solver_input%verbosity = verbosity
   end select
   call solver_maker(solver, solver_input,  error)

   ! SPD matrix
   amat = reshape([4.0_wp, 1.0_wp, 1.0_wp, &
   &1.0_wp, 3.0_wp, 2.0_wp, &
   &1.0_wp, 2.0_wp, 4.0_wp], [n, n])

   ! RHS vector
   xvec = [6.0_wp, 6.0_wp, 7.0_wp]

   ! Initial guess
   vrhs = [0.0_wp, 0.0_wp, 0.0_wp]

   ! Solve
   call solver%solve(amat=amat, xvec=xvec, vrhs=vrhs, error=error)
   if (allocated(error)) return

   ! Reference direct solver
   deallocate(solver_input)
   deallocate(solver)
   allocate(direct_input :: solver_input)
   select type(solver_input)
   type is (direct_input)
      solver_input%verbosity = verbosity
   end select
   call solver_maker(solver, solver_input,  error)
   expected = [1.0_wp, 1.0_wp, 1.0_wp]
   call solver%solve(amat=amat, xvec=xvec, vrhs=expected, error=error)
   if (allocated(error)) return

   ! Check solution
   if (any(abs(vrhs - expected) > thr)) then
      call test_failed(error, "CG solver failed for small SPD matrix")
      print'(a)', "Solution:"
      print'(3es21.14)', vrhs
      print'(a)', "Expected:"
      print'(3es21.14)', expected
   end if

end subroutine test_cg_spd_small

!> Test: Large SPD matrix (1000x1000)
subroutine test_cg_spd_large(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 1000
   real(wp), allocatable :: amat(:,:), xvec(:), vrhs(:)
   real(wp), allocatable :: expected(:), b(:)
   integer :: i, j

   ! Solver variables
   class(mchrg_solver_type), allocatable :: solver
   class(mchrg_solver_input), allocatable :: solver_input
   real(wp) :: tol = 1.0e-15_wp
   integer :: maxiter = 1000
   integer :: verbosity = 0

   allocate(cg_input :: solver_input)
   select type (solver_input)
   type is (cg_input)
      solver_input%cgtol = tol
      solver_input%cgmiter = maxiter
      solver_input%verbosity = verbosity
   end select
   call solver_maker(solver, solver_input,  error)

   allocate(amat(n,n), xvec(n), vrhs(n), expected(n), b(n))

   ! Create a diagonally dominant SPD matrix
   amat = 0.0_wp
   do i = 1, n
      amat(i,i) = real(n, wp) + real(i, wp)
      do j = 1, n
         if (i /= j) then
            amat(i,j) = 1.0_wp / (abs(i-j) + 1.0_wp)
         end if
      end do
   end do

   ! Create a solution vector
   expected = [(sin(real(i, wp) * 0.1_wp), i=1, n)]

   ! Compute RHS: b = A * expected
   b = 0.0_wp
   do i = 1, n
      do j = 1, n
         b(i) = b(i) + amat(i,j) * expected(j)
      end do
   end do

   ! Initial guess
   vrhs = 0.0_wp
   xvec = b

   ! Solve
   call solver%solve(amat=amat, xvec=xvec, vrhs=vrhs, error=error)
   if (allocated(error)) return

   ! Reference direct solver
   deallocate(solver_input)
   deallocate(solver)
   allocate(direct_input :: solver_input)
   select type(solver_input)
   type is (direct_input)
      solver_input%verbosity = verbosity
   end select
   call solver_maker(solver, solver_input,  error)

   ! Reference solution
   call solver%solve(amat=amat, xvec=xvec, vrhs=expected, error=error)
   if (allocated(error)) return

   ! Check solution
   if (any(abs(vrhs - expected) / max(1.0_wp, abs(expected)) > thr_rel)) then
      call test_failed(error, "CG solver failed for large SPD matrix")
      print'(a)', "Solution:"
      print'(3es21.14)', vrhs
      print'(a)', "Expected:"
      print'(3es21.14)', expected
   end if

end subroutine test_cg_spd_large

!> Test: Ill-conditioned matrix
subroutine test_cg_ill_conditioned(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 8
   real(wp) :: amat(n, n), xvec(n), vrhs(n)
   real(wp) :: expected(n), b(n)
   integer :: i

   ! Solver variables
   class(mchrg_solver_type), allocatable :: solver
   class(mchrg_solver_input), allocatable :: solver_input
   real(wp) :: tol = 1.0e-15_wp
   integer :: maxiter = 1000
   integer :: verbosity = 0

   allocate(cg_input :: solver_input)
   select type (solver_input)
   type is (cg_input)
      solver_input%cgtol = tol
      solver_input%cgmiter = maxiter
      solver_input%verbosity = verbosity
   end select
   call solver_maker(solver, solver_input,  error)

   ! Create an ill-conditioned diagonal matrix
   amat = 0.0_wp
   do i = 1, n
      amat(i,i) = 10.0_wp ** (-(i-1))
   end do

   ! Solution vector
   expected = [(1.0_wp, i=1, n)]

   ! Compute RHS
   b = 0.0_wp
   do i = 1, n
      b(i) = amat(i,i) * expected(i)
   end do

   ! Initial guess
   vrhs = 0.0_wp
   xvec = b

   ! Solve
   call solver%solve(amat=amat, xvec=xvec, vrhs=vrhs, error=error)
   if (allocated(error)) return

   ! Reference direct solver
   deallocate(solver_input)
   deallocate(solver)
   allocate(direct_input :: solver_input)
   select type(solver_input)
   type is (direct_input)
      solver_input%verbosity = verbosity
   end select
   call solver_maker(solver, solver_input,  error)

   ! Reference solution
   call solver%solve(amat=amat, xvec=xvec, vrhs=expected, error=error)
   if (allocated(error)) return

   ! Check solution
   if (any(abs(vrhs - expected) / max(1.0_wp, abs(expected)) > thr_rel)) then
      call test_failed(error, "CG solver failed for ill-conditioned matrix")
      print'(a)', "Solution:"
      print'(3es21.14)', vrhs
      print'(a)', "Expected:"
      print'(3es21.14)', expected
   end if

end subroutine test_cg_ill_conditioned

!> Test: Zero RHS vector
subroutine test_cg_zero_rhs(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 5
   real(wp) :: amat(n, n), xvec(n), vrhs(n)
   real(wp) :: expected(n)
   integer :: i

   ! Solver variables
   class(mchrg_solver_type), allocatable :: solver
   class(mchrg_solver_input), allocatable :: solver_input
   real(wp) :: tol = 1.0e-15_wp
   integer :: maxiter = 1000
   integer :: verbosity = 0

   allocate(cg_input :: solver_input)
   select type (solver_input)
   type is (cg_input)
      solver_input%cgtol = tol
      solver_input%cgmiter = maxiter
      solver_input%verbosity = verbosity
   end select
   call solver_maker(solver, solver_input,  error)

   ! Create a simple SPD matrix
   amat = 0.0_wp
   do i = 1, n
      amat(i,i) = real(i, wp)
   end do
   do i = 1, n - 1
      amat(i,i+1) = 0.1_wp
      amat(i+1,i) = 0.1_wp
   end do

   ! Zero RHS
   xvec = 0.0_wp

   ! Initial guess (non-zero to test convergence)
   vrhs = [1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp]

   ! Expected solution (all zeros)
   expected = 0.0_wp

   ! Solve
   call solver%solve(amat=amat, xvec=xvec, vrhs=vrhs, error=error)
   if (allocated(error)) return

   ! Reference direct solver
   deallocate(solver_input)
   deallocate(solver)
   allocate(direct_input :: solver_input)
   select type(solver_input)
   type is (direct_input)
      solver_input%verbosity = verbosity
   end select
   call solver_maker(solver, solver_input,  error)

   ! Reference solution
   call solver%solve(amat=amat, xvec=xvec, vrhs=expected, error=error)
   if (allocated(error)) return

   ! Check solution
   if (any(abs(vrhs - expected) / max(1.0_wp, abs(expected)) > thr_rel)) then
      call test_failed(error, "CG solver failed for zero right-hand side")
      print'(a)', "Solution:"
      print'(3es21.14)', vrhs
      print'(a)', "Expected:"
      print'(3es21.14)', expected
   end if

end subroutine test_cg_zero_rhs

!> Test: Random SPD matrix
subroutine test_cg_random_spd(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 100
   real(wp), allocatable :: amat(:,:), xvec(:), vrhs(:)
   real(wp), allocatable :: expected(:), b(:), temp(:,:)
   integer :: i, j, k, seed_size
   integer, allocatable :: seed(:)
   real(wp) :: r

   ! Solver variables
   class(mchrg_solver_type), allocatable :: solver
   class(mchrg_solver_input), allocatable :: solver_input
   real(wp) :: tol = 1.0e-15_wp
   integer :: maxiter = 1000
   integer :: verbosity = 0

   allocate(cg_input :: solver_input)
   select type (solver_input)
   type is (cg_input)
      solver_input%cgtol = tol
      solver_input%cgmiter = maxiter
      solver_input%verbosity = verbosity
   end select
   call solver_maker(solver, solver_input,  error)

   allocate(amat(n,n), xvec(n), vrhs(n), expected(n), b(n), temp(n,n))

   ! Initialize random seed
   call random_seed(size=seed_size)
   allocate(seed(seed_size))
   seed = 12345
   call random_seed(put=seed)

   ! Generate random matrix B
   temp = 0.0_wp
   do i = 1, n
      do j = 1, i
         call random_number(r)
         ! Random values in [-1, 1]
         temp(i,j) = 2.0_wp * r - 1.0_wp
         if (i /= j) then
            ! Make symmetric
            temp(j,i) = temp(i,j)
         end if
      end do
   end do

   ! Create SPD matrix: A = B^T * B + n*I (ensures positive definiteness)
   amat = 0.0_wp
   do i = 1, n
      do j = 1, n
         do k = 1, n
            amat(i,j) = amat(i,j) + temp(k,i) * temp(k,j)
         end do
      end do
      ! Add diagonal dominance
      amat(i,i) = amat(i,i) + real(n, wp)
   end do

   ! Generate random solution vector
   expected = 0.0_wp
   do i = 1, n
      call random_number(r)
      ! Random values in [-1, 1]
      expected(i) = 2.0_wp * r - 1.0_wp
   end do

   ! Compute RHS: b = A * expected
   b = 0.0_wp
   do i = 1, n
      do j = 1, n
         b(i) = b(i) + amat(i,j) * expected(j)
      end do
   end do

   ! Initial guess
   vrhs = 0.0_wp
   xvec = b

   ! Solve
   call solver%solve(amat=amat, xvec=xvec, vrhs=vrhs, error=error)
   if (allocated(error)) return

   ! Reference direct solver
   deallocate(solver_input)
   deallocate(solver)
   allocate(direct_input :: solver_input)
   select type(solver_input)
   type is (direct_input)
      solver_input%verbosity = verbosity
   end select
   call solver_maker(solver, solver_input,  error)

   ! Reference solution
   call solver%solve(amat=amat, xvec=xvec, vrhs=expected, error=error)
   if (allocated(error)) return

   ! Check solution
   if (any(abs(vrhs - expected) / max(1.0_wp, abs(expected)) > thr_rel)) then
      call test_failed(error, "CG solver failed for random SPD matrix")
      print'(a)', "Solution:"
      print'(3es21.14)', vrhs
      print'(a)', "Expected:"
      print'(3es21.14)', expected
   end if

end subroutine test_cg_random_spd


!> Test: Block CG with a dense matrix, inverting a random SPD matrix block-wise
subroutine test_block_cg_dense(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 150, nrhs = 16
   real(wp), allocatable :: amat(:, :), ainv(:, :), unity(:, :)
   real(wp), allocatable :: xmat(:, :)
   integer :: i, ivec, ncol
   type(cg_solver) :: solver

   call new_cg_solver(solver, cg_input(cgtol=1.0e-12_wp, cgmiter=1000, &
      & verbosity=0, use_nlist=.false.))

   call random_spd(n, amat)
   allocate(unity(n, n), source=0.0_wp)
   do i = 1, n
      unity(i, i) = 1.0_wp
   end do

   allocate(ainv(n, n))
   do ivec = 0, n - 1, nrhs
      ncol = min(nrhs, n - ivec)
      allocate(xmat(n, ncol), source=0.0_wp)
      call solver%solve_block(amat=amat, bmat=unity(:, ivec+1:ivec+ncol), &
         & xmat=xmat, error=error)
      if (allocated(error)) return
      ainv(:, ivec+1:ivec+ncol) = xmat
      deallocate(xmat)
   end do

   if (any(abs(matmul(amat, ainv) - unity) > thr_rel)) then
      call test_failed(error, "Block CG failed to invert a dense SPD matrix")
      print '(a, es12.4)', "Max deviation: ", &
         & maxval(abs(matmul(amat, ainv) - unity))
   end if

end subroutine test_block_cg_dense


!> Test: Block CG with an upper-triangular CSR matrix against the dense solve
subroutine test_block_cg_sparse(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 120, nrhs = 16, band = 5
   real(wp), allocatable :: amat(:, :), alist(:), bmat(:, :)
   real(wp), allocatable :: xdense(:, :), xsparse(:, :)
   integer :: iat, jat, nnz
   type(csr_list) :: list
   type(cg_solver) :: solver

   ! Banded diagonally dominant SPD matrix
   allocate(amat(n, n), source=0.0_wp)
   do iat = 1, n
      amat(iat, iat) = 2.0_wp * band + 1.0_wp + real(iat, wp) / n
      do jat = iat + 1, min(n, iat + band)
         amat(iat, jat) = -1.0_wp / real(jat - iat, wp)
         amat(jat, iat) = amat(iat, jat)
      end do
   end do

   ! Upper triangle in CSR format, diagonal element first in each row
   list%complete = .false.
   allocate(list%inl(n + 1), list%nlat(n * (band + 1)), alist(n * (band + 1)))
   nnz = 0
   do iat = 1, n
      list%inl(iat) = nnz + 1
      do jat = iat, min(n, iat + band)
         nnz = nnz + 1
         list%nlat(nnz) = jat
         alist(nnz) = amat(iat, jat)
      end do
   end do
   list%inl(n + 1) = nnz + 1

   allocate(bmat(n, nrhs))
   call random_number(bmat)
   allocate(xdense(n, nrhs), xsparse(n, nrhs), source=0.0_wp)

   call new_cg_solver(solver, cg_input(cgtol=1.0e-12_wp, cgmiter=1000, &
      & verbosity=0, use_nlist=.false.))
   call solver%solve_block(amat=amat, bmat=bmat, xmat=xdense, error=error)
   if (allocated(error)) return

   call new_cg_solver(solver, cg_input(cgtol=1.0e-12_wp, cgmiter=1000, &
      & verbosity=0, use_nlist=.true.))
   call solver%solve_block(alist=alist(:nnz), bmat=bmat, xmat=xsparse, &
      & list=list, error=error)
   if (allocated(error)) return

   if (any(abs(matmul(amat, xsparse) - bmat) > thr_rel)) then
      call test_failed(error, "Block CG failed for a sparse SPD matrix")
      return
   end if
   if (any(abs(xsparse - xdense) > thr_rel)) then
      call test_failed(error, "Sparse and dense block CG solutions differ")
   end if

end subroutine test_block_cg_sparse


!> Test: Block-CG inversion of a dense random SPD matrix
subroutine test_cg_invert_dense(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 150
   real(wp), allocatable :: amat(:, :), ainv(:, :)
   type(cg_solver) :: solver

   call new_cg_solver(solver, cg_input(cgtol=1.0e-14_wp, cgmiter=1000, &
      & verbosity=0, use_nlist=.false.))

   call random_spd(n, amat)
   allocate(ainv(n, n))
   call solver%invert(amat=amat, ainv=ainv, error=error)
   if (allocated(error)) return

   call check_inverse(error, amat, ainv)

end subroutine test_cg_invert_dense


!> Test: Block-CG inversion on the pattern of a complete list, for a matrix
!> of interleaved atom groups whose inverse has the same pattern
subroutine test_cg_invert_list(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 90, ngroup = 7
   real(wp), allocatable :: amat(:, :), alist(:), adense(:, :), ainvlist(:)
   type(csr_list) :: list
   type(cg_solver) :: solver
   integer :: iat
   integer(i8) :: kat
   real(wp) :: dev

   call grouped_spd(n, ngroup, amat)
   call dense_to_csr(amat, .true., alist, list)
   allocate(adense(n, n), ainvlist(size(list%nlat)))

   call new_cg_solver(solver, cg_input(cgtol=1.0e-14_wp, cgmiter=1000, &
      & verbosity=0, use_nlist=.false.))
   call solver%invert(amat=amat, ainv=adense, error=error)
   if (allocated(error)) return

   call new_cg_solver(solver, cg_input(cgtol=1.0e-14_wp, cgmiter=1000, &
      & verbosity=0, use_nlist=.true.))
   call solver%invert_list(alist, list, ainvlist, error=error)
   if (allocated(error)) return

   dev = 0.0_wp
   do iat = 1, n
      do kat = list%inl(iat), list%inl(iat + 1) - 1
         dev = max(dev, abs(ainvlist(kat) - adense(list%nlat(kat), iat)))
      end do
   end do
   if (dev > thr1) then
      call test_failed(error, "Inverse on the list differs from the dense inverse")
      print '(a, es12.4)', "Max deviation: ", dev
   end if

end subroutine test_cg_invert_list


!> Test: Block-CG inversion on an upper-triangular list reports an error
subroutine test_cg_invert_list_upper(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 20, band = 3
   real(wp), allocatable :: amat(:, :), alist(:), ainvlist(:)
   type(csr_list) :: list
   type(cg_solver) :: solver

   call banded_spd(n, band, amat)
   call dense_to_csr(amat, .false., alist, list)
   allocate(ainvlist(size(list%nlat)))

   call new_cg_solver(solver, cg_input(cgtol=1.0e-14_wp, cgmiter=1000, &
      & verbosity=0, use_nlist=.true.))
   call solver%invert_list(alist, list, ainvlist, error=error)
   if (.not. allocated(error)) then
      call test_failed(error, "Upper-triangular list was not reported")
      return
   end if
   deallocate(error)

end subroutine test_cg_invert_list_upper


!> Test: Block-CG inversion on the subsystems of spatial blocks for a chain
!> whose atom order is shuffled, the buffer reduces the subsystem error
subroutine test_cg_invert_list_local(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 300
   real(wp), parameter :: cutoff = 3.0_wp
   ! Deviation of the subsystems without buffer, the elements of the inverse
   ! decay to 1e-5 within the cutoff
   real(wp), parameter :: thr_sub = 1.0e-4_wp
   real(wp), allocatable :: xyz(:, :), amat(:, :), alist(:), adense(:, :), &
      & ainvlist(:)
   type(csr_list) :: list
   type(cg_solver) :: solver
   real(wp) :: dev0, devb

   call chain_spd(n, cutoff, xyz, amat)
   call dense_to_csr(amat, .true., alist, list)
   allocate(adense(n, n), ainvlist(size(list%nlat)))

   call new_cg_solver(solver, cg_input(cgtol=1.0e-14_wp, cgmiter=1000, &
      & verbosity=0, use_nlist=.false.))
   call solver%invert(amat=amat, ainv=adense, error=error)
   if (allocated(error)) return

   ! Subsystems of the rows of the block atoms
   call new_cg_solver(solver, cg_input(cgtol=1.0e-14_wp, cgmiter=1000, &
      & verbosity=0, use_nlist=.true.))
   call solver%invert_list(alist, list, ainvlist, xyz=xyz, error=error)
   if (allocated(error)) return
   dev0 = list_deviation(list, ainvlist, adense)

   ! Subsystems extended by a buffer of the cutoff
   call new_cg_solver(solver, cg_input(cgtol=1.0e-14_wp, cgmiter=1000, &
      & verbosity=0, use_nlist=.true., ainvbuf=cutoff))
   call solver%invert_list(alist, list, ainvlist, xyz=xyz, error=error)
   if (allocated(error)) return
   devb = list_deviation(list, ainvlist, adense)

   if (dev0 > thr_sub) then
      call test_failed(error, "Inverse on the subsystems differs from the dense inverse")
      print '(a, es12.4)', "Max deviation: ", dev0
      return
   end if
   if (devb > thr2 .or. devb > 1.0e-2_wp * dev0) then
      call test_failed(error, "Buffer of the subsystems does not reduce the deviation")
      print '(a, 2es12.4)', "Max deviation: ", dev0, devb
   end if

end subroutine test_cg_invert_list_local


!> Test: Block-CG inversion with a buffer of the subsystems reports missing
!> coordinates
subroutine test_cg_invert_list_buffer(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 30
   real(wp), allocatable :: xyz(:, :), amat(:, :), alist(:), ainvlist(:)
   type(csr_list) :: list
   type(cg_solver) :: solver

   call chain_spd(n, 3.0_wp, xyz, amat)
   call dense_to_csr(amat, .true., alist, list)
   allocate(ainvlist(size(list%nlat)))

   call new_cg_solver(solver, cg_input(cgtol=1.0e-14_wp, cgmiter=1000, &
      & verbosity=0, use_nlist=.true., ainvbuf=3.0_wp))
   call solver%invert_list(alist, list, ainvlist, error=error)
   if (.not. allocated(error)) then
      call test_failed(error, "Buffer without coordinates was not reported")
      return
   end if
   deallocate(error)

end subroutine test_cg_invert_list_buffer


!> Largest deviation of the inverse on the pattern of a list from the dense
!> inverse, relative to its largest element
function list_deviation(list, ainvlist, adense) result(dev)

   !> Complete neighborlist
   type(csr_list), intent(in) :: list

   !> Inverse on the pattern of the list
   real(wp), intent(in) :: ainvlist(:)

   !> Dense inverse
   real(wp), intent(in) :: adense(:, :)

   !> Relative deviation
   real(wp) :: dev

   integer :: iat
   integer(i8) :: kat

   dev = 0.0_wp
   do iat = 1, size(adense, 1)
      do kat = list%inl(iat), list%inl(iat + 1) - 1
         dev = max(dev, abs(ainvlist(kat) - adense(list%nlat(kat), iat)))
      end do
   end do
   dev = dev / maxval(abs(adense))

end function list_deviation


!> Check that a matrix is the symmetric inverse of an SPD matrix
subroutine check_inverse(error, amat, ainv)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> SPD matrix
   real(wp), intent(in) :: amat(:, :)

   !> Approximate inverse of the matrix
   real(wp), intent(in) :: ainv(:, :)

   real(wp), allocatable :: unity(:, :)
   integer :: i

   allocate(unity(size(amat, 1), size(amat, 1)), source=0.0_wp)
   do i = 1, size(amat, 1)
      unity(i, i) = 1.0_wp
   end do

   if (any(abs(matmul(amat, ainv) - unity) > thr1)) then
      call test_failed(error, "Block-CG inverse does not invert the matrix")
      print '(a, es12.4)', "Max deviation: ", &
         & maxval(abs(matmul(amat, ainv) - unity))
      return
   end if
   if (any(abs(ainv - transpose(ainv)) > thr1)) then
      call test_failed(error, "Block-CG inverse is not symmetric")
      print '(a, es12.4)', "Max asymmetry: ", maxval(abs(ainv - transpose(ainv)))
   end if

end subroutine check_inverse


!> Generate a banded, diagonally dominant SPD matrix
subroutine banded_spd(n, band, amat)

   !> Dimension of the matrix
   integer, intent(in) :: n

   !> Number of off-diagonals on each side
   integer, intent(in) :: band

   !> Banded SPD matrix
   real(wp), allocatable, intent(out) :: amat(:, :)

   integer :: iat, jat

   allocate(amat(n, n), source=0.0_wp)
   do iat = 1, n
      amat(iat, iat) = 2.0_wp * band + 1.0_wp + real(iat, wp) / n
      do jat = iat + 1, min(n, iat + band)
         amat(iat, jat) = -1.0_wp / real(jat - iat, wp)
         amat(jat, iat) = amat(iat, jat)
      end do
   end do

end subroutine banded_spd


!> Generate an SPD matrix for a chain of atoms with unit spacing, coupling
!> atoms within the cutoff by exp(-r). The atom order is shuffled along the
!> chain.
subroutine chain_spd(n, cutoff, xyz, amat)

   !> Number of atoms
   integer, intent(in) :: n

   !> Cutoff distance of the coupling
   real(wp), intent(in) :: cutoff

   !> Cartesian coordinates of the atoms
   real(wp), allocatable, intent(out) :: xyz(:, :)

   !> Chain SPD matrix
   real(wp), allocatable, intent(out) :: amat(:, :)

   integer :: iat, jat
   real(wp) :: dist

   allocate(xyz(3, n), source=0.0_wp)
   do iat = 1, n
      xyz(1, iat) = real(mod(7 * (iat - 1), n), wp)
   end do

   allocate(amat(n, n), source=0.0_wp)
   do iat = 1, n
      do jat = 1, n
         dist = abs(xyz(1, iat) - xyz(1, jat))
         if (jat == iat .or. dist > cutoff) cycle
         amat(jat, iat) = -exp(-dist)
      end do
      amat(iat, iat) = sum(abs(amat(:, iat))) + 10.0_wp
   end do

end subroutine chain_spd


!> Generate an SPD matrix coupling only atoms of the same group, the groups
!> interleave as the remainder of the atom index
subroutine grouped_spd(n, ngroup, amat)

   !> Dimension of the matrix
   integer, intent(in) :: n

   !> Number of groups
   integer, intent(in) :: ngroup

   !> Grouped SPD matrix
   real(wp), allocatable, intent(out) :: amat(:, :)

   integer :: iat, jat

   allocate(amat(n, n), source=0.0_wp)
   do iat = 1, n
      do jat = 1, n
         if (jat == iat .or. mod(iat, ngroup) /= mod(jat, ngroup)) cycle
         amat(jat, iat) = -1.0_wp / real(1 + abs(iat - jat), wp)
      end do
      amat(iat, iat) = sum(abs(amat(:, iat))) + 1.0_wp + real(iat, wp) / n
   end do

end subroutine grouped_spd


!> Compressed-row storage of the nonzero elements of a symmetric matrix with
!> the diagonal element first in each row
subroutine dense_to_csr(amat, complete, alist, list)

   !> Dense symmetric matrix
   real(wp), intent(in) :: amat(:, :)

   !> Store the complete matrix instead of its upper triangle
   logical, intent(in) :: complete

   !> Matrix values in compressed-row storage
   real(wp), allocatable, intent(out) :: alist(:)

   !> Compressed-row index of the matrix
   type(csr_list), intent(out) :: list

   integer :: n, iat, jat, jmin, nnz

   n = size(amat, 1)
   nnz = count(abs(amat) > 0.0_wp)
   list%complete = complete
   allocate(list%inl(n + 1), list%nlat(nnz), alist(nnz))
   nnz = 0
   do iat = 1, n
      list%inl(iat) = nnz + 1
      nnz = nnz + 1
      list%nlat(nnz) = iat
      alist(nnz) = amat(iat, iat)
      jmin = 1
      if (.not. complete) jmin = iat + 1
      do jat = jmin, n
         if (jat == iat .or. .not. abs(amat(jat, iat)) > 0.0_wp) cycle
         nnz = nnz + 1
         list%nlat(nnz) = jat
         alist(nnz) = amat(jat, iat)
      end do
   end do
   list%inl(n + 1) = nnz + 1

end subroutine dense_to_csr


!> Generate a well-conditioned random SPD matrix
subroutine random_spd(n, amat)

   !> Dimension of the matrix
   integer, intent(in) :: n

   !> Random SPD matrix
   real(wp), allocatable, intent(out) :: amat(:, :)

   integer :: i

   allocate(amat(n, n))
   call random_number(amat)
   amat = 0.5_wp * (amat + transpose(amat))
   do i = 1, n
      amat(i, i) = amat(i, i) + real(n, wp)
   end do

end subroutine random_spd

end module test_solver
