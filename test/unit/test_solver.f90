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
   use mctc_env, only: wp
   use mctc_env_testing, only: new_unittest, unittest_type, error_type, test_failed
   use mctc_io_structure, only: structure_type, new
   use mstore, only: get_structure
   use multicharge_model_type, only: mchrg_model_type
   use multicharge_model_eeqbc, only: eeqbc_model
   use multicharge_param, only: new_eeq2019_model, new_eeqbc2025_model
   use multicharge_charge, only: get_charges, get_eeq_charges, get_eeqbc_charges
   use multicharge_solver_type, only: mchrg_solver_type, mchrg_solver_input
   use multicharge_solver_direct, only : direct_solver, new_direct_solver, direct_input
   use multicharge_solver_cg, only : cg_solver, new_cg_solver, cg_input, get_blocks
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
   & new_unittest("block-cg-get-blocks", test_get_blocks), &
   & new_unittest("block-cg-dense", test_block_cg_dense), &
   & new_unittest("block-cg-sparse", test_block_cg_sparse) &
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


!> Test: Splitting a square matrix into column blocks
subroutine test_get_blocks(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 150
   real(wp), allocatable :: mat(:, :)
   real(wp), allocatable :: blocks(:, :, :)
   integer, allocatable :: ncol(:)
   integer :: iblk, ivec

   allocate(mat(n, n))
   call random_number(mat)
   call get_blocks(mat, blocks, ncol)

   if (any(shape(blocks) /= [n, 16, 10])) then
      call test_failed(error, "Wrong shape of the column blocks")
      return
   end if
   if (any(ncol /= [spread(16, 1, 9), 6])) then
      call test_failed(error, "Wrong number of columns per block")
      return
   end if
   do iblk = 1, size(ncol)
      ivec = (iblk - 1) * 16
      if (any(abs(blocks(:, :ncol(iblk), iblk) &
         & - mat(:, ivec+1:ivec+ncol(iblk))) > thr)) then
         call test_failed(error, "Column blocks do not match the matrix")
         return
      end if
   end do
   iblk = size(ncol)
   if (any(abs(blocks(:, ncol(iblk)+1:, iblk)) > thr)) then
      call test_failed(error, "Last column block is not zero-padded")
   end if

end subroutine test_get_blocks


!> Test: Block CG with a dense matrix, inverting a random SPD matrix block-wise
subroutine test_block_cg_dense(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 150
   real(wp), allocatable :: amat(:, :), ainv(:, :), unity(:, :)
   real(wp), allocatable :: blocks(:, :, :), xmat(:, :)
   integer, allocatable :: ncol(:)
   integer :: i, iblk, ivec
   type(cg_solver) :: solver

   call new_cg_solver(solver, cg_input(cgtol=1.0e-12_wp, cgmiter=1000, &
      & verbosity=0, use_nlist=.false.))

   call random_spd(n, amat)
   allocate(unity(n, n), source=0.0_wp)
   do i = 1, n
      unity(i, i) = 1.0_wp
   end do

   call get_blocks(unity, blocks, ncol)
   allocate(ainv(n, n))
   do iblk = 1, size(ncol)
      ivec = (iblk - 1) * size(blocks, 2)
      allocate(xmat(n, ncol(iblk)), source=0.0_wp)
      call solver%solve_block(amat=amat, bmat=blocks(:, :ncol(iblk), iblk), &
         & xmat=xmat, error=error)
      if (allocated(error)) return
      ainv(:, ivec+1:ivec+ncol(iblk)) = xmat
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
