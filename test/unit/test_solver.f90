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
   use iso_fortran_env, only : output_unit
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
   use multicharge_solver_cg, only : cg_solver, new_cg_solver, cg_input
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
      & new_unittest("cg-spd-medium", test_cg_spd_medium), &
      & new_unittest("cg-spd-large", test_cg_spd_large), &
      & new_unittest("cg-ill-conditioned", test_cg_ill_conditioned), &
      & new_unittest("cg-zero-rhs", test_cg_zero_rhs), &
      & new_unittest("cg-random-spd", test_cg_random_spd) &
     ! & new_unittest("time-scaling", test_cg_spd_time_scaling), &
     ! & new_unittest("time-scaling-tri-diag-matrix", test_cg_121_time_scaling) &
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
   real(wp) :: expected(n), residual

   ! Timer variables
   real(wp) :: start_cg, end_cg, start_direct, end_direct

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
   call cpu_time(start_cg)
   call solver%solve(amat=amat, xvec=xvec, vrhs=vrhs, error=error)
   if (allocated(error)) return
   call cpu_time(end_cg)

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
   call cpu_time(start_direct)
   call solver%solve(amat=amat, xvec=xvec, vrhs=expected,  error=error)
   call cpu_time(end_direct)
   
   ! Check solution
   if (any(abs(vrhs - expected) > thr)) then
      call test_failed(error, "CG solver failed for identity matrix")
      print'(a)', "Solution:"
      print'(3es21.14)', vrhs
      print'(a)', "Expected:"
      print'(3es21.14)', expected
   else
      print '("CG Solver CPU Time : ",f6.3," seconds.")',end_cg-start_cg
      print '("Direct Solver CPU Time : ",f6.3," seconds.")',end_direct-start_direct
      print '("CG Solver ime profit : ",f6.3)', (end_direct-start_direct)/(end_cg-start_cg)
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

   ! Timer variables
   real(wp) :: start_cg, end_cg, start_direct, end_direct

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
   
   ! Expected solution (x_i = 1/diag(i))
   expected = 1.0_wp / diag
   
   ! Solve
   call cpu_time(start_cg)
   call solver%solve(amat=amat, xvec=xvec, vrhs=vrhs, error=error)
   if (allocated(error)) return
   call cpu_time(end_cg)

   ! Reference direct solver
   deallocate(solver_input)
   deallocate(solver)
   allocate(direct_input :: solver_input)
   select type(solver_input)
      type is (direct_input)
      solver_input%verbosity = verbosity
   end select
   call solver_maker(solver, solver_input,  error)

   call cpu_time(start_direct)
   expected = [0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
   call solver%solve(amat=amat, xvec=xvec, vrhs=expected, error=error)
   call cpu_time(end_direct)
   
   ! Check solution
   if (any(abs(vrhs - expected) > thr)) then
      call test_failed(error, "CG solver failed for diagonal matrix")
      print'(a)', "Solution:"
      print'(3es21.14)', vrhs
      print'(a)', "Expected:"
      print'(3es21.14)', expected
   else
      print '("CG Solver CPU Time : ",f6.3," seconds.")',end_cg-start_cg
      print '("Direct Solver CPU Time : ",f6.3," seconds.")',end_direct-start_direct
      print '("CG Solver ime profit : ",f6.3)', (end_direct-start_direct)/(end_cg-start_cg)
   end if

end subroutine test_cg_diagonal_5x5

!> Test: Small SPD matrix
subroutine test_cg_spd_small(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 3
   real(wp) :: amat(n, n), xvec(n), vrhs(n)
   real(wp) :: expected(n), residual

   ! Timer variables
   real(wp) :: start_cg, end_cg, start_direct, end_direct

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

   ! SPD matrix: A = [4 1 1; 1 3 2; 1 2 4]
   amat = reshape([4.0_wp, 1.0_wp, 1.0_wp, &
                  &1.0_wp, 3.0_wp, 2.0_wp, &
                  &1.0_wp, 2.0_wp, 4.0_wp], [n, n])
   
   ! RHS vector
   xvec = [6.0_wp, 6.0_wp, 7.0_wp]
   
   ! Initial guess
   vrhs = [0.0_wp, 0.0_wp, 0.0_wp]
   
   ! Solve
   call cpu_time(start_cg)
   call solver%solve(amat=amat, xvec=xvec, vrhs=vrhs, error=error)
   if (allocated(error)) return
   call cpu_time(end_cg)

   ! Reference direct solver
   deallocate(solver_input)
   deallocate(solver)
   allocate(direct_input :: solver_input)
   select type(solver_input)
      type is (direct_input)
      solver_input%verbosity = verbosity
   end select
   call solver_maker(solver, solver_input,  error)

   ! Expected solution (precomputed)
   expected = [1.0_wp, 1.0_wp, 1.0_wp]
   call cpu_time(start_direct)
   call solver%solve(amat=amat, xvec=xvec, vrhs=expected, error=error)
   call cpu_time(end_direct)
   
   ! Check solution
   if (any(abs(vrhs - expected) > thr)) then
      call test_failed(error, "CG solver failed for small SPD matrix")
      print'(a)', "Solution:"
      print'(3es21.14)', vrhs
      print'(a)', "Expected:"
      print'(3es21.14)', expected
   else
      print '("CG Solver CPU Time : ",f6.3," seconds.")',end_cg-start_cg
      print '("Direct Solver CPU Time : ",f6.3," seconds.")',end_direct-start_direct
      print '("CG Solver ime profit : ",f6.3)', (end_direct-start_direct)/(end_cg-start_cg)
   end if

end subroutine test_cg_spd_small

!> Test: Medium SPD matrix
subroutine test_cg_spd_medium(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 40
   real(wp) :: amat(n, n), xvec(n), vrhs(n)
   real(wp) :: expected(n), b(n), residual
   integer :: i, j

   ! Timer variables
   real(wp) :: start_cg, end_cg, start_direct, end_direct

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

   ! Create a simple SPD matrix: A = I + 0.1*E where E is matrix of ones
   amat = 0.1_wp
   do i = 1, n
      amat(i,i) = 1.0_wp
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
   call cpu_time(start_cg)
   call solver%solve(amat=amat, xvec=xvec, vrhs=vrhs, error=error)
   if (allocated(error)) return
   call cpu_time(end_cg)

   ! Reference direct solver
   deallocate(solver_input)
   deallocate(solver)
   allocate(direct_input :: solver_input)
   select type(solver_input)
      type is (direct_input)
      solver_input%verbosity = verbosity
   end select
   call solver_maker(solver, solver_input,  error)

   ! Expected solution (precomputed)
   expected = [(real(i, wp), i=1, n)]
   call cpu_time(start_direct)
   call solver%solve(amat=amat, xvec=xvec, vrhs=expected, error=error)
   call cpu_time(end_direct)
   
   ! Check solution
   if (any(abs(vrhs - expected) / max(1.0_wp, abs(expected)) > thr_rel)) then
      call test_failed(error, "CG solver failed for medium SPD matrix")
      print'(a)', "Solution:"
      print'(3es21.14)', vrhs
      print'(a)', "Expected:"
      print'(3es21.14)', expected
   else
      print '("CG Solver CPU Time : ",f6.3," seconds.")',end_cg-start_cg
      print '("Direct Solver CPU Time : ",f6.3," seconds.")',end_direct-start_direct
      print '("CG Solver ime profit : ",f6.3)', (end_direct-start_direct)/(end_cg-start_cg)
   end if

end subroutine test_cg_spd_medium

!> Test: Large SPD matrix (1000x1000)
subroutine test_cg_spd_large(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 1000
   real(wp), allocatable :: amat(:,:), xvec(:), vrhs(:)
   real(wp), allocatable :: expected(:), b(:)
   integer :: i, j

  ! Timer variables
   real(wp) :: start_cg, end_cg, start_direct, end_direct

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
   call cpu_time(start_cg)
   call solver%solve(amat=amat, xvec=xvec, vrhs=vrhs, error=error)
   if (allocated(error)) return
   call cpu_time(end_cg)

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
   call cpu_time(start_direct)
   call solver%solve(amat=amat, xvec=xvec, vrhs=expected, error=error)
   call cpu_time(end_direct)
   
   ! Check solution
   if (any(abs(vrhs - expected) / max(1.0_wp, abs(expected)) > thr_rel)) then
      call test_failed(error, "CG solver failed for medium SPD matrix")
      print'(a)', "Solution:"
      print'(3es21.14)', vrhs
      print'(a)', "Expected:"
      print'(3es21.14)', expected
   else
      print '("CG Solver CPU Time : ",f6.3," seconds.")',end_cg-start_cg
      print '("Direct Solver CPU Time : ",f6.3," seconds.")',end_direct-start_direct
      print '("CG Solver ime profit : ",f6.3)', (end_direct-start_direct)/(end_cg-start_cg)
   end if

end subroutine test_cg_spd_large

!> Test: Ill-conditioned matrix
subroutine test_cg_ill_conditioned(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 8
   real(wp) :: amat(n, n), xvec(n), vrhs(n)
   real(wp) :: expected(n), b(n), residual
   integer :: i

   ! Timer variables
   real(wp) :: start_cg, end_cg, start_direct, end_direct

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
   call cpu_time(start_cg)
   call solver%solve(amat=amat, xvec=xvec, vrhs=vrhs, error=error)
   if (allocated(error)) return
   call cpu_time(end_cg)

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
   call cpu_time(start_direct)
   call solver%solve(amat=amat, xvec=xvec, vrhs=expected, error=error)
   call cpu_time(end_direct)
   
   ! Check solution
   if (any(abs(vrhs - expected) / max(1.0_wp, abs(expected)) > thr_rel)) then
      call test_failed(error, "CG solver failed for medium SPD matrix")
      print'(a)', "Solution:"
      print'(3es21.14)', vrhs
      print'(a)', "Expected:"
      print'(3es21.14)', expected
   else
      print '("CG Solver CPU Time : ",f6.3," seconds.")',end_cg-start_cg
      print '("Direct Solver CPU Time : ",f6.3," seconds.")',end_direct-start_direct
      print '("CG Solver ime profit : ",f6.3)', (end_direct-start_direct)/(end_cg-start_cg)
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

   ! Timer variables
   real(wp) :: start_cg, end_cg, start_direct, end_direct

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
      if (i > 1) amat(i,i-1) = 0.1_wp
      if (i < n) amat(i,i+1) = 0.1_wp
   end do
   
   ! Zero RHS
   xvec = 0.0_wp
   
   ! Initial guess (non-zero to test convergence)
   vrhs = [1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp]
   
   ! Expected solution (all zeros)
   expected = 0.0_wp
   
   ! Solve
   call cpu_time(start_cg)
   call solver%solve(amat=amat, xvec=xvec, vrhs=vrhs, error=error)
   if (allocated(error)) return
   call cpu_time(end_cg)

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
   call cpu_time(start_direct)
   call solver%solve(amat=amat, xvec=xvec, vrhs=expected, error=error)
   call cpu_time(end_direct)
   
   ! Check solution
   if (any(abs(vrhs - expected) / max(1.0_wp, abs(expected)) > thr_rel)) then
      call test_failed(error, "CG solver failed for medium SPD matrix")
      print'(a)', "Solution:"
      print'(3es21.14)', vrhs
      print'(a)', "Expected:"
      print'(3es21.14)', expected
   else
      print '("CG Solver CPU Time : ",f6.3," seconds.")',end_cg-start_cg
      print '("Direct Solver CPU Time : ",f6.3," seconds.")',end_direct-start_direct
      print '("CG Solver ime profit : ",f6.3)', (end_direct-start_direct)/(end_cg-start_cg)
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
   real(wp) :: r, max_rel_error

   ! Timer variables
   real(wp) :: start_cg, end_cg, start_direct, end_direct

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
   seed = 12345  ! Fixed seed for reproducibility
   call random_seed(put=seed)

   ! Generate random matrix B
   temp = 0.0_wp
   do i = 1, n
      do j = 1, i
         call random_number(r)
         temp(i,j) = 2.0_wp * r - 1.0_wp  ! Random values in [-1, 1]
         if (i /= j) then
            temp(j,i) = temp(i,j)  ! Make symmetric
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
      amat(i,i) = amat(i,i) + real(n, wp)  ! Add diagonal dominance
   end do
   
   ! Generate random solution vector
   expected = 0.0_wp
   do i = 1, n
      call random_number(r)
      expected(i) = 2.0_wp * r - 1.0_wp  ! Random values in [-1, 1]
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
   call cpu_time(start_cg)
   call solver%solve(amat=amat, xvec=xvec, vrhs=vrhs, error=error)
   if (allocated(error)) return
   call cpu_time(end_cg)

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
   call cpu_time(start_direct)
   call solver%solve(amat=amat, xvec=xvec, vrhs=expected, error=error)
   call cpu_time(end_direct)
   
   ! Check solution
   if (any(abs(vrhs - expected) / max(1.0_wp, abs(expected)) > thr_rel)) then
      call test_failed(error, "CG solver failed for medium SPD matrix")
      print'(a)', "Solution:"
      print'(3es21.14)', vrhs
      print'(a)', "Expected:"
      print'(3es21.14)', expected
   else
      print '("CG Solver CPU Time : ",f6.3," seconds.")',end_cg-start_cg
      print '("Direct Solver CPU Time : ",f6.3," seconds.")',end_direct-start_direct
      print '("CG Solver ime profit : ",f6.3)', (end_direct-start_direct)/(end_cg-start_cg)
   end if

end subroutine test_cg_random_spd

!> Test: Test time scaling of the CG solver
subroutine test_cg_spd_time_scaling(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: n 
   integer, parameter :: max_size=13
   real(wp), allocatable :: amat(:,:), xvec(:), vrhs(:)
   real(wp), allocatable :: expected(:), b(:)
   integer :: i, j, length

   ! Timer variables
   real(wp) :: start_cg, end_cg, start_direct, end_direct

   ! Scaling rows
   real(wp) :: iter_scal(max_size), dir_scal(max_size) 

   ! Solver variables
   class(mchrg_solver_type), allocatable :: solver_test, solver_ref
   class(mchrg_solver_input), allocatable :: solver_test_input, solver_ref_input
   real(wp) :: tol = 1.0e-15_wp
   integer :: maxiter = 1000
   integer :: verbosity = 0

   allocate(cg_input :: solver_test_input)
   select type (solver_test_input)
   type is (cg_input)
      solver_test_input%cgtol = tol
      solver_test_input%cgmiter = maxiter
      solver_test_input%verbosity = verbosity
   end select
   call solver_maker(solver_test, solver_test_input,  error)

   ! Reference solution
   allocate(direct_input :: solver_ref_input)
   select type(solver_ref_input)
      type is (direct_input)
      solver_ref_input%verbosity = verbosity
   end select
   call solver_maker(solver_ref, solver_ref_input,  error)

   do length = 1, max_size
      n = 2**length
      allocate(amat(n,n), xvec(n), vrhs(n), expected(n), b(n))
      write(*,*) "Size of the vector", n

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
      call cpu_time(start_cg)
      call solver_test%solve(amat=amat, xvec=xvec, vrhs=vrhs, error=error)
      if (allocated(error)) return
      call cpu_time(end_cg)
      
      ! Reference solution
      call cpu_time(start_direct)
      call solver_ref%solve(amat=amat, xvec=xvec, vrhs=expected, error=error)
      call cpu_time(end_direct)

      ! Scaling save
      iter_scal(length) = end_cg-start_cg
      dir_scal(length) = end_direct-start_direct
      
      ! Check solution
      if (any(abs(vrhs - expected) / max(1.0_wp, abs(expected)) > thr_rel)) then
         call test_failed(error, "CG solver failed for medium SPD matrix")
         print'(a)', "Solution:"
         print'(3es21.14)', vrhs
         print'(a)', "Expected:"
         print'(3es21.14)', expected
      else
         print '("CG Solver CPU Time : ",f6.3," seconds.")',end_cg-start_cg
         print '("Direct Solver CPU Time : ",f6.3," seconds.")',end_direct-start_direct
         ! Added check to prevent division by zero on very fast runs
         if (abs(end_cg-start_cg) > epsilon(1.0_wp)) then
            print '("CG Solver time profit : ",f6.3)', (end_direct-start_direct)/(end_cg-start_cg)
         end if
      end if

      deallocate(amat, xvec, vrhs, expected, b)
   enddo

   call write_vector(iter_scal, "CG Time Scaling")
   call write_vector(dir_scal, "DIRECT Time Scaling")

end subroutine test_cg_spd_time_scaling

!> Test: Test time scaling of the CG solver for the "-1 2 -1" matrix
subroutine test_cg_121_time_scaling(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: n 
   integer, parameter :: max_size=12
   real(wp), allocatable :: amat(:,:), xvec(:), vrhs(:)
   real(wp), allocatable :: expected(:), b(:)
   integer :: i, j, length

   ! Timer variables
   real(wp) :: start_cg, end_cg, start_direct, end_direct

   ! Scaling rows
   real(wp) :: iter_scal(max_size), dir_scal(max_size) 

   ! Solver variables
   class(mchrg_solver_type), allocatable :: solver_test, solver_ref
   class(mchrg_solver_input), allocatable :: solver_test_input, solver_ref_input
   real(wp) :: tol = 1.0e-15_wp
   integer :: maxiter = 1000
   integer :: verbosity = 0

   allocate(cg_input :: solver_test_input)
   select type (solver_test_input)
   type is (cg_input)
      solver_test_input%cgtol = tol
      solver_test_input%cgmiter = maxiter
      solver_test_input%verbosity = verbosity
   end select
   call solver_maker(solver_test, solver_test_input,  error)

   ! Reference solution
   allocate(direct_input :: solver_ref_input)
   select type(solver_ref_input)
      type is (direct_input)
      solver_ref_input%verbosity = verbosity
   end select
   call solver_maker(solver_ref, solver_ref_input,  error)

   do length = 1, max_size
      n = 2**length
      allocate(amat(n,n), xvec(n), vrhs(n), expected(n), b(n))
      write(*,*) "Size of the vector", n

      ! Create a diagonally dominant (-1 2 -1) matrix
      amat = 0.0_wp
      do i = 1, n
         amat(i,i) = 2.0_wp
         if ( i > 1) then
            amat(i-1,i) = -1.0_wp
         end if
         if ( i < n) then
            amat(i+1,i) = -1.0_wp
         end if
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
      call cpu_time(start_cg)
      call solver_test%solve(amat=amat, xvec=xvec, vrhs=vrhs, error=error)
      if (allocated(error)) return
      call cpu_time(end_cg)
      
      ! Reference solution
      call cpu_time(start_direct)
      call solver_ref%solve(amat=amat, xvec=xvec, vrhs=expected, error=error)
      call cpu_time(end_direct)

      ! Scaling save
      iter_scal(length) = end_cg-start_cg
      dir_scal(length) = end_direct-start_direct
      
      ! Check solution
      if (any(abs(vrhs - expected) / max(1.0_wp, abs(expected)) > thr_rel)) then
         call test_failed(error, "CG solver failed for medium SPD matrix")
         print'(a)', "Solution:"
         print'(3es21.14)', vrhs
         print'(a)', "Expected:"
         print'(3es21.14)', expected
      else
         print '("CG Solver CPU Time : ",f6.3," seconds.")',end_cg-start_cg
         print '("Direct Solver CPU Time : ",f6.3," seconds.")',end_direct-start_direct
         ! Added check to prevent division by zero on very fast runs
         if (abs(end_cg-start_cg) > epsilon(1.0_wp)) then
            print '("CG Solver time profit : ",f6.3)', (end_direct-start_direct)/(end_cg-start_cg)
         end if
      end if

      deallocate(amat, xvec, vrhs, expected, b)
   enddo

   call write_vector(iter_scal, "CG Time Scaling")
   call write_vector(dir_scal, "DIRECT Time Scaling")

end subroutine test_cg_121_time_scaling

! Additional subroutine for vector output

subroutine write_vector(vector, name, unit)
    implicit none
    real(wp),intent(in) :: vector(:)
    character(len=*),intent(in),optional :: name
    integer, intent(in),optional :: unit
    integer :: d
    integer :: i, j, k, l, istep, iunit

    d = size(vector, dim=1)

    if (present(unit)) then
        iunit = unit
    else
        iunit = output_unit
    end if

    if (present(name)) write(iunit,'(/,"vector printed:",1x,a)') name

    do j = 1, d
        write(iunit, '(i6)', advance='no') j
        write(iunit, '(1x,f15.10)', advance='no') vector(j)
        write(iunit, '(a)')
    end do

end subroutine write_vector

end module test_solver