# 1 "/home/zakharov/Documents/prog/mcharge/multicharge/test/unit/test_solver.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/home/zakharov/Documents/prog/mcharge/multicharge/cmake_build//"
# 1 "/home/zakharov/Documents/prog/mcharge/multicharge/test/unit/test_solver.f90"
module test_solver
   use mctc_env, only: wp
   use mctc_env_testing, only: new_unittest, unittest_type, error_type, test_failed
   use print_matrix, only: write_vector, write_matrix
   use mctc_io_structure, only: structure_type, new
   use mstore, only: get_structure
   use multicharge_model_type, only: mchrg_model_type
   use multicharge_model_eeqbc, only: eeqbc_model
   use multicharge_param, only: new_eeq2019_model, new_eeqbc2025_model
   use multicharge_model_cache, only: cache_container
   use multicharge_charge, only: get_charges, get_eeq_charges, get_eeqbc_charges
   use multicharge_solver_type, only: mchrg_solver_type, mchrg_solver_input
   use multicharge_solver_direct, only : mchrg_solver_direct, new_direct_solver, direct_input
   use multicharge_solver_cg, only : mchrg_solver_cg, new_cg_solver, cg_input
   use multicharge_solver, only: new_mchrg_solver, mchrg_solver_type, mchrg_solver_direct, &
      & mchrg_solver_cg, mchrg_solver_input, cg_input, direct_input
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
      & new_unittest("cg-preconditioned", test_cg_preconditioned) &
      !& new_unittest("time-scaling", test_cg_spd_time_scaling), &
      !& new_unittest("time-scaling-(-12-1)-matrix", test_cg_121_time_scaling) &
      & ]

end subroutine collect_solver

!> Test 1: Identity matrix 2x2
subroutine test_cg_identity_2x2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 2
   logical, parameter :: cpq = .false.
   real(wp) :: amat(n, n), xvec(n), vrhs(n), ainv(n, n)
   real(wp) :: expected(n), residual

   ! Timer variables
   real(wp) :: start_cg, end_cg, start_direct, end_direct


   ! Solver variables
   real(wp), allocatable :: tol
   integer, allocatable :: maxiter 
   class(mchrg_solver_type), allocatable :: solver
   class(mchrg_solver_input), allocatable :: solver_input

   allocate(cg_input :: solver_input)
   select type(solver_input)
   type is (cg_input)
      if (allocated(maxiter)) then
         solver_input%cgmiter = maxiter
      end if
      if (allocated(tol)) then
         solver_input%cgtol = tol
      end if
         block
             class(mchrg_solver_cg), allocatable :: tmp
             allocate(tmp)
             call new_cg_solver(tmp, solver_input)
             call move_alloc(tmp, solver)
         end block
   end select
   
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
   call solver%solve(amat, xvec, vrhs, ainv, cpq, error=error)
   if (allocated(error)) return
   call cpu_time(end_cg)

   ! Reference solution
   deallocate(solver_input)
   deallocate(solver)
   allocate(direct_input :: solver_input)
   allocate(mchrg_solver_direct :: solver)
   select type (solver_input)
   type is (direct_input)
      call new_direct_solver(solver, solver_input)
   end select

   expected = [0.0_wp, 0.0_wp]
   call cpu_time(start_direct)
   call solver%solve(amat, xvec, expected, ainv, cpq, error=error)
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

!> Test 2: Diagonal matrix 5x5
subroutine test_cg_diagonal_5x5(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 5
   logical, parameter :: cpq = .false.
   real(wp) :: amat(n, n), xvec(n), vrhs(n), ainv(n, n)
   real(wp) :: expected(n), diag(n)
   integer :: i

   ! Timer variables
   real(wp) :: start_cg, end_cg, start_direct, end_direct

   !> Solver variables
   real(wp), allocatable :: tol
   integer, allocatable :: maxiter 
   class(mchrg_solver_type), allocatable :: solver
   class(mchrg_solver_input), allocatable :: solver_input

   allocate(cg_input :: solver_input)
   select type(solver_input)
   type is (cg_input)
      if (allocated(maxiter)) then
         solver_input%cgmiter = maxiter
      end if
      if (allocated(tol)) then
         solver_input%cgtol = tol
      end if
         block
             class(mchrg_solver_cg), allocatable :: tmp
             allocate(tmp)
             call new_cg_solver(tmp, solver_input)
             call move_alloc(tmp, solver)
         end block
   end select

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
   call solver%solve(amat, xvec, vrhs, ainv, cpq, error=error)
   if (allocated(error)) return
   call cpu_time(end_cg)

   ! Reference solution
   deallocate(solver_input)
   deallocate(solver)
   allocate(direct_input :: solver_input)
   allocate(mchrg_solver_direct :: solver)
   select type (solver_input)
   type is (direct_input)
      call new_direct_solver(solver, solver_input)
   end select

   call cpu_time(start_direct)
   expected = [0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
   call solver%solve(amat, xvec, expected, ainv, cpq, error=error)
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

!> Test 3: Small SPD matrix
subroutine test_cg_spd_small(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 3
   logical, parameter :: cpq = .false.
   real(wp) :: amat(n, n), xvec(n), vrhs(n), ainv(n, n)
   real(wp) :: expected(n), residual

   ! Timer variables
   real(wp) :: start_cg, end_cg, start_direct, end_direct

   !> Solver variables
   real(wp), allocatable :: tol
   integer, allocatable :: maxiter 
   class(mchrg_solver_type), allocatable :: solver
   class(mchrg_solver_input), allocatable :: solver_input

   allocate(cg_input :: solver_input)
   select type(solver_input)
   type is (cg_input)
      if (allocated(maxiter)) then
         solver_input%cgmiter = maxiter
      end if
      if (allocated(tol)) then
         solver_input%cgtol = tol
      end if
         block
             class(mchrg_solver_cg), allocatable :: tmp
             allocate(tmp)
             call new_cg_solver(tmp, solver_input)
             call move_alloc(tmp, solver)
         end block
   end select

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
   call solver%solve(amat, xvec, vrhs, ainv, cpq, error=error)
   if (allocated(error)) return
   call cpu_time(end_cg)

   ! Reference solution
   deallocate(solver_input)
   deallocate(solver)
   allocate(direct_input :: solver_input)
   allocate(mchrg_solver_direct :: solver)
   select type (solver_input)
   type is (direct_input)
      call new_direct_solver(solver, solver_input)
   end select

   ! Expected solution (precomputed)
   expected = [1.0_wp, 1.0_wp, 1.0_wp]
   call cpu_time(start_direct)
   call solver%solve(amat, xvec, expected, ainv, cpq, error=error)
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

!> Test 4: Large SPD matrix (1000x1000)
subroutine test_cg_spd_large(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 1000
   logical, parameter :: cpq = .false.
   real(wp), allocatable :: amat(:,:), xvec(:), vrhs(:), ainv(:,:)
   real(wp), allocatable :: expected(:), b(:)
   integer :: i, j

  ! Timer variables
   real(wp) :: start_cg, end_cg, start_direct, end_direct

   !> Solver variables
   real(wp), allocatable :: tol
   integer, allocatable :: maxiter 
   class(mchrg_solver_type), allocatable :: solver
   class(mchrg_solver_input), allocatable :: solver_input

   allocate(cg_input :: solver_input)
   select type(solver_input)
   type is (cg_input)
      if (allocated(maxiter)) then
         solver_input%cgmiter = maxiter
      end if
      if (allocated(tol)) then
         solver_input%cgtol = tol
      end if
         block
             class(mchrg_solver_cg), allocatable :: tmp
             allocate(tmp)
             call new_cg_solver(tmp, solver_input)
             call move_alloc(tmp, solver)
         end block
   end select

   allocate(amat(n,n), xvec(n), vrhs(n), ainv(n,n), expected(n), b(n))

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
   call solver%solve(amat, xvec, vrhs, ainv, cpq, error=error)
   if (allocated(error)) return
   call cpu_time(end_cg)

   ! Reference solution
   deallocate(solver_input)
   deallocate(solver)
   allocate(direct_input :: solver_input)
   allocate(mchrg_solver_direct :: solver)
   select type (solver_input)
   type is (direct_input)
      call new_direct_solver(solver, solver_input)
   end select

   ! Reference solution
   call cpu_time(start_direct)
   call solver%solve(amat, xvec, expected, ainv, cpq, error=error)
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

!> Test 5: Ill-conditioned matrix
subroutine test_cg_ill_conditioned(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 8
   logical, parameter :: cpq = .false.
   real(wp) :: amat(n, n), xvec(n), vrhs(n), ainv(n, n)
   real(wp) :: expected(n), b(n), residual
   integer :: i

   ! Timer variables
   real(wp) :: start_cg, end_cg, start_direct, end_direct

   !> Solver variables
   real(wp), allocatable :: tol
   integer, allocatable :: maxiter 
   class(mchrg_solver_type), allocatable :: solver
   class(mchrg_solver_input), allocatable :: solver_input

   allocate(cg_input :: solver_input)
   select type(solver_input)
   type is (cg_input)
      if (allocated(maxiter)) then
         solver_input%cgmiter = maxiter
      end if
      if (allocated(tol)) then
         solver_input%cgtol = tol
      end if
         block
             class(mchrg_solver_cg), allocatable :: tmp
             allocate(tmp)
             call new_cg_solver(tmp, solver_input)
             call move_alloc(tmp, solver)
         end block
   end select

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
   call solver%solve(amat, xvec, vrhs, ainv, cpq, error=error)
   if (allocated(error)) return
   call cpu_time(end_cg)

   ! Reference solution
   deallocate(solver_input)
   deallocate(solver)
   allocate(direct_input :: solver_input)
   allocate(mchrg_solver_direct :: solver)
   select type (solver_input)
   type is (direct_input)
      call new_direct_solver(solver, solver_input)
   end select

   ! Reference solution
   call cpu_time(start_direct)
   call solver%solve(amat, xvec, expected, ainv, cpq, error=error)
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

!> Test 6: Zero RHS vector
subroutine test_cg_zero_rhs(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 5
   logical, parameter :: cpq = .false.
   real(wp) :: amat(n, n), xvec(n), vrhs(n), ainv(n, n)
   real(wp) :: expected(n)
   integer :: i

   ! Timer variables
   real(wp) :: start_cg, end_cg, start_direct, end_direct

   !> Solver variables
   real(wp), allocatable :: tol
   integer, allocatable :: maxiter 
   class(mchrg_solver_type), allocatable :: solver
   class(mchrg_solver_input), allocatable :: solver_input

   allocate(cg_input :: solver_input)
   select type(solver_input)
   type is (cg_input)
      if (allocated(maxiter)) then
         solver_input%cgmiter = maxiter
      end if
      if (allocated(tol)) then
         solver_input%cgtol = tol
      end if
         block
             class(mchrg_solver_cg), allocatable :: tmp
             allocate(tmp)
             call new_cg_solver(tmp, solver_input)
             call move_alloc(tmp, solver)
         end block
   end select

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
   call solver%solve(amat, xvec, vrhs, ainv, cpq, error=error)
   if (allocated(error)) return
   call cpu_time(end_cg)

   ! Reference solution
   deallocate(solver_input)
   deallocate(solver)
   allocate(direct_input :: solver_input)
   allocate(mchrg_solver_direct :: solver)
   select type (solver_input)
   type is (direct_input)
      call new_direct_solver(solver, solver_input)
   end select

   ! Reference solution
   call cpu_time(start_direct)
   call solver%solve(amat, xvec, expected, ainv, cpq, error=error)
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

!> Test 7: Random SPD matrix
subroutine test_cg_random_spd(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 100
   logical, parameter :: cpq = .false.
   real(wp), allocatable :: amat(:,:), xvec(:), vrhs(:), ainv(:,:)
   real(wp), allocatable :: expected(:), b(:), temp(:,:)
   integer :: i, j, k, seed_size
   integer, allocatable :: seed(:)
   real(wp) :: r, max_rel_error

   ! Timer variables
   real(wp) :: start_cg, end_cg, start_direct, end_direct

   !> Solver variables
   real(wp), allocatable :: tol
   integer, allocatable :: maxiter 
   class(mchrg_solver_type), allocatable :: solver
   class(mchrg_solver_input), allocatable :: solver_input

   allocate(cg_input :: solver_input)
   select type(solver_input)
   type is (cg_input)
      if (allocated(maxiter)) then
         solver_input%cgmiter = maxiter
      end if
      if (allocated(tol)) then
         solver_input%cgtol = tol
      end if
         block
             class(mchrg_solver_cg), allocatable :: tmp
             allocate(tmp)
             call new_cg_solver(tmp, solver_input)
             call move_alloc(tmp, solver)
         end block
   end select

   allocate(amat(n,n), xvec(n), vrhs(n), ainv(n,n), expected(n), b(n), temp(n,n))

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
   call solver%solve(amat, xvec, vrhs, ainv, cpq, error=error)
   if (allocated(error)) return
   call cpu_time(end_cg)

   ! Reference solution
   deallocate(solver_input)
   deallocate(solver)
   allocate(direct_input :: solver_input)
   allocate(mchrg_solver_direct :: solver)
   select type (solver_input)
   type is (direct_input)
      call new_direct_solver(solver, solver_input)
   end select

   ! Reference solution
   call cpu_time(start_direct)
   call solver%solve(amat, xvec, expected, ainv, cpq, error=error)
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


!> Test 8: Test with different preconditioner settings
subroutine test_cg_preconditioned(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 15
   logical, parameter :: cpq = .false.
   real(wp) :: amat(n, n), xvec(n), vrhs(n), ainv(n, n)
   real(wp) :: expected(n), b(n)
   integer :: i, j
   real(wp) :: residual_jacobi, residual_no_precond

   ! Timer variables
   real(wp) :: start_cg, end_cg, start_direct, end_direct

   !> Solver variables
   real(wp), allocatable :: tol
   integer, allocatable :: maxiter 
   class(mchrg_solver_type), allocatable :: solver
   class(mchrg_solver_input), allocatable :: solver_input

   allocate(cg_input :: solver_input)
   select type(solver_input)
   type is (cg_input)
      if (allocated(maxiter)) then
         solver_input%cgmiter = maxiter
      end if
      if (allocated(tol)) then
         solver_input%cgtol = tol
      end if
         block
             class(mchrg_solver_cg), allocatable :: tmp
             allocate(tmp)
             call new_cg_solver(tmp, solver_input)
             call move_alloc(tmp, solver)
         end block
   end select

   ! Create an SPD matrix with varying diagonal
   amat = 0.0_wp
   do i = 1, n
      amat(i,i) = 1000.0_wp ** ((i-1)/real(n-1, wp))  ! Diagonal from 1 to 1000
      do j = 1, n
         if (i /= j) then
            amat(i,j) = 0.1_wp / (abs(i-j) + 1.0_wp)
         end if
      end do
   end do
   
   ! Solution vector
   expected = [(sin(real(i, wp) * 0.5_wp), i=1, n)]
   
   ! Compute RHS
   b = 0.0_wp
   do i = 1, n
      do j = 1, n
         b(i) = b(i) + amat(i,j) * expected(j)
      end do
   end do
      
   vrhs = 0.0_wp
   xvec = b
   
   ! Solve
   call cpu_time(start_cg)
   call solver%solve(amat, xvec, vrhs, ainv, cpq, error=error)
   if (allocated(error)) return
   call cpu_time(end_cg)

   ! Reference solution
   deallocate(solver_input)
   deallocate(solver)
   allocate(direct_input :: solver_input)
   allocate(mchrg_solver_direct :: solver)
   select type (solver_input)
   type is (direct_input)
      call new_direct_solver(solver, solver_input)
   end select

   ! Reference solution
   call cpu_time(start_direct)
   call solver%solve(amat, xvec, expected, ainv, cpq, error=error)
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

end subroutine test_cg_preconditioned

!> Test 9: Test time scaling of the CG solver
subroutine test_cg_spd_time_scaling(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: n 
   integer, parameter :: max_size=13
   logical, parameter :: cpq = .false.
   real(wp), allocatable :: amat(:,:), xvec(:), vrhs(:), ainv(:,:)
   real(wp), allocatable :: expected(:), b(:)
   integer :: i, j, length

   ! Timer variables
   real(wp) :: start_cg, end_cg, start_direct, end_direct

   ! Scaling rows
   real(wp) :: iter_scal(max_size), dir_scal(max_size) 

   !> Solver variables
   real(wp), allocatable :: tol
   integer, allocatable :: maxiter 
   class(mchrg_solver_type), allocatable :: solver_test, solver_ref
   class(mchrg_solver_input), allocatable :: solver_test_input, solver_ref_input

   allocate(cg_input :: solver_test_input)
   select type(solver_test_input)
   type is (cg_input)
      if (allocated(maxiter)) then
         solver_test_input%cgmiter = maxiter
      end if
      if (allocated(tol)) then
         solver_test_input%cgtol = tol
      end if
         block
             class(mchrg_solver_cg), allocatable :: tmp
             allocate(tmp)
             call new_cg_solver(tmp, solver_test_input)
             call move_alloc(tmp, solver_test)
         end block
   end select

   ! Reference solution
   allocate(direct_input :: solver_ref_input)
   allocate(mchrg_solver_direct :: solver_ref)
   select type (solver_ref_input)
   type is (direct_input)
      call new_direct_solver(solver_ref, solver_ref_input)
   end select

   do length = 1, max_size
      n = 2**length
      allocate(amat(n,n), xvec(n), vrhs(n), ainv(n,n), expected(n), b(n))
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
      call solver_test%solve(amat, xvec, vrhs, ainv, cpq, error=error)
      if (allocated(error)) return
      call cpu_time(end_cg)
      
      ! Reference solution
      call cpu_time(start_direct)
      call solver_ref%solve(amat, xvec, expected, ainv, cpq, error=error)
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

      deallocate(amat, xvec, vrhs, ainv, expected, b)
   enddo

   call write_vector(iter_scal, "CG Time Scaling")
   call write_vector(dir_scal, "DIRECT Time Scaling")

end subroutine test_cg_spd_time_scaling

!> Test 10: Test time scaling of the CG solver for the "-1 2 -1" matrix
subroutine test_cg_121_time_scaling(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: n 
   integer, parameter :: max_size=12
   logical, parameter :: cpq = .false.
   real(wp), allocatable :: amat(:,:), xvec(:), vrhs(:), ainv(:,:)
   real(wp), allocatable :: expected(:), b(:)
   integer :: i, j, length

   ! Timer variables
   real(wp) :: start_cg, end_cg, start_direct, end_direct

   ! Scaling rows
   real(wp) :: iter_scal(max_size), dir_scal(max_size) 

   !> Solver variables
   real(wp), allocatable :: tol
   integer, allocatable :: maxiter 
   class(mchrg_solver_type), allocatable :: solver_test, solver_ref
   class(mchrg_solver_input), allocatable :: solver_test_input, solver_ref_input

   maxiter = 100000
   allocate(cg_input :: solver_test_input)
   select type(solver_test_input)
   type is (cg_input)
      if (allocated(maxiter)) then
         solver_test_input%cgmiter = maxiter
      end if
      if (allocated(tol)) then
         solver_test_input%cgtol = tol
      end if
         block
             class(mchrg_solver_cg), allocatable :: tmp
             allocate(tmp)
             call new_cg_solver(tmp, solver_test_input)
             call move_alloc(tmp, solver_test)
         end block
   end select

   ! Reference solution
   allocate(direct_input :: solver_ref_input)
   allocate(mchrg_solver_direct :: solver_ref)
   select type (solver_ref_input)
   type is (direct_input)
      call new_direct_solver(solver_ref, solver_ref_input)
   end select

   do length = 1, max_size
      n = 2**length
      allocate(amat(n,n), xvec(n), vrhs(n), ainv(n,n), expected(n), b(n))
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

      !call write_matrix(amat, 'Matrix')
      
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
      call solver_test%solve(amat, xvec, vrhs, ainv, cpq, error=error)
      if (allocated(error)) return
      call cpu_time(end_cg)
      
      ! Reference solution
      call cpu_time(start_direct)
      call solver_ref%solve(amat, xvec, expected, ainv, cpq, error=error)
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

      deallocate(amat, xvec, vrhs, ainv, expected, b)
   enddo

   call write_vector(iter_scal, "CG Time Scaling")
   call write_vector(dir_scal, "DIRECT Time Scaling")

end subroutine test_cg_121_time_scaling

end module test_solver
