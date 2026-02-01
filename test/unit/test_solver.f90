module test_solver
   use mctc_env, only: wp
   use mctc_env_testing, only: new_unittest, unittest_type, error_type, test_failed
   use mctc_io_structure, only: structure_type, new
   use mstore, only: get_structure
   use multicharge_model, only: mchrg_model_type
   use multicharge_model_eeqbc, only: eeqbc_model
   use multicharge_param, only: new_eeq2019_model, new_eeqbc2025_model
   use multicharge_model_cache, only: cache_container
   use multicharge_charge, only: get_charges, get_eeq_charges, get_eeqbc_charges
   use solver, only : mchrg_solver_type
   use solver_factory, only : new_mchrg_solver
   implicit none
   private

   public :: collect_solver

   real(wp), parameter :: thr = 1.0e-10_wp
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
      & new_unittest("cg-random-spd", test_cg_random_spd), &
      & new_unittest("cg-preconditioned", test_cg_preconditioned) &
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
   class(mchrg_solver_type), allocatable :: solver
   character(len=:), allocatable :: solver_choice
   integer, allocatable :: cgmiter
   real(wp), allocatable :: cgtol
   character(len=32), allocatable :: cgmode

   ! Setup solver
   solver_choice = "CG"
   cgmiter = 1000
   cgtol = 1.0e-12_wp
   cgmode = "default"
   
   solver = new_mchrg_solver(solver_choice, cgmiter, cgtol, cgmode)

   ! Identity matrix
   amat = 0.0_wp
   amat(1,1) = 1.0_wp
   amat(2,2) = 1.0_wp
   
   ! RHS vector
   xvec = [1.0_wp, 2.0_wp]
   
   ! Initial guess
   vrhs = [0.0_wp, 0.0_wp]
     
   ! Solve iteratively
   call solver%solve(amat, xvec, vrhs, ainv, cpq, error=error)
   if (allocated(error)) return

   ! Reference solution
   solver_choice = "DIRECT"
   expected = [0.0_wp, 0.0_wp]
   solver = new_mchrg_solver(solver_choice, cgmiter, cgtol, cgmode)
   call solver%solve(amat, xvec, expected, ainv, cpq, error=error)
   
   ! Check solution
   if (any(abs(vrhs - expected) > thr)) then
      call test_failed(error, "CG solver failed for identity matrix")
      print'(a)', "Solution:"
      print'(3es21.14)', vrhs
      print'(a)', "Expected:"
      print'(3es21.14)', expected
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
   class(mchrg_solver_type), allocatable :: solver
   character(len=:), allocatable :: solver_choice
   integer, allocatable :: cgmiter
   real(wp), allocatable :: cgtol
   character(len=32), allocatable :: cgmode
   integer :: i

   ! Setup solver
   solver_choice = "CG"
   cgmiter = 1000
   cgtol = 1.0e-12_wp
   cgmode = "default"
   
   solver = new_mchrg_solver(solver_choice, cgmiter, cgtol, cgmode)

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
   call solver%solve(amat, xvec, vrhs, ainv, cpq, error=error)
   if (allocated(error)) return

   ! Reference solution
   solver_choice = "DIRECT"
   expected = [0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
   solver = new_mchrg_solver(solver_choice, cgmiter, cgtol, cgmode)
   call solver%solve(amat, xvec, expected, ainv, cpq, error=error)
   
   ! Check solution
   if (any(abs(vrhs - expected) > thr)) then
      call test_failed(error, "CG solver failed for diagonal matrix")
      print'(a)', "Solution:"
      print'(3es21.14)', vrhs
      print'(a)', "Expected:"
      print'(3es21.14)', expected
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
   class(mchrg_solver_type), allocatable :: solver
   character(len=:), allocatable :: solver_choice
   integer, allocatable :: cgmiter
   real(wp), allocatable :: cgtol
   character(len=32), allocatable :: cgmode

   ! Setup solver
   solver_choice = "CG"
   cgmiter = 1000
   cgtol = 1.0e-12_wp
   cgmode = "default"
   
   solver = new_mchrg_solver(solver_choice, cgmiter, cgtol, cgmode)

   ! SPD matrix: A = [4 1 1; 1 3 2; 1 2 4]
   amat = reshape([4.0_wp, 1.0_wp, 1.0_wp, &
                  &1.0_wp, 3.0_wp, 2.0_wp, &
                  &1.0_wp, 2.0_wp, 4.0_wp], [n, n])
   
   ! RHS vector
   xvec = [6.0_wp, 6.0_wp, 7.0_wp]
   
   ! Initial guess
   vrhs = [0.0_wp, 0.0_wp, 0.0_wp]
   
   ! Expected solution (precomputed)
   expected = [1.0_wp, 1.0_wp, 1.0_wp]
   
   ! Solve
   call solver%solve(amat, xvec, vrhs, ainv, cpq, error=error)
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

!> Test 4: Medium SPD matrix
subroutine test_cg_spd_medium(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 10
   logical, parameter :: cpq = .false.
   real(wp) :: amat(n, n), xvec(n), vrhs(n), ainv(n, n)
   real(wp) :: expected(n), b(n), residual
   class(mchrg_solver_type), allocatable :: solver
   character(len=:), allocatable :: solver_choice
   integer, allocatable :: cgmiter
   real(wp), allocatable :: cgtol
   character(len=32), allocatable :: cgmode
   integer :: i, j

   ! Setup solver
   solver_choice = "CG"
   cgmiter = 2000
   cgtol = 1.0e-10_wp
   cgmode = "default"
   
   solver = new_mchrg_solver(solver_choice, cgmiter, cgtol, cgmode)

   ! Create a simple SPD matrix: A = I + 0.1*E where E is matrix of ones
   amat = 0.0_wp
   do i = 1, n
      amat(i,i) = 1.0_wp
      do j = 1, n
         if (i /= j) then
            amat(i,j) = 0.1_wp
         end if
      end do
   end do
   
   ! Create a solution vector
   expected = [(real(i, wp), i=1, n)]
   
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
   call solver%solve(amat, xvec, vrhs, ainv, cpq, error=error)
   if (allocated(error)) return
   
   ! Check solution
   if (any(abs(vrhs - expected) / max(1.0_wp, abs(expected)) > thr_rel)) then
      call test_failed(error, "CG solver failed for medium SPD matrix")
      print'(a)', "Relative error:"
      print'(3es21.14)', abs(vrhs - expected) / max(1.0_wp, abs(expected))
   end if

end subroutine test_cg_spd_medium

!> Test 5: Large SPD matrix (1000x1000)
subroutine test_cg_spd_large(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 1000
   logical, parameter :: cpq = .false.
   real(wp), allocatable :: amat(:,:), xvec(:), vrhs(:), ainv(:,:)
   real(wp), allocatable :: expected(:), b(:)
   class(mchrg_solver_type), allocatable :: solver
   character(len=:), allocatable :: solver_choice
   integer, allocatable :: cgmiter
   real(wp), allocatable :: cgtol
   character(len=32), allocatable :: cgmode
   integer :: i, j

   ! Setup solver
   solver_choice = "CG"
   cgmiter = 5000
   cgtol = 1.0e-8_wp
   cgmode = "default"
   
   solver = new_mchrg_solver(solver_choice, cgmiter, cgtol, cgmode)

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
   call solver%solve(amat, xvec, vrhs, ainv, cpq, error=error)
   if (allocated(error)) return
   
   ! Check solution with relative tolerance
   if (any(abs(vrhs - expected) / max(1.0_wp, abs(expected)) > thr_rel)) then
      call test_failed(error, "CG solver failed for large SPD matrix")
      print'(a,i0)', "Maximum relative error: ", &
           maxval(abs(vrhs - expected) / max(1.0_wp, abs(expected)))
   end if

end subroutine test_cg_spd_large

!> Test 6: Ill-conditioned matrix
subroutine test_cg_ill_conditioned(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 8
   logical, parameter :: cpq = .false.
   real(wp) :: amat(n, n), xvec(n), vrhs(n), ainv(n, n)
   real(wp) :: expected(n), b(n), residual
   class(mchrg_solver_type), allocatable :: solver
   character(len=:), allocatable :: solver_choice
   integer, allocatable :: cgmiter
   real(wp), allocatable :: cgtol
   character(len=32), allocatable :: cgmode
   integer :: i

   ! Setup solver
   solver_choice = "CG"
   cgmiter = 5000
   cgtol = 1.0e-10_wp
   cgmode = "default"
   
   solver = new_mchrg_solver(solver_choice, cgmiter, cgtol, cgmode)

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
   call solver%solve(amat, xvec, vrhs, ainv, cpq, error=error)
   if (allocated(error)) return
   
   ! Check solution with relative tolerance
   if (any(abs(vrhs - expected) / max(1.0_wp, abs(expected)) > thr_rel)) then
      call test_failed(error, "CG solver failed for ill-conditioned matrix")
      print'(a)', "Condition numbers range from 1 to 10^7"
      print'(a,es10.2)', "Max relative error: ", &
           maxval(abs(vrhs - expected) / max(1.0_wp, abs(expected)))
   end if

end subroutine test_cg_ill_conditioned

!> Test 7: Zero RHS vector
subroutine test_cg_zero_rhs(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 5
   logical, parameter :: cpq = .false.
   real(wp) :: amat(n, n), xvec(n), vrhs(n), ainv(n, n)
   real(wp) :: expected(n)
   class(mchrg_solver_type), allocatable :: solver
   character(len=:), allocatable :: solver_choice
   integer, allocatable :: cgmiter
   real(wp), allocatable :: cgtol
   character(len=32), allocatable :: cgmode
   integer :: i

   ! Setup solver
   solver_choice = "CG"
   cgmiter = 1000
   cgtol = 1.0e-12_wp
   cgmode = "default"
   
   solver = new_mchrg_solver(solver_choice, cgmiter, cgtol, cgmode)

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
   call solver%solve(amat, xvec, vrhs, ainv, cpq, error=error)
   if (allocated(error)) return
   
   ! Check solution
   if (any(abs(vrhs) > thr)) then
      call test_failed(error, "CG solver failed for zero RHS vector")
      print'(a)', "Solution should be zero vector:"
      print'(3es21.14)', vrhs
   end if

end subroutine test_cg_zero_rhs

!> Test 8: Random SPD matrix
subroutine test_cg_random_spd(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 1000
   logical, parameter :: cpq = .false.
   real(wp), allocatable :: amat(:,:), xvec(:), vrhs(:), ainv(:,:)
   real(wp), allocatable :: expected(:), b(:), temp(:,:)
   class(mchrg_solver_type), allocatable :: solver
   character(len=:), allocatable :: solver_choice
   integer, allocatable :: cgmiter
   real(wp), allocatable :: cgtol
   character(len=32), allocatable :: cgmode
   integer :: i, j, k, seed_size
   integer, allocatable :: seed(:)
   real(wp) :: r, max_rel_error

   ! Setup solver
   solver_choice = "CG"
   cgmiter = 3000
   cgtol = 1.0e-8_wp
   cgmode = "default"
   
   solver = new_mchrg_solver(solver_choice, cgmiter, cgtol, cgmode)

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
   call solver%solve(amat, xvec, vrhs, ainv, cpq, error=error)
   if (allocated(error)) return
   
   ! Check solution with relative tolerance
   max_rel_error = maxval(abs(vrhs - expected) / max(1.0_wp, abs(expected)))
   if (max_rel_error > thr_rel) then
      call test_failed(error, "CG solver failed for random SPD matrix")
      print'(a,es10.2)', "Maximum relative error: ", max_rel_error
   end if

end subroutine test_cg_random_spd


!> Test 9: Test with different preconditioner settings
subroutine test_cg_preconditioned(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 15
   logical, parameter :: cpq = .false.
   real(wp) :: amat(n, n), xvec(n), vrhs(n), ainv(n, n)
   real(wp) :: expected(n), b(n)
   class(mchrg_solver_type), allocatable :: solver
   character(len=:), allocatable :: solver_choice
   integer, allocatable :: cgmiter
   real(wp), allocatable :: cgtol
   character(len=32), allocatable :: cgmode
   integer :: i, j
   real(wp) :: residual_jacobi, residual_no_precond

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
   
   ! Test 1: With Jacobi preconditioner (default)
   solver_choice = "CG"
   cgmiter = 1000
   cgtol = 1.0e-10_wp
   cgmode = "default"
   
   solver = new_mchrg_solver(solver_choice, cgmiter, cgtol, cgmode)
   
   vrhs = 0.0_wp
   xvec = b
   
   call solver%solve(amat, xvec, vrhs, ainv, cpq, error=error)
   if (allocated(error)) return
   
   residual_jacobi = sqrt(sum((vrhs - expected)**2)) / sqrt(sum(expected**2))
   
   ! Note: The current implementation always uses Jacobi preconditioning
   ! This test verifies it works correctly
   if (residual_jacobi > thr_rel) then
      call test_failed(error, "CG solver with Jacobi preconditioner failed")
      print'(a,es10.2)', "Relative residual with Jacobi: ", residual_jacobi
   end if

end subroutine test_cg_preconditioned

end module test_solver