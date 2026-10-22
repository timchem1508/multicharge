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

!> Unit tests for CSR-list based multicharge calculations
module test_csrlist
   use iso_fortran_env, only : output_unit
   use mctc_env, only : wp, i8
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, test_failed
   use mctc_cutoff, only : get_lattice_points
   use mctc_io_structure, only : structure_type, new
   use mctc_csrlist, only : csr_list, new_csr_list
   use mstore, only : get_structure
   use multicharge_blas, only : gemv
   use multicharge_model_type, only : mchrg_model_type, mchrg_hessian, get_dqdr_row, &
      & new_partial_derivs
   use multicharge_model_eeqbc, only : eeqbc_model
   use multicharge_param, only : new_eeq2019_model, new_eeqbc2025_model
   use multicharge_model_cache, only : mchrg_cache
   use mctc_ncoord, only : ncoord_type
   use multicharge_charge, only : get_charges, get_eeq_charges, get_eeqbc_charges
   use multicharge_solver_type, only : mchrg_solver_type, mchrg_solver_input
   use multicharge_solver_direct, only : direct_solver, new_direct_solver, direct_input
   use multicharge_solver_cg, only : cg_solver, new_cg_solver, cg_input
   implicit none
   private

   public :: collect_csrlist

   !> Tight tolerance for direct numerical comparisons
   real(wp), parameter :: thr = 100 * epsilon(1.0_wp)

   !> Tolerance for neighborlist and dense-matrix comparisons
   real(wp), parameter :: thr1 = 1.0e5_wp*epsilon(1.0_wp)

   !> Tolerance for numerical derivatives
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

   !> neighborlist cutoff used by shared test helpers
   real(wp), parameter :: cutoff = 29.0_wp

   !> Cutoff of a Hessian list holding all pairs of the test structures
   real(wp), parameter :: allpairs = 200.0_wp


contains


!> Collect all exported unit tests
subroutine collect_csrlist(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
   & new_unittest("eeqbc-components-mb01", test_eeqbc_components_mb01), &
   & new_unittest("eeqbc-charges-mb01", test_eeqbc_q_mb01), &
   & new_unittest("eeqbc-charges-mb02", test_eeqbc_q_mb02), &
   & new_unittest("eeqbc-charges-actinides", test_eeqbc_q_actinides), &
   & new_unittest("eeqbc-energy-mb03", test_eeqbc_e_mb03), &
   & new_unittest("eeqbc-energy-mb04", test_eeqbc_e_mb04), &
   & new_unittest("eeqbc-energy-co2", test_eeqbc_e_co2),  &
   & new_unittest("eeqbc-energy-ice", test_eeqbc_e_ice), &
   & new_unittest("eeqbc-energy-ice-supercell", test_eeqbc_e_ice222), &
   & new_unittest("eeqbc-gradient-mb05", test_eeqbc_g_mb05), &
   & new_unittest("eeqbc-gradient-mb06", test_eeqbc_g_mb06), &
   & new_unittest("eeqbc-gradient-co2", test_eeqbc_g_co2), &
   & new_unittest("eeqbc-gradient-ice", test_eeqbc_g_ice), &
   & new_unittest("eeqbc-gradient-ice-supercell", test_eeqbc_g_ice222), &
   & new_unittest("eeqbc-api-mb01", test_eeqbc_api_mb01), &
   & new_unittest("eeqbc-api-mb02", test_eeqbc_api_mb02), &
   & new_unittest("eeqbc-api-actinides", test_eeqbc_api_actinides), &
   & new_unittest("partial-derivs-mb01", test_partial_derivs_mb01), &
   & new_unittest("partial-derivs-co2", test_partial_derivs_co2), &
   & new_unittest("partial-derivs-ice", test_partial_derivs_ice), &
   & new_unittest("d2cdr2-mb01", test_d2cdr2_mb01), &
   & new_unittest("d2cdr2-mb05", test_d2cdr2_mb05), &
   & new_unittest("explicit-hessian-coulomb-mb01", test_hexpl_coulomb_mb01), &
   & new_unittest("explicit-hessian-en-mb01", test_hexpl_en_mb01), &
   & new_unittest("explicit-hessian-mb05", test_hexpl_mb05), &
   & new_unittest("explicit-hessian-dfdq-mb06", test_hexpl_dfdq_mb06), &
   & new_unittest("explicit-hessian-znooh", test_hexpl_znooh), &
   & new_unittest("explicit-hessian-sparse-upu23", test_hexpl_sparse_upu23), &
   & new_unittest("explicit-hessian-sparse-dfdq-upu23", test_hexpl_sparse_dfdq_upu23), &
   & new_unittest("implicit-hessian-list-mb05", test_himpl_list_mb05), &
   & new_unittest("implicit-hessian-sparse-list-upu23", test_himpl_sparse_list_upu23), &
   & new_unittest("hessian-list-mb05", test_hess_list_mb05), &
   & new_unittest("hessian-list-znooh", test_hess_list_znooh), &
   & new_unittest("hessian-pattern-list-upu23", test_hess_pattern_list_upu23), &
   & new_unittest("external-hessian-list-mb05", test_hext_list_mb05), &
   & new_unittest("external-hessian-list-znooh", test_hext_list_znooh), &
   & new_unittest("external-hessian-sparse-list-upu23", test_hext_sparse_list_upu23) &
   & ]

end subroutine collect_csrlist

!> Construct a solver from its polymorphic input configuration
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

!> Replicate a structure along each lattice-vector direction
subroutine make_supercell(mol, rep)

   !> Structure to replicate
   type(structure_type), intent(inout) :: mol

   !> Replication factors along the three lattice vectors
   integer, intent(in) :: rep(3)

   real(wp), allocatable :: xyz(:, :), lattice(:, :)
   integer, allocatable :: num(:)
   integer :: i, j, k, c

   num = reshape(spread([mol%num(mol%id)], 2, product(rep)), [product(rep)*mol%nat])
   lattice = reshape(&
      [rep(1)*mol%lattice(:, 1), rep(2)*mol%lattice(:, 2), rep(3)*mol%lattice(:, 3)], &
      shape(mol%lattice))
   allocate(xyz(3, product(rep)*mol%nat))
   c = 0
   do i = 0, rep(1)-1
      do j = 0, rep(2)-1
         do k = 0, rep(3)-1
            xyz(:, c+1:c+mol%nat) = mol%xyz &
            & + spread(matmul(mol%lattice, [real(wp):: i, j, k]), 2, mol%nat)
            c = c + mol%nat
         end do
      end do
   end do

   call new(mol, num, xyz, lattice=lattice)
end subroutine make_supercell

subroutine gen_test_molecular(error, mol, model, qref, eref)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Reference charges
   real(wp), intent(in), optional :: qref(:)

   !> Reference energies
   real(wp), intent(in), optional :: eref(:)

   type(mchrg_cache), allocatable :: cache

   ! Solver variables
   class(mchrg_solver_type), allocatable :: solver
   class(mchrg_solver_input), allocatable :: solver_input
   real(wp) :: tol = 1.0e-15_wp
   integer :: maxiter = 1000
   integer :: verbosity = 0

   type(csr_list), allocatable :: list

   real(wp) :: trans(3, 1) = 0.0_wp
   real(wp), allocatable :: energy(:)
   real(wp), allocatable :: qvec(:)

   allocate(cg_input :: solver_input)
   select type (solver_input)
      type is (cg_input)
      solver_input%cgtol = tol
      solver_input%cgmiter = maxiter
      solver_input%verbosity = verbosity
      solver_input%use_nlist = .true.
   end select
   call solver_maker(solver, solver_input, error)

   allocate(cache)

   if (present(eref)) then
      allocate (energy(mol%nat))
      energy(:) = 0.0_wp
   end if
   if (present(qref)) then
      allocate (qvec(mol%nat))
   end if

   ! Build adjacency list
   allocate(list)
   call new_csr_list(list, mol, error, cutoff=cutoff)

   call model%update(mol, cache, trans, grad=.false., list=list)
   call model%solve(mol, solver, cache, error, energy=energy, qvec=qvec, list=list, &
   & unit=output_unit)
   if (allocated(error)) return

   if (present(qref)) then
      if (any(abs(qvec - qref) > thr)) then
         call test_failed(error, "Partial charges do not match")
         print'(a)', "Charges:"
         print'(3es21.14)', qvec
         print'(a)', "diff:"
         print'(3es21.14)', qvec - qref
      end if
   end if
   if (allocated(error)) return

   if (present(eref)) then
      if (any(abs(energy - eref) > thr)) then
         call test_failed(error, "Energies do not match")
         print'(a)', "Energy:"
         print'(3es21.14)', energy
         print'(a)', "diff:"
         print'(3es21.14)', energy - eref
      end if
   end if

end subroutine gen_test_molecular

subroutine test_components(error, mol, model)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model


   type(mchrg_cache), allocatable :: cache1, cache2

   ! Solver variables
   class(mchrg_solver_type), allocatable :: solver
   class(mchrg_solver_input), allocatable :: solver_input
   real(wp) :: tol = 1.0e-15_wp
   integer :: maxiter = 1000
   integer :: verbosity = 0

   type(csr_list), allocatable :: list

   real(wp) :: trans(3, 1) = 0.0_wp
   real(wp), allocatable :: cmat(:, :), amat(:, :)

   integer :: iat, jat, ndim
   integer(i8) :: kat

   allocate(cg_input :: solver_input)
   select type (solver_input)
      type is (cg_input)
      solver_input%cgtol = tol
      solver_input%cgmiter = maxiter
      solver_input%verbosity = verbosity
      solver_input%use_nlist = .true.
      ndim = mol%nat
   end select
   call solver_maker(solver, solver_input, error)

   allocate(cache1)

   call model%update(mol, cache1, trans, grad=.true.)
   call model%solve(mol, solver, cache1, error, unit=output_unit)
   if (allocated(error)) return

   allocate(cache2)
   allocate(list)

   call new_csr_list(list, mol, error, cutoff=cutoff)
   call model%update(mol, cache2, trans, grad=.true., list=list)

   call model%solve(mol, solver, cache2, error, list=list, unit=output_unit)
   if (allocated(error)) return

   allocate(cmat(mol%nat, mol%nat), source=0.0_wp)
   allocate(amat(mol%nat, mol%nat), source=0.0_wp)
   do iat = 1, mol%nat
      do kat = list%inl(iat), list%inl(iat+1) - 1
         jat = list%nlat(kat)
         cmat(iat, jat) = cache2%clist(kat)
         cmat(jat, iat) = cache2%clist(kat)
         amat(iat, jat) = cache2%alist(kat)
         amat(jat, iat) = cache2%alist(kat)
      end do
   end do

   if (any(abs(cmat(:, :) - cache1%cmat(:, :)) > thr1)) then
      call test_failed(error, "C-matrix does not match")
      print'(a)', "C-matrix:"
      print'(3es21.14)', cmat
      print'(a)', "numcmat:"
      print'(3es21.14)', cache1%cmat
      print'(a)', "diff:"
      print'(3es21.14)', cmat - cache1%cmat
   end if

   if (any(abs(cache1%xvec(:) - cache2%xvec(:)) > thr1)) then
      call test_failed(error, "x-vector does not match")
      print'(a)', "x-vector:"
      print'(3es21.14)', cache1%xvec
      print'(a)', "numxvec:"
      print'(3es21.14)', cache2%xvec
      print'(a)', "diff:"
      print'(3es21.14)', cache1%xvec - cache2%xvec
   end if

   if (any(abs(amat(:, :) - cache1%amat(:, :)) > thr1)) then
      call test_failed(error, "A-matrix does not match")
      print'(a)', "A-matrix:"
      print'(3es21.14)', amat
      print'(a)', "numamat:"
      print'(3es21.14)', cache1%amat
      print'(a)', "diff:"
      print'(3es21.14)', amat - cache1%amat
   end if


end subroutine test_components

subroutine gen_test_periodic(error, mol, model)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Reference charges
   real(wp), allocatable :: qref(:)

   !> Reference energies
   real(wp), allocatable :: eref(:)

   type(mchrg_cache), allocatable :: cache
   type(csr_list), allocatable :: list

   ! Solver variables
   class(mchrg_solver_type), allocatable :: solver
   class(mchrg_solver_input), allocatable :: solver_input
   real(wp) :: tol = 1.0e-15_wp
   integer :: maxiter = 1000
   integer :: verbosity = 0

   integer :: ndim
   real(wp), allocatable :: trans(:, :)
   real(wp), allocatable :: energy(:)
   real(wp), allocatable :: qvec(:)
   real(wp), allocatable :: amat_dir(:,:)
   real(wp), allocatable :: amat_list(:,:)

   allocate(cg_input :: solver_input)
   select type (solver_input)
      type is (cg_input)
      solver_input%cgtol = tol
      solver_input%cgmiter = maxiter
      solver_input%verbosity = verbosity
      ndim = mol%nat
      solver_input%use_nlist = .true.
   end select
   call solver_maker(solver, solver_input, error)

   allocate(cache)

   call get_lattice_points(mol%periodic, mol%lattice, model%ncoord%cutoff, trans)

   allocate(eref(mol%nat))
   eref(:) = 0.0_wp
   allocate(energy(mol%nat))
   energy(:) = 0.0_wp
   allocate (qref(mol%nat))
   allocate (qvec(mol%nat))
   allocate(amat_dir(mol%nat, mol%nat), amat_list(mol%nat, mol%nat))

   call model%update(mol, cache, trans, grad=.false.)
   call model%solve(mol, solver, cache, error, energy=eref, qvec=qref, unit=output_unit)


   if (allocated(error)) return

   ! Build adjacency list
   deallocate(cache)
   allocate(cache)
   allocate(list)
   call new_csr_list(list, mol, error, cache%wsc, 29.0_wp)
   call model%update(mol, cache, trans, grad=.false., list=list)
   call model%solve(mol, solver, cache, error, energy=energy, qvec=qvec, list=list, &
   & unit=output_unit)

   if (any(abs(qvec - qref) > thr1)) then
      call test_failed(error, "Partial charges do not match")
      print'(a)', "Charges:"
      print'(3es21.14)', qvec
      print'(a)', "diff:"
      print'(3es21.14)', qvec - qref
   end if

   if (any(abs(energy - eref) > thr1)) then
      call test_failed(error, "Energies do not match")
      print'(a)', "Energy:"
      print'(3es21.14)', energy
      print'(a)', "diff:"
      print'(3es21.14)', energy - eref
   end if
   if (allocated(error)) return


end subroutine gen_test_periodic

subroutine test_grad(error, mol, model)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   type(mchrg_cache), allocatable :: cache1, cache2

   ! Solver variables
   class(mchrg_solver_type), allocatable :: solver
   class(mchrg_solver_input), allocatable :: solver_input
   real(wp) :: tol = 1.0e-15_wp
   integer :: maxiter = 1000
   integer :: verbosity = 0

   type(csr_list), allocatable :: list

   integer :: ndim
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), allocatable :: energy(:), gradient(:, :), sigma(:, :)
   real(wp), allocatable :: numgrad(:, :), numsigma(:, :)

   allocate(cg_input :: solver_input)
   select type (solver_input)
      type is (cg_input)
      solver_input%cgtol = tol
      solver_input%cgmiter = maxiter
      solver_input%verbosity = verbosity
      solver_input%use_nlist = .true.
      ndim = mol%nat
   end select
   call solver_maker(solver, solver_input, error)

   allocate(cache1)

   allocate (energy(mol%nat), gradient(3, mol%nat), sigma(3, 3))
   energy(:) = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   ! Build adjacency list
   allocate(numsigma(3, 3), source=0.0_wp)
   allocate(numgrad(3, mol%nat), source=0.0_wp)

   call model%update(mol, cache1, trans, grad=.true.)
   call model%solve(mol, solver, cache1, error, &
   & gradient=numgrad, sigma=numsigma, unit=output_unit)
   if (allocated(error)) return

   allocate(cache2)
   allocate(list)
   gradient = 0.0_wp
   sigma(:, :) = 0.0_wp
   call new_csr_list(list, mol, error, cutoff=cutoff)
   call model%update(mol, cache2, trans, grad=.true., list=list)

   call model%solve(mol, solver, cache2, error, &
   & gradient=gradient, sigma=sigma, list=list, unit=output_unit)
   if (allocated(error)) return


   if (any(abs(gradient(:, :) - numgrad(:, :)) > thr2)) then
      call test_failed(error, "Derivative of energy does not match")
      print'(a)', "Energy gradient:"
      print'(3es21.14)', gradient
      print'(a)', "numgrad:"
      print'(3es21.14)', numgrad
      print'(a)', "diff:"
      print'(3es21.14)', gradient - numgrad
   end if

   if (any(abs(sigma(:, :) - numsigma(:, :)) > thr2)) then
      call test_failed(error, "Derivative of energy does not match")
      print'(a)', "Energy sigma:"
      print'(3es21.14)', sigma
      print'(a)', "numsigma:"
      print'(3es21.14)', numsigma
      print'(a)', "diff:"
      print'(3es21.14)', sigma - numsigma
   end if

end subroutine test_grad

subroutine test_grad_periodic(error, mol, model)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   type(mchrg_cache), allocatable :: cache1, cache2

   ! Solver variables
   class(mchrg_solver_type), allocatable :: solver
   class(mchrg_solver_input), allocatable :: solver_input
   real(wp) :: tol = 1.0e-15_wp
   integer :: maxiter = 1000
   integer :: verbosity = 0

   type(csr_list), allocatable :: list

   integer :: iat, ndim
   real(wp), allocatable :: energy(:), gradient(:, :), sigma(:, :)
   real(wp), allocatable :: numgrad(:, :), numsigma(:, :)
   real(wp), allocatable :: trans(:, :)
   real(wp), allocatable :: dcdrdiag(:, :)

   allocate(cg_input :: solver_input)
   select type (solver_input)
      type is (cg_input)
      solver_input%cgtol = tol
      solver_input%cgmiter = maxiter
      solver_input%verbosity = verbosity
      solver_input%use_nlist = .true.
      ndim = mol%nat
   end select
   call solver_maker(solver, solver_input, error)

   allocate(cache1)

   allocate (energy(mol%nat), gradient(3, mol%nat), sigma(3, 3))
   energy(:) = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   ! Build adjacency list
   allocate(numsigma(3, 3), source=0.0_wp)
   allocate(numgrad(3, mol%nat), source=0.0_wp)

   call get_lattice_points(mol%periodic, mol%lattice, model%ncoord%cutoff, trans)

   call model%update(mol, cache1, trans, grad=.true.)
   call model%solve(mol, solver, cache1, error, &
   & gradient=numgrad, sigma=numsigma, unit=output_unit)
   allocate(dcdrdiag(3, mol%nat), source = 0.0_wp)
   do iat = 1, mol%nat
      dcdrdiag(:, iat) = cache1%dcdr(:, iat, iat)
   end do

   if (allocated(error)) return

   allocate(cache2)
   allocate(list)
   gradient = 0.0_wp
   sigma(:, :) = 0.0_wp
   call new_csr_list(list, mol, error, cache2%wsc, cutoff=29.0_wp)
   call model%update(mol, cache2, trans, grad=.true., list=list)

   call model%solve(mol, solver, cache2, error, &
   & gradient=gradient, sigma=sigma, list=list, unit=output_unit)
   if (allocated(error)) return

   if (any(abs(gradient(:, :) - numgrad(:, :)) > thr2)) then
      call test_failed(error, "Derivative of energy does not match")
      print'(a)', "Energy gradient:"
      print'(3es21.14)', gradient
      print'(a)', "numgrad:"
      print'(3es21.14)', numgrad
      print'(a)', "diff:"
      print'(3es21.14)', gradient - numgrad
   end if

   if (any(abs(sigma(:, :) - numsigma(:, :)) > thr2)) then
      call test_failed(error, "Derivative of energy does not match")
      print'(a)', "Energy sigma:"
      print'(3es21.14)', sigma
      print'(a)', "numsigma:"
      print'(3es21.14)', numsigma
      print'(a)', "diff:"
      print'(3es21.14)', sigma - numsigma
   end if

end subroutine test_grad_periodic

subroutine test_api(error, mol, model, qref)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Reference charges
   real(wp), intent(in) :: qref(:)

   type(csr_list), allocatable :: list
   real(wp), allocatable :: qvec(:)

   allocate(list)
   call new_csr_list(list, mol, error, cutoff=cutoff)

   ! Check wrapper functions
   allocate (qvec(mol%nat), source=0.0_wp)
   call get_charges(model, mol, error, qvec, list=list)
   if (allocated(error)) return

   if (any(abs(qvec - qref) > thr)) then
      call test_failed(error, "Partial charges do not match")
      print'(a)', "Charges:"
      print'(3es21.14)', qvec
      print'(a)', "diff:"
      print'(3es21.14)', qvec - qref
   end if
   if (allocated(error)) return

   qvec = 0.0_wp
   call get_eeqbc_charges(mol, error, qvec, list=list)
   if (allocated(error)) return

   if (any(abs(qvec - qref) > thr)) then
      call test_failed(error, "Partial charges do not match")
      print'(a)', "Charges:"
      print'(3es21.14)', qvec
      print'(a)', "diff:"
      print'(3es21.14)', qvec - qref
   end if
   if (allocated(error)) return

end subroutine test_api

subroutine test_eeqbc_components_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "01")

   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_components(error, mol, model)

end subroutine test_eeqbc_components_mb01

subroutine test_eeqbc_q_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   real(wp), parameter :: ref(16) = [&
   & 6.32665177635486E-1_wp, -7.79000021394748E-3_wp, -7.28560267968198E-1_wp, &
   &-4.10042373033769E-2_wp, -4.60026363641469E-1_wp,  2.20004260276381E-1_wp, &
   &-4.51954034871547E-2_wp, -7.57357565689225E-1_wp, -5.27004217539908E-1_wp, &
   & 3.09566865491112E-1_wp,  2.27494681048596E-1_wp, -3.45550645258975E-1_wp, &
   &-3.99729823174914E-2_wp,  8.48742713551116E-1_wp, -4.49898375805073E-1_wp, &
   & 1.16388636122213E+0_wp]

   call get_structure(mol, "MB16-43", "01")

   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call gen_test_molecular(error, mol, model, qref=ref)

end subroutine test_eeqbc_q_mb01

subroutine test_eeqbc_q_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   real(wp), parameter :: ref(16) = [&
   &-1.68706023409942E-2_wp, -4.31585754567721E-1_wp, -5.63457264808067E-3_wp, &
   &-6.88598679593678E-1_wp,  8.95507266279121E-1_wp,  2.70560711967719E-1_wp, &
   & 1.06704514101413E-2_wp, -1.87359068857452E-3_wp,  2.75916162155502E-1_wp, &
   & 2.68351655726932E-1_wp,  4.83228968763888E-3_wp,  5.84779537866739E-1_wp, &
   &-6.15231036823021E-1_wp,  1.47042131652275E-2_wp, -1.44439879995398E-2_wp, &
   &-5.51084063597411E-1_wp]

   call get_structure(mol, "MB16-43", "02")

   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call gen_test_molecular(error, mol, model, qref=ref)

end subroutine test_eeqbc_q_mb02

subroutine test_eeqbc_q_actinides(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   real(wp), parameter :: ref(17) = [&
   & 4.80669255571553E-2_wp, -2.23681443476404E-1_wp,  2.95879870118081E-1_wp, &
   & 1.40528551895345E-1_wp,  1.54496956558730E-1_wp, -1.31765078691778E-1_wp, &
   & 6.43024168002473E-2_wp, -3.08519563494370E-1_wp, -2.76459245598841E-1_wp, &
   & 1.78293689441428E-1_wp, -2.11018657500951E-1_wp, -1.03628773361279E-1_wp, &
   &-1.71248308078648E-1_wp,  2.54400229067594E-1_wp, -5.83023049918706E-2_wp, &
   & 2.01328580342047E-1_wp,  1.47326155413513E-1_wp]

   call get_structure(mol, "f-block", "Fr_to_Lr")

   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call gen_test_molecular(error, mol, model, qref=ref)

end subroutine test_eeqbc_q_actinides

subroutine test_eeqbc_e_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   real(wp), parameter :: ref(16) = [&
   &-6.85383909094970E-2_wp, -2.19074341545678E+0_wp, -4.74601764501974E-3_wp, &
   &-7.87937923890688E-1_wp, -3.03527925157768E+0_wp, -2.82220139737265E-1_wp, &
   &-5.40497002004504E-1_wp, -5.41240617940733E-3_wp, -8.05028961511096E-3_wp, &
   &-1.75953554576722E-2_wp, -1.37304680376243E+0_wp, -3.53988490287759E-1_wp, &
   &-9.94227733791934E-1_wp, -1.58035207823740E-1_wp, -1.51308961030806E+0_wp, &
   &-4.01533041771599E-3_wp]

   call get_structure(mol, "MB16-43", "03")

   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call gen_test_molecular(error, mol, model, eref=ref)

end subroutine test_eeqbc_e_mb03

subroutine test_eeqbc_e_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   real(wp), parameter :: ref(16) = [&
   &-3.11006035589390E-3_wp, -5.28599433581519E-2_wp, -9.35709515842972E-5_wp, &
   &-6.16264514697450E-1_wp, -1.87462481346083E+0_wp, -2.21662455395278E-3_wp, &
   &-2.51039203171614E-2_wp, -8.93995363870797E-3_wp, -1.62270805124076E-4_wp, &
   &-9.27331161591396E-4_wp, -4.55713551016470E-1_wp, -1.49805355046533E+0_wp, &
   &-1.06116580759546E-2_wp, -1.45309833473282E+0_wp, -2.30657216723778E-2_wp, &
   & 1.53766210354731E-6_wp]

   call get_structure(mol, "MB16-43", "04")

   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call gen_test_molecular(error, mol, model, eref=ref)

end subroutine test_eeqbc_e_mb04

subroutine test_eeqbc_g_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "05")

   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_grad(error, mol, model)

end subroutine test_eeqbc_g_mb05

subroutine test_eeqbc_g_mb06(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "06")

   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_grad(error, mol, model)

end subroutine test_eeqbc_g_mb06

subroutine test_eeqbc_g_co2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "X23", "CO2")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_grad_periodic(error, mol, model)

end subroutine test_eeqbc_g_co2

subroutine test_eeqbc_g_ice(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "ICE10", "vi")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_grad_periodic(error, mol, model)

end subroutine test_eeqbc_g_ice

subroutine test_eeqbc_e_ice(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "ICE10", "vi")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call gen_test_periodic(error, mol, model)

end subroutine test_eeqbc_e_ice

subroutine test_eeqbc_e_co2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "X23", "CO2")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call gen_test_periodic(error, mol, model)

end subroutine test_eeqbc_e_co2

subroutine test_eeqbc_e_ice222(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   integer, parameter :: supercell(*) = [2, 2, 2]

   call get_structure(mol, "ICE10", "vi")
   call make_supercell(mol, supercell)
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call gen_test_periodic(error, mol, model)

end subroutine test_eeqbc_e_ice222

subroutine test_eeqbc_g_ice222(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   integer, parameter :: supercell(*) = [2, 2, 2]

   call get_structure(mol, "ICE10", "vi")
   call make_supercell(mol, supercell)
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_grad_periodic(error, mol, model)

end subroutine test_eeqbc_g_ice222

subroutine test_eeqbc_api_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   real(wp), parameter :: ref(16) = [&
   & 6.32665177635486E-1_wp, -7.79000021394748E-3_wp, -7.28560267968198E-1_wp, &
   &-4.10042373033769E-2_wp, -4.60026363641469E-1_wp,  2.20004260276381E-1_wp, &
   &-4.51954034871547E-2_wp, -7.57357565689225E-1_wp, -5.27004217539908E-1_wp, &
   & 3.09566865491112E-1_wp,  2.27494681048596E-1_wp, -3.45550645258975E-1_wp, &
   &-3.99729823174914E-2_wp,  8.48742713551116E-1_wp, -4.49898375805073E-1_wp, &
   & 1.16388636122213E+0_wp]

   call get_structure(mol, "MB16-43", "01")

   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_api(error, mol, model, qref=ref)

end subroutine test_eeqbc_api_mb01

subroutine test_eeqbc_api_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   real(wp), parameter :: ref(16) = [&
   &-1.68706023409942E-2_wp, -4.31585754567721E-1_wp, -5.63457264808067E-3_wp, &
   &-6.88598679593678E-1_wp,  8.95507266279121E-1_wp,  2.70560711967719E-1_wp, &
   & 1.06704514101413E-2_wp, -1.87359068857452E-3_wp,  2.75916162155502E-1_wp, &
   & 2.68351655726932E-1_wp,  4.83228968763888E-3_wp,  5.84779537866739E-1_wp, &
   &-6.15231036823021E-1_wp,  1.47042131652275E-2_wp, -1.44439879995398E-2_wp, &
   &-5.51084063597411E-1_wp]

   call get_structure(mol, "MB16-43", "02")

   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_api(error, mol, model, qref=ref)

end subroutine test_eeqbc_api_mb02

subroutine test_eeqbc_api_actinides(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   real(wp), parameter :: ref(17) = [&
   & 4.80669255571553E-2_wp, -2.23681443476404E-1_wp,  2.95879870118081E-1_wp, &
   & 1.40528551895345E-1_wp,  1.54496956558730E-1_wp, -1.31765078691778E-1_wp, &
   & 6.43024168002473E-2_wp, -3.08519563494370E-1_wp, -2.76459245598841E-1_wp, &
   & 1.78293689441428E-1_wp, -2.11018657500951E-1_wp, -1.03628773361279E-1_wp, &
   &-1.71248308078648E-1_wp,  2.54400229067594E-1_wp, -5.83023049918706E-2_wp, &
   & 2.01328580342047E-1_wp,  1.47326155413513E-1_wp]

   call get_structure(mol, "f-block", "Fr_to_Lr")

   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_api(error, mol, model, qref=ref)

end subroutine test_eeqbc_api_actinides

!> Expand the half neighbour list storage of the partial derivatives
subroutine expand_partial_derivs(mol, list, dlist_ut, dlist_lt, dabdr)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Half CSR neighbour list
   type(csr_list), intent(in) :: list

   !> Diagonal and upper triangle of the partial derivatives
   real(wp), intent(in) :: dlist_ut(:, :)

   !> Lower triangle of the partial derivatives
   real(wp), intent(in) :: dlist_lt(:, :)

   !> Dense partial derivatives
   real(wp), intent(out) :: dabdr(:, :, :)

   integer :: iat
   integer(i8) :: kat

   ! The entry of a pair holds both components in the two triangles
   dabdr(:, :, :) = 0.0_wp
   do iat = 1, mol%nat
      do kat = list%inl(iat), list%inl(iat + 1) - 1
         dabdr(:, list%nlat(kat), iat) = dlist_ut(:, kat)
         if (list%nlat(kat) /= iat) &
            & dabdr(:, iat, list%nlat(kat)) = dlist_lt(:, kat)
      end do
   end do

end subroutine expand_partial_derivs


!> Compare the partial derivatives and charge derivatives obtained with the
!> half neighbour list against the dense evaluation with the block CG and the
!> direct solver
subroutine test_partial_derivs(error, mol, model)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   class(mchrg_solver_type), allocatable :: solver, direct
   class(mchrg_solver_input), allocatable :: solver_input, direct_in
   type(mchrg_cache), allocatable :: cache1, cache2, cache3, cache4
   type(csr_list), allocatable :: list
   integer :: ndim, iscale, jat
   real(wp), parameter :: scales(2, 3) = reshape([1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, &
      & -1.0_wp, 1.0_wp], [2, 3])
   real(wp), allocatable :: trans(:, :), dabdr(:, :, :), pvec(:)
   real(wp), allocatable :: dqdr1(:, :, :), dqdL1(:, :, :)
   real(wp), allocatable :: dqdr2(:, :, :), dqdL2(:, :, :)
   real(wp), allocatable :: dqdr3(:, :, :), dqdL3(:, :, :)
   real(wp), allocatable :: dqrow(:, :), dqdL4(:, :, :)

   allocate(cg_input :: solver_input)
   select type (solver_input)
   type is (cg_input)
      solver_input%cgtol = 1.0e-15_wp
      solver_input%cgmiter = 1000
      solver_input%verbosity = 0
      solver_input%use_nlist = .true.
   end select
   ndim = mol%nat
   call solver_maker(solver, solver_input, error)
   if (allocated(error)) return

   allocate(direct_input :: direct_in)
   select type (direct_in)
   type is (direct_input)
      direct_in%verbosity = 0
   end select
   call solver_maker(direct, direct_in, error)
   if (allocated(error)) return

   call get_lattice_points(mol%periodic, mol%lattice, model%ncoord%cutoff, trans)
   allocate(dabdr(3, mol%nat, ndim))
   allocate(dqdr1(3, mol%nat, mol%nat), dqdL1(3, 3, mol%nat))
   allocate(dqdr2(3, mol%nat, mol%nat), dqdL2(3, 3, mol%nat))
   allocate(dqdr3(3, mol%nat, mol%nat), dqdL3(3, 3, mol%nat))
   allocate(dqrow(3, mol%nat), dqdL4(3, 3, mol%nat))

   ! Dense reference with the block CG solver
   allocate(cache1)
   call model%update(mol, cache1, trans, grad=.true., dqdr=.true.)
   call model%solve(mol, solver, cache1, error, dqdr=dqdr1, dqdL=dqdL1, &
      & unit=output_unit)
   if (allocated(error)) return

   ! Dense reference with the direct solver
   allocate(cache3)
   call model%update(mol, cache3, trans, grad=.true., dqdr=.true.)
   call model%solve(mol, direct, cache3, error, dqdr=dqdr3, dqdL=dqdL3, &
      & unit=output_unit)
   if (allocated(error)) return

   ! Half neighbour list, inverse on the pattern of the list
   allocate(cache2, list)
   if (any(mol%periodic)) then
      call new_csr_list(list, mol, error, cache2%wsc, cutoff=cutoff, complete=.false.)
   else
      call new_csr_list(list, mol, error, cutoff=cutoff, complete=.false.)
   end if
   if (allocated(error)) return
   call model%update(mol, cache2, trans, grad=.true., list=list, dqdr=.true.)
   call model%solve(mol, solver, cache2, error, dqdr=dqdr2, dqdL=dqdL2, list=list, &
      & unit=output_unit)
   if (allocated(error)) return

   if (any(abs(dqdr2 - dqdr1) > thr2) .or. any(abs(dqdr2 - dqdr3) > thr2)) then
      call test_failed(error, "Charge derivatives w.r.t. positions do not match")
      print'(2es21.14)', maxval(abs(dqdr2 - dqdr1)), maxval(abs(dqdr2 - dqdr3))
      return
   end if
   if (any(abs(dqdL2 - dqdL1) > thr2) .or. any(abs(dqdL2 - dqdL3) > thr2)) then
      call test_failed(error, "Charge derivatives w.r.t. strain do not match")
      print'(2es21.14)', maxval(abs(dqdL2 - dqdL1)), maxval(abs(dqdL2 - dqdL3))
      return
   end if

   ! Position derivatives only in compressed storage, without a dense dqdr,
   ! fresh cache keeps the Wigner-Seitz cell of the neighbour list
   allocate(cache4)
   cache4%wsc = cache2%wsc
   call move_alloc(cache4, cache2)
   call model%update(mol, cache2, trans, grad=.true., list=list, dqdr=.true.)
   call model%solve(mol, solver, cache2, error, dqdL=dqdL4, list=list, &
      & unit=output_unit)
   if (allocated(error)) return

   do jat = 1, mol%nat
      call get_dqdr_row(cache2, list, jat, dqrow)
      if (any(abs(dqrow - dqdr2(:, :, jat)) > thr1)) then
         call test_failed(error, "Compressed charge derivatives do not match")
         print'(es21.14)', maxval(abs(dqrow - dqdr2(:, :, jat)))
         return
      end if
   end do
   if (any(abs(dqdL4 - dqdL2) > thr1)) then
      call test_failed(error, "Charge derivatives w.r.t. strain without dqdr differ")
      print'(es21.14)', maxval(abs(dqdL4 - dqdL2))
      return
   end if

   ! Coulomb matrix and electronegativity parts separately and added to the
   ! same partial derivatives, same charges
   cache2%vrhs(:) = cache1%vrhs
   do iscale = 1, 3
      call new_partial_derivs(mol, ndim, cache1)
      call model%get_xvec_derivs(mol, ndim, cache1, scales(2, iscale))
      call model%get_coulomb_derivs(mol, ndim, cache1, scales(1, iscale), &
         & cache1%vrhs(:mol%nat))
      call new_partial_derivs(mol, ndim, cache2, list)
      call model%get_xvec_derivs(mol, ndim, cache2, scales(2, iscale), list=list)
      call model%get_coulomb_derivs(mol, ndim, cache2, scales(1, iscale), &
         & cache2%vrhs(:mol%nat), list=list)
      call expand_partial_derivs(mol, list, cache2%dabdrlist_ut, &
         & cache2%dabdrlist_lt, dabdr)

      if (any(abs(dabdr - cache1%dabdr) > thr2)) then
         call test_failed(error, "Partial derivatives w.r.t. positions do not match")
         print'(es21.14)', maxval(abs(dabdr - cache1%dabdr))
         return
      end if
      if (any(abs(cache2%dabdL - cache1%dabdL) > thr2)) then
         call test_failed(error, "Partial derivatives w.r.t. strain do not match")
         print'(es21.14)', maxval(abs(cache2%dabdL - cache1%dabdL))
         return
      end if
   end do

   ! Coulomb matrix derivatives contracted with an external vector, the partial
   ! derivatives of the last combination stay as they are
   allocate(pvec(mol%nat))
   do jat = 1, mol%nat
      pvec(jat) = cos(real(jat, wp)) - 0.3_wp
   end do
   call new_partial_derivs(mol, ndim, cache1, dafdr=.true.)
   call model%get_coulomb_derivs(mol, ndim, cache1, 0.5_wp, pvec, dafdr=.true.)
   call new_partial_derivs(mol, ndim, cache2, list, dafdr=.true.)
   call model%get_coulomb_derivs(mol, ndim, cache2, 0.5_wp, pvec, list=list, &
      & dafdr=.true.)
   call expand_partial_derivs(mol, list, cache2%dafdrlist_ut, cache2%dafdrlist_lt, &
      & dabdr)
   if (any(abs(dabdr - cache1%dafdr) > thr1)) then
      call test_failed(error, "Coulomb matrix derivatives times p do not match")
      print'(es21.14)', maxval(abs(dabdr - cache1%dafdr))
      return
   end if
   call expand_partial_derivs(mol, list, cache2%dabdrlist_ut, cache2%dabdrlist_lt, &
      & dabdr)
   if (any(abs(dabdr - cache1%dabdr) > thr2)) then
      call test_failed(error, "Partial derivatives changed by dA/dR*p")
      print'(es21.14)', maxval(abs(dabdr - cache1%dabdr))
      return
   end if

end subroutine test_partial_derivs


subroutine test_partial_derivs_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "01")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_partial_derivs(error, mol, model)

end subroutine test_partial_derivs_mb01


subroutine test_partial_derivs_co2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "X23", "CO2")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_partial_derivs(error, mol, model)

end subroutine test_partial_derivs_co2


subroutine test_partial_derivs_ice(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "ICE10", "vi")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_partial_derivs(error, mol, model)

end subroutine test_partial_derivs_ice


!> Check the diagonal blocks of the second derivative of the capacitance matrix
!> from the neighbour list against finite differences of the diagonal of its
!> first derivative
subroutine test_d2cdr2(error, mol, model)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   type(mchrg_cache), allocatable :: cache
   type(csr_list), allocatable :: list
   integer :: iat, ic, ndim
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), parameter :: step = 1.0e-6_wp
   real(wp) :: dcdrr1(3), dcdrr2(3), dcdrl1(3), dcdrl2(3)
   real(wp), allocatable :: numd2cdr2(:, :, :)

   allocate(cache, list)
   allocate(numd2cdr2(3, 3, mol%nat), source=0.0_wp)
   ndim = mol%nat

   call new_csr_list(list, mol, error, cutoff=cutoff, complete=.false.)
   if (allocated(error)) return

   ! The capacitance matrix depends on the positions only
   call model%update(mol, cache, trans, grad=.true., list=list)

   do iat = 1, mol%nat
      do ic = 1, 3
         ! First right-hand side (x+h)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call model%get_capacitance_matrix(mol, ndim, cache, list=list)
         dcdrr1 = cache%dcdrdiag(:, iat)

         ! Second right-hand side (x+2h)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call model%get_capacitance_matrix(mol, ndim, cache, list=list)
         dcdrr2 = cache%dcdrdiag(:, iat)

         ! First left-hand side (x-h)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 3*step
         call model%get_capacitance_matrix(mol, ndim, cache, list=list)
         dcdrl1 = cache%dcdrdiag(:, iat)

         ! Second left-hand side (x-2h)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - step
         call model%get_capacitance_matrix(mol, ndim, cache, list=list)
         dcdrl2 = cache%dcdrdiag(:, iat)

         ! Return to original position
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + 2*step

         ! f'(x) ≈ [f(x-2h) - 8f(x-h) + 8f(x+h) - f(x+2h)]/(12h)
         numd2cdr2(:, ic, iat) = (dcdrl2 - 8.0_wp*dcdrl1 + 8.0_wp*dcdrr1 &
            & - dcdrr2)/(12.0_wp*step)
      end do
   end do

   ! Analytical second derivative
   allocate(cache%d2cdr2diag(3, 3, mol%nat))
   select type (model)
   type is (eeqbc_model)
      call model%get_d2cmat_0d_list(mol, list, cache)
   class default
      call test_failed(error, "Capacitance matrix requires the EEQBC model")
      return
   end select

   if (any(abs(cache%d2cdr2diag - numd2cdr2) > thr2)) then
      call test_failed(error, "Second derivative of the C matrix does not match")
      print'(a)', "d2cdr2diag:"
      print'(3es21.12)', cache%d2cdr2diag
      print'(a)', "numd2cdr2:"
      print'(3es21.12)', numd2cdr2
      print'(a)', "diff:"
      print'(3es21.12)', cache%d2cdr2diag - numd2cdr2
   end if

end subroutine test_d2cdr2


subroutine test_d2cdr2_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "01")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_d2cdr2(error, mol, model)

end subroutine test_d2cdr2_mb01


subroutine test_d2cdr2_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "05")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_d2cdr2(error, mol, model)

end subroutine test_d2cdr2_mb05


!> Check the explicit Hessian, the second derivative of alpha * p^T A q
!> + beta * p^T X at fixed p and charges q, with a half neighbour list against
!> central differences of the gradient with the same list, against the dense
!> evaluation if the list holds all pairs, and its symmetry and translational
!> invariance
subroutine test_explicit_hessian(error, mol, model, alpha, beta, nlcutoff, pvec)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Prefactor of the Coulomb matrix term
   real(wp), intent(in) :: alpha

   !> Prefactor of the electronegativity term
   real(wp), intent(in) :: beta

   !> Cutoff of the neighbour list
   real(wp), intent(in) :: nlcutoff

   !> Contraction vector, the charges are used if absent
   real(wp), intent(in), optional :: pvec(:)

   class(mchrg_solver_type), allocatable :: solver
   class(mchrg_solver_input), allocatable :: solver_input
   type(mchrg_cache), allocatable :: cache
   type(csr_list), allocatable :: list
   type(mchrg_hessian), allocatable :: sphess
   integer :: iat, ic, ndim
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), parameter :: step = 1.0e-4_wp
   real(wp) :: xorig
   real(wp), allocatable :: qvec(:), p(:), hessian(:, :), refhess(:, :), numhess(:, :)
   real(wp), allocatable :: gr1(:, :), gr2(:, :), gl1(:, :), gl2(:, :)

   allocate(cg_input :: solver_input)
   select type (solver_input)
   type is (cg_input)
      solver_input%cgtol = 1.0e-15_wp
      solver_input%cgmiter = 1000
      solver_input%verbosity = 0
   end select
   call solver_maker(solver, solver_input, error)
   if (allocated(error)) return

   ndim = 3*mol%nat
   allocate(qvec(mol%nat), hessian(ndim, ndim), refhess(ndim, ndim), &
      & numhess(ndim, ndim))
   allocate(gr1(3, mol%nat), gr2(3, mol%nat), gl1(3, mol%nat), gl2(3, mol%nat))

   ! Charges of the reference structure without a neighbour list
   allocate(cache)
   call model%update(mol, cache, trans, grad=.true., hess=.true.)
   call model%solve(mol, solver, cache, error, qvec=qvec, unit=output_unit)
   if (allocated(error)) return
   if (present(pvec)) then
      p = pvec
   else
      p = qvec
   end if

   allocate(list)
   call new_csr_list(list, mol, error, cutoff=nlcutoff, complete=.false.)
   if (allocated(error)) return

   refhess(:, :) = 0.0_wp
   hessian(:, :) = 0.0_wp
   select type (model)
   type is (eeqbc_model)
      ! Dense evaluation from the cache of the charges
      call model%get_explicit_hessian(mol, cache, p, refhess, error, &
         & alpha=alpha, beta=beta)
      if (allocated(error)) return

      ! Neighbour list evaluation at the same charges, the list of the Hessian
      ! holds all pairs to keep the products of the CN derivatives with the
      ! second neighbours
      call fixed_charge_cache(mol, model, list, qvec, cache)
      call new_allpairs_hessian(mol, sphess, error)
      if (allocated(error)) return
      call model%get_explicit_hessian(mol, cache, p, sphess, error, &
         & alpha=alpha, beta=beta, list=list)
      if (allocated(error)) return
      call expand_hessian(mol, sphess, hessian)
   class default
      call test_failed(error, "Explicit Hessian requires the EEQBC model")
      return
   end select

   ! Without all pairs in the list the dense evaluation is a different model
   if (size(list%nlat) == mol%nat * (mol%nat + 1) / 2) then
      if (any(abs(hessian - refhess) > thr1)) then
         call test_failed(error, "Explicit Hessian does not match the dense one")
         print'(a)', "diff:"
         print'(3es21.12)', hessian - refhess
         return
      end if
   end if

   ! Central differences of the analytical gradient with the same list
   do iat = 1, mol%nat
      do ic = 1, 3
         xorig = mol%xyz(ic, iat)
         mol%xyz(ic, iat) = xorig + step
         call fixed_charge_gradient(mol, model, list, qvec, p, alpha, beta, gr1)
         mol%xyz(ic, iat) = xorig + 2*step
         call fixed_charge_gradient(mol, model, list, qvec, p, alpha, beta, gr2)
         mol%xyz(ic, iat) = xorig - step
         call fixed_charge_gradient(mol, model, list, qvec, p, alpha, beta, gl1)
         mol%xyz(ic, iat) = xorig - 2*step
         call fixed_charge_gradient(mol, model, list, qvec, p, alpha, beta, gl2)
         mol%xyz(ic, iat) = xorig

         ! f'(x) ≈ [f(x-2h) - 8f(x-h) + 8f(x+h) - f(x+2h)]/(12h)
         numhess(:, 3*(iat - 1) + ic) = reshape((gl2 - 8.0_wp*gl1 + 8.0_wp*gr1 &
            & - gr2)/(12.0_wp*step), [ndim])
      end do
   end do

   if (any(abs(hessian - numhess) > thr2)) then
      call test_failed(error, "Explicit Hessian does not match")
      print'(a)', "hessian:"
      print'(3es21.12)', hessian
      print'(a)', "numhess:"
      print'(3es21.12)', numhess
      print'(a)', "diff:"
      print'(3es21.12)', hessian - numhess
      return
   end if

   if (any(abs(hessian - transpose(hessian)) > thr1)) then
      call test_failed(error, "Explicit Hessian is not symmetric")
      print'(a)', "asymmetry:"
      print'(3es21.12)', hessian - transpose(hessian)
      return
   end if

   ! A rigid translation of all atoms leaves the expression unchanged
   do ic = 1, 3
      if (any(abs(sum(hessian(ic::3, :), dim=1)) > thr1)) then
         call test_failed(error, "Explicit Hessian is not translationally invariant")
         print'(a)', "sum over atoms:"
         print'(3es21.12)', sum(hessian(ic::3, :), dim=1)
         return
      end if
   end do

end subroutine test_explicit_hessian


!> Cache for the neighbour list evaluation at fixed charges, including the CN
!> and local charge derivatives of the Hessian
subroutine fixed_charge_cache(mol, model, list, qvec, cache)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Half neighbour list
   type(csr_list), intent(in) :: list

   !> Fixed atomic partial charges
   real(wp), intent(in) :: qvec(:)

   !> Multicharge cache
   type(mchrg_cache), intent(out) :: cache

   real(wp), parameter :: trans(3, 1) = 0.0_wp

   call model%update(mol, cache, trans, grad=.true., list=list, hess=.true.)
   call model%get_capacitance_matrix(mol, mol%nat, cache, list=list)
   call model%get_xvec(mol, mol%nat, cache, list=list)
   cache%vrhs = qvec

end subroutine fixed_charge_cache


!> Analytical gradient of alpha * p^T A q + beta * p^T X at fixed p and charges
!> q for a displaced structure using a fresh cache and the given neighbour list
subroutine fixed_charge_gradient(mol, model, list, qvec, p, alpha, beta, gradient)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Half neighbour list
   type(csr_list), intent(in) :: list

   !> Fixed atomic partial charges
   real(wp), intent(in) :: qvec(:)

   !> Contraction vector
   real(wp), intent(in) :: p(:)

   !> Prefactor of the Coulomb matrix term
   real(wp), intent(in) :: alpha

   !> Prefactor of the electronegativity term
   real(wp), intent(in) :: beta

   !> Gradient of the expression
   real(wp), intent(out) :: gradient(:, :)

   type(mchrg_cache), allocatable :: cache
   real(wp) :: sigma(3, 3)

   allocate(cache)
   call fixed_charge_cache(mol, model, list, qvec, cache)

   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp
   call model%get_grad(mol, cache, p, gradient, sigma, alpha=alpha, beta=beta, &
      & list=list)

end subroutine fixed_charge_gradient


!> Check the implicit Hessian with a half neighbour list against the dense
!> product of the expanded partial and charge derivatives of the same list,
!> and against the dense evaluation if the list holds all pairs
subroutine test_implicit_hessian(error, mol, model, nlcutoff)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Cutoff of the neighbour list
   real(wp), intent(in) :: nlcutoff

   class(mchrg_solver_type), allocatable :: solver
   class(mchrg_solver_input), allocatable :: solver_input
   type(mchrg_cache), allocatable :: cache
   type(csr_list), allocatable :: list
   type(mchrg_hessian), allocatable :: sphess
   integer :: jat, ndim, xat, yat, ii, jj
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp) :: blk(3, 3)
   real(wp), allocatable :: dabdr(:, :, :), dqdr(:, :, :), dqdL(:, :, :)
   real(wp), allocatable :: hessian(:, :), refhess(:, :)

   allocate(cg_input :: solver_input)
   select type (solver_input)
   type is (cg_input)
      solver_input%cgtol = 1.0e-15_wp
      solver_input%cgmiter = 1000
      solver_input%verbosity = 0
      solver_input%use_nlist = .true.
   end select
   call solver_maker(solver, solver_input, error)
   if (allocated(error)) return

   ndim = 3*mol%nat
   allocate(dabdr(3, mol%nat, mol%nat), dqdr(3, mol%nat, mol%nat), &
      & dqdL(3, 3, mol%nat), hessian(ndim, ndim), refhess(ndim, ndim))

   allocate(cache, list)
   call new_csr_list(list, mol, error, cutoff=nlcutoff, complete=.false.)
   if (allocated(error)) return
   call model%update(mol, cache, trans, grad=.false., list=list, dqdr=.true.)
   call model%solve(mol, solver, cache, error, dqdL=dqdL, list=list, &
      & unit=output_unit)
   if (allocated(error)) return

   call new_allpairs_hessian(mol, sphess, error)
   if (allocated(error)) return
   select type (model)
   type is (eeqbc_model)
      call model%get_implicit_hessian(mol, cache, sphess, list=list)
   class default
      call test_failed(error, "Implicit Hessian requires the EEQBC model")
      return
   end select
   call expand_hessian(mol, sphess, hessian)

   ! The list stores the rows of each atom in its columns, so the reference is
   ! formed in the same order, -(dq/dR)^T (dX/dR - dA/dR*q)
   call expand_partial_derivs(mol, list, cache%dabdrlist_ut, cache%dabdrlist_lt, &
      & dabdr)
   do jat = 1, mol%nat
      call get_dqdr_row(cache, list, jat, dqdr(:, :, jat))
   end do
   refhess(:, :) = -matmul(reshape(dqdr, [ndim, mol%nat]), &
      & transpose(reshape(dabdr, [ndim, mol%nat])))

   ! With the cutoff the product is not symmetric. The Hessian holds its blocks
   ! in the rows of atom x and the columns of the atoms y > x, their transposes
   ! below the diagonal and the symmetrized diagonal blocks.
   do xat = 1, mol%nat
      ii = 3 * (xat - 1)
      blk(:, :) = refhess(ii+1:ii+3, ii+1:ii+3)
      refhess(ii+1:ii+3, ii+1:ii+3) = 0.5_wp * (blk + transpose(blk))
      do yat = xat + 1, mol%nat
         jj = 3 * (yat - 1)
         refhess(jj+1:jj+3, ii+1:ii+3) = transpose(refhess(ii+1:ii+3, jj+1:jj+3))
      end do
   end do

   if (any(abs(hessian - refhess) > thr1)) then
      call test_failed(error, "Implicit Hessian does not match the expanded product")
      print'(a)', "diff:"
      print'(3es21.12)', hessian - refhess
      return
   end if

   ! Without all pairs in the list the dense evaluation is a different model
   if (size(list%nlat) == mol%nat * (mol%nat + 1) / 2) then
      deallocate(cache)
      allocate(cache)
      call model%update(mol, cache, trans, grad=.false., dqdr=.true.)
      call model%solve(mol, solver, cache, error, dqdr=dqdr, dqdL=dqdL, &
         & unit=output_unit)
      if (allocated(error)) return

      refhess(:, :) = 0.0_wp
      select type (model)
      type is (eeqbc_model)
         call model%get_implicit_hessian(mol, cache, refhess, dqdr=dqdr)
      end select

      if (any(abs(hessian - refhess) > thr2)) then
         call test_failed(error, "Implicit Hessian does not match the dense one")
         print'(a)', "diff:"
         print'(3es21.12)', hessian - refhess
      end if
   end if

end subroutine test_implicit_hessian


!> Check the analytical Hessian with a half neighbour list holding all pairs
!> against the dense evaluation, and its symmetry
subroutine test_hessian(error, mol, model)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   class(mchrg_solver_type), allocatable :: solver
   class(mchrg_solver_input), allocatable :: solver_input
   type(mchrg_cache), allocatable :: cache
   type(csr_list), allocatable :: list
   type(mchrg_hessian), allocatable :: sphess
   integer :: ndim
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), allocatable :: hessian(:, :), refhess(:, :)

   allocate(cg_input :: solver_input)
   select type (solver_input)
   type is (cg_input)
      solver_input%cgtol = 1.0e-15_wp
      solver_input%cgmiter = 1000
      solver_input%verbosity = 0
      solver_input%use_nlist = .true.
   end select
   call solver_maker(solver, solver_input, error)
   if (allocated(error)) return

   ndim = 3*mol%nat
   allocate(hessian(ndim, ndim), refhess(ndim, ndim))

   allocate(list)
   call new_csr_list(list, mol, error, cutoff=cutoff, complete=.false.)
   if (allocated(error)) return
   if (size(list%nlat) /= mol%nat * (mol%nat + 1) / 2) then
      call test_failed(error, "Neighbour list does not hold all pairs")
      return
   end if

   select type (model)
   type is (eeqbc_model)
      allocate(cache)
      call model%update(mol, cache, trans, grad=.false., hess=.true.)
      call model%solve(mol, solver, cache, error)
      if (allocated(error)) return
      call model%get_hessian(mol, solver, cache, error, refhess)
      if (allocated(error)) return

      ! The list of the Hessian is a copy of the neighbour list
      deallocate(cache)
      allocate(cache, sphess)
      call model%update(mol, cache, trans, grad=.false., hess=.true., list=list)
      call model%solve(mol, solver, cache, error, list=list)
      if (allocated(error)) return
      call model%get_hessian(mol, solver, cache, error, sphess, list=list)
      if (allocated(error)) return
   class default
      call test_failed(error, "Analytical Hessian requires the EEQBC model")
      return
   end select
   call expand_hessian(mol, sphess, hessian)

   if (any(abs(hessian - refhess) > thr2)) then
      call test_failed(error, "Hessian does not match the dense one")
      print'(a)', "diff:"
      print'(3es21.12)', hessian - refhess
      return
   end if

   if (any(abs(hessian - transpose(hessian)) > thr2)) then
      call test_failed(error, "Hessian is not symmetric")
      print'(a)', "asymmetry:"
      print'(3es21.12)', hessian - transpose(hessian)
      return
   end if

end subroutine test_hessian


!> Check that the compressed Hessian on a short list of its own holds the
!> blocks of the Hessian on a list with all pairs, and the same factors of the
!> charge constraint term
subroutine test_hessian_pattern(error, mol, model, nlcutoff, hcutoff)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Cutoff of the neighbour list
   real(wp), intent(in) :: nlcutoff

   !> Cutoff of the short list of the Hessian
   real(wp), intent(in) :: hcutoff

   class(mchrg_solver_type), allocatable :: solver
   class(mchrg_solver_input), allocatable :: solver_input
   type(mchrg_cache), allocatable :: cache
   type(csr_list), allocatable :: list
   type(mchrg_hessian), allocatable :: refhess, sphess
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   integer :: xat
   integer(i8) :: hat, kat
   logical :: found

   allocate(cg_input :: solver_input)
   select type (solver_input)
   type is (cg_input)
      solver_input%cgtol = 1.0e-15_wp
      solver_input%cgmiter = 1000
      solver_input%verbosity = 0
      solver_input%use_nlist = .true.
   end select
   call solver_maker(solver, solver_input, error)
   if (allocated(error)) return

   allocate(list)
   call new_csr_list(list, mol, error, cutoff=nlcutoff, complete=.false.)
   if (allocated(error)) return

   allocate(sphess)
   call new_csr_list(sphess%list, mol, error, cutoff=hcutoff, complete=.false.)
   if (allocated(error)) return
   if (size(sphess%list%nlat) >= mol%nat * (mol%nat + 1) / 2) then
      call test_failed(error, "Short list of the Hessian holds all pairs")
      return
   end if

   select type (model)
   type is (eeqbc_model)
      allocate(cache)
      call new_allpairs_hessian(mol, refhess, error)
      if (allocated(error)) return
      call model%update(mol, cache, trans, grad=.false., hess=.true., list=list)
      call model%solve(mol, solver, cache, error, list=list)
      if (allocated(error)) return
      call model%get_hessian(mol, solver, cache, error, refhess, list)
      if (allocated(error)) return

      deallocate(cache)
      allocate(cache)
      call model%update(mol, cache, trans, grad=.false., hess=.true., list=list)
      call model%solve(mol, solver, cache, error, list=list)
      if (allocated(error)) return
      call model%get_hessian(mol, solver, cache, error, sphess, list)
      if (allocated(error)) return
   class default
      call test_failed(error, "Analytical Hessian requires the EEQBC model")
      return
   end select

   if (abs(sphess%weight - refhess%weight) > thr1 &
      & .or. any(abs(sphess%dlambda - refhess%dlambda) > thr1)) then
      call test_failed(error, "Charge constraint term depends on the list of the Hessian")
      return
   end if

   do xat = 1, mol%nat
      do hat = sphess%list%inl(xat), sphess%list%inl(xat + 1) - 1
         found = .false.
         do kat = refhess%list%inl(xat), refhess%list%inl(xat + 1) - 1
            if (refhess%list%nlat(kat) /= sphess%list%nlat(hat)) cycle
            found = .true.
            if (any(abs(sphess%hesslist(:, :, hat) - refhess%hesslist(:, :, kat)) &
               & > thr1)) then
               call test_failed(error, "Hessian block does not match the full list")
               print'(a, 2i6)', "block", xat, sphess%list%nlat(hat)
               print'(3es21.12)', sphess%hesslist(:, :, hat) - refhess%hesslist(:, :, kat)
               return
            end if
         end do
         if (.not. found) then
            call test_failed(error, "Pair of the short list missing in the full list")
            return
         end if
      end do
   end do

end subroutine test_hessian_pattern


!> Check the Hessian of an external function with a half neighbour list holding
!> all pairs against the dense evaluation. The second and the mixed derivatives
!> of the external function are given on shorter half lists of their own.
subroutine test_external_hessian(error, mol, model)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   class(mchrg_solver_type), allocatable :: solver
   class(mchrg_solver_input), allocatable :: solver_input
   type(mchrg_cache), allocatable :: cache
   type(csr_list), allocatable :: list, qlist, rlist
   type(mchrg_hessian), allocatable :: sphess
   integer :: ndim
   real(wp), allocatable :: dfdq(:), d2fdq2(:, :), dfdrdq(:, :, :), qval(:)
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), allocatable :: rval_ut(:, :), rval_lt(:, :), hessian(:, :), refhess(:, :)

   allocate(cg_input :: solver_input)
   select type (solver_input)
   type is (cg_input)
      solver_input%cgtol = 1.0e-15_wp
      solver_input%cgmiter = 1000
      solver_input%verbosity = 0
      solver_input%use_nlist = .true.
   end select
   call solver_maker(solver, solver_input, error)
   if (allocated(error)) return

   allocate(list, qlist, rlist)
   call new_csr_list(list, mol, error, cutoff=cutoff, complete=.false.)
   if (allocated(error)) return
   if (size(list%nlat) /= mol%nat * (mol%nat + 1) / 2) then
      call test_failed(error, "Neighbour list does not hold all pairs")
      return
   end if
   call new_csr_list(qlist, mol, error, cutoff=6.0_wp, complete=.false.)
   if (allocated(error)) return
   call new_csr_list(rlist, mol, error, cutoff=8.0_wp, complete=.false.)
   if (allocated(error)) return
   call get_external_derivs(mol, qlist, rlist, dfdq, d2fdq2, dfdrdq, qval, &
      & rval_ut, rval_lt)

   ndim = 3*mol%nat
   allocate(hessian(ndim, ndim), refhess(ndim, ndim))
   select type (model)
   type is (eeqbc_model)
      allocate(cache)
      call model%update(mol, cache, trans, grad=.false., hess=.true.)
      call model%solve(mol, solver, cache, error)
      if (allocated(error)) return
      call model%get_external_hessian(mol, solver, cache, error, dfdq, d2fdq2, &
         & dfdrdq, refhess)
      if (allocated(error)) return

      ! The list of the Hessian is a copy of the neighbour list
      deallocate(cache)
      allocate(cache, sphess)
      call model%update(mol, cache, trans, grad=.false., hess=.true., list=list)
      call model%solve(mol, solver, cache, error, list=list)
      if (allocated(error)) return
      call model%get_external_hessian(mol, solver, cache, error, dfdq, qval, &
         & qlist, rval_ut, rval_lt, rlist, sphess, list)
      if (allocated(error)) return
   class default
      call test_failed(error, "Analytical Hessian requires the EEQBC model")
      return
   end select
   call expand_hessian(mol, sphess, hessian)

   if (any(abs(hessian - refhess) > thr2)) then
      call test_failed(error, "External Hessian does not match the dense one")
      print'(a)', "diff:"
      print'(3es21.12)', hessian - refhess
   end if

end subroutine test_external_hessian


!> Check the Hessian of an external function with a half neighbour list holding
!> part of the pairs against its evaluation from the inverse, the partial
!> derivatives and the charge derivatives on the same list. Each atom x only
!> couples through its neighbours b, c and d,
!> dp_b/dR_x = sum_c S_bc (G_cx + sum_d F_cd dq_d/dR_x), and the mixed
!> derivatives G of the pairs outside of the list are dropped. Without df/dq the
!> adjoint vector vanishes and with it the explicit part, so the blocks only
!> hold the response to the derivatives of the external function.
subroutine test_external_hessian_sparse(error, mol, model, nlcutoff, ecutoff)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Cutoff of the neighbour list
   real(wp), intent(in) :: nlcutoff

   !> Cutoff of the lists of the derivatives of the external function
   real(wp), intent(in) :: ecutoff

   class(mchrg_solver_type), allocatable :: solver
   class(mchrg_solver_input), allocatable :: solver_input
   type(mchrg_cache), allocatable :: cache
   type(csr_list), allocatable :: list, qlist, rlist
   type(mchrg_hessian), allocatable :: sphess
   integer :: xat, yat, bat, cat, ic
   integer(i8) :: kat, hat
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp) :: blk(3, 3)
   logical, allocatable :: nbr(:, :)
   real(wp), allocatable :: dfdq(:), d2fdq2(:, :), dfdrdq(:, :, :), qval(:)
   real(wp), allocatable :: rval_ut(:, :), rval_lt(:, :), dqdL(:, :, :)
   real(wp), allocatable :: ainv(:, :), dabdr(:, :, :), dqdr(:, :, :)
   real(wp), allocatable :: mvec(:, :), pvec(:, :), fuvec(:), tvec(:), wvec(:, :)

   allocate(cg_input :: solver_input)
   select type (solver_input)
   type is (cg_input)
      solver_input%cgtol = 1.0e-15_wp
      solver_input%cgmiter = 1000
      solver_input%verbosity = 0
      solver_input%use_nlist = .true.
   end select
   call solver_maker(solver, solver_input, error)
   if (allocated(error)) return

   allocate(list, qlist, rlist)
   call new_csr_list(list, mol, error, cutoff=nlcutoff, complete=.false.)
   if (allocated(error)) return
   call new_csr_list(qlist, mol, error, cutoff=ecutoff, complete=.false.)
   if (allocated(error)) return
   call new_csr_list(rlist, mol, error, cutoff=ecutoff, complete=.false.)
   if (allocated(error)) return
   if (size(rlist%nlat) <= size(list%nlat)) then
      call test_failed(error, "Lists of the external derivatives are not longer")
      return
   end if
   call get_external_derivs(mol, qlist, rlist, dfdq, d2fdq2, dfdrdq, qval, &
      & rval_ut, rval_lt)
   dfdq(:) = 0.0_wp

   select type (model)
   type is (eeqbc_model)
      allocate(cache)
      call new_allpairs_hessian(mol, sphess, error)
      if (allocated(error)) return
      call model%update(mol, cache, trans, grad=.false., hess=.true., list=list)
      call model%solve(mol, solver, cache, error, list=list)
      if (allocated(error)) return
      call model%get_external_hessian(mol, solver, cache, error, dfdq, qval, &
         & qlist, rval_ut, rval_lt, rlist, sphess, list)
      if (allocated(error)) return
   class default
      call test_failed(error, "Analytical Hessian requires the EEQBC model")
      return
   end select

   ! Inverse, partial derivatives and charge derivatives without the charge
   ! constraint term on the list
   deallocate(cache)
   allocate(cache)
   allocate(dqdL(3, 3, mol%nat))
   call model%update(mol, cache, trans, grad=.false., list=list, dqdr=.true.)
   call model%solve(mol, solver, cache, error, dqdL=dqdL, list=list, &
      & unit=output_unit)
   if (allocated(error)) return

   allocate(nbr(mol%nat, mol%nat), source=.false.)
   allocate(ainv(mol%nat, mol%nat), source=0.0_wp)
   do xat = 1, mol%nat
      do kat = list%inl(xat), list%inl(xat + 1) - 1
         bat = list%nlat(kat)
         nbr(bat, xat) = .true.
         nbr(xat, bat) = .true.
         ainv(bat, xat) = cache%ainvlist(kat)
         ainv(xat, bat) = cache%ainvlist(kat)
      end do
   end do
   allocate(dabdr(3, mol%nat, mol%nat), dqdr(3, mol%nat, mol%nat))
   call expand_partial_derivs(mol, list, cache%dabdrlist_ut, cache%dabdrlist_lt, &
      & dabdr)
   call expand_partial_derivs(mol, list, cache%dqdrlist_ut, cache%dqdrlist_lt, &
      & dqdr)
   do xat = 1, mol%nat
      do bat = 1, mol%nat
         if (.not. nbr(bat, xat)) dfdrdq(:, xat, bat) = 0.0_wp
      end do
   end do

   ! Factor w = G^T u + B^T A^-1 F u of the charge constraint term
   allocate(fuvec(mol%nat), tvec(mol%nat), wvec(3, mol%nat))
   fuvec(:) = matmul(d2fdq2, cache%uvec)
   tvec(:) = 0.0_wp
   call solver%solve(alist=cache%alist, xvec=fuvec, vrhs=tvec, list=list, &
      & error=error)
   if (allocated(error)) return
   do xat = 1, mol%nat
      wvec(:, xat) = matmul(dfdrdq(:, xat, :), cache%uvec) &
         & + matmul(dabdr(:, xat, :), tvec)
   end do
   if (abs(sphess%weight - dot_product(cache%uvec, fuvec)) > thr1 &
      & .or. any(abs(sphess%wvec - wvec) > thr1)) then
      call test_failed(error, "Charge constraint term of the external Hessian "// &
         & "does not match")
      print'(a)', "diff:"
      print'(3es21.12)', sphess%wvec - wvec
      return
   end if

   ! dp_b/dR_x for the neighbours b of x and the blocks (x, y) of y >= x
   allocate(mvec(3, mol%nat), pvec(3, mol%nat))
   do xat = 1, mol%nat
      mvec(:, :) = 0.0_wp
      pvec(:, :) = 0.0_wp
      do cat = 1, mol%nat
         if (.not. nbr(cat, xat)) cycle
         mvec(:, cat) = dfdrdq(:, xat, cat)
         do bat = 1, mol%nat
            if (nbr(bat, xat)) mvec(:, cat) = mvec(:, cat) &
               & + d2fdq2(cat, bat) * dqdr(:, xat, bat)
         end do
      end do
      do bat = 1, mol%nat
         if (.not. nbr(bat, xat)) cycle
         do cat = 1, mol%nat
            if (nbr(cat, xat)) pvec(:, bat) = pvec(:, bat) &
               & + ainv(bat, cat) * mvec(:, cat)
         end do
      end do

      do hat = sphess%list%inl(xat), sphess%list%inl(xat + 1) - 1
         yat = sphess%list%nlat(hat)
         blk(:, :) = 0.0_wp
         do bat = 1, mol%nat
            if (.not. nbr(bat, xat)) cycle
            do ic = 1, 3
               blk(:, ic) = blk(:, ic) + dfdrdq(:, yat, bat) * dqdr(ic, xat, bat) &
                  & + dabdr(:, yat, bat) * pvec(ic, bat)
            end do
         end do
         if (yat == xat) blk(:, :) = 0.5_wp * (blk + transpose(blk))
         if (any(abs(sphess%hesslist(:, :, hat) - blk) > thr1)) then
            call test_failed(error, "External Hessian block does not match the "// &
               & "restriction to the neighbours")
            print'(a, 2i6)', "block", xat, yat
            print'(3es21.12)', sphess%hesslist(:, :, hat) - blk
            return
         end if
      end do
   end do

end subroutine test_external_hessian_sparse


!> Partial derivatives of an external function on half lists of their own and
!> in dense storage: df/dq, the second derivatives d2f/dq2 = 0.1 / (1 + r_ab)
!> on the pairs of qlist and the mixed derivatives d2f/dR_x dq_b on the pairs
!> of rlist, which are not symmetric in x and b
subroutine get_external_derivs(mol, qlist, rlist, dfdq, d2fdq2, dfdrdq, qval, &
   & rval_ut, rval_lt)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Half list of the second derivatives w.r.t. the charges
   type(csr_list), intent(in) :: qlist

   !> Half list of the mixed derivatives
   type(csr_list), intent(in) :: rlist

   !> Derivative w.r.t. the charges
   real(wp), allocatable, intent(out) :: dfdq(:)

   !> Second derivative w.r.t. the charges
   real(wp), allocatable, intent(out) :: d2fdq2(:, :)

   !> Mixed derivative, dfdrdq(:, x, b) = d2f/dR_x dq_b
   real(wp), allocatable, intent(out) :: dfdrdq(:, :, :)

   !> Upper triangle of the second derivative w.r.t. the charges on qlist
   real(wp), allocatable, intent(out) :: qval(:)

   !> Diagonal and upper triangle of the mixed derivative on rlist
   real(wp), allocatable, intent(out) :: rval_ut(:, :)

   !> Lower triangle of the mixed derivative on rlist
   real(wp), allocatable, intent(out) :: rval_lt(:, :)

   real(wp), parameter :: efield(3) = [0.1_wp, -0.2_wp, 0.3_wp]
   real(wp), parameter :: kpot = 0.05_wp
   integer :: iat, jat
   integer(i8) :: kat
   real(wp) :: vec(3), r2

   allocate(dfdq(mol%nat), d2fdq2(mol%nat, mol%nat), dfdrdq(3, mol%nat, mol%nat), &
      & qval(size(qlist%nlat)), rval_ut(3, size(rlist%nlat)), &
      & rval_lt(3, size(rlist%nlat)))

   do iat = 1, mol%nat
      dfdq(iat) = dot_product(efield, mol%xyz(:, iat)) &
         & + 0.5_wp * kpot * sum(mol%xyz(:, iat)**2)
   end do

   d2fdq2(:, :) = 0.0_wp
   do iat = 1, mol%nat
      do kat = qlist%inl(iat), qlist%inl(iat + 1) - 1
         jat = qlist%nlat(kat)
         qval(kat) = 0.1_wp / (1.0_wp + norm2(mol%xyz(:, jat) - mol%xyz(:, iat)))
         d2fdq2(jat, iat) = qval(kat)
         d2fdq2(iat, jat) = qval(kat)
      end do
   end do

   ! The entry of atom x in the row of atom b holds d2f/dR_x dq_b in rval_ut and
   ! d2f/dR_b dq_x in rval_lt
   do iat = 1, mol%nat
      do kat = rlist%inl(iat), rlist%inl(iat + 1) - 1
         jat = rlist%nlat(kat)
         vec(:) = mol%xyz(:, jat) - mol%xyz(:, iat)
         r2 = sum(vec**2)
         if (jat == iat) then
            rval_ut(:, kat) = efield + kpot * mol%xyz(:, iat)
            rval_lt(:, kat) = 0.0_wp
         else
            rval_ut(:, kat) = 0.02_wp * (1.0_wp + 0.1_wp * real(iat, wp)) * vec &
               & / (1.0_wp + r2)
            rval_lt(:, kat) = -0.03_wp * (1.0_wp + 0.1_wp * real(jat, wp)) * vec &
               & / (1.0_wp + r2)
         end if
      end do
   end do
   call expand_partial_derivs(mol, rlist, rval_ut, rval_lt, dfdrdq)

end subroutine get_external_derivs


!> Compressed Hessian with zero blocks on a half list holding all pairs
subroutine new_allpairs_hessian(mol, hessian, error)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Compressed Hessian
   type(mchrg_hessian), allocatable, intent(out) :: hessian

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   allocate(hessian)
   call new_csr_list(hessian%list, mol, error, cutoff=allpairs, complete=.false.)
   if (allocated(error)) return
   if (size(hessian%list%nlat) /= mol%nat * (mol%nat + 1) / 2) then
      call test_failed(error, "List of the Hessian does not hold all pairs")
      return
   end if
   allocate(hessian%hesslist(3, 3, size(hessian%list%nlat, kind=i8)), &
      & source=0.0_wp)

end subroutine new_allpairs_hessian


!> Expand the compressed Hessian and its factored charge constraint term, if
!> present, into dense storage
subroutine expand_hessian(mol, hessian, dense)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Compressed Hessian
   type(mchrg_hessian), intent(in) :: hessian

   !> Dense Hessian
   real(wp), intent(out) :: dense(:, :)

   integer :: xat, yat, ic, ii, jj
   integer(i8) :: hat

   dense(:, :) = 0.0_wp
   if (allocated(hessian%dlambda)) then
      do xat = 1, mol%nat
         ii = 3 * (xat - 1)
         do yat = 1, mol%nat
            jj = 3 * (yat - 1)
            do ic = 1, 3
               dense(jj+1:jj+3, ii+ic) = hessian%weight * hessian%dlambda(ic, xat) &
                  & * hessian%dlambda(:, yat)
               if (allocated(hessian%wvec)) dense(jj+1:jj+3, ii+ic) = &
                  & dense(jj+1:jj+3, ii+ic) &
                  & - hessian%wvec(ic, xat) * hessian%dlambda(:, yat) &
                  & - hessian%dlambda(ic, xat) * hessian%wvec(:, yat)
            end do
         end do
      end do
   end if
   do xat = 1, mol%nat
      ii = 3 * (xat - 1)
      do hat = hessian%list%inl(xat), hessian%list%inl(xat + 1) - 1
         jj = 3 * (hessian%list%nlat(hat) - 1)
         dense(jj+1:jj+3, ii+1:ii+3) = dense(jj+1:jj+3, ii+1:ii+3) &
            & + hessian%hesslist(:, :, hat)
         ! The half list keeps one block of each pair
         if (jj == ii) cycle
         dense(ii+1:ii+3, jj+1:jj+3) = dense(ii+1:ii+3, jj+1:jj+3) &
            & + transpose(hessian%hesslist(:, :, hat))
      end do
   end do

end subroutine expand_hessian


subroutine test_hexpl_coulomb_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "01")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_explicit_hessian(error, mol, model, alpha=1.0_wp, beta=0.0_wp, &
      & nlcutoff=cutoff)

end subroutine test_hexpl_coulomb_mb01


subroutine test_hexpl_en_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "01")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_explicit_hessian(error, mol, model, alpha=0.0_wp, beta=1.0_wp, &
      & nlcutoff=cutoff)

end subroutine test_hexpl_en_mb01


subroutine test_hexpl_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "05")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   ! Prefactors of the electrostatic energy
   call test_explicit_hessian(error, mol, model, alpha=0.5_wp, beta=-1.0_wp, &
      & nlcutoff=cutoff)

end subroutine test_hexpl_mb05


subroutine test_hexpl_dfdq_mb06(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   real(wp), allocatable :: pvec(:)

   integer :: iat

   call get_structure(mol, "MB16-43", "06")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   allocate(pvec(mol%nat))
   do iat = 1, mol%nat
      pvec(iat) = 0.1_wp * sum(mol%xyz(:, iat))
   end do
   ! Prefactors of an external function with the adjoint vector pvec
   call test_explicit_hessian(error, mol, model, alpha=-1.0_wp, beta=1.0_wp, &
      & nlcutoff=cutoff, pvec=pvec)

end subroutine test_hexpl_dfdq_mb06


subroutine test_hexpl_znooh(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   integer, parameter :: nat = 4
   real(wp), parameter :: charge = -1.0_wp
   integer, parameter :: num(nat) = [30, 8, 8, 1]
   real(wp), parameter :: xyz(3, nat) = reshape([ &
   & -0.30631629283878_wp, -1.11507514203552_wp, +0.00000000000000_wp, &
   & -0.06543072660074_wp, -4.32862093666082_wp, +0.00000000000000_wp, &
   & -0.64012239724097_wp, +2.34966763895920_wp, +0.00000000000000_wp, &
   & +1.01186941668051_wp, +3.09402843973713_wp, +0.00000000000000_wp],&
   & [3, nat])

   call new(mol, num, xyz, charge)
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_explicit_hessian(error, mol, model, alpha=0.5_wp, beta=-1.0_wp, &
      & nlcutoff=cutoff)

end subroutine test_hexpl_znooh


subroutine test_hexpl_sparse_upu23(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "UPU23", "0a")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   ! Short cutoff, the list only holds part of the pairs
   call test_explicit_hessian(error, mol, model, alpha=0.5_wp, beta=-1.0_wp, &
      & nlcutoff=12.0_wp)

end subroutine test_hexpl_sparse_upu23


subroutine test_hexpl_sparse_dfdq_upu23(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   real(wp), allocatable :: pvec(:)

   integer :: iat

   call get_structure(mol, "UPU23", "0a")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   allocate(pvec(mol%nat))
   do iat = 1, mol%nat
      pvec(iat) = 0.1_wp * sum(mol%xyz(:, iat))
   end do
   ! Short cutoff, the list only holds part of the pairs
   call test_explicit_hessian(error, mol, model, alpha=-1.0_wp, beta=1.0_wp, &
      & nlcutoff=12.0_wp, pvec=pvec)

end subroutine test_hexpl_sparse_dfdq_upu23


subroutine test_himpl_list_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "05")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_implicit_hessian(error, mol, model, nlcutoff=cutoff)

end subroutine test_himpl_list_mb05


subroutine test_himpl_sparse_list_upu23(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "UPU23", "0a")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   ! Short cutoff, the half list only holds part of the pairs
   call test_implicit_hessian(error, mol, model, nlcutoff=12.0_wp)

end subroutine test_himpl_sparse_list_upu23


subroutine test_hess_list_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "05")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_hessian(error, mol, model)

end subroutine test_hess_list_mb05


subroutine test_hess_pattern_list_upu23(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "UPU23", "0a")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   ! The half list of the Hessian is shorter than the half neighbour list
   call test_hessian_pattern(error, mol, model, nlcutoff=12.0_wp, hcutoff=8.0_wp)

end subroutine test_hess_pattern_list_upu23


subroutine test_hess_list_znooh(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_znooh(mol)
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_hessian(error, mol, model)

end subroutine test_hess_list_znooh


subroutine test_hext_list_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "05")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_external_hessian(error, mol, model)

end subroutine test_hext_list_mb05


subroutine test_hext_list_znooh(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_znooh(mol)
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_external_hessian(error, mol, model)

end subroutine test_hext_list_znooh


subroutine test_hext_sparse_list_upu23(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "UPU23", "0a")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   ! The lists of the external derivatives are longer than the neighbour list
   call test_external_hessian_sparse(error, mol, model, nlcutoff=12.0_wp, &
      & ecutoff=16.0_wp)

end subroutine test_hext_sparse_list_upu23


!> Charged ZnOOH molecule
subroutine get_znooh(mol)

   !> Molecular structure data
   type(structure_type), intent(out) :: mol

   integer, parameter :: nat = 4
   real(wp), parameter :: charge = -1.0_wp
   integer, parameter :: num(nat) = [30, 8, 8, 1]
   real(wp), parameter :: xyz(3, nat) = reshape([ &
   & -0.30631629283878_wp, -1.11507514203552_wp, +0.00000000000000_wp, &
   & -0.06543072660074_wp, -4.32862093666082_wp, +0.00000000000000_wp, &
   & -0.64012239724097_wp, +2.34966763895920_wp, +0.00000000000000_wp, &
   & +1.01186941668051_wp, +3.09402843973713_wp, +0.00000000000000_wp],&
   & [3, nat])

   call new(mol, num, xyz, charge)

end subroutine get_znooh


end module test_csrlist
