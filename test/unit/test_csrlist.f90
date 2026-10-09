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
   use multicharge_model_type, only : mchrg_model_type
   use multicharge_model_eeqbc, only : eeqbc_model
   use multicharge_param, only : new_eeq2019_model, new_eeqbc2025_model
   use multicharge_model_cache, only : mchrg_cache
   use multicharge_ncoord, only : get_pair_derivs
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
   & new_unittest("pair-derivs-mb01", test_pair_derivs_mb01), &
   & new_unittest("pair-derivs-charged", test_pair_derivs_charged), &
   & new_unittest("pair-derivs-co2", test_pair_derivs_co2), &
   & new_unittest("pair-derivs-ice", test_pair_derivs_ice), &
   & new_unittest("partial-derivs-mb01", test_partial_derivs_mb01), &
   & new_unittest("partial-derivs-co2", test_partial_derivs_co2), &
   & new_unittest("partial-derivs-ice", test_partial_derivs_ice) &
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


!> Rebuild dense CN derivatives from the pair derivatives and compare them
subroutine check_pair_derivs(error, mol, ncoord, trans, cn, dcndr, dcndL)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Coordination number container
   class(ncoord_type), intent(in) :: ncoord

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Coordination numbers including the maximum CN cutoff
   real(wp), intent(in) :: cn(:)

   !> Reference derivative w.r.t. the Cartesian coordinates
   real(wp), intent(in) :: dcndr(:, :, :)

   !> Reference derivative w.r.t. strain deformations
   real(wp), intent(in) :: dcndL(:, :, :)

   integer :: iat, jat, itr
   real(wp) :: vec(3), dpair(3, 2), spair(3, 3, 2)
   real(wp), allocatable :: dcndr_pair(:, :, :), dcndL_pair(:, :, :)

   allocate(dcndr_pair(3, mol%nat, mol%nat), dcndL_pair(3, 3, mol%nat), &
      & source=0.0_wp)

   do iat = 1, mol%nat
      do jat = 1, iat
         do itr = 1, size(trans, 2)
            vec = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, itr))
            call get_pair_derivs(ncoord, mol, iat, jat, vec, cn(iat), cn(jat), &
               & dpair, spair)

            dcndr_pair(:, iat, iat) = dcndr_pair(:, iat, iat) + dpair(:, 1)
            dcndr_pair(:, jat, iat) = dcndr_pair(:, jat, iat) - dpair(:, 1)
            dcndr_pair(:, iat, jat) = dcndr_pair(:, iat, jat) + dpair(:, 2)
            dcndr_pair(:, jat, jat) = dcndr_pair(:, jat, jat) - dpair(:, 2)

            dcndL_pair(:, :, iat) = dcndL_pair(:, :, iat) + spair(:, :, 1)
            dcndL_pair(:, :, jat) = dcndL_pair(:, :, jat) + spair(:, :, 2)
         end do
      end do
   end do

   if (any(abs(dcndr_pair - dcndr) > thr1)) then
      call test_failed(error, "Pair derivatives w.r.t. positions do not match")
      print'(es21.14)', maxval(abs(dcndr_pair - dcndr))
      return
   end if

   if (any(abs(dcndL_pair - dcndL) > thr1)) then
      call test_failed(error, "Pair derivatives w.r.t. strain do not match")
      print'(es21.14)', maxval(abs(dcndL_pair - dcndL))
      return
   end if

end subroutine check_pair_derivs


!> Check pair derivatives of the CN and the local charges of a model
subroutine test_pair_derivs(error, mol, model)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   real(wp), allocatable :: trans(:, :)
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: qloc(:), dqlocdr(:, :, :), dqlocdL(:, :, :)

   allocate(cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat))
   allocate(qloc(mol%nat), dqlocdr(3, mol%nat, mol%nat), dqlocdL(3, 3, mol%nat))

   call get_lattice_points(mol%periodic, mol%lattice, model%ncoord%cutoff, trans)

   call model%ncoord%get_coordination_number(mol, trans, cn, dcndr=dcndr, &
      & dcndL=dcndL)
   call check_pair_derivs(error, mol, model%ncoord, trans, cn, dcndr, dcndL)
   if (allocated(error)) return

   ! The maximum CN cutoff acts on the local charges before the total charge
   ! is distributed
   call model%local_charge(mol, trans, qloc, dqlocdr=dqlocdr, dqlocdL=dqlocdL)
   call check_pair_derivs(error, mol, model%ncoord_en, trans, &
      & qloc - mol%charge / real(mol%nat, wp), dqlocdr, dqlocdL)

end subroutine test_pair_derivs


subroutine test_pair_derivs_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "01")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_pair_derivs(error, mol, model)

end subroutine test_pair_derivs_mb01


subroutine test_pair_derivs_charged(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "02")
   mol%charge = -2.0_wp
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_pair_derivs(error, mol, model)

end subroutine test_pair_derivs_charged


subroutine test_pair_derivs_co2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "X23", "CO2")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_pair_derivs(error, mol, model)

end subroutine test_pair_derivs_co2


subroutine test_pair_derivs_ice(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "ICE10", "vi")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_pair_derivs(error, mol, model)

end subroutine test_pair_derivs_ice


!> Expand the complete neighbour list storage of the partial derivatives
subroutine expand_partial_derivs(mol, list, cache, dabdr)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Complete CSR neighbour list
   type(csr_list), intent(in) :: list

   !> Cache holding the partial derivatives in neighbour list storage
   type(mchrg_cache), intent(in) :: cache

   !> Dense partial derivatives
   real(wp), intent(out) :: dabdr(:, :, :)

   integer :: iat
   integer(i8) :: kat

   dabdr(:, :, :) = 0.0_wp
   do iat = 1, mol%nat
      do kat = list%inl(iat), list%inl(iat + 1) - 1
         dabdr(:, list%nlat(kat), iat) = cache%dabdrlist(:, kat)
      end do
   end do

end subroutine expand_partial_derivs


!> Compare the partial derivatives and charge derivatives obtained with the
!> complete neighbour list against the dense evaluation with the block CG and
!> the direct solver
subroutine test_partial_derivs(error, mol, model)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   class(mchrg_solver_type), allocatable :: solver, direct
   class(mchrg_solver_input), allocatable :: solver_input, direct_in
   type(mchrg_cache), allocatable :: cache1, cache2, cache3
   type(csr_list), allocatable :: list
   integer :: ndim, iscale, jat
   real(wp), parameter :: scales(2, 2) = reshape([1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & [2, 2])
   real(wp), allocatable :: trans(:, :), dabdr(:, :, :)
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
   call model%update(mol, cache1, trans, grad=.true.)
   call model%solve(mol, solver, cache1, error, dqdr=dqdr1, dqdL=dqdL1, &
      & unit=output_unit)
   if (allocated(error)) return

   ! Dense reference with the direct solver
   allocate(cache3)
   call model%update(mol, cache3, trans, grad=.true.)
   call model%solve(mol, direct, cache3, error, dqdr=dqdr3, dqdL=dqdL3, &
      & unit=output_unit)
   if (allocated(error)) return

   ! Complete neighbour list, inverse on the pattern of the list
   allocate(cache2, list)
   if (any(mol%periodic)) then
      call new_csr_list(list, mol, error, cache2%wsc, cutoff=cutoff, complete=.true.)
   else
      call new_csr_list(list, mol, error, cutoff=cutoff, complete=.true.)
   end if
   if (allocated(error)) return
   call model%update(mol, cache2, trans, grad=.true., list=list)
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

   ! Position derivatives only in compressed storage, without a dense dqdr
   call model%solve(mol, solver, cache2, error, dqdL=dqdL4, list=list, &
      & unit=output_unit)
   if (allocated(error)) return

   do jat = 1, mol%nat
      call cache2%get_dqdr_row(jat, dqrow)
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

   ! Coulomb matrix and electronegativity parts separately, same charges
   cache2%vrhs(:) = cache1%vrhs
   do iscale = 1, 2
      call model%get_partial_derivs(mol, ndim, cache1, scales(1, iscale), &
         & scales(2, iscale))
      call model%get_partial_derivs(mol, ndim, cache2, scales(1, iscale), &
         & scales(2, iscale), list=list)
      call expand_partial_derivs(mol, list, cache2, dabdr)

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


end module test_csrlist
