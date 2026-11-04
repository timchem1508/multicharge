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

!> Unit tests for the CSR neighbour list representation of the charge models,
!> validated against the dense representation of the same model
module test_csrlist
   use mctc_csrlist, only : csr_list, new_csr_list
   use mctc_cutoff, only : get_lattice_points
   use mctc_env, only : wp, i8
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, test_failed
   use mctc_io_structure, only : structure_type, new
   use mstore, only : get_structure
   use multicharge_charge, only : get_charges, get_eeqbc_charges
   use multicharge_model_cache, only : mchrg_cache
   use multicharge_model_type, only : mchrg_model_type
   use multicharge_param, only : new_eeq2019_model, new_eeqbc2025_model
   use multicharge_solver_cg, only : cg_solver, new_cg_solver, cg_input
   use multicharge_solver_direct, only : direct_solver, new_direct_solver, &
      & direct_input
   implicit none
   private

   public :: collect_csrlist

   !> Tolerance for the matrices and vectors of the linear system
   real(wp), parameter :: thr = 1.0e5_wp*epsilon(1.0_wp)

   !> Tolerance for the solution and its derivatives
   real(wp), parameter :: thr1 = 1.0e-10_wp

   !> Neighbour list cutoff
   real(wp), parameter :: cutoff = 29.0_wp


contains


!> Collect all exported unit tests
subroutine collect_csrlist(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("eeqbc-components-mb01", test_eeqbc_components_mb01), &
      & new_unittest("eeqbc-components-co2", test_eeqbc_components_co2), &
      & new_unittest("eeqbc-energy-mb03", test_eeqbc_e_mb03), &
      & new_unittest("eeqbc-energy-mb04", test_eeqbc_e_mb04), &
      & new_unittest("eeqbc-energy-actinides", test_eeqbc_e_actinides), &
      & new_unittest("eeqbc-energy-co2", test_eeqbc_e_co2), &
      & new_unittest("eeqbc-energy-ice", test_eeqbc_e_ice), &
      & new_unittest("eeqbc-energy-ice-supercell", test_eeqbc_e_ice222), &
      & new_unittest("eeqbc-gradient-mb05", test_eeqbc_g_mb05), &
      & new_unittest("eeqbc-gradient-mb06", test_eeqbc_g_mb06), &
      & new_unittest("eeqbc-gradient-co2", test_eeqbc_g_co2), &
      & new_unittest("eeqbc-gradient-ice", test_eeqbc_g_ice), &
      & new_unittest("eeqbc-external-gradient-mb07", test_eeqbc_eg_mb07), &
      & new_unittest("eeqbc-external-gradient-co2", test_eeqbc_eg_co2), &
      & new_unittest("eeqbc-api-mb01", test_eeqbc_api_mb01), &
      & new_unittest("eeqbc-api-actinides", test_eeqbc_api_actinides), &
      & new_unittest("eeq-list-unsupported", test_eeq_list_error, &
      & should_fail=.true.), &
      & new_unittest("direct-list-unsupported", test_direct_list_error, &
      & should_fail=.true.), &
      & new_unittest("complete-list-unsupported", test_complete_list_error, &
      & should_fail=.true.), &
      & new_unittest("periodic-list-without-wsc", test_periodic_wsc_error, &
      & should_fail=.true.), &
      & new_unittest("api-list-charge-derivatives", test_api_dqdr_error, &
      & should_fail=.true.) &
      & ]

end subroutine collect_csrlist


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
      & [rep(1)*mol%lattice(:, 1), rep(2)*mol%lattice(:, 2), rep(3)*mol%lattice(:, 3)], &
      & shape(mol%lattice))
   allocate(xyz(3, product(rep)*mol%nat))
   c = 0
   do i = 0, rep(1) - 1
      do j = 0, rep(2) - 1
         do k = 0, rep(3) - 1
            xyz(:, c+1:c+mol%nat) = mol%xyz &
               & + spread(matmul(mol%lattice, [real(wp) :: i, j, k]), 2, mol%nat)
            c = c + mol%nat
         end do
      end do
   end do

   call new(mol, num, xyz, lattice=lattice)

end subroutine make_supercell


!> Build the neighbour list, periodic images are resolved with the
!> Wigner-Seitz cell stored in the cache
subroutine build_list(error, mol, cache, list)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Cache receiving the Wigner-Seitz cell of the neighbour list
   type(mchrg_cache), intent(inout) :: cache

   !> CSR neighbour list
   type(csr_list), intent(out) :: list

   if (any(mol%periodic)) then
      call new_csr_list(list, mol, error, cache%wsc, cutoff)
   else
      call new_csr_list(list, mol, error, cutoff=cutoff)
   end if

end subroutine build_list


!> Update the cache and solve the model with the CG solver
subroutine solve_model(error, mol, model, cache, energy, qvec, gradient, sigma, list)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Charge model
   class(mchrg_model_type), intent(in) :: model

   !> Cache of the model
   type(mchrg_cache), intent(inout) :: cache

   !> Atom-resolved energy
   real(wp), intent(inout), optional :: energy(:)

   !> Atomic partial charges
   real(wp), intent(out), optional :: qvec(:)

   !> Energy gradient
   real(wp), intent(inout), optional :: gradient(:, :)

   !> Strain derivative
   real(wp), intent(inout), optional :: sigma(:, :)

   !> Optional CSR neighbour list
   type(csr_list), intent(in), optional :: list

   type(cg_input) :: input
   type(cg_solver) :: solver
   real(wp), allocatable :: trans(:, :)

   call new_cg_solver(solver, input)
   call get_lattice_points(mol%periodic, mol%lattice, model%ncoord%cutoff, trans)
   call model%update(mol, cache, trans, present(gradient), list)
   call model%solve(mol, solver, cache, error, energy=energy, qvec=qvec, &
      & gradient=gradient, sigma=sigma, list=list)

end subroutine solve_model


!> Compare two arrays and report the largest deviation
subroutine check_array(error, label, actual, expected, tol)

   !> Error handling
   type(error_type), allocatable, intent(inout) :: error

   !> Name of the compared quantity
   character(len=*), intent(in) :: label

   !> Neighbour list result
   real(wp), intent(in) :: actual(:)

   !> Dense reference result
   real(wp), intent(in) :: expected(:)

   !> Tolerance of the comparison
   real(wp), intent(in) :: tol

   if (allocated(error)) return
   if (any(abs(actual - expected) > tol)) then
      call test_failed(error, label//" does not match the dense representation")
      print '(a, es12.4)', "Largest deviation:", maxval(abs(actual - expected))
   end if

end subroutine check_array


!> Expand a matrix in the CSR list layout to its dense symmetric form
subroutine expand_list(nat, list, mlist, mat)

   !> Number of atoms
   integer, intent(in) :: nat

   !> CSR neighbour list
   type(csr_list), intent(in) :: list

   !> Matrix in the CSR list layout
   real(wp), intent(in) :: mlist(:)

   !> Dense matrix
   real(wp), allocatable, intent(out) :: mat(:, :)

   integer :: iat, jat
   integer(i8) :: kat

   allocate(mat(nat, nat), source=0.0_wp)
   do iat = 1, nat
      do kat = list%inl(iat), list%inl(iat + 1) - 1
         jat = list%nlat(kat)
         mat(jat, iat) = mlist(kat)
         mat(iat, jat) = mlist(kat)
      end do
   end do

end subroutine expand_list


!> Compare the capacitance matrix, Coulomb matrix and right-hand side
subroutine gen_test_components(error, mol)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   class(mchrg_model_type), allocatable :: model
   type(mchrg_cache) :: dense, sparse
   type(csr_list) :: list
   real(wp), allocatable :: cmat(:, :), amat(:, :)

   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return

   call solve_model(error, mol, model, dense)
   if (allocated(error)) return

   call build_list(error, mol, sparse, list)
   if (allocated(error)) return
   call solve_model(error, mol, model, sparse, list=list)
   if (allocated(error)) return

   call expand_list(mol%nat, list, sparse%clist, cmat)
   call expand_list(mol%nat, list, sparse%alist, amat)

   call check_array(error, "Capacitance matrix", reshape(cmat, [size(cmat)]), &
      & reshape(dense%cmat, [size(dense%cmat)]), thr)
   call check_array(error, "Coulomb matrix", reshape(amat, [size(amat)]), &
      & reshape(dense%amat, [size(dense%amat)]), thr)
   call check_array(error, "Electronegativity vector", sparse%xvec, dense%xvec, thr)

end subroutine gen_test_components


!> Compare partial charges and energies
subroutine gen_test_energy(error, mol)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   class(mchrg_model_type), allocatable :: model
   type(mchrg_cache) :: dense, sparse
   type(csr_list) :: list
   real(wp), allocatable :: eref(:), energy(:), qref(:), qvec(:)

   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return

   allocate(eref(mol%nat), energy(mol%nat), source=0.0_wp)
   allocate(qref(mol%nat), qvec(mol%nat))

   call solve_model(error, mol, model, dense, energy=eref, qvec=qref)
   if (allocated(error)) return

   call build_list(error, mol, sparse, list)
   if (allocated(error)) return
   call solve_model(error, mol, model, sparse, energy=energy, qvec=qvec, list=list)
   if (allocated(error)) return

   call check_array(error, "Partial charges", qvec, qref, thr1)
   call check_array(error, "Energies", energy, eref, thr1)

end subroutine gen_test_energy


!> Compare the energy gradient and strain derivative
subroutine gen_test_gradient(error, mol)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   class(mchrg_model_type), allocatable :: model
   type(mchrg_cache) :: dense, sparse
   type(csr_list) :: list
   real(wp), allocatable :: energy(:), gradient(:, :), sigma(:, :)
   real(wp), allocatable :: gref(:, :), sref(:, :)

   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return

   allocate(energy(mol%nat), source=0.0_wp)
   allocate(gref(3, mol%nat), gradient(3, mol%nat), sref(3, 3), sigma(3, 3), &
      & source=0.0_wp)

   call solve_model(error, mol, model, dense, energy=energy, gradient=gref, sigma=sref)
   if (allocated(error)) return

   call build_list(error, mol, sparse, list)
   if (allocated(error)) return
   call solve_model(error, mol, model, sparse, energy=energy, gradient=gradient, &
      & sigma=sigma, list=list)
   if (allocated(error)) return

   call check_array(error, "Energy gradient", reshape(gradient, [size(gradient)]), &
      & reshape(gref, [size(gref)]), thr1)
   call check_array(error, "Strain derivative", reshape(sigma, [size(sigma)]), &
      & reshape(sref, [size(sref)]), thr1)

end subroutine gen_test_gradient


!> Compare the adjoint gradient of an external function of the charges
subroutine gen_test_external_gradient(error, mol)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   class(mchrg_model_type), allocatable :: model
   type(mchrg_cache) :: dense, sparse
   type(csr_list) :: list
   type(cg_input) :: input
   type(cg_solver) :: solver
   integer :: iat
   real(wp), allocatable :: dfdq(:), qvec(:), gradient(:, :), sigma(:, :)
   real(wp), allocatable :: gref(:, :), sref(:, :)

   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call new_cg_solver(solver, input)

   allocate(qvec(mol%nat), dfdq(mol%nat))
   allocate(gref(3, mol%nat), gradient(3, mol%nat), sref(3, 3), sigma(3, 3), &
      & source=0.0_wp)
   dfdq(:) = [(0.1_wp * real(modulo(3 * iat, 7) - 3, wp), iat = 1, mol%nat)]

   call solve_model(error, mol, model, dense, qvec=qvec, gradient=gref, sigma=sref)
   if (allocated(error)) return
   gref(:, :) = 0.0_wp
   sref(:, :) = 0.0_wp
   call model%get_external_gradient(mol, solver, dense, error, dfdq, gref, sref)
   if (allocated(error)) return

   call build_list(error, mol, sparse, list)
   if (allocated(error)) return
   call solve_model(error, mol, model, sparse, qvec=qvec, gradient=gradient, &
      & sigma=sigma, list=list)
   if (allocated(error)) return
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp
   call model%get_external_gradient(mol, solver, sparse, error, dfdq, gradient, &
      & sigma, list=list)
   if (allocated(error)) return

   call check_array(error, "External gradient", reshape(gradient, [size(gradient)]), &
      & reshape(gref, [size(gref)]), thr1)
   call check_array(error, "External strain derivative", reshape(sigma, [size(sigma)]), &
      & reshape(sref, [size(sref)]), thr1)

end subroutine gen_test_external_gradient


!> Compare the charges of the high-level interface
subroutine gen_test_api(error, mol)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   class(mchrg_model_type), allocatable :: model
   type(csr_list) :: list
   real(wp), allocatable :: qref(:), qvec(:)

   allocate(qref(mol%nat), qvec(mol%nat))
   call new_csr_list(list, mol, error, cutoff=cutoff)
   if (allocated(error)) return

   call get_eeqbc_charges(mol, error, qref)
   if (allocated(error)) return

   call get_eeqbc_charges(mol, error, qvec, list=list)
   if (allocated(error)) return
   call check_array(error, "Partial charges", qvec, qref, thr1)
   if (allocated(error)) return

   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   qvec(:) = 0.0_wp
   call get_charges(model, mol, error, qvec, list=list)
   if (allocated(error)) return
   call check_array(error, "Partial charges", qvec, qref, thr1)

end subroutine gen_test_api


subroutine test_eeqbc_components_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "01")
   call gen_test_components(error, mol)

end subroutine test_eeqbc_components_mb01


subroutine test_eeqbc_components_co2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "CO2")
   call gen_test_components(error, mol)

end subroutine test_eeqbc_components_co2


subroutine test_eeqbc_e_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "03")
   call gen_test_energy(error, mol)

end subroutine test_eeqbc_e_mb03


subroutine test_eeqbc_e_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "04")
   call gen_test_energy(error, mol)

end subroutine test_eeqbc_e_mb04


subroutine test_eeqbc_e_actinides(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "f-block", "Fr_to_Lr")
   call gen_test_energy(error, mol)

end subroutine test_eeqbc_e_actinides


subroutine test_eeqbc_e_co2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "CO2")
   call gen_test_energy(error, mol)

end subroutine test_eeqbc_e_co2


subroutine test_eeqbc_e_ice(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "ICE10", "vi")
   call gen_test_energy(error, mol)

end subroutine test_eeqbc_e_ice


subroutine test_eeqbc_e_ice222(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "ICE10", "vi")
   call make_supercell(mol, [2, 2, 2])
   call gen_test_energy(error, mol)

end subroutine test_eeqbc_e_ice222


subroutine test_eeqbc_g_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "05")
   call gen_test_gradient(error, mol)

end subroutine test_eeqbc_g_mb05


subroutine test_eeqbc_g_mb06(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "06")
   call gen_test_gradient(error, mol)

end subroutine test_eeqbc_g_mb06


subroutine test_eeqbc_g_co2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "CO2")
   call gen_test_gradient(error, mol)

end subroutine test_eeqbc_g_co2


subroutine test_eeqbc_g_ice(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "ICE10", "vi")
   call gen_test_gradient(error, mol)

end subroutine test_eeqbc_g_ice


subroutine test_eeqbc_eg_mb07(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "07")
   call gen_test_external_gradient(error, mol)

end subroutine test_eeqbc_eg_mb07


subroutine test_eeqbc_eg_co2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "CO2")
   call gen_test_external_gradient(error, mol)

end subroutine test_eeqbc_eg_co2


subroutine test_eeqbc_api_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "01")
   call gen_test_api(error, mol)

end subroutine test_eeqbc_api_mb01


subroutine test_eeqbc_api_actinides(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "f-block", "Fr_to_Lr")
   call gen_test_api(error, mol)

end subroutine test_eeqbc_api_actinides


!> The long-range EEQ model does not support the neighbour list
subroutine test_eeq_list_error(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   type(mchrg_cache) :: cache
   type(csr_list) :: list
   real(wp), allocatable :: qvec(:)

   call get_structure(mol, "MB16-43", "01")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   allocate(qvec(mol%nat))

   call build_list(error, mol, cache, list)
   if (allocated(error)) return
   call solve_model(error, mol, model, cache, qvec=qvec, list=list)

end subroutine test_eeq_list_error


!> The direct solver requires the dense representation
subroutine test_direct_list_error(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   type(mchrg_cache) :: cache
   type(csr_list) :: list
   type(direct_input) :: input
   type(direct_solver) :: solver
   real(wp), allocatable :: qvec(:), trans(:, :)

   call get_structure(mol, "MB16-43", "01")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   allocate(qvec(mol%nat))
   call new_direct_solver(solver, input)

   call build_list(error, mol, cache, list)
   if (allocated(error)) return
   call get_lattice_points(mol%periodic, mol%lattice, model%ncoord%cutoff, trans)
   call model%update(mol, cache, trans, .false., list)
   call model%solve(mol, solver, cache, error, qvec=qvec, list=list)

end subroutine test_direct_list_error


!> The model routines require the upper triangular neighbour list
subroutine test_complete_list_error(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   type(mchrg_cache) :: cache
   type(csr_list) :: list
   real(wp), allocatable :: qvec(:)

   call get_structure(mol, "MB16-43", "01")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   allocate(qvec(mol%nat))

   call new_csr_list(list, mol, error, cutoff=cutoff, complete=.true.)
   if (allocated(error)) return
   call solve_model(error, mol, model, cache, qvec=qvec, list=list)

end subroutine test_complete_list_error


!> Periodic neighbour lists must provide the Wigner-Seitz images of the cache
subroutine test_periodic_wsc_error(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   type(mchrg_cache) :: cache
   type(csr_list) :: list
   real(wp), allocatable :: qvec(:)

   call get_structure(mol, "X23", "CO2")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   allocate(qvec(mol%nat))

   call new_csr_list(list, mol, error, cutoff=cutoff)
   if (allocated(error)) return
   call solve_model(error, mol, model, cache, qvec=qvec, list=list)

end subroutine test_periodic_wsc_error


!> Charge derivatives are not available from the neighbour list interface
subroutine test_api_dqdr_error(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(csr_list) :: list
   real(wp), allocatable :: qvec(:), dqdr(:, :, :), dqdL(:, :, :)

   call get_structure(mol, "MB16-43", "01")
   allocate(qvec(mol%nat), dqdr(3, mol%nat, mol%nat), dqdL(3, 3, mol%nat))

   call new_csr_list(list, mol, error, cutoff=cutoff)
   if (allocated(error)) return
   call get_eeqbc_charges(mol, error, qvec, dqdr, dqdL, list=list)

end subroutine test_api_dqdr_error


end module test_csrlist
