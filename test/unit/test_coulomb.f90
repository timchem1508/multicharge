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

!> Unit tests for the particle mesh Ewald operator of the EEQ model and the
!> projected CG solver, validated against the dense Coulomb matrix
module test_coulomb
   use mctc_csrlist, only : csr_list, new_csr_list
   use mctc_cutoff, only : get_lattice_points
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, test_failed
   use mctc_io_structure, only : structure_type, new
   use mstore, only : get_structure
   use multicharge_coulomb, only : coulomb_input, coulomb_method
   use multicharge_model_cache, only : mchrg_cache
   use multicharge_model_type, only : mchrg_model_type
   use multicharge_param, only : new_eeq2019_model, new_eeqbc2025_model
   use multicharge_solver, only : mchrg_solver_type, new_mchrg_solver, cg_input, &
      & direct_input
   implicit none
   private

   public :: collect_coulomb

   !> Tolerance for results of a tightly converged operator
   real(wp), parameter :: thr = 1.0e-7_wp

   !> Tolerance for results of the operator with the default accuracy
   real(wp), parameter :: thr2 = 1.0e-5_wp

   !> Accuracy of the tightly converged particle mesh Ewald operator
   real(wp), parameter :: tight = 1.0e-10_wp

   !> Convergence threshold of the CG solver
   real(wp), parameter :: cgtol = 1.0e-12_wp


contains


!> Collect all exported unit tests
subroutine collect_coulomb(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("pme-mb01", test_pme_mb01), &
      & new_unittest("pme-urea", test_pme_urea), &
      & new_unittest("pme-charged-co2", test_pme_co2_charged), &
      & new_unittest("pme-amylose", test_pme_amylose), &
      & new_unittest("pme-urea-supercell", test_pme_urea222), &
      & new_unittest("pme-ice-hexagonal", test_pme_ice), &
      & new_unittest("pme-default-mb05", test_pme_default_mb05), &
      & new_unittest("projected-dense-mb02", test_projected_mb02), &
      & new_unittest("projected-dense-urea", test_projected_urea), &
      & new_unittest("projected-eeqbc-list", test_projected_eeqbc), &
      & new_unittest("pme-external-gradient", test_pme_external_gradient), &
      & new_unittest("operator-direct-solver", test_operator_direct, &
      & should_fail=.true.), &
      & new_unittest("operator-invalid-accuracy", test_invalid_accuracy, &
      & should_fail=.true.) &
      & ]

end subroutine collect_coulomb


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


!> Create a direct or a CG solver
subroutine get_solver(error, solver, iterative, projected)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Solver instance
   class(mchrg_solver_type), allocatable, intent(out) :: solver

   !> Use the iterative CG solver
   logical, intent(in) :: iterative

   !> Enforce the total charge by projection in the CG solver
   logical, intent(in) :: projected

   type(cg_input) :: cg
   type(direct_input) :: direct

   if (iterative) then
      cg%cgtol = cgtol
      cg%projected = projected
      call new_mchrg_solver(solver, cg, error)
   else
      call new_mchrg_solver(solver, direct, error)
   end if

end subroutine get_solver


!> Update the cache and solve the model
subroutine solve_model(error, mol, model, solver, energy, qvec, gradient, sigma, list)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Charge model
   class(mchrg_model_type), intent(in) :: model

   !> Solver instance
   class(mchrg_solver_type), intent(in) :: solver

   !> Atom-resolved energy
   real(wp), intent(inout), contiguous :: energy(:)

   !> Atomic partial charges
   real(wp), intent(out), contiguous :: qvec(:)

   !> Energy gradient
   real(wp), intent(inout), contiguous, optional :: gradient(:, :)

   !> Strain derivative
   real(wp), intent(inout), contiguous, optional :: sigma(:, :)

   !> Optional CSR neighbour list
   type(csr_list), intent(in), optional :: list

   type(mchrg_cache) :: cache
   real(wp), allocatable :: trans(:, :)

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

   !> Result of the tested representation
   real(wp), intent(in) :: actual(:)

   !> Reference result
   real(wp), intent(in) :: expected(:)

   !> Tolerance of the comparison
   real(wp), intent(in) :: tol

   if (allocated(error)) return
   if (any(abs(actual - expected) > tol)) then
      call test_failed(error, label//" does not match the reference")
      print "(a, es12.4)", "Largest deviation:", maxval(abs(actual - expected))
   end if

end subroutine check_array


!> Compare charges, energies and optionally the derivatives of a matrix-free
!> operator with the dense Coulomb matrix and the direct solver
subroutine gen_test_operator(error, mol, input, tol, grad)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Representation of the Coulomb matrix
   type(coulomb_input), intent(in) :: input

   !> Tolerance of the comparison
   real(wp), intent(in) :: tol

   !> Compare the gradient and the strain derivative
   logical, intent(in) :: grad

   class(mchrg_model_type), allocatable :: model
   class(mchrg_solver_type), allocatable :: solver
   real(wp), allocatable :: eref(:), energy(:), qref(:), qvec(:)
   real(wp), allocatable :: gref(:, :), gradient(:, :), sref(:, :), sigma(:, :)

   allocate(eref(mol%nat), energy(mol%nat), source=0.0_wp)
   allocate(qref(mol%nat), qvec(mol%nat))
   allocate(gref(3, mol%nat), gradient(3, mol%nat), sref(3, 3), sigma(3, 3), &
      & source=0.0_wp)

   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call get_solver(error, solver, .false., .false.)
   if (allocated(error)) return
   if (grad) then
      call solve_model(error, mol, model, solver, eref, qref, gref, sref)
   else
      call solve_model(error, mol, model, solver, eref, qref)
   end if
   if (allocated(error)) return

   call new_eeq2019_model(mol, model, error, coulomb=input)
   if (allocated(error)) return
   call get_solver(error, solver, .true., .true.)
   if (allocated(error)) return
   if (grad) then
      call solve_model(error, mol, model, solver, energy, qvec, gradient, sigma)
   else
      call solve_model(error, mol, model, solver, energy, qvec)
   end if
   if (allocated(error)) return

   call check_array(error, "Partial charges", qvec, qref, tol)
   call check_array(error, "Total energy", [sum(energy)], [sum(eref)], tol * mol%nat)
   if (grad) then
      call check_array(error, "Gradient", reshape(gradient, [size(gradient)]), &
         & reshape(gref, [size(gref)]), tol)
      call check_array(error, "Strain derivative", reshape(sigma, [9]), &
         & reshape(sref, [9]), tol * mol%nat)
   end if

end subroutine gen_test_operator


!> Compare the projected CG solver with the direct solver
subroutine gen_test_projected(error, mol, eeqbc)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Use the EEQBC model with the CSR neighbour list
   logical, intent(in) :: eeqbc

   class(mchrg_model_type), allocatable :: model
   class(mchrg_solver_type), allocatable :: solver
   type(csr_list) :: list
   real(wp), allocatable :: eref(:), energy(:), qref(:), qvec(:)

   allocate(eref(mol%nat), energy(mol%nat), source=0.0_wp)
   allocate(qref(mol%nat), qvec(mol%nat))

   if (eeqbc) then
      call new_eeqbc2025_model(mol, model, error)
      if (allocated(error)) return
      call new_csr_list(list, mol, error, cutoff=29.0_wp)
      if (allocated(error)) return
      call get_solver(error, solver, .true., .false.)
      if (allocated(error)) return
      call solve_model(error, mol, model, solver, eref, qref, list=list)
      if (allocated(error)) return
      call get_solver(error, solver, .true., .true.)
      if (allocated(error)) return
      call solve_model(error, mol, model, solver, energy, qvec, list=list)
   else
      call new_eeq2019_model(mol, model, error)
      if (allocated(error)) return
      call get_solver(error, solver, .false., .false.)
      if (allocated(error)) return
      call solve_model(error, mol, model, solver, eref, qref)
      if (allocated(error)) return
      call get_solver(error, solver, .true., .true.)
      if (allocated(error)) return
      call solve_model(error, mol, model, solver, energy, qvec)
   end if
   if (allocated(error)) return

   call check_array(error, "Partial charges", qvec, qref, thr)
   call check_array(error, "Energies", energy, eref, thr)

end subroutine gen_test_projected


subroutine test_pme_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(coulomb_input) :: input

   call get_structure(mol, "MB16-43", "01")
   input%method = coulomb_method%pme
   input%accuracy = tight
   call gen_test_operator(error, mol, input, thr, .true.)

end subroutine test_pme_mb01


subroutine test_pme_urea(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(coulomb_input) :: input

   call get_structure(mol, "X23", "urea")
   input%method = coulomb_method%pme
   input%accuracy = tight
   call gen_test_operator(error, mol, input, thr, .true.)

end subroutine test_pme_urea


subroutine test_pme_co2_charged(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(coulomb_input) :: input

   call get_structure(mol, "X23", "CO2")
   mol%charge = 1.0_wp
   input%method = coulomb_method%pme
   input%accuracy = tight
   call gen_test_operator(error, mol, input, thr, .true.)

end subroutine test_pme_co2_charged


subroutine test_pme_amylose(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(coulomb_input) :: input

   call get_structure(mol, "AMYLOSE", "Amylose8")
   input%method = coulomb_method%pme
   input%accuracy = tight
   call gen_test_operator(error, mol, input, thr, .true.)

end subroutine test_pme_amylose


subroutine test_pme_urea222(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(coulomb_input) :: input

   call get_structure(mol, "X23", "urea")
   call make_supercell(mol, [2, 2, 2])
   input%method = coulomb_method%pme
   input%accuracy = tight
   call gen_test_operator(error, mol, input, thr, .true.)

end subroutine test_pme_urea222


subroutine test_pme_ice(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(coulomb_input) :: input

   ! Hexagonal lattice with non-orthogonal grid axes
   call get_structure(mol, "ICE10", "ih")
   input%method = coulomb_method%pme
   input%accuracy = tight
   call gen_test_operator(error, mol, input, thr, .true.)

end subroutine test_pme_ice


subroutine test_pme_default_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(coulomb_input) :: input

   call get_structure(mol, "MB16-43", "05")
   input%method = coulomb_method%pme
   call gen_test_operator(error, mol, input, thr2, .true.)

end subroutine test_pme_default_mb05


subroutine test_projected_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "02")
   mol%charge = 1.0_wp
   call gen_test_projected(error, mol, .false.)

end subroutine test_projected_mb02


subroutine test_projected_urea(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "urea")
   call gen_test_projected(error, mol, .false.)

end subroutine test_projected_urea


subroutine test_projected_eeqbc(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "03")
   mol%charge = -1.0_wp
   call gen_test_projected(error, mol, .true.)

end subroutine test_projected_eeqbc


subroutine test_pme_external_gradient(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(coulomb_input) :: input
   class(mchrg_model_type), allocatable :: model
   class(mchrg_solver_type), allocatable :: solver
   type(mchrg_cache) :: cache
   real(wp), allocatable :: trans(:, :), qvec(:), dfdq(:)
   real(wp), allocatable :: gref(:, :), gradient(:, :), sref(:, :), sigma(:, :)
   integer :: iat

   call get_structure(mol, "MB16-43", "04")
   allocate(qvec(mol%nat), dfdq(mol%nat))
   allocate(gref(3, mol%nat), gradient(3, mol%nat), sref(3, 3), sigma(3, 3), &
      & source=0.0_wp)
   do iat = 1, mol%nat
      dfdq(iat) = 0.1_wp * real(mod(7 * iat, 5) - 2, wp)
   end do

   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call get_solver(error, solver, .false., .false.)
   if (allocated(error)) return
   call get_lattice_points(mol%periodic, mol%lattice, model%ncoord%cutoff, trans)
   call model%update(mol, cache, trans, .true.)
   call model%solve(mol, solver, cache, error, qvec=qvec)
   if (allocated(error)) return
   call model%get_external_gradient(mol, solver, cache, error, dfdq, gref, sref)
   if (allocated(error)) return

   input%method = coulomb_method%pme
   input%accuracy = tight
   call new_eeq2019_model(mol, model, error, coulomb=input)
   if (allocated(error)) return
   call get_solver(error, solver, .true., .true.)
   if (allocated(error)) return
   call model%update(mol, cache, trans, .true.)
   call model%solve(mol, solver, cache, error, qvec=qvec)
   if (allocated(error)) return
   call model%get_external_gradient(mol, solver, cache, error, dfdq, gradient, sigma)
   if (allocated(error)) return

   call check_array(error, "External gradient", reshape(gradient, [size(gradient)]), &
      & reshape(gref, [size(gref)]), thr)
   call check_array(error, "External strain derivative", reshape(sigma, [9]), &
      & reshape(sref, [9]), thr)

end subroutine test_pme_external_gradient


subroutine test_operator_direct(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(coulomb_input) :: input
   class(mchrg_model_type), allocatable :: model
   class(mchrg_solver_type), allocatable :: solver
   real(wp), allocatable :: energy(:), qvec(:)

   call get_structure(mol, "MB16-43", "01")
   allocate(energy(mol%nat), qvec(mol%nat), source=0.0_wp)
   input%method = coulomb_method%pme
   call new_eeq2019_model(mol, model, error, coulomb=input)
   if (allocated(error)) return
   call get_solver(error, solver, .false., .false.)
   if (allocated(error)) return
   call solve_model(error, mol, model, solver, energy, qvec)

end subroutine test_operator_direct


subroutine test_invalid_accuracy(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(coulomb_input) :: input
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "01")
   input%method = coulomb_method%pme
   input%accuracy = 0.0_wp
   call new_eeq2019_model(mol, model, error, coulomb=input)

end subroutine test_invalid_accuracy

end module test_coulomb
