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

!> Tests for the weighted first charge response from the get_charges cache
module test_charge_gradient
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, test_failed
   use mctc_io_structure, only : structure_type, new
   use mstore, only : get_structure
   use multicharge, only : get_charges, add_charge_gradient, add_charge_hessian, &
      & mchrg_cache, mchrg_model_type, new_eeq2019_model, new_eeqbc2025_model
   use multicharge_solver_direct, only : direct_solver, new_direct_solver, &
      & direct_input
   implicit none
   private

   public :: collect_charge_gradient

   real(wp), parameter :: thr1 = 1.0e5_wp * epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

   !> Skew lattice for the compact periodic test cell
   real(wp), parameter :: skew_lattice(3, 3) = reshape([ &
      & 5.1_wp, 0.0_wp, 0.0_wp, 0.4_wp, 6.2_wp, 0.0_wp, &
      & 0.2_wp, 0.7_wp, 7.1_wp], [3, 3])

   !> Positions of the Zn-O-H test system
   real(wp), parameter :: znoh_xyz(3, 3) = reshape([ &
      & 0.0_wp, 0.0_wp, 0.0_wp, 3.0_wp, 0.1_wp, 0.2_wp, &
      & 0.5_wp, 2.0_wp, 0.3_wp], [3, 3])


contains


!> Collect all exported unit tests
subroutine collect_charge_gradient(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("eeq-molecule", test_eeq_molecule), &
      & new_unittest("eeqbc-molecule", test_eeqbc_molecule), &
      & new_unittest("eeq-charged", test_eeq_charged), &
      & new_unittest("eeqbc-charged", test_eeqbc_charged), &
      & new_unittest("eeq-periodic", test_eeq_periodic), &
      & new_unittest("eeqbc-periodic", test_eeqbc_periodic), &
      & new_unittest("eeq-numerical", test_eeq_numerical), &
      & new_unittest("eeqbc-numerical", test_eeqbc_numerical), &
      & new_unittest("accumulate", test_accumulate), &
      & new_unittest("charge-conservation", test_conservation), &
      & new_unittest("cache-reuse", test_cache_reuse), &
      & new_unittest("external-gradient", test_external_gradient) &
      & ]

end subroutine collect_charge_gradient


!> Fixed, non-uniform weights for the atomic charges
pure function get_weights(nat) result(weights)

   !> Number of atoms
   integer, intent(in) :: nat

   !> Weights of the atomic partial charges
   real(wp) :: weights(nat)

   integer :: iat

   weights = [(sin(real(2 * iat, wp)), iat = 1, nat)]

end function get_weights


!> Compare the weighted response with the contraction of dq/dR and dq/dL
subroutine test_contraction_gen(error, mol, model)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Charge model
   class(mchrg_model_type), intent(in) :: model

   integer :: iat
   type(mchrg_cache) :: cache
   real(wp) :: sigma(3, 3), sigma_ref(3, 3)
   real(wp), allocatable :: qvec(:), dqdr(:, :, :), dqdL(:, :, :), weights(:)
   real(wp), allocatable :: gradient(:, :), gradient_ref(:, :)

   allocate(qvec(mol%nat), dqdr(3, mol%nat, mol%nat), dqdL(3, 3, mol%nat), &
      & gradient(3, mol%nat), gradient_ref(3, mol%nat))
   weights = get_weights(mol%nat)

   call get_charges(model, mol, error, qvec, dqdr, dqdL, cache=cache)
   if (allocated(error)) return

   gradient_ref(:, :) = 0.0_wp
   sigma_ref(:, :) = 0.0_wp
   do iat = 1, mol%nat
      gradient_ref(:, :) = gradient_ref + weights(iat) * dqdr(:, :, iat)
      sigma_ref(:, :) = sigma_ref + weights(iat) * dqdL(:, :, iat)
   end do

   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp
   call add_charge_gradient(model, mol, error, weights, gradient, sigma, cache)
   if (allocated(error)) return

   if (any(abs(gradient - gradient_ref) > thr1)) then
      call test_failed(error, "Weighted charge gradient does not match dq/dR")
      print "(3es21.12)", gradient - gradient_ref
      return
   end if
   if (any(abs(sigma - sigma_ref) > thr1)) then
      call test_failed(error, "Weighted charge strain derivative does not match dq/dL")
      print "(3es21.12)", sigma - sigma_ref
   end if

end subroutine test_contraction_gen

!> Compare the weighted response with central differences of the weighted charges
subroutine test_numerical_gen(error, mol, model)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Charge model
   class(mchrg_model_type), intent(in) :: model

   real(wp), parameter :: step = 1.0e-5_wp

   integer :: iat, ic
   type(mchrg_cache) :: cache
   real(wp) :: sigma(3, 3), qr, ql
   real(wp), allocatable :: qvec(:), weights(:), gradient(:, :), numgrad(:, :)

   allocate(qvec(mol%nat), gradient(3, mol%nat), numgrad(3, mol%nat))
   weights = get_weights(mol%nat)

   do iat = 1, mol%nat
      do ic = 1, 3
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call get_charges(model, mol, error, qvec)
         if (allocated(error)) return
         qr = dot_product(weights, qvec)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2 * step
         call get_charges(model, mol, error, qvec)
         if (allocated(error)) return
         ql = dot_product(weights, qvec)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numgrad(ic, iat) = 0.5_wp * (qr - ql) / step
      end do
   end do

   call get_charges(model, mol, error, qvec, cache=cache)
   if (allocated(error)) return
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp
   call add_charge_gradient(model, mol, error, weights, gradient, sigma, cache)
   if (allocated(error)) return

   if (any(abs(gradient - numgrad) > thr2)) then
      call test_failed(error, &
         & "Weighted charge gradient differs from finite differences")
      print "(3es21.12)", gradient - numgrad
   end if

end subroutine test_numerical_gen

subroutine test_eeq_molecule(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "01")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_contraction_gen(error, mol, model)

end subroutine test_eeq_molecule

subroutine test_eeqbc_molecule(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "01")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_contraction_gen(error, mol, model)

end subroutine test_eeqbc_molecule

subroutine test_eeq_charged(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call new(mol, [30, 8, 1], znoh_xyz, charge=1.0_wp)
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_contraction_gen(error, mol, model)

end subroutine test_eeq_charged

subroutine test_eeqbc_charged(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call new(mol, [30, 8, 1], znoh_xyz, charge=1.0_wp)
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_contraction_gen(error, mol, model)

end subroutine test_eeqbc_charged

subroutine test_eeq_periodic(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "X23", "urea")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_contraction_gen(error, mol, model)
   if (allocated(error)) return

   call new(mol, [30, 8, 1], znoh_xyz, charge=1.0_wp, lattice=skew_lattice, &
      & periodic=[.true., .true., .true.])
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_contraction_gen(error, mol, model)

end subroutine test_eeq_periodic

subroutine test_eeqbc_periodic(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "X23", "urea")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_contraction_gen(error, mol, model)
   if (allocated(error)) return

   call new(mol, [30, 8, 1], znoh_xyz, charge=1.0_wp, lattice=skew_lattice, &
      & periodic=[.true., .true., .true.])
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_contraction_gen(error, mol, model)

end subroutine test_eeqbc_periodic

subroutine test_eeq_numerical(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "01")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_numerical_gen(error, mol, model)

end subroutine test_eeq_numerical

subroutine test_eeqbc_numerical(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "01")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_numerical_gen(error, mol, model)

end subroutine test_eeqbc_numerical

!> The response is added to existing gradient and strain derivative arrays
subroutine test_accumulate(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   type(mchrg_cache) :: cache
   real(wp) :: sigma(3, 3), sigma0(3, 3), offset(3, 3)
   real(wp), allocatable :: qvec(:), weights(:), gradient(:, :), gradient0(:, :)

   call get_structure(mol, "X23", "urea")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   allocate(qvec(mol%nat), gradient(3, mol%nat), gradient0(3, mol%nat))
   weights = get_weights(mol%nat)

   call get_charges(model, mol, error, qvec, cache=cache)
   if (allocated(error)) return
   gradient0(:, :) = 0.0_wp
   sigma0(:, :) = 0.0_wp
   call add_charge_gradient(model, mol, error, weights, gradient0, sigma0, cache)
   if (allocated(error)) return

   offset(:, :) = reshape([0.1_wp, -0.2_wp, 0.3_wp, -0.4_wp, 0.5_wp, -0.6_wp, &
      & 0.7_wp, -0.8_wp, 0.9_wp], [3, 3])
   gradient(:, :) = 0.25_wp
   sigma(:, :) = offset
   call add_charge_gradient(model, mol, error, weights, gradient, sigma, cache)
   if (allocated(error)) return

   if (any(abs(gradient - 0.25_wp - gradient0) > thr1) &
      & .or. any(abs(sigma - offset - sigma0) > thr1)) then
      call test_failed(error, "Charge gradient does not accumulate")
   end if

end subroutine test_accumulate

!> Uniform weights measure the total charge, which does not respond
subroutine test_conservation(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   type(mchrg_cache) :: cache
   real(wp) :: sigma(3, 3)
   real(wp), allocatable :: qvec(:), weights(:), gradient(:, :)

   call new(mol, [30, 8, 1], znoh_xyz, charge=1.0_wp, lattice=skew_lattice, &
      & periodic=[.true., .true., .true.])
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   allocate(qvec(mol%nat), weights(mol%nat), gradient(3, mol%nat))
   weights(:) = 1.0_wp

   call get_charges(model, mol, error, qvec, cache=cache)
   if (allocated(error)) return
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp
   call add_charge_gradient(model, mol, error, weights, gradient, sigma, cache)
   if (allocated(error)) return

   if (any(abs(gradient) > thr1) .or. any(abs(sigma) > thr1)) then
      call test_failed(error, "First response violates total charge conservation")
   end if

end subroutine test_conservation

!> A charge Hessian contraction leaves the response cache unchanged
subroutine test_cache_reuse(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   type(mchrg_cache) :: cache
   real(wp) :: sigma(3, 3), sigma0(3, 3)
   real(wp), allocatable :: qvec(:), weights(:), gradient(:, :), gradient0(:, :)
   real(wp), allocatable :: hessian(:, :)

   call new(mol, [30, 8, 1], znoh_xyz, charge=1.0_wp)
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   allocate(qvec(mol%nat), gradient(3, mol%nat), gradient0(3, mol%nat), &
      & hessian(3 * mol%nat, 3 * mol%nat))
   weights = get_weights(mol%nat)

   call get_charges(model, mol, error, qvec, cache=cache)
   if (allocated(error)) return
   gradient0(:, :) = 0.0_wp
   sigma0(:, :) = 0.0_wp
   call add_charge_gradient(model, mol, error, weights, gradient0, sigma0, cache)
   if (allocated(error)) return

   hessian(:, :) = 0.0_wp
   call add_charge_hessian(model, mol, error, weights, hessian, cache)
   if (allocated(error)) return

   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp
   call add_charge_gradient(model, mol, error, weights, gradient, sigma, cache)
   if (allocated(error)) return

   if (any(abs(gradient - gradient0) > thr1) &
      & .or. any(abs(sigma - sigma0) > thr1)) then
      call test_failed(error, "Charge Hessian changed the response cache")
   end if

end subroutine test_cache_reuse

!> The adjoint gradient with a direct solver agrees on the same cache
subroutine test_external_gradient(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   type(mchrg_cache) :: cache
   type(direct_input) :: solver_input
   type(direct_solver) :: solver
   real(wp) :: sigma(3, 3), sigma_ref(3, 3)
   real(wp), allocatable :: qvec(:), weights(:), gradient(:, :)
   real(wp), allocatable :: gradient_ref(:, :)

   call get_structure(mol, "X23", "urea")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   allocate(qvec(mol%nat), gradient(3, mol%nat), gradient_ref(3, mol%nat))
   weights = get_weights(mol%nat)
   call new_direct_solver(solver, solver_input)

   call get_charges(model, mol, error, qvec, cache=cache)
   if (allocated(error)) return
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp
   call add_charge_gradient(model, mol, error, weights, gradient, sigma, cache)
   if (allocated(error)) return

   gradient_ref(:, :) = 0.0_wp
   sigma_ref(:, :) = 0.0_wp
   call model%get_external_gradient(mol, solver, cache, error, weights, &
      & gradient_ref, sigma_ref)
   if (allocated(error)) return

   if (any(abs(gradient - gradient_ref) > thr1) &
      & .or. any(abs(sigma - sigma_ref) > thr1)) then
      call test_failed(error, &
         & "Cached charge gradient differs from the adjoint external gradient")
   end if

end subroutine test_external_gradient

end module test_charge_gradient
