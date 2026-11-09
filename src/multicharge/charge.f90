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

!> @file multicharge/charge.f90
!> Contains functions to calculate the partial charges with and without
!> separate charge model setup

!> Interface to the charge models
module multicharge_charge
   use iso_fortran_env, only : output_unit
   use mctc_env, only : error_type, fatal_error, wp
   use mctc_io, only : structure_type
   use mctc_cutoff, only : get_lattice_points
   use mctc_csrlist, only : csr_list
   use multicharge_model_type, only : mchrg_model_type
   use multicharge_model_cache, only : mchrg_cache
   use multicharge_solver_type, only : mchrg_solver_type, mchrg_solver_input
   use multicharge_solver_cg, only : cg_solver, new_cg_solver, cg_input
   use multicharge_solver_direct, only : direct_solver, new_direct_solver, direct_input
   use multicharge_param, only : new_eeq2019_model, new_eeqbc2025_model
   implicit none
   private

   public :: get_charges, get_eeq_charges, get_eeqbc_charges, add_charge_hessian


contains


!> Classical electronegativity equilibration charges
subroutine get_charges(mchrg_model, mol, error, qvec, dqdr, dqdL, list, efield, &
   & cache)

   !> Multicharge model
   class(mchrg_model_type), intent(in) :: mchrg_model

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Atomic partial charges
   real(wp), intent(out), contiguous :: qvec(:)

   !> Derivative of the partial charges w.r.t. the Cartesian coordinates
   real(wp), intent(out), contiguous, optional :: dqdr(:, :, :)

   !> Derivative of the partial charges w.r.t. strain deformations
   real(wp), intent(out), contiguous, optional :: dqdL(:, :, :)

   !> Sorted complete neighborlist, see sort_neighborlist
   type(csr_list), intent(in), optional :: list

   !> Optional external electric field
   real(wp), intent(in), contiguous, optional :: efield(:)

   !> Optional cache retaining the solve and the CN derivatives for
   !> add_charge_hessian, recompute after any change of the structure or model
   type(mchrg_cache), intent(inout), optional :: cache

   type(mchrg_cache), allocatable :: local_cache
   type(mchrg_cache) :: empty_cache
   logical :: grad, hess
   real(wp), allocatable :: trans(:, :)

   class(cg_solver), allocatable :: solver
   class(cg_input), allocatable :: solver_input

   allocate(solver_input)
   solver_input%use_nlist = present(list)
   allocate(solver)
   call new_cg_solver(solver, solver_input)

   ! With a neighborlist the position derivatives can stay in the cache
   grad = present(dqdL) .and. present(dqdr)
   hess = present(cache)
   if (present(list)) then
      if ((grad .or. hess) .and. .not. list%complete) then
         call fatal_error(error, "Neighborlist must be sorted and complete "// &
            & "when computing charge gradients or Hessian")
         return
      end if
   end if

   call get_lattice_points(mol%periodic, mol%lattice, mchrg_model%ncoord%cutoff, trans)
   if (present(cache)) then
      ! Start from an empty cache, the Hessian checks for allocated data
      cache = empty_cache
      call mchrg_model%update(mol, cache, trans, grad, list=list, dqdr=grad, &
         & hess=.true.)
      call mchrg_model%solve(mol, solver, cache, error, &
         & qvec=qvec, dqdr=dqdr, dqdL=dqdL, list=list, &
         & efield=efield, unit=output_unit)
   else
      allocate(local_cache)
      call mchrg_model%update(mol, local_cache, trans, grad, list=list, dqdr=grad)
      call mchrg_model%solve(mol, solver, local_cache, error, &
         & qvec=qvec, dqdr=dqdr, dqdL=dqdL, list=list, &
         & efield=efield, unit=output_unit)
   end if

end subroutine get_charges


!> Add the Hessian of the charges contracted with fixed atomic weights,
!> sum_i weights(i) d2q_i/dR2, at a fixed lattice. The cache has to come from
!> get_charges with the same structure and model and without a neighborlist.
subroutine add_charge_hessian(mchrg_model, mol, error, weights, hessian, cache)

   !> Multicharge model
   class(mchrg_model_type), intent(in) :: mchrg_model

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Weights of the atomic partial charges
   real(wp), intent(in) :: weights(:)

   !> Nuclear Hessian of shape (3*nat, 3*nat), the contraction is added
   real(wp), intent(inout) :: hessian(:, :)

   !> Cache retained by get_charges, holds partial derivatives afterwards
   type(mchrg_cache), intent(inout) :: cache

   real(wp), allocatable :: whess(:, :)
   class(cg_solver), allocatable :: solver
   class(cg_input), allocatable :: solver_input

   if (any(shape(hessian) /= 3 * mol%nat)) then
      call fatal_error(error, "Invalid dimensions for charge Hessian contraction")
      return
   end if

   ! Same solver as in get_charges, the constraint response of its solve is
   ! kept in the cache
   allocate(solver_input)
   allocate(solver)
   call new_cg_solver(solver, solver_input)

   ! The weighted sum of the charges is linear in the charges, its Hessian is
   ! the external Hessian with df/dq = weights and no further derivatives
   allocate(whess(3 * mol%nat, 3 * mol%nat))
   call mchrg_model%get_external_hessian(mol, solver, cache, error, weights, &
      & hessian=whess)
   if (allocated(error)) return
   hessian(:, :) = hessian + whess

end subroutine add_charge_hessian


!> Obtain charges from electronegativity equilibration model
subroutine get_eeq_charges(mol, error, qvec, dqdr, dqdL, efield)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Atomic partial charges
   real(wp), intent(out), contiguous :: qvec(:)

   !> Derivative of the partial charges w.r.t. the Cartesian coordinates
   real(wp), intent(out), contiguous, optional :: dqdr(:, :, :)

   !> Derivative of the partial charges w.r.t. strain deformations
   real(wp), intent(out), contiguous, optional :: dqdL(:, :, :)

   !> Optional external electric field
   real(wp), intent(in), contiguous, optional :: efield(:)

   class(mchrg_model_type), allocatable :: eeq_model

   call new_eeq2019_model(mol, eeq_model, error)

   call get_charges(eeq_model, mol, error, qvec, dqdr, dqdL, efield=efield)

end subroutine get_eeq_charges


!> Obtain charges from bond capacity electronegativity equilibration model
subroutine get_eeqbc_charges(mol, error, qvec, dqdr, dqdL, list, efield)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Atomic partial charges
   real(wp), intent(out), contiguous :: qvec(:)

   !> Derivative of the partial charges w.r.t. the Cartesian coordinates
   real(wp), intent(out), contiguous, optional :: dqdr(:, :, :)

   !> Derivative of the partial charges w.r.t. strain deformations
   real(wp), intent(out), contiguous, optional :: dqdL(:, :, :)

   !> neighborlist
   type(csr_list), intent(in), optional :: list

   !> Optional external electric field
   real(wp), intent(in), contiguous, optional :: efield(:)

   class(mchrg_model_type), allocatable :: eeqbc_model

   call new_eeqbc2025_model(mol, eeqbc_model, error)

   call get_charges(eeqbc_model, mol, error, qvec, dqdr, dqdL, list=list, efield=efield)

end subroutine get_eeqbc_charges


end module multicharge_charge
