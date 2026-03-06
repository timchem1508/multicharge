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
   use mctc_env, only : error_type, wp
   use mctc_io, only : structure_type
   use mctc_cutoff, only : get_lattice_points
   use multicharge_model_type, only : mchrg_model_type
   use multicharge_solver_type, only : mchrg_solver_type, mchrg_solver_input
   use multicharge_solver_direct, only : direct_solver, new_direct_solver, direct_input, direct_cache
   use multicharge_solver_cg, only : cg_solver, new_cg_solver, cg_input, cg_cache
   use multicharge_solver_cache, only: mchrg_solver_cache
   use multicharge_param, only : new_eeq2019_model, new_eeqbc2025_model

   implicit none
   private

   public :: get_charges, get_eeq_charges, get_eeqbc_charges

contains


!> Classical electronegativity equilibration charges
subroutine get_charges(mchrg_model, mol, error, qvec, dqdr, dqdL)
   
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

   logical :: grad
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: qloc(:), dqlocdr(:, :, :), dqlocdL(:, :, :)
   real(wp), allocatable :: trans(:, :)

   class(mchrg_solver_type), allocatable :: solver
   class(mchrg_solver_input), allocatable :: solver_input

   allocate(direct_input :: solver_input)
   select type(solver_input)
   type is (direct_input)
      block
          class(direct_solver), allocatable :: tmp
          allocate(tmp)
          call new_direct_solver(tmp, solver_input)
          call move_alloc(tmp, solver)
      end block
   end select

   grad = present(dqdr) .and. present(dqdL)

   allocate(cn(mol%nat), qloc(mol%nat))
   if (grad) then
      allocate(dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat))
      allocate (dqlocdr(3, mol%nat, mol%nat), dqlocdL(3, 3, mol%nat))
   end if

   call get_lattice_points(mol%periodic, mol%lattice, mchrg_model%ncoord%cutoff, trans)
   call mchrg_model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)
   call mchrg_model%local_charge(mol, trans, qloc, dqlocdr, dqlocdL)
   call mchrg_model%solve(mol, solver, error, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL, &
      & qvec=qvec, dqdr=dqdr, dqdL=dqdL, unit=output_unit)

end subroutine get_charges


!> Obtain charges from electronegativity equilibration model
subroutine get_eeq_charges(mol, error, qvec, dqdr, dqdL)

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

   class(mchrg_model_type), allocatable :: eeq_model

   call new_eeq2019_model(mol, eeq_model, error)

   call get_charges(eeq_model, mol, error, qvec, dqdr, dqdL)

end subroutine get_eeq_charges


!> Obtain charges from bond capacity electronegativity equilibration model
subroutine get_eeqbc_charges(mol, error, qvec, dqdr, dqdL)

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

   class(mchrg_model_type), allocatable :: eeqbc_model

   call new_eeqbc2025_model(mol, eeqbc_model, error)

   call get_charges(eeqbc_model, mol, error, qvec, dqdr, dqdL)

end subroutine get_eeqbc_charges


end module multicharge_charge
