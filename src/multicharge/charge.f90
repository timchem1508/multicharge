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
   use multicharge_model_cache, only: mchrg_cache
   use multicharge_solver_type, only : mchrg_solver_type, mchrg_solver_input
   use multicharge_solver_direct, only : direct_solver, new_direct_solver, direct_input
   use multicharge_solver_cg, only : cg_solver, new_cg_solver, cg_input
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

   type(mchrg_cache), allocatable :: cache
   logical :: grad
   real(wp), allocatable :: trans(:, :)

   class(direct_solver), allocatable :: solver
   class(direct_input), allocatable :: solver_input

   allocate(solver_input)
   allocate(solver)
   call new_direct_solver(solver, solver_input)

   grad = present(dqdr) .and. present(dqdL)

   allocate(cache)
   call get_lattice_points(mol%periodic, mol%lattice, mchrg_model%ncoord%cutoff, trans)
   call mchrg_model%update(mol, cache, trans, grad)
   call mchrg_model%solve(mol, solver, cache, error, &
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
