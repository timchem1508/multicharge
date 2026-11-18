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

!> @file multicharge/model/type.F90
!> Provides a general base class for charge models

#ifndef IK
#define IK i4
#endif

!> Abstract base type and shared operations for charge models
module multicharge_model_type
   use iso_fortran_env, only : output_unit
   use mctc_env, only : timer_type, format_time, error_type, fatal_error, wp, i8, &
      & ik => IK
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use mctc_io_math, only : matinv_3x3
   use mctc_io_structure_atomsort, only : get_spatial_order
   use mctc_cutoff, only : get_lattice_points
   use mctc_ncoord, only : ncoord_type
   use mctc_csrlist, only : csr_list, find_lower_start, spsymv_csr, spmm_csr
   use mctc_wignerseitz, only : wignerseitz_cell
   use multicharge_blas, only : gemv, symv, gemm
   use multicharge_linal, only : max_group, add_transposed_product, &
      & gather_group_list, multiply_group_list
   use multicharge_model_cache, only : mchrg_cache, clear_solve_derivs
   use multicharge_solver_type, only : mchrg_solver_type
   use multicharge_solver_cg, only : cg_solver

   implicit none
   private

   public :: mchrg_model_type, mchrg_hessian, get_dir_trans, get_rec_trans, &
      & hess_index, get_dqdr_row, new_partial_derivs, new_hessian_blocks

   ! Used by the Hessian submodule, gfortran does not link private procedures
   public :: get_inverse_list, get_constraint_factors, get_group_dqdr, &
      & get_max_neighbors, get_adjoint_vector, print_stage_time, &
      & ensure_partial_derivs, get_dqdL_list


   !> Molecular nuclear Hessian in compressed storage on a sorted list
   type :: mchrg_hessian

      !> Sorted list of the blocks, the half list of the update or the complete list
      type(csr_list), allocatable :: list

      !> Hessian blocks on the upper triangle of the list, the entry of atom j >= i
      !> in the row of atom i holds the block (j, i)
      real(wp), allocatable :: hesslist(:, :, :)

      !> Position derivatives of the Lagrange multiplier
      real(wp), allocatable :: dlambda(:, :)

      !> Second factor of the charge constraint term of an external function
      real(wp), allocatable :: wvec(:, :)

      !> Weight of the charge constraint term
      real(wp) :: weight = 0.0_wp

   end type mchrg_hessian


   !> Abstract multicharge model type
   type, abstract :: mchrg_model_type

      !> Electronegativity
      real(wp), allocatable :: chi(:)

      !> Charge width
      real(wp), allocatable :: rad(:)

      !> Chemical hardness
      real(wp), allocatable :: eta(:)

      !> CN scaling factor for electronegativity
      real(wp), allocatable :: kcnchi(:)

      !> Local charge scaling factor for electronegativity
      real(wp), allocatable :: kqchi(:)

      !> Local charge scaling factor for chemical hardness
      real(wp), allocatable :: kqeta(:)

      !> CN scaling factor for charge width
      real(wp), allocatable :: kcnrad(:)

      !> Scaling factor for the external electric field
      real(wp) :: efield_scale = 1.0_wp

      !> Coordination number
      class(ncoord_type), allocatable :: ncoord

      !> Electronegativity weighted CN for local charge
      class(ncoord_type), allocatable :: ncoord_en
   contains

      !> Solve linear equations for the charge model
      procedure :: solve

      !> Calculate the gradient of the electrostatic energy
      procedure :: get_electrostatic_gradient

      !> Get external gradient
      procedure :: get_external_gradient

      !> Calculate the derivatives of the charges w.r.t. positions and strain
      procedure :: get_charge_derivs

      !> Calculate the analytical Hessian of an external function (molecular)
      procedure :: get_external_hessian_dense

      !> Calculate the analytical Hessian of an external function on a list
      procedure :: get_external_hessian_list

      !> Calculate the analytical Hessian of an external function of the charges
      generic :: get_external_hessian => get_external_hessian_dense, &
         & get_external_hessian_list

      !> Calculate the analytical Hessian of the electrostatic energy (molecular)
      procedure :: get_electrostatic_hessian_dense

      !> Calculate the analytical Hessian of the electrostatic energy on a list
      procedure :: get_electrostatic_hessian_list

      !> Calculate the analytical Hessian of the electrostatic energy
      generic :: get_electrostatic_hessian => get_electrostatic_hessian_dense, &
         & get_electrostatic_hessian_list

      !> Calculate the implicit Hessian (molecular)
      procedure :: get_implicit_hessian_dense

      !> Calculate the implicit Hessian on a list (molecular)
      procedure :: get_implicit_hessian_list

      !> Calculate the implicit Hessian from the response of the charges
      generic :: get_implicit_hessian => get_implicit_hessian_dense, &
         & get_implicit_hessian_list

      !> Calculate semi-numerical Hessian in packed lower-triangle storage
      procedure, private :: get_numhess_packed

      !> Calculate semi-numerical Hessian in dense symmetric storage
      procedure, private :: get_numhess_dense

      !> Calculate semi-numerical Hessian and pressure tensor
      generic :: get_numhess => get_numhess_packed, get_numhess_dense

      !> Update cache
      procedure(update), deferred :: update

      !> Calculate capacitance matrix
      procedure(get_capacitance_matrix), deferred :: get_capacitance_matrix

      !> Calculate the derivatives of the capacitance matrix
      procedure(get_capacitance_derivs), deferred :: get_capacitance_derivs

      !> Calculate right-hand side (electronegativity)
      procedure(get_xvec), deferred :: get_xvec

      !> Add the electronegativity vector derivatives, beta * dX/dR
      procedure(get_xvec_derivs), deferred :: get_xvec_derivs

      !> Calculate Coulomb matrix
      procedure(get_coulomb_matrix), deferred :: get_coulomb_matrix

      !> Add the Coulomb matrix derivatives contracted with a vector
      procedure(get_coulomb_derivs), deferred :: get_coulomb_derivs

      !> Calculate capacitance-corrected electronegativity derivatives
      procedure(get_grad), deferred :: get_grad

      !> Calculate the explicit Hessian at fixed p and q (molecular)
      procedure(get_explicit_hessian_dense), deferred :: get_explicit_hessian_dense

      !> Calculate the explicit Hessian at fixed p and q on a list (molecular)
      procedure(get_explicit_hessian_list), deferred :: get_explicit_hessian_list

      !> Calculate the explicit Hessian at fixed p and q
      generic :: get_explicit_hessian => get_explicit_hessian_dense, &
         & get_explicit_hessian_list

   end type mchrg_model_type

   abstract interface
      !> Update the cache, qcp adds the missing CN derivatives and keeps the rest
      subroutine update(self, mol, cache, trans, list, qcp)
         import :: mchrg_model_type, structure_type, mchrg_cache, csr_list, wp

         !> Multicharge model type
         class(mchrg_model_type), intent(in) :: self

         !> Structure type
         type(structure_type), intent(in) :: mol

         !> Multicharge cache containing CN, local charges, and a Wigner-Seitz cell
         type(mchrg_cache), intent(inout) :: cache

         !> Lattice vectors
         real(wp), intent(in) :: trans(:, :)

         !> Multicharge neighborlist type
         type(csr_list), intent(in), optional :: list

         !> Flag for the CN derivatives of the charge response
         logical, intent(in), optional :: qcp
      end subroutine update

      !> Construct a capacitance matrix using cached CN and charge data
      subroutine get_capacitance_matrix(self, mol, ndim, cache, list)
         import :: mchrg_model_type, structure_type, mchrg_cache, csr_list, wp

         !> Multicharge model type
         class(mchrg_model_type), intent(in) :: self

         !> Structure type
         type(structure_type), intent(in) :: mol

         !> System size
         integer, intent(in) :: ndim

         !> Multicharge cache holding the capacitance matrix and optional derivatives
         type(mchrg_cache), intent(inout) :: cache

         !> Multicharge neighborlist type
         type(csr_list), intent(in), optional :: list
      end subroutine get_capacitance_matrix

      !> Missing derivatives of the capacitance matrix
      subroutine get_capacitance_derivs(self, mol, cache, list)
         import :: mchrg_model_type, structure_type, mchrg_cache, csr_list

         !> Multicharge model type
         class(mchrg_model_type), intent(in) :: self

         !> Structure type
         type(structure_type), intent(in) :: mol

         !> Multicharge cache holding the capacitance matrix and its derivatives
         type(mchrg_cache), intent(inout) :: cache

         !> Multicharge neighborlist type of the update
         type(csr_list), intent(in), optional :: list
      end subroutine get_capacitance_derivs

      !> Coulomb interaction matrix (A-matrix) construction
      subroutine get_coulomb_matrix(self, mol, ndim, cache, list)
         import :: mchrg_model_type, structure_type, mchrg_cache, csr_list, wp

         !> Multicharge model type
         class(mchrg_model_type), intent(in) :: self

         !> Structure type
         type(structure_type), intent(in) :: mol

         !> System size
         integer, intent(in) :: ndim

         !> Multicharge cache holding the Coulomb matrix
         type(mchrg_cache), intent(inout) :: cache

         !> Multicharge neighborlist type
         type(csr_list), intent(in), optional :: list
      end subroutine get_coulomb_matrix

      !> Coulomb matrix derivatives contracted with a vector, alpha * dA/dR*p
      subroutine get_coulomb_derivs(self, mol, ndim, cache, alpha, p, &
         & list, dadr, dadrlist, dadL)
         import :: mchrg_model_type, structure_type, mchrg_cache, csr_list, wp

         !> Multicharge model type
         class(mchrg_model_type), intent(in) :: self

         !> Structure type
         type(structure_type), intent(in) :: mol

         !> System size
         integer, intent(in) :: ndim

         !> Multicharge cache holding the Coulomb matrix
         type(mchrg_cache), intent(inout) :: cache

         !> dA/dR*p multiplier
         real(wp), intent(in) :: alpha

         !> Vector contracted with the Coulomb matrix derivatives
         real(wp), intent(in) :: p(:)

         !> Derivatives w.r.t. positions, dadr(:, x, b) of component b w.r.t. atom x
         real(wp), intent(inout), optional :: dadr(:, :, :)

         !> Compressed geometrical derivatives of the Coulomb matrix multiplied by p
         real(wp), intent(inout), optional :: dadrlist(:, :)

         !> Lattice derivatives of the Coulomb matrix multiplied p
         real(wp), intent(inout) :: dadL(:, :, :)

         !> Multicharge neighborlist type
         type(csr_list), intent(in), optional :: list
      end subroutine get_coulomb_derivs

      !> Electronegativity vector derivatives, beta * dX/dR
      subroutine get_xvec_derivs(self, mol, ndim, cache, beta, list)
         import :: mchrg_model_type, structure_type, mchrg_cache, csr_list, wp

         !> Multicharge model type
         class(mchrg_model_type), intent(in) :: self

         !> Structure type
         type(structure_type), intent(in) :: mol

         !> System size
         integer, intent(in) :: ndim

         !> Multicharge cache holding the partial derivatives
         type(mchrg_cache), intent(inout) :: cache

         !> Multiplier of the electronegativity vector derivatives
         real(wp), intent(in) :: beta

         !> Multicharge neighborlist type (sorted complete)
         type(csr_list), intent(in), optional :: list
      end subroutine get_xvec_derivs

      !> Electronegativity vector construction
      subroutine get_xvec(self, mol, ndim, cache, list, efield)
         import :: mchrg_model_type, mchrg_cache, structure_type, csr_list, wp

         !> Multicharge model type
         class(mchrg_model_type), intent(in) :: self

         !> Structure type
         type(structure_type), intent(in) :: mol

         !> System size
         integer, intent(in) :: ndim

         !> Multicharge cache holding the electronegativity vector and workspace
         type(mchrg_cache), intent(inout) :: cache

         !> Multicharge neighborlist type
         type(csr_list), intent(in), optional :: list

         !> External electric field
         real(wp), intent(in), optional :: efield(:)
      end subroutine get_xvec

      !> Calculate capacitance-corrected electronegativity derivatives
      subroutine get_grad(self, mol, cache, p, gradient, sigma, alpha, beta, list)
         import :: mchrg_model_type, structure_type, mchrg_cache, csr_list, wp

         !> Multicharge model type
         class(mchrg_model_type), intent(in) :: self

         !> Structure type
         type(structure_type), intent(in) :: mol

         !> Multicharge cache
         type(mchrg_cache), intent(in) :: cache

         !> Charge-like contraction vector
         real(wp), intent(in) :: p(:)

         !> Energy gradient
         real(wp), intent(inout) :: gradient(:, :)

         !> Stress tensor
         real(wp), intent(inout) :: sigma(:, :)

         !> Gradient scaling factor
         real(wp), intent(in), optional :: alpha

         !> Stress scaling factor
         real(wp), intent(in), optional :: beta

         !> Multicharge neighborlist type (sorted complete)
         type(csr_list), optional, intent(in) :: list
      end subroutine get_grad

      !> Explicit part of the molecular nuclear Hessian at fixed p and charges q
      subroutine get_explicit_hessian_dense(self, mol, cache, p, hessian, elst, &
         & error, alpha, beta)
         import :: mchrg_model_type, structure_type, mchrg_cache, error_type, wp

         !> Multicharge model type
         class(mchrg_model_type), intent(in) :: self

         !> Structure type
         type(structure_type), intent(in) :: mol

         !> Multicharge cache
         type(mchrg_cache), intent(inout) :: cache

         !> Charge-like contraction vector
         real(wp), intent(in) :: p(:)

         !> Molecular nuclear Hessian
         real(wp), intent(inout) :: hessian(:, :)

         !> Strain Hessian x -> (1 + eps) x, element (a, b, c, d) is eps_ab, eps_cd
         real(wp), intent(inout) :: elst(:, :, :, :)

         !> Error handling
         type(error_type), allocatable, intent(out) :: error

         !> Optional Coulomb matrix prefactor
         real(wp), intent(in), optional :: alpha

         !> Optional electronegativity prefactor
         real(wp), intent(in), optional :: beta
      end subroutine get_explicit_hessian_dense

      !> Explicit part of the molecular nuclear Hessian on a sorted complete list
      subroutine get_explicit_hessian_list(self, mol, cache, p, hessian, elst, &
         & error, list, alpha, beta)
         import :: mchrg_model_type, structure_type, mchrg_cache, csr_list, &
            & mchrg_hessian, error_type, wp

         !> Multicharge model type
         class(mchrg_model_type), intent(in) :: self

         !> Structure type
         type(structure_type), intent(in) :: mol

         !> Multicharge cache
         type(mchrg_cache), intent(inout) :: cache

         !> Charge-like contraction vector
         real(wp), intent(in) :: p(:)

         !> Molecular nuclear Hessian in compressed storage
         type(mchrg_hessian), intent(inout) :: hessian

         !> Strain Hessian x -> (1 + eps) x, element (a, b, c, d) is eps_ab, eps_cd
         real(wp), intent(inout) :: elst(:, :, :, :)

         !> Error handling
         type(error_type), allocatable, intent(out) :: error

         !> Multicharge neighborlist type (sorted complete)
         type(csr_list), intent(in) :: list

         !> Optional Coulomb matrix prefactor
         real(wp), intent(in), optional :: alpha

         !> Optional electronegativity prefactor
         real(wp), intent(in), optional :: beta
      end subroutine get_explicit_hessian_list

   end interface

   interface
      !> Compute the nuclear Hessian of the electrostatic energy
      module subroutine get_electrostatic_hessian_dense(self, mol, solver, cache, &
         & error, hessian, elst, gradient, sigma, unit, verbosity)

         !> Electronegativity-equilibration model
         class(mchrg_model_type), intent(in) :: self

         !> Molecular structure data
         type(structure_type), intent(in) :: mol

         !> Solver instance
         class(mchrg_solver_type), intent(in) :: solver

         !> Cache handling
         type(mchrg_cache), intent(inout) :: cache

         !> Error handling
         type(error_type), allocatable, intent(out) :: error

         !> Molecular nuclear Hessian
         real(wp), intent(out) :: hessian(:, :)

         !> Strain Hessian x -> (1 + eps) x, element (a, b, c, d) is eps_ab, eps_cd
         real(wp), intent(inout) :: elst(:, :, :, :)

         !> Optional gradient of the energy, added to
         real(wp), intent(inout), optional :: gradient(:, :)

         !> Optional virial of the energy, added to
         real(wp), intent(inout), optional :: sigma(:, :)

         !> Output unit
         integer, intent(in), optional :: unit

         !> Verbosity level
         integer, intent(in), optional :: verbosity
      end subroutine get_electrostatic_hessian_dense

      !> Compute the molecular nuclear Hessian of the electrostatic energy on a list
      module subroutine get_electrostatic_hessian_list(self, mol, solver, cache, &
         & error, hessian, elst, list, gradient, sigma, unit, verbosity)

         !> Electronegativity-equilibration model
         class(mchrg_model_type), intent(in) :: self

         !> Molecular structure data
         type(structure_type), intent(in) :: mol

         !> Solver instance
         class(mchrg_solver_type), intent(in) :: solver

         !> Cache handling
         type(mchrg_cache), intent(inout) :: cache

         !> Error handling
         type(error_type), allocatable, intent(out) :: error

         !> Molecular nuclear Hessian in compressed storage, its blocks are formed on
         !> its list, the sorted half list of the update, or on the complete list
         type(mchrg_hessian), intent(inout) :: hessian

         !> Strain Hessian x -> (1 + eps) x, element (a, b, c, d) is eps_ab, eps_cd
         real(wp), intent(inout) :: elst(:, :, :, :)

         !> Sorted complete neighborlist, walked along the complete rows
         type(csr_list), intent(in) :: list

         !> Optional gradient of the energy, added to
         real(wp), intent(inout), optional :: gradient(:, :)

         !> Optional virial of the energy, added to
         real(wp), intent(inout), optional :: sigma(:, :)

         !> Output unit
         integer, intent(in), optional :: unit

         !> Verbosity level
         integer, intent(in), optional :: verbosity
      end subroutine get_electrostatic_hessian_list

      !> Allocate the zero blocks of a compressed Hessian on the upper triangle
      !> of its list
      module subroutine new_hessian_blocks(hessian, list)

         !> Molecular nuclear Hessian in compressed storage
         type(mchrg_hessian), intent(inout) :: hessian

         !> Sorted list of the blocks, absent to keep the list of the Hessian
         type(csr_list), intent(in), optional :: list
      end subroutine new_hessian_blocks

      !> Compute the implicit Hessian of the electrostatic energy
      module subroutine get_implicit_hessian_dense(self, mol, cache, hessian, &
         & elst, dqdr, dqdL)

         !> Electronegativity-equilibration model
         class(mchrg_model_type), intent(in) :: self

         !> Molecular structure data
         type(structure_type), intent(in) :: mol

         !> Cache handling
         type(mchrg_cache), intent(in) :: cache

         !> Molecular nuclear Hessian
         real(wp), intent(inout) :: hessian(:, :)

         !> Strain Hessian x -> (1 + eps) x, element (a, b, c, d) is eps_ab, eps_cd
         real(wp), intent(inout) :: elst(:, :, :, :)

         !> Derivatives of the charges w.r.t. positions, dqdr(:, x, b) = dq_b/dR_x
         real(wp), intent(in) :: dqdr(:, :, :)

         !> Derivatives of the charges w.r.t. strain
         real(wp), intent(in) :: dqdL(:, :, :)
      end subroutine get_implicit_hessian_dense

      !> Compute the implicit Hessian of the electrostatic energy on a list
      module subroutine get_implicit_hessian_list(self, mol, cache, hessian, &
         & elst, list, times)

         !> Electronegativity-equilibration model
         class(mchrg_model_type), intent(in) :: self

         !> Molecular structure data
         type(structure_type), intent(in) :: mol

         !> Cache handling
         type(mchrg_cache), intent(in) :: cache

         !> Molecular nuclear Hessian in compressed storage
         type(mchrg_hessian), intent(inout) :: hessian

         !> Strain Hessian x -> (1 + eps) x, element (a, b, c, d) is eps_ab, eps_cd
         real(wp), intent(inout) :: elst(:, :, :, :)

         !> Multicharge neighborlist type (sorted complete)
         type(csr_list), intent(in) :: list

         !> Times of the group derivatives and their contraction
         real(wp), intent(out), optional :: times(2)
      end subroutine get_implicit_hessian_list

      !> Compute the nuclear Hessian of an external function of the charges
      module subroutine get_external_hessian_dense(self, mol, solver, cache, &
         & error, dfdq, d2fdq2, dfdrdq, hessian, elst, gradient, sigma, unit, &
         & verbosity)

         !> Electronegativity-equilibration model
         class(mchrg_model_type), intent(in) :: self

         !> Molecular structure data
         type(structure_type), intent(in) :: mol

         !> Solver instance
         class(mchrg_solver_type), intent(in) :: solver

         !> Cache handling
         type(mchrg_cache), intent(inout) :: cache

         !> Error handling
         type(error_type), allocatable, intent(out) :: error

         !> Derivative of the external function w.r.t. atomic charges
         real(wp), intent(in) :: dfdq(:)

         !> Second derivative of the external function w.r.t. atomic charges
         real(wp), intent(in) :: d2fdq2(:, :)

         !> Mixed derivative of the external function w.r.t. positions and charges
         real(wp), intent(in) :: dfdrdq(:, :, :)

         !> Nuclear Hessian of the external function without d2f/dR2
         real(wp), intent(out) :: hessian(:, :)

         !> Strain Hessian x -> (1 + eps) x, element (a, b, c, d) is eps_ab, eps_cd
         real(wp), intent(inout) :: elst(:, :, :, :)

         !> Optional gradient of the external function, added to
         real(wp), intent(inout), optional :: gradient(:, :)

         !> Optional virial of the external function, added to
         real(wp), intent(inout), optional :: sigma(:, :)

         !> Output unit
         integer, intent(in), optional :: unit

         !> Verbosity level
         integer, intent(in), optional :: verbosity
      end subroutine get_external_hessian_dense

      !> Compute the nuclear Hessian of an external function on a list
      module subroutine get_external_hessian_list(self, mol, solver, cache, &
         & error, dfdq, d2fdq2, d2fdq2_list, dfdrdq, dfdrdq_list, hessian, &
         & elst, list, gradient, sigma, unit, verbosity)

         !> Electronegativity-equilibration model
         class(mchrg_model_type), intent(in) :: self

         !> Molecular structure data
         type(structure_type), intent(in) :: mol

         !> Solver instance
         class(mchrg_solver_type), intent(in) :: solver

         !> Cache handling
         type(mchrg_cache), intent(inout) :: cache

         !> Error handling
         type(error_type), allocatable, intent(out) :: error

         !> Derivative of the external function w.r.t. atomic charges
         real(wp), intent(in) :: dfdq(:)

         !> Second derivative of the external function w.r.t. atomic charges
         real(wp), intent(in), contiguous :: d2fdq2(:)

         !> Sorted complete list of the second derivative w.r.t. atomic charges
         type(csr_list), intent(in) :: d2fdq2_list

         !> Mixed derivative of the external function w.r.t. positions and charges
         real(wp), intent(in) :: dfdrdq(:, :)

         !> Sorted complete list of the mixed derivative w.r.t. positions and charges
         type(csr_list), intent(in) :: dfdrdq_list

         !> Molecular nuclear Hessian in compressed storage
         type(mchrg_hessian), intent(inout) :: hessian

         !> Strain Hessian x -> (1 + eps) x, element (a, b, c, d) is eps_ab, eps_cd
         real(wp), intent(inout) :: elst(:, :, :, :)

         !> Neighborlist type
         type(csr_list), intent(in) :: list

         !> Optional gradient of the external function, added to
         real(wp), intent(inout), optional :: gradient(:, :)

         !> Optional virial of the external function, added to
         real(wp), intent(inout), optional :: sigma(:, :)

         !> Output unit
         integer, intent(in), optional :: unit

         !> Verbosity level
         integer, intent(in), optional :: verbosity
      end subroutine get_external_hessian_list

      !> Semi-numerical Hessian and pressure tensor
      module subroutine get_numhess_packed(self, mol, solver, cache, error, &
         & qvec, energy, grad, sigma, hess, press, list, unit, verbosity)

         !> Electronegativity-equilibration model
         class(mchrg_model_type), intent(in) :: self

         !> Molecular structure data
         type(structure_type), intent(in) :: mol

         !> The solver instance
         class(mchrg_solver_type), intent(in) :: solver

         !> Cache handling for the unperturbed system
         type(mchrg_cache), intent(inout) :: cache

         !> Error handling
         type(error_type), allocatable, intent(out) :: error

         !> Atomic partial charges of the unperturbed system
         real(wp), intent(out), contiguous :: qvec(:)

         !> Electrostatic energy of the unperturbed system
         real(wp), intent(inout), contiguous :: energy(:)

         !> Energy gradient of the unperturbed system
         real(wp), intent(inout), contiguous :: grad(:, :)

         !> Virial of the unperturbed system
         real(wp), intent(inout), contiguous :: sigma(:, :)

         !> Hessian matrix d2E/dR2 in packed lower-triangle storage
         real(wp), intent(out) :: hess(:)

         !> Virial derivatives w.r.t. strain (3, 3, 3, 3)
         real(wp), intent(out) :: press(:, :, :, :)

         !> neighborlist optional type
         type(csr_list), intent(in), optional :: list

         !> Output unit
         integer, intent(in), optional :: unit

         !> Verbosity level
         integer, intent(in), optional :: verbosity
      end subroutine get_numhess_packed

      !> Semi-numerical Hessian in dense symmetric storage
      module subroutine get_numhess_dense(self, mol, solver, cache, error, &
         & qvec, energy, grad, sigma, hess, press, list, unit, verbosity)

         !> Electronegativity-equilibration model
         class(mchrg_model_type), intent(in) :: self

         !> Molecular structure data
         type(structure_type), intent(in) :: mol

         !> The solver instance
         class(mchrg_solver_type), intent(in) :: solver

         !> Cache handling for the unperturbed system
         type(mchrg_cache), intent(inout) :: cache

         !> Error handling
         type(error_type), allocatable, intent(out) :: error

         !> Atomic partial charges of the unperturbed system
         real(wp), intent(out), contiguous :: qvec(:)

         !> Electrostatic energy of the unperturbed system
         real(wp), intent(inout), contiguous :: energy(:)

         !> Energy gradient of the unperturbed system
         real(wp), intent(inout), contiguous :: grad(:, :)

         !> Virial of the unperturbed system
         real(wp), intent(inout), contiguous :: sigma(:, :)

         !> Hessian matrix d2E/dR2
         real(wp), intent(out) :: hess(:, :)

         !> Virial derivatives w.r.t. positions (3, 3, 3, nat) or strain (3, 3, 3, 3)
         real(wp), intent(out) :: press(:, :, :, :)

         !> neighborlist optional type
         type(csr_list), intent(in), optional :: list

         !> Output unit
         integer, intent(in), optional :: unit

         !> Verbosity level
         integer, intent(in), optional :: verbosity
      end subroutine get_numhess_dense

      !> Position of a symmetric-matrix element in packed lower-triangle storage
      elemental module function hess_index(idx, jdx) result(ij)

         !> Row index
         integer, intent(in) :: idx

         !> Column index
         integer, intent(in) :: jdx

         !> Position in the packed lower triangle
         integer :: ij
      end function hess_index
   end interface

   !> Twice pi
   real(wp), parameter :: twopi = 2.0_wp * pi

   !> Smallest positive working-precision number
   real(wp), parameter :: eps = tiny(1.0_wp)


contains


!> Generate direct lattice translation vectors for a periodic structure
subroutine get_dir_trans(mol, trans, cutoff)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Translation vectors
   !> Shape: (3, ntrans)
   real(wp), allocatable, intent(out) :: trans(:, :)

   !> Optional lattice-vector cutoff
   real(wp), intent(in), optional :: cutoff

   integer, parameter :: rep(3) = [2, 2, 2]

   if (present(cutoff)) then
      call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans)
   else
      call get_lattice_points(mol%lattice, rep, .true., trans)
   end if

end subroutine get_dir_trans


!> Generate reciprocal lattice translation vectors for a periodic structure
subroutine get_rec_trans(mol, trans, cutoff)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Reciprocal translation vectors
   !> Shape: (3, ntrans)
   real(wp), allocatable, intent(out) :: trans(:, :)

   !> Optional reciprocal-vector cutoff
   real(wp), intent(in), optional :: cutoff

   integer, parameter :: rep(3) = [2, 2, 2]
   real(wp) :: rec_lat(3, 3)

   rec_lat = twopi * transpose(matinv_3x3(mol%lattice))
   if (present(cutoff)) then
      call get_lattice_points(mol%periodic, rec_lat, cutoff, trans)
   else
      call get_lattice_points(rec_lat, rep, .false., trans)
   end if

end subroutine get_rec_trans


!> Solve the charge model for the charges and the electrostatic energy
subroutine solve(self, mol, solver, cache, error, energy, qvec, list, efield, &
   & verbosity, unit)

   !> Electronegativity-equilibration model
   class(mchrg_model_type), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> The solver instance
   class(mchrg_solver_type), intent(in) :: solver

   !> Cache handling
   type(mchrg_cache), intent(inout) :: cache

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Optional electrostatic energy result
   real(wp), intent(inout), contiguous, optional :: energy(:)

   !> Optional atomic partial charges result
   real(wp), intent(out), contiguous, optional :: qvec(:)

   !> Optional sorted neighborlist of the update, see sort_neighborlist
   type(csr_list), intent(in), optional :: list

   !> Optional external electric field
   real(wp), intent(in), contiguous, optional :: efield(:)

   !> Optional print verbosity
   integer, intent(in), optional :: verbosity

   !> Output unit
   integer, intent(in), optional :: unit

   integer :: iat, ndim, print_unit, verbosity_solve
   real(wp) :: uvecsum, vvecsum, lambda
   real(wp), allocatable :: unitvec(:), vvec(:), avec(:)
   logical :: add_lagr
   type(timer_type) :: timer

   ! The augmented system of the direct solver is only set up as dense matrix
   if (.not. solver%need_pos_def .and. present(list)) then
      call fatal_error(error, "The direct solver does not support a neighborlist")
      return
   end if

   verbosity_solve = 0
   if (present(verbosity)) verbosity_solve = verbosity
   print_unit = output_unit
   if (present(unit)) print_unit = unit

   ! The CG solver requires a positive-definite system
   add_lagr = .not. solver%need_pos_def
   ndim = mol%nat
   if (add_lagr) ndim = mol%nat + 1

   ! New matrices and charges invalidate their inverse and derivatives
   call clear_solve_derivs(cache)

   call timer%push("total")
   call timer%push("setup")
   call self%get_capacitance_matrix(mol, ndim, cache, list)
   call self%get_coulomb_matrix(mol, ndim, cache, list)
   call self%get_xvec(mol, ndim, cache, list, efield)
   if (.not. allocated(cache%vrhs)) allocate(cache%vrhs(mol%nat + 1))
   call timer%pop
   call print_solve_header(print_unit, verbosity_solve, timer%get("setup"))

   if (add_lagr) then
      cache%vrhs = cache%xvec
      call solver%solve(amat=cache%amat, xvec=cache%xvec, vrhs=cache%vrhs, &
         & new_unit=print_unit, error=error)
      if (allocated(error)) return
   else
      if (.not. allocated(cache%uvec)) allocate(cache%uvec(mol%nat))
      allocate(unitvec(mol%nat), source=1.0_wp)
      allocate(vvec(mol%nat))
      ! Initial guess from the diagonal
      if (present(list)) then
         do iat = 1, mol%nat
            cache%uvec(iat) = 1.0_wp / (cache%alist(list%inl(iat)) + eps)
            vvec(iat) = - cache%xvec(iat) / (cache%alist(list%inl(iat)) + eps)
         end do
      else
         do iat = 1, mol%nat
            cache%uvec(iat) = 1.0_wp / (cache%amat(iat, iat) + eps)
            vvec(iat) = - cache%xvec(iat) / (cache%amat(iat, iat) + eps)
         end do
      end if

      ! Constrained response: A*uvec = 1
      call print_constrained_system_message(print_unit, verbosity_solve, 'u')
      call solver%solve(amat=cache%amat, alist=cache%alist, xvec=unitvec, &
         & vrhs=cache%uvec, list=list, new_unit=print_unit, error=error)
      if (allocated(error)) return

      ! Unconstrained response: A*vvec = -xvec
      call print_constrained_system_message(print_unit, verbosity_solve, 'v')
      call solver%solve(amat=cache%amat, alist=cache%alist, xvec=-cache%xvec, &
         & vrhs=vvec, list=list, new_unit=print_unit, error=error)
      if (allocated(error)) return

      ! Lagrange multiplier and projection of uvec on vvec
      uvecsum = sum(cache%uvec)
      vvecsum = sum(vvec)
      lambda = - (mol%charge + vvecsum) / (uvecsum + eps)
      cache%vrhs(:mol%nat) = -vvec - lambda * cache%uvec
      cache%vrhs(mol%nat + 1) = lambda
   end if

   if (present(qvec)) qvec(:) = cache%vrhs(:mol%nat)

   if (present(energy)) then
      call timer%push("energy")
      if (present(list)) then
         allocate(avec(mol%nat))
         call spsymv_csr(mol%nat, cache%alist, list%inl, list%nlat, cache%vrhs, &
            & avec, uplo='U')
         cache%xvec(:mol%nat) = 0.5_wp * avec - cache%xvec(:mol%nat)
      else
         call symv(cache%amat, cache%vrhs, cache%xvec(:mol%nat), &
            & alpha=0.5_wp, beta=-1.0_wp, uplo='l')
      end if
      ! Exclude the constraint term
      if (ndim > mol%nat) then
         cache%xvec(:mol%nat) = cache%xvec(:mol%nat) - 0.5_wp * cache%vrhs(mol%nat + 1)
      end if
      energy(:) = energy(:) + cache%vrhs(:mol%nat) * cache%xvec(:mol%nat)
      call timer%pop
      call print_energy_time(print_unit, verbosity_solve, timer%get("energy"))
   end if

   call timer%pop
   call print_total_time(print_unit, verbosity_solve, timer%get("total"))

end subroutine solve


!> Gradient and virial of the electrostatic energy at the charges of a solve
subroutine get_electrostatic_gradient(self, mol, cache, error, gradient, sigma, &
   & list, unit, verbosity)

   !> Electronegativity-equilibration model
   class(mchrg_model_type), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Cache handling
   type(mchrg_cache), intent(inout) :: cache

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Gradient of the electrostatic energy, added to
   real(wp), intent(inout), contiguous :: gradient(:, :)

   !> Virial of the electrostatic energy, added to
   real(wp), intent(inout), contiguous :: sigma(:, :)

   !> Optional sorted neighborlist of the update
   type(csr_list), intent(in), optional :: list

   !> Output unit
   integer, intent(in), optional :: unit

   !> Verbosity level
   integer, intent(in), optional :: verbosity

   integer :: print_unit, verbosity_grad
   type(timer_type) :: timer

   if (.not. allocated(cache%vrhs)) then
      call fatal_error(error, "Charges are not in the cache, solve first")
      return
   end if

   verbosity_grad = 0
   if (present(verbosity)) verbosity_grad = verbosity
   print_unit = output_unit
   if (present(unit)) print_unit = unit

   call timer%push("gradient")
   call self%get_capacitance_derivs(mol, cache, list)
   call self%get_grad(mol, cache, cache%vrhs(:mol%nat), gradient, sigma, &
      & alpha=0.5_wp, beta=-1.0_wp, list=list)
   call timer%pop
   call print_gradient_time(print_unit, verbosity_grad, timer%get("gradient"))

end subroutine get_electrostatic_gradient


!> Derivatives of the charges w.r.t. atomic positions and lattice vectors
subroutine get_charge_derivs(self, mol, solver, cache, error, dqdr, dqdL, &
   & dqdrlist, dlambdadr, list, unit, verbosity)

   !> Electronegativity-equilibration model
   class(mchrg_model_type), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> The solver instance
   class(mchrg_solver_type), intent(in) :: solver

   !> Cache handling
   type(mchrg_cache), intent(inout) :: cache

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Optional charge derivatives w.r.t. positions, dense matrices only
   real(wp), intent(out), contiguous, optional :: dqdr(:, :, :)

   !> Optional charge derivatives w.r.t. lattice vectors
   real(wp), intent(out), contiguous, optional :: dqdL(:, :, :)

   !> Optional charge derivatives w.r.t. positions on the list
   real(wp), allocatable, intent(out), optional :: dqdrlist(:, :)

   !> Optional position derivatives of the Lagrange multiplier on the list
   real(wp), allocatable, intent(out), optional :: dlambdadr(:, :)

   !> Optional sorted complete neighborlist of the update
   type(csr_list), intent(in), optional :: list

   !> Output unit
   integer, intent(in), optional :: unit

   !> Verbosity level
   integer, intent(in), optional :: verbosity

   integer :: jat, print_unit, verbosity_q
   real(wp) :: scale, uvecsum
   real(wp), allocatable :: xvec(:), vtmp(:), dlam(:, :)
   type(timer_type) :: timer

   call check_charge_derivs(mol, solver, cache, error, list)
   if (allocated(error)) return
   if (present(list) .and. present(dqdr)) then
      call fatal_error(error, "Dense charge derivatives require dense matrices, "// &
         & "use dqdrlist on a neighborlist")
      return
   end if

   verbosity_q = 0
   if (present(verbosity)) verbosity_q = verbosity
   print_unit = output_unit
   if (present(unit)) print_unit = unit

   ! Right-hand sides of the response equations, dX/dR - dA/dR*q
   call timer%push("total")
   call timer%push("partial")
   call ensure_partial_derivs(self, mol, cache, list)
   call timer%pop
   call print_gradient_header(print_unit, verbosity_q, timer%get("partial"))

   if (present(list)) then
      call timer%push("inverse")
      select type (solver)
      class is (cg_solver)
         call get_inverse_list(mol, solver, cache, list, error, unit)
      end select
      call timer%pop
      if (allocated(error)) return

      call timer%push("product")
      call get_dqdr_list(mol, cache, list)
      call get_constraint_factors(mol, cache, list, dlam)
      call move_alloc(dlam, cache%dlambdadr)
      if (present(dqdL)) call get_dqdL_list(mol, cache, list, dqdL)
      call timer%pop
      if (present(dqdrlist)) dqdrlist = cache%dqdrlist
      if (present(dlambdadr)) dlambdadr = cache%dlambdadr

      if (verbosity_q > 0) then
         call print_stage_time(print_unit, "Inverse on the neighborlist", &
            & timer%get("inverse"))
         call print_stage_time(print_unit, "Charge derivatives S*dB/dR", &
            & timer%get("product"))
         write(print_unit, '(a)') ''
      end if
   else
      ! Projected inverse P J^-1 in ainv(:nat, :nat) for both solvers
      if (.not. allocated(cache%ainv)) then
         select type (solver)
         class is (cg_solver)
            allocate(cache%ainv(mol%nat, mol%nat))
            call solver%invert(amat=cache%amat, ainv=cache%ainv, new_unit=unit, &
               & error=error)
            if (.not. allocated(error)) then
               ! (P J^-1)^T = J^-1 - u*u^T/sum(u)
               uvecsum = sum(cache%uvec)
               !$omp parallel do default(none) schedule(runtime) &
               !$omp shared(mol, cache, uvecsum) private(jat, scale)
               do jat = 1, mol%nat
                  scale = cache%uvec(jat) / (uvecsum + eps)
                  cache%ainv(:, jat) = cache%ainv(:, jat) - scale * cache%uvec
               end do
            end if
         class default
            ! Inverse of the augmented matrix
            allocate(cache%ainv(size(cache%amat, 1), size(cache%amat, 2)))
            allocate(xvec(size(cache%amat, 1)), source=0.0_wp)
            allocate(vtmp(size(cache%amat, 1)))
            call solver%solve(amat=cache%amat, xvec=xvec, vrhs=vtmp, &
               & ainv=cache%ainv, cpq=.true., new_unit=unit, error=error)
         end select
         if (allocated(error)) then
            deallocate(cache%ainv)
            return
         end if
      end if

      if (present(dqdr)) call gemm(cache%dabdr, cache%ainv(:mol%nat, :mol%nat), dqdr)
      if (present(dqdL)) call gemm(cache%dabdL, cache%ainv(:mol%nat, :mol%nat), dqdL)
   end if

   call timer%pop
   call print_gradient_time(print_unit, verbosity_q, timer%get("total"))

end subroutine get_charge_derivs


!> Check the cache and the setup of the charge derivatives
subroutine check_charge_derivs(mol, solver, cache, error, list)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> The solver instance
   class(mchrg_solver_type), intent(in) :: solver

   !> Cache handling
   type(mchrg_cache), intent(in) :: cache

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Optional sorted complete neighborlist of the update
   type(csr_list), intent(in), optional :: list

   logical :: same_list

   if (.not. allocated(cache%vrhs)) then
      call fatal_error(error, "Charges are not in the cache, solve first")
      return
   end if

   if (.not. present(list)) then
      if (.not. allocated(cache%dcndr)) then
         call fatal_error(error, "Charge derivatives require the CN derivatives, "// &
            & "update the cache with qcp")
      end if
      return
   end if

   if (.not. list%complete) then
      call fatal_error(error, "The charge derivatives require a sorted "// &
         & "complete list, see sort_neighborlist")
      return
   end if
   select type (solver)
   class is (cg_solver)
   class default
      call fatal_error(error, "Charge derivatives on a neighborlist require "// &
         & "the CG solver")
      return
   end select
   if (.not. allocated(cache%dcndrlist)) then
      call fatal_error(error, "Charge derivatives require the CN derivatives "// &
         & "on the neighborlist, update the cache with qcp")
      return
   end if
   same_list = size(cache%hinl) == size(list%inl)
   if (same_list) same_list = all(cache%hinl == list%inl)
   if (.not. same_list) then
      call fatal_error(error, "The cache was updated with a different neighborlist")
      return
   end if
   if (size(list%inl) /= mol%nat + 1) then
      call fatal_error(error, "The neighborlist does not match the structure")
   end if

end subroutine check_charge_derivs


!> Partial derivatives dX/dR - dA/dR*q, formed only if not in the cache
subroutine ensure_partial_derivs(self, mol, cache, list)

   !> Electronegativity-equilibration model
   class(mchrg_model_type), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Cache with the charges and the CN derivatives
   type(mchrg_cache), intent(inout) :: cache

   !> Optional sorted complete neighborlist
   type(csr_list), intent(in), optional :: list

   if (present(list)) then
      if (allocated(cache%dabdrlist)) return
   else if (allocated(cache%dabdr)) then
      return
   end if

   call self%get_capacitance_derivs(mol, cache, list)
   call new_partial_derivs(mol, mol%nat, cache, list)
   call self%get_xvec_derivs(mol, mol%nat, cache, beta=1.0_wp, list=list)
   call self%get_coulomb_derivs(mol, mol%nat, cache, alpha=-1.0_wp, &
      & p=cache%vrhs(:mol%nat), list=list, dadr=cache%dabdr, &
      & dadrlist=cache%dabdrlist, dadL=cache%dabdL)

end subroutine ensure_partial_derivs


!> Inverse of the coefficient matrix on a sorted neighborlist, kept if in the cache
subroutine get_inverse_list(mol, solver, cache, list, error, unit, half)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Conjugate-gradient solver with the tolerance of the inverse
   class(cg_solver), intent(in) :: solver

   !> Cache with the coefficient matrix on the list, receives the inverse
   type(mchrg_cache), intent(inout) :: cache

   !> Sorted complete neighborlist
   type(csr_list), intent(in) :: list

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Output unit
   integer, intent(in), optional :: unit

   !> Sorted list of the update holding the matrix and the inverse on its upper
   !> triangle, absent for the complete list
   type(csr_list), intent(in), optional :: half

   integer(i8) :: nnz

   if (present(half)) then
      nnz = size(half%nlat, kind=i8)
   else
      nnz = size(list%nlat, kind=i8)
   end if
   if (allocated(cache%ainvlist)) then
      if (size(cache%ainvlist, kind=i8) == nnz) return
      deallocate(cache%ainvlist)
   end if
   allocate(cache%ainvlist(nnz))
   call solver%invert_list(cache%alist, list, cache%ainvlist, xyz=mol%xyz, &
      & new_unit=unit, error=error, half=half)

end subroutine get_inverse_list


!> Position derivatives of the charges on a sorted complete neighborlist
subroutine get_dqdr_list(mol, cache, list)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Cache with the inverse and the partial derivatives on the list
   type(mchrg_cache), intent(inout) :: cache

   !> Sorted complete neighborlist
   type(csr_list), intent(in) :: list

   integer :: igrp, ngrp, nbmax, x0, ng, nb, n, i
   integer(i8) :: kat
   integer, allocatable :: order(:), bloc(:), blist(:)
   integer(i8), allocatable :: ufin(:), qpos(:, :)
   real(wp), allocatable :: pan(:, :), qcol(:, :)

   call get_spatial_order(mol%xyz, order)
   ngrp = (mol%nat + max_group - 1) / max_group
   nbmax = min(mol%nat, max_group * get_max_neighbors(mol%nat, list))

   allocate(ufin(mol%nat))
   do i = 1, mol%nat
      ufin(i) = find_lower_start(list, i) - 1_i8
   end do

   if (allocated(cache%dqdrlist)) deallocate(cache%dqdrlist)
   allocate(cache%dqdrlist(3, size(list%nlat, kind=i8)), source=0.0_wp)

   ! Only the group of atom x writes the derivatives w.r.t. R_x in the row of x
   !$omp parallel default(none) &
   !$omp shared(mol, cache, list, ufin, order, ngrp, nbmax) &
   !$omp private(igrp, x0, ng, nb, n, i, kat, bloc, blist, qpos, pan, qcol)
   allocate(bloc(mol%nat), source=0)
   allocate(blist(nbmax), qpos(max_group, nbmax), pan(3 * max_group, nbmax), &
      & qcol(3 * max_group, nbmax))
   !$omp do schedule(runtime)
   do igrp = 1, ngrp
      x0 = max_group * (igrp - 1)
      ng = min(max_group, mol%nat - x0)
      call get_group_dqdr(cache, list, list, ufin, order(x0 + 1:x0 + ng), nb, &
         & bloc, blist, qpos, pan, qcol)
      do n = 1, nb
         do i = 1, ng
            kat = qpos(i, n)
            if (kat > 0) cache%dqdrlist(:, kat) = qcol(3 * i - 2:3 * i, n)
         end do
         bloc(blist(n)) = 0
      end do
   end do
   !$omp end do
   deallocate(bloc, blist, qpos, pan, qcol)
   !$omp end parallel

end subroutine get_dqdr_list


!> Charge constraint factors of the position derivatives of the charges
subroutine get_constraint_factors(mol, cache, list, dlambdadr)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Cache with the partial derivatives and the constraint response
   type(mchrg_cache), intent(in) :: cache

   !> Sorted complete neighborlist
   type(csr_list), intent(in) :: list

   !> Position derivatives of the Lagrange multiplier
   real(wp), allocatable, intent(out) :: dlambdadr(:, :)

   allocate(dlambdadr(3, mol%nat), source=0.0_wp)
   call add_transposed_product(list, cache%dabdrlist, cache%uvec, dlambdadr)
   dlambdadr(:, :) = dlambdadr / (sum(cache%uvec) + eps)

end subroutine get_constraint_factors


!> Derivatives of the charges w.r.t. the lattice vectors on a list
subroutine get_dqdL_list(mol, cache, list, dqdL)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Cache with the inverse and the strain derivatives
   type(mchrg_cache), intent(in) :: cache

   !> Sorted complete neighborlist
   type(csr_list), intent(in) :: list

   !> Derivative of the atomic partial charges w.r.t. lattice vectors
   real(wp), intent(out) :: dqdL(:, :, :)

   integer :: iat, jat
   real(wp) :: uvecsum, scale
   real(wp) :: wlat(3, 3)
   real(wp), allocatable :: dabdL(:, :), sdabdL(:, :)

   uvecsum = sum(cache%uvec)
   wlat(:, :) = 0.0_wp
   do iat = 1, mol%nat
      wlat(:, :) = wlat + cache%dabdL(:, :, iat) * cache%uvec(iat)
   end do

   ! The nine strain components of the atoms are the columns of the product
   dabdL = transpose(reshape(cache%dabdL(:, :, :mol%nat), [9, mol%nat]))
   allocate(sdabdL(mol%nat, 9))
   call spmm_csr('N', mol%nat, 9, mol%nat, 1.0_wp, 'SU', cache%ainvlist, list%nlat, &
      & list%inl(:mol%nat), list%inl(2:), dabdL, mol%nat, 0.0_wp, sdabdL, mol%nat)
   do jat = 1, mol%nat
      scale = cache%uvec(jat) / (uvecsum + eps)
      dqdL(:, :, jat) = reshape(sdabdL(jat, :), [3, 3]) - scale * wlat
   end do

end subroutine get_dqdL_list


!> Charge derivatives of a group of atoms for the union of their neighbors
subroutine get_group_dqdr(cache, list, half, hfin, group, nb, bloc, blist, qpos, &
   & pan, qcol)

   !> Cache with the inverse and the partial derivatives on the list
   type(mchrg_cache), intent(in) :: cache

   !> Sorted complete neighborlist of the partial derivatives
   type(csr_list), intent(in) :: list

   !> Sorted list of the update with the inverse on its upper triangle
   type(csr_list), intent(in) :: half

   !> End of the upper triangle in each row of the list of the inverse
   integer(i8), intent(in), contiguous :: hfin(:)

   !> Atoms of the group
   integer, intent(in), contiguous :: group(:)

   !> Number of atoms in the union of the neighbors of the group
   integer, intent(out) :: nb

   !> Index of each atom in the union, zero outside
   integer, intent(inout), contiguous :: bloc(:)

   !> Atoms of the union
   integer, intent(inout), contiguous :: blist(:)

   !> Position of the entry of b in the row of x, zero if not a neighbor
   integer(i8), intent(inout), contiguous :: qpos(:, :)

   !> Partial derivatives dB_c/dR_x of the union in the layout of qcol
   real(wp), intent(inout), contiguous :: pan(:, :)

   !> Charge derivatives dq_b/dR_x of the union
   real(wp), intent(inout), contiguous :: qcol(:, :)

   integer :: n, i

   call gather_group_list(list, cache%dabdrlist, group, nb, bloc, blist, qpos, pan)
   call multiply_group_list(half, hfin, cache%ainvlist, nb, bloc, blist, pan, qcol)
   do n = 1, nb
      do i = 1, size(group)
         if (qpos(i, n) == 0) qcol(3 * i - 2:3 * i, n) = 0.0_wp
      end do
   end do

end subroutine get_group_dqdr


!> Largest number of entries of a row of a complete neighborlist
function get_max_neighbors(nat, list) result(kmax)

   !> Number of atoms
   integer, intent(in) :: nat

   !> Complete neighborlist
   type(csr_list), intent(in) :: list

   !> Largest number of neighbors
   integer :: kmax

   kmax = 0
   if (nat > 0) kmax = int(maxval(list%inl(2:nat + 1) - list%inl(:nat)))

end function get_max_neighbors


!> Allocate the partial derivatives in the cache and set them to zero
subroutine new_partial_derivs(mol, ndim, cache, list, dafdr)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> System size
   integer, intent(in) :: ndim

   !> Multicharge cache
   type(mchrg_cache), intent(inout) :: cache

   !> Multicharge neighborlist type (sorted complete)
   type(csr_list), intent(in), optional :: list

   !> Allocate dafdr for an external vector instead of the partial derivatives
   logical, intent(in), optional :: dafdr

   logical :: ext
   integer(i8) :: nnz

   ext = .false.
   if (present(dafdr)) ext = dafdr

   if (present(list)) then
      nnz = size(list%nlat, kind=i8)
      if (ext) then
         if (allocated(cache%dafdrlist)) deallocate(cache%dafdrlist)
         allocate(cache%dafdrlist(3, nnz), source=0.0_wp)
      else
         if (allocated(cache%dabdrlist)) deallocate(cache%dabdrlist)
         allocate(cache%dabdrlist(3, nnz), source=0.0_wp)
      end if
   else if (ext) then
      if (allocated(cache%dafdr)) deallocate(cache%dafdr)
      allocate(cache%dafdr(3, mol%nat, mol%nat), source=0.0_wp)
   else
      if (allocated(cache%dabdr)) deallocate(cache%dabdr)
      allocate(cache%dabdr(3, mol%nat, ndim), source=0.0_wp)
   end if
   if (ext) then
      if (allocated(cache%dafdL)) deallocate(cache%dafdL)
      allocate(cache%dafdL(3, 3, ndim), source=0.0_wp)
   else
      if (allocated(cache%dabdL)) deallocate(cache%dabdL)
      allocate(cache%dabdL(3, 3, ndim), source=0.0_wp)
   end if

end subroutine new_partial_derivs


!> Position derivatives of the charge of one atom
subroutine get_dqdr_row(cache, list, jat, dqdr)

   !> Cache with the compressed position derivatives
   type(mchrg_cache), intent(in) :: cache

   !> Sorted complete neighborlist of the charge derivatives
   type(csr_list), intent(in) :: list

   !> Atom of the charge
   integer, intent(in) :: jat

   !> Derivatives dq_j/dR_i for all atoms i, dimension (3, nat)
   real(wp), intent(out) :: dqdr(:, :)

   integer :: iat
   integer(i8) :: kat, pos

   do iat = 1, size(dqdr, 2)
      dqdr(:, iat) = -cache%dlambdadr(:, iat) * cache%uvec(jat)
   end do

   ! dq_j/dR_i is held at the entry of jat in the row of each neighbor iat
   do kat = list%inl(jat), list%inl(jat + 1) - 1
      iat = list%nlat(kat)
      do pos = list%inl(iat), list%inl(iat + 1) - 1
         if (list%nlat(pos) == jat) then
            dqdr(:, iat) = dqdr(:, iat) + cache%dqdrlist(:, pos)
            exit
         end if
      end do
   end do

end subroutine get_dqdr_row


!> External gradient calculation using the adjoint state method
subroutine get_external_gradient(self, mol, solver, cache, error, &
   & dfdq, dfdr, dfdL, list, unit, verbosity)

   !> Electronegativity-equilibration model
   class(mchrg_model_type), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Solver instance
   class(mchrg_solver_type), intent(in) :: solver

   !> Cache handling
   type(mchrg_cache), intent(inout) :: cache

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Derivative of the objective w.r.t. atomic partial charges
   real(wp), intent(in) :: dfdq(:)

   !> External gradient w.r.t. positions
   real(wp), intent(inout) :: dfdr(:, :)

   !> External gradient w.r.t. lattice vectors
   real(wp), intent(inout) :: dfdL(:, :)

   !> Neighborlist optional type
   type(csr_list), intent(in), optional :: list

   !> Output unit
   integer, intent(in), optional :: unit

   !> Verbosity level
   integer, intent(in), optional :: verbosity

   real(wp), allocatable :: padj(:)
   integer :: print_unit, verbosity_solve
   type(timer_type) :: timer

   verbosity_solve = 0
   if (present(verbosity)) verbosity_solve = verbosity
   print_unit = output_unit
   if (present(unit)) print_unit = unit

   if (size(dfdq) > mol%nat) then
      call fatal_error(error, "External partial derivative is wrong size")
      return
   end if

   call timer%push("setup_external")
   call timer%pop

   call print_gradient_header(print_unit, verbosity_solve, &
      & timer%get("setup_external"))
   call timer%push("external_gradient")

   if (.not. allocated(cache%vrhs)) then
      call fatal_error(error, "Charges are not in the cache, solve first")
      return
   end if
   if (.not. (allocated(cache%amat) .or. allocated(cache%alist))) then
      call fatal_error(error, "Coefficient matrix is not in the cache, solve again")
      return
   end if
   if (.not. allocated(cache%uvec) .and. solver%need_pos_def) then
      call fatal_error(error, "Constrained response J*uvec = 1 is not allocated")
      return
   end if

   call get_adjoint_vector(mol, solver, cache, error, dfdq, padj, list, &
      & print_unit, verbosity_solve)
   if (allocated(error)) return

   ! External gradient by adjoint contraction
   call self%get_capacitance_derivs(mol, cache, list)
   call self%get_grad(mol, cache, padj(:mol%nat), dfdr, dfdL, alpha=-1.0_wp, &
      & beta=1.0_wp, list=list)

   ! pop dfdr timer
   call timer%pop
   call print_gradient_time(print_unit, verbosity_solve, &
      & timer%get("external_gradient"))

end subroutine get_external_gradient


!> Adjoint vector p = S df/dq of an external function of the charges.
subroutine get_adjoint_vector(mol, solver, cache, error, dfdq, padj, list, unit, &
   & verbosity)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Solver instance
   class(mchrg_solver_type), intent(in) :: solver

   !> Cache with the coefficient matrix and the constraint response
   type(mchrg_cache), intent(in) :: cache

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Derivative of the external function w.r.t. atomic charges
   real(wp), intent(in) :: dfdq(:)

   !> Adjoint vector, of the dimension of the linear system
   real(wp), allocatable, intent(out) :: padj(:)

   !> Multicharge neighborlist type (sorted complete)
   type(csr_list), intent(in), optional :: list

   !> Output unit
   integer, intent(in), optional :: unit

   !> Verbosity level
   integer, intent(in), optional :: verbosity

   integer :: iat, print_unit, verbosity_solve
   real(wp) :: scale
   real(wp), allocatable :: dfdq_loc(:)

   verbosity_solve = 0
   if (present(verbosity)) verbosity_solve = verbosity
   print_unit = output_unit
   if (present(unit)) print_unit = unit

   call print_adjoint_message(print_unit, verbosity_solve)
   if (solver%need_pos_def) then
      ! Unconstrained response J*y = df/dq from the diagonal guess
      allocate(padj(mol%nat))
      if (present(list)) then
         do iat = 1, mol%nat
            padj(iat) = dfdq(iat) / (cache%alist(list%inl(iat)) + eps)
         end do
      else
         do iat = 1, mol%nat
            padj(iat) = dfdq(iat) / (cache%amat(iat, iat) + eps)
         end do
      end if
      call solver%solve(amat=cache%amat, alist=cache%alist, xvec=dfdq, &
         & vrhs=padj, list=list, new_unit=print_unit, error=error)
      if (allocated(error)) return

      ! Projection of the constraint response u on y
      scale = sum(padj) / (sum(cache%uvec) + eps)
      padj(:) = padj - scale * cache%uvec
   else
      allocate(dfdq_loc(mol%nat + 1), source=0.0_wp)
      dfdq_loc(:mol%nat) = dfdq
      allocate(padj(mol%nat + 1))
      call solver%solve(amat=cache%amat, xvec=dfdq_loc, vrhs=padj, &
         & new_unit=print_unit, error=error)
   end if

end subroutine get_adjoint_vector


!> Print header for charge equilibration solver
subroutine print_solve_header(unit, verbosity, timer)
   !> Output unit
   integer, intent(in) :: unit

   !> Verbosity level
   integer, intent(in) :: verbosity

   !> Elapsed setup time
   real(wp), intent(in) :: timer

   if (verbosity > 0) then
      write(unit, '(54("-"))')
      write(unit, '(13x, a)') "Charge equilibration solver"
      write(unit, '(54("-"))')
      write(unit, '(a)') ''
      if (verbosity > 1) then
         write(unit, '(a, 1x, a)') "Setup time : ", format_time(timer)
         write(unit, '(a)') ''
      end if
   end if
end subroutine print_solve_header


!> Print header for gradient calculations
subroutine print_gradient_header(unit, verbosity, timer)
   !> Output unit
   integer, intent(in) :: unit

   !> Verbosity level
   integer, intent(in) :: verbosity

   !> Elapsed gradient setup time
   real(wp), intent(in) :: timer

   if (verbosity > 0) then
      write(unit, '(54("-"))')
      write(unit, '(17x, a)') "Gradient Calculations"
      write(unit, '(54("-"))')
      write(unit, '(a)') ''
      if (verbosity > 1) then
         write(unit, '(a, 1x, a)') "Gradient setup time : ", format_time(timer)
         write(unit, '(a)') ''
      end if
   end if
end subroutine print_gradient_header


!> Print message for constrained system solves
subroutine print_constrained_system_message(unit, verbosity, vector)
   !> Output unit
   integer, intent(in) :: unit

   !> Verbosity level
   integer, intent(in) :: verbosity

   !> Constrained-system identifier
   character, intent(in) :: vector

   if (verbosity > 0) then
      if (vector == 'u') then
         write(unit, '(a)') 'Solving constrained system: J*u = 1'
      else if (vector == 'v') then
         write(unit, '(a)') 'Solving unconstrained system: J*v = chi'
      else if (vector == 'y') then
         write(unit, '(a)') 'Solving derivative unconstrained system: J*y = df/dq'
      end if
      write(unit, '(a)') ''
   end if
end subroutine print_constrained_system_message


!> Print message for adjoint system solve
subroutine print_adjoint_message(unit, verbosity)
   !> Output unit
   integer, intent(in) :: unit

   !> Verbosity level
   integer, intent(in) :: verbosity

   if (verbosity > 0) then
      write(unit, '(a)') 'Solving adjoint system: J*y = dfdq'
      write(unit, '(a)') ''
   end if
end subroutine print_adjoint_message


!> Print gradient calculation time
subroutine print_gradient_time(unit, verbosity, timer)
   !> Output unit
   integer, intent(in) :: unit

   !> Verbosity level
   integer, intent(in) :: verbosity

   !> Elapsed gradient calculation time
   real(wp), intent(in) :: timer

   if (verbosity > 1) then
      write(unit, '(a, 1x, a)') "Gradient calculation time : ", format_time(timer)
      write(unit, '(a)') ''
   end if
end subroutine print_gradient_time


!> Print energy calculation time
subroutine print_energy_time(unit, verbosity, timer)
   !> Output unit
   integer, intent(in) :: unit

   !> Verbosity level
   integer, intent(in) :: verbosity

   !> Elapsed energy calculation time
   real(wp), intent(in) :: timer

   if (verbosity > 1) then
      write(unit, '(a, 1x, a)') "Energy calculation time : ", format_time(timer)
      write(unit, '(a)') ''
   end if
end subroutine print_energy_time


!> Print total solve time
subroutine print_total_time(unit, verbosity, timer)
   !> Output unit
   integer, intent(in) :: unit

   !> Verbosity level
   integer, intent(in) :: verbosity

   !> Elapsed total solve time
   real(wp), intent(in) :: timer

   if (verbosity > 1) then
      write(unit, '(a, 1x, a)') "Total solve time : ", format_time(timer)
      write(unit, '(a)') ''
   end if
end subroutine print_total_time


!> Print the elapsed time of a stage of a calculation
subroutine print_stage_time(unit, label, time)
   !> Output unit
   integer, intent(in) :: unit

   !> Name of the stage
   character(len=*), intent(in) :: label

   !> Elapsed time in seconds
   real(wp), intent(in) :: time

   write(unit, '(a, t34, ":", a)') label, format_time(time)
end subroutine print_stage_time


end module multicharge_model_type
