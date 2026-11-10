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
   use mctc_csrlist, only : csr_list, new_csr_list, distance_screening, &
      & find_lower_start, sort_neighborlist, spsymv_csr, spmm_csr
   use mctc_wignerseitz, only : wignerseitz_cell
   use multicharge_blas, only : gemv, symv, gemm
   use multicharge_model_cache, only : mchrg_cache
   use multicharge_solver_type, only : mchrg_solver_type
   use multicharge_solver_cg, only : cg_solver

   implicit none
   private

   public :: mchrg_model_type, mchrg_hessian, get_dir_trans, get_rec_trans, &
      & hess_index, get_dqdr_row, new_partial_derivs, new_hessian_blocks


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

      !> Get external gradient
      procedure :: get_external_gradient

      !> Calculate the analytical Hessian of an external function (molecular)
      procedure :: get_external_hessian_dense

      !> Calculate the analytical Hessian of an external function on a list
      procedure :: get_external_hessian_list

      !> Calculate the analytical Hessian of an external function of the charges
      generic :: get_external_hessian => get_external_hessian_dense, &
         & get_external_hessian_list

      !> Calculate the analytical Hessian of the electrostatic energy (molecular)
      procedure :: get_hessian_dense

      !> Calculate the analytical Hessian of the electrostatic energy on a list
      procedure :: get_hessian_list

      !> Calculate the analytical Hessian of the electrostatic energy
      generic :: get_hessian => get_hessian_dense, get_hessian_list

      !> Calculate the implicit Hessian (molecular)
      procedure :: get_implicit_hessian_dense

      !> Calculate the implicit Hessian on a list (molecular)
      procedure :: get_implicit_hessian_list

      !> Calculate the implicit Hessian from the response of the charges
      generic :: get_implicit_hessian => get_implicit_hessian_dense, &
         & get_implicit_hessian_list

      !> Calculate local charges from electronegativity weighted CN
      procedure :: local_charge

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
      !> Update model-dependent quantities and cache
      subroutine update(self, mol, cache, trans, grad, list, dqdr, hess)
         import :: mchrg_model_type, structure_type, mchrg_cache, csr_list, wp

         !> Multicharge model type
         class(mchrg_model_type), intent(in) :: self

         !> Structure type
         type(structure_type), intent(in) :: mol

         !> Multicharge neighborlist type
         type(csr_list), intent(in), optional :: list

         !> Multicharge cache containing CN, local charges, and a Wigner-Seitz cell
         type(mchrg_cache), intent(inout) :: cache

         !> Lattice vectors
         real(wp), intent(in) :: trans(:, :)

         !> Flag to compute derivatives
         logical, intent(in) :: grad

         !> Flag for the charge derivatives w.r.t. positions and strain
         logical, intent(in), optional :: dqdr

         !> Flag for the Hessian
         logical, intent(in), optional :: hess
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
      subroutine get_coulomb_derivs(self, mol, ndim, cache, alpha, p, list, dafdr)
         import :: mchrg_model_type, structure_type, mchrg_cache, csr_list, wp

         !> Multicharge model type
         class(mchrg_model_type), intent(in) :: self

         !> Structure type
         type(structure_type), intent(in) :: mol

         !> System size
         integer, intent(in) :: ndim

         !> Multicharge cache holding the partial derivatives
         type(mchrg_cache), intent(inout) :: cache

         !> Multiplier of the Coulomb matrix derivatives
         real(wp), intent(in) :: alpha

         !> Vector contracted with the Coulomb matrix derivatives
         real(wp), intent(in) :: p(:)

         !> Multicharge neighborlist type (sorted complete)
         type(csr_list), intent(in), optional :: list

         !> Add to dafdr for an external vector instead
         logical, intent(in), optional :: dafdr
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
      subroutine get_explicit_hessian_dense(self, mol, cache, p, hessian, error, &
         & alpha, beta)
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

         !> Error handling
         type(error_type), allocatable, intent(out) :: error

         !> Optional Coulomb matrix prefactor
         real(wp), intent(in), optional :: alpha

         !> Optional electronegativity prefactor
         real(wp), intent(in), optional :: beta
      end subroutine get_explicit_hessian_dense

      !> Explicit part of the molecular nuclear Hessian on a sorted complete list
      subroutine get_explicit_hessian_list(self, mol, cache, p, hessian, error, &
         & list, alpha, beta)
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

   !> Twice pi
   real(wp), parameter :: twopi = 2.0_wp * pi

   !> Smallest positive working-precision number
   real(wp), parameter :: eps = tiny(1.0_wp)

   !> Number of atoms in a group of the charge derivatives on a neighborlist
   integer, parameter :: max_group = 8


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


!> Top-level solve routine with optional persistent cache
subroutine solve(self, mol, solver, cache, error, &
   & energy, gradient, sigma, qvec, dqdr, dqdL, list, efield, verbosity, unit)

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

   !> Optional atomic partial charges result
   real(wp), intent(out), contiguous, optional :: qvec(:)

   !> Optional electrostatic energy result
   real(wp), intent(inout), contiguous, optional :: energy(:)

   !> Optional gradient for electrostatic energy
   real(wp), intent(inout), contiguous, optional :: gradient(:, :)

   !> Optional stress tensor for electrostatic energy
   real(wp), intent(inout), contiguous, optional :: sigma(:, :)

   !> Optional derivative of the atomic partial charges w.r.t. atomic positions
   real(wp), intent(out), contiguous, optional :: dqdr(:, :, :)

   !> Optional derivative of the atomic partial charges w.r.t. lattice vectors
   real(wp), intent(out), contiguous, optional :: dqdL(:, :, :)

   !> Optional sorted complete neighborlist, see sort_neighborlist
   type(csr_list), intent(in), optional :: list

   !> Optional external electric field
   real(wp), intent(in), contiguous, optional :: efield(:)

   !> Optional print verbossity number input flag
   integer, intent(in), optional :: verbosity

   !> Output unit
   integer, intent(in), optional :: unit

   integer :: iat, ndim

   real(wp), allocatable :: unitvec(:)
   real(wp), allocatable :: vvec(:), avec(:)
   real(wp) :: uvecsum
   real(wp) :: vvecsum
   real(wp) :: lambda

   logical :: grad, cpq, rebuilt, same_list
   logical :: add_lagr = .true.
   type(timer_type) :: timer
   integer :: print_unit, verbosity_solve

   ! Calculate gradient if the respective arrays are present
   grad = present(gradient) .and. present(sigma)
   cpq = present(dqdL) .and. (present(dqdr) .or. present(list))

   ! The augmented system of the direct solver is only set up as dense matrix
   if (.not. solver%need_pos_def .and. present(list)) then
      call fatal_error(error, "The direct solver does not support a neighborlist")
      return
   end if

   ! Charges, energy and gradient only use the upper triangle of a sorted list,
   ! the charge derivatives walk the complete rows
   if (present(list) .and. cpq) then
      if (.not. list%complete) then
         call fatal_error(error, "The charge derivatives require a sorted "// &
            & "complete list, see sort_neighborlist")
         return
      end if
   end if

   ! The CN derivatives of the charge derivatives are stored by the update
   if (cpq) then
      if (present(list)) then
         if (.not. allocated(cache%dcndrlist)) then
            call fatal_error(error, "Charge derivatives require the CN derivatives "// &
               & "on the neighborlist, update the cache with dqdr")
            return
         end if
         ! The pair arrays of the cache are formed on the list of the update
         same_list = size(cache%hinl) == size(list%inl)
         if (same_list) same_list = all(cache%hinl == list%inl)
         if (.not. same_list) then
            call fatal_error(error, "The cache was updated with a different "// &
               & "neighborlist")
            return
         end if
      else if (.not. allocated(cache%dcndr)) then
         call fatal_error(error, "Charge derivatives require the CN derivatives, "// &
            & "update the cache with dqdr")
         return
      end if
   end if

   if (.not. present(verbosity)) then
      verbosity_solve = 0
   else
      verbosity_solve = verbosity
   end if

   if (present(unit)) then
      print_unit = unit
   else
      print_unit = output_unit
   end if

   ! The CG solver requires a positive-definite system
   if (solver%need_pos_def) then
      ndim = mol%nat
      add_lagr = .false.
   else
      ndim = mol%nat + 1
      add_lagr = .true.
   end if

   call timer%push("total")
   call timer%push("setup")

   ! Setup the system matrices and vectors
   call self%get_capacitance_matrix(mol, ndim, cache, list)
   call self%get_coulomb_matrix(mol, ndim, cache, list)
   call self%get_xvec(mol, ndim, cache, list, efield)
   if (.not. allocated(cache%vrhs)) then
      allocate(cache%vrhs(mol%nat + 1))
   end if

   ! pop setup timer
   call timer%pop

   ! Print header
   call print_solve_header(print_unit, verbosity_solve, timer%get("setup"))

   if (add_lagr) then
      if (.not. allocated(cache%ainv)) then
         allocate(cache%ainv(ndim, ndim))
      end if
      cache%vrhs = cache%xvec
      cache%ainv = cache%amat
      call solver%solve(amat=cache%amat, xvec=cache%xvec, &
         & vrhs=cache%vrhs, ainv=cache%ainv, cpq=cpq, &
         & new_unit=print_unit, error=error)

   else
      if (.not. allocated(cache%uvec)) then
         allocate(cache%uvec(mol%nat))
      end if
      allocate(unitvec(mol%nat))
      allocate(vvec(mol%nat))
      ! Initial guess
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

      unitvec = 1.0_wp

      call print_constrained_system_message(print_unit, verbosity_solve, 'u')
      ! Constrained response: A*uvec = 1
      call solver%solve(amat=cache%amat, alist=cache%alist, xvec=unitvec, &
         & vrhs=cache%uvec, list=list, new_unit=print_unit, error=error)
      call print_constrained_system_message(print_unit, verbosity_solve, 'v')

      ! Unconstrained response: A*uvec = -xvec
      call solver%solve(amat=cache%amat, alist=cache%alist, xvec=-cache%xvec, &
         & vrhs=vvec, list=list, new_unit=print_unit, error=error)
      uvecsum = sum(cache%uvec)
      vvecsum = sum(vvec)
      ! Lagrangian multiplier
      lambda = - (mol%charge + vvecsum) / (uvecsum + eps)

      ! Projection of uvec on vvec
      cache%vrhs(:mol%nat) = -vvec - lambda * cache%uvec
      cache%vrhs(mol%nat + 1) = lambda

   end if

   ! Partial charges if present
   if (present(qvec)) then
      qvec(:) = cache%vrhs(:mol%nat)
   end if

   ! Electrostatic energy if present
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
      if (ndim > mol%nat) then
         ! Correct xvec to exclude constraint term
         cache%xvec(:mol%nat) = cache%xvec(:mol%nat) - 0.5_wp * cache%vrhs(mol%nat + 1)
      end if
      energy(:) = energy(:) + cache%vrhs(:mol%nat) * cache%xvec(:mol%nat)
      call timer%pop
      call print_energy_time(print_unit, verbosity_solve, timer%get("energy"))
   end if

   ! Calculate gradients if requested
   if (grad) then
      call timer%push("gradient")

      call self%get_grad(mol, cache, cache%vrhs(:mol%nat), gradient, sigma, &
         & alpha=0.5_wp, beta=-1.0_wp, list=list)

      ! pop gradient timer
      call timer%pop
      call print_gradient_time(print_unit, verbosity_solve, timer%get("gradient"))
   end if

   ! Calculate charge derivatives if requested
   if (cpq) then
      call timer%push("setup_gradient")

      ! Right-hand sides of the response equations, dX/dR - dA/dR*q
      call timer%push("dabdr_setup")
      call new_partial_derivs(mol, ndim, cache, list)
      call self%get_xvec_derivs(mol, ndim, cache, beta=1.0_wp, list=list)
      call self%get_coulomb_derivs(mol, ndim, cache, alpha=-1.0_wp, &
         & p=cache%vrhs(:mol%nat), list=list)
      call timer%pop
      if (verbosity_solve > 1) then
         write(output_unit, '(a, 1x, a)') &
            & "Partial derivatives setup time : ", &
            & format_time(timer%get("dabdr_setup"))
         write(output_unit, '(a)') ''
      end if

      ! pop gradient setup
      call timer%pop
      call print_gradient_header(print_unit, verbosity_solve, &
         & timer%get("setup_gradient"))

      call timer%push("cpq")
      call get_q_derivs(mol, solver, cache, error, ndim, dqdr, dqdL, list, &
         & unit=print_unit, verbosity=verbosity_solve)
      ! pop cpq timer
      call timer%pop
      call print_gradient_time(print_unit, verbosity_solve, timer%get("cpq"))
   end if

   ! pop total solve timer
   call timer%pop
   call print_total_time(print_unit, verbosity_solve, timer%get("total"))

end subroutine solve


!> Derivatives of the partial charges w.r.t. atomic positions and lattice vectors
subroutine get_q_derivs(mol, solver, cache, error, ndim, dqdr, dqdL, list, unit, &
   & verbosity)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> The solver instance
   class(mchrg_solver_type), intent(in) :: solver

   !> Cache handling
   type(mchrg_cache), intent(inout) :: cache

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Dimension of the linear system
   integer, intent(in) :: ndim

   !> Derivative of the atomic partial charges w.r.t. atomic positions
   real(wp), intent(out), optional :: dqdr(:, :, :)

   !> Derivative of the atomic partial charges w.r.t. lattice vectors
   real(wp), intent(out) :: dqdL(:, :, :)

   !> neighborlist optional type
   type(csr_list), intent(in), optional :: list

   !> Output unit
   integer, intent(in), optional :: unit

   !> Verbosity level
   integer, intent(in), optional :: verbosity

   real(wp), allocatable :: daqxdr(:, :, :), daqxdL(:, :, :), ainv(:, :)
   real(wp), allocatable :: dlambdadr(:, :)
   real(wp) :: scale, uvecsum
   integer :: jat, xat, bat, print_unit
   integer(i8) :: kat
   type(timer_type) :: timer

   ! Without a neighborlist the position derivatives are only formed densely
   if (.not. (present(list) .or. present(dqdr))) then
      call fatal_error(error, "Charge derivatives without a neighborlist require dqdr")
      return
   end if

   select type (solver)
   class is (cg_solver)
      if (present(list)) then
         call timer%push("inverse")
         call get_inverse_list(mol, solver, cache, list, error, unit)
         call timer%pop
         if (allocated(error)) return

         call timer%push("product")
         call get_dqdr_list(mol, cache, list)
         call get_constraint_factors(mol, cache, list, dlambdadr)
         call move_alloc(dlambdadr, cache%dlambdadr)
         call get_dqdL_list(mol, cache, list, dqdL)
         call timer%pop

         if (present(verbosity)) then
            if (verbosity > 0) then
               print_unit = output_unit
               if (present(unit)) print_unit = unit
               call print_stage_time(print_unit, "Inverse on the neighborlist", &
                  & timer%get("inverse"))
               call print_stage_time(print_unit, "Charge derivatives S*dB/dR", &
                  & timer%get("product"))
               write(print_unit, '(a)') ''
            end if
         end if

         ! dq_b/dR_x is held at the entry of b in the row of x
         if (present(dqdr)) then
            !$omp parallel do default(none) schedule(runtime) &
            !$omp shared(mol, cache, list, dqdr) private(xat, bat, kat)
            do xat = 1, mol%nat
               do bat = 1, mol%nat
                  dqdr(:, xat, bat) = -cache%dlambdadr(:, xat) * cache%uvec(bat)
               end do
               do kat = list%inl(xat), list%inl(xat + 1) - 1
                  bat = list%nlat(kat)
                  dqdr(:, xat, bat) = dqdr(:, xat, bat) + cache%dqdrlist(:, kat)
               end do
            end do
         end if
         return
      end if

      allocate(ainv(mol%nat, mol%nat))
      call solver%invert(amat=cache%amat, ainv=ainv, new_unit=unit, &
         & error=error)
      if (allocated(error)) return

      ! Projection onto the charge constraint, (P J^-1)^T = J^-1 - u*u^T/sum(u)
      uvecsum = sum(cache%uvec)
      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(mol, cache, ainv, uvecsum) private(jat, scale)
      do jat = 1, mol%nat
         scale = cache%uvec(jat) / (uvecsum + eps)
         ainv(:, jat) = ainv(:, jat) - scale * cache%uvec
      end do

      call gemm(cache%dabdr(:, :, :mol%nat), ainv, dqdr)
      call gemm(cache%dabdL(:, :, :mol%nat), ainv, dqdL)

   class default
      ! Non-iterative solve using the inverse of the augmented matrix
      if (.not. allocated(cache%ainv) .or. .not. present(dqdr)) then
         call fatal_error(error, &
            & "Charge derivatives require the inverse of the augmented matrix")
         return
      end if
      allocate(daqxdr(3, mol%nat, ndim), source=0.0_wp)
      allocate(daqxdL(3, 3, ndim), source=0.0_wp)
      daqxdr(:, :, :mol%nat) = cache%dabdr(:, :, :mol%nat)
      daqxdL(:, :, :mol%nat) = cache%dabdL(:, :, :mol%nat)
      call gemm(daqxdr, cache%ainv(:, :mol%nat), dqdr, alpha=1.0_wp)
      call gemm(daqxdL, cache%ainv(:, :mol%nat), dqdL, alpha=1.0_wp)
   end select

end subroutine get_q_derivs


!> Inverse of the coefficient matrix on a sorted complete neighborlist
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
      if (size(cache%ainvlist, kind=i8) /= nnz) deallocate(cache%ainvlist)
   end if
   if (.not. allocated(cache%ainvlist)) allocate(cache%ainvlist(nnz))
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


!> Add the transposed product of partial derivatives on a list with a vector
subroutine add_transposed_product(list, dlist, vec, wvec)

   !> Sorted complete neighborlist
   type(csr_list), intent(in) :: list

   !> Partial derivatives on the list
   real(wp), intent(in) :: dlist(:, :)

   !> Vector contracted with the components
   real(wp), intent(in) :: vec(:)

   !> Product for each atom, dimension (3, nat)
   real(wp), intent(inout) :: wvec(:, :)

   integer :: xat
   integer(i8) :: kat
   real(wp) :: wx(3)

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(list, dlist, vec, wvec) private(xat, kat, wx)
   do xat = 1, size(wvec, 2)
      wx(:) = 0.0_wp
      do kat = list%inl(xat), list%inl(xat + 1) - 1
         wx(:) = wx + dlist(:, kat) * vec(list%nlat(kat))
      end do
      wvec(:, xat) = wvec(:, xat) + wx
   end do

end subroutine add_transposed_product


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


!> Derivatives of the adjoint vector for a group of atoms
subroutine get_group_dpdr(cache, half, hfin, d2fdq2, d2fdq2_list, ufin_d2f, &
   & ng, nb, bloc, blist, qpos, qcol, pan, pcol)

   !> Cache with the inverse and the mixed derivatives on the list
   type(mchrg_cache), intent(in) :: cache

   !> Sorted list of the update with the inverse on its upper triangle
   type(csr_list), intent(in) :: half

   !> End of the upper triangle in each row of the list of the inverse
   integer(i8), intent(in), contiguous :: hfin(:)

   !> Second derivatives w.r.t. the charges on the upper triangle of their list
   real(wp), intent(in), contiguous :: d2fdq2(:)

   !> Sorted complete list of the second derivatives w.r.t. the charges
   type(csr_list), intent(in) :: d2fdq2_list

   !> End of the upper triangle in each row of d2fdq2_list
   integer(i8), intent(in), contiguous :: ufin_d2f(:)

   !> Number of atoms in the group
   integer, intent(in) :: ng

   !> Number of atoms in the union of the neighbors of the group
   integer, intent(in) :: nb

   !> Index of each atom in the union, zero outside
   integer, intent(in), contiguous :: bloc(:)

   !> Atoms of the union
   integer, intent(in), contiguous :: blist(:)

   !> Position of the entry of each union atom in the row of x
   integer(i8), intent(in), contiguous :: qpos(:, :)

   !> Charge derivatives dq_b/dR_x of the union
   real(wp), intent(in), contiguous :: qcol(:, :)

   !> Workspace for G + d2f/dq2 dq/dR of the union
   real(wp), intent(inout), contiguous :: pan(:, :)

   !> Derivatives dp_b/dR_x of the adjoint vector of the union
   real(wp), intent(inout), contiguous :: pcol(:, :)

   integer :: n, i
   integer(i8) :: kat

   ! G + d2f/dq2 dq/dR on the neighbors of x
   call multiply_group_list(d2fdq2_list, ufin_d2f, d2fdq2, nb, bloc, blist, &
      & qcol, pan)
   do n = 1, nb
      do i = 1, ng
         kat = qpos(i, n)
         if (kat > 0) then
            pan(3 * i - 2:3 * i, n) = pan(3 * i - 2:3 * i, n) &
               & + cache%dafdrlist(:, kat)
         else
            pan(3 * i - 2:3 * i, n) = 0.0_wp
         end if
      end do
   end do
   ! Product with the inverse on the neighbors of x
   call multiply_group_list(half, hfin, cache%ainvlist, nb, bloc, blist, pan, &
      & pcol)
   do n = 1, nb
      do i = 1, ng
         if (qpos(i, n) == 0) pcol(3 * i - 2:3 * i, n) = 0.0_wp
      end do
   end do

end subroutine get_group_dpdr


!> Union of the neighbors of a group of atoms on a sorted complete list
subroutine gather_group_list(list, dlist, group, nb, bloc, blist, qpos, pan)

   !> Sorted complete neighborlist
   type(csr_list), intent(in) :: list

   !> Partial derivatives on the list, the row of x holds the column x
   real(wp), intent(in), contiguous :: dlist(:, :)

   !> Atoms of the group
   integer, intent(in), contiguous :: group(:)

   !> Number of atoms in the union of the neighbors of the group
   integer, intent(out) :: nb

   !> Index of each atom in the union, zero outside
   integer, intent(inout), contiguous :: bloc(:)

   !> Atoms of the union, sorted in ascending order
   integer, intent(inout), contiguous :: blist(:)

   !> Position of the entry of each union atom in the row of x
   integer(i8), intent(inout), contiguous :: qpos(:, :)

   !> Partial derivatives dB_c/dR_x of the union
   real(wp), intent(inout), contiguous :: pan(:, :)

   integer :: i, ic, xat, bat, n
   integer(i8) :: kat

   nb = 0
   do i = 1, size(group)
      xat = group(i)
      ic = 3 * (i - 1)
      do kat = list%inl(xat), list%inl(xat + 1) - 1
         bat = list%nlat(kat)
         n = bloc(bat)
         if (n == 0) then
            nb = nb + 1
            n = nb
            bloc(bat) = n
            blist(n) = bat
            qpos(:, n) = 0_i8
            pan(:, n) = 0.0_wp
         end if
         pan(ic + 1:ic + 3, n) = dlist(:, kat)
         qpos(i, n) = kat
      end do
   end do

end subroutine gather_group_list


!> Product of a symmetric matrix on a list with the derivatives of a union
pure subroutine multiply_group_list(list, ufin, smat, nb, bloc, blist, pan, qcol)

   !> Sorted complete list of the symmetric matrix
   type(csr_list), intent(in) :: list

   !> End of the upper triangle in each row of the list
   integer(i8), intent(in), contiguous :: ufin(:)

   !> Symmetric matrix on the upper triangle of the list
   real(wp), intent(in), contiguous :: smat(:)

   !> Number of atoms in the union
   integer, intent(in) :: nb

   !> Index of each atom in the union, zero outside
   integer, intent(in), contiguous :: bloc(:)

   !> Atoms of the union
   integer, intent(in), contiguous :: blist(:)

   !> Derivatives of the union with fixed leading dimension 3 * max_group
   real(wp), intent(in) :: pan(3 * max_group, *)

   !> Product with the symmetric matrix on the union
   real(wp), intent(inout) :: qcol(3 * max_group, *)

   integer, parameter :: nvec = 3 * max_group
   integer :: n, m, bat
   integer(i8) :: hat
   real(wp) :: s, pn(nvec), acc(nvec)

   qcol(:, :nb) = 0.0_wp
   do n = 1, nb
      bat = blist(n)
      pn(:) = pan(:, n)
      acc(:) = qcol(:, n) + smat(list%inl(bat)) * pn(:)
      do hat = list%inl(bat) + 1, ufin(bat)
         m = bloc(list%nlat(hat))
         if (m == 0) cycle
         s = smat(hat)
         acc(:) = acc(:) + s * pan(:, m)
         qcol(:, m) = qcol(:, m) + s * pn(:)
      end do
      qcol(:, n) = acc(:)
   end do

end subroutine multiply_group_list


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
   if (.not. ext) then
      if (allocated(cache%dabdL)) deallocate(cache%dabdL)
      allocate(cache%dabdL(3, 3, ndim), source=0.0_wp)
   end if

end subroutine new_partial_derivs


!> Add partial derivatives on another sorted complete list
subroutine add_list_derivs(nat, list, dlist, plist, pval)

   !> Number of atoms
   integer, intent(in) :: nat

   !> Sorted complete neighborlist
   type(csr_list), intent(in) :: list

   !> Partial derivatives on the neighborlist
   real(wp), intent(inout) :: dlist(:, :)

   !> Sorted complete list of the added partial derivatives
   type(csr_list), intent(in) :: plist

   !> Added partial derivatives on their list
   real(wp), intent(in) :: pval(:, :)

   integer :: xat
   integer(i8) :: kat, pat
   integer(i8), allocatable :: loc(:)

   ! loc holds the position of each atom in the row of the neighborlist
   !$omp parallel default(none) shared(nat, list, dlist, plist, pval) &
   !$omp private(xat, kat, pat, loc)
   allocate(loc(nat), source=0_i8)
   !$omp do schedule(runtime)
   do xat = 1, nat
      do kat = list%inl(xat), list%inl(xat + 1) - 1
         loc(list%nlat(kat)) = kat
      end do
      do pat = plist%inl(xat), plist%inl(xat + 1) - 1
         kat = loc(plist%nlat(pat))
         if (kat == 0) cycle
         dlist(:, kat) = dlist(:, kat) + pval(:, pat)
      end do
      do kat = list%inl(xat), list%inl(xat + 1) - 1
         loc(list%nlat(kat)) = 0
      end do
   end do
   !$omp end do
   deallocate(loc)
   !$omp end parallel

end subroutine add_list_derivs


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

   ! Get variables from the model cache
   if (.not. allocated(cache%amat)) then
      call fatal_error(error, "J-matrix is not allocated")
      return
   end if
   if (.not. allocated(cache%uvec) .and. solver%need_pos_def) then
      call fatal_error(error, "Constrained response J*uvec = 1 is not allocated")
      return
   end if

   call get_adjoint_vector(mol, solver, cache, error, dfdq, padj, list, &
      & print_unit, verbosity_solve)
   if (allocated(error)) return

   ! Evaluate external gradients via adjoint contraction

   call self%get_grad(mol, cache, padj(:mol%nat), dfdr, dfdL, alpha=-1.0_wp, &
      & beta=1.0_wp, list=list)

   ! pop dfdr timer
   call timer%pop
   call print_gradient_time(print_unit, verbosity_solve, &
      & timer%get("external_gradient"))

end subroutine get_external_gradient


!> Compute the molecular nuclear Hessian of the electrostatic energy
subroutine get_hessian_dense(self, mol, solver, cache, error, hessian, &
   & unit, verbosity)

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

   !> Output unit
   integer, intent(in), optional :: unit

   !> Verbosity level
   integer, intent(in), optional :: verbosity

   real(wp), allocatable :: dqdr(:, :, :), dqdL(:, :, :)

   if (any(mol%periodic)) then
      call fatal_error(error, &
         & "Analytical Hessian is not available for periodic systems")
      return
   end if

   ! Get variables from the model cache
   if (.not. allocated(cache%vrhs) .or. .not. allocated(cache%dcndr)) then
      call fatal_error(error, "Charges and CN derivatives are not in the cache, "// &
         & "update the cache with hess and solve first")
      return
   end if

   ! Charge derivatives
   call new_partial_derivs(mol, mol%nat, cache)
   call self%get_xvec_derivs(mol, mol%nat, cache, beta=1.0_wp)
   call self%get_coulomb_derivs(mol, mol%nat, cache, alpha=-1.0_wp, &
      & p=cache%vrhs(:mol%nat))
   allocate(dqdL(3, 3, mol%nat), dqdr(3, mol%nat, mol%nat))
   call get_q_derivs(mol, solver, cache, error, mol%nat, dqdr, dqdL, &
      & unit=unit, verbosity=verbosity)
   if (allocated(error)) return

   ! Second derivatives of 0.5 * q^T A q - q^T X at fixed charges
   hessian(:, :) = 0.0_wp
   call self%get_explicit_hessian(mol, cache, cache%vrhs(:mol%nat), hessian, &
      & error, alpha=0.5_wp, beta=-1.0_wp)
   if (allocated(error)) return
   call self%get_implicit_hessian(mol, cache, hessian, dqdr=dqdr)

end subroutine get_hessian_dense


!> Compute the molecular nuclear Hessian of the electrostatic energy on a list
subroutine get_hessian_list(self, mol, solver, cache, error, hessian, list, &
   & unit, verbosity)

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

   !> Sorted complete neighborlist, walked along the complete rows
   type(csr_list), intent(in) :: list

   !> Output unit
   integer, intent(in), optional :: unit

   !> Verbosity level
   integer, intent(in), optional :: verbosity

   real(wp) :: tgroup(2)
   type(timer_type) :: timer
   integer :: print_unit, verbosity_hess

   ! Without a list of its own the blocks stay on the complete list
   if (.not. allocated(hessian%list)) hessian%list = list
   call check_hessian_setup(self, mol, solver, cache, hessian, list, error)
   if (allocated(error)) return

   print_unit = output_unit
   if (present(unit)) print_unit = unit
   verbosity_hess = 0
   if (present(verbosity)) verbosity_hess = verbosity

   call timer%push("total")

   ! Partial derivatives dX/dR - dA/dR*q and the inverse on the list
   call print_hessian_header(print_unit, verbosity_hess)
   call timer%push("dabdr")
   call new_partial_derivs(mol, mol%nat, cache, list)
   call self%get_xvec_derivs(mol, mol%nat, cache, beta=1.0_wp, list=list)
   call self%get_coulomb_derivs(mol, mol%nat, cache, alpha=-1.0_wp, &
      & p=cache%vrhs(:mol%nat), list=list)
   call timer%pop
   if (allocated(cache%tpos)) deallocate(cache%tpos)
   call timer%push("inverse")
   select type (solver)
   class is (cg_solver)
      call get_inverse_list(mol, solver, cache, list, error, unit, hessian%list)
   end select
   call timer%pop
   if (allocated(error)) return
   deallocate(cache%alist)

   call new_hessian_blocks(hessian)

   ! Response of the charges, -(dX/dR - dA/dR*q)^T dq/dR
   call timer%push("implicit")
   call self%get_implicit_hessian(mol, cache, hessian, list, times=tgroup)
   call timer%pop
   deallocate(cache%ainvlist, cache%dabdrlist)

   ! Second derivatives of 0.5 * q^T A q - q^T X at fixed charges
   call timer%push("explicit")
   call self%get_explicit_hessian(mol, cache, cache%vrhs(:mol%nat), hessian, &
      & error, list, alpha=0.5_wp, beta=-1.0_wp)
   call timer%pop
   call timer%pop
   if (allocated(error)) return

   if (verbosity_hess > 0) call print_hessian_stages(print_unit, timer, tgroup)

end subroutine get_hessian_list


!> Allocate the zero blocks of a compressed Hessian on the upper triangle of its list
subroutine new_hessian_blocks(hessian, list)

   !> Molecular nuclear Hessian in compressed storage
   type(mchrg_hessian), intent(inout) :: hessian

   !> Sorted list of the blocks, absent to keep the list of the Hessian
   type(csr_list), intent(in), optional :: list

   if (present(list)) hessian%list = list
   if (allocated(hessian%hesslist)) deallocate(hessian%hesslist)
   allocate(hessian%hesslist(3, 3, size(hessian%list%nlat, kind=i8)), source=0.0_wp)

end subroutine new_hessian_blocks


!> Check the setup of the Hessian on a neighborlist
subroutine check_hessian_setup(self, mol, solver, cache, hessian, list, error)

   !> Electronegativity-equilibration model
   class(mchrg_model_type), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Solver instance
   class(mchrg_solver_type), intent(in) :: solver

   !> Cache updated with the list of the Hessian
   type(mchrg_cache), intent(in) :: cache

   !> Molecular nuclear Hessian with its list
   type(mchrg_hessian), intent(in) :: hessian

   !> Sorted complete neighborlist
   type(csr_list), intent(in) :: list

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iat
   integer(i8) :: nup
   logical :: same

   if (any(mol%periodic)) then
      call fatal_error(error, &
         & "Analytical Hessian is not available for periodic systems")
      return
   end if
   if (.not. list%complete) then
      call fatal_error(error, "The Hessian on a neighborlist requires a sorted "// &
         & "complete list, see sort_neighborlist")
      return
   end if
   select type (solver)
   class is (cg_solver)
   class default
      call fatal_error(error, "The Hessian on a neighborlist requires the CG solver")
      return
   end select
   if (.not. allocated(cache%vrhs) .or. .not. allocated(cache%dcndrlist)) then
      call fatal_error(error, "Charges and CN derivatives on the list are not "// &
         & "in the cache, update the cache with hess and solve first")
      return
   end if

   ! The mirrored CN derivatives of the pair arrays hold only without CN cap
   if (self%ncoord%cut > 0.0_wp) then
      call fatal_error(error, "The Hessian on a neighborlist requires the "// &
         & "coordination numbers without cap")
      return
   end if
   if (allocated(self%ncoord_en)) then
      if (self%ncoord_en%cut > 0.0_wp) then
         call fatal_error(error, "The Hessian on a neighborlist requires the "// &
            & "local charges without cap")
         return
      end if
   end if

   ! The cache holds the pair arrays on the rows of the list of the Hessian
   if (size(list%inl) /= mol%nat + 1 .or. .not. allocated(cache%hinl)) then
      call fatal_error(error, "The cache was updated with a different neighborlist")
      return
   end if
   if (size(hessian%list%inl) /= mol%nat + 1 .or. size(cache%hinl) /= mol%nat + 1) then
      call fatal_error(error, "The cache was updated with a different neighborlist")
      return
   end if
   if (any(cache%hinl /= hessian%list%inl) &
      & .or. size(cache%dcndrlist, 2, kind=i8) /= size(hessian%list%nlat, kind=i8)) &
      & then
      call fatal_error(error, "The cache was updated with a different neighborlist")
      return
   end if

   ! Each row of the list of the Hessian starts with the upper triangle of the
   ! same row of the complete list, a half list holds nothing else
   same = .true.
   !$omp parallel do default(none) schedule(runtime) reduction(.and.:same) &
   !$omp shared(mol, list, hessian) private(iat, nup)
   do iat = 1, mol%nat
      nup = find_lower_start(list, iat) - list%inl(iat)
      if (find_lower_start(hessian%list, iat) - hessian%list%inl(iat) /= nup) then
         same = .false.
      else
         same = same .and. all(hessian%list%nlat(hessian%list%inl(iat): &
            & hessian%list%inl(iat) + nup - 1) &
            & == list%nlat(list%inl(iat):list%inl(iat) + nup - 1))
      end if
   end do
   if (.not. same) then
      call fatal_error(error, "The list of the Hessian is not the upper triangle "// &
         & "of the complete neighborlist")
      return
   end if

end subroutine check_hessian_setup


!> Compute the implicit Hessian of the electrostatic energy
subroutine get_implicit_hessian_dense(self, mol, cache, hessian, dqdr)

   !> Electronegativity-equilibration model
   class(mchrg_model_type), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Cache handling
   type(mchrg_cache), intent(in) :: cache

   !> Molecular nuclear Hessian
   real(wp), intent(inout) :: hessian(:, :)

   !> Derivatives of the charges w.r.t. positions, dqdr(:, x, b) = dq_b/dR_x
   real(wp), intent(in) :: dqdr(:, :, :)

   call gemm(reshape(cache%dabdr(:, :, :mol%nat), [3 * mol%nat, mol%nat]), &
      & reshape(dqdr, [3 * mol%nat, mol%nat]), hessian, transb='t', &
      & alpha=-1.0_wp, beta=1.0_wp)

end subroutine get_implicit_hessian_dense


!> Compute the implicit Hessian of the electrostatic energy on a list
subroutine get_implicit_hessian_list(self, mol, cache, hessian, list, times)

   !> Electronegativity-equilibration model
   class(mchrg_model_type), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Cache handling
   type(mchrg_cache), intent(in) :: cache

   !> Molecular nuclear Hessian in compressed storage
   type(mchrg_hessian), intent(inout) :: hessian

   !> Multicharge neighborlist type (sorted complete)
   type(csr_list), intent(in) :: list

   !> Times of the group derivatives and their contraction
   real(wp), intent(out), optional :: times(2)

   real(wp) :: tgroup(2)
   integer, allocatable :: order(:)

   call get_spatial_order(mol%xyz, order)
   call add_implicit_blocks_list(mol, cache, hessian, list, order, tgroup)
   if (present(times)) times(:) = tgroup

   call get_constraint_factors(mol, cache, list, hessian%dlambda)
   hessian%weight = sum(cache%uvec)
   if (allocated(hessian%wvec)) deallocate(hessian%wvec)

end subroutine get_implicit_hessian_list


!> Add the implicit blocks for groups of atoms on a sorted complete list
subroutine add_implicit_blocks_list(mol, cache, hessian, list, order, times, &
   & d2fdq2, d2fdq2_list)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Cache with the inverse, partial and mixed derivatives
   type(mchrg_cache), intent(in) :: cache

   !> Molecular nuclear Hessian with allocated blocks on its list, the list of
   !> the update with the inverse on its upper triangle
   type(mchrg_hessian), intent(inout) :: hessian

   !> Sorted complete neighborlist of the partial derivatives
   type(csr_list), intent(in) :: list

   !> Atoms in the spatial order of the groups
   integer, intent(in), contiguous :: order(:)

   !> Times of the group derivatives and their contraction
   real(wp), intent(out) :: times(2)

   !> Second derivatives of an external function w.r.t. the charges
   real(wp), intent(in), contiguous, optional :: d2fdq2(:)

   !> Sorted complete list of the second derivatives w.r.t. the charges
   type(csr_list), intent(in), optional :: d2fdq2_list

   integer :: igrp, ngrp, nbmax, numax, x0, ng, i, a, na, xat, yat, nb, nu, n, u
   integer :: act(max_group)
   integer(i8) :: hat, tick(3), rate
   real(wp) :: scale, tpan, tcon
   real(wp) :: blk(3, 3, max_group)
   logical :: ext
   integer, allocatable :: bloc(:), hloc(:), blist(:), ulist(:)
   integer(i8), allocatable :: hfin(:), ufin_d2f(:), qpos(:, :), jmap(:, :)
   real(wp), allocatable :: pan(:, :), qcol(:, :), pcol(:, :)

   times(:) = 0.0_wp
   ngrp = (mol%nat + max_group - 1) / max_group
   if (ngrp < 1) return
   nbmax = min(mol%nat, max_group * get_max_neighbors(mol%nat, list))
   numax = nbmax
   ext = present(d2fdq2) .and. present(d2fdq2_list)
   scale = -1.0_wp
   if (ext) scale = 1.0_wp

   ! Rows of the blocks and of the inverse
   allocate(hfin(mol%nat))
   do i = 1, mol%nat
      hfin(i) = find_lower_start(hessian%list, i) - 1_i8
   end do
   if (ext) then
      allocate(ufin_d2f(mol%nat))
      do i = 1, mol%nat
         ufin_d2f(i) = find_lower_start(d2fdq2_list, i) - 1_i8
      end do
   end if

   ! bloc and hloc hold the position of each atom in the unions of the group
   tpan = 0.0_wp
   tcon = 0.0_wp
   !$omp parallel default(none) &
   !$omp shared(mol, cache, hessian, list, hfin, ufin_d2f, order, ngrp, nbmax) &
   !$omp shared(numax, ext, scale, d2fdq2, d2fdq2_list) &
   !$omp private(igrp, x0, ng, i, a, na, act, xat, yat, nb, nu, n, u, hat, tick, rate) &
   !$omp private(bloc, hloc, blist, ulist, qpos, jmap, pan, qcol, pcol, blk) &
   !$omp reduction(+:tpan, tcon)
   allocate(bloc(mol%nat), hloc(mol%nat), source=0)
   allocate(blist(nbmax), ulist(numax), qpos(max_group, nbmax), &
      & jmap(max_group, numax), pan(3 * max_group, nbmax), &
      & qcol(3 * max_group, nbmax))
   if (ext) allocate(pcol(3 * max_group, nbmax))
   !$omp do schedule(runtime)
   do igrp = 1, ngrp
      call system_clock(tick(1), rate)
      x0 = max_group * (igrp - 1)
      ng = min(max_group, mol%nat - x0)

      call get_group_dqdr(cache, list, hessian%list, hfin, order(x0 + 1:x0 + ng), &
         & nb, bloc, blist, qpos, pan, qcol)
      if (ext) call get_group_dpdr(cache, hessian%list, hfin, d2fdq2, d2fdq2_list, &
         & ufin_d2f, ng, nb, bloc, blist, qpos, qcol, pan, pcol)
      call system_clock(tick(2))

      nu = 0
      do i = 1, ng
         xat = order(x0 + i)
         do hat = hessian%list%inl(xat), hfin(xat)
            yat = hessian%list%nlat(hat)
            if (hloc(yat) == 0) then
               nu = nu + 1
               hloc(yat) = nu
               ulist(nu) = yat
               jmap(:ng, nu) = 0
            end if
            jmap(i, hloc(yat)) = hat
         end do
      end do

      do u = 1, nu
         yat = ulist(u)
         hloc(yat) = 0
         na = 0
         do i = 1, ng
            hat = jmap(i, u)
            if (hat == 0) cycle
            na = na + 1
            act(na) = i
         end do
         if (ext) then
            call add_column_products_ext(list, cache%dabdrlist, cache%dafdrlist, &
               & yat, bloc, pcol, qcol, na, act, blk)
         else
            call add_column_products(list, cache%dabdrlist, yat, bloc, qcol, &
               & na, act, blk)
         end if
         do a = 1, na
            i = act(a)
            hat = jmap(i, u)
            if (yat == order(x0 + i)) then
               hessian%hesslist(:, :, hat) = hessian%hesslist(:, :, hat) &
                  & + 0.5_wp * scale * (blk(:, :, a) + transpose(blk(:, :, a)))
            else
               hessian%hesslist(:, :, hat) = hessian%hesslist(:, :, hat) &
                  & + scale * blk(:, :, a)
            end if
         end do
      end do
      do n = 1, nb
         bloc(blist(n)) = 0
      end do
      call system_clock(tick(3))
      tpan = tpan + real(tick(2) - tick(1), wp) / real(rate, wp)
      tcon = tcon + real(tick(3) - tick(2), wp) / real(rate, wp)
   end do
   !$omp end do
   deallocate(bloc, hloc, blist, ulist, qpos, jmap, pan, qcol)
   if (ext) deallocate(pcol)
   !$omp end parallel
   times(:) = [tpan, tcon]

end subroutine add_implicit_blocks_list


!> Compute the products of a column of partial derivatives with group derivatives
pure subroutine add_column_products(list, dlist, yat, bloc, col, na, act, blk)

   !> Sorted complete neighborlist
   type(csr_list), intent(in) :: list

   !> Partial derivatives on the list, the row of y holds the column y
   real(wp), intent(in), contiguous :: dlist(:, :)

   !> Atom of the column
   integer, intent(in) :: yat

   !> Index of each atom in the union, zero outside
   integer, intent(in), contiguous :: bloc(:)

   !> Derivatives of the group for the union with leading dimension 3 * max_group
   real(wp), intent(in) :: col(3 * max_group, *)

   !> Number of group atoms with an upper-triangle block with y
   integer, intent(in) :: na

   !> Indices of the active group atoms
   integer, intent(in) :: act(max_group)

   !> 3x3 blocks for the active group atoms
   real(wp), intent(out) :: blk(3, 3, max_group)

   integer :: a, ic, n
   integer(i8) :: kat
   real(wp) :: dy(3)

   blk(:, :, :na) = 0.0_wp
   if (na == max_group) then
      do kat = list%inl(yat), list%inl(yat + 1) - 1
         n = bloc(list%nlat(kat))
         if (n == 0) cycle
         dy(:) = dlist(:, kat)
         do a = 1, max_group
            ic = 3 * (a - 1)
            blk(:, 1, a) = blk(:, 1, a) + dy * col(ic + 1, n)
            blk(:, 2, a) = blk(:, 2, a) + dy * col(ic + 2, n)
            blk(:, 3, a) = blk(:, 3, a) + dy * col(ic + 3, n)
         end do
      end do
   else
      do kat = list%inl(yat), list%inl(yat + 1) - 1
         n = bloc(list%nlat(kat))
         if (n == 0) cycle
         dy(:) = dlist(:, kat)
         do a = 1, na
            ic = 3 * (act(a) - 1)
            blk(:, 1, a) = blk(:, 1, a) + dy * col(ic + 1, n)
            blk(:, 2, a) = blk(:, 2, a) + dy * col(ic + 2, n)
            blk(:, 3, a) = blk(:, 3, a) + dy * col(ic + 3, n)
         end do
      end do
   end if

end subroutine add_column_products


!> Compute the fused products of two columns of partial derivatives with group derivatives
pure subroutine add_column_products_ext(list, dlist1, dlist2, yat, bloc, col1, &
   & col2, na, act, blk)

   !> Sorted complete neighborlist
   type(csr_list), intent(in) :: list

   !> First partial derivatives on the list
   real(wp), intent(in), contiguous :: dlist1(:, :)

   !> Second partial derivatives on the list
   real(wp), intent(in), contiguous :: dlist2(:, :)

   !> Atom of the column
   integer, intent(in) :: yat

   !> Index of each atom in the union, zero outside
   integer, intent(in), contiguous :: bloc(:)

   !> First derivatives of the group for the union
   real(wp), intent(in) :: col1(3 * max_group, *)

   !> Second derivatives of the group for the union
   real(wp), intent(in) :: col2(3 * max_group, *)

   !> Number of group atoms with an upper-triangle block with y
   integer, intent(in) :: na

   !> Indices of the active group atoms
   integer, intent(in) :: act(max_group)

   !> 3x3 blocks for the active group atoms
   real(wp), intent(out) :: blk(3, 3, max_group)

   integer :: a, ic, n
   integer(i8) :: kat
   real(wp) :: dy1(3), dy2(3)

   blk(:, :, :na) = 0.0_wp
   do kat = list%inl(yat), list%inl(yat + 1) - 1
      n = bloc(list%nlat(kat))
      if (n == 0) cycle
      dy1(:) = dlist1(:, kat)
      dy2(:) = dlist2(:, kat)
      do a = 1, na
         ic = 3 * (act(a) - 1)
         blk(:, 1, a) = blk(:, 1, a) + dy1 * col1(ic + 1, n) + dy2 * col2(ic + 1, n)
         blk(:, 2, a) = blk(:, 2, a) + dy1 * col1(ic + 2, n) + dy2 * col2(ic + 2, n)
         blk(:, 3, a) = blk(:, 3, a) + dy1 * col1(ic + 3, n) + dy2 * col2(ic + 3, n)
      end do
   end do

end subroutine add_column_products_ext


!> Compute the nuclear Hessian of an external function of the charges
subroutine get_external_hessian_dense(self, mol, solver, cache, error, dfdq, &
   & d2fdq2, dfdrdq, hessian, unit, verbosity)

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

   !> Output unit
   integer, intent(in), optional :: unit

   !> Verbosity level
   integer, intent(in), optional :: verbosity

   real(wp), allocatable :: dqdr(:, :, :), dqdL(:, :, :), padj(:)
   integer :: ndim

   ndim = 3 * mol%nat
   if (size(dfdq) /= mol%nat .or. any(shape(d2fdq2) /= [mol%nat, mol%nat]) &
      & .or. any(shape(dfdrdq) /= [3, mol%nat, mol%nat]) &
      & .or. any(shape(hessian) /= [ndim, ndim])) then
      call fatal_error(error, "External partial derivatives are wrong size")
      return
   end if

   ! Get variables from the model cache
   if (.not. allocated(cache%vrhs) .or. .not. allocated(cache%dcndr)) then
      call fatal_error(error, "Charges and CN derivatives are not in the cache, "// &
         & "update the cache with hess and solve first")
      return
   end if

   ! Charge derivatives
   call new_partial_derivs(mol, mol%nat, cache)
   call self%get_xvec_derivs(mol, mol%nat, cache, beta=1.0_wp)
   call self%get_coulomb_derivs(mol, mol%nat, cache, alpha=-1.0_wp, &
      & p=cache%vrhs(:mol%nat))
   allocate(dqdL(3, 3, mol%nat), dqdr(3, mol%nat, mol%nat))
   call get_q_derivs(mol, solver, cache, error, mol%nat, dqdr, dqdL, &
      & unit=unit, verbosity=verbosity)
   if (allocated(error)) return

   call get_adjoint_vector(mol, solver, cache, error, dfdq, padj, unit=unit, &
      & verbosity=verbosity)
   if (allocated(error)) return

   ! Mixed derivatives of the Lagrangian, G^T = d2f/dRdq - dA/dR*p
   call new_partial_derivs(mol, size(padj), cache, dafdr=.true.)
   call self%get_coulomb_derivs(mol, size(padj), cache, alpha=-1.0_wp, &
      & p=padj(:mol%nat), dafdr=.true.)
   cache%dafdr = cache%dafdr + dfdrdq

   call gemm(reshape(cache%dafdr, [ndim, mol%nat]), &
      & reshape(dqdr, [ndim, mol%nat]), hessian, transb='t')
   call gemm(dqdr, d2fdq2, cache%dafdr, beta=1.0_wp)
   call gemm(reshape(dqdr, [ndim, mol%nat]), &
      & reshape(cache%dafdr, [ndim, mol%nat]), hessian, transb='t', &
      & beta=1.0_wp)
   deallocate(cache%dafdr)

   ! Second derivatives of p^T (X - A*q) at fixed charges and adjoint vector
   call self%get_explicit_hessian(mol, cache, padj(:mol%nat), hessian, &
      & error, alpha=-1.0_wp, beta=1.0_wp)
   if (allocated(error)) return

end subroutine get_external_hessian_dense


!> Compute the nuclear Hessian of an external function on a list
subroutine get_external_hessian_list(self, mol, solver, cache, error, dfdq, &
   & d2fdq2, d2fdq2_list, dfdrdq, dfdrdq_list, hessian, list, unit, verbosity)

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

   !> Neighborlist type
   type(csr_list), intent(in) :: list

   !> Output unit
   integer, intent(in), optional :: unit

   !> Verbosity level
   integer, intent(in), optional :: verbosity

   integer :: print_unit, verbosity_hess
   real(wp) :: tgroup(2)
   integer, allocatable :: order(:)
   real(wp), allocatable :: padj(:), fuvec(:), tvec(:)
   type(timer_type) :: timer

   ! Without a list of its own the blocks stay on the complete list
   if (.not. allocated(hessian%list)) hessian%list = list
   call check_hessian_setup(self, mol, solver, cache, hessian, list, error)
   if (allocated(error)) return
   call check_external_derivs(mol, dfdq, d2fdq2, d2fdq2_list, dfdrdq, &
      & dfdrdq_list, error)
   if (allocated(error)) return

   print_unit = output_unit
   if (present(unit)) print_unit = unit
   verbosity_hess = 0
   if (present(verbosity)) verbosity_hess = verbosity

   call timer%push("total")

   ! Adjoint vector p = S df/dq and the response of the charge constraint
   call timer%push("adjoint")
   ! The solves with J only read its upper triangle on the list of the update
   call get_adjoint_vector(mol, solver, cache, error, dfdq, padj, hessian%list, &
      & unit, verbosity)
   if (.not. allocated(error)) call get_constraint_response(mol, solver, cache, &
      & hessian%list, d2fdq2, d2fdq2_list, fuvec, tvec, print_unit, error)
   call timer%pop
   if (allocated(error)) return

   call print_hessian_header(print_unit, verbosity_hess)
   call timer%push("dabdr")
   call new_partial_derivs(mol, mol%nat, cache, list)
   call self%get_xvec_derivs(mol, mol%nat, cache, beta=1.0_wp, list=list)
   call self%get_coulomb_derivs(mol, mol%nat, cache, alpha=-1.0_wp, &
      & p=cache%vrhs(:mol%nat), list=list)
   call new_partial_derivs(mol, mol%nat, cache, list, dafdr=.true.)
   call self%get_coulomb_derivs(mol, mol%nat, cache, alpha=-1.0_wp, p=padj, &
      & list=list, dafdr=.true.)
   call add_list_derivs(mol%nat, list, cache%dafdrlist, dfdrdq_list, dfdrdq)
   call timer%pop
   if (allocated(cache%tpos)) deallocate(cache%tpos)
   call timer%push("inverse")
   select type (solver)
   class is (cg_solver)
      call get_inverse_list(mol, solver, cache, list, error, unit, hessian%list)
   end select
   call timer%pop
   if (allocated(error)) return
   deallocate(cache%alist)

   call new_hessian_blocks(hessian)
   call get_external_constraint_factors(mol, cache, list, fuvec, tvec, hessian)

   call timer%push("implicit")
   call get_spatial_order(mol%xyz, order)
   call add_implicit_blocks_list(mol, cache, hessian, list, order, tgroup, &
      & d2fdq2, d2fdq2_list)
   call timer%pop
   deallocate(cache%ainvlist, cache%dabdrlist, cache%dafdrlist)

   ! Second derivatives of p^T (X - A*q) at fixed charges and adjoint vector
   call timer%push("explicit")
   call self%get_explicit_hessian(mol, cache, padj, hessian, error, list, &
      & alpha=-1.0_wp, beta=1.0_wp)
   call timer%pop
   call timer%pop
   if (allocated(error)) return

   if (verbosity_hess > 0) call print_hessian_stages(print_unit, timer, tgroup, &
      & adjoint=.true.)

end subroutine get_external_hessian_list


!> Check the derivatives of an external function on their lists
subroutine check_external_derivs(mol, dfdq, d2fdq2, d2fdq2_list, dfdrdq, &
   & dfdrdq_list, error)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Derivative of the external function w.r.t. atomic charges
   real(wp), intent(in) :: dfdq(:)

   !> Second derivative w.r.t. atomic charges on the upper triangle of its list
   real(wp), intent(in) :: d2fdq2(:)

   !> Sorted complete list of the second derivative w.r.t. atomic charges
   type(csr_list), intent(in) :: d2fdq2_list

   !> Mixed derivative w.r.t. positions and charges on its list
   real(wp), intent(in) :: dfdrdq(:, :)

   !> Sorted complete list of the mixed derivative w.r.t. positions and charges
   type(csr_list), intent(in) :: dfdrdq_list

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   if (.not. (d2fdq2_list%complete .and. dfdrdq_list%complete)) then
      call fatal_error(error, "The external Hessian on a neighborlist requires "// &
         & "sorted complete lists")
      return
   end if

   if (size(dfdq) /= mol%nat &
      & .or. size(d2fdq2, kind=i8) /= size(d2fdq2_list%nlat, kind=i8) &
      & .or. size(dfdrdq, 1) /= 3 &
      & .or. size(dfdrdq, 2, kind=i8) /= size(dfdrdq_list%nlat, kind=i8)) then
      call fatal_error(error, "External partial derivatives are wrong size")
      return
   end if

end subroutine check_external_derivs


!> Response of the charge constraint to an external function
subroutine get_constraint_response(mol, solver, cache, list, d2fdq2, &
   & d2fdq2_list, fuvec, tvec, unit, error)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Solver instance
   class(mchrg_solver_type), intent(in) :: solver

   !> Cache with the coefficient matrix on the list and the constraint response
   type(mchrg_cache), intent(in) :: cache

   !> Sorted list of the update with the coefficient matrix on its upper triangle
   type(csr_list), intent(in) :: list

   !> Second derivative w.r.t. atomic charges on the upper triangle of its list
   real(wp), intent(in), contiguous :: d2fdq2(:)

   !> Sorted complete list of the second derivative w.r.t. atomic charges
   type(csr_list), intent(in) :: d2fdq2_list

   !> Product F u
   real(wp), allocatable, intent(out) :: fuvec(:)

   !> Response t = A^-1 F u
   real(wp), allocatable, intent(out) :: tvec(:)

   !> Output unit
   integer, intent(in) :: unit

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iat

   allocate(fuvec(mol%nat), tvec(mol%nat))
   call spsymv_csr(mol%nat, d2fdq2, d2fdq2_list%inl, d2fdq2_list%nlat, cache%uvec, &
      & fuvec, uplo='U')

   ! Diagonal guess of the CG solve
   do iat = 1, mol%nat
      tvec(iat) = fuvec(iat) / (cache%alist(list%inl(iat)) + eps)
   end do
   call solver%solve(alist=cache%alist, xvec=fuvec, vrhs=tvec, list=list, &
      & new_unit=unit, error=error)

end subroutine get_constraint_response


!> Charge constraint factors of the Hessian of an external function
subroutine get_external_constraint_factors(mol, cache, list, fuvec, tvec, &
   & hessian)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Cache with the partial and mixed derivatives
   type(mchrg_cache), intent(in) :: cache

   !> Sorted complete neighborlist
   type(csr_list), intent(in) :: list

   !> Product F u of the second derivatives w.r.t. the charges
   real(wp), intent(in) :: fuvec(:)

   !> Response t = A^-1 F u
   real(wp), intent(in) :: tvec(:)

   !> Molecular nuclear Hessian, receives the factors
   type(mchrg_hessian), intent(inout) :: hessian

   call get_constraint_factors(mol, cache, list, hessian%dlambda)
   if (allocated(hessian%wvec)) deallocate(hessian%wvec)
   allocate(hessian%wvec(3, mol%nat), source=0.0_wp)
   call add_transposed_product(list, cache%dafdrlist, cache%uvec, hessian%wvec)
   call add_transposed_product(list, cache%dabdrlist, tvec, hessian%wvec)
   hessian%weight = dot_product(cache%uvec, fuvec)

end subroutine get_external_constraint_factors


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


!> Semi-numerical Hessian and pressure tensor
subroutine get_numhess_packed(self, mol, solver, cache, error, qvec, energy, grad, &
   & sigma, hess, press, list, unit, verbosity)

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

   real(wp), parameter :: step = 1.0e-6_wp
   type(structure_type) :: mol_work
   real(wp), allocatable :: trans(:, :)
   real(wp), allocatable :: g_plus(:, :), g_minus(:, :)
   real(wp) :: s_plus(3, 3), s_minus(3, 3)
   real(wp), allocatable :: xyz_orig(:, :)
   real(wp) :: lattice_orig(3, 3), eps_mat(3, 3)
   integer :: ic, jc, kc, lc, jat
   logical :: periodic

   hess(:) = 0.0_wp
   press(:, :, :, :) = 0.0_wp

   ! Mutable local copy of the structure
   periodic = any(mol%periodic)
   mol_work = mol
   allocate(xyz_orig(3, mol%nat))
   xyz_orig(:, :) = mol%xyz
   lattice_orig(:, :) = 0.0_wp
   if (periodic) lattice_orig(:, :) = mol%lattice

   ! Evaluate the unperturbed system
   call get_lattice_points(mol%periodic, mol%lattice, self%ncoord%cutoff, trans)
   call self%update(mol, cache, trans, .true., list)
   call self%solve(mol, solver, cache, error, qvec=qvec, energy=energy, &
      & gradient=grad, sigma=sigma, list=list, unit=unit, verbosity=verbosity)
   if (allocated(error)) return

   allocate(g_plus(3, mol%nat), g_minus(3, mol%nat))

   ! Cartesian displacements
   do jat = 1, mol%nat
      do jc = 1, 3
         mol_work%xyz(jc, jat) = xyz_orig(jc, jat) + step
         call get_displaced_gradient(self, mol_work, solver, error, g_plus, s_plus, &
            & list, unit, verbosity)
         if (allocated(error)) return

         mol_work%xyz(jc, jat) = xyz_orig(jc, jat) - step
         call get_displaced_gradient(self, mol_work, solver, error, g_minus, s_minus, &
            & list, unit, verbosity)
         if (allocated(error)) return

         mol_work%xyz(jc, jat) = xyz_orig(jc, jat)

         call add_hess_column(3 * (jat - 1) + jc, &
            & 0.5_wp * (g_plus - g_minus) / step, hess)
         if (size(press, 4) == mol%nat) then
            press(:, :, jc, jat) = 0.5_wp * (s_plus - s_minus) / step
         end if
      end do
   end do

   ! Strain deformations
   if (size(press, 4) == 3) then
      eps_mat(:, :) = 0.0_wp
      do ic = 1, 3
         eps_mat(ic, ic) = 1.0_wp
      end do

      do kc = 1, 3
         do lc = 1, 3
            eps_mat(lc, kc) = eps_mat(lc, kc) + step
            mol_work%xyz(:, :) = matmul(eps_mat, xyz_orig)
            if (periodic) mol_work%lattice(:, :) = matmul(eps_mat, lattice_orig)
            call get_displaced_gradient(self, mol_work, solver, error, g_plus, s_plus, &
               & list, unit, verbosity)
            if (allocated(error)) return

            eps_mat(lc, kc) = eps_mat(lc, kc) - 2.0_wp * step
            mol_work%xyz(:, :) = matmul(eps_mat, xyz_orig)
            if (periodic) mol_work%lattice(:, :) = matmul(eps_mat, lattice_orig)
            call get_displaced_gradient(self, mol_work, solver, error, g_minus, s_minus, &
               & list, unit, verbosity)
            if (allocated(error)) return

            eps_mat(lc, kc) = eps_mat(lc, kc) + step
            mol_work%xyz(:, :) = xyz_orig
            if (periodic) mol_work%lattice(:, :) = lattice_orig

            press(:, :, lc, kc) = 0.5_wp * (s_plus - s_minus) / step
         end do
      end do
   end if

end subroutine get_numhess_packed


!> Semi-numerical Hessian in dense symmetric storage
subroutine get_numhess_dense(self, mol, solver, cache, error, qvec, energy, grad, &
   & sigma, hess, press, list, unit, verbosity)

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

   real(wp), allocatable :: hess_packed(:)
   integer :: idx, jdx, ij, ndim

   hess(:, :) = 0.0_wp

   ndim = 3 * mol%nat
   allocate(hess_packed(ndim * (ndim + 1) / 2))
   call self%get_numhess_packed(mol, solver, cache, error, qvec, energy, grad, sigma, &
      & hess_packed, press, list, unit, verbosity)
   if (allocated(error)) return

   ij = 0
   do idx = 1, ndim
      do jdx = 1, idx
         ij = ij + 1
         hess(jdx, idx) = hess_packed(ij)
         hess(idx, jdx) = hess_packed(ij)
      end do
   end do

end subroutine get_numhess_dense


!> Accumulate one column of the Hessian in packed storage
subroutine add_hess_column(jdx, dgdr, hess)

   !> Index of the displaced Cartesian coordinate
   integer, intent(in) :: jdx

   !> Derivative of the energy gradient w.r.t. the displaced coordinate
   real(wp), intent(in) :: dgdr(:, :)

   !> Packed lower-triangle Hessian
   real(wp), intent(inout) :: hess(:)

   integer :: iat, ic, idx, ij

   do iat = 1, size(dgdr, 2)
      do ic = 1, 3
         idx = 3 * (iat - 1) + ic
         ij = hess_index(idx, jdx)
         if (idx == jdx) then
            hess(ij) = dgdr(ic, iat)
         else
            hess(ij) = hess(ij) + 0.5_wp * dgdr(ic, iat)
         end if
      end do
   end do

end subroutine add_hess_column


!> Position of a symmetric-matrix element in packed lower-triangle storage
elemental function hess_index(idx, jdx) result(ij)

   !> Row index
   integer, intent(in) :: idx

   !> Column index
   integer, intent(in) :: jdx

   !> Position in the packed lower triangle
   integer :: ij

   ij = max(idx, jdx) * (max(idx, jdx) - 1) / 2 + min(idx, jdx)

end function hess_index


!> Energy gradient and virial for a displaced structure
subroutine get_displaced_gradient(self, mol, solver, error, gradient, sigma, &
   & list, unit, verbosity)

   !> Electronegativity-equilibration model
   class(mchrg_model_type), intent(in) :: self

   !> Displaced molecular structure data
   type(structure_type), intent(in) :: mol

   !> The solver instance
   class(mchrg_solver_type), intent(in) :: solver

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Energy gradient
   real(wp), intent(out), contiguous :: gradient(:, :)

   !> Virial
   real(wp), intent(out), contiguous :: sigma(:, :)

   !> Reference neighborlist
   type(csr_list), intent(in), optional :: list

   !> Output unit
   integer, intent(in), optional :: unit

   !> Verbosity level
   integer, intent(in), optional :: verbosity

   type(mchrg_cache) :: cache
   type(csr_list), allocatable :: list_work
   type(distance_screening) :: screen
   real(wp), allocatable :: trans(:, :)

   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   call get_lattice_points(mol%periodic, mol%lattice, self%ncoord%cutoff, trans)

   if (present(list)) then
      allocate(list_work)
      if (any(mol%periodic)) then
         call new_csr_list(list_work, screen, mol, error, cache%wsc, list%cutoff, &
            & complete=.true.)
         if (allocated(error)) return
         call sort_neighborlist(list_work, cache%wsc)
      else
         call new_csr_list(list_work, screen, mol, error, cutoff=list%cutoff, &
            & complete=.true.)
         if (allocated(error)) return
         call sort_neighborlist(list_work)
      end if
   end if

   call self%update(mol, cache, trans, .true., list_work)
   call self%solve(mol, solver, cache, error, gradient=gradient, sigma=sigma, &
      & list=list_work, unit=unit, verbosity=verbosity)

end subroutine get_displaced_gradient


!> Local charges calculation
subroutine local_charge(self, mol, trans, qloc, dqlocdr, dqlocdL, &
   & list, dqlocdrlist)
   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Lattice translation vectors
   real(wp), intent(in) :: trans(:, :)

   !> Local atomic partial charges
   real(wp), intent(out) :: qloc(:)

   !> Optional derivative of local atomic partial charges w.r.t. atomic positions
   real(wp), intent(out), optional :: dqlocdr(3, mol%nat, mol%nat)

   !> Optional derivative of local atomic partial charges w.r.t. lattice vectors
   real(wp), intent(out), optional :: dqlocdL(3, 3, mol%nat)

   !> Multicharge neighborlist type
   type(csr_list), intent(in), optional :: list

   !> Optional derivative of local charges w.r.t. positions on the list
   real(wp), intent(out), optional :: dqlocdrlist(:, :)

   qloc = 0.0_wp
   if (present(dqlocdr) .and. present(dqlocdL)) then
      dqlocdr = 0.0_wp
      dqlocdL = 0.0_wp
   end if
   if (present(list) .and. present(dqlocdrlist) .and. present(dqlocdL)) then
      dqlocdrlist = 0.0_wp
      dqlocdL = 0.0_wp
   end if
   ! Get the electronegativity weighted CN for local charge
   if (allocated(self%ncoord_en)) then
      call self%ncoord_en%get_coordination_number(mol, trans, qloc, &
         & dcndr=dqlocdr, dcndrlist=dqlocdrlist, dcndL=dqlocdL, list=list)
   end if

   ! Distribute the total charge equally
   qloc = qloc + mol%charge / real(mol%nat, wp)

end subroutine local_charge


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


!> Print header for the analytical Hessian
subroutine print_hessian_header(unit, verbosity)
   !> Output unit
   integer, intent(in) :: unit

   !> Verbosity level
   integer, intent(in) :: verbosity

   if (verbosity > 0) then
      write(unit, '(54("-"))')
      write(unit, '(18x, a)') "Analytical Hessian"
      write(unit, '(54("-"))')
      write(unit, '(a)') ''
   end if
end subroutine print_hessian_header


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


!> Print the elapsed time of the stages of the Hessian on a neighborlist
subroutine print_hessian_stages(unit, timer, tgroup, adjoint)
   !> Output unit
   integer, intent(in) :: unit

   !> Timer of the stages
   type(timer_type), intent(in) :: timer

   !> Times of the group derivatives and their contraction
   real(wp), intent(in) :: tgroup(2)

   !> Print the stage of the adjoint solves of an external function
   logical, intent(in), optional :: adjoint

   if (present(adjoint)) then
      if (adjoint) call print_stage_time(unit, "Adjoint vector", &
         & timer%get("adjoint"))
   end if
   call print_stage_time(unit, "Partial derivatives dB/dR", timer%get("dabdr"))
   call print_stage_time(unit, "Inverse on the neighborlist", &
      & timer%get("inverse"))
   call print_stage_time(unit, "Implicit Hessian", timer%get("implicit"))
   write(unit, '(a, t34, ":", 2(1x, f0.3))') &
      & "  S*dB/dR/contraction (thread s)", tgroup
   call print_stage_time(unit, "Explicit Hessian", timer%get("explicit"))
   call print_stage_time(unit, "Total Hessian", timer%get("total"))
   write(unit, '(a)') ''
end subroutine print_hessian_stages

end module multicharge_model_type
