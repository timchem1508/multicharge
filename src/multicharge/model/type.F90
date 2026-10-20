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
   use mctc_cutoff, only : get_lattice_points
   use mctc_ncoord, only : ncoord_type
   use mctc_csrlist, only : csr_list, spgemv_csr, spsymv_csr, new_csr_list
   use multicharge_blas, only : gemv, symv, gemm
   use multicharge_model_cache, only : mchrg_cache
   use multicharge_solver_type, only : mchrg_solver_type
   use multicharge_solver_cg, only : cg_solver
   use multicharge_lowerlist, only : lower_list, count_neighbors, get_neighbor

   implicit none
   private

   public :: mchrg_model_type, mchrg_hessian, get_dir_trans, get_rec_trans, &
      & hess_index, get_dqdr_row, new_partial_derivs


   !> Molecular nuclear Hessian in compressed-row storage on a half neighborlist
   !> of its own, which keeps the blocks of the atoms y >= x in the row of x, the
   !> block of y and x is the transpose. The charge constraint couples all atoms
   !> by the dense rank-1 term sum(u) dlambda dlambda^T, which is kept in its
   !> factors, the block of the atoms x and y is
   !> hesslist(:, :, (x, y)) + uvecsum * dlambda(:, x) dlambda(:, y)^T.
   type :: mchrg_hessian

      !> Half neighborlist of the atom pairs with a stored block, the first entry
      !> of each row is the atom itself
      type(csr_list) :: list

      !> Blocks without the charge constraint term on the pattern of the list,
      !> the entry of atom y in the row of atom x holds d2E/dR_y,jc dR_x,ic in
      !> (jc, ic)
      real(wp), allocatable :: hesslist(:, :, :)

      !> Position derivatives of the Lagrange multiplier, the factors of the
      !> charge constraint term
      real(wp), allocatable :: dlambda(:, :)

      !> Sum of the constraint response, sum(u) = 1^T A^-1 1, the weight of the
      !> charge constraint term
      real(wp) :: uvecsum = 0.0_wp

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

      !> Calculate the analytical Hessian of an external function of the
      !> charges (molecular)
      procedure :: get_external_hessian

      !> Calculate the analytical Hessian of the electrostatic energy (molecular)
      procedure :: get_hessian_dense

      !> Calculate the analytical Hessian of the electrostatic energy on a half
      !> neighborlist in compressed storage (molecular)
      procedure :: get_hessian_list

      !> Calculate the analytical Hessian of the electrostatic energy
      generic :: get_hessian => get_hessian_dense, get_hessian_list

      !> Calculate the implicit Hessian -(dX/dR - dA/dR*q)^T dq/dR from the
      !> response of the charges (molecular)
      procedure :: get_implicit_hessian_dense

      !> Calculate the implicit Hessian from the inverse and the partial
      !> derivatives on a half neighborlist (molecular)
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

      !> Add the Coulomb matrix derivatives contracted with a vector,
      !> alpha * dA/dR*p
      procedure(get_coulomb_derivs), deferred :: get_coulomb_derivs

      !> Calculate capacitance-corrected electronegativity derivatives
      procedure(get_grad), deferred :: get_grad

      !> Calculate the explicit Hessian alpha * p^T d2A/dR2 q + beta * p^T d2X/dR2
      !> at fixed p and q (molecular)
      procedure(get_explicit_hessian_dense), deferred :: get_explicit_hessian_dense

      !> Calculate the explicit Hessian at fixed p and q on a half neighborlist
      !> in compressed storage (molecular)
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

         !> Flag for the charge derivatives w.r.t. positions and strain, stores
         !> the CN (and local charge) derivatives in the cache
         logical, intent(in), optional :: dqdr

         !> Flag for the Hessian, stores the CN (and local charge) derivatives in
         !> the cache
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

      !> Coulomb matrix derivatives contracted with a vector, alpha * dA/dR*p,
      !> w.r.t. positions and strain, added to the partial derivatives in the
      !> cache, dabdr (dabdrlist_ut and dabdrlist_lt with a neighborlist) and
      !> dabdL, allocated with new_partial_derivs
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

         !> Multicharge neighborlist type (half)
         type(csr_list), intent(in), optional :: list

         !> Add to dafdr (dafdrlist_ut and dafdrlist_lt with a neighborlist)
         !> for an external vector instead, without strain derivatives
         logical, intent(in), optional :: dafdr
      end subroutine get_coulomb_derivs

      !> Electronegativity vector derivatives, beta * dX/dR, w.r.t. positions and
      !> strain, added to the partial derivatives in the cache, dabdr
      !> (dabdrlist_ut and dabdrlist_lt with a neighborlist) and dabdL, allocated
      !> with new_partial_derivs
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

         !> Multicharge neighborlist type (half)
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

         !> neighborlist (each unordered pair appears once)
         type(csr_list), optional, intent(in) :: list
      end subroutine get_grad

      !> Explicit part of the molecular nuclear Hessian at fixed p and charges q,
      !> alpha * p^T d2A/dR2 q + beta * p^T d2X/dR2, added to the Hessian
      subroutine get_explicit_hessian_dense(self, mol, cache, p, hessian, error, &
         & alpha, beta)
         import :: mchrg_model_type, structure_type, mchrg_cache, error_type, wp

         !> Multicharge model type
         class(mchrg_model_type), intent(in) :: self

         !> Structure type
         type(structure_type), intent(in) :: mol

         !> Multicharge cache with the position derivatives of the CN and the
         !> local charges from the update for the Hessian
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

      !> Explicit part of the molecular nuclear Hessian at fixed p and charges q
      !> on a half neighborlist, added to the blocks of the compressed Hessian on
      !> the pattern of its own half list
      subroutine get_explicit_hessian_list(self, mol, cache, p, hessian, error, &
         & list, alpha, beta)
         import :: mchrg_model_type, structure_type, mchrg_cache, csr_list, &
            & mchrg_hessian, error_type, wp

         !> Multicharge model type
         class(mchrg_model_type), intent(in) :: self

         !> Structure type
         type(structure_type), intent(in) :: mol

         !> Multicharge cache with the position derivatives of the CN and the
         !> local charges on the neighborlist from the update for the Hessian
         type(mchrg_cache), intent(inout) :: cache

         !> Charge-like contraction vector
         real(wp), intent(in) :: p(:)

         !> Molecular nuclear Hessian in compressed storage with its list and
         !> allocated blocks
         type(mchrg_hessian), intent(inout) :: hessian

         !> Error handling
         type(error_type), allocatable, intent(out) :: error

         !> Multicharge neighborlist type (half)
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

   !> Largest number of atoms in a group of the charge derivatives on a
   !> neighborlist
   integer, parameter :: max_group = 8

   !> Largest distance of an atom of a group to its first atom, relative to the
   !> cutoff of the list
   real(wp), parameter :: group_scale = 0.2_wp


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

   !> Optional derivative of the atomic partial charges w.r.t. atomic positions,
   !> with a neighborlist they are kept in compressed-row storage in the cache
   !> and only expanded into dqdr if present
   real(wp), intent(out), contiguous, optional :: dqdr(:, :, :)

   !> Optional derivative of the atomic partial charges w.r.t. lattice vectors,
   !> requests the charge derivatives together with dqdr or a neighborlist
   real(wp), intent(out), contiguous, optional :: dqdL(:, :, :)

   !> neighborlist optional type
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

   logical :: grad, cpq
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

   ! Charge derivatives are evaluated on a half neighborlist with its lower
   ! triangle
   if (cpq .and. present(list)) then
      if (list%complete) then
         call fatal_error(error, "Charge derivatives on a neighborlist require "// &
            & "a half neighborlist")
         return
      end if
      if (.not. allocated(cache%lowlist%inl)) then
         call fatal_error(error, "Charge derivatives on a half neighborlist "// &
            & "require its lower triangle, update the cache with dqdr")
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
         if (list%complete) then
            call spgemv_csr(mol%nat, cache%alist, list%inl, list%nlat, &
               & cache%vrhs, avec)
         else
            call spsymv_csr(mol%nat, cache%alist, list%inl, list%nlat, &
               & cache%vrhs, avec)
         end if
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

   !> Derivative of the atomic partial charges w.r.t. atomic positions, with a
   !> neighborlist only expanded from the compressed storage if present
   real(wp), intent(out), optional :: dqdr(:, :, :)

   !> Derivative of the atomic partial charges w.r.t. lattice vectors
   real(wp), intent(out) :: dqdL(:, :, :)

   !> neighborlist optional type
   type(csr_list), intent(in), optional :: list

   !> Output unit
   integer, intent(in), optional :: unit

   !> Verbosity level, the times of the inverse and the product on a
   !> neighborlist are printed from 1
   integer, intent(in), optional :: verbosity

   real(wp), allocatable :: daqxdr(:, :, :), daqxdL(:, :, :), ainv(:, :)
   real(wp), allocatable :: dlambdadr(:, :)
   real(wp) :: scale, uvecsum
   integer :: jat, print_unit
   type(timer_type) :: timer

   ! Without a neighborlist the position derivatives are only formed densely
   if (.not. (present(list) .or. present(dqdr))) then
      call fatal_error(error, "Charge derivatives without a neighborlist require dqdr")
      return
   end if

   select type (solver)
   class is (cg_solver)
      if (present(list)) then
         ! J^-1 decays fast with the distance and is kept on the pattern of the
         ! list, the dense part of (P J^-1)^T = J^-1 - u*u^T/sum(u) is the
         ! rank-1 term of the charge constraint, which is rebuilt from uvec
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

         if (present(dqdr)) then
            !$omp parallel do default(none) schedule(runtime) &
            !$omp shared(mol, cache, list, dqdr) private(jat)
            do jat = 1, mol%nat
               call get_dqdr_row(cache, list, jat, dqdr(:, :, jat))
            end do
         end if
         return
      end if

      ! Inverse of J by block CG instead of 3N+9 solves with the columns of
      ! dB/dR - dA/dR*q, then dq/dR = (dB/dR - dA/dR*q)^T (P J^-1)^T
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


!> Inverse S = J^-1 of the coefficient matrix on the pattern of a half
!> neighborlist in cache%ainvlist with the lower triangle of the list from the
!> cache
subroutine get_inverse_list(mol, solver, cache, list, error, unit)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Conjugate-gradient solver with the tolerance of the inverse
   class(cg_solver), intent(in) :: solver

   !> Cache with the coefficient matrix on the list, receives the inverse
   type(mchrg_cache), intent(inout) :: cache

   !> Half neighborlist
   type(csr_list), intent(in) :: list

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Output unit
   integer, intent(in), optional :: unit

   if (allocated(cache%ainvlist)) then
      if (size(cache%ainvlist, kind=i8) /= size(list%nlat, kind=i8)) &
         & deallocate(cache%ainvlist)
   end if
   if (.not. allocated(cache%ainvlist)) &
      & allocate(cache%ainvlist(size(list%nlat, kind=i8)))
   call solver%invert_list(cache%alist, list, cache%ainvlist, xyz=mol%xyz, &
      & low=cache%lowlist, new_unit=unit, error=error)

end subroutine get_inverse_list


!> Position derivatives of the charges from the inverse S = J^-1 on the pattern
!> of a half neighborlist, dq_j/dR_i = (S B)_ji - u_j w_i / sum(u) with
!> B = dB/dR - dA/dR*q and w = B^T u. The sparse product S B is kept in
!> compressed storage on the pattern of the list, dq_j/dR_i for i >= j in
!> cache%dqdrlist_ut at the entry of i in the row of j and for i < j in
!> cache%dqdrlist_lt at the entry of j in the row of i. The columns of S B are
!> formed for groups of consecutive atoms, the elements outside of the pattern
!> are of the order of the inverse at the cutoff and are dropped. The factors
!> w / sum(u) of the charge constraint term follow from get_constraint_factors.
subroutine get_dqdr_list(mol, cache, list)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Cache with the inverse and the partial derivatives on the list, receives
   !> the compressed position derivatives
   type(mchrg_cache), intent(inout) :: cache

   !> Half neighborlist
   type(csr_list), intent(in) :: list

   integer :: igrp, ngrp, nbmax, x0, ng, nb, n, i
   integer(i8) :: kat
   integer, allocatable :: gstart(:), bloc(:), blist(:)
   integer(i8), allocatable :: qpos(:, :)
   real(wp), allocatable :: pan(:, :), qcol(:, :)

   call get_list_groups(mol, list, gstart)
   ngrp = size(gstart) - 1
   nbmax = min(mol%nat, max_group * get_max_neighbors(mol%nat, list, cache%lowlist))

   if (allocated(cache%dqdrlist_ut)) deallocate(cache%dqdrlist_ut)
   if (allocated(cache%dqdrlist_lt)) deallocate(cache%dqdrlist_lt)
   allocate(cache%dqdrlist_ut(3, size(list%nlat, kind=i8)), &
      & cache%dqdrlist_lt(3, size(list%nlat, kind=i8)), source=0.0_wp)

   ! Only the group of atom x writes the derivatives w.r.t. R_x
   !$omp parallel default(none) &
   !$omp shared(mol, cache, list, gstart, ngrp, nbmax) &
   !$omp private(igrp, x0, ng, nb, n, i, kat, bloc, blist, qpos, pan, qcol)
   allocate(bloc(mol%nat), source=0)
   allocate(blist(nbmax), qpos(max_group, nbmax), pan(3 * max_group, nbmax), &
      & qcol(3 * max_group, nbmax))
   !$omp do schedule(runtime)
   do igrp = 1, ngrp
      x0 = gstart(igrp)
      ng = gstart(igrp + 1) - x0
      call get_group_dqdr(cache, list, x0, ng, nb, bloc, blist, qpos, pan, qcol)
      do n = 1, nb
         do i = 1, ng
            kat = qpos(i, n)
            if (kat > 0) then
               cache%dqdrlist_ut(:, kat) = qcol(3 * i - 2:3 * i, n)
            else if (kat < 0) then
               cache%dqdrlist_lt(:, -kat) = qcol(3 * i - 2:3 * i, n)
            end if
         end do
         bloc(blist(n)) = 0
      end do
   end do
   !$omp end do
   deallocate(bloc, blist, qpos, pan, qcol)
   !$omp end parallel

end subroutine get_dqdr_list


!> Factors w / sum(u) of the charge constraint term of the position derivatives
!> of the charges, w = B^T u with the partial derivatives B = dB/dR - dA/dR*q on
!> a half neighborlist and the constraint response u
subroutine get_constraint_factors(mol, cache, list, dlambdadr)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Cache with the partial derivatives on the list and the constraint
   !> response
   type(mchrg_cache), intent(in) :: cache

   !> Half neighborlist
   type(csr_list), intent(in) :: list

   !> Factors of the charge constraint term, the position derivatives of the
   !> Lagrange multiplier
   real(wp), allocatable, intent(out) :: dlambdadr(:, :)

   integer :: iat, jat
   integer(i8) :: kat
   real(wp) :: uvecsum

   uvecsum = sum(cache%uvec)
   allocate(dlambdadr(3, mol%nat), source=0.0_wp)
   do iat = 1, mol%nat
      do kat = list%inl(iat), list%inl(iat + 1) - 1
         jat = list%nlat(kat)
         dlambdadr(:, jat) = dlambdadr(:, jat) &
            & + cache%dabdrlist_ut(:, kat) * cache%uvec(iat)
         if (jat /= iat) dlambdadr(:, iat) = dlambdadr(:, iat) &
            & + cache%dabdrlist_lt(:, kat) * cache%uvec(jat)
      end do
   end do
   dlambdadr(:, :) = dlambdadr / (uvecsum + eps)

end subroutine get_constraint_factors


!> Derivatives of the charges w.r.t. the lattice vectors, B_L^T S with the
!> inverse S on a half neighborlist and the strain derivatives B_L of the
!> partial derivatives, and their charge constraint term
subroutine get_dqdL_list(mol, cache, list, dqdL)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Cache with the inverse on the list, the strain derivatives of the partial
   !> derivatives and the lower triangle of the list
   type(mchrg_cache), intent(in) :: cache

   !> Half neighborlist
   type(csr_list), intent(in) :: list

   !> Derivative of the atomic partial charges w.r.t. lattice vectors
   real(wp), intent(out) :: dqdL(:, :, :)

   integer :: iat, jat, ient
   integer(i8) :: kat
   real(wp) :: uvecsum, scale
   real(wp) :: wlat(3, 3)

   uvecsum = sum(cache%uvec)
   wlat(:, :) = 0.0_wp
   do iat = 1, mol%nat
      wlat(:, :) = wlat + cache%dabdL(:, :, iat) * cache%uvec(iat)
   end do

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(mol, cache, list, uvecsum, wlat, dqdL) &
   !$omp private(iat, jat, ient, kat, scale)
   do jat = 1, mol%nat
      scale = cache%uvec(jat) / (uvecsum + eps)
      dqdL(:, :, jat) = -scale * wlat
      do ient = 1, count_neighbors(list, cache%lowlist, jat)
         call get_neighbor(list, cache%lowlist, jat, ient, iat, kat)
         dqdL(:, :, jat) = dqdL(:, :, jat) + &
            & cache%ainvlist(kat) * cache%dabdL(:, :, iat)
      end do
   end do

end subroutine get_dqdL_list


!> Charge derivatives dq_b/dR_x of a group of consecutive atoms x for the
!> union of their neighbors b on a half neighborlist, in
!> qcol(3*(i-1)+1:3*i, n) for atom i of the group and atom n of the union.
!> The product of the inverse S on the list with the columns x of the partial
!> derivatives is only formed on the rows of the union, the elements of the
!> atoms b outside of the list of x are dropped. The caller resets bloc for
!> the atoms of the union.
subroutine get_group_dqdr(cache, list, x0, ng, nb, bloc, blist, qpos, pan, qcol)

   !> Cache with the inverse and the partial derivatives on the list and the
   !> lower triangle of the list
   type(mchrg_cache), intent(in) :: cache

   !> Half neighborlist
   type(csr_list), intent(in) :: list

   !> First atom of the group
   integer, intent(in) :: x0

   !> Number of atoms of the group
   integer, intent(in) :: ng

   !> Number of atoms in the union of the neighbors of the group
   integer, intent(out) :: nb

   !> Index of each atom in the union, zero outside
   integer, intent(inout), contiguous :: bloc(:)

   !> Atoms of the union
   integer, intent(inout), contiguous :: blist(:)

   !> Position of dq_b/dR_x in the compressed storage on the list, zero if b is
   !> not a neighbor of x, negative for the lower triangle
   integer(i8), intent(inout), contiguous :: qpos(:, :)

   !> Partial derivatives dB_c/dR_x of the union in the layout of qcol
   real(wp), intent(inout), contiguous :: pan(:, :)

   !> Charge derivatives dq_b/dR_x of the union
   real(wp), intent(inout), contiguous :: qcol(:, :)

   integer :: n, i

   call gather_group_list(list, cache%lowlist, cache%dabdrlist_ut, &
      & cache%dabdrlist_lt, x0, ng, nb, bloc, blist, qpos, pan)
   call multiply_group_list(list, cache%ainvlist, 3 * ng, nb, bloc, blist, pan, &
      & qcol)
   do n = 1, nb
      do i = 1, ng
         if (qpos(i, n) == 0) qcol(3 * i - 2:3 * i, n) = 0.0_wp
      end do
   end do

end subroutine get_group_dqdr


!> Union of the neighbors of a group of consecutive atoms on a half list with
!> the columns x of the partial derivatives of the group atoms. The row of x in
!> the list holds the diagonal and B(c, x) of c > x in the lower triangle, the
!> atoms c < x come from the lower triangle of the list and B(c, x) from the
!> upper triangle in the row of c.
subroutine gather_group_list(list, low, but, blt, x0, ng, nb, bloc, blist, qpos, &
   & pan)

   !> Half neighborlist
   type(csr_list), intent(in) :: list

   !> Lower triangle of the list
   type(lower_list), intent(in) :: low

   !> Diagonal and upper triangle of the partial derivatives on the list
   real(wp), intent(in), contiguous :: but(:, :)

   !> Lower triangle of the partial derivatives on the list
   real(wp), intent(in), contiguous :: blt(:, :)

   !> First atom of the group
   integer, intent(in) :: x0

   !> Number of atoms of the group
   integer, intent(in) :: ng

   !> Number of atoms in the union of the neighbors of the group
   integer, intent(out) :: nb

   !> Index of each atom in the union, zero outside
   integer, intent(inout), contiguous :: bloc(:)

   !> Atoms of the union
   integer, intent(inout), contiguous :: blist(:)

   !> Position of the pair of x and each atom of the union in the list, zero if
   !> the atom is not a neighbor of x, negative for the lower triangle
   integer(i8), intent(inout), contiguous :: qpos(:, :)

   !> Partial derivatives dB_c/dR_x of the union
   real(wp), intent(inout), contiguous :: pan(:, :)

   integer :: i, ic, xat, n
   integer(i8) :: kat, lat

   nb = 0
   do i = 1, ng
      xat = x0 + i - 1
      ic = 3 * (i - 1)
      kat = list%inl(xat)
      call add_union_atom(xat, nb, bloc, blist, qpos, pan, n)
      pan(ic + 1:ic + 3, n) = but(:, kat)
      qpos(i, n) = kat
      do kat = list%inl(xat) + 1, list%inl(xat + 1) - 1
         call add_union_atom(list%nlat(kat), nb, bloc, blist, qpos, pan, n)
         pan(ic + 1:ic + 3, n) = blt(:, kat)
         qpos(i, n) = -kat
      end do
      do lat = low%inl(xat), low%inl(xat + 1) - 1
         kat = low%pos(lat)
         call add_union_atom(low%nlat(lat), nb, bloc, blist, qpos, pan, n)
         pan(ic + 1:ic + 3, n) = but(:, kat)
         qpos(i, n) = kat
      end do
   end do

end subroutine gather_group_list


!> Index of an atom in the union of the neighbors of a group, a new atom is
!> appended without derivatives
pure subroutine add_union_atom(bat, nb, bloc, blist, qpos, pan, n)

   !> Atom to add
   integer, intent(in) :: bat

   !> Number of atoms in the union
   integer, intent(inout) :: nb

   !> Index of each atom in the union, zero outside
   integer, intent(inout), contiguous :: bloc(:)

   !> Atoms of the union
   integer, intent(inout), contiguous :: blist(:)

   !> Positions of the derivatives of the union
   integer(i8), intent(inout), contiguous :: qpos(:, :)

   !> Partial derivatives of the union
   real(wp), intent(inout), contiguous :: pan(:, :)

   !> Index of the atom in the union
   integer, intent(out) :: n

   n = bloc(bat)
   if (n > 0) return
   nb = nb + 1
   n = nb
   bloc(bat) = n
   blist(n) = bat
   qpos(:, n) = 0_i8
   pan(:, n) = 0.0_wp

end subroutine add_union_atom


!> Product of the inverse S on a half list with the partial derivatives of the
!> union U of a group, qcol = S_UU pan. The rows of the union are read from the
!> list, each pair of the upper triangle also adds the mirrored element of the
!> symmetric inverse to the row of its other atom, so the lower triangle of the
!> list is not needed.
subroutine multiply_group_list(list, ainv, nw, nb, bloc, blist, pan, qcol)

   !> Half neighborlist
   type(csr_list), intent(in) :: list

   !> Upper triangle of the inverse on the pattern of the list
   real(wp), intent(in), contiguous :: ainv(:)

   !> Number of derivatives of each atom of the union, three per group atom
   integer, intent(in) :: nw

   !> Number of atoms in the union
   integer, intent(in) :: nb

   !> Index of each atom in the union, zero outside
   integer, intent(in), contiguous :: bloc(:)

   !> Atoms of the union
   integer, intent(in), contiguous :: blist(:)

   !> Partial derivatives of the union
   real(wp), intent(in), contiguous :: pan(:, :)

   !> Product with the inverse on the union
   real(wp), intent(inout), contiguous :: qcol(:, :)

   integer :: n, m, bat
   integer(i8) :: kat
   real(wp) :: s, acc(3 * max_group)

   qcol(:nw, :nb) = 0.0_wp
   do n = 1, nb
      bat = blist(n)
      acc(:nw) = 0.0_wp
      do kat = list%inl(bat), list%inl(bat + 1) - 1
         m = bloc(list%nlat(kat))
         if (m == 0) cycle
         s = ainv(kat)
         acc(:nw) = acc(:nw) + s * pan(:nw, m)
         if (m /= n) qcol(:nw, m) = qcol(:nw, m) + s * pan(:nw, n)
      end do
      qcol(:nw, n) = qcol(:nw, n) + acc(:nw)
   end do

end subroutine multiply_group_list


!> Groups of consecutive atoms close to each other for the charge derivatives
!> on a neighborlist, which share most of their neighbors
subroutine get_list_groups(mol, list, gstart)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Half neighborlist
   type(csr_list), intent(in) :: list

   !> First atom of each group, the last entry is the number of atoms plus one
   integer, allocatable, intent(out) :: gstart(:)

   real(wp) :: radius

   radius = 0.0_wp
   if (allocated(list%cutoff)) radius = group_scale * list%cutoff
   call get_atom_groups(mol, radius, max_group, gstart)

end subroutine get_list_groups


!> Largest number of neighbors of an atom in a half list with its lower
!> triangle, including the atom itself
function get_max_neighbors(nat, list, low) result(kmax)

   !> Number of atoms
   integer, intent(in) :: nat

   !> Half neighborlist
   type(csr_list), intent(in) :: list

   !> Lower triangle of the list
   type(lower_list), intent(in) :: low

   !> Largest number of neighbors
   integer :: kmax

   integer :: iat

   kmax = 0
   do iat = 1, nat
      kmax = max(kmax, count_neighbors(list, low, iat))
   end do

end function get_max_neighbors


!> Allocate the partial derivatives in the cache and set them to zero, which
!> get_xvec_derivs and get_coulomb_derivs add to: dabdr (dabdrlist_ut and
!> dabdrlist_lt with a neighborlist) and dabdL, or with dafdr the Coulomb
!> matrix derivatives contracted with an external vector, dafdr (dafdrlist_ut
!> and dafdrlist_lt with a neighborlist)
subroutine new_partial_derivs(mol, ndim, cache, list, dafdr)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> System size
   integer, intent(in) :: ndim

   !> Multicharge cache
   type(mchrg_cache), intent(inout) :: cache

   !> Multicharge neighborlist type (half)
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
         if (allocated(cache%dafdrlist_ut)) deallocate(cache%dafdrlist_ut)
         if (allocated(cache%dafdrlist_lt)) deallocate(cache%dafdrlist_lt)
         allocate(cache%dafdrlist_ut(3, nnz), cache%dafdrlist_lt(3, nnz), &
            & source=0.0_wp)
      else
         if (allocated(cache%dabdrlist_ut)) deallocate(cache%dabdrlist_ut)
         if (allocated(cache%dabdrlist_lt)) deallocate(cache%dabdrlist_lt)
         allocate(cache%dabdrlist_ut(3, nnz), cache%dabdrlist_lt(3, nnz), &
            & source=0.0_wp)
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


!> Position derivatives of the charge of one atom w.r.t. all atomic positions,
!> expanded from the compressed storage on the pattern of the half
!> neighborlist of the solve and the charge constraint term
subroutine get_dqdr_row(cache, list, jat, dqdr)

   !> Cache with the compressed position derivatives and the lower triangle of
   !> the list
   type(mchrg_cache), intent(in) :: cache

   !> Half neighborlist of the charge derivatives
   type(csr_list), intent(in) :: list

   !> Atom of the charge
   integer, intent(in) :: jat

   !> Derivatives dq_j/dR_i for all atoms i, dimension (3, nat)
   real(wp), intent(out) :: dqdr(:, :)

   integer :: iat
   integer(i8) :: kat

   do iat = 1, size(dqdr, 2)
      dqdr(:, iat) = -cache%dlambdadr(:, iat) * cache%uvec(jat)
   end do

   ! The atoms of lower index are in the lower triangle of the list
   do kat = list%inl(jat), list%inl(jat + 1) - 1
      iat = list%nlat(kat)
      dqdr(:, iat) = dqdr(:, iat) + cache%dqdrlist_ut(:, kat)
   end do
   do kat = cache%lowlist%inl(jat), cache%lowlist%inl(jat + 1) - 1
      iat = cache%lowlist%nlat(kat)
      dqdr(:, iat) = dqdr(:, iat) + cache%dqdrlist_lt(:, cache%lowlist%pos(kat))
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

   !> neighborlist optional type
   type(csr_list), intent(in), optional :: list

   !> Output unit
   integer, intent(in), optional :: unit

   !> Verbosity level
   integer, intent(in), optional :: verbosity

   integer :: iat
   integer :: ndim
   real(wp), allocatable :: yvec(:)
   real(wp), allocatable :: padj(:)
   real(wp), allocatable :: dfdq_loc(:)
   real(wp) :: uvecsum
   real(wp) :: yvecsum
   real(wp) :: scale
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

   if (solver%need_pos_def) then
      ndim = mol%nat
   else
      ndim = mol%nat + 1
      allocate(dfdq_loc(ndim), source=0.0_wp)
      dfdq_loc(:mol%nat) = dfdq
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

   if (solver%need_pos_def) then
      allocate(yvec(ndim))
      do iat = 1, mol%nat
         yvec(iat) = dfdq(iat) / cache%amat(iat, iat)
      end do

      ! Unconstrained response: J*yvec = dfdq
      call print_adjoint_message(print_unit, verbosity_solve)
      call solver%solve(amat=cache%amat, alist=cache%alist, &
         & xvec=dfdq, vrhs=yvec, list=list, error=error)
      if (allocated(error)) return

      ! Projection of uvec on yvec
      yvecsum = sum(yvec)
      uvecsum = sum(cache%uvec)
      scale = yvecsum / (uvecsum + eps)
      allocate(padj(mol%nat))
      padj = yvec - scale * cache%uvec
   else
      allocate(padj(ndim))

      ! Direct solution: J*yvec = dfdq
      call print_adjoint_message(print_unit, verbosity_solve)
      call solver%solve(amat=cache%amat, xvec=dfdq_loc, vrhs=padj, &
         & new_unit=print_unit, error=error)
      if (allocated(error)) return
   end if

   ! Evaluate external gradients via adjoint contraction:
   ! dfdr = p^T * (db/dr - dA/dr X q)

   call self%get_grad(mol, cache, padj(:mol%nat), dfdr, dfdL, alpha=-1.0_wp, &
      & beta=1.0_wp, list=list)

   ! pop dfdr timer
   call timer%pop
   call print_gradient_time(print_unit, verbosity_solve, &
      & timer%get("external_gradient"))

end subroutine get_external_gradient


!> Compute the molecular nuclear Hessian of the electrostatic energy, the
!> explicit part at fixed charges and the implicit part from the response of
!> the charges. The charges, energy, gradient and virial of the same solve are
!> optional.
subroutine get_hessian_dense(self, mol, solver, cache, error, hessian, qvec, &
   & energy, gradient, sigma, unit, verbosity)

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

   !> Optional atomic partial charges
   real(wp), intent(out), contiguous, optional :: qvec(:)

   !> Optional electrostatic energy
   real(wp), intent(inout), contiguous, optional :: energy(:)

   !> Optional energy gradient
   real(wp), intent(inout), contiguous, optional :: gradient(:, :)

   !> Optional virial
   real(wp), intent(inout), contiguous, optional :: sigma(:, :)

   !> Output unit
   integer, intent(in), optional :: unit

   !> Verbosity level
   integer, intent(in), optional :: verbosity

   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), allocatable :: dqdr(:, :, :), dqdL(:, :, :)

   if (any(mol%periodic)) then
      call fatal_error(error, &
         & "Analytical Hessian is not available for periodic systems")
      return
   end if

   ! Charges and their derivatives
   call self%update(mol, cache, trans, grad=.false., hess=.true.)
   allocate(dqdL(3, 3, mol%nat), dqdr(3, mol%nat, mol%nat))
   call self%solve(mol, solver, cache, error, energy=energy, gradient=gradient, &
      & sigma=sigma, qvec=qvec, dqdr=dqdr, dqdL=dqdL, verbosity=verbosity, &
      & unit=unit)
   if (allocated(error)) return

   ! Second derivatives of 0.5 * q^T A q - q^T X at fixed charges
   hessian(:, :) = 0.0_wp
   call self%get_explicit_hessian(mol, cache, cache%vrhs(:mol%nat), hessian, &
      & error, alpha=0.5_wp, beta=-1.0_wp)
   if (allocated(error)) return
   call self%get_implicit_hessian(mol, cache, hessian, dqdr=dqdr)

end subroutine get_hessian_dense


!> Compute the molecular nuclear Hessian of the electrostatic energy on a half
!> neighborlist, which requires the CG solver with neighborlist support. The
!> Hessian is kept in compressed storage on the pattern of its own half list,
!> which keeps the blocks of the atoms y >= x, and the dense charge constraint
!> term in its factors. Without a neighborlist the list of the Hessian serves as
!> neighborlist, with a neighborlist the list of the Hessian is a copy of it
!> unless it is set up on entry. The charges, energy, gradient and virial of
!> the same solve are optional.
subroutine get_hessian_list(self, mol, solver, cache, error, hessian, list, qvec, &
   & energy, gradient, sigma, unit, verbosity)

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

   !> Molecular nuclear Hessian in compressed storage, the pattern of an
   !> allocated list is kept
   type(mchrg_hessian), intent(inout) :: hessian

   !> Multicharge neighborlist type (half), the list of the Hessian is used if
   !> absent
   type(csr_list), intent(in), optional :: list

   !> Optional atomic partial charges
   real(wp), intent(out), contiguous, optional :: qvec(:)

   !> Optional electrostatic energy
   real(wp), intent(inout), contiguous, optional :: energy(:)

   !> Optional energy gradient
   real(wp), intent(inout), contiguous, optional :: gradient(:, :)

   !> Optional virial
   real(wp), intent(inout), contiguous, optional :: sigma(:, :)

   !> Output unit
   integer, intent(in), optional :: unit

   !> Verbosity level
   integer, intent(in), optional :: verbosity

   if (any(mol%periodic)) then
      call fatal_error(error, &
         & "Analytical Hessian is not available for periodic systems")
      return
   end if

   if (present(list)) then
      if (.not. allocated(hessian%list%inl)) hessian%list = list
      call get_hessian_csr(self, mol, solver, cache, error, hessian, list, &
         & qvec, energy, gradient, sigma, unit, verbosity)
   else
      if (.not. allocated(hessian%list%inl)) then
         call fatal_error(error, "The Hessian on a neighborlist requires a "// &
            & "neighborlist or the list of the Hessian")
         return
      end if
      ! The list of the Hessian is only read, its blocks are written
      call get_hessian_csr(self, mol, solver, cache, error, hessian, &
         & hessian%list, qvec, energy, gradient, sigma, unit, verbosity)
   end if

end subroutine get_hessian_list


!> Compressed molecular nuclear Hessian on a half neighborlist with the
!> allocated half list of the Hessian. The implicit part is formed first from
!> the inverse and the partial derivatives, which are freed before the
!> explicit part.
subroutine get_hessian_csr(self, mol, solver, cache, error, hessian, list, qvec, &
   & energy, gradient, sigma, unit, verbosity)

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

   !> Molecular nuclear Hessian in compressed storage with its list
   type(mchrg_hessian), intent(inout) :: hessian

   !> Multicharge neighborlist type (half)
   type(csr_list), intent(in) :: list

   !> Optional atomic partial charges
   real(wp), intent(out), contiguous, optional :: qvec(:)

   !> Optional electrostatic energy
   real(wp), intent(inout), contiguous, optional :: energy(:)

   !> Optional energy gradient
   real(wp), intent(inout), contiguous, optional :: gradient(:, :)

   !> Optional virial
   real(wp), intent(inout), contiguous, optional :: sigma(:, :)

   !> Output unit
   integer, intent(in), optional :: unit

   !> Verbosity level
   integer, intent(in), optional :: verbosity

   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp) :: tgroup(2)
   type(timer_type) :: timer
   integer :: print_unit, verbosity_hess

   print_unit = output_unit
   if (present(unit)) print_unit = unit
   verbosity_hess = 0
   if (present(verbosity)) verbosity_hess = verbosity

   if (list%complete .or. hessian%list%complete) then
      call fatal_error(error, "The Hessian on a neighborlist requires a half "// &
         & "neighborlist and a half list of the Hessian")
      return
   end if
   select type (solver)
   class is (cg_solver)
   class default
      call fatal_error(error, "The Hessian on a neighborlist requires the CG solver")
      return
   end select

   call timer%push("total")

   ! CN and local charge derivatives, and the lower triangle of the list
   call timer%push("update")
   call self%update(mol, cache, trans, grad=.false., list=list, hess=.true.)
   call timer%pop

   ! Charges, energy and gradient
   call timer%push("charges")
   call self%solve(mol, solver, cache, error, energy=energy, gradient=gradient, &
      & sigma=sigma, qvec=qvec, list=list, verbosity=verbosity, unit=unit)
   call timer%pop
   if (allocated(error)) return

   ! Partial derivatives dX/dR - dA/dR*q and the inverse on the list
   call print_hessian_header(print_unit, verbosity_hess)
   call timer%push("dabdr")
   call new_partial_derivs(mol, mol%nat, cache, list)
   call self%get_xvec_derivs(mol, mol%nat, cache, beta=1.0_wp, list=list)
   call self%get_coulomb_derivs(mol, mol%nat, cache, alpha=-1.0_wp, &
      & p=cache%vrhs(:mol%nat), list=list)
   call timer%pop
   call timer%push("inverse")
   select type (solver)
   class is (cg_solver)
      call get_inverse_list(mol, solver, cache, list, error, unit)
   end select
   call timer%pop
   if (allocated(error)) return
   deallocate(cache%alist)

   if (allocated(hessian%hesslist)) deallocate(hessian%hesslist)
   allocate(hessian%hesslist(3, 3, size(hessian%list%nlat, kind=i8)), &
      & source=0.0_wp)

   ! Response of the charges, -(dX/dR - dA/dR*q)^T dq/dR
   call timer%push("implicit")
   call self%get_implicit_hessian(mol, cache, hessian, list, times=tgroup)
   call timer%pop
   deallocate(cache%ainvlist, cache%dabdrlist_ut, cache%dabdrlist_lt)

   ! Second derivatives of 0.5 * q^T A q - q^T X at fixed charges
   call timer%push("explicit")
   call self%get_explicit_hessian(mol, cache, cache%vrhs(:mol%nat), hessian, &
      & error, list, alpha=0.5_wp, beta=-1.0_wp)
   call timer%pop
   call timer%pop
   if (allocated(error)) return

   if (verbosity_hess > 0) then
      call print_stage_time(print_unit, "CN and local charge derivatives", &
         & timer%get("update"))
      call print_stage_time(print_unit, "Charges and energy", timer%get("charges"))
      call print_stage_time(print_unit, "Partial derivatives dB/dR", &
         & timer%get("dabdr"))
      call print_stage_time(print_unit, "Inverse on the neighborlist", &
         & timer%get("inverse"))
      call print_stage_time(print_unit, "Implicit Hessian", timer%get("implicit"))
      write(print_unit, '(a, t34, ":", 2(1x, f0.3))') &
         & "  S*dB/dR/contraction (thread s)", tgroup
      call print_stage_time(print_unit, "Explicit Hessian", timer%get("explicit"))
      call print_stage_time(print_unit, "Total Hessian", timer%get("total"))
      write(print_unit, '(a)') ''
   end if

end subroutine get_hessian_csr


!> Compute the implicit part of the molecular nuclear Hessian of the
!> electrostatic energy from the response of the charges,
!> -(dX/dR - dA/dR*q)^T dq/dR. The cache has to hold the partial derivatives
!> dX/dR - dA/dR*q of the dense charge derivatives.
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


!> Compute the implicit part of the molecular nuclear Hessian of the
!> electrostatic energy on a half neighborlist from the inverse S = J^-1 on
!> the pattern of the list and the partial derivatives B = dX/dR - dA/dR*q in
!> the cache. The block of the atoms x <= y is
!> -sum_b (dX_b/dR_y - (dA/dR_y*q)_b) (dq_b/dR_x)^T over the atoms b in the
!> lists of x and y, without the constraint term of the charge derivatives.
!> The charge derivatives dq_b/dR_x = (S B)_bx on the pattern of the list are
!> formed for groups of consecutive atoms x and used at once, they are not
!> stored. The block of y and x is the transpose and the diagonal blocks are
!> symmetrized, so this part is symmetric. The blocks are added on the pattern
!> of the half list of the Hessian, which only keeps the blocks of y >= x. The
!> charge constraint gives the dense rank-1 term
!> sum(u) * vec(dlambda) vec(dlambda)^T, which is kept in its factors.
subroutine get_implicit_hessian_list(self, mol, cache, hessian, list, times)

   !> Electronegativity-equilibration model
   class(mchrg_model_type), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Cache handling
   type(mchrg_cache), intent(in) :: cache

   !> Molecular nuclear Hessian in compressed storage with its list and
   !> allocated blocks
   type(mchrg_hessian), intent(inout) :: hessian

   !> Multicharge neighborlist type (half)
   type(csr_list), intent(in) :: list

   !> Time of the charge derivatives of the groups and of their contraction
   !> in seconds, summed over the threads
   real(wp), intent(out), optional :: times(2)

   real(wp) :: tgroup(2)
   integer, allocatable :: gstart(:)

   call get_list_groups(mol, list, gstart)
   call add_implicit_blocks_list(mol, cache, hessian, list, gstart, tgroup)
   if (present(times)) times(:) = tgroup

   ! Charge constraint term, sum_b (dX_b/dR_x - (dA/dR_x*q)_b) * uvec(b) is
   ! sum(u) * dlambda(:, x), where sum(u) = 1^T A^-1 1 is positive
   call get_constraint_factors(mol, cache, list, hessian%dlambda)
   hessian%uvecsum = sum(cache%uvec)

end subroutine get_implicit_hessian_list


!> Add the implicit blocks of the atoms x <= y on a half neighborlist with its
!> lower triangle in the cache for groups of consecutive atoms. The charge
!> derivatives of a group are contracted with the column y of the partial
!> derivatives at once, so each column is read once per group. The column y of
!> the partial derivatives holds the diagonal and the lower triangle in the row
!> of y and the upper triangle in the rows of the neighbors of lower index,
!> which is gathered once in the order of the lower triangle of the list. The
!> Hessian keeps the blocks of y >= x on a half list, only the group of atom x
!> writes the blocks (x, y), so the groups are independent.
subroutine add_implicit_blocks_list(mol, cache, hessian, list, gstart, times)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Cache with the inverse and the partial derivatives on the half list and
   !> its lower triangle
   type(mchrg_cache), intent(in) :: cache

   !> Molecular nuclear Hessian on a half list with allocated blocks
   type(mchrg_hessian), intent(inout) :: hessian

   !> Half neighborlist
   type(csr_list), intent(in) :: list

   !> First atom of each group, the last entry is the number of atoms plus one
   integer, intent(in) :: gstart(:)

   !> Time of the charge derivatives of the groups and of their contraction
   !> in seconds, summed over the threads
   real(wp), intent(out) :: times(2)

   integer :: igrp, ngrp, nbmax, numax, x0, ng, i, ic, xat, yat, nb, nu, n, u
   integer(i8) :: kat, lat, hat, tick(3), rate
   real(wp) :: dy(3), tpan, tcon
   integer, allocatable :: bloc(:), hloc(:), blist(:), ulist(:)
   integer(i8), allocatable :: qpos(:, :), jmap(:, :)
   real(wp), allocatable :: dlow(:, :), pan(:, :), qcol(:, :), acc(:, :, :)

   times(:) = 0.0_wp
   ngrp = size(gstart) - 1
   if (ngrp < 1) return
   nbmax = min(mol%nat, max_group * get_max_neighbors(mol%nat, list, cache%lowlist))
   numax = min(mol%nat, max_group * int(maxval(hessian%list%inl(2:mol%nat + 1) &
      & - hessian%list%inl(:mol%nat))))

   ! Upper triangle of the columns of the partial derivatives in the order of
   ! the lower triangle of the list
   allocate(dlow(3, size(cache%lowlist%nlat, kind=i8)))
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(mol, cache, dlow) private(xat, lat)
   do xat = 1, mol%nat
      do lat = cache%lowlist%inl(xat), cache%lowlist%inl(xat + 1) - 1
         dlow(:, lat) = cache%dabdrlist_ut(:, cache%lowlist%pos(lat))
      end do
   end do

   ! bloc and hloc hold the position of each atom in the unions of the group
   tpan = 0.0_wp
   tcon = 0.0_wp
   !$omp parallel default(none) &
   !$omp shared(mol, cache, hessian, list, dlow, gstart, ngrp, nbmax) &
   !$omp shared(numax) &
   !$omp private(igrp, x0, ng, i, ic, xat, yat, nb, nu, n, u, kat, lat, hat) &
   !$omp private(tick, rate, dy, bloc, hloc, blist, ulist, qpos, jmap, pan, qcol) &
   !$omp private(acc) &
   !$omp reduction(+:tpan, tcon)
   allocate(bloc(mol%nat), hloc(mol%nat), source=0)
   allocate(blist(nbmax), ulist(numax), qpos(max_group, nbmax), &
      & jmap(max_group, numax), pan(3 * max_group, nbmax), &
      & qcol(3 * max_group, nbmax), acc(3, 3, max_group))
   !$omp do schedule(runtime)
   do igrp = 1, ngrp
      call system_clock(tick(1), rate)
      x0 = gstart(igrp)
      ng = gstart(igrp + 1) - x0

      ! Charge derivatives dq_b/dR_x of the atoms x of the group for the union
      ! of their neighbors b, zero if b is not a neighbor of x
      call get_group_dqdr(cache, list, x0, ng, nb, bloc, blist, qpos, pan, qcol)
      call system_clock(tick(2))

      ! Union of the atoms y >= x in the rows of the Hessian of the group with
      ! the positions of the blocks (x, y), zero if not present
      nu = 0
      do i = 1, ng
         xat = x0 + i - 1
         do hat = hessian%list%inl(xat), hessian%list%inl(xat + 1) - 1
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

      ! acc(jc, ic, i) = sum_b (dX_b/dR_y,jc - (dA/dR_y,jc*q)_b) dq_b/dR_x,ic
      ! for the atoms x of the group, the block (x, y) is -acc
      do u = 1, nu
         yat = ulist(u)
         hloc(yat) = 0
         acc(:, :, :ng) = 0.0_wp
         kat = list%inl(yat)
         n = bloc(yat)
         if (n > 0) then
            dy = cache%dabdrlist_ut(:, kat)
            do i = 1, ng
               ic = 3 * (i - 1)
               acc(:, 1, i) = acc(:, 1, i) + dy * qcol(ic + 1, n)
               acc(:, 2, i) = acc(:, 2, i) + dy * qcol(ic + 2, n)
               acc(:, 3, i) = acc(:, 3, i) + dy * qcol(ic + 3, n)
            end do
         end if
         do kat = list%inl(yat) + 1, list%inl(yat + 1) - 1
            n = bloc(list%nlat(kat))
            if (n == 0) cycle
            dy = cache%dabdrlist_lt(:, kat)
            do i = 1, ng
               ic = 3 * (i - 1)
               acc(:, 1, i) = acc(:, 1, i) + dy * qcol(ic + 1, n)
               acc(:, 2, i) = acc(:, 2, i) + dy * qcol(ic + 2, n)
               acc(:, 3, i) = acc(:, 3, i) + dy * qcol(ic + 3, n)
            end do
         end do
         do lat = cache%lowlist%inl(yat), cache%lowlist%inl(yat + 1) - 1
            n = bloc(cache%lowlist%nlat(lat))
            if (n == 0) cycle
            dy = dlow(:, lat)
            do i = 1, ng
               ic = 3 * (i - 1)
               acc(:, 1, i) = acc(:, 1, i) + dy * qcol(ic + 1, n)
               acc(:, 2, i) = acc(:, 2, i) + dy * qcol(ic + 2, n)
               acc(:, 3, i) = acc(:, 3, i) + dy * qcol(ic + 3, n)
            end do
         end do
         do i = 1, ng
            hat = jmap(i, u)
            if (hat == 0) cycle
            if (yat == x0 + i - 1) then
               hessian%hesslist(:, :, hat) = hessian%hesslist(:, :, hat) &
                  & - 0.5_wp * (acc(:, :, i) + transpose(acc(:, :, i)))
            else
               hessian%hesslist(:, :, hat) = hessian%hesslist(:, :, hat) &
                  & - acc(:, :, i)
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
   deallocate(bloc, hloc, blist, ulist, qpos, jmap, pan, qcol, acc)
   !$omp end parallel
   times(:) = [tpan, tcon]

end subroutine add_implicit_blocks_list


!> Groups of consecutive atoms within a distance of the first atom of the group
subroutine get_atom_groups(mol, radius, max_group, gstart)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Largest distance of an atom of a group to its first atom
   real(wp), intent(in) :: radius

   !> Largest number of atoms in a group
   integer, intent(in) :: max_group

   !> First atom of each group, the last entry is the number of atoms plus one
   integer, allocatable, intent(out) :: gstart(:)

   integer :: ngrp, iat, ng
   integer, allocatable :: tmp(:)

   allocate(tmp(mol%nat + 1))
   ngrp = 0
   iat = 1
   do while (iat <= mol%nat)
      ngrp = ngrp + 1
      tmp(ngrp) = iat
      ng = 1
      do while (iat + ng <= mol%nat .and. ng < max_group)
         if (norm2(mol%xyz(:, iat + ng) - mol%xyz(:, iat)) > radius) exit
         ng = ng + 1
      end do
      iat = iat + ng
   end do
   tmp(ngrp + 1) = mol%nat + 1
   gstart = tmp(:ngrp + 1)

end subroutine get_atom_groups


!> Compute the molecular nuclear Hessian of an external function f(q, R) of the
!> charges with the adjoint vector p = S df/dq, the explicit part
!> p^T (d2X/dR2 - d2A/dR2*q) at fixed charges and the implicit part
!> G^T dq/dR + (dq/dR)^T G + (dq/dR)^T d2f/dq2 dq/dR from the response of the
!> charges with G = d2f/dRdq - dA/dR*p. The derivatives of f are taken at the
!> charges of the structure, the second derivative d2f/dR2 at fixed charges is
!> added by the caller. The charges, energy, gradient and virial of the
!> electrostatic energy of the same solve are optional.
subroutine get_external_hessian(self, mol, solver, cache, error, dfdq, d2fdq2, &
   & dfdrdq, hessian, qvec, energy, gradient, sigma, unit, verbosity)

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

   !> Mixed derivative of the external function w.r.t. positions and charges,
   !> dfdrdq(:, x, b) = d2f/dR_x dq_b
   real(wp), intent(in) :: dfdrdq(:, :, :)

   !> Molecular nuclear Hessian of the external function without d2f/dR2
   real(wp), intent(out) :: hessian(:, :)

   !> Optional atomic partial charges
   real(wp), intent(out), contiguous, optional :: qvec(:)

   !> Optional electrostatic energy
   real(wp), intent(inout), contiguous, optional :: energy(:)

   !> Optional energy gradient
   real(wp), intent(inout), contiguous, optional :: gradient(:, :)

   !> Optional virial
   real(wp), intent(inout), contiguous, optional :: sigma(:, :)

   !> Output unit
   integer, intent(in), optional :: unit

   !> Verbosity level
   integer, intent(in), optional :: verbosity

   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), allocatable :: dqdr(:, :, :), dqdL(:, :, :), dtmpdr(:, :, :)
   real(wp), allocatable :: yvec(:), padj(:), dfdq_loc(:)
   real(wp) :: uvecsum, yvecsum, scale
   integer :: iat, ndim, print_unit, verbosity_solve

   if (any(mol%periodic)) then
      call fatal_error(error, &
         & "Analytical Hessian is not available for periodic systems")
      return
   end if

   if (size(dfdq) /= mol%nat .or. any(shape(d2fdq2) /= [mol%nat, mol%nat]) &
      & .or. any(shape(dfdrdq) /= [3, mol%nat, mol%nat])) then
      call fatal_error(error, "External partial derivatives are wrong size")
      return
   end if

   verbosity_solve = 0
   if (present(verbosity)) verbosity_solve = verbosity
   print_unit = output_unit
   if (present(unit)) print_unit = unit

   ! Charges and their derivatives
   call self%update(mol, cache, trans, grad=.false., hess=.true.)
   allocate(dqdL(3, 3, mol%nat), dqdr(3, mol%nat, mol%nat))
   call self%solve(mol, solver, cache, error, energy=energy, gradient=gradient, &
      & sigma=sigma, qvec=qvec, dqdr=dqdr, dqdL=dqdL, verbosity=verbosity, &
      & unit=unit)
   if (allocated(error)) return

   ! Adjoint vector p = S df/dq
   if (solver%need_pos_def) then
      ndim = mol%nat
      allocate(yvec(ndim))
      do iat = 1, mol%nat
         yvec(iat) = dfdq(iat) / cache%amat(iat, iat)
      end do

      ! Unconstrained response: J*yvec = dfdq
      call print_adjoint_message(print_unit, verbosity_solve)
      call solver%solve(amat=cache%amat, xvec=dfdq, vrhs=yvec, &
         & new_unit=print_unit, error=error)
      if (allocated(error)) return

      ! Projection of uvec on yvec
      yvecsum = sum(yvec)
      uvecsum = sum(cache%uvec)
      scale = yvecsum / (uvecsum + eps)
      allocate(padj(mol%nat))
      padj = yvec - scale * cache%uvec
   else
      ndim = mol%nat + 1
      allocate(dfdq_loc(ndim), source=0.0_wp)
      dfdq_loc(:mol%nat) = dfdq
      allocate(padj(ndim))

      ! Direct solution: J*yvec = dfdq
      call print_adjoint_message(print_unit, verbosity_solve)
      call solver%solve(amat=cache%amat, xvec=dfdq_loc, vrhs=padj, &
         & new_unit=print_unit, error=error)
      if (allocated(error)) return
   end if

   hessian(:, :) = 0.0_wp

   ! Mixed derivatives of the Lagrangian, G^T = d2f/dRdq - dA/dR*p
   call new_partial_derivs(mol, ndim, cache, dafdr=.true.)
   call self%get_coulomb_derivs(mol, ndim, cache, alpha=-1.0_wp, &
      & p=padj(:mol%nat), dafdr=.true.)
   cache%dafdr = cache%dafdr + dfdrdq

   ! Implicit part G^T dq/dR + (dq/dR)^T G
   call gemm(reshape(cache%dafdr, [3 * mol%nat, mol%nat]), &
      & reshape(dqdr, [3 * mol%nat, mol%nat]), hessian, transb='t', &
      & alpha=+1.0_wp, beta=1.0_wp)

   call gemm(reshape(dqdr, [3 * mol%nat, mol%nat]), &
      & reshape(cache%dafdr, [3 * mol%nat, mol%nat]), hessian, transb='t', &
      & alpha=+1.0_wp, beta=1.0_wp)

   ! Implicit part (dq/dR)^T d2f/dq2 dq/dR
   allocate(dtmpdr(3, mol%nat, mol%nat))
   call gemm(dqdr, d2fdq2, dtmpdr)
   call gemm(reshape(dtmpdr, [3 * mol%nat, mol%nat]), &
      & reshape(dqdr, [3 * mol%nat, mol%nat]), hessian, transb='t', &
      & alpha=+1.0_wp, beta=1.0_wp)

   ! Second derivatives of p^T (X - A*q) at fixed charges and adjoint vector
   call self%get_explicit_hessian(mol, cache, padj(:mol%nat), hessian, &
      & error, alpha=-1.0_wp, beta=1.0_wp)
   if (allocated(error)) return

end subroutine get_external_hessian


!> Semi-numerical Hessian and pressure tensor from central differences of
!> the analytical energy gradient and virial
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


!> Semi-numerical Hessian in dense symmetric storage, obtained by expanding the
!> packed lower triangle
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


!> Accumulate one column of the Cartesian Hessian into packed lower-triangle
!> storage. Every off-diagonal element is reached once from each of the two
!> displacements that contribute to it, so both estimates are averaged.
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


!> Energy gradient and virial for a displaced structure using a fresh cache
!> and, if requested, a rebuilt neighborlist
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

   !> Reference neighborlist, only its settings are reused
   type(csr_list), intent(in), optional :: list

   !> Output unit
   integer, intent(in), optional :: unit

   !> Verbosity level
   integer, intent(in), optional :: verbosity

   type(mchrg_cache) :: cache
   type(csr_list), allocatable :: list_work
   real(wp), allocatable :: trans(:, :)

   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   call get_lattice_points(mol%periodic, mol%lattice, self%ncoord%cutoff, trans)

   if (present(list)) then
      allocate(list_work)
      if (any(mol%periodic)) then
         call new_csr_list(list_work, mol, error, cache%wsc, list%cutoff)
      else
         call new_csr_list(list_work, mol, error, cutoff=list%cutoff)
      end if
      if (allocated(error)) return
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

   !> Optional derivative of local atomic partial charges w.r.t. atomic positions
   !> in the CSR list layout: d(qloc_j)/d(r_i) for each list entry, with
   !> d(qloc_i)/d(r_i) stored in the diagonal entry list%inl(i)
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

end module multicharge_model_type
