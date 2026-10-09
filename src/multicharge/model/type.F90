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
   use mctc_csrlist, only : csr_list, spgemv_csr, spsymv_csr, spmm_csr, new_csr_list
   use multicharge_blas, only : gemv, symv, gemm
   use multicharge_model_cache, only : mchrg_cache
   use multicharge_solver_type, only : mchrg_solver_type
   use multicharge_solver_cg, only : cg_solver

   implicit none
   private

   public :: mchrg_model_type, get_dir_trans, get_rec_trans, hess_index

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


      !> Calculate Coulomb matrix
      procedure(get_coulomb_matrix), deferred :: get_coulomb_matrix

      !> Calculate alpha * dA/dR*q + beta * dX/dR
      procedure(get_partial_derivs), deferred :: get_partial_derivs

      !> Calculate capacitance-corrected electronegativity derivatives
      procedure(get_grad), deferred :: get_grad

   end type mchrg_model_type

   abstract interface
      !> Update model-dependent quantities and cache
      subroutine update(self, mol, cache, trans, grad, list)
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

      !> Linear combination of the Coulomb matrix derivatives contracted with the
      !> charges and the electronegativity vector derivatives,
      !> alpha * dA/dR*q + beta * dX/dR, w.r.t. positions and strain
      subroutine get_partial_derivs(self, mol, ndim, cache, alpha, beta, list)
         import :: mchrg_model_type, structure_type, mchrg_cache, csr_list, wp

         !> Multicharge model type
         class(mchrg_model_type), intent(in) :: self

         !> Structure type
         type(structure_type), intent(in) :: mol

         !> System size
         integer, intent(in) :: ndim

         !> Multicharge cache holding the partial derivatives
         type(mchrg_cache), intent(inout) :: cache

         !> Multiplier of the Coulomb matrix derivatives, dA/dR*q
         real(wp), intent(in) :: alpha

         !> Multiplier of the electronegativity vector derivatives, dX/dR
         real(wp), intent(in) :: beta

         !> Multicharge neighborlist type (complete)
         type(csr_list), intent(in), optional :: list
      end subroutine get_partial_derivs

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

   ! Charge derivatives are evaluated row by row on a complete neighborlist
   if (cpq .and. present(list)) then
      if (.not. list%complete) then
         call fatal_error(error, "Charge derivatives require a complete neighborlist")
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
      call self%get_partial_derivs(mol, ndim, cache, alpha=-1.0_wp, beta=1.0_wp, &
         & list=list)
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
         & unit=print_unit)
      ! pop cpq timer
      call timer%pop
      call print_gradient_time(print_unit, verbosity_solve, timer%get("cpq"))
   end if

   ! pop total solve timer
   call timer%pop
   call print_total_time(print_unit, verbosity_solve, timer%get("total"))

end subroutine solve


!> Derivatives of the partial charges w.r.t. atomic positions and lattice vectors
subroutine get_q_derivs(mol, solver, cache, error, ndim, dqdr, dqdL, list, unit)

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

   real(wp), allocatable :: daqxdr(:, :, :), daqxdL(:, :, :), ainv(:, :)
   real(wp) :: scale, uvecsum
   integer :: jat

   ! Compressed position derivatives of an earlier call are outdated
   if (allocated(cache%dqdrlist)) deallocate(cache%dqdrlist)
   if (allocated(cache%dqdrscal)) deallocate(cache%dqdrscal)

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
         if (allocated(cache%ainvlist)) then
            if (size(cache%ainvlist, kind=i8) /= size(list%nlat, kind=i8)) &
               & deallocate(cache%ainvlist)
         end if
         if (.not. allocated(cache%ainvlist)) then
            allocate(cache%ainvlist(size(list%nlat, kind=i8)))
         end if
         call solver%invert_list(cache%alist, list, cache%ainvlist, &
            & xyz=mol%xyz, new_unit=unit, error=error)
         if (allocated(error)) return

         call get_q_derivs_list(mol, solver, cache, list, dqdL, error, unit)
         if (allocated(error)) return

         if (present(dqdr)) then
            !$omp parallel do default(none) schedule(runtime) &
            !$omp shared(mol, cache, dqdr) private(jat)
            do jat = 1, mol%nat
               call cache%get_dqdr_row(jat, dqdr(:, :, jat))
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


!> Charge derivatives from the inverse S = J^-1 on the pattern of the complete
!> neighborlist. The position derivatives dq_j/dR_i = (S B)_ji - u_j w_i / sum(u)
!> with B = dB/dR - dA/dR*q and w = B^T u are kept in compressed-row storage:
!> the sparse product S B in cache%dqdrlist on the pattern cache%dqdrpat and
!> the factors w / sum(u) of the charge constraint term in cache%dqdrscal.
subroutine get_q_derivs_list(mol, solver, cache, list, dqdL, error, unit)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Conjugate-gradient solver with the threshold of the inverse
   class(cg_solver), intent(in) :: solver

   !> Cache with the inverse on the neighborlist, receives the compressed
   !> position derivatives
   type(mchrg_cache), intent(inout) :: cache

   !> Complete neighborlist
   type(csr_list), intent(in) :: list

   !> Derivative of the atomic partial charges w.r.t. lattice vectors
   real(wp), intent(out) :: dqdL(:, :, :)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Output unit
   integer, intent(in), optional :: unit

   integer :: iat, jat, ic, info
   integer(i8) :: kat, lat, nkeep, nnzq
   real(wp) :: thr, uvecsum, scale
   real(wp) :: wlat(3, 3)
   ! Inverse on its pattern without the elements below the threshold
   integer(i8), allocatable :: sinl(:)
   integer, allocatable :: snlat(:)
   real(wp), allocatable :: sval(:)
   ! Placeholders for the product values not referenced by the pattern stage
   integer :: jdum(0)
   real(wp) :: cdum(0)

   ! Elements of the inverse below the threshold are dropped
   thr = 0.0_wp
   do iat = 1, mol%nat
      thr = max(thr, cache%ainvlist(list%inl(iat)))
   end do
   thr = solver%ainvthr * thr
   allocate(sinl(mol%nat + 1))
   sinl(1) = 1
   do iat = 1, mol%nat
      sinl(iat + 1) = sinl(iat)
      do kat = list%inl(iat), list%inl(iat + 1) - 1
         if (abs(cache%ainvlist(kat)) > thr) sinl(iat + 1) = sinl(iat + 1) + 1
      end do
   end do
   nkeep = sinl(mol%nat + 1) - 1
   allocate(snlat(nkeep), sval(nkeep))
   lat = 1
   do iat = 1, mol%nat
      do kat = list%inl(iat), list%inl(iat + 1) - 1
         if (abs(cache%ainvlist(kat)) <= thr) cycle
         snlat(lat) = list%nlat(kat)
         sval(lat) = cache%ainvlist(kat)
         lat = lat + 1
      end do
   end do

   ! Sparse product S B, B holds dB_k/dR_i in row k, the pattern is formed once
   ! and the values for each Cartesian component
   if (allocated(cache%dqdrpat%inl)) deallocate(cache%dqdrpat%inl)
   if (allocated(cache%dqdrpat%nlat)) deallocate(cache%dqdrpat%nlat)
   cache%dqdrpat%complete = .true.
   allocate(cache%dqdrpat%inl(mol%nat + 1))
   call spmm_csr("N", 1, 0, mol%nat, mol%nat, mol%nat, sval, snlat, sinl, &
      & cache%dabdrlist(1, :), list%nlat, list%inl, cdum, jdum, &
      & cache%dqdrpat%inl, 0_i8, info)
   nnzq = cache%dqdrpat%inl(mol%nat + 1) - 1
   allocate(cache%dqdrpat%nlat(nnzq), cache%dqdrlist(3, nnzq))
   do ic = 1, 3
      call spmm_csr("N", 2, 0, mol%nat, mol%nat, mol%nat, sval, snlat, sinl, &
         & cache%dabdrlist(ic, :), list%nlat, list%inl, cache%dqdrlist(ic, :), &
         & cache%dqdrpat%nlat, cache%dqdrpat%inl, nnzq, info)
      if (info /= 0) then
         call fatal_error(error, "Sparse product of the charge derivatives failed")
         return
      end if
   end do

   if (present(unit) .and. solver%verbosity > 0) then
      write(unit, '(a, 1x, i0, a, i0)') "Inverse elements on neighborlist :", &
         & nkeep, " of ", size(list%nlat, kind=i8)
      write(unit, '(a, 1x, i0, a, i0)') "Charge derivative elements       :", &
         & nnzq, " of ", int(mol%nat, i8)**2
      write(unit, '(a)') ''
   end if

   ! Contractions with the constraint response, w = B^T u
   uvecsum = sum(cache%uvec)
   allocate(cache%dqdrscal(3, mol%nat), source=0.0_wp)
   wlat(:, :) = 0.0_wp
   do iat = 1, mol%nat
      do kat = list%inl(iat), list%inl(iat + 1) - 1
         jat = list%nlat(kat)
         cache%dqdrscal(:, jat) = cache%dqdrscal(:, jat) &
            & + cache%dabdrlist(:, kat) * cache%uvec(iat)
      end do
      wlat(:, :) = wlat + cache%dabdL(:, :, iat) * cache%uvec(iat)
   end do
   cache%dqdrscal(:, :) = cache%dqdrscal / (uvecsum + eps)

   ! Lattice derivatives B_L^T S and their charge constraint term
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(mol, cache, sinl, snlat, sval, uvecsum, wlat, dqdL) &
   !$omp private(iat, jat, kat, scale)
   do jat = 1, mol%nat
      scale = cache%uvec(jat) / (uvecsum + eps)
      dqdL(:, :, jat) = -scale * wlat
      do kat = sinl(jat), sinl(jat + 1) - 1
         iat = snlat(kat)
         dqdL(:, :, jat) = dqdL(:, :, jat) + sval(kat) * cache%dabdL(:, :, iat)
      end do
   end do

end subroutine get_q_derivs_list


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

   call self%get_grad(mol, cache, padj, dfdr, dfdL, alpha=-1.0_wp, &
      & beta=1.0_wp, list=list)

   ! pop dfdr timer
   call timer%pop
   call print_gradient_time(print_unit, verbosity_solve, &
      & timer%get("external_gradient"))

end subroutine get_external_gradient


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

   !> Virial derivatives w.r.t. positions (3, 3, 3, nat) or strain (3, 3, 3, 3)
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

end module multicharge_model_type
