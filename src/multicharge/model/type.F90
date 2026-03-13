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

!> @file multicharge/model/type.f90
!> Provides a general base class for the charge models

#ifndef IK
#define IK i4
#endif

!> General charge model
module multicharge_model_type
   use iso_fortran_env, only : output_unit
   use mctc_env, only: timer_type, format_time, error_type, fatal_error, wp,  ik => IK
   use mctc_io, only: structure_type
   use mctc_io_constants, only: pi
   use mctc_io_math, only: matinv_3x3
   use mctc_cutoff, only: get_lattice_points
   use mctc_ncoord, only: ncoord_type
   use multicharge_blas, only: gemv, symv, gemm
   use multicharge_lapack, only: sytrf, sytrs
   use multicharge_wignerseitz, only: wignerseitz_cell_type, new_wignerseitz_cell
   use multicharge_model_cache, only: model_cache, cache_container
   use multicharge_solver_type, only: mchrg_solver_type
   
   implicit none
   private

   public :: mchrg_model_type, get_dir_trans, get_rec_trans

   !> Abstract multicharge model type
   type, abstract :: mchrg_model_type
      !> Size of a system
      integer, allocatable :: ndim
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
      real(wp), allocatable :: kcnrad
      !> Coordination number
      class(ncoord_type), allocatable :: ncoord
      !> Electronegativity weighted CN for local charge
      class(ncoord_type), allocatable :: ncoord_en
   contains
      !> Solve linear equations for the charge model
      procedure :: solve
      !> Calculate local charges from electronegativity weighted CN
      procedure :: local_charge
      !> Update cache
      procedure(update), deferred :: update
      !> Calculate right-hand side (electronegativity)
      procedure(get_xvec), deferred :: get_xvec
      !> Calculate xvec Gradients
      procedure(get_xvec_derivs), deferred :: get_xvec_derivs
      !> Calculate Coulomb matrix
      procedure(get_coulomb_matrix), deferred :: get_coulomb_matrix
      !> Calculate Coulomb matrix derivatives
      procedure(get_coulomb_derivs), deferred :: get_coulomb_derivs
   end type mchrg_model_type

   abstract interface
      subroutine update(self, mol, cache, ndim, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL)
         import :: mchrg_model_type, structure_type, cache_container, wp
         class(mchrg_model_type), intent(in) :: self
         type(structure_type), intent(in) :: mol
         type(cache_container), intent(inout) :: cache
         integer, intent(in) :: ndim   
         real(wp), intent(in) :: cn(:)
         real(wp), intent(in), optional :: qloc(:)
         real(wp), intent(in), optional :: dcndr(:, :, :)
         real(wp), intent(in), optional :: dcndL(:, :, :)
         real(wp), intent(in), optional :: dqlocdr(:, :, :)
         real(wp), intent(in), optional :: dqlocdL(:, :, :)
      end subroutine update

      subroutine get_coulomb_matrix(self, mol, cache, amat)
         import :: mchrg_model_type, structure_type, cache_container, wp
         class(mchrg_model_type), intent(in) :: self
         type(structure_type), intent(in) :: mol
         type(cache_container), intent(inout) :: cache
         real(wp), intent(out) :: amat(:, :)
      end subroutine get_coulomb_matrix

      subroutine get_coulomb_derivs(self, mol, cache, qvec, dadr, dadL, atrace)
         import :: mchrg_model_type, structure_type, cache_container, wp
         class(mchrg_model_type), intent(in) :: self
         type(structure_type), intent(in) :: mol
         type(cache_container), intent(inout) :: cache
         real(wp), intent(in) :: qvec(:)
         real(wp), intent(out) :: dadr(:, :, :), dadL(:, :, :), atrace(:, :)
      end subroutine get_coulomb_derivs

      subroutine get_xvec(self, mol, cache, xvec)
         import :: mchrg_model_type, cache_container, structure_type, wp
         class(mchrg_model_type), intent(in) :: self
         type(structure_type), intent(in) :: mol
         type(cache_container), intent(inout) :: cache
         real(wp), intent(out) :: xvec(:)
      end subroutine get_xvec

      subroutine get_xvec_derivs(self, mol, cache, dxdr, dxdL)
         import :: mchrg_model_type, structure_type, cache_container, wp
         class(mchrg_model_type), intent(in) :: self
         type(structure_type), intent(in) :: mol
         type(cache_container), intent(inout) :: cache
         real(wp), intent(out), contiguous :: dxdr(:, :, :)
         real(wp), intent(out), contiguous :: dxdL(:, :, :)
      end subroutine get_xvec_derivs    
   end interface

   real(wp), parameter :: twopi = 2 * pi
   real(wp), parameter :: eps = tiny(1.0_wp)

contains

subroutine get_dir_trans(lattice, trans)
   real(wp), intent(in) :: lattice(:, :)
   real(wp), allocatable, intent(out) :: trans(:, :)
   integer, parameter :: rep(3) = [2, 2, 2]

   call get_lattice_points(lattice, rep, .true., trans)

end subroutine get_dir_trans

subroutine get_rec_trans(lattice, trans)
   real(wp), intent(in) :: lattice(:, :)
   real(wp), allocatable, intent(out) :: trans(:, :)
   integer, parameter :: rep(3) = [2, 2, 2]
   real(wp) :: rec_lat(3, 3)

   rec_lat = twopi * transpose(matinv_3x3(lattice))
   call get_lattice_points(rec_lat, rep, .false., trans)

end subroutine get_rec_trans

subroutine solve(self, mol, solver, error, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL, &
   & energy, gradient, sigma, qvec, dqdr, dqdL, dfdq, verbosity, unit)
   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in):: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> The solver instance
   class(mchrg_solver_type), intent(in) :: solver
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Coordination number
   real(wp), intent(in), contiguous :: cn(:)
   !> Local atomic partial charges
   real(wp), intent(in), contiguous :: qloc(:)
   !> Optional derivative of the coordination number w.r.t. atomic positions
   real(wp), intent(in), contiguous, optional :: dcndr(:, :, :)
   !> Optional derivative of the coordination number w.r.t. lattice vectors
   real(wp), intent(in), contiguous, optional :: dcndL(:, :, :)
   !> Optional derivative of the local atomic partial charges w.r.t. atomic positions
   real(wp), intent(in), contiguous, optional :: dqlocdr(:, :, :)
   !> Optional derivative of the local atomic partial charges w.r.t. lattice vectors
   real(wp), intent(in), contiguous, optional :: dqlocdL(:, :, :)
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
   !> Optional derivative of the electrostatic energy w.r.t. atomic partial charges
   real(wp), intent(in), contiguous, optional :: dfdq(:)
   !> Optional print verbossity number input flag
   integer, intent(in), optional :: verbosity
   !> Output unit
   integer, intent(in), optional :: unit

   type(cache_container), allocatable :: cache
   logical :: main_mode, dfdr_mode

   allocate(cache)

   ! Determine if we need the forward solution (charges, energy, gradients, charge derivatives)
   main_mode = present(qvec) .or. present(energy) .or. present(gradient) .or. &
               present(sigma) .or. present(dqdr) .or. present(dqdL)
   dfdr_mode = present(dfdq)

   if (main_mode) then
      ! First phase: compute energy, charges, and possibly forward gradients/derivatives
      call solve_main(self, mol, solver, error, cn, qloc, dcndr, dcndL, &
           dqlocdr, dqlocdL, energy, gradient, sigma, qvec, dqdr, dqdL, &
           verbosity, unit, cache)
      if (allocated(error)) return
   end if

   if (dfdr_mode) then
      ! Second phase: compute adjoint gradient (requires cache from forward solve)
      if (.not. main_mode) then
         call fatal_error(error, "Adjoint gradients requested but forward solve was not performed; cache not available")
         return
      end if
      call solve_dfdr(self, mol, solver, error, dfdq, gradient, unit, verbosity, cache)
   end if

end subroutine solve

!> Main solve: obtain charges, energy, gradients, and store data for adjoint
subroutine solve_main(self, mol, solver, error, cn, qloc, dcndr, dcndL, &
      & dqlocdr, dqlocdL, energy, gradient, sigma, qvec, dqdr, dqdL, &
      & verbosity, unit, cache)
   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Solver instance
   class(mchrg_solver_type), intent(in) :: solver
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Coordination number
   real(wp), intent(in) :: cn(:)
   !> Local atomic partial charges
   real(wp), intent(in) :: qloc(:)
   !> Derivatives of coordination number w.r.t. positions
   real(wp), intent(in), optional :: dcndr(:, :, :)
   !> Derivatives of coordination number w.r.t. lattice vectors
   real(wp), intent(in), optional :: dcndL(:, :, :)
   !> Derivatives of local charges w.r.t. positions
   real(wp), intent(in), optional :: dqlocdr(:, :, :)
   !> Derivatives of local charges w.r.t. lattice vectors
   real(wp), intent(in), optional :: dqlocdL(:, :, :)
   !> Electrostatic energy (incremented)
   real(wp), intent(inout), optional :: energy(:)
   !> Gradient of electrostatic energy (incremented)
   real(wp), intent(inout), optional :: gradient(:, :)
   !> Stress tensor of electrostatic energy (incremented)
   real(wp), intent(inout), optional :: sigma(:, :)
   !> Atomic partial charges
   real(wp), intent(out), optional :: qvec(:)
   !> Derivatives of charges w.r.t. positions
   real(wp), intent(out), optional :: dqdr(:, :, :)
   !> Derivatives of charges w.r.t. lattice vectors
   real(wp), intent(out), optional :: dqdL(:, :, :)
   !> Verbosity level
   integer, intent(in), optional :: verbosity
   !> Output unit
   integer, intent(in), optional :: unit
   !> Cache to store data for later adjoint calculation
   type(cache_container), intent(inout) :: cache

   integer :: ndim, iat, ic, jc
   logical :: add_lagr, grad, cpq, dcn
   real(wp), allocatable :: xvec(:), vrhs(:), amat(:, :), diag(:)
   real(wp), allocatable :: ainv(:, :), jmat(:, :)
   real(wp), allocatable :: unitvec(:), uvec(:), vvec(:)
   real(wp) :: uvecsum, vvecsum, lambda
   real(wp), allocatable :: dadr(:, :, :), dadL(:, :, :), atrace(:, :)
   real(wp), allocatable :: dxdr(:, :, :), dxdL(:, :, :)
   real(wp), allocatable :: derivsum(:, :, :)
   type(timer_type) :: timer
   integer :: verbosity_solve, print_unit

   dcn = present(dcndr) .and. present(dcndL)
   grad = present(gradient) .and. present(sigma) .and. dcn
   cpq = present(dqdr) .and. present(dqdL) .and. dcn

   verbosity_solve = 0
   if (present(verbosity)) verbosity_solve = verbosity
   print_unit = output_unit
   if (present(unit)) print_unit = unit

   ! Determine system size and constraint handling
   if (solver%need_pos_def) then
      ndim = mol%nat
      add_lagr = .false.
   else
      ndim = mol%nat + 1
      add_lagr = .true.
   end if

   ! Update model cache
   call self%update(mol, cache, ndim, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL)

   call timer%push("total")
   call timer%push("setup")

   ! Coulomb matrix
   allocate(amat(ndim, ndim))
   call self%get_coulomb_matrix(mol, cache, amat)
   ! Store for later adjoint use
   cache%amat = amat

   ! Right‑hand side (electronegativity vector)
   allocate(xvec(ndim))
   call self%get_xvec(mol, cache, xvec)
   allocate(vrhs(mol%nat+1))

   call timer%pop
   call print_solve_header(print_unit, verbosity_solve, timer)

   ! Solve the linear system
   if (add_lagr) then
      vrhs = xvec
      ainv = amat
      call solver%solve(amat, xvec, vrhs, ainv=ainv, cpq=cpq, &
         & new_unit=print_unit, error=error)
      jmat = amat(:mol%nat, :mol%nat)
   else
      ! Unconstrained case: build response vectors
      allocate(vvec(mol%nat), uvec(mol%nat), unitvec(mol%nat), diag(mol%nat))
      do iat = 1, mol%nat
         diag(iat) = amat(iat, iat)
      end do
      unitvec = 1.0_wp
      uvec = 1.0_wp / (diag + eps)
      vvec = -xvec(:mol%nat) / (diag + eps)

      call print_constrained_system_message(print_unit, verbosity_solve, 'u')
      call solver%solve(amat, unitvec, uvec, new_unit=print_unit, error=error)
      call print_constrained_system_message(print_unit, verbosity_solve, 'v')
      call solver%solve(amat, -xvec(:mol%nat), vvec, new_unit=print_unit, error=error)

      uvecsum = sum(uvec)
      vvecsum = sum(vvec)
      lambda = -(mol%charge + vvecsum) / (uvecsum + eps)
      vrhs(:mol%nat) = -vvec - lambda * uvec
      vrhs(mol%nat+1) = lambda

      jmat = amat(:mol%nat, :mol%nat)
      ! Store uvec for later adjoint
      cache%uvec = uvec
   end if

   ! Store charges if requested
   if (present(qvec)) qvec(:) = vrhs(:mol%nat)

   ! Electrostatic energy
   if (present(energy)) then
      call symv(jmat, vrhs(:mol%nat), xvec(:mol%nat), &
         & alpha=0.5_wp, beta=-1.0_wp, uplo='l')
      energy(:) = energy(:) + vrhs(:mol%nat) * xvec(:mol%nat)
   end if

   if (present(gradient) .or. present(dqdr)) then
      call print_gradient_header(print_unit, verbosity_solve)
      allocate(dadr(3, mol%nat, mol%nat), dadL(3, 3, mol%nat), atrace(3, mol%nat))
      allocate(dxdr(3, mol%nat, mol%nat), dxdL(3, 3, mol%nat))
      call self%get_xvec_derivs(mol, cache, dxdr, dxdL)
      call self%get_coulomb_derivs(mol, cache, vrhs, dadr, dadL, atrace)
      do iat = 1, mol%nat
         dadr(:, iat, iat) = atrace(:, iat) + dadr(:, iat, iat)
      end do
      ! Build and store derivsum = dxdr - dadr
      allocate(derivsum(3, mol%nat, mol%nat))
      do iat = 1, mol%nat
         derivsum(:, :, iat) = dxdr(:, :, iat) - dadr(:, :, iat)
      end do
      cache%derivsum = derivsum
   end if

   ! Forward gradients if requested
   if (grad) then
      call timer%push("gradient")
      call gemv(dadr, vrhs(:mol%nat), gradient, beta=1.0_wp, alpha=0.5_wp)
      call gemv(dxdr, vrhs(:mol%nat), gradient, beta=1.0_wp, alpha=-1.0_wp)
      call gemv(dadL, vrhs, sigma, beta=1.0_wp, alpha=0.5_wp)
      call gemv(dxdL, vrhs, sigma, beta=1.0_wp, alpha=-1.0_wp)
      call timer%pop
      call print_gradient_time(print_unit, verbosity_solve, timer)
   end if

   ! Charge derivatives if requested
   if (cpq) then
      call print_gradient_header(print_unit, verbosity_solve)
      ! Constrained case: use inverse from solver
      do iat = 1, mol%nat
         dadr(:, :, iat) = -dxdr(:, :, iat) + dadr(:, :, iat)
         dadL(:, :, iat) = -dxdL(:, :, iat) + dadL(:, :, iat)
      end do
      call gemm(dadr, ainv(:, :mol%nat), dqdr, alpha=-1.0_wp)
      call gemm(dadL, ainv(:, :mol%nat), dqdL, alpha=-1.0_wp)

      ! pop gradient
      call timer%pop
      call print_gradient_time(print_unit, verbosity_solve, timer)
   end if

   call timer%pop
   call print_total_time(print_unit, verbosity_solve, timer)

end subroutine solve_main

!> Adjoint gradient calculation using cached data from solve_forward
subroutine solve_dfdr(self, mol, solver, error, dfdq, gradient, unit, verbosity, cache)
   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Solver instance
   class(mchrg_solver_type), intent(in) :: solver
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Derivative of the objective w.r.t. atomic partial charges
   real(wp), intent(in) :: dfdq(:)
   !> Gradient of the objective (incremented)
   real(wp), intent(inout) :: gradient(:, :)
   !> Output unit
   integer, intent(in), optional :: unit
   !> Verbosity level
   integer, intent(in), optional :: verbosity
   !> Cache containing amat, derivsum, and (if allocated) uvec
   type(cache_container), intent(in) :: cache

   integer :: iat
   real(wp), allocatable :: jmat(:, :), diag(:), yvec(:), uvec(:), unitvec(:), padj(:)
   real(wp) :: uvecsum, yvecsum, scale
   integer :: print_unit, verbosity_solve
   type(timer_type) :: timer

   verbosity_solve = 0
   if (present(verbosity)) verbosity_solve = verbosity
   print_unit = output_unit
   if (present(unit)) print_unit = unit

   call print_gradient_header(print_unit, verbosity_solve)
   call timer%push("gradient")

   ! Retrieve stored data
   if (.not. allocated(cache%amat) .or. .not. allocated(cache%derivsum)) then
      call fatal_error(error, "Cache does not contain required arrays (amat, derivsum)")
      return
   end if

   jmat = cache%amat(:mol%nat, :mol%nat) 

   if (solver%need_pos_def) then
      allocate(diag(mol%nat))
      do iat = 1, mol%nat
         diag(iat) = jmat(iat, iat)
      end do
      ! Unconstrained case: use stored uvec
      if (.not. allocated(cache%uvec)) then
         write(print_unit, '(a)') "[WARN] Cache does not have unit vector response, calculat it on-the-fly."
         allocate(unitvec(mol%nat), source = 1.0_wp)
         allocate(uvec(mol%nat), source = 1.0_wp / (diag + eps))
         call solver%solve(jmat, unitvec, uvec, new_unit=print_unit, error=error)
      else 
         uvec = cache%uvec
      end if
      
      allocate(yvec(mol%nat))
      yvec = dfdq
      call solver%solve(jmat, dfdq, yvec, error=error)
      if (allocated(error)) return

      yvecsum = sum(yvec)
      uvecsum = sum(uvec)
      scale = yvecsum / (uvecsum + eps)
      allocate(padj(mol%nat))
      padj = yvec - scale * uvec

      call gemv(cache%derivsum, padj, gradient, alpha=1.0_wp, beta=1.0_wp)

   else
      ! Constrained case: solve for uvec and yvec using jmat
      allocate(unitvec(mol%nat), source = 1.0_wp)
      allocate(uvec(mol%nat), source = unitvec)
      call solver%solve(jmat, unitvec, uvec, new_unit=print_unit, error=error)
      if (allocated(error)) return
      uvecsum = sum(uvec)

      yvec = dfdq
      call solver%solve(jmat, dfdq, yvec, new_unit=print_unit, error=error)
      if (allocated(error)) return
      yvecsum = sum(yvec)

      scale = yvecsum / (uvecsum + eps)
      allocate(padj(mol%nat))
      padj = yvec - scale * uvec

      call gemv(cache%derivsum, padj, gradient, alpha=1.0_wp, beta=1.0_wp)
   end if

   ! pop gradient
   call timer%pop
   call print_gradient_time(print_unit, verbosity_solve, timer)

end subroutine solve_dfdr

!> Local charges calculation
subroutine local_charge(self, mol, trans, qloc, dqlocdr, dqlocdL)
   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points
   real(wp), intent(in) :: trans(:, :)
   !> Local atomic partial charges
   real(wp), intent(out) :: qloc(:)
   !> Optional derivative of local atomic partial charges w.r.t. atomic positions
   real(wp), intent(out), optional :: dqlocdr(3, mol%nat, mol%nat)
   !> Optional derivative of local atomic partial charges w.r.t. lattice vectors
   real(wp), intent(out), optional :: dqlocdL(3, 3, mol%nat)

   qloc = 0.0_wp
   if (present(dqlocdr) .and. present(dqlocdL)) then
      dqlocdr = 0.0_wp
      dqlocdL = 0.0_wp
   end if
   ! Get the electronegativity weighted CN for local charge
   ! Derivatives depend only in this CN
   if (allocated(self%ncoord_en)) then
      call self%ncoord_en%get_coordination_number(mol, trans, qloc, dqlocdr, dqlocdL)
   end if

   ! Distribute the total charge equally
   qloc = qloc + mol%charge / real(mol%nat, wp)

end subroutine local_charge

!> Print header for charge equilibration solver
subroutine print_solve_header(unit, verbosity, timer)
   integer, intent(in) :: unit, verbosity
   type(timer_type), intent(in) :: timer

   if (verbosity > 0) then
      write(unit, '(54("-"))')
      write(unit, '(13x, a)') "Charge equilibration solver"
      write(unit, '(54("-"))')
      write(unit, '(a)') ''
      if (verbosity > 1) then
         write(unit, '(a, 1x, a)') "Setup time : ", &
            format_time(timer%get("setup"))
         write(unit, '(a)') ''
      end if
   end if
end subroutine print_solve_header

!> Print header for gradient calculations
subroutine print_gradient_header(unit, verbosity)
   integer, intent(in) :: unit, verbosity

   if (verbosity > 0) then
      write(unit, '(54("-"))')
      write(unit, '(17x, a)') "Gradient Calculations"
      write(unit, '(54("-"))')
      write(unit, '(a)') ''
   end if
end subroutine print_gradient_header

!> Print message for constrained system solves
subroutine print_constrained_system_message(unit, verbosity, vector)
   integer, intent(in) :: unit, verbosity
   character, intent(in) :: vector

   if (verbosity > 0) then
      if (vector == 'u') then
         write(unit, '(a)') 'Solving constrained system: J*u = 1'
      else if (vector == 'v') then
         write(unit, '(a)') 'Solving constrained system: J*v = chi'
      else if (vector == 'y') then
         write(unit, '(a)') 'Solving derivative constrained system: J*y = df/dq'
      end if
      write(unit, '(a)') ''
   end if
end subroutine print_constrained_system_message

!> Print message for adjoint system solve
subroutine print_adjoint_message(unit, verbosity)
   integer, intent(in) :: unit, verbosity

   if (verbosity > 0) then
      write(unit, '(a)') 'Solving adjoint system: J*y = dfdq'
      write(unit, '(a)') ''
   end if
end subroutine print_adjoint_message

!> Print gradient calculation time
subroutine print_gradient_time(unit, verbosity, timer)
   integer, intent(in) :: unit, verbosity
   type(timer_type), intent(in) :: timer

   if (verbosity > 1) then
      write(unit, '(a, 1x, a)') "Gradient calculation time : ", &
         format_time(timer%get("gradient"))
      write(unit, '(a)') ''
   end if
end subroutine print_gradient_time

!> Print total solve time
subroutine print_total_time(unit, verbosity, timer)
   integer, intent(in) :: unit, verbosity
   type(timer_type), intent(in) :: timer

   if (verbosity > 1) then
      write(unit, '(a, 1x, a)') "Total solve time : ", format_time(timer%get("total"))
      write(unit, '(a)') ''
   end if
end subroutine print_total_time


end module multicharge_model_type