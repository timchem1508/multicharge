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
   use mctc_env, only: timer_type, format_time, error_type, fatal_error, wp, ik => IK
   use mctc_io, only: structure_type
   use mctc_io_constants, only: pi
   use mctc_io_math, only: matinv_3x3
   use mctc_cutoff, only: get_lattice_points
   use mctc_ncoord, only: ncoord_type
   use multicharge_blas, only: gemv, symv, gemm
   use multicharge_lapack, only: sytrf, sytrs
   use multicharge_wignerseitz, only: wignerseitz_cell_type, new_wignerseitz_cell
   use multicharge_model_cache, only: mchrg_cache
   use multicharge_solver_type, only: mchrg_solver_type
   
   implicit none
   private

   public :: mchrg_model_type, get_dir_trans, get_rec_trans

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
      real(wp), allocatable :: kcnrad
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
      !> Update cache
      procedure(update), deferred :: update
      !> Calculate capacitance matrix
      procedure(get_capacitance_matrix), deferred :: get_capacitance_matrix
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
      subroutine update(self, mol, cache, trans, dcndr, dcndL)
         import :: mchrg_model_type, structure_type, mchrg_cache, wp
         class(mchrg_model_type), intent(in) :: self
         type(structure_type), intent(in) :: mol
         type(mchrg_cache), intent(inout) :: cache
         real(wp), intent(in) :: trans(:, :)
         real(wp), intent(inout), contiguous, optional :: dcndr(:, :, :)
         real(wp), intent(inout), contiguous, optional :: dcndL(:, :, :)
      end subroutine update

      subroutine get_capacitance_matrix(self, mol, ndim, cache)
         import :: mchrg_model_type, structure_type, mchrg_cache, wp
         class(mchrg_model_type), intent(in) :: self
         type(structure_type), intent(in) :: mol
         integer, intent(in) :: ndim
         type(mchrg_cache), intent(inout) :: cache
      end subroutine get_capacitance_matrix

      subroutine get_coulomb_matrix(self, mol, ndim, cache)
         import :: mchrg_model_type, structure_type, mchrg_cache, wp
         class(mchrg_model_type), intent(in) :: self
         type(structure_type), intent(in) :: mol
         integer, intent(in) :: ndim
         type(mchrg_cache), intent(inout) :: cache
      end subroutine get_coulomb_matrix

      subroutine get_coulomb_derivs(self, mol, ndim, cache)
         import :: mchrg_model_type, structure_type, mchrg_cache, wp
         class(mchrg_model_type), intent(in) :: self
         type(structure_type), intent(in) :: mol
         integer, intent(in) :: ndim
         type(mchrg_cache), intent(inout) :: cache
      end subroutine get_coulomb_derivs

      subroutine get_xvec(self, mol, ndim, cache)
         import :: mchrg_model_type, mchrg_cache, structure_type, wp
         class(mchrg_model_type), intent(in) :: self
         type(structure_type), intent(in) :: mol
         integer, intent(in) :: ndim
         type(mchrg_cache), intent(inout) :: cache
      end subroutine get_xvec

      subroutine get_xvec_derivs(self, mol, ndim, cache)
         import :: mchrg_model_type, structure_type, mchrg_cache, wp
         class(mchrg_model_type), intent(in) :: self
         type(structure_type), intent(in) :: mol
         integer, intent(in) :: ndim
         type(mchrg_cache), intent(inout) :: cache
      end subroutine get_xvec_derivs    
   end interface

   real(wp), parameter :: twopi = 2 * pi
   real(wp), parameter :: eps = tiny(1.0_wp)
   real(wp), parameter :: geom_tol = 1.0e-10_wp

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

!> Top-level solve routine with optional persistent cache
subroutine solve(self, mol, solver, cache, error, &
   & energy, gradient, sigma, qvec, dqdr, dqdL, verbosity, unit)
   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in):: self
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
   !> Optional print verbossity number input flag
   integer, intent(in), optional :: verbosity
   !> Output unit
   integer, intent(in), optional :: unit

   integer :: iat, ndim

   real(wp), allocatable :: unitvec(:)
   real(wp), allocatable :: vvec(:)
   real(wp) :: uvecsum
   real(wp) :: vvecsum
   real(wp) :: lambda 
   real(wp), allocatable :: daqxdr(:,:,:)
   real(wp), allocatable :: daqxdL(:,:,:)

   logical :: grad, cpq, dcn
   logical :: add_lagr = .true.  
   type(timer_type) :: timer
   integer :: print_unit, verbosity_solve


   ! Calculate gradient if the respective arrays are present
   dcn = allocated(cache%dcndr) .and. allocated(cache%dcndL)
   grad = present(gradient) .and. present(sigma) .and. dcn
   cpq = present(dqdr) .and. present(dqdL) .and. dcn 

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

   ! The cg_solver requires postive definite system
   if (solver%need_pos_def .eqv. .true.) then
      ndim = mol%nat 
      add_lagr = .false.
   else
      ndim = mol%nat + 1
      add_lagr = .true.
   end if

   call timer%push("total")
   call timer%push("setup")

   call self%get_capacitance_matrix(mol, ndim, cache)

   ! Setup the Coulomb matrix 
   call self%get_coulomb_matrix(mol, ndim, cache)

   ! Get RHS of ES equation
   call self%get_xvec(mol, ndim, cache)
   if (.not. allocated(cache%vrhs)) then
      allocate(cache%vrhs(mol%nat + 1))
   end if

   ! pop setup timer
   call timer%pop 

   ! Print header
   call print_solve_header(print_unit, verbosity_solve, timer)

   if (add_lagr .eqv. .true.) then
      if (.not. allocated(cache%ainv)) then
         allocate(cache%ainv(ndim, ndim))
      end if
      cache%vrhs = cache%xvec
      cache%ainv = cache%amat
      call solver%solve(cache%amat, cache%xvec, cache%vrhs, ainv=cache%ainv, &
         & cpq=cpq, new_unit=print_unit, error=error)

   else
      if (.not. allocated(cache%uvec)) then
         allocate(cache%uvec(mol%nat))  
      end if
      allocate(unitvec(mol%nat))
      allocate(vvec(mol%nat))
      ! Initial guess
      do iat = 1, mol%nat
         cache%uvec(iat) = 1.0_wp / (cache%amat(iat, iat) + eps)
         vvec(iat) = - cache%xvec(iat) / (cache%amat(iat, iat) + eps)
      end do

      unitvec = 1.0_wp

      call print_constrained_system_message(print_unit, verbosity_solve, 'u')
      ! Constrained response: A*uvec = 1
      call solver%solve(amat=cache%amat, xvec=unitvec, vrhs=cache%uvec, new_unit=print_unit, error=error)
      call print_constrained_system_message(print_unit, verbosity_solve, 'v')
      ! Constrained response: A*uvec = chi
      call solver%solve(amat=cache%amat, xvec=-cache%xvec, vrhs=vvec, new_unit=print_unit, error=error)
      uvecsum = sum(cache%uvec)
      vvecsum = sum(vvec)
      ! Lagrangian multiplier
      lambda = - (mol%charge + vvecsum) / (uvecsum + eps)

      ! Reconstruct the full VRHS and JMAT
      cache%vrhs(:mol%nat) = -vvec - lambda * cache%uvec
      cache%vrhs(mol%nat + 1) = lambda

   end if

   ! Partial charges if present
   if (present(qvec)) then
      qvec(:) = cache%vrhs(:mol%nat)
   end if

   ! Electrostatic energy if present
   if (present(energy)) then
      call symv(cache%amat, cache%vrhs, cache%xvec(:mol%nat), &
         & alpha=0.5_wp, beta=-1.0_wp, uplo='l')
      if (ndim > mol%nat) then
         ! Reconstructing xvec for a proper energy calculation
         cache%xvec(:mol%nat) = cache%xvec(:mol%nat) - 0.5_wp * cache%vrhs(mol%nat + 1)
      end if
      energy(:) = energy(:) + cache%vrhs(:mol%nat) * cache%xvec(:mol%nat)
   end if

   ! Allocate and get amat derivatives
   if (dcn) then
      call timer%push("setup")
      call self%get_xvec_derivs(mol, ndim, cache)
      call self%get_coulomb_derivs(mol, ndim, cache)
      ! pop gradient setup
      call timer%pop
      call print_gradient_header(print_unit, verbosity_solve, timer)
   end if

   ! Calculate gradients if requested
   if (grad) then
      call timer%push("gradient") 
      allocate(daqxdr(3, mol%nat, ndim), source=0.0_wp)
      allocate(daqxdL(3, 3, ndim), source=0.0_wp)
      do iat = 1, mol%nat
         daqxdr(:, :, iat) = - cache%dxdr(:, :, iat) + 0.5_wp * cache%dadr(:, :, iat)
         daqxdL(:, :, iat) = - cache%dxdL(:, :, iat) + 0.5_wp * cache%dadL(:, :, iat)
      end do
      call gemv(daqxdr(:, :, :mol%nat), cache%vrhs(:mol%nat), gradient, beta=1.0_wp, alpha=1.0_wp)
      call gemv(daqxdL, cache%vrhs, sigma, beta=1.0_wp, alpha=1.0_wp)
      ! pop gradient timer
      call timer%pop 
      call print_gradient_time(print_unit, verbosity_solve, timer)
   end if

   ! Calculate charge derivatives if requested
   if (cpq) then
      call timer%push("gradient")
      allocate(daqxdr(3, mol%nat, ndim), source=0.0_wp)
      allocate(daqxdL(3, 3, ndim), source=0.0_wp)
      do iat = 1, mol%nat
         daqxdr(:, :, iat) = cache%dxdr(:, :, iat) - cache%dadr(:, :, iat)
         daqxdL(:, :, iat) = cache%dxdL(:, :, iat) - cache%dadL(:, :, iat)
      end do
      call gemm(daqxdr, cache%ainv(:, :mol%nat), dqdr, alpha=1.0_wp)
      call gemm(daqxdL, cache%ainv(:, :mol%nat), dqdL, alpha=1.0_wp)
      ! pop gradient timer
      call timer%pop 
      call print_gradient_time(print_unit, verbosity_solve, timer)
   end if

   ! pop total solve timer
   call timer%pop 
   call print_total_time(print_unit, verbosity_solve, timer)

end subroutine solve

!> Adjoint gradient calculation using cached data
subroutine get_external_gradient(self, mol, solver, cache, error, dfdq, dfdr, dfdL, unit, verbosity)
   !> Electronegativity equilibration model
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
   real(wp), intent(inout) :: dfdL(:,:)
   !> Output unit
   integer, intent(in), optional :: unit
   !> Verbosity level
   integer, intent(in), optional :: verbosity

   integer :: iat
   integer :: ndim
   real(wp), allocatable :: daqxdr(:,:,:)
   real(wp), allocatable :: daqxdL(:,:,:)
   real(wp), allocatable :: yvec(:)
   real(wp), allocatable :: unitvec(:)
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

   call timer%push("setup")
   call timer%pop

   call print_gradient_header(print_unit, verbosity_solve, timer)
   call timer%push("gradient")

   ! Get variables from the model cache
   if (allocated(cache%dadr) .and. allocated(cache%dadL) &
   & .and. allocated(cache%dxdr) .and. allocated(cache%dxdL)) then
      allocate(daqxdr(3, mol%nat, mol%nat), &
         & source=cache%dxdr(:, :, :mol%nat) - cache%dadr(:, :, :mol%nat))
      allocate(daqxdL(3, 3, mol%nat), &
         & source=cache%dxdL(:, :, :mol%nat) - cache%dadL(:, :, :mol%nat))
   else
      call fatal_error(error, "J-matrix and electronegativity derivatives are not allocated")
      return
   end if
   if (.not. allocated(cache%amat)) then
      call fatal_error(error, "J-matrix is not allocated")
      return
   end if
   if (.not. allocated(cache%uvec) .and. solver%need_pos_def) then
      call fatal_error(error, "Constraint response J*uvec = 1 is not allocated")
      return
   end if

   if (solver%need_pos_def) then
      allocate(yvec(ndim))
      do iat = 1, mol%nat
         yvec(iat) = dfdq(iat) / cache%amat(iat, iat)
      end do

      ! Constrained response: J*yvec = dfdq
      call print_adjoint_message(print_unit, verbosity_solve)
      call solver%solve(cache%amat, dfdq, yvec, error=error)
      if (allocated(error)) return

      ! Projection of 
      yvecsum = sum(yvec)
      uvecsum = sum(cache%uvec)
      scale = yvecsum / (uvecsum + eps)
      allocate(padj(mol%nat))
      padj = yvec - scale * cache%uvec
      
      ! Calculation of external derivatives using the adjoint method
      ! dfdr = p^T * (db/dr - dA/dr X q)
      call gemv(daqxdr, padj, dfdr, alpha=1.0_wp, beta=0.0_wp)
      ! dfdL = p^T * (db/dL - dA/dL X q)
      call gemv(daqxdL, padj, dfdL, alpha=1.0_wp, beta=0.0_wp)

   else
      allocate(padj(ndim))

      ! Direct solution: J*yvec = dfdq
      call print_adjoint_message(print_unit, verbosity_solve)
      call solver%solve(cache%amat, dfdq_loc, padj, new_unit=print_unit, error=error)
      if (allocated(error)) return

      ! Calculation of external derivatives using the adjoint method
      ! dfdr = p^T * (db/dr - dA/dr X q)
      call gemv(daqxdr, padj, dfdr, alpha=1.0_wp, beta=0.0_wp)
      ! dfdL = p^T * (db/dL - dA/dL X q)
      call gemv(daqxdL, padj, dfdL, alpha=1.0_wp, beta=0.0_wp)

   end if

   ! pop dfdr
   call timer%pop
   call print_gradient_time(print_unit, verbosity_solve, timer)

end subroutine get_external_gradient

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
subroutine print_gradient_header(unit, verbosity, timer)
   integer, intent(in) :: unit, verbosity
   type(timer_type), intent(in) :: timer

   if (verbosity > 0) then
      write(unit, '(54("-"))')
      write(unit, '(17x, a)') "Gradient Calculations"
      write(unit, '(54("-"))')
      write(unit, '(a)') ''
      if (verbosity > 1) then
         write(unit, '(a, 1x, a)') "Setup time : ", &
            format_time(timer%get("setup"))
         write(unit, '(a)') ''
      end if
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