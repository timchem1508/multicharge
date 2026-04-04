! This file is part of mctc-lib.
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
   use mctc_ncoord, only: adjacency_list
   use multicharge_blas, only: gemv, symv, gemm
   use multicharge_blascomp, only: gemv_cmp
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
      !> Update model-dependent quantities and cache
      subroutine update(self, mol, cache, trans, grad, list)
         import :: mchrg_model_type, structure_type, mchrg_cache, adjacency_list, wp
         !> Multicharge model type
         class(mchrg_model_type), intent(in) :: self
         !> Structure type
         type(structure_type), intent(in) :: mol
         !> Multicharge neighbourlist type
         type(adjacency_list), intent(in), optional :: list
         !> Multicharge cache 
         !> Allocation: cn, qloc, wsc
         type(mchrg_cache), intent(inout) :: cache
         !> Lattice vectors
         real(wp), intent(in) :: trans(:, :)
         !> Flag to compute derivatives (dcndr, dcndL, dqlocdr, dqlocdL)
         logical, intent(in) :: grad
      end subroutine update

      !> Capacitance matrix construction using cached CN/charge data (only for the EEQBC model)
      subroutine get_capacitance_matrix(self, mol, ndim, cache, list)
         import :: mchrg_model_type, structure_type, mchrg_cache, adjacency_list, wp
         !> Multicharge model type
         class(mchrg_model_type), intent(in) :: self
         !> Structure type
         type(structure_type), intent(in) :: mol
         !> System size
         integer, intent(in) :: ndim
         !> Multicharge cache 
         !> Allocation: cmat, (dcdr, dcdL if cache%dcndr/L and cache%dqlocdr/L allocated)
         type(mchrg_cache), intent(inout) :: cache
         !> Multicharge neighbourlist type
         type(adjacency_list), intent(in), optional :: list
      end subroutine get_capacitance_matrix

      !> Coulomb interaction matrix (A-matrix) construction
      subroutine get_coulomb_matrix(self, mol, ndim, cache, list)
         import :: mchrg_model_type, structure_type, mchrg_cache, adjacency_list, wp
         !> Multicharge model type
         class(mchrg_model_type), intent(in) :: self
         !> Structure type
         type(structure_type), intent(in) :: mol
         !> System size
         integer, intent(in) :: ndim
         !> Multicharge cache 
         !> Allocation: amat
         type(mchrg_cache), intent(inout) :: cache
         !> Multicharge neighbourlist type
         type(adjacency_list), intent(in), optional :: list
      end subroutine get_coulomb_matrix


      !> Coulomb matrix derivatives contracted with charges
      subroutine get_coulomb_derivs(self, mol, ndim, cache, list)
         import :: mchrg_model_type, structure_type, mchrg_cache, adjacency_list, wp
         !> Multicharge model type
         class(mchrg_model_type), intent(in) :: self
         !> Structure type
         type(structure_type), intent(in) :: mol
         !> System size
         integer, intent(in) :: ndim
         !> Multicharge cache 
         !> Allocation: dadr, dadL
         type(mchrg_cache), intent(inout) :: cache
         !> Multicharge neighbourlist type
         type(adjacency_list), intent(in), optional :: list
      end subroutine get_coulomb_derivs

      !> Electronegativity vector construction

      subroutine get_xvec(self, mol, ndim, cache, list)
         import :: mchrg_model_type, mchrg_cache, structure_type, adjacency_list, wp
         !> Multicharge model type
         class(mchrg_model_type), intent(in) :: self
         !> Structure type
         type(structure_type), intent(in) :: mol
         !> System size
         integer, intent(in) :: ndim
         !> Multicharge cache 
         !> Allocation: xvec, xtmp
         type(mchrg_cache), intent(inout) :: cache
         !> Multicharge neighbourlist type
         type(adjacency_list), intent(in), optional :: list
      end subroutine get_xvec

      !> Derivatives of electronegativity vector
      subroutine get_xvec_derivs(self, mol, ndim, cache, list)
         import :: mchrg_model_type, structure_type, mchrg_cache, adjacency_list, wp
         !> Multicharge model type
         class(mchrg_model_type), intent(in) :: self
         !> Structure type
         type(structure_type), intent(in) :: mol
         !> System size
         integer, intent(in) :: ndim
         !> Multicharge cache 
         !> Allocation: dxdr, dxdL
         type(mchrg_cache), intent(inout) :: cache
         !> Multicharge neighbourlist type
         type(adjacency_list), intent(in), optional :: list
      end subroutine get_xvec_derivs    
   end interface

   real(wp), parameter :: twopi = 2 * pi
   real(wp), parameter :: eps = tiny(1.0_wp)
   real(wp), parameter :: geom_tol = 1.0e-10_wp

contains

!> Generate direct lattice translation vectors within a supercell of 2×2×2 repetitions.
subroutine get_dir_trans(lattice, trans)
   !> Lattice parameters (3×3 matrix)
   real(wp), intent(in) :: lattice(:, :)
   !> Output translation vectors (3 × N) where N = 2×2×2 = 8
   real(wp), allocatable, intent(out) :: trans(:, :)
   integer, parameter :: rep(3) = [2, 2, 2]

   call get_lattice_points(lattice, rep, .true., trans)

end subroutine get_dir_trans

!> Generate reciprocal lattice translation vectors within a supercell of 2×2×2 repetitions.
subroutine get_rec_trans(lattice, trans)
   !> Lattice parameters (3×3 matrix)
   real(wp), intent(in) :: lattice(:, :)
   !> Output translation vectors in reciprocal space (3 × N) where N = 2×2×2 = 8
   real(wp), allocatable, intent(out) :: trans(:, :)
   integer, parameter :: rep(3) = [2, 2, 2]
   real(wp) :: rec_lat(3, 3)

   rec_lat = twopi * transpose(matinv_3x3(lattice))
   call get_lattice_points(rec_lat, rep, .false., trans)

end subroutine get_rec_trans

!> Top-level solve routine with optional persistent cache
subroutine solve(self, mol, solver, cache, error, &
   & energy, gradient, sigma, qvec, dqdr, dqdL, list, verbosity, unit)
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
   !> Neighbour list optional type
   type(adjacency_list), intent(in), optional :: list
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

   call self%get_capacitance_matrix(mol, ndim, cache, list)

   ! Setup the Coulomb matrix 
   call self%get_coulomb_matrix(mol, ndim, cache, list)

   ! Get RHS of ES equation
   call self%get_xvec(mol, ndim, cache, list)
   if (.not. allocated(cache%vrhs)) then
      allocate(cache%vrhs(mol%nat + 1))
   end if

   ! pop setup timer
   call timer%pop 

   ! Print header
   call print_solve_header(print_unit, verbosity_solve, timer%get("setup"))

   if (add_lagr .eqv. .true.) then
      if (.not. allocated(cache%ainv)) then
         allocate(cache%ainv(ndim, ndim))
      end if
      cache%vrhs = cache%xvec
      cache%ainv = cache%amat
      call solver%solve(amat=cache%amat, xvec=cache%xvec, vrhs=cache%vrhs, ainv=cache%ainv, &
         & cpq=cpq, new_unit=print_unit, error=error)

   else
      if (.not. allocated(cache%uvec)) then
         allocate(cache%uvec(mol%nat))  
      end if
      allocate(unitvec(mol%nat))
      allocate(vvec(mol%nat))
      ! Initial guess
      if (present(list)) then
         cache%uvec(:) = 1.0_wp / (cache%adiag(:) + eps)
         vvec(:) = - cache%xvec(:) / (cache%adiag(:) + eps)
      else
         do iat = 1, mol%nat
            cache%uvec(iat) = 1.0_wp / (cache%amat(iat, iat) + eps)
            vvec(iat) = - cache%xvec(iat) / (cache%amat(iat, iat) + eps)
         end do
      end if

      unitvec = 1.0_wp

      call print_constrained_system_message(print_unit, verbosity_solve, 'u')
      ! Constrained response: A*uvec = 1
      call solver%solve(amat=cache%amat, alist=cache%alist, adiag=cache%adiag, xvec=unitvec, &
            & vrhs=cache%uvec, list=list, new_unit=print_unit, error=error)
      call print_constrained_system_message(print_unit, verbosity_solve, 'v')
      ! Constrained response: A*uvec = chi
      call solver%solve(amat=cache%amat, alist=cache%alist, adiag=cache%adiag, xvec=-cache%xvec, &
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
         call gemv_cmp(list, cache%alist, cache%adiag, cache%vrhs, cache%xvec(:mol%nat), &
            & alpha=0.5_wp, beta=-1.0_wp)
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

   ! Allocate and get amat derivatives
   if (dcn) then
      call timer%push("setup_gradient")
      call self%get_xvec_derivs(mol, ndim, cache, list)
      call self%get_coulomb_derivs(mol, ndim, cache, list)
      allocate(daqxdr(3, mol%nat, ndim), source=0.0_wp)
      allocate(daqxdL(3, 3, ndim), source=0.0_wp)
      ! pop gradient setup
      call timer%pop
      call print_gradient_header(print_unit, verbosity_solve, timer%get("setup_gradient"))
   end if

   ! Calculate gradients if requested
   if (grad) then
      call timer%push("gradient") 
      do iat = 1, mol%nat
         daqxdr(:, :, iat) = - cache%dxdr(:, :, iat) + 0.5_wp * cache%dadr(:, :, iat)
         daqxdL(:, :, iat) = - cache%dxdL(:, :, iat) + 0.5_wp * cache%dadL(:, :, iat)
      end do
      call gemv(daqxdr(:, :, :mol%nat), cache%vrhs(:mol%nat), gradient, beta=1.0_wp, alpha=1.0_wp)
      call gemv(daqxdL, cache%vrhs, sigma, beta=1.0_wp, alpha=1.0_wp)
      ! pop gradient timer
      call timer%pop 
      call print_gradient_time(print_unit, verbosity_solve, timer%get("gradient"))
   end if

   ! Calculate charge derivatives if requested
   if (cpq) then
      call timer%push("cpq")
      do iat = 1, mol%nat
         daqxdr(:, :, iat) = cache%dxdr(:, :, iat) - cache%dadr(:, :, iat)
         daqxdL(:, :, iat) = cache%dxdL(:, :, iat) - cache%dadL(:, :, iat)
      end do
      call gemm(daqxdr, cache%ainv(:, :mol%nat), dqdr, alpha=1.0_wp)
      call gemm(daqxdL, cache%ainv(:, :mol%nat), dqdL, alpha=1.0_wp)
      ! pop cpq timer
      call timer%pop 
      call print_gradient_time(print_unit, verbosity_solve, timer%get("cpq"))
   end if

   ! pop total solve timer
   call timer%pop 
   call print_total_time(print_unit, verbosity_solve, timer%get("total"))

end subroutine solve

!> Adjoint external gradient calculation using cached data 
!
!> This routine evaluates dF/dR and dF/dL from the derivative of the
!> objective w.r.t. charges (dF/dq), avoiding explicit differentiation
!> of the charge solution by solving an adjoint system.
subroutine get_external_gradient(self, mol, solver, cache, error, dfdq, dfdr, dfdL, list, unit, verbosity)
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
   !> Neighbour list optional type
   type(adjacency_list), intent(in), optional :: list
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

   call timer%push("setup_external")
   call timer%pop

   call print_gradient_header(print_unit, verbosity_solve, timer%get("setup_external"))
   call timer%push("external_gradient")

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
      call solver%solve(amat=cache%amat, alist=cache%alist, adiag=cache%adiag, &
         & xvec=dfdq, vrhs=yvec, list=list, error=error)
      if (allocated(error)) return

      ! Projection of uvec on yvec
      yvecsum = sum(yvec)
      uvecsum = sum(cache%uvec)
      scale = yvecsum / (uvecsum + eps)
      allocate(padj(mol%nat))
      padj = yvec - scale * cache%uvec
      
      ! Evaluate external gradients via adjoint contraction:
      ! dfdr = p^T * (db/dr - dA/dr X q)
      call gemv(daqxdr, padj, dfdr, alpha=1.0_wp, beta=0.0_wp)
      ! dfdL = p^T * (db/dL - dA/dL X q)
      call gemv(daqxdL, padj, dfdL, alpha=1.0_wp, beta=0.0_wp)

   else
      allocate(padj(ndim))

      ! Direct solution: J*yvec = dfdq
      call print_adjoint_message(print_unit, verbosity_solve)
      call solver%solve(amat=cache%amat, xvec=dfdq_loc, vrhs=padj, new_unit=print_unit, error=error)
      if (allocated(error)) return

      ! Evaluate external gradients via adjoint contraction:
      ! dfdr = p^T * (db/dr - dA/dr X q)
      call gemv(daqxdr, padj, dfdr, alpha=1.0_wp, beta=0.0_wp)
      ! dfdL = p^T * (db/dL - dA/dL X q)
      call gemv(daqxdL, padj, dfdL, alpha=1.0_wp, beta=0.0_wp)

   end if

   ! pop dfdr timer
   call timer%pop
   call print_gradient_time(print_unit, verbosity_solve, timer%get("external_gradient"))

end subroutine get_external_gradient

!> Local charges calculation
subroutine local_charge(self, mol, trans, qloc, dqlocdr, dqlocdL, &
   & list, dqlocdrlist, dqlocdrdiag)
   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: trans(:, :)
   !> Local atomic partial charges
   real(wp), intent(out) :: qloc(:)
   !> Optional derivative of local atomic partial charges w.r.t. atomic positions
   real(wp), intent(out), optional :: dqlocdr(3, mol%nat, mol%nat)
   !> Optional derivative of local atomic partial charges w.r.t. lattice vectors
   real(wp), intent(out), optional :: dqlocdL(3, 3, mol%nat)
   !> Lattice points
   type(adjacency_list), intent(in), optional :: list
   !> Optional derivative of local atomic partial charges w.r.t. atomic positions
   real(wp), intent(out), optional :: dqlocdrlist(:, :), dqlocdrdiag(:, :)

   qloc = 0.0_wp
   if (present(dqlocdr) .and. present(dqlocdL)) then
      dqlocdr = 0.0_wp
      dqlocdL = 0.0_wp
   end if
   if (present(list) .and. present(dqlocdrlist) .and. present(dqlocdrdiag) .and. present(dqlocdL)) then
      dqlocdrlist = 0.0_wp
      dqlocdrdiag = 0.0_wp
      dqlocdL = 0.0_wp
   end if
   ! Get the electronegativity weighted CN for local charge
   ! Derivatives depend only in this CN
   if (allocated(self%ncoord_en)) then
      call self%ncoord_en%get_coordination_number(mol, trans, qloc, dcndr=dqlocdr, &
         & dcndrlist=dqlocdrlist, dcndrdiag=dqlocdrdiag, dcndL=dqlocdL, list=list)
   end if

   ! Distribute the total charge equally
   qloc = qloc + mol%charge / real(mol%nat, wp)

end subroutine local_charge

!> Print header for charge equilibration solver
subroutine print_solve_header(unit, verbosity, timer)
   integer, intent(in) :: unit, verbosity
   real(wp):: timer

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
   integer, intent(in) :: unit, verbosity
   real(wp):: timer

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
   real(wp):: timer

   if (verbosity > 1) then
      write(unit, '(a, 1x, a)') "Gradient calculation time : ", format_time(timer)
      write(unit, '(a)') ''
   end if
end subroutine print_gradient_time

!> Print gradient calculation time
subroutine print_energy_time(unit, verbosity, timer)
   integer, intent(in) :: unit, verbosity
   real(wp):: timer

   if (verbosity > 1) then
      write(unit, '(a, 1x, a)') "Energy calculation time : ", format_time(timer)
      write(unit, '(a)') ''
   end if
end subroutine print_energy_time

!> Print total solve time
subroutine print_total_time(unit, verbosity, timer)
   integer, intent(in) :: unit, verbosity
   real(wp):: timer

   if (verbosity > 1) then
      write(unit, '(a, 1x, a)') "Total solve time : ", format_time(timer)
      write(unit, '(a)') ''
   end if
end subroutine print_total_time

end module multicharge_model_type