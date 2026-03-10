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
   use multicharge_lapack, only: sytrf, sytrs, sytri
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
   ! Optional derivative of the electrostatic energy w.r.t. atomic partial charges
   real(wp), intent(in), contiguous, optional :: dfdq(:)
   ! Optional print verbossity number input flag
   integer, intent(in), optional :: verbosity
   ! Output unit
   integer, intent(in), optional :: unit

   integer :: ic, jc, iat, jat, ndim
   logical :: grad, cpq, dcn
   logical :: add_lagr = .true.   

   ! Local variable for print verbosity
   integer :: verbosity_solve
   integer :: print_unit
   ! Variables for solving ES equation
   real(wp), allocatable :: xvec(:), vrhs(:), amat(:, :)
   real(wp), allocatable :: ainv(:, :), jmat(:, :)
   ! Gradients
   ! dadr, dadl and atrace already includes multiplication by q
   real(wp), allocatable :: dadr(:, :, :), dadL(:, :, :), atrace(:, :)
   real(wp), allocatable :: dxdr(:, :, :), dxdL(:, :, :)
   type(cache_container), allocatable :: cache
   real(wp), allocatable :: trans(:, :)

   ! Response vectors: vvec = electronegativity response (Jv=chi)
   ! uvec = constraint response (Ju=1)
   real(wp), allocatable :: unitvec(:)
   real(wp), allocatable :: vvec(:), uvec(:)
   ! Sums of the v and u vector elements
   real(wp) :: uvecsum, vvecsum
   real(wp), allocatable :: jinv(:, :)
   ! Lagrangian factor for constraint
   real(wp) :: lambda 

   ! Gradient solver
   logical :: adj_grad
   ! Derivative response J*y=dfdq
   real(wp), allocatable :: yvec(:)
   real(wp) :: yvecsum, factor
   ! Adjoint vector p
   real(wp), allocatable :: padj(:)
   ! -dA/dr*q + db/dr :: dxdr - dadr 
   real(wp), allocatable :: derivsum(:,:,:)

   ! Temporary arrays for dqdr and dqdL solution
   real(wp), allocatable :: rhs(:), sol(:), diag(:)
   real(wp) :: scale

   ! Timer
   type(timer_type) :: timer

   ! Calculate gradient if the respective arrays are present
   dcn = present(dcndr) .and. present(dcndL)
   grad = present(gradient) .and. present(sigma) .and. dcn .and. .not. present(dfdq)
   cpq = present(dqdr) .and. present(dqdL) .and. dcn .and. .not. present(dfdq)
   adj_grad = present(dfdq)

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

   ! The cg_solver requires postive definite system,
   ! so the Lagrangian constraint is handled separately.
   if (solver%need_pos_def .eqv. .true.) then
      ndim = mol%nat 
      add_lagr = .false.
   else
      ndim = mol%nat + 1
      add_lagr = .true.
   end if

   ! Update cache
   allocate(cache)
   call self%update(mol, cache, ndim, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL)

   call timer%push("total")
   call timer%push("setup")
   ! Setup the Coulomb matrix 
   allocate(amat(ndim, ndim))
   call self%get_coulomb_matrix(mol, cache, amat)

   ! Get RHS of ES equation
   allocate(xvec(ndim))
   call self%get_xvec(mol, cache, xvec)

   ! stop setup timer
   call timer%pop 

   ! Print header
   call print_solve_header(print_unit, verbosity_solve, timer)

   allocate(jmat(mol%nat, mol%nat))

   if (add_lagr .eqv. .true.) then

      vrhs = xvec
      ainv = amat

      ! Solving the linear system
      call solver%solve(amat, xvec, vrhs, ainv=ainv, cpq=cpq, new_unit=print_unit, error=error)

      jmat = amat(:mol%nat, :mol%nat)

   else
      ! Constrained system
      allocate(vvec(mol%nat))  
      allocate(uvec(mol%nat))
      allocate(unitvec(mol%nat))

      ! Initial guess for vvec: v = chi / diag(J)
      do ic = 1, mol%nat
         uvec(ic)= 1.0_wp/(amat(ic, ic) + eps)
         vvec(ic)= -xvec(ic)/(amat(ic, ic) + eps)
      end do
      unitvec = 1.0_wp
      ainv = amat

      call print_constrained_system_message(print_unit, verbosity_solve, 'u')
      ! Constrained response: J*u = 1
      call solver%solve(amat=amat, xvec=unitvec, vrhs=uvec, new_unit=print_unit, error=error)

      call print_constrained_system_message(print_unit, verbosity_solve, 'v')
      ! Constrained response: J*u = chi
      call solver%solve(amat=amat, xvec=-xvec, vrhs=vvec, new_unit=print_unit, error=error)

      uvecsum = sum(uvec)
      vvecsum = sum(vvec)

      ! Lagrangian multiplier
      lambda = - (mol%charge + vvecsum) / (uvecsum + eps)

      ! Reconstruct the full VRHS for gradient calculations
      allocate(vrhs(mol%nat+1))

      vrhs(:mol%nat) = -vvec - lambda * uvec
      vrhs(mol%nat+1) = lambda

      jmat = amat(:mol%nat, :mol%nat)

   end if

   ! Partial charges
   if (present(qvec)) then
      qvec(:) = vrhs(:mol%nat)
   end if

   ! Electrostatic energy
   if (present(energy)) then
      call symv(jmat, vrhs(:mol%nat), xvec(:mol%nat), &
         & alpha=0.5_wp, beta=-1.0_wp, uplo='l')
      energy(:) = energy(:) + vrhs(:mol%nat) * xvec(:mol%nat)
   end if

   ! Allocate and get amat derivatives
   if (grad .or. cpq .or. adj_grad) then
      allocate(dadr(3, mol%nat, mol%nat), dadL(3, 3, mol%nat), atrace(3, mol%nat))
      allocate(dxdr(3, mol%nat, mol%nat), dxdL(3, 3, mol%nat))
      call self%get_xvec_derivs(mol, cache, dxdr, dxdL)
      call self%get_coulomb_derivs(mol, cache, vrhs, dadr, dadL, atrace)
      do iat = 1, mol%nat
         dadr(:, iat, iat) = atrace(:, iat) + dadr(:, iat, iat)
      end do
   end if

   if (adj_grad) then
      ! Adjoint gradient calculation
      call print_adjoint_message(print_unit, verbosity_solve)
      call timer%push("gradient") 
      allocate(yvec(mol%nat))
      allocate(padj(mol%nat))
      allocate(derivsum(3, mol%nat, mol%nat), source=0.0_wp)

      if (solver%need_pos_def .eqv. .false.) then
         ! Constrained response: J*u = 1
         allocate(uvec(mol%nat))
         allocate(unitvec(mol%nat))
         unitvec = 1.0_wp
         call solver%solve(amat=jmat, xvec=unitvec, vrhs=uvec, &
            & new_unit=print_unit, error=error)
         uvecsum = sum(uvec)

         yvec = dfdq
         ! Derivative response: J*y = dfdq
         call solver%solve(amat=jmat, xvec=dfdq, vrhs=yvec, &
               & new_unit=print_unit, error=error)
         yvecsum = sum(yvec)
         ! Project out the component of yvec along the constraint direction uvec
         padj = yvec - (yvecsum / (uvecsum + eps)) * uvec

         ! Gradient: df/dr = padj^T * (- dA/dr * q + dx/dr)
         do iat = 1, mol%nat
            derivsum(:, :, iat) = dxdr(:, :, iat) - dadr(:, :, iat)
         end do
         call gemv(derivsum(:, :, :mol%nat), padj(:mol%nat), gradient, beta=0.0_wp, alpha=1.0_wp)

      else 
         ! Initial guess for yvec: y = dfdq / diag(J)
         do ic = 1, mol%nat
            yvec(ic)= dfdq(ic)/(amat(ic, ic) + eps)
         end do

         ! Derivative response: J*y = dfdq
         call solver%solve(amat=amat, xvec=dfdq, vrhs=yvec, &
               & new_unit=print_unit, error=error)
         yvecsum = sum(yvec)

         ! Project out the component of yvec along the constraint direction uvec
         padj = yvec - (yvecsum / (uvecsum + eps)) * uvec

         ! Gradient: df/dr = padj^T * (- dA/dr * q + dx/dr)
         do iat = 1, mol%nat
            derivsum(:, :, iat) = dxdr(:, :, iat) - dadr(:, :, iat)
         end do
         call gemv(derivsum(:, :, :mol%nat), padj(:mol%nat), gradient, beta=0.0_wp, alpha=1.0_wp)

      end if

      ! Stop gradient timer
      call timer%pop 
      call print_gradient_time(print_unit, verbosity_solve, timer)
   end if 

   ! Calculate gradients if requested
   if (grad) then
      call timer%push("gradient") 

      call gemv(dadr(:, :, :mol%nat), vrhs(:mol%nat), gradient, beta=1.0_wp, alpha=0.5_wp)
      call gemv(dxdr(:, :, :mol%nat), vrhs(:mol%nat), gradient, beta=1.0_wp, alpha=-1.0_wp)
      call gemv(dadL, vrhs, sigma, beta=1.0_wp, alpha=0.5_wp)
      call gemv(dxdL, vrhs, sigma, beta=1.0_wp, alpha=-1.0_wp)

      ! stop gradient timer
      call timer%pop 
      call print_gradient_time(print_unit, verbosity_solve, timer)
   end if

   ! Calculate charge derivatives if requested
   if (cpq) then
      call timer%push("gradient") 
      if (solver%need_pos_def) then
         ! Diagonal for initial guess
         allocate(diag(mol%nat))
         do iat = 1, mol%nat
            diag(iat) = amat(iat, iat)
         end do

         !$omp parallel default(none) &
         !$omp shared(mol, solver, dqdr, dqdL, amat, dadr, dxdr, dadL, dxdL, &
         !$omp        uvec, uvecsum, diag, print_unit, error) &
         !$omp private(iat, ic, jc, rhs, sol, scale)
         allocate(rhs(mol%nat), sol(mol%nat))
         ! Position derivatives (dq/dR) 
         !$omp do schedule(runtime)
         do iat = 1, mol%nat
            do ic = 1, 3
               ! RHS = db/dR - dA/dR * q
               rhs(:) = dxdr(ic, iat, :) - dadr(ic, iat, :)
               ! Initial guess: rhs / diag
               sol(:) = rhs(:) / diag(:)
               ! Solve J * m = rhs
               call solver%solve(amat=amat, xvec=rhs, vrhs=sol, &
                  & new_unit=print_unit, error=error)
               ! Projection factor
               scale = sum(sol) / uvecsum
               ! dq/dR = m - scale * u
               dqdr(ic, iat, :) = sol(:) - scale * uvec(:)
            end do
         end do
         !$omp end do
         ! Charge virial (dq/dL) 
         !$omp do schedule(runtime)
         do ic = 1, 3
            do jc = 1, 3
               rhs(:) = dxdL(ic, jc, :) - dadL(ic, jc, :)
               sol(:) = rhs(:) / diag(:)
               call solver%solve(amat=amat, xvec=rhs, vrhs=sol, &
                  & new_unit=print_unit, error=error)
               scale = sum(sol) / uvecsum
               dqdL(ic, jc, :) = sol(:) - scale * uvec(:)
            end do
         end do
         !$omp end do
         deallocate(rhs, sol)
         !$omp end parallel
      else
         ! Original inverse‑based method for augmented matrix
         do iat = 1, mol%nat
            dadr(:, :, iat) = -dxdr(:, :, iat) + dadr(:, :, iat)
            dadL(:, :, iat) = -dxdL(:, :, iat) + dadL(:, :, iat)
         end do
         call gemm(dadr, ainv(:, :mol%nat), dqdr, alpha=-1.0_wp)
         call gemm(dadL, ainv(:, :mol%nat), dqdL, alpha=-1.0_wp)
      end if
      
      ! Stop gradient timer
      call timer%pop 
      call print_gradient_time(print_unit, verbosity_solve, timer)
   end if

   ! stop total solve timer
   
   call timer%pop 

   call print_total_time(print_unit, verbosity_solve, timer)

end subroutine solve

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
      write(unit, '(a)') "             Charge equilibration solver            "
      write(unit, '(54("-"))')
      write(unit, '(a)') ''
      if (verbosity > 1) then
         write(unit, '(a, 1x, a)') "Setup time : ", &
            format_time(timer%get("setup"))
         write(unit, '(a)') ''
      end if
   end if
end subroutine print_solve_header

!> Print message for constrained system solves
subroutine print_constrained_system_message(unit, verbosity, which)
   integer, intent(in) :: unit, verbosity
   character, intent(in) :: which

   if (verbosity > 0) then
      if (which == 'u') then
         write(unit, '(a)') 'Solving constrained system: J*u = 1'
      else if (which == 'v') then
         write(unit, '(a)') 'Solving constrained system: J*v = chi'
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