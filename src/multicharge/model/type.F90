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
   use solver_type, only: mchrg_solver_type, mchrg_solver_input
   
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
      !> Solver for the linear equations
      class(mchrg_solver_type), allocatable :: solver
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

subroutine solve(self, mol, slv, error, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL, &
   & energy, gradient, sigma, qvec, dqdr, dqdL, dfdq, verbose, new_unit)
   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in), target :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> The solver instance
   class(mchrg_solver_type), intent(in) :: slv
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
   real(wp), intent(in), optional :: dfdq(:)
   ! Optional print verbossity number input flag
   integer, intent(in), optional :: verbose
   ! Output unit
   integer, intent(in), optional :: new_unit

   integer :: ic, jc, iat, ndim
   logical :: grad, cpq, dcn
   logical :: add_lagr = .true.   

   ! Local variable for print verbosity
   integer :: verbose_solve
   integer :: out
   ! Variables for solving ES equation
   real(wp), allocatable :: xvec(:), vrhs(:), amat(:, :)
   real(wp), allocatable :: ainv(:, :), jmat(:, :)
   ! Gradients
   real(wp), allocatable :: dadr(:, :, :), dadL(:, :, :), atrace(:, :)
   real(wp), allocatable :: dxdr(:, :, :), dxdL(:, :, :)
   type(cache_container), allocatable :: cache
   real(wp), allocatable :: trans(:, :)

   ! Response vectors: vvec = electronegativity response (Jv=chi)
   ! uvec = constraint response (Ju=1)
   real(wp), allocatable :: chivec(:), unitvec(:)
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
   real(wp) :: yvecsum
   ! Adjoint vector p
   real(wp), allocatable :: padj(:)

   !> Timer
   type(timer_type) :: timer

   ! Calculate gradient if the respective arrays are present
   dcn = present(dcndr) .and. present(dcndL)
   grad = present(gradient) .and. present(sigma) .and. dcn
   cpq = present(dqdr) .and. present(dqdL) .and. dcn
   adj_grad = present(dfdq) .and. slv%need_pos_def .eqv. .true. 

   if (.not. present(verbose)) then
      verbose_solve = 0
   else
      verbose_solve = verbose
   end if 

   if (present(new_unit)) then
      out = new_unit
   else
      out = output_unit
   end if

   ! The cg_solver requires postive definite system,
   ! so the Lagrangian constraint is handled separately.
   if (slv%need_pos_def .eqv. .true.) then
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

   ! Setup X-vector and A^-1 for the linear system
   allocate(xvec(ndim))
   call self%get_xvec(mol, cache, xvec)

   call timer%pop ! stop setup timer

   ! Print header
   call print_solve_header(out, verbose_solve, timer)

   allocate(jmat(mol%nat, mol%nat))

   if (add_lagr .eqv. .true.) then

      vrhs = xvec
      ainv = amat

      ! Solving the linear system
      call slv%solve(amat, xvec, vrhs, ainv, cpq, new_unit=out, error=error)

      ! Partial charges
      if (present(qvec)) then
         qvec(:) = vrhs(:mol%nat)
      end if

      jmat = amat(:mol%nat, :mol%nat)
      ! Electrostatic energy
      if (present(energy)) then
         call symv(jmat, vrhs(:mol%nat), xvec(:mol%nat), &
            & alpha=0.5_wp, beta=-1.0_wp, uplo='l')
         energy(:) = energy(:) + vrhs(:mol%nat) * xvec(:mol%nat)
      end if

   else
      ! Constrained system
      allocate(vvec(mol%nat))  
      allocate(uvec(mol%nat))
      allocate(unitvec(mol%nat))
      allocate(chivec(mol%nat))
      allocate(jinv(mol%nat, mol%nat))

      ! J matrix and chi vector for the constrained system
      jmat = amat(:mol%nat, :mol%nat)
      chivec = -xvec(:mol%nat)

      ! Initial guess for vvec: v = chi / diag(J)
      !$omp parallel do default(none) shared(mol,jmat,uvec,vvec,chivec)
      do ic = 1, mol%nat
         uvec(ic)= 1.0_wp/jmat(ic, ic) + tiny(1.0_wp)
         vvec(ic)= chivec(ic)/jmat(ic, ic) + tiny(1.0_wp)
      end do
      !$omp end parallel do
      unitvec = 1.0_wp

      call print_constrained_system_message(out, verbose_solve, 'u')
      ! Constrained response: J*u = 1
      call slv%solve(amat=jmat, xvec=unitvec, vrhs=uvec, ainv=jinv, cpq=cpq, new_unit=out, error=error)

      call print_constrained_system_message(out, verbose_solve, 'v')
      ! Constrained response: J*u = chi
      call slv%solve(amat=jmat, xvec=chivec, vrhs=vvec, ainv=jinv, cpq=cpq, new_unit=out, error=error)

      uvecsum = sum(uvec)
      vvecsum = sum(vvec)

      ! Lagrangian multiplier
      lambda = - (mol%charge + vvecsum) / (uvecsum + tiny(1.0_wp))

      ! Reconstruct the full VRHS for gradient calculations
      if (allocated(vrhs)) deallocate(vrhs)
      allocate(vrhs(mol%nat+1))

      vrhs(:mol%nat) = -vvec - lambda * uvec
      vrhs(mol%nat+1) = lambda

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
   end if

   ! Allocate and get amat derivatives
   if (grad .or. cpq) then
      allocate(dadr(3, mol%nat, mol%nat), dadL(3, 3, mol%nat), atrace(3, mol%nat))
      allocate(dxdr(3, mol%nat, mol%nat), dxdL(3, 3, mol%nat))
      call self%get_xvec_derivs(mol, cache, dxdr, dxdL)
      call self%get_coulomb_derivs(mol, cache, vrhs, dadr, dadL, atrace)
      do iat = 1, mol%nat
         dadr(:, iat, iat) = atrace(:, iat) + dadr(:, iat, iat)
      end do
   end if

   ! Calculate gradients if requested
   if (grad) then
      call timer%push("gradient")
      call print_gradient_message(out, verbose_solve)

      if (adj_grad) then
         ! Adjoint gradient calculation
         if (add_lagr .eqv. .true.) then
            ! If the constraint response have not been calculated before
            if (allocated(uvec)) deallocate(uvec)
            allocate(uvec(mol%nat))
            if (allocated(unitvec)) deallocate(unitvec)
            allocate(unitvec(mol%nat))
            if (allocated(jinv)) deallocate(jinv)
            allocate(jinv(mol%nat, mol%nat))
            ! Initial guess for uvec: u = 1 / diag(J)
            do ic = 1, mol%nat
               uvec(ic)= 1.0_wp/jmat(ic, ic) + tiny(1.0_wp)
            end do
            unitvec = 1.0_wp
            ! Constrained response: J*u = 1
            call slv%solve(amat=jmat, xvec=unitvec, vrhs=uvec, &
               & ainv=jinv, cpq=cpq, new_unit=out, error=error)
            uvecsum = sum(uvec)
         end if 

         call print_adjoint_message(out, verbose_solve)

         ! Solving the adjoint system J*y = dfdq
         allocate(yvec(mol%nat))
         allocate(padj(mol%nat))

         ! Initial guess for yvec: y = dfdq / diag(J)
         do ic = 1, mol%nat
            yvec(ic)= dfdq(ic)/jmat(ic, ic) + tiny(1.0_wp)
         end do

         ! Derivative response: J*y = dfdq
         call slv%solve(amat=jmat, xvec=dfdq, vrhs=yvec, &
               & ainv=jinv, cpq=cpq, new_unit=out, error=error)
         yvecsum = sum(yvec)

         ! Project out the component of yvec along the constraint direction uvec
         padj = yvec - (yvecsum / (uvecsum + tiny(1.0_wp))) * uvec

         ! Gradient: df/dr = -padj^T * (dA/dr * q + dx/dr)
         gradient = 0.0_wp

         !$omp parallel do collapse(2) default(none) &
         !$omp shared(mol,gradient,dadr,dxdr,vrhs,padj) private(iat,ic) schedule(static)
         do iat = 1, mol%nat
            do ic = 1, mol%nat
               gradient(1,iat) = gradient(1,iat) &
                  - padj(ic) * (dadr(1,iat,ic) * vrhs(ic) + dxdr(1,iat,ic))
               gradient(2,iat) = gradient(2,iat) &
                  - padj(ic) * (dadr(2,iat,ic) * vrhs(ic) + dxdr(2,iat,ic))
               gradient(3,iat) = gradient(3,iat) &
                  - padj(ic) * (dadr(3,iat,ic) * vrhs(ic) + dxdr(3,iat,ic))
            end do
         end do
         !$omp end parallel do
         
      else
         ! If no input dfdq present
         gradient = 0.0_wp
         call gemv(dadr(:, :, :mol%nat), vrhs(:mol%nat), gradient, beta=1.0_wp, alpha=0.5_wp)
         call gemv(dxdr(:, :, :mol%nat), vrhs(:mol%nat), gradient, beta=1.0_wp, alpha=-1.0_wp)
         call gemv(dadL, vrhs, sigma, beta=1.0_wp, alpha=0.5_wp)
         call gemv(dxdL, vrhs, sigma, beta=1.0_wp, alpha=-1.0_wp)
      end if

      call timer%pop ! stop gradient timer
      call print_gradient_time(out, verbose_solve, timer)
   end if

   ! Calculate charge derivatives if requested
   if (cpq) then
      do iat = 1, mol%nat
         dadr(:, :, iat) = -dxdr(:, :, iat) + dadr(:, :, iat)
         dadL(:, :, iat) = -dxdL(:, :, iat) + dadL(:, :, iat)
      end do
      call gemm(dadr, ainv(:, :mol%nat), dqdr, alpha=-1.0_wp)
      call gemm(dadL, ainv(:, :mol%nat), dqdL, alpha=-1.0_wp)
   end if     
   
   call timer%pop ! stop total solve timer

   call print_total_time(out, verbose_solve, timer)

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
subroutine print_solve_header(unit, verbose, timer)
   integer, intent(in) :: unit, verbose
   type(timer_type), intent(in) :: timer

   if (verbose > 0) then
      write(unit, '(a)') ''
      write(unit, '(54("-"))')
      write(unit, '(a)') "             Charge equilibration solver            "
      write(unit, '(54("-"))')
      write(unit, '(a)') ''
      if (verbose > 1) then
         write(unit, '(a, 1x, a)') "A-matrix and X-vector formation time : ", &
            format_time(timer%get("setup"))
      end if
   end if
end subroutine print_solve_header

!> Print message for constrained system solves
subroutine print_constrained_system_message(unit, verbose, which)
   integer, intent(in) :: unit, verbose
   character, intent(in) :: which

   if (verbose > 0) then
      write(unit, '(a)') ''
      if (which == 'u') then
         write(unit, '(a)') 'Solving constrained system: J*u = 1'
      else if (which == 'v') then
         write(unit, '(a)') 'Solving constrained system: J*v = chi'
      end if
   end if
end subroutine print_constrained_system_message

!> Print message for gradient calculation start
subroutine print_gradient_message(unit, verbose)
   integer, intent(in) :: unit, verbose

   if (verbose > 0) then
      write(unit, '(a)') ''
      write(unit, '(a)') 'Calculating gradients'
   end if
end subroutine print_gradient_message

!> Print message for adjoint system solve
subroutine print_adjoint_message(unit, verbose)
   integer, intent(in) :: unit, verbose

   if (verbose > 0) then
      write(unit, '(a)') ''
      write(unit, '(a)') 'Solving adjoint system: J*y = dfdq'
   end if
end subroutine print_adjoint_message

!> Print gradient calculation time
subroutine print_gradient_time(unit, verbose, timer)
   integer, intent(in) :: unit, verbose
   type(timer_type), intent(in) :: timer

   if (verbose > 1) then
      write(unit, '(a)') ''
      write(unit, '(a, 1x, a)') "Gradient calculation time : ", &
         format_time(timer%get("gradient"))
      write(unit, '(a)') ''
   end if
end subroutine print_gradient_time

!> Print total solve time
subroutine print_total_time(unit, verbose, timer)
   integer, intent(in) :: unit, verbose
   type(timer_type), intent(in) :: timer

   if (verbose > 1) then
      write(unit, '(a, 1x, a)') "Total solve time : ", format_time(timer%get("total"))
      write(unit, '(a)') ''
      write(unit, '(a)') ''
   end if
end subroutine print_total_time


end module multicharge_model_type