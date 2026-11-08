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

!> @file multicharge/solver/type.f90
!> Provides a general base class for the linear system solvers

!> Base types for the linear system solvers
module multicharge_solver_type
   use mctc_csrlist, only : csr_list
   use mctc_env, only : error_type, wp
   implicit none
   private

   public :: mchrg_solver_type, mchrg_solver_input, mchrg_operator_type


   !> Abstract linear operator for matrix-free iterative solvers
   type, abstract :: mchrg_operator_type
   contains

      !> Apply the operator to a vector
      procedure(apply_operator), deferred :: apply

      !> Diagonal of the operator used for preconditioning
      procedure(get_operator_diagonal), deferred :: get_diagonal

   end type mchrg_operator_type


   !> Abstract base type for multicharge solvers
   type, abstract :: mchrg_solver_type

      !> Whether the solver requires a positive-definite matrix
      logical :: need_pos_def = .false.

      !> Whether the solver enforces the total charge by projection
      logical :: projected = .false.
   contains

      !> Solve the linear system
      procedure(solve), deferred :: solve

   end type mchrg_solver_type


   abstract interface
      !> Product of the linear operator with a vector
      subroutine apply_operator(self, xvec, yvec)
         import :: mchrg_operator_type, wp

         !> Linear operator
         class(mchrg_operator_type), intent(in) :: self

         !> Input vector
         real(wp), intent(in) :: xvec(:)

         !> Operator-vector product
         real(wp), intent(out) :: yvec(:)

      end subroutine apply_operator

      !> Diagonal elements of the linear operator
      subroutine get_operator_diagonal(self, diag)
         import :: mchrg_operator_type, wp

         !> Linear operator
         class(mchrg_operator_type), intent(in) :: self

         !> Diagonal of the operator
         real(wp), intent(out) :: diag(:)

      end subroutine get_operator_diagonal

      !> Solve the linear system of equations
      subroutine solve(self, amat, xvec, vrhs, ainv, alist, list, op, qtot, &
         & lambda, new_unit, error)
         import :: mchrg_solver_type, mchrg_operator_type, error_type, wp, csr_list

         !> Solver instance
         class(mchrg_solver_type), intent(in) :: self

         !> Dense coefficient matrix of the linear system
         real(wp), intent(in), optional :: amat(:, :)

         !> Right-hand side vector
         real(wp), intent(in) :: xvec(:)

         !> On input: initial guess; on output: solution
         real(wp), intent(inout), contiguous :: vrhs(:)

         !> Inverse matrix, requested by passing the argument
         real(wp), intent(out), optional :: ainv(:, :)

         !> Coefficient matrix values in the sparsity pattern of the CSR list
         real(wp), intent(in), optional :: alist(:)

         !> CSR neighbour list, selects the sparse representation of the matrix
         type(csr_list), intent(in), optional :: list

         !> Matrix-free linear operator, selects the operator representation
         class(mchrg_operator_type), intent(in), optional :: op

         !> Total sum of the solution, solves the system with a Lagrange
         !> multiplier enforcing the constraint by projection
         real(wp), intent(in), optional :: qtot

         !> Lagrange multiplier of the projected system
         real(wp), intent(out), optional :: lambda

         !> Output unit
         integer, intent(in), optional :: new_unit

         !> Error handling
         type(error_type), allocatable, intent(out) :: error

      end subroutine solve
   end interface


   !> Abstract base type for solver configuration
   type, abstract :: mchrg_solver_input
   end type mchrg_solver_input


end module multicharge_solver_type
