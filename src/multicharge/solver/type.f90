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
!> Provides a general base class for the linear system solvers.

module multicharge_solver_type
   use mctc_env, only: error_type, wp
   use mctc_csrlist, only: csr_list
   implicit none
   private

   public :: mchrg_solver_type, mchrg_solver_input

   !> Abstract base type for Multi-charge solvers
   type, abstract :: mchrg_solver_type
      !> Type of matrix availiable for solver
      !> CG solver can use only positive definite matrices
      !> Direct solver can use either type of matrix
      logical, allocatable :: need_pos_def
   contains
      procedure(solve), deferred :: solve
   end type mchrg_solver_type


   abstract interface
      subroutine solve(self, amat, alist, xvec, vrhs, ainv, cpq, list, new_unit, error)
         import :: mchrg_solver_type, error_type, wp, csr_list
         class(mchrg_solver_type), intent(in) :: self
         !> A matrix of Ax=b system
         real(wp), intent(in), optional  :: amat(:, :)
         !> Off-diagonall elements of matrix for in compressed
         real(wp), intent(in), optional  :: alist(:)
         !> Right-hand side vector
         real(wp), intent(in)  :: xvec(:)
         !> On input: initial guess; on output: solution
         real(wp), intent(inout), contiguous :: vrhs(:)
         !> Inverse matrix
         real(wp), intent(out), optional :: ainv(:, :)
         !> Flag for coupled-perturbed equations
         logical, intent(in), optional :: cpq
         !> Neighbour list optional type
         type(csr_list), intent(in), optional :: list
         !> Output unit
         integer, intent(in), optional :: new_unit
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine solve
   end interface

   !> Solver input abstract type
   type, abstract, public :: mchrg_solver_input
   end type mchrg_solver_input



end module multicharge_solver_type
