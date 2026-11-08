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

!> @dir multicharge/coulomb
!> Contains the matrix-free representations of the Coulomb matrix

!> @file multicharge/coulomb.f90
!> Provides a reexport of the matrix-free Coulomb operators

!> Proxy module to construct the matrix-free Coulomb operators
module multicharge_coulomb
   use mctc_env, only : error_type, fatal_error, wp
   use mctc_io, only : structure_type
   use multicharge_coulomb_pme, only : pme_operator, new_pme_operator
   use multicharge_coulomb_type, only : coulomb_operator, coulomb_input, &
      & coulomb_method
   implicit none
   private

   public :: coulomb_operator, coulomb_input, coulomb_method
   public :: pme_operator
   public :: new_coulomb_operator


contains


!> Construct a matrix-free Coulomb operator from its input configuration
subroutine new_coulomb_operator(op, mol, error, input, rad, eta, alpha_ref)
   !> Coulomb operator to construct
   class(coulomb_operator), allocatable, intent(out) :: op
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Input configuration
   type(coulomb_input), intent(in) :: input
   !> Species-wise charge widths
   real(wp), intent(in) :: rad(:)
   !> Species-wise chemical hardness
   real(wp), intent(in) :: eta(:)
   !> Ewald parameter of the reference dense matrix in periodic systems
   real(wp), intent(in) :: alpha_ref

   if (.not. input%accuracy > 0.0_wp) then
      call fatal_error(error, "Accuracy of the Coulomb operator must be positive")
      return
   end if

   select case(input%method)
   case(coulomb_method%pme)
      block
         type(pme_operator), allocatable :: tmp
         allocate(tmp)
         call new_pme_operator(tmp, mol, error, input, rad, eta, alpha_ref)
         if (allocated(error)) return
         call move_alloc(tmp, op)
      end block
   case default
      call fatal_error(error, "The dense Coulomb matrix has no matrix-free operator")
   end select

end subroutine new_coulomb_operator

end module multicharge_coulomb
