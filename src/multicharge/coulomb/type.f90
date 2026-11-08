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

!> @file multicharge/coulomb/type.f90
!> Provides the base class for matrix-free Coulomb operators

!> Matrix-free representation of the Coulomb matrix of Gaussian charges.
!>
!> The operator represents A_ij = erf(gam_ij r_ij) / r_ij for i /= j, summed
!> over all lattice images for periodic systems, and eta_i + sqrt(2/pi)/rad_i
!> on the diagonal, with gam_ij = 1/sqrt(rad_i^2 + rad_j^2).
module multicharge_coulomb_type
   use mctc_env, only : error_type, wp
   use multicharge_solver_type, only : mchrg_operator_type
   implicit none
   private

   public :: coulomb_operator, coulomb_input, coulomb_method


   !> Available representations of the Coulomb matrix
   type :: enum_coulomb_method
      !> Dense Coulomb matrix
      integer :: dense = 0

      !> Particle mesh Ewald summation
      integer :: pme = 1
   end type enum_coulomb_method

   !> Actual enumerator of the Coulomb matrix representations
   type(enum_coulomb_method), parameter :: coulomb_method = enum_coulomb_method()

   !> Input configuration of the Coulomb matrix representation
   type :: coulomb_input
      !> Representation of the Coulomb matrix
      integer :: method = coulomb_method%dense

      !> Target accuracy of the operator, per pair term in Hartree/e^2
      real(wp) :: accuracy = 1.0e-7_wp

      !> Real-space cutoff of the particle mesh Ewald summation in Bohr
      real(wp), allocatable :: cutoff

      !> Maximal grid spacing of the particle mesh Ewald summation in Bohr
      real(wp), allocatable :: spacing

      !> Order of the cardinal B-splines of the particle mesh Ewald summation
      integer, allocatable :: order
   end type coulomb_input

   !> Abstract matrix-free Coulomb operator
   type, abstract, extends(mchrg_operator_type) :: coulomb_operator
      !> Number of atoms
      integer :: nat = 0
   contains
      !> Contract the operator derivatives with two charge-like vectors
      procedure(add_gradient), deferred :: add_gradient
      !> Write a summary of the operator setup
      procedure(write_info), deferred :: write_info
   end type coulomb_operator


   abstract interface
      !> Add scale * sum_ij p_i q_j dA_ij/dR to the gradient and the strain
      !> derivative to sigma
      subroutine add_gradient(self, pvec, qvec, gradient, sigma, scale, error)
         import :: coulomb_operator, wp, error_type

         !> Coulomb operator
         class(coulomb_operator), intent(in) :: self

         !> First charge-like vector
         real(wp), intent(in) :: pvec(:)

         !> Second charge-like vector
         real(wp), intent(in) :: qvec(:)

         !> Cartesian gradient
         real(wp), intent(inout) :: gradient(:, :)

         !> Strain derivative
         real(wp), intent(inout) :: sigma(:, :)

         !> Scaling factor of the contraction
         real(wp), intent(in) :: scale

         !> Error handling
         type(error_type), allocatable, intent(out) :: error

      end subroutine add_gradient

      !> Write a summary of the operator parameters
      subroutine write_info(self, unit)
         import :: coulomb_operator

         !> Coulomb operator
         class(coulomb_operator), intent(in) :: self

         !> Output unit
         integer, intent(in) :: unit

      end subroutine write_info
   end interface

end module multicharge_coulomb_type
