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

!> @file multicharge/cache.f90
!> Contains the cache baseclass for the charge models and a container for mutable cache data

!> Cache for charge models
module multicharge_model_cache
   use mctc_env, only: wp
   use mctc_io, only: structure_type
   use multicharge_wignerseitz, only: wignerseitz_cell_type
   implicit none
   private

   !> Cache for the charge model
   type, public :: mchrg_cache
      !> Coordination number array
      real(wp), allocatable :: cn(:)
      !> Ewald separation parameter
      real(wp) :: alpha
      !> Wigner-Seitz cell
      type(wignerseitz_cell_type) :: wsc
      !> Translation matrix
      real(wp), allocatable :: trans(:, :)
      !> Local charges
      real(wp), allocatable :: qloc(:)
      !> Full Maxwell capacitance matrix
      real(wp), allocatable :: cmat(:, :)
      !> Compressed version of the C-matrix
      real(wp), allocatable :: clist(:)
      !> Store tmp array from xvec calculation for reuse
      real(wp), allocatable :: xtmp(:)
      !> Electronegativity vector
      real(wp), allocatable :: xvec(:)
      !> Coulomb matrix
      real(wp), allocatable :: amat(:, :)
      !> Compressed version of the A-matrix
      real(wp), allocatable :: alist(:)
      !> Inversed amat
      real(wp), allocatable :: ainv(:,:)
      !> Solution of the ES equation
      real(wp), allocatable :: vrhs(:)
      !> Constraint response: jamt*uvec=1
      real(wp), allocatable :: uvec(:)
      !> Coordination number gradient w.r.t the positions
      real(wp), allocatable :: dcndr(:, :, :)
      !> Coordination number gradient w.r.t the lattice vectors
      real(wp), allocatable :: dcndL(:, :, :)
      !> Local charge derivatives w.r.t positions
      real(wp), allocatable :: dqlocdr(:, :, :)
      !> Local charge derivatives w.r.t lattice vectors
      real(wp), allocatable :: dqlocdL(:, :, :)
      !> Derivative of Maxwell capacitance matrix w.r.t positions
      real(wp), allocatable :: dcdr(:, :, :)
      !> Derivative of Maxwell capacitance matrix w.r.t positions in compressed format
      real(wp), allocatable :: dcdrdiag(:, :)
      !> Derivative of Maxwell capacitance matrix w.r.t lattice vectors
      real(wp), allocatable :: dcdL(:, :, :)
      !> Coulomb matrix derivatives w.r.t positions
      real(wp), allocatable :: dadr(:, :, :)
      !> Coulomb matrix derivatives w.r.t lattice vectors
      real(wp), allocatable :: dadL(:, :, :)
      !> Electronegativity derivatives w.r.t positions
      real(wp), allocatable :: dxdr(:, :, :)
      !> Electronegativity derivatives w.r.t lattice vectors
      real(wp), allocatable :: dxdL(:, :, :)
      !> Logical flag for gradient calculation
      logical :: grad

   end type mchrg_cache

end module multicharge_model_cache
