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

!> @file multicharge/model/cache.f90
!> Contains the cache baseclass for the charge models and a container for
!> mutable cache data

!> Cache for charge models
module multicharge_model_cache
   use mctc_env, only : wp, i8
   use mctc_io, only : structure_type
   use mctc_wignerseitz, only : wignerseitz_cell
   use multicharge_lowerlist, only : lower_list
   implicit none
   private

   public :: mchrg_cache

   !> Cache for the charge model
   type :: mchrg_cache

      !> Coordination number array
      real(wp), allocatable :: cn(:)

      !> Coordination number derivatives w.r.t. positions
      real(wp), allocatable :: dcndr(:, :, :)

      !> Coordination number derivatives w.r.t. lattice vectors
      real(wp), allocatable :: dcndL(:, :, :)

      !> Coordination number derivatives w.r.t. positions on the pattern of the
      !> neighborlist, the entry of atom j in the row of atom i holds dCN_j/dR_i.
      !> Without a maximum CN cutoff dCN_i/dR_j = -dCN_j/dR_i, so a half list
      !> only stores the upper triangle.
      real(wp), allocatable :: dcndrlist(:, :)

      !> Ewald separation parameter
      real(wp) :: alpha

      !> Wigner-Seitz cell
      type(wignerseitz_cell) :: wsc

      !> Translation matrix
      real(wp), allocatable :: trans(:, :)

      !> Local charges
      real(wp), allocatable :: qloc(:)

      !> Local charge derivatives w.r.t. positions
      real(wp), allocatable :: dqlocdr(:, :, :)

      !> Local charge derivatives w.r.t. lattice vectors
      real(wp), allocatable :: dqlocdL(:, :, :)

      !> Local charge derivatives w.r.t. positions on the pattern of the
      !> neighborlist, the entry of atom j in the row of atom i holds dqloc_j/dR_i.
      !> Without a maximum CN cutoff dqloc_i/dR_j = dqloc_j/dR_i, so a half list
      !> only stores the upper triangle.
      real(wp), allocatable :: dqlocdrlist(:, :)

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
      real(wp), allocatable :: ainv(:, :)

      !> Inverse of amat on the pattern of the neighborlist, the charge
      !> constraint part -uvec*uvec^T/sum(uvec) is not stored. A half list only
      !> stores the upper triangle of the symmetric inverse.
      real(wp), allocatable :: ainvlist(:)

      !> Solution of the ES equation
      real(wp), allocatable :: vrhs(:)

      !> Constraint response: jamt*uvec=1
      real(wp), allocatable :: uvec(:)

      !> Derivative of Maxwell capacitance matrix w.r.t positions
      real(wp), allocatable :: dcdr(:, :, :)

      !> Derivative of Maxwell capacitance matrix w.r.t positions in compressed format
      real(wp), allocatable :: dcdrdiag(:, :)

      !> Derivative of Maxwell capacitance matrix w.r.t lattice vectors
      real(wp), allocatable :: dcdL(:, :, :)

      !> Linear combination alpha * dA/dR*q + beta * dX/dR w.r.t. positions,
      !> dabdr(:, x, b) is the derivative of component b w.r.t. atom x
      real(wp), allocatable :: dabdr(:, :, :)

      !> Linear combination alpha * dA/dL*q + beta * dX/dL w.r.t. strain
      real(wp), allocatable :: dabdL(:, :, :)

      !> Diagonal and upper triangle of dabdr on a half neighborlist, the entry
      !> of atom x in the row of atom i holds the derivative of component i
      !> w.r.t. atom x
      real(wp), allocatable :: dabdrlist_ut(:, :)

      !> Lower triangle of dabdr on a half neighborlist, the entry of atom x in
      !> the row of atom i holds the derivative of component x w.r.t. atom i,
      !> the diagonal entries are not used
      real(wp), allocatable :: dabdrlist_lt(:, :)

      !> Coulomb matrix derivatives contracted with an external vector,
      !> alpha * dA/dR*p w.r.t. positions, dafdr(:, x, b) is the derivative of
      !> component b w.r.t. atom x
      real(wp), allocatable :: dafdr(:, :, :)

      !> Diagonal and upper triangle of dafdr on a half neighborlist, the entry
      !> of atom x in the row of atom i holds the derivative of component i
      !> w.r.t. atom x
      real(wp), allocatable :: dafdrlist_ut(:, :)

      !> Lower triangle of dafdr on a half neighborlist, the entry of atom x in
      !> the row of atom i holds the derivative of component x w.r.t. atom i,
      !> the diagonal entries are not used
      real(wp), allocatable :: dafdrlist_lt(:, :)

      !> Lower triangle of the half neighborlist of the charge derivatives and
      !> the Hessian, row x holds the atoms b < x with the position of the entry
      !> of x in row b
      type(lower_list) :: lowlist

      !> Diagonal and upper triangle of the position derivatives of the charges
      !> on the pattern of the half neighborlist of the solve without the charge
      !> constraint term, the entry of atom i in row j holds dq_j/dR_i
      real(wp), allocatable :: dqdrlist_ut(:, :)

      !> Lower triangle of the position derivatives of the charges on the
      !> pattern of the half neighborlist of the solve without the charge
      !> constraint term, the entry of atom i in row j holds dq_i/dR_j, the
      !> diagonal entries are not used
      real(wp), allocatable :: dqdrlist_lt(:, :)

      !> Factors of the charge constraint term of the position derivatives,
      !> dq_j/dR_i = dqdrlist_ut(:, (j, i)) - dlambdadr(:, i) * uvec(j) for i >= j
      real(wp), allocatable :: dlambdadr(:, :)

      !> Diagonal part of the second derivative of the Maxwell capacitance
      !> matrix w.r.t. positions of (3, 3, mol%nat) size.
      real(wp), allocatable :: d2cdr2diag(:, :, :)

      !> Logical flag for gradient calculation
      logical :: grad

      !> Logical flag for Hessian calculation
      logical :: hess
   end type mchrg_cache

end module multicharge_model_cache
