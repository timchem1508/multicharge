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
   use mctc_csrlist, only : csr_list
   use mctc_io, only : structure_type
   use mctc_wignerseitz, only : wignerseitz_cell
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
      !> neighborlist, the entry of atom j in the row of atom i holds dCN_j/dR_i
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
      !> neighborlist, the entry of atom j in the row of atom i holds dqloc_j/dR_i
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

      !> Inverse of amat on the pattern of the complete neighborlist, the
      !> charge constraint part -uvec*uvec^T/sum(uvec) is not stored
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

      !> Complete neighborlist storage of dabdr, the entry of a neighbor in the
      !> row of atom i holds the derivative of component i w.r.t. that neighbor
      real(wp), allocatable :: dabdrlist(:, :)

      !> Compressed-row pattern of the position derivatives of the charges, row j
      !> holds the atoms whose displacement changes the charge of atom j
      type(csr_list) :: dqdrpat

      !> Position derivatives of the charges on dqdrpat without the charge
      !> constraint term, the entry of atom i in row j holds dq_j/dR_i
      real(wp), allocatable :: dqdrlist(:, :)

      !> Factors of the charge constraint term of the position derivatives,
      !> dq_j/dR_i = dqdrlist(:, (j, i)) - dlambdadr(:, i) * uvec(j)
      real(wp), allocatable :: dlambdadr(:, :)

      !> Diagonal part of the second derivative of the Maxwell capacitance
      !> matrix w.r.t. positions of (3, 3, mol%nat) size.
      real(wp), allocatable :: d2cdr2diag(:, :, :)

      !> Sparse CSR storage of the hessian
      real(wp), allocatable :: hesslist(:, :, :)

      !> Logical flag for gradient calculation
      logical :: grad

      !> Logical flag for Hessian calculation
      logical :: hess
   end type mchrg_cache

end module multicharge_model_cache
