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
   implicit none
   private

   public :: mchrg_cache, clear_cache_derivs, clear_solve_derivs

   !> Cache for the charge model
   type :: mchrg_cache

      !> Coordination number array
      real(wp), allocatable :: cn(:)

      !> Coordination number derivatives w.r.t. positions
      real(wp), allocatable :: dcndr(:, :, :)

      !> Coordination number derivatives w.r.t. lattice vectors
      real(wp), allocatable :: dcndL(:, :, :)

      !> Coordination number derivatives
      real(wp), allocatable :: dcndrlist(:, :)

      !> Column of each entry of the list of the update in dcndrlist and
      !> dqlocdrlist, 0 for entries without CN or local charge derivatives
      integer, allocatable :: cnmap(:)

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

      !> Local charge derivatives on the list
      real(wp), allocatable :: dqlocdrlist(:, :)

      !> Row offsets of the list of the update, which holds the pair arrays clist,
      !> alist, ainvlist, dcndrlist and dqlocdrlist. Row i of a sorted half list is
      !> the upper triangle of row i of the sorted complete list.
      integer(i8), allocatable :: hinl(:)

      !> Transposed positions of the entries of the sorted complete neighborlist
      integer(i8), allocatable :: tpos(:)

      !> Offsets of the rows of the CN and local-charge derivative entries on the list
      integer(i8), allocatable :: cninl(:)

      !> Positions of the CN and local-charge derivative entries in the sorted
      !> complete list
      integer(i8), allocatable :: cnpos(:)

      !> Full Maxwell capacitance matrix
      real(wp), allocatable :: cmat(:, :)

      !> Capacitance matrix on the upper triangle of the list of the update
      real(wp), allocatable :: clist(:)

      !> Store tmp array from xvec calculation for reuse
      real(wp), allocatable :: xtmp(:)

      !> Electronegativity vector
      real(wp), allocatable :: xvec(:)

      !> Coulomb matrix
      real(wp), allocatable :: amat(:, :)

      !> Coulomb matrix on the upper triangle of the list of the update
      real(wp), allocatable :: alist(:)

      !> Inversed amat
      real(wp), allocatable :: ainv(:, :)

      !> Inverse of amat on the upper triangle of the list of the update
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

      !> Partial derivatives alpha * dA/dR*q + beta * dX/dR
      real(wp), allocatable :: dabdr(:, :, :)

      !> Linear combination alpha * dA/dL*q + beta * dX/dL w.r.t. strain
      real(wp), allocatable :: dabdL(:, :, :)

      !> Partial derivatives dabdr on the sorted complete neighborlist
      real(wp), allocatable :: dabdrlist(:, :)

      !> Coulomb matrix derivatives contracted with an external vector
      real(wp), allocatable :: dafdr(:, :, :)

      !> Strain derivative of the Coulomb matrix contracted with an external
      !> vector, set together with dafdr
      real(wp), allocatable :: dafdL(:, :, :)

      !> Derivatives dafdr on the sorted complete neighborlist
      real(wp), allocatable :: dafdrlist(:, :)

      !> Position derivatives of the charges on the neighborlist
      real(wp), allocatable :: dqdrlist(:, :)

      !> Charge constraint factors of the position derivatives
      real(wp), allocatable :: dlambdadr(:, :)

      !> Diagonal blocks of the second derivative of the capacitance matrix
      real(wp), allocatable :: d2cdr2diag(:, :, :)

      !> External electric field of the last solve
      real(wp), allocatable :: efield(:)
   end type mchrg_cache


contains


!> Drop all derivatives of the previous structure
subroutine clear_cache_derivs(cache)

   !> Multicharge cache
   type(mchrg_cache), intent(inout) :: cache

   call clear_solve_derivs(cache)

   if (allocated(cache%dcndr)) deallocate(cache%dcndr)
   if (allocated(cache%dcndL)) deallocate(cache%dcndL)
   if (allocated(cache%dcndrlist)) deallocate(cache%dcndrlist)
   if (allocated(cache%cnmap)) deallocate(cache%cnmap)
   if (allocated(cache%dqlocdr)) deallocate(cache%dqlocdr)
   if (allocated(cache%dqlocdL)) deallocate(cache%dqlocdL)
   if (allocated(cache%dqlocdrlist)) deallocate(cache%dqlocdrlist)
   if (allocated(cache%tpos)) deallocate(cache%tpos)
   if (allocated(cache%cninl)) deallocate(cache%cninl)
   if (allocated(cache%cnpos)) deallocate(cache%cnpos)
   if (allocated(cache%dcdr)) deallocate(cache%dcdr)
   if (allocated(cache%dcdrdiag)) deallocate(cache%dcdrdiag)
   if (allocated(cache%dcdL)) deallocate(cache%dcdL)
   if (allocated(cache%d2cdr2diag)) deallocate(cache%d2cdr2diag)

end subroutine clear_cache_derivs


!> Drop the inverse and derivatives of the previous solve
subroutine clear_solve_derivs(cache)

   !> Multicharge cache
   type(mchrg_cache), intent(inout) :: cache

   if (allocated(cache%ainv)) deallocate(cache%ainv)
   if (allocated(cache%ainvlist)) deallocate(cache%ainvlist)
   if (allocated(cache%dabdr)) deallocate(cache%dabdr)
   if (allocated(cache%dabdL)) deallocate(cache%dabdL)
   if (allocated(cache%dabdrlist)) deallocate(cache%dabdrlist)
   if (allocated(cache%dafdr)) deallocate(cache%dafdr)
   if (allocated(cache%dafdL)) deallocate(cache%dafdL)
   if (allocated(cache%dafdrlist)) deallocate(cache%dafdrlist)
   if (allocated(cache%dqdrlist)) deallocate(cache%dqdrlist)
   if (allocated(cache%dlambdadr)) deallocate(cache%dlambdadr)

end subroutine clear_solve_derivs


end module multicharge_model_cache
