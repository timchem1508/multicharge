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

!> @file multicharge/atomsort.f90
!> Atom sorting utilities for multicharge calculations

module multicharge_atomsort
   use mctc_env, only : wp, i8
   use mctc_io, only : structure_type
   implicit none
   private

   public :: get_spatial_order

   !> Number of atoms in a group of the charge derivatives on a neighborlist,
   !> the last group of the spatial order may hold fewer atoms
   integer, parameter :: max_group = 8

contains

!> Order of the atoms from a recursive bisection of the structure, which is
!> independent of the input order of the atoms. Each max_group consecutive
!> atoms of the order form a group of the charge derivatives on a
!> neighborlist, whose atoms are close to each other and share most of their
!> neighbors.
subroutine get_spatial_order(mol, order)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Atoms in the spatial order, order(i) is the atom at position i
   integer, allocatable, intent(out) :: order(:)

   integer :: iat

   order = [(iat, iat = 1, mol%nat)]
   call bisect_atoms(mol%xyz, order)

end subroutine get_spatial_order


!> Split the atoms at the median of the longest side of their bounding box
!> until each part holds at most max_group atoms. The first part holds half
!> of the groups of max_group atoms, so only the last part of the order holds
!> fewer atoms.
recursive subroutine bisect_atoms(xyz, order)

   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)

   !> Atoms to split, reordered in place
   integer, intent(inout), contiguous :: order(:)

   integer :: iat, nleft, dir
   real(wp) :: lower(3), upper(3)

   if (size(order) <= max_group) return

   lower(:) = xyz(:, order(1))
   upper(:) = lower
   do iat = 2, size(order)
      lower(:) = min(lower, xyz(:, order(iat)))
      upper(:) = max(upper, xyz(:, order(iat)))
   end do
   dir = maxloc(upper - lower, 1)

   nleft = max_group * ((size(order) + max_group - 1) / max_group / 2)
   call select_atoms(xyz(dir, :), order, nleft)
   call bisect_atoms(xyz, order(:nleft))
   call bisect_atoms(xyz, order(nleft + 1:))

end subroutine bisect_atoms


!> Partial sort of the atoms by a coordinate with the selection algorithm of
!> Hoare, the first k atoms of the order have no larger coordinate than the
!> others
subroutine select_atoms(coord, order, k)

   !> Coordinate of each atom
   real(wp), intent(in) :: coord(:)

   !> Atoms to reorder
   integer, intent(inout), contiguous :: order(:)

   !> Number of atoms in the first part
   integer, intent(in) :: k

   integer :: lo, hi, i, j, tmp
   real(wp) :: pivot, first, last

   lo = 1
   hi = size(order)
   do while (lo < hi)
      ! Median of the first, k-th and last coordinate of the range
      first = coord(order(lo))
      last = coord(order(hi))
      pivot = max(min(first, last), min(max(first, last), coord(order(k))))
      i = lo
      j = hi
      do while (i <= j)
         do while (coord(order(i)) < pivot)
            i = i + 1
         end do
         do while (pivot < coord(order(j)))
            j = j - 1
         end do
         if (i <= j) then
            tmp = order(i)
            order(i) = order(j)
            order(j) = tmp
            i = i + 1
            j = j - 1
         end if
      end do
      ! The atoms between j and i have the coordinate of the pivot
      if (j < k) lo = i
      if (k < i) hi = j
   end do

end subroutine select_atoms

end module multicharge_atomsort