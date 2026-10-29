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
!> Half storage of pair arrays on sorted neighborlists

!> Sorted complete neighborlists and the half storage of pair arrays on them.
!> A sorted complete neighborlist keeps the atom itself first in each row,
!> followed by its neighbors in descending order of their index,
!>
!>   row i = i | neighbors j > i, descending | neighbors j < i, descending
!>
!> The first part of each row, the atom and its neighbors of higher index, is
!> the upper triangle of the list. Arrays of symmetric or antisymmetric pair
!> quantities only keep the upper triangle in the half storage, the upper
!> triangle of row i is held at hinl(i):hinl(i+1)-1 in the order of the row.
!> The pair of an atom with a neighbor of lower index is found in the upper
!> triangle of the row of the neighbor by a binary search. Arrays of
!> non-symmetric pair quantities keep both entries of a pair on the complete
!> list. The rows of a half list are its upper triangles, so the routines of
!> the half storage also take a half list, which suffices for the charges and
!> their energy gradient.
module multicharge_atomsort
   use mctc_env, only : wp, i8, error_type
   use mctc_io, only : structure_type
   use mctc_csrlist, only : csr_list, new_csr_list, distance_screening, &
      & sort_neighborlist, find_lower_start
   use mctc_wignerseitz, only : wignerseitz_cell
   implicit none
   private

   public :: new_sorted_list, get_half_offsets, get_half_list, find_neighbor, &
      & get_pair_positions, symv_half

   !> Product of a symmetric matrix in the half storage of a sorted complete
   !> neighborlist with a vector or with a vector for each atom
   interface symv_half
      module procedure :: symv_half_1d
      module procedure :: symv_half_2d
   end interface symv_half

contains

!> Create a complete neighborlist with sorted rows, the atom itself first and
!> its neighbors in descending order of their index
subroutine new_sorted_list(list, mol, error, cutoff, wsc)

   !> Sorted complete neighborlist
   type(csr_list), intent(out) :: list

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Realspace cutoff of the neighborlist
   real(wp), intent(in), optional :: cutoff

   !> Wigner-Seitz cell, receives the images of the entries of a periodic list
   type(wignerseitz_cell), intent(inout), optional :: wsc

   !> Distance screening, the threshold is set from the cutoff by the list
   type(distance_screening) :: screen

   call new_csr_list(list, screen, mol, error, wsc, cutoff, complete=.true.)
   if (allocated(error)) return
   call sort_neighborlist(list, wsc)

end subroutine new_sorted_list


!> Offsets of the upper triangles of the rows of a sorted complete
!> neighborlist in the half storage of symmetric pair arrays. The upper
!> triangle of row i, the atom itself and its neighbors of higher index, is
!> held at hinl(i):hinl(i+1)-1, hinl(nat+1) - 1 is the size of the storage.
!> The rows of a half list are its upper triangles, its offsets are kept.
subroutine get_half_offsets(list, hinl)

   !> Sorted complete neighborlist or half list
   type(csr_list), intent(in) :: list

   !> Offset of the upper triangle of each row in the half storage
   integer(i8), allocatable, intent(out) :: hinl(:)

   integer :: iat, nat

   if (.not. list%complete) then
      hinl = list%inl
      return
   end if

   nat = size(list%inl) - 1
   allocate(hinl(nat + 1))
   hinl(1) = 1_i8

   ! The upper triangle ends before the first neighbor of lower index
   !$omp parallel do default(none) schedule(static) shared(list, hinl, nat) &
   !$omp private(iat)
   do iat = 1, nat
      hinl(iat + 1) = find_lower_start(list, iat) - list%inl(iat)
   end do
   !$omp end parallel do

   do iat = 1, nat
      hinl(iat + 1) = hinl(iat) + hinl(iat + 1)
   end do

end subroutine get_half_offsets


!> Half neighborlist of the upper triangles of a sorted complete neighborlist,
!> whose entries follow the half storage of the complete list. It serves the
!> routines of half lists, e.g. the coordination numbers with their position
!> derivatives in the half storage.
subroutine get_half_list(list, hinl, half)

   !> Sorted complete neighborlist
   type(csr_list), intent(in) :: list

   !> Offset of the upper triangle of each row in the half storage
   integer(i8), intent(in) :: hinl(:)

   !> Half neighborlist of the upper triangles
   type(csr_list), intent(out) :: half

   integer :: iat, nat
   integer(i8) :: koff

   nat = size(list%inl) - 1
   half%complete = .false.
   if (allocated(list%cutoff)) half%cutoff = list%cutoff
   if (allocated(list%trans)) half%trans = list%trans
   half%inl = hinl(:nat + 1)
   allocate(half%nlat(hinl(nat + 1) - 1))
   if (allocated(list%nltr)) allocate(half%nltr(hinl(nat + 1) - 1))

   !$omp parallel do default(none) schedule(static) shared(list, hinl, half, nat) &
   !$omp private(iat, koff)
   do iat = 1, nat
      koff = list%inl(iat) - hinl(iat)
      half%nlat(hinl(iat):hinl(iat + 1) - 1) = &
         & list%nlat(hinl(iat) + koff:hinl(iat + 1) - 1 + koff)
      if (allocated(list%nltr)) half%nltr(hinl(iat):hinl(iat + 1) - 1) = &
         & list%nltr(hinl(iat) + koff:hinl(iat + 1) - 1 + koff)
   end do
   !$omp end parallel do

end subroutine get_half_list


!> Position of an atom in the row of another atom of a sorted complete
!> neighborlist by a binary search, zero if the atoms are not neighbors
pure function find_neighbor(list, iat, jat) result(kat)

   !> Sorted complete neighborlist
   type(csr_list), intent(in) :: list

   !> Atom of the row
   integer, intent(in) :: iat

   !> Atom to find in the row
   integer, intent(in) :: jat

   !> Position of the entry of jat in the row of iat
   integer(i8) :: kat

   integer(i8) :: base, hi, len, half

   kat = list%inl(iat)
   if (jat == iat) return

   ! First neighbor not above jat, the neighbors are in descending order. The
   ! search keeps the range [base, base + len] of this position and halves it
   ! without branches on the comparison.
   base = kat + 1
   hi = list%inl(iat + 1) - 1
   len = hi - base + 1
   kat = 0_i8
   if (len < 1) return
   do while (len > 1)
      half = len / 2
      if (list%nlat(base + half - 1) > jat) base = base + half
      len = len - half
   end do
   if (list%nlat(base) > jat) base = base + 1
   if (base <= hi) then
      if (list%nlat(base) == jat) kat = base
   end if

end function find_neighbor


!> Positions of the pairs of an atom with the entries of its row in a sorted
!> complete neighborlist, in the half storage of symmetric pair arrays and, if
!> requested, as the positions of the atom in the rows of its neighbors. The
!> pairs with neighbors of lower index are found by a binary search.
pure subroutine get_pair_positions(list, hinl, iat, hpos, tpos)

   !> Sorted complete neighborlist
   type(csr_list), intent(in) :: list

   !> Offset of the upper triangle of each row in the half storage
   integer(i8), intent(in) :: hinl(:)

   !> Atom of the row
   integer, intent(in) :: iat

   !> Position of the pair of each entry of the row in the half storage
   integer(i8), intent(out) :: hpos(:)

   !> Position of the atom in the row of each entry of the row
   integer(i8), intent(out), optional :: tpos(:)

   integer :: jat, ient
   integer(i8) :: ist, kat

   ist = list%inl(iat)
   do ient = 1, int(list%inl(iat + 1) - ist)
      jat = list%nlat(ist + ient - 1)
      if (jat >= iat) then
         hpos(ient) = hinl(iat) + ient - 1
         if (present(tpos)) tpos(ient) = find_neighbor(list, jat, iat)
      else
         kat = find_neighbor(list, jat, iat)
         hpos(ient) = hinl(jat) + kat - list%inl(jat)
         if (present(tpos)) tpos(ient) = kat
      end if
   end do

end subroutine get_pair_positions


!> Product of a symmetric matrix in the half storage of a sorted complete
!> neighborlist or a half list with a vector, y = A x. Each pair of an upper
!> triangle also adds the mirrored element to the row of its neighbor, which
!> the threads collect in private vectors.
subroutine symv_half_1d(list, hinl, amat, xvec, yvec)

   !> Sorted complete neighborlist or half list
   type(csr_list), intent(in) :: list

   !> Offset of the upper triangle of each row in the half storage
   integer(i8), intent(in) :: hinl(:)

   !> Upper triangle of the symmetric matrix in the half storage
   real(wp), intent(in) :: amat(:)

   !> Input vector
   real(wp), intent(in) :: xvec(:)

   !> Output vector, its elements of the atoms are overwritten
   real(wp), intent(inout) :: yvec(:)

   integer :: iat, jat, nat
   integer(i8) :: hat, koff
   real(wp) :: aij, yi
   real(wp), allocatable :: ypriv(:)

   nat = size(list%inl) - 1
   yvec(:nat) = 0.0_wp

   !$omp parallel default(none) shared(list, hinl, amat, xvec, yvec, nat) &
   !$omp private(iat, jat, hat, koff, aij, yi, ypriv)
   allocate(ypriv(nat), source=0.0_wp)
   !$omp do schedule(guided)
   do iat = 1, nat
      koff = list%inl(iat) - hinl(iat)
      yi = amat(hinl(iat)) * xvec(iat)
      do hat = hinl(iat) + 1, hinl(iat + 1) - 1
         jat = list%nlat(hat + koff)
         aij = amat(hat)
         yi = yi + aij * xvec(jat)
         ypriv(jat) = ypriv(jat) + aij * xvec(iat)
      end do
      ypriv(iat) = ypriv(iat) + yi
   end do
   !$omp end do
   !$omp critical (symv_half_1d_)
   yvec(:nat) = yvec(:nat) + ypriv
   !$omp end critical (symv_half_1d_)
   deallocate(ypriv)
   !$omp end parallel

end subroutine symv_half_1d


!> Product of a symmetric matrix in the half storage of a sorted complete
!> neighborlist or a half list with a vector for each atom,
!> y(:, i) = sum_j A_ij x(:, j)
subroutine symv_half_2d(list, hinl, amat, xvec, yvec)

   !> Sorted complete neighborlist or half list
   type(csr_list), intent(in) :: list

   !> Offset of the upper triangle of each row in the half storage
   integer(i8), intent(in) :: hinl(:)

   !> Upper triangle of the symmetric matrix in the half storage
   real(wp), intent(in) :: amat(:)

   !> Input vectors of the atoms
   real(wp), intent(in) :: xvec(:, :)

   !> Output vectors, the ones of the atoms are overwritten
   real(wp), intent(inout) :: yvec(:, :)

   integer :: iat, jat, nat, nv
   integer(i8) :: hat, koff
   real(wp) :: aij
   real(wp), allocatable :: ypriv(:, :), yi(:)

   nat = size(list%inl) - 1
   nv = size(xvec, 1)
   yvec(:, :nat) = 0.0_wp

   !$omp parallel default(none) shared(list, hinl, amat, xvec, yvec, nat, nv) &
   !$omp private(iat, jat, hat, koff, aij, yi, ypriv)
   allocate(ypriv(nv, nat), source=0.0_wp)
   allocate(yi(nv))
   !$omp do schedule(guided)
   do iat = 1, nat
      koff = list%inl(iat) - hinl(iat)
      yi(:) = amat(hinl(iat)) * xvec(:, iat)
      do hat = hinl(iat) + 1, hinl(iat + 1) - 1
         jat = list%nlat(hat + koff)
         aij = amat(hat)
         yi(:) = yi + aij * xvec(:, jat)
         ypriv(:, jat) = ypriv(:, jat) + aij * xvec(:, iat)
      end do
      ypriv(:, iat) = ypriv(:, iat) + yi
   end do
   !$omp end do
   !$omp critical (symv_half_2d_)
   yvec(:, :nat) = yvec(:, :nat) + ypriv
   !$omp end critical (symv_half_2d_)
   deallocate(ypriv, yi)
   !$omp end parallel

end subroutine symv_half_2d

end module multicharge_atomsort
